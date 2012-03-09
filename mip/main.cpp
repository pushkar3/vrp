#include "gurobi_c++.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
#include <map>
#include <string>
using namespace std;

class Node {
public:
	int value;
	int robot;
	string name;
	vector<string> neighbor;

	Node() {}

	Node(string _name) {
		name = _name;
		value = 0;
		robot = 0;
	}

	void setValue(int _value) {
		value = _value;
	}

	void setRobot(int _r) {
		robot = _r;
	}

	int isRobot() {
		return robot;
	}

	void setNeighbor(Node* n) {
		for (int i = 0; i < neighbor.size(); i++)
			if(n->name.compare(neighbor[i]) == 0)
				return;
		neighbor.push_back(n->name);
	}

	void listAllNeighbors() {
		int i = 0;
		for (; i < neighbor.size() - 1; i++) {
			cout << neighbor[i] << ", ";
		}
		cout << neighbor[i];
	}
	~Node() {}
};

class Edge {
public:
	string name;
	string nameFull;
	string i, j;
	int cost;

	Edge() {}

	Edge(string _name, string _nameFull, int _cost) {
		name = _name;
		nameFull = _nameFull;
		cost = _cost;
	}

	void setNodes(Node* _i, Node* _j) {
		i = _i->name;
		j = _j->name;
	}
	~Edge() {}
};

class Graph {
public:
	map<string, Node> node;
	map<string, Edge> edge;
	int nSteinerNodes;

	Graph() {
		nSteinerNodes = 0;
	}

	Graph(int n_nodes) {
		for (int i = 0; i < n_nodes; i++) {
			ostringstream nodeName;
			nodeName << "y" << i;
			node.insert(pair<string, Node> (nodeName.str(), Node(nodeName.str())));
		}
		nSteinerNodes = node.size();
	}

	Node* getNode(int i) {
		ostringstream index;
		index << "y" << i;
		return &(node[index.str()]);
	}

	Edge* getEdge(int i){
		ostringstream index;
		index << "e" << i;
		return &(edge[index.str()]);
	}

	Edge* getEdge(string fullName) {
		map<string, Edge>::iterator it;
		for (it = edge.begin(); it != edge.end(); ++it) {
			if(fullName.compare((*it).second.nameFull) == 0)
				return &((*it).second);
		}
		return NULL;
	}

	void setNodeAsRobot(string nodeIndex) {
		node[nodeIndex].setRobot(1);
		nSteinerNodes--;
	}

	void resetNode(string nodeIndex) {
		node[nodeIndex].setRobot(0);
		nSteinerNodes++;
	}

	void addNode(int c) {
		ostringstream nodeName;
		nodeName << "y" << node.size();
		Node n(nodeName.str());
		n.setValue(c);
		node.insert(pair<string, Node>(nodeName.str(), n));
		nSteinerNodes++;
	}

	void addEdge(int i, int j, int c) {
		if (i == j) return;
		if (j > i) {
			int t = i; i = j; j = t;
		}
		if(i < 0 || j < 0 || i >= node.size()
				|| j >= node.size()) return;

		Node* nodei = getNode(i);
		Node* nodej = getNode(j);

		nodei->setNeighbor(nodej);
		nodej->setNeighbor(nodei);

		ostringstream edgeName, edgeNameFull;
		edgeName << "e" << edge.size();
		edgeNameFull << "x" << j << "." << i;

		if(getEdge(edgeNameFull.str()) == NULL) {
			Edge e(edgeName.str(), edgeNameFull.str(), c);
			e.setNodes(nodei, nodej);
			edge.insert(pair<string, Edge> (edgeName.str(), e));
		}
	}

	void listAllNodes() {
		map<string, Node>::iterator it;
		for (it = node.begin(); it != node.end(); ++it) {
			cout << (*it).first << " => ";
			(*it).second.listAllNeighbors();
			cout << endl;
		}
	}

	void listAllEdges() {
		map<string, Edge>::iterator it;
		for (it = edge.begin(); it != edge.end(); ++it) {
			cout << (*it).first << "\t " << (*it).second.nameFull << endl;
		}
	}

	int getTotalNodes() {
		return node.size();
	}

	int getTotalSteinerNodes() {
		return nSteinerNodes;
	}

	int getTotalRobotNodes() {
		return (node.size() - nSteinerNodes);
	}

	int getTotalEdges() {
		return edge.size();
	}

	// Assume nodes are unique
	vector<int> getEdgesFromSubset(vector<int> nodes, int subTour) {
		vector<int> _edgesSubTour;
		vector<int> _edgesCut;
		vector<string> _nodes;

		for (int i = 0; i < nodes.size(); i++) {
			Node* n = getNode(nodes[i]);
			_nodes.push_back(n->name);
		}

		map<string, Edge>::iterator it;
		int j = 0;
		for (it = edge.begin(); it != edge.end(); ++it, j++) {
			Edge e = (*it).second;
			int vertices = 0;
			for (int i = 0; i < _nodes.size(); i++) {
				if(e.i.compare(_nodes[i]) == 0) vertices++;
				if(e.j.compare(_nodes[i]) == 0) vertices++;
			}
			if(vertices == 1) _edgesCut.push_back(j);
			if(vertices == 2) _edgesSubTour.push_back(j);
		}

		if(subTour) return _edgesSubTour;
		else return _edgesCut;
	}
};

int gridSize = 3;
int gridPos(int r, int c) {
	return c+gridSize*r;
}

int main(int argc, char *argv[]) {
	GRBEnv* env = 0;
	GRBVar* x = 0;
	GRBVar* y = 0;

	int nRobots = 3;
	Graph g;

	vector<int> terminal_x, terminal_y;
	terminal_x.push_back(0);
	terminal_y.push_back(0);

	for(int i = 0; i < gridSize; i++) {
		for(int j = 0; j < gridSize; j++) {
			int c = 0;
			for (int k = 0; k < terminal_x.size(); k++) {
				c += abs(terminal_x[k]-i)+abs(terminal_y[k]-j);
			}
			if (c == 0)
				c = 100; // robot
			g.addNode(c);
		}
	}

	for(int i = 0; i < gridSize; i++) {
		for(int j = 0; j < gridSize; j++) {
			int c = gridPos(i, j);

			if (j != 0) {
				int l = gridPos(i, j-1);
				g.addEdge(c, l, 1);
			}

			if (j != gridSize-1) {
				int r = gridPos(i, j+1);
				g.addEdge(c, r, 1);
			}

			if (i != 0) {
				int u = gridPos(i-1, j);
				g.addEdge(c, u, 1);
			}

			if (i != gridSize-1) {
				int d = gridPos(i+1, j);
				g.addEdge(c, d, 1);
			}
		}
	}

	g.listAllNodes();
	g.listAllEdges();

	try {

		env = new GRBEnv();
		GRBModel model = GRBModel(*env);
		model.set(GRB_StringAttr_ModelName, "MST");

		int i = 0;
		map<string, Edge>::iterator itE;
		map<string, Node>::iterator itN;

		//
		x = model.addVars(g.getTotalEdges(), GRB_BINARY);
		model.update();

		for (i = 0, itE = g.edge.begin(); itE != g.edge.end(); ++itE, i++) {
			string edgeName = (*itE).first;
			int edgeCost = (*itE).second.cost;
			x[i].set(GRB_DoubleAttr_Obj, edgeCost);
			x[i].set(GRB_StringAttr_VarName, edgeName);
		}

		y = model.addVars(g.getTotalNodes(), GRB_BINARY);
		model.update();

		for (i = 0, itN = g.node.begin(); itN != g.node.end(); ++itN, i++) {
			string nodeName = (*itN).first;
			int nodeValue = (*itN).second.value;
			y[i].set(GRB_DoubleAttr_Obj, nodeValue);
			y[i].set(GRB_StringAttr_VarName, nodeName);
		}

		model.set(GRB_IntAttr_ModelSense, 1);
		model.update();

		// Constraints
		// Total number of robots
		ostringstream nodeSumString;
		GRBLinExpr nodeSumExpr = 0;
		for (i = 0; i < g.getTotalNodes(); i++) {
			nodeSumExpr += y[i];
		}
		model.addConstr(nodeSumExpr == nRobots, nodeSumString.str());

		// Sum of all edges = n-1
		ostringstream edgeSumString;
		GRBLinExpr edgeSumExpr = 0;
		edgeSumString << "edge sum";
		for (i = 0; i < g.getTotalEdges(); i++) {
			edgeSumExpr += x[i];
		}
		model.addConstr(edgeSumExpr == g.getTotalNodes() - 1, edgeSumString.str());

		// Subtour Elimination Constraints
		vector<int> nodeSubset;
		nodeSubset.push_back(0);
		nodeSubset.push_back(1);
		nodeSubset.push_back(2);
		nodeSubset.push_back(3);
		nodeSubset.push_back(4);

		int c_n = 2;
		while (c_n < nodeSubset.size()) {
			int c[c_n + 1];

			for (int j = 0; j < c_n; j++)
				c[j] = 0;

			c[c_n] = g.getTotalNodes();

			int j = 0;
			while (1) {

				int skip = 0;
				for (int i = 0; i < c_n; i++) {
					if (c[i] == c[i + 1])
						skip = 1;
				}

				if (!skip) {
					ostringstream subTourString;
					GRBLinExpr subTourExpr;

					vector<int> _nodes;
					for (int i = 0; i < c_n; i++) _nodes.push_back(c[i]);
					vector<int> _edges = g.getEdgesFromSubset(_nodes, 1);

					subTourString << "s";
					for(int i = 0; i < _edges.size(); i++) {
						subTourString << "." << _edges[i];
						subTourExpr += x[_edges[i]];
					}

					if(_edges.size() > 0)
						model.addConstr(subTourExpr <= c_n - 1, subTourString.str());
				}

				j = 0;
				while (c[j] == c[j + 1]) {
					c[j] = 0;
					j++;
				}
				c[j]++;
				if (j == c_n)
					break;
			}
			c_n++;
		}

		model.update();

		// Cutset Constraints
		c_n = 1;
		while (c_n < nodeSubset.size()) {
			int c[c_n + 1];

			for (int j = 0; j < c_n; j++)
				c[j] = 0;

			c[c_n] = g.getTotalNodes();

			int j = 0;
			while (1) {

				int skip = 0;
				for (int i = 0; i < c_n; i++) {
					if (c[i] == c[i + 1])
						skip = 1;
				}

				if (!skip) {
					ostringstream cutSetString;
					GRBLinExpr cutSetExpr;

					vector<int> _nodes;
					for (int i = 0; i < c_n; i++) _nodes.push_back(c[i]);
					vector<int> _edges = g.getEdgesFromSubset(_nodes, 0);

					cutSetString << "c";
					for(int i = 0; i < _edges.size(); i++) {
						cutSetString << "." << _edges[i];
						cutSetExpr += x[_edges[i]];
					}

					if(_edges.size() > 0)
						model.addConstr(cutSetExpr >= 1, cutSetString.str());
				}

				j = 0;
				while (c[j] == c[j + 1]) {
					c[j] = 0;
					j++;
				}
				c[j]++;
				if (j == c_n)
					break;
			}
			c_n++;
		}

		model.update();


		model.write("prob.lp");

		// Initialization
		//		// Guess at the starting point: close the plant with the highest
		//		// fixed costs; open all others
		//		// First, open all plants
		//		for (int p = 0; p < n; ++p) {
		//			y[p].set(GRB_DoubleAttr_Start, 1.0);
		//		}
		//
		//		// Now close the plant with the highest fixed cost
		//		cout << "Initial guess:" << endl;
		//		double maxFixed = -GRB_INFINITY;
		//		for (int p = 0; p < nPlants; ++p) {
		//			if (FixedCosts[p] > maxFixed) {
		//				maxFixed = FixedCosts[p];
		//			}
		//		}
		//		for (int p = 0; p < nPlants; ++p) {
		//			if (FixedCosts[p] == maxFixed) {
		//				open[3].set(GRB_DoubleAttr_Start, 0.0);
		//				cout << "Closing plant " << p << endl << endl;
		//				break;
		//			}
		//		}

		// Use barrier to solve root relaxation
		model.getEnv().set(GRB_IntParam_Method, GRB_METHOD_BARRIER);

		// Solve
		model.optimize();

		// Print solution
		cout << "\nTOTAL COSTS: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
		cout << "Solution:" << endl;
		cout << "Robot Positions:" << endl;
		for (int i = 0; i < g.getTotalNodes(); i++) {
			if (y[i].get(GRB_DoubleAttr_X) > 0)
				cout << y[i].get(GRB_StringAttr_VarName) << endl;
		}
		cout << "Edge Positions:" << endl;
		for (int i = 0; i < g.getTotalEdges(); i++) {
			if (x[i].get(GRB_DoubleAttr_X) > 0)
				cout << x[i].get(GRB_StringAttr_VarName) << endl;
		}

	} catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	} catch (...) {
		cout << "Exception during optimization" << endl;
	}

	delete env;
	return 0;
}

