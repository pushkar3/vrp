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
	int index;
	int value;
	int robot;
	vector<int> neighbor;
	vector<int> edge;

	Node() { }

	Node(int _index) {
		index = _index;
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

	void setNeighbor(Node* n, int edgeIndex) {
		for (int i = 0; i < neighbor.size(); i++)
			if(n->index == neighbor[i])	return;
		neighbor.push_back(n->index);
		edge.push_back(edgeIndex);
	}

	void listAllNeighbors() {
		int i = 0;
		for (; i < neighbor.size() - 1; i++) {
			cout << "y" << neighbor[i] << ", ";
		}
		cout << "y" << neighbor[i];

		cout << " (";
		for (i = 0; i < edge.size() - 1; i++) {
			cout << "e" << edge[i] << ", ";
		}
		cout << "e" << edge[i] << ")";
	}
	~Node() {}
};

class Edge {
public:
	string name;
	int index;
	int i, j;
	int cost;

	Edge() {}

	Edge(string _name, int _index, int _cost) {
		index = _index;
		name = _name;
		cost = _cost;
	}

	void setNodes(Node* _i, Node* _j) {
		i = _i->index;
		j = _j->index;
	}
	~Edge() {}
};

class Graph {
public:
	map<int, Node> node;
	map<int, Edge> edge;
	int nSteinerNodes;

	Graph() {
		nSteinerNodes = 0;
	}

	Graph(int n_nodes) {
		for (int i = 0; i < n_nodes; i++) {
			node.insert(pair<int, Node> (i, Node(i)) );
		}
		nSteinerNodes = node.size();
	}

	Node* getNode(int i) {
		return &(node[i]);
	}

	Edge* getEdge(int i){
		return &(edge[i]);
	}

	int getEdge(int i, int j) {
		int r = 0;
		map<int, Edge>::iterator it;
		for (it = edge.begin(); it != edge.end(); ++it) {
			int ei = (*it).second.i;
			int ej = (*it).second.j;
			if ( (i==ei && j==ej) || (i==ej && j==ei) ) {
				return r;
			}
			r++;
		}
		return -1;
	}

	void setNodeAsRobot(int nodeIndex) {
		node[nodeIndex].setRobot(1);
		nSteinerNodes--;
	}

	void resetNode(int nodeIndex) {
		node[nodeIndex].setRobot(0);
		nSteinerNodes++;
	}

	void addNode(int c) {
		Node n(node.size());
		n.setValue(c);
		node.insert(pair<int, Node>(node.size(), n));
		nSteinerNodes++;
	}

	void addEdge(int i, int j, int c) {
		if (i == j) return;
		if (i > j) {
			int t = i; i = j; j = t;
		}
		if(i < 0 || j < 0 || i >= node.size()
				|| j >= node.size()) return;

		Node* nodei = getNode(i);
		Node* nodej = getNode(j);

		ostringstream edgeName;
		edgeName << "x" << i << "." << j;

		nodei->setNeighbor(nodej, edge.size());
		nodej->setNeighbor(nodei, edge.size());

		if(getEdge(i, j) == -1) {
			Edge e(edgeName.str(), edge.size(), c);
			e.setNodes(nodei, nodej);
			edge.insert(pair<int, Edge> (edge.size(), e));
		}
	}

	void listAllNodes() {
		map<int, Node>::iterator it;
		for (it = node.begin(); it != node.end(); ++it) {
			cout << (*it).first << " [" << (*it).second.value << "] => ";
			(*it).second.listAllNeighbors();
			cout << endl;
		}
	}

	void listAllEdges() {
		map<int, Edge>::iterator it;
		for (it = edge.begin(); it != edge.end(); ++it) {
			cout << "e" << (*it).first << "\t " << (*it).second.name;
			cout << " => " << (*it).second.i << "," << (*it).second.j << endl;
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

		map<int, Edge>::iterator it;
		for (it = edge.begin(); it != edge.end(); ++it) {
			int eIndex = (*it).first;
			Edge e = (*it).second;
			int vertices = 0;
			for (int i = 0; i < nodes.size(); i++) {
				if(e.i == nodes[i]) vertices++;
				if(e.j == nodes[i]) vertices++;
			}
			if(vertices == 1) _edgesCut.push_back(eIndex);
			if(vertices == 2) _edgesSubTour.push_back(eIndex);
		}

		if(subTour) return _edgesSubTour;
		else return _edgesCut;
	}
};

int gridSize = 3;
int nRobots = 6;
int gridPos(int r, int c) {
	return c+gridSize*r;
}

int main(int argc, char *argv[]) {
	GRBEnv* env = 0;
	GRBVar* x = 0;  // edges
	GRBVar* y = 0;  // robots
	GRBVar* z = 0;  // holes

	Graph g;

	vector<int> terminal, terminal_x, terminal_y;
	terminal_x.push_back(1);
	terminal_y.push_back(1);

	for (int i = 0; i < terminal_x.size(); i++)
		terminal.push_back(gridPos(terminal_x[i], terminal_y[i]));

	for(int i = 0; i < gridSize; i++) {
		for(int j = 0; j < gridSize; j++) {
			int c = 0;
			for (int k = 0; k < terminal_x.size(); k++) {
				c += (abs(terminal_x[k]-i)+abs(terminal_y[k]-j));
			}
			g.addNode(c);
			cout << c << " ";
		}
		cout << endl;
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
		map<int, Edge>::iterator itE;
		map<int, Node>::iterator itN;

		cout << "Objective" << endl;
		// Objective function
		x = model.addVars(g.getTotalEdges(), GRB_BINARY);
		model.update();

		for (i = 0, itE = g.edge.begin(); itE != g.edge.end(); ++itE, i++) {
			int edgeIndex = (*itE).first;
			int edgeCost = (*itE).second.cost;
			edgeCost = 0;
			ostringstream edgeName;
			edgeName << "e" << edgeIndex;
			x[i].set(GRB_DoubleAttr_Obj, edgeCost);
			x[i].set(GRB_StringAttr_VarName, edgeName.str());
		}

		y = model.addVars(g.getTotalNodes(), GRB_BINARY);
		model.update();

		for (i = 0, itN = g.node.begin(); itN != g.node.end(); ++itN, i++) {
			int nodeIndex = (*itN).first;
			int nodeValue = (*itN).second.value;
			nodeValue = -10;
			ostringstream nodeName;
			nodeName << "y" << nodeIndex;
			y[i].set(GRB_DoubleAttr_Obj, nodeValue);
			y[i].set(GRB_StringAttr_VarName, nodeName.str());
		}

		z = model.addVars(g.getTotalNodes(), GRB_BINARY);
		model.update();

		for (i = 0, itN = g.node.begin(); itN != g.node.end(); ++itN, i++) {
			int nodeIndex = (*itN).first;
			int nodeValue = (*itN).second.value;
			ostringstream nodeName;
			nodeName << "z" << nodeIndex;
			z[i].set(GRB_DoubleAttr_Obj, nodeValue);
			z[i].set(GRB_StringAttr_VarName, nodeName.str());
		}

		model.set(GRB_IntAttr_ModelSense, 1);
		model.update();

		cout << "Constraints.";
		// Constraints
		// Total number of robots
		ostringstream nodeSumString;
		nodeSumString << "ysum." << nRobots;
		GRBLinExpr nodeSumExpr = 0;
		for (i = 0; i < g.getTotalNodes(); i++) {
			nodeSumExpr += y[i];
		}
		model.addConstr(nodeSumExpr == nRobots, nodeSumString.str());

		model.update();

		cout << ".";
		// Select each node as a robot or a hole
		for (i = 0; i < g.getTotalNodes(); i++) {
			ostringstream nodeConditionString;
			GRBLinExpr nodeConditionExpr = 0;
			nodeConditionString << "nc" << i;
			nodeConditionExpr = y[i] + z[i];
			model.addConstr(nodeConditionExpr <= 1, nodeConditionString.str());
		}

		model.update();

		cout << ".";
		// Sum of all edges = n-1
		ostringstream edgeSumString;
		GRBLinExpr edgeSumExpr = 0;
		GRBLinExpr vertexSumExpr = 0;
		edgeSumString << "edge sum";
		for (i = 0; i < g.getTotalEdges(); i++)
			edgeSumExpr += x[i];

		for (i = 0; i < g.getTotalNodes(); i++)
			vertexSumExpr += (y[i]+z[i]);

		model.addConstr(edgeSumExpr == vertexSumExpr - 1, edgeSumString.str());

		model.update();

		cout << "." << endl;

		// Subtour Elimination Constraints
#if 1
		for (int c_n = 2; c_n < g.getTotalNodes(); c_n++) {
			int c[c_n + 1];

			for (int j = 0; j < c_n+1; j++)
				c[j] = 0;

			c[c_n] = g.getTotalNodes();

			int j = 0;
			cout << "" << c_n << "/" <<  g.getTotalNodes() << endl;
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

					GRBLinExpr c_nExpr = 0;
					for(int i = 0; i < _nodes.size(); i++)
						c_nExpr += (y[_nodes[i]]+z[_nodes[i]]);

					if(_edges.size() > 0)
						model.addConstr(c_n*subTourExpr <= (c_n-1)*c_nExpr, subTourString.str());
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
		}

		model.update();
#else

		// Cutset Constraints
		int c_n = 1;
		while (c_n < g.getTotalNodes()) {
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
#endif

		cout << ".";
		// New cutset conditions
		for (int i = 0; i < g.getTotalNodes(); i++) {
			ostringstream ncutSetString_1, ncutSetString_2;
			GRBLinExpr ncutSetExpr_1, ncutSetExpr_2;

			vector<int> nodes;
			nodes.push_back(i);
			vector<int> edges = g.getEdgesFromSubset(nodes, 0);

			ncutSetString_1 << "nc1." << i;
			ncutSetString_2 << "nc2." << i;
			for(int j = 0; j < edges.size(); j++) {
				ncutSetExpr_1 += x[edges[j]];
				ncutSetExpr_2 += x[edges[j]];
			}

			if(edges.size() > 0) {
				model.addConstr(ncutSetExpr_1 >= y[i] + 2*z[i], ncutSetString_1.str());
				model.addConstr(ncutSetExpr_2 <= y[i] + 4*z[i], ncutSetString_2.str());
			}
		}

		model.update();

		cout << "." << endl;
		// Set Terminal Positions
		for (int k = 0; k < terminal.size(); k++) {
			ostringstream terminalString;
			GRBLinExpr terminalExpr = z[terminal[k]];
			terminalString << "terminal" << terminal[k];
			model.addConstr(terminalExpr == 1, terminalString.str());
		}

		model.update();
		model.write("prob.lp");

		// Initialization

		cout << "Solving" << endl;
		// Use barrier to solve root relaxation
		model.getEnv().set(GRB_IntParam_Method, GRB_METHOD_BARRIER);

		// Solve
		model.optimize();

		// Print solution
		cout << "\nTOTAL COSTS: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
		cout << "Solution" << endl;
		cout << endl;
		cout << "Robot/Hole Positions:" << endl;
		for (int i = 0; i < g.getTotalNodes(); i++) {
			if (y[i].get(GRB_DoubleAttr_X) > 0)
				cout << y[i].get(GRB_StringAttr_VarName) << " => " << y[i].get(GRB_DoubleAttr_X) << endl;
		}
		cout << endl;
		for (int i = 0; i < g.getTotalNodes(); i++) {
			if (z[i].get(GRB_DoubleAttr_X) > 0)
				cout << z[i].get(GRB_StringAttr_VarName) << " => " << z[i].get(GRB_DoubleAttr_X) << endl;
		}
		cout << endl;
		cout << "Edge Positions:" << endl;
		for (int i = 0; i < g.getTotalEdges(); i++) {
			if (x[i].get(GRB_DoubleAttr_X) > 0)
				cout << x[i].get(GRB_StringAttr_VarName) << " => " << x[i].get(GRB_DoubleAttr_X) << endl;
		}

		cout << endl;
		cout << endl;

		int c = 0;
		for (int i = 0; i < gridSize; i++) {
			for (int j = 0; j < gridSize; j++) {
				int n = gridPos(i, j);
				if (y[n].get(GRB_DoubleAttr_X) > 0)	cout << "R";
				else if (z[n].get(GRB_DoubleAttr_X) > 0) {
					int is_t = 0;
					for (int k = 0; k < terminal.size(); k++)
						if(n == terminal[k]) {
							is_t = 1;
						}
					if (is_t) cout << "T";
					else cout << "O";
				}
				else cout << " ";
				cout << " ";
				int c = g.getEdge(n, n+1);
				if (c >= 0 && x[c].get(GRB_DoubleAttr_X) > 0) cout << "--";
				else cout << "  ";
				cout << " ";
			}
			cout << endl;
			int n = gridPos(i, 0);
			for (int k = 0; k < gridSize; k++) {
				int c = g.getEdge(n+k, n+k+gridSize);
				if (c >= 0 && x[c].get(GRB_DoubleAttr_X) > 0) cout << "|";
				else cout << " ";
				cout << "    ";
			}
			cout << endl;
		}

		cout << endl;
		cout << endl;

		ofstream ofs;
		ostringstream ofsName;
		ofsName << "latex_" << gridSize << "_" << nRobots << ".txt";
		ofs.open(ofsName.str().c_str());
		ofs << "\\begin{tikzpicture}[scale=.5]\\footnotesize\n";
		ofs << " \\begin{scope}<+->;\n";
		ofs << "  \\draw[step=1cm,gray,very thin] (0, 0) grid (" << gridSize << ", " << gridSize << ");\n";
		ofs << " \\end{scope}\n";
		ofs << "\n";
		ofs << " \\begin{scope}[thin, black]\n";
		for (int e = 0; e < g.getTotalEdges(); e++) {
			if (x[e].get(GRB_DoubleAttr_X) > 0) {
				int i = g.getEdge(e)->i;
				int j = g.getEdge(e)->j;
				float xi = float(i/gridSize) + 0.5;
				float yi = float(i%gridSize) + 0.5;
				float xj = float(j/gridSize) + 0.5;
				float yj = float(j%gridSize) + 0.5;
				ofs << "  \\draw (" << xi << ", " << yi << ") node {}  -- (" << xj << ", " << yj << ") node {};\n";
			}
		}
		ofs << " \\end{scope}\n";
		ofs << "\n";
		ofs << " \\begin{scope}[thick,red]\n";
		for (int n = 0; n < g.getTotalNodes(); n++) {
			int i = n/gridSize;
			int j = n%gridSize;
			if (y[n].get(GRB_DoubleAttr_X) > 0) {
				ofs << "  \\filldraw[thin,blue,opacity=.2] ("<< i <<", "<< j <<") rectangle ("<< i+1 <<", "<< j+1 <<");\n";
			}
			else if (z[n].get(GRB_DoubleAttr_X) > 0) {
				int is_t = 0;
				for (int k = 0; k < terminal.size(); k++)
					if(n == terminal[k]) {
						is_t = 1;
					}
				if (is_t) ofs << "  \\filldraw[thin,red,opacity=.5] ("<< i <<", "<< j <<") rectangle ("<< i+1 <<", "<< j+1 <<");\n";
				else ofs << "  \\filldraw[thin,green,opacity=.3] ("<< i <<", "<< j <<") rectangle ("<< i+1 <<", "<< j+1 <<");\n";
			}
		}
		ofs << " \\end{scope}\n";
		ofs << "\\end{tikzpicture}\n";
		ofs.close();

	} catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	} catch (...) {
		cout << "Exception during optimization" << endl;
	}

	delete env;
	return 0;
}

