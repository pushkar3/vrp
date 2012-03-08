#include "gurobi_c++.h"
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
	vector<Node*> neighbor;

	Node(int _value) {
		value = _value;
	}

	void setNeighbor(Node* n) {
		neighbor.push_back(n);
	}

	void listAllNeighbors() {
		int i = 0;
		for (; i < neighbor.size()-1; i++) {
			cout << neighbor[i]->value << ", ";
		}
		cout << neighbor[i]->value;
	}
};

class Graph {
public:
	vector<Node> node;
	map<string, int> edge;
	vector<string> edge_counter;

	Graph(int n_nodes) {
		for (int i = 0; i < n_nodes; i++)
			node.push_back(Node(i));
	}

	void addEdge(int i, int j, int c) {
		if (i == j) return;
		if (j > i) {
			int t = i; i = j; j = t;
		}
		node[i].setNeighbor(&node[j]);
		node[j].setNeighbor(&node[i]);
		ostringstream edgeName;
		edgeName << "x" << j << "." << i;
		edge.insert(pair<string, int>(edgeName.str(), c));
	}

	void listAllNodes() {
		for (int i = 0; i < node.size(); i++) {
			cout << node[i].value << " => " << "[";
			node[i].listAllNeighbors();
			cout << "]" << endl;
		}
	}

	void listAllEdges() {
		edge_counter.clear();
		map<string, int>::iterator it;
		for (it = edge.begin(); it != edge.end(); ++it) {
			edge_counter.push_back((*it).first);
			cout << (*it).first << "\t " << (*it).second << endl;
		}
	}

	int findEdgeIndexByName(string name) {
		int i = 0;
		while (i < edge_counter.size()) {
			if (edge_counter[i].compare(name) == 0)
				return i;
		}
		return -1;
	}

	int findEdgeIndexByNode(int i, int j) {
		if (j > i) {
			int t = i; i = j; j = t;
		}
		ostringstream edgeName;
		edgeName << "x" << j << "." << i;
		return findEdgeIndexByName(edgeName.str());
	}

	string findEdgeByIndex(int i) {
		return edge_counter[i];
	}

	int findEdgeValue(string name) {
		return (edge[name]);
	}

	int findEdgeValue(int i) {
		return findEdgeValue(findEdgeByIndex(i));
	}

	int getTotalNodes() {
		return node.size();
	}

	int getTotalEdges() {
		return edge.size();
	}
};

int main(int argc, char *argv[]) {
	GRBEnv* env = 0;
	GRBVar* x = 0;

	int c[6];

	Graph g(6);
	g.addEdge(0, 2, 1);
	g.addEdge(0, 1, 2);
	g.addEdge(2, 3, 1);
	g.addEdge(3, 1, 1);
	g.addEdge(3, 4, 1);
	g.addEdge(4, 5, 1);

	g.listAllNodes();
	g.listAllEdges();

	int n = g.getTotalNodes();
	int m = g.getTotalEdges();

	try {

		env = new GRBEnv();
		GRBModel model = GRBModel(*env);
		model.set(GRB_StringAttr_ModelName, "MST");

		x = model.addVars(m, GRB_BINARY);
		model.update();

		for (int i = 0; i < m; i++) {
			int edgeCost = g.findEdgeValue(i);
			string edgeName = g.findEdgeByIndex(i);
			x[i].set(GRB_DoubleAttr_Obj, edgeCost);
			x[i].set(GRB_StringAttr_VarName, edgeName);
		}

		model.set(GRB_IntAttr_ModelSense, 1);
		model.update();


		// Constraints
		// Sum of all edges = n-1
		ostringstream edgeSumString;
		GRBLinExpr edgeSumExpr = 0;
		edgeSumString << "edge sum";
		for (int i = 0; i < m; i++) {
			edgeSumExpr += x[i];
		}
		model.addConstr(edgeSumExpr == n-1, edgeSumString.str());

		// Subtour Elimination Constraints
		int c_n = 2;
		while (c_n < n) {
			int c[c_n + 1];

			for (int j = 0; j < c_n; j++)
				c[j] = 0;

			c[c_n] = n;

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
					subTourString << "c";
					for (int i = 0; i < c_n; i++) {
						subTourString << "." << i;
						subTourExpr += x[i];
					}
					model.addConstr(subTourExpr <= c_n-1, subTourString.str());
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
		for (int i = 0; i < m; i++) {
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

