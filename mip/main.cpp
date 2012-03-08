#include "gurobi_c++.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
using namespace std;

int main(int argc, char *argv[]) {
	GRBEnv* env = 0;
	GRBVar** x = 0;

	int c[6];

	try {
		const int n = 6;
		double c[][n] = {
			{ 0, 0, 0, 0, 0, 0 },
			{ 0, 0, 0, 0, 0, 0 },
			{ 0, 0, 0, 0, 0, 0 },
			{ 0, 0, 0, 0, 0, 0 },
			{ 0, 0, 0, 0, 0, 0 },
			{ 0, 0, 0, 0, 0, 0 }
		};
		int x_adj[][n] = {
			{ 0, 0, 0, 0, 0, 0 },
			{ 0, 0, 0, 0, 0, 0 },
			{ 0, 0, 0, 0, 0, 0 },
			{ 0, 0, 0, 0, 0, 0 },
			{ 0, 0, 0, 0, 0, 0 },
			{ 0, 0, 0, 0, 0, 0 }
		};

		x_adj[1-1][3-1] = 1;
		x_adj[3-1][4-1] = 1;
		x_adj[2-1][4-1] = 1;
		x_adj[4-1][5-1] = 1;
		x_adj[5-1][6-1] = 1;

		for (int i = 0; i < n; i++) {
			for (int j = i; j < n; j++) {
				if (x_adj[i][j] == 1 || x_adj[j][i] == 1) {
					x_adj[i][j] = x_adj[j][i] = 1;
					c[i][j] = c[j][i] = 1;
				}
				else {
					c[i][j] = c[j][i] = 100;
				}
			}
		}

		env = new GRBEnv();
		GRBModel model = GRBModel(*env);
		model.set(GRB_StringAttr_ModelName, "MST");

		x = new GRBVar*[n];
		for (int i = 0; i < n; i++) {
			x[i] = model.addVars(n, GRB_BINARY);
			model.update();
			for (int j = i+1; j < n; j++) {
				ostringstream xname;
				xname << "x" << i+1 << "." << j+1;
				x[i][j].set(GRB_DoubleAttr_Obj, c[i][j]);
				x[i][j].set(GRB_StringAttr_VarName, xname.str());
			}
		}

		model.set(GRB_IntAttr_ModelSense, 1);

		model.update();

		// Constraints
		// Sum of all edges = n-1
		ostringstream x_allString;
		GRBLinExpr x_allExpr = 0;
		x_allString << "x_all";
		for (int i = 0; i < n; i++) {
			for (int j = i+1; j < n; j++) {
				x_allExpr += x[i][j];
			}
		}
		model.addConstr(x_allExpr == n-1, x_allString.str());

		// Subtour Elimination Constraints
		ostringstream subTourString;
		GRBLinExpr subTourExpr;
		subTourString << "subtour1";

		subTourExpr += x[0][2];
		subTourExpr += x[2][3];
		//model.addConstr(subTourExpr <= 1, subTourString.str());
//		for (int j = 0; j < n; j++)
//			c[j] = 0;
//
//		c[0] = 1;
//
//		while (1) {
//			ostringstream subTourString;
//			GRBLinExpr subTourExpr;
//			for (int i = 0; i < n; i++) {
//				printf("%d ", c[i]);
//			}
//			printf("\n");
//
//			int j = 0;
//			while (c[j] == n) {
//				c[j] = 0;
//				j++;
//			}
//
//			c[j]++;
//
//			if (j == n-1)
//				break;
//
//		}


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
		cout << endl;
		for (int i = 0; i < n; i++) {
			for (int j = i; j < n; j++) {
				if (x[i][j].get(GRB_DoubleAttr_X) > 0)
					cout << i+1 << ", " << j+1  << " is " << x[i][j].get(GRB_DoubleAttr_X) << endl;
			}
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

