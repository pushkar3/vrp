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
	GRBVar* y = 0;
	GRBVar* z = 0;
	GRBVar** xc = 0;
	GRBVar* yc = 0;
	GRBVar* zc = 0;

	try {
		const int n = 6;
		double c[] = { 1, 0, 1, 1, 2, 2 };
		int x_adj[][n] = {
			{	0, 1, 0, 0, 0, 0 },
			{	1, 0, 1, 1, 0, 0 },
			{	0, 1, 0, 0, 0, 0 },
			{	0, 1, 0, 0, 1, 1 },
			{	0, 0, 0, 1, 0, 0 },
			{	0, 0, 0, 1, 0, 0 }
		};

		env = new GRBEnv();
		GRBModel model = GRBModel(*env);
		model.set(GRB_StringAttr_ModelName, "facility");

		y = model.addVars(n, GRB_BINARY);
		z = model.addVars(n, GRB_BINARY);
		model.update();
		for (int p = 0; p < n; ++p) {
			ostringstream yname, zname;
			yname << "y" << p;
			zname << "z" << p;
			y[p].set(GRB_DoubleAttr_Obj, c[p]);
			y[p].set(GRB_StringAttr_VarName, yname.str());
			z[p].set(GRB_DoubleAttr_Obj, 0);
			z[p].set(GRB_StringAttr_VarName, zname.str());
		}

		x = new GRBVar*[n];
		for (int w = 0; w < n; ++w) {
			x[w] = model.addVars(n, GRB_BINARY);
			model.update();
			for (int p = 0; p < n; ++p) {
				ostringstream xname;
				xname << "x" << p << "." << w;
				x[w][p].set(GRB_DoubleAttr_Obj, 0);
				x[w][p].set(GRB_StringAttr_VarName, xname.str());
			}
		}

		model.set(GRB_IntAttr_ModelSense, 1);

		model.update();

		// Constraints
		ostringstream cname;
		cname << "TotalRobots" << n;
		GRBLinExpr y_total = 0;
		for (int p = 0; p < n; ++p) {
			y_total += y[p];
		}
		model.addConstr(y_total >= 4, cname.str());

		for (int i = 0; i < n; i++) {
			for (int j = i+1; j < n; j++) {
				ostringstream x_dir_s;
				GRBLinExpr x_dir = 0;
				x_dir = x[i][j] - x[j][i];
				x_dir_s << "xdir" << i << j;
				//model.addConstr(x_dir == 0, x_dir_s.str());
			}
		}

		for (int i = 0; i < n; i++) {
			ostringstream yz;
			yz << "yz" << i << i;
			model.addConstr(y[i]+z[i] == 1, yz.str());
		}

		for (int i = 0; i < n; i++) {
			for (int j = i+1; j < n; j++) {
				if(x_adj[i][j] != 0) {
					GRBLinExpr yz1 = 0;
					GRBLinExpr yz2 = 0;

					yz1 = y[i] + z[j] - x[j][i];
					yz2 = -z[i] + x[i][j];

					ostringstream yz1_s, yz2_s;
					yz1_s << "yz1 " << i << i << j;
					yz2_s << "yz2 " << i << i << j;
					model.addConstr(yz1 >= 0, yz1_s.str());
					model.addConstr(yz2 == 0, yz2_s.str());

				}
			}
		}

		// Initialization
		model.addConstr(z[1] == 1, "z1");

		model.update();
		model.write("prob.lp");


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
		cout << "SOLUTION:" << endl;
		for (int p = 0; p < n; ++p) {
			cout << p << "\t" << y[p].get(GRB_DoubleAttr_X) << "  " << z[p].get(GRB_DoubleAttr_X) << endl;
		}

		cout << endl;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				cout << x[i][j].get(GRB_DoubleAttr_X) << "\t";
			}
			cout << endl;
		}


	} catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	} catch (...) {
		cout << "Exception during optimization" << endl;
	}

	delete[] y;
	delete[] z;
	delete env;
	return 0;
}

