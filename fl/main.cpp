#include "gurobi_c++.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>
using namespace std;

struct point_t {
	double x, y;
	double c;
	point_t() {
		x = 0.0f;
		y = 0.0f;
		c = 0.0f;
	}

	point_t(double _x, double _y) {
		x = _x;
		y = _y;
		c = 0.0f;
	}

	double distance_to(point_t p) {
		double x = p.x - x;
		double y = p.y - y;
		return sqrt(x*x + y*y);
	}

	void print() {
		cout << x << " " << y << " cost: " << c << endl;
	}
};

struct line_t {
	point_t p[2];
	double m;
	double b;

	line_t(point_t p1, point_t p2) {
		p[0] = p1;
		p[1] = p2;
		double _x = p[0].x - p[1].x;
		double _y = p[0].y - p[1].y;
		if (_x != 0) m = fabs(_y/_x);
		else m = 0.0f;
		b = p[0].y - m*p[0].x;
	}

	double distance_to(point_t p) {
		double nr = fabs(p.y - m*p.x - b);
		double dr = sqrt(m*m + 1);
		return (nr/dr);
	}

	void print() {
		cout << "Line y = " << m << " x + " << b << endl;
	}
};

double r = 50.0f;
vector<point_t> points;
vector<line_t> lines;


int main(int argc, char *argv[]) {
	double x_a[100], y_a[100], x_b[100], y_b[100];
	int n = 0;
	cout << "opening " << argv[1] << endl;
	fstream f(argv[1], fstream::in);

	int i = 0;
	if(f.is_open()) {
		while(f.good()) {
			f >> x_a[i] >> y_a[i] >> x_b[i] >> y_b[i];
			i++;
		}
	}
	else
		cout << "cannot open" << endl;

	n = i-1;

	for (int i = 0; i < n; i++) {
		lines.push_back(line_t(point_t(x_a[i], y_a[i]), point_t(x_b[i], y_b[i])));
	}

	for (int i = 0; i < lines.size(); i++) {
		lines[i].print();
	}

	f.close();
	return 0;
}

