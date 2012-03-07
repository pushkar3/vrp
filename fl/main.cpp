#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>

#include <opencv/cv.h>
#include <opencv/highgui.h>

using namespace std;
using namespace cv;

struct line_t {
	Point2f p[2];
	double m;
	double b;
	double length;

	line_t(Point2f p1, Point2f p2) {
		p[0] = p1;
		p[1] = p2;
		double _x = p[0].x - p[1].x;
		double _y = p[0].y - p[1].y;
		if (_x != 0) m = fabs(_y/_x);
		else m = 0.0f;
		b = p[0].y - m*p[0].x;
		length = norm(p1-p2);
	}

	double distance_to(Point2f p) {
		double nr = fabs(p.y - m*p.x - b);
		double dr = sqrt(m*m + 1);
		return (nr/dr);
	}

	void print() {
		cout << "Line y = " << m << " x + " << b << endl;
	}
};

double r = 50.0f;
vector<Point2f> points;
vector<line_t> lines;


int main(int argc, char *argv[]) {
	double x_a[100], y_a[100], x_b[100], y_b[100];
	int n = 0;
	cout << "opening " << argv[1] << endl;
	fstream f(argv[1], fstream::in);

	Mat display(600, 600, CV_8UC1, Scalar(255));

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
		Point2f start(x_a[i], y_a[i]);
		Point2f end(x_b[i], y_b[i]);
		lines.push_back(line_t(start, end));
		line(display, start, end, Scalar(0));
	}

	for (int i = 0; i < lines.size(); i++) {
		lines[i].print();
	}

	f.close();

	double r = 50;
	cout << "r is " << r << endl;

	double h = r*cos(M_PI/3);

	vector<Point2f> pnts;
	for (int i = 0; i < lines.size(); i++) {
		while(r < lines[i].length) {
			//pnts.push_back(lines.p[0]+);
		}
	}

	imshow("display", display);
	waitKey(0);

	return 0;
}


