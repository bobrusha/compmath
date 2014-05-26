#include <iostream>
#include <math.h>
#include <vector>

#define M_PI (3,14159265)
using namespace std;

const double A = 0.04;
const double B = 0.05;
const double xi = 1. / 18.;

double f1(vector <double>);
double f2(vector <double>);

vector <double> solve();

void print(const vector <double>);
int main(){
	print(solve(,,B*M_PI, A*M_PI, 2));

	system("pause");
	return 0;
}

double f1(vector <double> y){
	return A*y[1];
}

double f2(vector <double> y){
	return -B*y[0];
}

vector <double> solve( double l, double r, const double y0, const double y1, const int n){
	double h = (r - l) / n;
	double c2 = xi;
	double a21 = c2;
	double b2 = 1 / (2 * c2);
	double b1 = 1 - 1 / (2 * c2);

	vector <double> y;
	y.resize(2);
	vector <double> tmp;
	tmp.resize(2);

	y[0] = y0;
	y[1] = y1;

	for (int i = 0; i < n; ++i){
		double k11 = h * f1(y);
		double k12 = h * f2(y);

		tmp[0] = y[0] + a21 * k11;
		tmp[1] = y[1] + a21 * k12;

		double k21 = h * f1(tmp);
		double k22 = h * f2(tmp);

		y[0] += b1 * k11 + b2 * k21;
		y[1] += b1 * k12 + b2 * k22;
	}
	return y;
}



void print(const vector <double> x){
	for (int i = 0; i < x.size(); ++i){
		cout << x[i] << endl;
	}
	return;
}