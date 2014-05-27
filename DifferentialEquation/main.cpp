#include <iostream>
#include <math.h>
#include <vector>

#define M_PI (3.14159265)
using namespace std;

const double A = 0.04;
const double B = 0.05;
const double xi = 1. / 18.;

double f1(vector <double>);
double f2(vector <double>);

vector <double> solve(double, double, const double, const double, const int);
vector <double> solveWithOpponent(double, double, double, double, const int);
vector <double> getSolutionRunge(double, double, double, double, const int);
void print(const vector <double>);
void copy(vector <double>&, const vector <double>);
int main(){
	print(getSolutionRunge(0, M_PI, B*M_PI, A*M_PI, 2));
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

	double k11, k12, k21, k22;
	for (int i = 0; i < n; ++i){
		k11 = h * f1(y);
		k12 = h * f2(y);

		tmp[0] = y[0] + a21 * k11;
		tmp[1] = y[1] + a21 * k12;

		k21 = h * f1(tmp);
		k22 = h * f2(tmp);

		y[0] += b1 * k11 + b2 * k21;
		y[1] += b1 * k12 + b2 * k22;
	}
	return y;
}

vector <double> solveWithOpponent(double l, double r, double y0, double y1, const int n){
	double h = (r - l) / n;
	double c2 = xi;
	double a21 = c2;
	double b2 = 1 / (2 * c2);
	double b1 = 1 - 1 / (2 * c2);

	vector <double> y;
	y.push_back(y0);
	y.push_back(y1);

	double k11 = h*f1(y);
	double k12 = h*f2(y);

	vector <double> tmp;
	tmp.push_back(y[0]);
	tmp.push_back(y[1]);
	for (int i = 0; i < n; ++i){
		tmp[0] = y[0] + .5 * k11;
		tmp[1] = y[1] + .5 * k12;

		double k21 = h*f1(tmp);
		double k22 = h*f2(tmp);

		tmp[0] += 2 * k21 - k11;
		tmp[1] += 2 * k22 - k12;

		double k31 = h*f1(tmp);
		double k32 = h*f2(tmp);

		y[0] += (1. / 6.)*(k11 + 4 * k21 + k31);
		y[1] += (1. / 6.)*(k12 + 4 * k22 + k32);
	}
	return y;
}

vector <double> getSolutionRunge(const double l, const double r, double y0, double y1, const int n){
	const double eps = 10e-9;
	double err1, err2;
	int i = 1;
	int steps = n;
	vector <double> v1, v2;
	v1.resize(2);
	v2.resize(2);
	copy(v2, solve(l, r, y0, y1, steps));
	do{
		steps *= 2.;
		copy(v1, v2);
		copy(v2, solve(l, r, y0, y1, steps));
		cout << i << endl;
		print(v2);
		err1 = (v2[0] - v1[0])/3.;
		err2 = (v2[1] - v1[1])/3.;
		cout << err1 << "  " << err2 << endl<<endl;
		++i;
	} while (fabs(err1) > eps || fabs(err2) > eps);
	return v2;
}
/*
vector <double> getSolutionAutoStep(const double l, const double r, const double y0, const double y1 ){
	const double eps = 10e-5;
	vector <double> y;
	y.push_back(y0);
	y.push_back(y1);

	double delta;
	if (l > r)
		delta = pow(1 / l, 3) + pow(sqrt(pow(f1(y), 2) + pow(f2(y), 2)), 3);
	else
		delta = pow(1 / r, 3) + pow(sqrt(pow(f1(y), 2) + pow(f2(y), 2)), 3);
	double h = pow(eps / delta,1./3.);

}
*/


void print(const vector <double> x){
	for (int i = 0; i < x.size(); ++i){
		cout << x[i] << endl;
	}
	return;
}
void copy(vector <double>& y1, const vector <double> y2){
	for (int i = 0; i < y1.size(); ++i){
		y1[i] = y2[i];
	}
	return;
}