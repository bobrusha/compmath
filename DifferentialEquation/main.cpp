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
vector <double> getSolutionRungeBase(double, double, double, double, const int);
vector <double> getSolutionRungeOpponent(double, double, double, double, const int);
vector <double> getSolutionAutoStepBase(double, double, double, double);
vector <double> getSolutionAutoStepOpponent(double, double, double, double);

void print(const vector <double>);
void copy(vector <double>&, const vector <double>);
int main(){
	cout.precision(9);
	cout << "Base with Runge rule:" << endl;
	print(getSolutionRungeBase(0, M_PI, B*M_PI, A*M_PI, 2));
	cout << "Opponent with Runge rule:" << endl;
	print(getSolutionRungeOpponent(0, M_PI, B*M_PI, A*M_PI, 2));
	cout << "Base with auto step" << endl;
	print(getSolutionAutoStepBase(0, M_PI, B*M_PI, A*M_PI));
	cout << "Opponent with auto step" << endl;
	print(getSolutionAutoStepOpponent(0, M_PI, B*M_PI, A*M_PI));

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

	for (int i = 0; i < n; ++i){
		double k11 = h*f1(y);
		double k12 = h*f2(y);

		vector <double> tmp;
		tmp.push_back(y[0]);
		tmp.push_back(y[1]);

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
vector <double> getSolutionRungeBase(const double l, const double r, double y0, double y1, const int n){
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
		//cout << i << endl;
		//print(v2);
		err1 = (v2[0] - v1[0]) / 3.;
		err2 = (v2[1] - v1[1]) / 3.;
		//cout << err1 << "  " << err2 << endl << endl;
		++i;
	} while (fabs(err1) > eps || fabs(err2) > eps);
	cout << "Numbers of iterations = " << i << endl;
	return v2;
}
vector <double> getSolutionRungeOpponent(const double l, const double r, double y0, double y1, const int n){
	const double eps = 10e-9;
	double err1, err2;
	int i = 1;
	int steps = n;
	vector <double> v1, v2;
	v1.resize(2);
	v2.resize(2);
	copy(v2, solveWithOpponent(l, r, y0, y1, steps));
	do{
		steps *= 2.;
		copy(v1, v2);
		copy(v2, solveWithOpponent(l, r, y0, y1, steps));
		//cout << i << endl;
		//print(v2);
		err1 = (v2[0] - v1[0]) / 3.;
		err2 = (v2[1] - v1[1]) / 3.;
		//cout << err1 << "  " << err2 << endl<<endl;
		++i;
	} while (fabs(err1) > eps || fabs(err2) > eps);
	cout << "Numbers of iterations = " << i << endl;
	return v2;
}

vector <double> getSolutionAutoStepBase(const double l, const double r, const double y0, const double y1 ){
	const int k = 2;
	const double eps = 10e-9;
	int i = 0;

	vector <double> y;
	y.push_back(y0);
	y.push_back(y1);

	double delta;
	if (fabs(l) > fabs(r))
		delta = pow(1 / l, 3) + pow(sqrt(pow(f1(y), 2) + pow(f2(y), 2)), 3);
	else
		delta = pow(1 / r, 3) + pow(sqrt(pow(f1(y), 2) + pow(f2(y), 2)), 3);

	double h = pow(eps / delta, 1. / 3.);
	double xj = l;
	int steps = (int)((r - l) / h) + 1;
	h = (r - l) / steps;

	vector <double> v1, v2, vj;
	v1.resize(2);
	v2.resize(2);

	vj.push_back(y0);
	vj.push_back(y1);

	double err, err1, err2;

	while (xj < r){
		++i;
		copy(v1, solve(xj, xj + h, vj[0], vj[1], 1));
		copy(v2, solve(xj, xj + h / 2., vj[0], vj[1], 1));
		copy(v2, solve(xj + h / 2., xj + h, v2[0], v2[1], 1));

		err1 = (v2[0] - v1[0]) / (1 - pow(2, -2));
		err2 = (v2[1] - v1[1]) / (1 - pow(2, -2));
		err = sqrt(pow(err1, 2) + pow(err2, 2));

		if (fabs(err) >(eps * pow(2, k))){
			h /= 2.;
		}
		else if (err <= (eps * pow(2, k)) && err > eps){
			h /= 2.;
			copy(vj, v2);
		}
		else if (err <= eps && err >= (eps / pow(2, k + 1)))
			copy(vj, v1);
		else{
			h *= 2.;
			copy(vj, v1);
		}
		xj += h;
	}
	cout << "Numbers of iterations = " << i << endl;
	return vj;
}
vector <double> getSolutionAutoStepOpponent(const double l, const double r, const double y0, const double y1){
	const int k = 2;
	const double eps = 10e-9;
	int i = 0;

	vector <double> y;
	y.push_back(y0);
	y.push_back(y1);

	double delta;
	if (fabs(l) > fabs(r))
		delta = pow(1 / l, 3) + pow(sqrt(pow(f1(y), 2) + pow(f2(y), 2)), 3);
	else
		delta = pow(1 / r, 3) + pow(sqrt(pow(f1(y), 2) + pow(f2(y), 2)), 3);

	double h = pow(eps / delta, 1. / 3.);
	double xj = l;
	int steps = (int)((r - l) / h) + 1;
	h = (r - l) / steps;

	vector <double> v1, v2, vj;
	v1.resize(2);
	v2.resize(2);

	vj.push_back(y0);
	vj.push_back(y1);

	double err, err1, err2;

	while (xj < r){
		++i;
		copy(v1, solveWithOpponent(xj, xj + h, vj[0], vj[1], 1));
		copy(v2, solveWithOpponent(xj, xj + h/2., vj[0], vj[1], 1));
		copy(v2, solveWithOpponent(xj + h / 2., xj + h, v2[0], v2[1], 1));

		err1 = (v2[0] - v1[0]) / (1 - pow(2,-2));
		err2 = (v2[1] - v1[1]) / (1 - pow(2,-2));
		err = sqrt(pow(err1, 2) + pow(err2, 2));

		if (fabs(err) > (eps * pow(2, k))){
			h /= 2.;
		}
		else if (err <= (eps * pow(2, k)) && err > eps){
			h /= 2.;
			copy(vj, v2);
		}
		else if (err <= eps && err >= (eps / pow(2, k + 1)))
			copy(vj, v1);
		else{
			h *= 2.;
			copy(vj, v1);
		}
		xj += h;
	}
	cout << "Numbers of iterations = " << i << endl;
	return vj;
}
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