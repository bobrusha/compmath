#include <iostream>
#include "matrix.h"

#define pi (3.141592653589793)
//result = 3.57886

double w_f (double, vector < double > );
double w_d(double, vector < double >);
void MethodNewton (const double , const vector <double>, vector <double>& );

double Gauss(const double, const double);
double NewtonCotes(const double , const double );
double calculateIntegralG(const double, const double, const double);
double calculateIntegralNC(const double, const double, const double);

double AitkenG(const double, const double);
double AitkenNC(const double, const double);

void findRootsCardano(const vector<double>, vector <double> &);

int main(){
	//cout<<Gauss(0.1, 2.3);
	cout << calculateIntegralG(0.1, 2.3, 3);
	//cout << NewtonCotes(0.1, 2.3);
	//cout << calculateIntegralNC(0.1, 2.3, AitkenNC(0.1, 2.3)) << endl;
	system("pause");
	return 0;
}
double f ( double x ){
	double res = 2.5*cos(2 * x)*exp(2 * x / 3) + 4 * sin(3.5*x)*exp(-3 * x) + 3 * x;
	if (fabs(res) < EPS){
		res = 0.0;
	}
	return res;
}

double p (double x, double a, double b){
	double alfa = 0.2, beta = 0;
	double res = pow((x - a), (-alfa))*pow((b - x), (-beta));
	if (fabs(res) < EPS){
		res = 0.0;
	}
	return res;
}
vector <double> calculateMomentsG(double a, double b){
	vector <double> moments;
	moments.resize(6);
	moments[0] = 1.250 *(pow(b - 0.10, 0.80) - pow(a - 0.10, 0.80));
	moments[1] = 0.0110062 * ((1 + 8 * b) *pow(-1 + 10 * b, 0.8) -  (1 + 8 * a) * pow(-1 + 10 * a, 0.8));
	moments[2] = pow(-0.1 + b, 0.8) * ((0.00496032 + (0.0396825 + 0.357143 * b) * b) -
				 (0.00496032 + (0.0396825 + 0.357143 * a) * a));
	moments[3] = 0.0000620651 * pow(-1 + 10 * b, 0.8) *((1 + 8 * b + 72 * pow(b, 2) + 672 * pow(b, 3)) - 
		(1 + 8 * a + 72 * pow(a, 2) + 672 * pow(a, 3)));
	moments[4] = 0.00012413 * pow(-1 + 10 * b, 0.8) * ((0.0416667 + b/3 + 3*pow(b, 2) + 28*pow( b , 3) + 266 *pow(b, 4)) -
		 (0.0416667 + a / 3 + 3 * pow(a, 2) + 28 * pow(a, 3) + 266 * pow(a, 4)));
	moments[5] = 0.0000299624 * pow( - 1.0 + 10.0 * b, 0.80) * ((0.014881 + 0.119048*b + 1.07143 * pow(b, 2) +
		10 * pow(b, 3) + 95 * pow(b, 4) + 912 * pow(b, 5)) - (0.014881 + 0.119048*a + 1.07143 * pow(a, 2) +
		10 * pow(a, 3) + 95 * pow(a, 4) + 912 * pow(a, 5)));

	return moments;
}
double Gauss (const double l, const double r){
	//моменты
	vector <double> M;
	M.resize(6);
	copy(M, calculateMomentsG(l, r));

	//решаем слау
	vector <double> mns;
	for (int i = 0; i < 3; ++i){
		mns.push_back(-M[i + 3]);
	}

	vector <vector <double> > mis;
	mis.resize(3);
	for (int i = 0; i < 3; ++i){
		for (int j = 0; j < 3; ++j){
			mis[i].push_back(M[i+j]);
		}
	}
	vector <vector <double>> P, L, U;
	P.resize(3);
	L.resize(3);
	U.resize(3);

	for (int i = 0; i < 3; ++i){
		P[i].resize(3);
		L[i].resize(3);
		U[i].resize(3);
	}

	PLU(mis, P, L, U);

	vector <double> a;
	a.resize(3);
	copy(a, solveLinerSystem(P, L, U, mns));

	vector < double > x;
	MethodNewton( (l+r)/2, a, x);
	//findRootsCardano(a, x);

	vector <double> A;
	A.resize(3);

	vector <vector <double>> X;
	X.resize(3);
	for (int i = 0; i < 3; ++i){
		for (int j = 0; j < 3; ++j){
			X[i].push_back(pow(x[j], i));
		}
	}

	PLU(X, P, L, U);
	vector <double> ms;
	for (int i = 0; i < 3; ++i){
		ms.push_back(M[i]);
	}
	copy(A, solveLinerSystem(P,L, U, ms));

	double res = A[0] * f(x[0]) + A[1] * f(x[1]) + A[2] * f(x[2]);
	if (fabs(res) < EPS)
		res = 0.0;
	return res;
}

double calculateIntegralWithStepG(const double l, const double r, const int h){
	double res = 0;
	for (int i = 0; i < h; ++i){
		res += Gauss(l + i * ((r - l) / h), l + (i + 1)*((r - l) / h));
	}
	return res;
}
double calculateIntegralG(const double l, const double r, const double n)
{
	const double L = 2.0;
	double sh1 = 0.0;
	double sh2 = 0.0;
	double R = 23.0;
	unsigned int i = 1;
	do{
		sh1 = calculateIntegralWithStepG(l, r, i);
		sh2 = calculateIntegralWithStepG(l, r, 2 * i);
		R = (sh2 - sh1) / (pow(L, n) - 1.0);
		cout << sh2 << " - " << sh1 << endl;
		++i;
	} while (fabs(R) > EPS);
	cout << "Numbers of iterations: " << i - 1 << endl;
	return sh2;
}
vector<double> calculateMomentsForCotes(const double z0, const double z1){
	vector <double> res;
	res.resize(3);

	res[0] = 10 * (pow(z1 - 0.1, 0.8) - pow(z0 - 0.1, 0.8)) / (8);
	res[1] = 10 * ((pow(z1 - 0.1, 1.8) - pow(z0 - 0.1, 1.8)) / (18)) 
		+ 0.1*res[0];
	res[2] = 10*((pow(z1 - 0.1, 2.8) - pow(z0 - 0.1, 2.8)) / (28))
		+ 0.2 * res[1] - 0.01 * res[0];
	
	return res;
}
double NewtonCotes(const double a, const double b){
	unsigned int m = 1;
	unsigned int k = 2*m; //количество отрезков
	double H = (b - a) / k;

	vector <double> moments;
	moments.resize(k+1);

	vector <double> z;
	for (int i = 0; i < k + 1; ++i){
		z.push_back(a+i*H);
	}

	copy(moments, calculateMomentsForCotes( z[0] , z[2]));
	checkEps(moments);
	
	//cout << "Matrix moments:" << endl;
	//print(moments);

	vector <double> A;
	A.push_back((moments[2] - moments[1]*(z[1] + z[2]) + moments[0]*z[2]*z[1])/
		((z[1] - z[0])*(z[2] - z[0])));
	A.push_back((- moments[2] + moments[1] * (z[0] + z[2]) - moments[0] * z[0] * z[2]) /
		((z[1] - z[0])*(z[2] - z[1])));
	A.push_back((moments[2] - moments[1] * (z[0] + z[1]) + moments[0] * z[0] * z[1]) /
		((z[2] - z[1])*(z[2] - z[0])));
	checkEps(A);
	//cout << "A: " << endl;
	//print(A);
	
	//cout << f(z[0]) << "   " << f(z[1]) << "  " << f(z[2]) << endl;
	double res;

	res = A[0] * f(z[0]) + pow(A[2], m)*f(z[k]);
	for (int i = 1; i <= m; ++i){
		res+= pow(A[1],i) *f(z[2*i-1]);
	}
	return res;
}

//находит корни на отрезке от l до r
void MethodNewton(const double x0, const vector <double> a, vector <double>& x){
	double prev = x0;
	double x1 = x0;
	do {
		x1 -= w_f(prev, a) / w_d(prev, a);
		prev = x1;

	} while (fabs(w_f(x1, a) - w_f(prev, a)) > EPS);

	x.push_back(x1);

	double x2, x3;
	double p = a[2] + x1;
	double q = a[1] + x1*(a[2] + x0);
	double D = pow(p,2) - 4*q;
	x2 = ( -p - sqrt(D)) / 2.0;
	x3 = ( -p + sqrt(D)) / 2.0;
	x.push_back(x2);
	x.push_back(x3);
	return;
}
double calculateIntegralWithStepNC(const double l, const double r, const double h){
	double res = 0;
	for (int i = 0; i < h; ++i)
		res += NewtonCotes(l + i * ((r - l) / h), l + (i + 1)*((r - l) / h));
	return res;
}

double calculateIntegralNC(const double l, const double r, const double n){
	const double L = 2.0;
	double sh1 = 0.0;
	double sh2 = 0.0;
	double R = 0.0;
	int i = 1;

	do{
		sh1 = calculateIntegralWithStepNC(l, r, i);
		sh2 = calculateIntegralWithStepNC(l, r, 2*i);
		R = (sh2 - sh1) / (pow(L, n) - 1.0);
		//cout << sh2 << " - " << sh1 << endl;
		++i;
	} while (fabs(R) > EPS);
	cout << "Numbers oi iterations: " << i << endl;
	return sh2;
}
void findRootsCardano(const vector <double> a, vector<double> & x){
	double Q = (pow(a[0], 2) - 3 * a[1]) / 9;
	double R = (2 * pow(a[0], 3) - 9 * a[0] * a[1] + 27 * a[2]) / 54;
	double t = acos(R / sqrt(pow(Q, 3))) / 3;
	x.push_back(-2 * sqrt(Q)*cos(t) - a[0] / 3),
	x.push_back(-2 * sqrt(Q)*cos(t + (2 * pi / 3)) - a[0] / 3);
	x.push_back(-2 * sqrt(Q)*cos(t - (2 * pi / 3)) - a[0] / 3);
	return;
}

double AitkenNC(const double l, const double r){

	const int h = 8;
	double sh1 = calculateIntegralWithStepNC(l,r, h);
	double sh2 = calculateIntegralWithStepNC(l, r, h / 2);
	double sh3 = calculateIntegralWithStepNC(l, r, h / 4);

	double cm = (sh1 - sh2) / (sh2 - sh3);
	double m = log(cm) / log(l);

	return m;
}
double AitkenG(const double l, const double r){
	const int h = 8;
	double sh1 = calculateIntegralWithStepG(l, r, h);
	double sh2 = calculateIntegralWithStepG(l, r, h / 2);
	double sh3 = calculateIntegralWithStepG(l, r, h / 4);

	double cm = (sh1 - sh2) / (sh2 - sh3);
	double m = log(cm) / log(l);

	return m;
}
double w_f ( const double x, const vector <double> a){
	double res = pow(x, 3) + a[2] * pow(x, 2) + a[1] * x + a[0];

	return res;
}

double w_d ( double x, vector <double> a){
	double res = 3 * pow(x, 2) + 2 * a[0] * x + a[1];
	if (fabs(res) < EPS)
		res = 0.0;
	return res;
}