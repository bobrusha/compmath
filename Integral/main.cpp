#include <iostream>
#include <iomanip>
#include "matrix.h"

#define pi (3.141592653589793)
#define M_2PI (2.*pi)
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
	
	cout << setprecision(12)<< calculateIntegralG(0.1, 2.3, AitkenG(0.1, 2.3));
	//cout << NewtonCotes(0.1, 2.3);
	//cout << calculateIntegralNC(0.1, 2.3, AitkenNC(0.1, 2.3)) << endl;
	system("pause");
	return 0;
}
double f ( const double x ){
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
	moments[2] = 0.357143 * pow(-0.1 + b, 0.80)*(0.0138889 + b / 9 + pow(b, 2)) - 0.357143 * pow(-0.1 + a, 0.80)*(0.0138889 + a / 9 + pow(a, 2));
	moments[3] = 0.0726174 * pow(-0.5 + 5. * b, 0.80)* (0.0014881 + 0.0119048*b + 0.107143 * pow(b, 2) + 1.*pow(b, 3)) -
		0.0726174 * pow(-0.5 + 5. * a, 0.80)* (0.0014881 + 0.0119048*a + 0.107143 * pow(a, 2) + 1.*pow(a, 3));
	moments[4] = 0.0574887*pow(-0.5 + 5.* b, 0.80) *(0.000156642 + 0.00125313 * b + 0.0112782*pow(b, 2) + 0.105263*pow(b, 3) + 1. * pow(b, 4)) -
		0.0574887*pow(-0.5 + 5.* a, 0.80) *(0.000156642 + 0.00125313 * a + 0.0112782*pow(a, 2) + 0.105263*pow(a, 3) + 1. * pow(a, 4));
	moments[5] = 0.0475769 * pow(-0.5 + 5. * b, 0.80) * (0.0000163168 + 0.000130535 * b + 0.00117481*pow(b, 2) + 0.0109649*pow(b, 3) + 0.104167*pow(b, 4) + 1.*pow(b, 5)) -
		0.0475769 * pow(-0.5 + 5. * a, 0.80) * (0.0000163168 + 0.000130535 * a + 0.00117481*pow(a, 2) + 0.0109649*pow(a, 3) + 0.104167*pow(a, 4) + 1.*pow(a, 5));
	return moments;
}
double Gauss (const double l, const double r){
	//моменты
	vector <double> M;
	M.resize(6);
	copy(M, calculateMomentsG(l, r));
	checkEps(M);
	//cout << "Moments: " << endl;
	//print(M);
	//решаем слау
	vector <double> mns;
	for (int i = 0; i < 3; ++i){
		mns.push_back(-M[i + 3]);
	}
	checkEps(mns);
	vector <vector <double> > mis;
	mis.resize(3);
	for (int i = 0; i < 3; ++i){
		for (int j = 0; j < 3; ++j){
			mis[i].push_back(M[i+j]);
		}
	}
	checkEps(mis);

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
	checkEps(a);
	//cout << "vector a: " << endl;
	//print(a);

	vector < double > x;
	x.resize(3);
	MethodNewton( (l+r)/2, a, x);
	//findRootsCardano(a, x);
	checkEps(x);
	cout << "x: " << endl;
	print(x);
	vector <double> A;
	A.resize(3);

	vector <vector <double>> X;
	X.resize(3);
	for (int i = 0; i < 3; ++i){
		for (int j = 0; j < 3; ++j){
			X[i].push_back(pow(x[j], i));
		}
	}
	checkEps(X);

	PLU(X, P, L, U);
	vector <double> ms;
	for (int i = 0; i < 3; ++i){
		ms.push_back(M[i]);
	}
	copy(A, solveLinerSystem(P,L, U, ms));
	checkEps(A);
	double res = 0.0;
	res = A[0] * f(x[0]) + A[1] * f(x[1]) + A[2] * f(x[2]);
	/*
	if (fabs(A[0]) > EPS){
		res += A[0] * f(x[0]);
	}
	if (fabs(A[1]) > EPS){
		res += A[1] * f(x[1]);
	}
	if ( fabs (A[2]) > EPS){
		res += A[2] * f(x[2]);
	}
	if (fabs (res) < EPS)
		res = 0.0;
	*/
	//cout << "res" << res << endl;
	return res;
}

double calculateIntegralWithStepG(const double l, const double r, const int h){
	double res = 0.0;
	for (int i = 0; i < h; ++i){
		if (i == 12 && h == 14){
			//cout << "QQ" << endl;
		}
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
		cout << "SH1: " << setprecision(12) << sh1 << endl;
		cout << " == " << endl;
		sh2 = calculateIntegralWithStepG(l, r, 2 * i);
		cout << "SH2: " << setprecision(12) << sh2 << endl;
		cout << " == " << endl;
		R = (sh2 - sh1) / (pow(L, n) - 1.0);
		cout << " ======= "<< endl;
		++i;
		if (i == 7){
			cout << "QQ" << endl;
		}
	} while (fabs(R) > EPS);
	cout << "Numbers of iterations: " << i - 1 << endl;
	return sh1;
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

	} while (fabs(w_f(x1, a)) > EPS);

	x[0] = x1;

	double x2, x3;
	double p = a[2] + x1;
	double q = a[1] + x1*(a[2] + x1);
	double D = p*p - 4*q;
	x[1] = ( -p - sqrt(D)) / 2.0;
	x[2] = ( -p + sqrt(D)) / 2.0;
	
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
void findRootsCardano(const vector <double> a_v, vector<double> & x){
	double a, b, c;
	a = a_v[2];
	b = a_v[1];
	c = a_v[0];

	double q, r, r2, q3;
	q = (a*a - 3.*b) / 9.; r = (a*(2.*a*a - 9.*b) + 27.*c) / 54.;
	r2 = r*r; q3 = q*q*q;
	if (r2 < q3) {
		double t = acos(r / sqrt(q3));
		a /= 3.; 
		q = -2.*sqrt(q);
		x[0] = q*cos(t / 3.) - a;
		x[1] = q*cos((t - M_2PI) / 3.) - a;
		x[2] = q*cos((t + M_2PI) / 3.) - a;
		return;
	}
	else {
		cout << "Error! Complex roots" << endl;
	}
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
	double res = 3 * pow(x, 2) + 2 * a[2] * x + a[1];
	return res;
}