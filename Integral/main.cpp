#include <iostream>
#include "matrix.h"
//result = 3.57886

double w_f (double, vector < double > );
double w_d(double, vector < double >);
double MethodNewton (const double , const double , const vector <double>);

double gauss(const double, const double);
double NewtonCotes(const double , const double );
double calculateIntegralG(const double, const double, const int);
double calculateIntegralNC(const double, const double, const int);

int main(){
	//cout<<gauss();
	//cout << NewtonCotes(0.1, 2.3);
	cout << calculateIntegralNC(0.1, 2.3, 5) << endl;
	system("pause");
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
vector <double> calculateMoments(double a, double b){
	vector <double> moments;
	moments.resize(6);
	moments[0] = 1.25 * pow(b - 0.1, 0.8) - 1.25 * pow(a - 0.1, 0.8);
	moments[1] = 0.0110062 * (1 + 8 * b) *pow(-1 + 10 * b, 0.8) - 0.0110062 * (1 + 8 * a) * pow(-1 + 10 * a, 0.8);
	moments[2] = pow(-0.1 + b, 0.8) * (0.00496032 + (0.0396825 + 0.357143 * b) * b) -
				pow(-0.1 + a, 0.8) * (0.00496032 + (0.0396825 + 0.357143 * a) * a);
	moments[3] = 0.0000620651 * pow(-1 + 10 * b, 0.8) *(1 + 8 * b + 72 * pow(b, 2) + 672 * pow(b, 3));
	moments[4] = 0.00012413 * pow(-1 + 10 * b, 0.8) * (0.0416667 + b/3 + 3*pow(b, 2) + 28*pow( b , 3) + 266 *pow(b, 4)) -
		0.00012413 * pow(-1 + 10 * a, 0.8) * (0.0416667 + a / 3 + 3 * pow(a, 2) + 28 * pow(a, 3) + 266 * pow(a, 4));
	moments[5] = 0.0000299624 * (-1 + 10 * b, 0.8) * (0.014881 + 0.119048*b + 1.07143 * pow(b, 2) +
		10 * pow(b, 3) + 95 * pow(b, 4) + 912 * pow(b, 5)) - 0.0000299624 * (-1 + 10 * a, 0.8) * (0.014881 + 0.119048*a + 1.07143 * pow(a, 2) +
		10 * pow(a, 3) + 95 * pow(a, 4) + 912 * pow(a, 5));

	return moments;
}
double gauss (const double l, const double r){
	//1. Считаем моменты весовой функции

	vector <double> moments;
	moments.resize(6);
	copy(moments, calculateMoments(l, r));
	checkEps(moments);
	cout << "Moments: " << endl;
	print(moments);

	//2. Решаем СЛАУ
	vector < vector < double> > M;
	M.resize(3);

	for (int i = 0; i < 3; ++i){
		for (int j = 0; j < 3; ++j){
			M[i].push_back(moments[i + j]);
		}
	}

	vector <double> mju;
	for (int s = 0; s < 3; ++s){
		mju.push_back(-moments[3+s]);
	}
	
	vector < vector<double> >P;
	vector < vector<double> >L;
	vector < vector<double> >U;

	P.resize(3);
	L.resize(3);
	U.resize(3);
	
	for (int i = 0; i < 3; ++i){
		P[i].resize(3);
		L[i].resize(3);
		U[i].resize(3);
	}

	PLU(M, P, L, U);

	vector <double> a;
	a.resize(3);
	copy( a, solveLinerSystem(P, L, U, mju));
	cout << "Print a" << endl;
	print(a);
	//3. Находим корни узлового многочлена
	vector <double> x;
	
	x.push_back(MethodNewton(-58.0, -50.0, a));
	x.push_back(MethodNewton ( 0.0, 0.6, a));
	x.push_back(MethodNewton(1.0, 2.0, a));
	cout << "x" << endl;
	print(x);

	// 4. Находим Ai
	vector <vector <double> > X;
	X.resize(3);
	for (int i = 0; i < 3; ++i){
		for (int j = 0; j < 3; ++j){
			X[i].push_back( pow (x[j], i));
		}
	}

	cout << "Matrix X:" << endl;
	print(X);

	PLU( X, P, L, U );
	
	vector < double > A;
	A.resize(3);
	
	for (int i = 0; i < 3; ++i){
		mju[i] = moments[i];
	}
	cout << "mju " << endl;
	print(mju);
	copy(A, solveLinerSystem(P, L, U, mju));
	cout << "Matrix A: " << endl;
	print(A);

	double res = A[0] * x[0] + A[1] * x[1] + A[2] * x[2];
	cout << res << endl;
	if (fabs(res) < EPS)
		res = 0.0;
	return res;
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
double MethodNewton(const double l, const double r, const vector <double> a){
	double res = (r+l)/2;

	while ( fabs(w_f(res, a)) > EPS){
		res -= w_f(res, a) / w_d(res, a);
	}
	/*
	if (fabs(res) < EPS)
		res = 0.0;
	*/
	return res;
}
double calculateIntegralWithStepNC(const double l, const double r, const double h){
	double res = 0;
	for (int i = 0; i < h; ++i)
		res += NewtonCotes(l + i * ((r - l) / h), l + (i + 1)*((r - l) / h));
	return res;
}

double calculateIntegralNC(const double l, const double r, const int n){
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
double w_f ( const double x, const vector <double> a){
	double res = pow(x, 3) + a[0] * pow(x, 2) + a[1] * x + a[2];
	if (fabs(res) < EPS)
		res = 0.0;
	return res;
}

double w_d ( double x, vector <double> a){
	double res = 3 * pow(x, 2) + 2 * a[0] * x + a[1];
	if (fabs(res) < EPS)
		res = 0.0;
	return res;
}