#include "matrix.h"
using namespace std;

const double ACCURACY = 10e-12;
//-------------------------------
double f1(double);
double d1(double);

double MethodNuton(double, double);
void MethodNuton(double, double, double, double, double&, double&);

double f2(double, double);
double f3(double, double);

vector <double> calculateMatrix(const vector <double> &);
vector <vector <double> > calculateMatrixJacobi(const vector <double> &);
vector <double> MethodNuton(const vector <double>&, const bool, const int);

//-------------------------------
int main(){
	/*
	double solution = MethodNuton(0.0, 1.0);
	cout << "1) Solution of nonlinear equation:" << endl;
	cout << solution << endl << endl;
	cout << "Check:" << endl;
	cout << "f(" << solution << ") = " << f1(solution) << endl << endl;
	
	double xn = 0, yn = 0 ;

	MethodNuton(3.0, 4.0, 1.0, 2.0, xn, yn);
	cout << "2)Solution of nonlinear system:" << endl;
	cout << xn << " " << yn << endl;

	cout << "Check: " << endl;
	cout << f2(xn, yn) << " " << f3(xn, yn) << endl;
	
	
	cout << "===============" << endl << endl;
	cout << "Solution of task2:" << endl;
	*/
	vector <double> v;
	v.push_back(0.5);
	v.push_back(0.5);
	v.push_back(1.5);
	v.push_back(-1.0);
	v.push_back(-0.5);
	v.push_back(1.5);
	v.push_back(0.5);
	v.push_back(-0.5);
	v.push_back(1.5);
	v.push_back(-1.5);

	print( MethodNuton(v, true, 5 ));

	system("pause");
	return 0;
}
// f(x) = sqrt(x) - cos(x)
double f1(double x){

	if ((sqrt(x) - cos(x)) < EPS) return 0.0;
	return (sqrt(x) - cos(x));
}

// f'(x) = - 1/(2*sqrt(x)) + sin(x)
double d1(double x){
	if ((1.0 / (2.0*sqrt(x)) + sin(x)) < EPS) return 0.0;
	return (1.0 / (2.0*sqrt(x)) + sin(x));
}

double MethodNuton(const double a, const double b){
	double x = b - f1(b) / d1(b);
	if (fabs(x - b) > ACCURACY){
		if (fabs(x - a) > EPS) x = MethodNuton(a, x);
	}
	if (x < EPS){
		x = 0.0;
	}
	return x;
}

double f2(double x, double y){
	double res = cos(x - 1) + y - 0.5;
	if (fabs(res) < EPS)
		return 0.0;
	else
		return res;
}

double f3(double x, double y){
	double res = x - cos(y) - 3.0;
	if (fabs(res) < EPS)
		return 0.0;
	else
		return res;
}

void MethodNuton(double a0, double a1, double c0, double c1, double& xn, double& yn){
	vector <double> F;
	F.resize(2);

	F[0] = f2 ( a1, c1);
	F[1] = f3 ( a1, c1);
	//cout << "Matrix F" << F[0] << " " << F[1] << endl;

	vector <vector <double> > J;
	vector <vector <double> > P;                                              
	vector <vector <double> > L;
	vector <vector <double> > U;
	J.resize(2);
	P.resize(2);
	L.resize(2);
	U.resize(2);

	for (int i = 0; i < 2; ++i){
		J[i].resize(2);
		P[i].resize(2);
		L[i].resize(2);
		U[i].resize(2);
	}

	J[0][0] = -sin(a1 - 1);
	J[0][1] = 1;
	J[1][0] = 1;
	J[1][1] = sin(c1);
	checkEps(J);

	//cout << "Matrix Jacobi:" << endl;
	//print(J);
	cout << endl;

	PLU(J, P, L, U);

	vector <vector <double> > InvJ;
	InvJ.resize(2);
	for (int i = 0; i < 2; ++i)
		InvJ[i].resize(2);

	copy( InvJ, getInverseMatrix(P, L, U));
	//cout << "Inverse Matrix:" << endl;
	//print(InvJ);

	vector<double> M;
	M.resize(F.size());
	M = multiplication(InvJ, F);
	xn = a1 - M[0];
	yn = c1 - M[1];

	if (fabs(xn) < EPS) xn = 0.0;
	if (fabs(yn) < EPS) yn = 0.0;

	//cout << "============ " << endl;
	//cout << xn << " " << yn << endl;

	if (fabs(a1 - xn) > ACCURACY && fabs(c1 - yn) > ACCURACY){
		MethodNuton(a0, xn, c0, yn, xn, yn);
	}
	return;
}

vector <double> calculateMatrix (const vector <double> & X){
	vector <double> res;
	res.resize(X.size());

	res[0] = -(cos(X[0] * X[1]) - exp(-3 * X[2]) + X[3] * pow(X[4], 2) - X[5] - sinh(2 * X[7]) * X[8] + 2 * X[9] + 2.0004339741653854440);
	res[1] = -(sin(X[0] * X[1]) + X[2] * X[8] * X[6] - exp(-X[9] + X[5]) + 3 * pow(X[4], 2) - X[5] * (X[7] + 1) + 10.886272036407019994);
	res[2] = -(X[0] - X[1] + X[2] - X[3] + X[4] - X[5] + X[6] - X[7] + X[8] - X[9] - 3.1361904761904761904);
	res[3] = -(2 * cos(-X[8] + X[3]) + X[4] / (X[2] + X[0]) - sin(pow(X[1], 2)) + pow(cos(X[6] * X[9]), 2) - X[7] - 0.1707472705022304757);
	res[4] = -(sin(X[4]) + 2 * X[7] * (X[2] + X[0]) - exp(-X[6] * (-X[9] + X[5])) + 2 * cos(X[1]) - 1 / (X[3] - X[8]) - 0.3685896273101277862);
	res[5] = -(exp(X[0] - X[3] - X[8]) + pow(X[4], 2) / X[7] + cos(3 * X[9] * X[1]) / 2 - X[5] * X[2] + 2.0491086016771875115);
	res[6] = -(pow(X[1], 3) * X[6] - sin(X[9] / X[4] + X[7]) + (X[0] - X[5]) * cos(X[3]) + X[2] - 0.7380430076202798014);
	res[7] = -(X[4] * pow(X[0] - 2 * X[5], 2) - 2 * sin(-X[8] + X[2]) + 1.5 * X[3] - exp(X[1] * X[6] + X[9]) + 3.5668321989693809040);
	res[8] = -(7 / X[5] + exp(X[4] + X[3]) - 2 * X[1] * X[7] * X[9] * X[6] + 3 * X[8] - 3 * X[0] - 8.4394734508383257499);
	res[9] = -(X[9] * X[0] + X[8] * X[1] - X[7] * X[2] + sin(X[3] + X[4] + X[5]) * X[6] - 0.78238095238095238096);

	checkEps(res);
	return res;
}

vector <vector <double> > calculateMatrixJacobi(const vector <double> & X){
	vector <vector <double> > res;
	res.resize(10);
	for (int i = 0; i < 10; ++i){
		res[i].resize(10);
	}
	res[0][0] = -sin(X[0] * X[1]) * X[1];
	res[0][1] = -sin(X[0] * X[1]) * X[0];
	res[0][2] = 3.0 * exp(-3.0 * X[2]);
	res[0][3] = X[4] * X[4];
	res[0][4] = 2.0 * X[3] * X[4];
	res[0][5] = -1.0;
	res[0][6] = 0.0;
	res[0][7] = -2.0 * cosh(2.0 * X[7]) * X[8];
	res[0][8] = -sinh(2.0 * X[7]);
	res[0][9] = 2.0;
	res[1][0] = cos(X[0] * X[1]) * X[1];
	res[1][1] = cos(X[0] * X[1]) * X[0];
	res[1][2] = X[8] * X[6];
	res[1][3] = 0.0;
	res[1][4] = 6.0 * X[4];
	res[1][5] = -exp(-X[9] + X[5]) - X[7] - 1.0;
	res[1][6] = X[2] * X[8];
	res[1][7] = -X[5];
	res[1][8] = X[2] * X[6];
	res[1][9] = exp(-X[9] + X[5]);
	res[2][0] = 1.0;
	res[2][1] = -1.0;
	res[2][2] = 1.0;
	res[2][3] = -1.0;
	res[2][4] = 1.0;
	res[2][5] = -1.0;
	res[2][6] = 1.0;
	res[2][7] = -1.0;
	res[2][8] = 1.0;
	res[2][9] = -1.0;
	res[3][0] = -X[4] / pow(X[2] + X[0], 2.0);
	res[3][1] = -2.0 * cos(X[1] * X[1]) * X[1];
	res[3][2] = -X[4] / pow(X[2] + X[0], 2.0);
	res[3][3] = -2.0 * sin(-X[8] + X[3]);
	res[3][4] = pow(X[2] + X[0], -1.0);
	res[3][5] = 0.0;
	res[3][6] = -2.0 * cos(X[6] * X[9]) * sin(X[6] * X[9]) * X[9];
	res[3][7] = -1.0;
	res[3][8] = 2.0 * sin(-X[8] + X[3]);
	res[3][9] = -2.0 * cos(X[6] * X[9]) * sin(X[6] * X[9]) * X[6];
	res[4][0] = 2.0 * X[7];
	res[4][1] = -2.0 * sin(X[1]);
	res[4][2] = 2.0 * X[7];
	res[4][3] = pow(-X[8] + X[3], -2.0);
	res[4][4] = cos(X[4]);
	res[4][5] = X[6] * exp(-X[6] * (-X[9] + X[5]));
	res[4][6] = -(X[9] - X[5]) * exp(-X[6] * (-X[9] + X[5]));
	res[4][7] = 2.0 * X[2] + 2.0 * X[0];
	res[4][8] = -pow(-X[8] + X[3], -2);
	res[4][9] = -X[6] * exp(-X[6] * (-X[9] + X[5]));
	res[5][0] = exp(X[0] - X[3] - X[8]);
	res[5][1] = -3.0 / 2.0 * sin(3.0 * X[9] * X[1]) * X[9];
	res[5][2] = -X[5];
	res[5][3] = -exp(X[0] - X[3] - X[8]);
	res[5][4] = 2.0 * X[4] / X[7];
	res[5][5] = -X[2];
	res[5][6] = 0.0;
	res[5][7] = -pow(X[4], 2.0) / pow(X[7], 2.0);
	res[5][8] = -exp(X[0] - X[3] - X[8]);
	res[5][9] = -3.0 / 2.0 * sin(3.0 * X[9] * X[1]) * X[1];
	res[6][0] = cos(X[3]);
	res[6][1] = 3.0 * X[1] * X[1] * X[6];
	res[6][2] = 1.0;
	res[6][3] = -(X[0] - X[5]) * sin(X[3]);
	res[6][4] = cos(X[9] / X[4] + X[7]) * X[9] * pow(X[4], -2.0);
	res[6][5] = -cos(X[3]);
	res[6][6] = pow(X[1], 3.0);
	res[6][7] = -cos(X[9] / X[4] + X[7]);
	res[6][8] = 0.0;
	res[6][9] = -cos(X[9] / X[4] + X[7]) / X[4];
	res[7][0] = 2.0 * X[4] * (X[0] - 2.0 * X[5]);
	res[7][1] = -X[6] * exp(X[1] * X[6] + X[9]);
	res[7][2] = -2.0 * cos(-X[8] + X[2]);
	res[7][3] = 1.5;
	res[7][4] = pow(X[0] - 2.0 * X[5], 2.0);
	res[7][5] = -4.0 * X[4] * (X[0] - 2.0 * X[5]);
	res[7][6] = -X[1] * exp(X[1] * X[6] + X[9]);
	res[7][7] = 0.0;
	res[7][8] = 2.0 * cos(-X[8] + X[2]);
	res[7][9] = -exp(X[1] * X[6] + X[9]);
	res[8][0] = -3.0;
	res[8][1] = -2.0 * X[7] * X[9] * X[6];
	res[8][2] = 0.0;
	res[8][3] = exp(X[4] + X[3]);
	res[8][4] = exp(X[4] + X[3]);
	res[8][5] = -7.0 * pow(X[5], -2.0);
	res[8][6] = -2.0 * X[1] * X[7] * X[9];
	res[8][7] = -2.0 * X[1] * X[9] * X[6];
	res[8][8] = 3.0;
	res[8][9] = -2.0 * X[1] * (X[7] * X[6]);
	res[9][0] = X[9];
	res[9][1] = X[8];
	res[9][2] = -X[7];
	res[9][3] = cos(X[3] + X[4] + X[5]) * X[6];
	res[9][4] = cos(X[3] + X[4] + X[5]) * X[6];
	res[9][5] = cos(X[3] + X[4] + X[5]) * X[6];
	res[9][6] = sin(X[3] + X[4] + X[5]);
	res[9][7] = -X[2];
	res[9][8] = X[1];
	res[9][9] = X[0];

	checkEps(res); 

	return res;
}

vector <double> MethodNuton( const vector <double>& v, const bool is_modified, const int index ){
	int iter = 0;
	
	vector < double> res;
	res.resize(10);
	copy(res, v);

	vector <double> F;
	F.resize(10);
	
	vector <vector <double> > J;
	J.resize(10);
	for (int i = 0; i < 10; ++i){
		J[i].resize(10);
	}
	zero(J);
	
	vector <vector <double> > P;
	vector <vector <double> > L;
	vector <vector <double> > U;

	P.resize(10);
	L.resize(10);
	U.resize(10);

	for (int i = 0; i < 10; ++i){
		P[i].resize(10);
		L[i].resize(10);
		U[i].resize(10);
	}

	vector <double> d; //d =x^(k+1) - x^(k) 
	d.resize(10);
	bool first = true;
	unsigned int k = 0;
	do{
		copy(F, calculateMatrix(res));

		if ((!is_modified) || (is_modified && k < index)){			//k % index == 0 )){
			first = false;
			copy(J, calculateMatrixJacobi(res));
			PLU(J, P, L, U);
		}

		copy(d, solveLinerSystem(P, L, U, F));

		checkEps(d);
		for (int i = 0; i < v.size(); ++i){
			res[i] += d[i];
		}
		++k;
		++iter;
	} while ( norm(d) > EPS );
	cout << endl;
	cout << k;
	cout << endl;
	return res;
}

