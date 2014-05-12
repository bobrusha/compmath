#define _MATRIX_H_
#ifdef _MATRIX_H_

#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

const double EPS = 10e-9;
void print(const vector <double> &);
void print(const vector < vector <double> >&);
void copy(vector <double> &, const vector <double> &);
void copy(vector < vector <double> >&, const vector <vector <double>> &);
void checkEps(vector <double> &);
void checkEps(vector <vector <double> > &);
void zero(vector <vector <double> > &);

vector <double> multiplication(const vector <vector <double> >& l, const vector < double >& r){
	vector <double> res;
	res.resize(r.size());

	for (int i = 0; i < r.size(); ++i){
		res[i] = 0;
		for (int j = 0; j < r.size(); ++j){
			res[i] += l[i][j] * r[j];
		}
	}
	checkEps(res);
	return res;
}
vector < vector <double> > multiplication(const vector < vector <double> >& l, const vector < vector <double> >&  r)
{
	vector <vector <double> > res;
	if (l.size() == r[0].size()){


		res.resize(r.size());
		for (int i = 0; i < l.size(); ++i){
			res[i].resize(l[0].size());
		}

		for (int i = 0; i < l[0].size(); ++i)
		{
			for (int j = 0; j < l[0].size(); ++j)
			{
				for (int k = 0; k < l[0].size(); ++k)
				{
					res[i][j] += l[i][k] * r[k][j];
				}
				if (fabs(res[i][j]) < EPS) { res[i][j] = 0.0; }
			}
		}
		checkEps(res);
		return res;
	}
	return res;
}

vector < vector <double> > addition(const vector < vector <double> >& l, const vector < vector <double> >&  r)
{
	vector <vector <double> > res;
	res.resize(r.size());
	for (int i = 0; i < l.size(); ++i){
		res[i].resize(l[0].size());
	}

	for (int i = 0; i < l.size(); ++i){
		for (int j = 0; j < l.size(); ++j){
			res[i][j] = l[i][j] + r[i][j];
		}
	}
	return res;
}

vector < vector <double> > subtraction(vector < vector <double> >& l, vector < vector <double> >&  r)
{
	vector <vector <double> > res;
	res.resize(r.size());
	for (int i = 0; i < l.size(); ++i){
		res[i].resize(l[0].size());
	}

	for (int i = 0; i < l.size(); ++i){
		for (int j = 0; j < l.size(); ++j){
			res[i][j] = l[i][j] - r[i][j];
		}
	}
	return res;
}

int calculateRank(vector < vector <double> >& x)
{
	int r = 0;

	for (int i = 0; i < x[0].size(); ++i)
	{
		for (int j = 0; j < x.size(); ++j)
		{
			if (x[i][j]) {
				++r;
				break;
			}
		}
	}
	return r;
}

vector < vector <double> > identityMatrix(const int n){
	vector < vector <double> > res;
	res.resize(n);
	for (int i = 0; i < n; ++i){
		res[i].resize(n);
	}
	for (int i = 0; i < n; ++i){
		for (int j = 0; j < n; ++j){
			res[i][j] = 0;
		}
		res[i][i] = 1;
	}
	return res;
}

void swapRows(vector < vector <double> >& x, const int r1, const int r2){
	int sz = x.size();

	if (r1 >= sz || r2 >= sz){
		cout << "Error in swapRows: one of numbers bigger than numbers of rows" << endl;
		return;
	}
	double tmp;

	for (int i = 0; i < x.size(); ++i){
		tmp = x[i][r1];
		x[i][r1] = x[i][r2];
		x[i][r2] = tmp;
	}
	return;
}

void swapColumns(vector < vector <double> >& x, const int r1, const int r2){
	int sz = x[0].size();

	if (r1 >= sz || r2 >= sz)
	{
		cout << "Error in swapRows: one of numbers bigger than numbers of rows" << endl;
		return;
	}
	double tmp;

	for (int i = 0; i < x.size(); ++i)
	{
		tmp = x[r1][i];
		x[r1][i] = x[r2][i];
		x[r2][i] = tmp;
	}
	return;
}

void PLU(const vector < vector <double> > &A, vector < vector <double> >& P, vector < vector <double> > &L, vector < vector <double> > &U) {
	//n - размерность исходной матрицы

	const unsigned int n = A.size();

	vector <int> pArray;
	pArray.resize(n);

	double pivot;

	unsigned int k1 = 0;
	unsigned int i;

	double detP = 1;
	double detU = 1;

	vector <vector <double> >C;
	C.resize(n);

	for (i = 0; i < n; ++i){
		pArray[i] = i;
		C[i].resize(n);
	}

	copy(C, A);

	//загружаем в матрицу P единичную матрицу
	zero(L);
	zero(U);

	copy(P, identityMatrix(A.size()));

	for (int k = 0; k < n; ++k) {
		//поиск опорного элемента
		pivot = 0;
		for (i = k; i < n; ++i) {
			if (fabs(C[i][k]) > pivot) {
				pivot = fabs(C[i][k]);
				k1 = i;
				detP = -detP;
			}
		}

		if (fabs(pivot) < EPS) {
			cout << "Error: Matrix is singular" << endl;
			cout << k << endl;
			continue;
		}

		int tmp = pArray[k];
		pArray[k] = pArray[k1];
		pArray[k1] = tmp;

		for (i = 0; i < n; ++i){
			double tmp1 = C[k][i];
			C[k][i] = C[k1][i];
			C[k1][i] = tmp1;
		}

		for (int i = k + 1; i < n; ++i) {
			C[i][k] /= C[k][k];
			for (int j = k + 1; j < n; ++j)
				C[i][j] -= C[i][k] * C[k][j];
		}
	}
	//--------------------------------------------------------------------
	for (i = 0; i < n; ++i){
		for (int j = 0; j < n; ++j){
			P[i][j] = 0;
			if (i < j){
				U[i][j] = C[i][j];
				L[i][j] = 0;
			}
			else if (i > j){
				L[i][j] = C[i][j];
				U[i][j] = 0;
			}
			else if (i == j){
				U[i][j] = C[i][j];
				detU *= C[i][j];
				L[i][j] = 1;
			}
		}
	}
	for (i = 0; i < n; ++i){
		P[i][pArray[i]] = 1;
	}
	
	return;
}

int GaussJordanElimination(vector <vector < double > >& A, vector < double > & b){
	vector <vector < double > > TMP;
	TMP.resize(A.size());

	for (int i = 0; i < A.size(); ++i)
		TMP[i].resize(A[0].size() + 1);

	int sz = A[0].size();

	for (int i = 0; i < sz; ++i){
		for (int j = 0; j < sz; ++j){
			TMP[i][j] = A[i][j];
		}
	}

	for (int i = 0; i < sz; i++){
		TMP[i][sz] = b[i];
	}

	std::vector<int> where(sz, -1);

	const int n = A[0].size();
	const int m = A.size();

	for (int col = 0, row = 0; col < m && row < n; ++col){
		int sel = row;
		for (int i = row; i < n; i++){
			if (fabs(TMP[i][col]) > fabs(TMP[sel][col])){
				sel = i;
			}
		}
		if (fabs(TMP[sel][col]) < EPS)
			continue;
		for (int i = col; i <= m; i++){
			swapColumns(TMP, sel, i);
		}

		where[col] = row;

		for (int i = 0; i < n; i++){
			if (i != row){
				double tmp = TMP[i][col] / TMP[row][col];
				for (int j = col; j <= m; ++j){
					TMP[i][j] -= TMP[row][j] * tmp;
				}
			}
		}
		++row;
	}
	for (int i = 0; i < b.size(); ++i){
		b[i] = 0;
	}

	for (int i = 0; i < m; i++){
		if (where[i] != -1){
			b[i] = TMP[where[i]][m] / TMP[where[i]][i];
		}
	}

	for (int i = 0; i < A[0].size(); i++){
		for (int j = 0; j < A.size(); j++){
			A[i][j] = TMP[i][j];
		}
	}

	for (int i = 0; i < n; ++i) {
		double sum = 0.0;
		for (int j = 0; j < m; ++j)
			sum += b[j] * TMP[i][j];
		if (fabs(sum - TMP[i][m]) > EPS){
			return 0;
		}
	}

	for (int i = 0; i < m; i++){
		if (where[i] == -1){
			return 2;
		}
	}
	return 1;
}

vector <double> solveLinerSystem(const vector < vector <double> > & P, const vector <vector <double> > & L, const vector <vector <double> > & U, vector <double>& b){
	const unsigned int n = P.size();
	double det = 1;
	for (int i = 0; i < n; ++i){
		det *= U[i][i];
	}

	if (det){
		vector <double> x, y, p;
		p.resize(n);

		double sum = 0;
		for (int i = 0; i < n; ++i){
			x.push_back(0);
			y.push_back(0);
			for (int j = 0; j < n; ++j){
				if (P[i][j] == 1){
					p[i] = j;
				}
			}
		}
		for (int i = 0; i < n; ++i){
			sum = 0;
			for (int j = 0; j < i; ++j){
				sum += L[i][j] * y[j];
			}
			sum = b[p[i]] - sum;
			y[i] = sum;
		}
		for (int i = n - 1; i >= 0; --i){
			sum = 0;
			for (int j = i + 1; j < n; ++j){
				sum += U[i][j] * x[j];
			}
			sum = (y[i] - sum) / U[i][i];
			x[i] = sum;
		}
		return x;
	}
	else{
		cout << "ERROR!!!" << endl;
		return b;
	}
}

vector <vector <double> > getInverseMatrix(const vector <vector <double> > & P, const vector < vector <double> >& L, const vector <vector <double>> & U)
{
	vector <vector <double> > Inv;
	Inv.resize(L[0].size());
	for (int i = 0; i < L[0].size(); ++i){
		Inv[i].resize(L[0].size());
	}
	zero(Inv);

	for (int i = 0; i < L[0].size(); ++i)
	{
		vector <double> Ei;
		Ei.resize(L[0].size());
		for (int j = 0; j < L[0].size(); ++j){
			Ei[j] = 0;
		}
		Ei[i] = 1.0;

		vector <double> Ai;

		Ai.resize(L[0].size());

		for (int j = 0; j < L[0].size(); ++j){
			Ai[j] = 0;
		}

		copy(Ai, solveLinerSystem(P, L, U, Ei));

		for (int j = 0; j < L.size(); ++j){
			Inv[i][j] = Ai[j];
			if (fabs(Inv[i][j]) < EPS) Inv[i][j] = 0.0;
		}
	}

	return multiplication(P, Inv);
}


void print(const vector <double> & x){
	for (int i = 0; i < x.size(); ++i){
		cout << x[i] << endl;
	}
	cout << endl;
	return;
}

void print(const vector < vector <double> >& x){
	for (int i = 0; i < x.size(); ++i){
		for (int j = 0; j < x[0].size(); ++j){
			cout << x[i][j] << " ";
		}
		cout << endl;
	}
	return;
}

void copy(vector < vector <double> >& l, const vector <vector <double>> & r){
	for (int i = 0; i < l[0].size(); ++i){
		for (int j = 0; j < l.size(); ++j){
			l[i][j] = r[i][j];
		}
	}
}
void copy(vector <double> & l, const vector <double> & r){
	for (int j = 0; j < l.size(); ++j){
		l[j] = r[j];
	}
}
void checkEps(vector <double> & l){
	for (int i = 0; i < l.size(); ++i)
	if (fabs(l[i]) < EPS)
		l[i] = 0;
	return;
}
void checkEps(vector <vector <double> > & l){
	for (int i = 0; i < l[0].size(); ++i){
		for (int j = 0; j < l.size(); ++j){
			if (fabs(l[i][j]) < EPS)
				l[i][j] = 0;
		}
	}
	return;
}
void zero(vector <vector <double> > & l){
	for (int i = 0; i < l[0].size(); ++i){
		for (int j = 0; j < l.size(); ++j){
			l[i][j] = 0;
		}
	}
	return;
}

double norm(vector <double> & A){
	double res = fabs(A[0]);
	for (int i = 0; i < A.size(); ++i){
		if (fabs(A[i]) > fabs(res)) res = fabs(A[i]);
	}
	return res;
}

#endif	