#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

//stale
double c, alfa, q, to;
double s = 2;
double k = 50;
//zmienne
double l;
int iw;
double de;
double *ttab;
double *ptab;
double **htab;
double **hptab;
//struktura elementu
struct element {
	double cc;
	double hloc[2][2];
	element() {
		cc = s*k / de;
		hloc[0][0] = cc;
		hloc[0][1] = -cc;
		hloc[1][0] = -cc;
		hloc[1][1] = cc;
	}
};

void wypelnijStale(){
	//takie jak w pdf
	c = 40;
	alfa = 10;
	q = -150;
	to = 400;
}

double dlugoscElementu() {
	return l / (iw-1);
}

void wypelnijPusteTablice() {
	//nowy wektor nieznanych temperatur
	ttab = new double[iw];
	for (int i = 0; i < iw; i++)
		ttab[i] = 0;
	//nowy wektor p i jego wypelnienie
	ptab = new double[iw];
	for (int i = 0; i < iw; i++)
		ptab[i] = 0;
	//nowa macierz h
	htab = new double *[iw];
	for (int i = 0; i<iw; i++)
		htab[i] = new double[iw];
	for (int i = 0; i < iw; i++) {
		for (int j = 0; j < iw; j++)
			htab[i][j] = 0;
	}
}

void wypelnijPtab() {
	ptab[0] = q*s;
	ptab[iw - 1] = to*s*alfa*(-1);
}

void wypelnijHtab() {
	element *etab = new element[iw - 1];
	
	for (int i = 0; i < iw - 1; i++) {
		htab[i][i] += etab[i].hloc[0][0];
		htab[i][i+1] += etab[i].hloc[0][1];
		htab[i+1][i] += etab[i].hloc[1][0];
		htab[i+1][i+1] += etab[i].hloc[1][1];
	}
	htab[iw-1][iw-1] += alfa*s;
}

void zlaczHiP() {
	//nowa macierz hp
	hptab = new double *[iw];
	for (int i = 0; i<iw+1; i++)
		hptab[i] = new double[iw];
	for (int i = 0; i < iw; i++) {
		for (int j = 0; j < iw+1; j++)
			hptab[i][j] = 0;
	}

	//dodanie h tab
	for (int i = 0; i < iw; i++) {
		for (int j = 0; j < iw; j++)
			hptab[i][j] = htab[i][j];
	}

	//dodanie p tab
	for (int i = 0; i < iw; i++)
		hptab[i][iw] = ptab[i];
}

vector<double> gauss() {
	int n = iw;

	for (int i = 0; i<n; i++) {
		// Search for maximum in this column
		double maxEl = abs(hptab[i][i]);
		int maxRow = i;
		for (int k = i + 1; k<n; k++) {
			if (abs(hptab[k][i]) > maxEl) {
				maxEl = abs(hptab[k][i]);
				maxRow = k;
			}
		}

		// Swap maximum row with current row (column by column)
		for (int k = i; k<n + 1; k++) {
			double tmp = hptab[maxRow][k];
			hptab[maxRow][k] = hptab[i][k];
			hptab[i][k] = tmp;
		}

		// Make all rows below this one 0 in current column
		for (int k = i + 1; k<n; k++) {
			double c = -hptab[k][i] / hptab[i][i];
			for (int j = i; j<n + 1; j++) {
				if (i == j) {
					hptab[k][j] = 0;
				}
				else {
					hptab[k][j] += c * hptab[i][j];
				}
			}
		}
	}

	// Solve equation Ax=b for an upper triangular matrix A
	vector<double> x(n);
	for (int i = n - 1; i >= 0; i--) {
		x[i] = -(hptab[i][n] / hptab[i][i]);
		for (int k = i - 1; k >= 0; k--) {
			hptab[k][n] += hptab[k][i] * x[i];
		}
	}
	return x;
}

void wypiszTemp() {
	cout << "Wyznaczone temperatury:" << endl;
	for (int i = 0; i < iw; i++)
		cout << "t[" << i + 1 << "]= " << ttab[i] << endl;
}

void test1() {
	cout << endl;
	cout << "Wektor P:" << endl;
	for (int i = 0; i < iw; i++)
		cout <<"p["<<i+1<<"]="<<ptab[i] << " ";
	cout << endl;
}

void test2() {
	cout << endl;
	cout << "Macierz Globalna H:" << endl;
	for (int i = 0; i < iw; i++) {
		for (int j = 0; j < iw; j++)
			cout << htab[i][j]<< " ";
		cout << endl;
	}	
	cout << endl;
}

void test3() {
	cout << "Macierz Globalna HP:" << endl;
	for (int i = 0; i < iw; i++) {
		for (int j = 0; j < iw+1; j++)
			cout << hptab[i][j] << " ";
		cout << endl;
	}
	cout << endl;
}


int main() {
	//l=5, iw=3 430 422,5 415
	cout << "Dlugosc preta: ";
	cin >> l;
	cout << "Ilosc wezlow: ";
	cin >> iw;

	de = dlugoscElementu();
	wypelnijStale();
	wypelnijPusteTablice();
	wypelnijPtab();
	wypelnijHtab();
	zlaczHiP();

	cout << endl;
	test1();
	test2();
	test3();


	
	vector<double> tempe(iw);
	tempe= gauss();
	cout << "Wyznaczone temperatury:" << endl;
	for (int i = 0; i < iw; i++)
		cout << "t[" << i + 1 << "]= " << tempe[i] << endl;
	
	system("pause");
}
