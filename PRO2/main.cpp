#include <iostream>
//#include <cmath>
#include <vector>
#include <fstream>
using namespace std;

//stale z pliku
int iw;
double rmax, taldelta, talmax, c, ro, k, alfa, t0, talfa;
//zmienne
double rdelta;
// int iteracje (taumax/taudelta)+1
fstream plik;
double *ttab;
double *ptab;
double **htab;
double **hptab;
//kwadratura
double w1 = 1, w2 = 1;
double e1 = -0.5773502692, e2 = 0.5773502692;
double n11 = 0.5*(1 - e1);
double n12 = 0.5*(1 - e2);
double n21 = 0.5*(1 + e1);
double n22 = 0.5*(1 + e2);

struct element {
	double hloc[2][2];
	double ploc[2];

	element() {
		hloc[0][0] = 0;
		hloc[0][1] = 0;
		hloc[1][0] = 0;
		hloc[1][1] = 0;
		ploc[0] = 0;
		ploc[1] = 0;
	}

	element(double r1, double r2, double temp0, double temp1) {
		double rp1 = n11*r1 + n21*r2;
		double rp2 = n12*r1 + n22*r2;
	
		hloc[0][0] = k / rdelta*(rp1 + rp2) + c*ro*rdelta / taldelta*(n11*n11*rp1 + n12*n12*rp2);
		hloc[0][1] = -k / rdelta*(rp1 + rp2) + c*ro*rdelta / taldelta*(n11*n21*rp1 + n12*n22*rp2);
		hloc[1][0] = hloc[0][1];
		hloc[1][1] = k / rdelta*(rp1 + rp2) + c*ro*rdelta / taldelta*(n21*n21*rp1 + n22*n22*rp2);
		//cout << endl << hloc[0][0] << " " << hloc[0][1] << endl << hloc[1][0] << " " << hloc[1][1] << endl;
		
		double tp1 = n11*temp0 + n21*temp1;
		double tp2 = n12*temp0 + n22*temp1;

		ploc[0] = c*ro*rdelta / taldelta*(tp1*rp1*n11 + tp2*rp2*n12);
		ploc[1] = c*ro*rdelta / taldelta*(tp1*rp1*n21 + tp2*rp2*n22);
		cout << endl << ploc[0] << endl;
	}
};

element *etab;

void wczytajDaneWejsciowe() {
	plik.open("dane.txt", ios::in);
	if (plik.good()) {
		plik >> iw >> rmax >> taldelta >> talmax>> c >> ro >> k >> alfa >> t0 >> talfa;
		plik.close();
		rdelta = rmax / (iw - 1);
		cout << "Dane wejsciowe zostaly poprawnie odczytane\n";
	}
	else
		cout << "Blad odczytu pliku 'dane.txt'\n";
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

void wypelnijHtab() {
	etab = new element[iw - 1];
	
	for (int i = 0; i < iw - 1; i++)
		etab[i] = element(i*rdelta, (i+1)*rdelta, t0, talfa);

	for (int i = 0; i < iw - 1; i++) {
		htab[i][i] += etab[i].hloc[0][0];
		htab[i][i + 1] += etab[i].hloc[0][1];
		htab[i + 1][i] += etab[i].hloc[1][0];
		htab[i + 1][i + 1] += etab[i].hloc[1][1];
		//added
		ptab[i] += etab[i].ploc[0];
		ptab[i + 1] += etab[i].ploc[1];
	}
	htab[iw - 1][iw - 1] += 2*alfa*rmax;
	//added
	ptab[iw - 1] -= 2 * alfa*rmax*talfa;
}

void zlaczHiP() {
	//nowa macierz hp
	hptab = new double *[iw];
	for (int i = 0; i<iw + 1; i++)
		hptab[i] = new double[iw];
	for (int i = 0; i < iw; i++) {
		for (int j = 0; j < iw + 1; j++)
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
		// szuka max
		double maxEl = abs(hptab[i][i]);
		int maxRow = i;
		for (int k = i + 1; k<n; k++) {
			if (abs(hptab[k][i]) > maxEl) {
				maxEl = abs(hptab[k][i]);
				maxRow = k;
			}
		}

		// zmiana wierszy max z obecnym
		for (int k = i; k<n + 1; k++) {
			double tmp = hptab[maxRow][k];
			hptab[maxRow][k] = hptab[i][k];
			hptab[i][k] = tmp;
		}

		// zerowanie
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

	// wyznacz szukany wektor
	vector<double> x(n);
	for (int i = n - 1; i >= 0; i--) {
		x[i] = -(hptab[i][n] / hptab[i][i]);
		for (int k = i - 1; k >= 0; k--) {
			hptab[k][n] += hptab[k][i] * x[i];
		}
	}
	return x;
}

void wypiszTempiZapisz() {
	vector<double> tempe(iw);

	tempe = gauss();

	cout << "Wyniki:" << endl;
	for (int i = 0; i < iw; i++)
		cout << "t[" << i + 1 << "]=" << tempe[i]<<"   ";
	cout << endl;

	ofstream plik2;
	plik2.open("wyniki.txt");
	if (plik2.good()) {
		for (int i = 0; i < iw; i++)
			plik2 << "t[" << i + 1 << "]=" << tempe[i] << "   ";
		plik2 << endl;
		plik2.close();
	}
	else
		cout << "\nBlad zapisu pliku 'wyniki.txt'\n\n";
}

void test0() {
	cout << endl;
	//iw rmax, taldelta, talmax, c, ro, k, alfa, t0, talfa
	printf("iw=%d\nrmax=%lf\ntaldelta=%lf\ntalmax=%lf\nc=%lf\nro=%lf\nk=%lf\nalfa=%lf\nt0=%lf\ntalfa=%lf\n", iw, rmax, taldelta, talmax, c, ro, k, alfa, t0, talfa);
	
	cout << endl;
}

void test1() {
	cout << endl;
	cout << "Macierz Globalna H:" << endl;
	for (int i = 0; i < iw; i++) {
		for (int j = 0; j < iw; j++)
			cout << htab[i][j] << " ";
		cout << endl;
	}
	cout << endl;
}

void test2() {
	cout << endl;
	cout << "Wektor P:" << endl;
	for (int i = 0; i < iw; i++)
		cout << "p[" << i + 1 << "]=" << ptab[i] << " ";
	cout << endl<<endl<<endl;
}

void test3() {
	cout << "Macierz Globalna HP:" << endl;
	for (int i = 0; i < iw; i++) {
		for (int j = 0; j < iw + 1; j++)
			cout << hptab[i][j] << " ";
		cout << endl;
	}
	cout << endl;
}

int main() {
	wczytajDaneWejsciowe();
	wypelnijPusteTablice();
	//test0();
	wypelnijHtab();
	//test1();
	test2();
	zlaczHiP();
	//test3();
	wypiszTempiZapisz();
	system("pause");
}

/*
jak ma wyglad format wyjsciowy?
popraw wypisztempizapisz()
czy wszystkie dane wejsciowe?
jak z proejktem do tego zadania?
htab[iw - 1][iw - 1] jest ~~470 a ma byc 527
*/
