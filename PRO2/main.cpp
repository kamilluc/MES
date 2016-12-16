#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
using namespace std;

//stale z pliku
int iw;
double rmax, rdelta, taldelta, c, ro, k, alfa, t0, talfa;
fstream plik;

void wczytajDaneWejsciowe() {
	plik.open("dane.txt", ios::in);
	if (plik.good()) {
		plik >> iw >> rmax >> rdelta >> taldelta >> c >> ro >> k >> alfa >> t0 >> talfa;
		plik.close();
		cout << "Dane wejsciowe zostaly poprawnie odczytane\n";
	}
	else
		cout << "Blad odczytu pliku 'dane.txt'\n";
}

double *ttab;
double *ptab;
double **htab;
double **hptab;

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
		cout << "t[" << i + 1 << "]= " << tempe[i] << endl;

	ofstream plik2;
	plik2.open("wyniki.txt");
	if (plik2.good()) {
		for (int i = 0; i < iw; i++)
			plik2 << "t[" << i + 1 << "]= " << tempe[i] << endl;
		plik2.close();
	}
	else
		cout << "\nBlad zapisu pliku\n\n";
}

int main() {
	wczytajDaneWejsciowe();

	system("pause");
}

/*
jak ma wyglad format wyjsciowy popraw wypisztempizapisz()
czy wszystkei dane wejsciowe?
*/
