// finite diference method.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <chrono>

using namespace std;

//napisane vo Visual Studio 2015 - (kompilátor: Visual C++ 14.0)

vector<double> Gauss_Jordan(const int N, vector<vector<double>> matrix, vector<double> vec)
{
	vector<double> swap, x;
	double help;
	int pivot_pos;
	double pivot;
	int i, j, k, m, n;

	swap.resize(N);
	x.resize(N);

	//gauss jordan elimination
	for (k = 0; k < N; k++)
	{
		//pivotovanie
		pivot = matrix[k][k];
		pivot_pos = k;
		for (i = k; i < N; i++)
		{
			if (abs(matrix[i][k]) > abs(pivot))
			{
				pivot = matrix[i][k];
				pivot_pos = i;
			}
		}

		//prehodenie riadkov matice
		for (j = 0; j < N; j++)
		{
			swap[j] = matrix[k][j];
			matrix[k][j] = matrix[pivot_pos][j];
			matrix[pivot_pos][j] = swap[j];
		}

		//prehodenie prvkov vektora
		help = vec[k];
		vec[k] = vec[pivot_pos];
		vec[pivot_pos] = help;

		//koniec pivotovania

		//eliminácia
		for (i = 0; i < N; i++)
		{

			if (i != k)
			{
				help = matrix[i][k] / matrix[k][k];
				for (j = 0; j < N; j++)
				{
					matrix[i][j] -= matrix[k][j] * help;
					if (abs(matrix[i][j]) < 10E-14) matrix[i][j] = 0; //predchadzanie nestabilitam vznikajucim nepresostou pri deleni

				}
				vec[i] -= vec[k] * help;
			}
		}

	}

	for (j = 0; j < N; j++)
	{
		x[j] = vec[j] / matrix[j][j];
	}

	return x;
}

vector<double> TDMA(const int N, vector<vector<double>> matrix, vector<double> vec) //pre tridiagonalne matice (Thomas algorithm)
{
	int i;
	vector<double> x;
	double help;
	x.resize(N);

	for (i = 1; i < N; i++)
	{
		help = matrix[i][i - 1] / matrix[i - 1][i - 1];
		matrix[i][i - 1] -= matrix[i - 1][i - 1] * help;
		matrix[i][i] -= matrix[i - 1][i] * help;;
		vec[i] -= vec[i-1] * help;
	}

	x[N - 1] = vec[N - 1] / matrix[N - 1][N - 1];

	for (i = N - 2; i > -1; i--)
		x[i] = (vec[i] - x[i + 1] * matrix[i][i + 1]) / matrix[i][i];

		return x;
}

int main()
{
	const double pi = 3.14159265358979323;
	vector<vector<double>> matrix, riesenie;
	vector<double> vec;
	int i, j;
	const int n = 302;
	double h;
	bool sw;

	//vo¾ba metody vypoètu
	cout << "Akou metodou chces riesit sustavu rovnic?" << endl;
	cout << "Gauss-Jordanova metoda (0)" << endl;
	cout << "TDMA - Thomas algorithm (1)" << endl;
	cin >> sw;
	cout << endl;

	//nastavenie rozmerov poli
	riesenie.resize(4);
	vec.resize(n);
	matrix.resize(n);
	for (i = 0; i < 4; i++) riesenie[i].resize(n);
	for (i = 0; i < n; i++) matrix[i].resize(n);
	
	//vypocet priestoroveho kroku h
	h = (pi - 0) / double(n - 1);

	//zapis diskretizacie prietoru
	for (i = 0; i < n; i++)
		riesenie[0][i] = i*h;

	//inicializacia merania casu, vyuziva chrono kniznicu
	auto zaciatok_central = chrono::high_resolution_clock::now();

	//metoda vyuzivajuca centralnu diferenciu - naplnenie matice a vektora pravej strany
	for (i = 1; i < n - 1; i++)
	{
		matrix[i][i - 1] = -1.0 - (h*cos(i*h)) / 2.0;
		matrix[i][i] = 2.0 + (sin(i*h) + 1)*pow(h, 2.0);
		matrix[i][i + 1] = -1.0 + (h*cos(i*h)) / 2.0;

		vec[i] = -pow(h, 2.0)*(-3 * sin(i*h) - 2);
	}

	//okrajove podmienky
	matrix[0][0] = 1;
	matrix[n - 1][n - 1] = 1;

	vec[0] = 1.0;
	vec[n-1] = 1.0;

	//volanie funkcie ktora riesi sustavu rovnic a nasledne kopirovanie do pola rieseni
	if (sw == 0) riesenie[1].swap(Gauss_Jordan(n, matrix, vec));
	else riesenie[1].swap(TDMA(n, matrix, vec));
	

	auto koniec_central = chrono::high_resolution_clock::now();
	auto zaciatok_upw_ups = chrono::high_resolution_clock::now();

	//upwind/upstream schema - naplnenie matice a vektora pravej strany
	for (i = 1; i < n - 1; i++)
	{
		if (cos(i*h) < 0)
		{
			matrix[i][i - 1] = -1.0;
			matrix[i][i] = 2.0 + (sin(i*h) + 1)*pow(h, 2.0) - h*cos(i*h);
			matrix[i][i + 1] = -1.0 + (h*cos(i*h));
		}
		if (cos(i*h) >= 0)
		{
			matrix[i][i - 1] = -1.0 - h*cos(i*h);
			matrix[i][i] = 2.0 + (sin(i*h) + 1)*pow(h, 2.0) + h*cos(i*h);
			matrix[i][i + 1] = -1.0;
		}


		vec[i] = -pow(h, 2.0)*(-3 * sin(i*h) - 2);
	}

	//okrajove podmienky
	matrix[0][0] = 1;
	matrix[n - 1][n - 1] = 1;

	vec[0] = 1.0;
	vec[n - 1] = 1.0;

	//volanie funkcie ktora riesi sustavu rovnic a nasledne kopirovanie do pola rieseni
	if (sw == 0) riesenie[2].swap(Gauss_Jordan(n, matrix, vec));
	else riesenie[2].swap(TDMA(n, matrix, vec));

	//presne riesenie y=sin(x) + 1

	for (i = 0; i < n ; i++)
	{
		riesenie[3][i] = sin(i*h) + 1;
	}

	auto koniec_upw_ups = chrono::high_resolution_clock::now();

	//vypis dlzky vypoctu jednotlivych schem
	if (sw == 0) cout << "Na vypocet bola pouzita Gauss-Jordanova metoda" << endl;
	else cout << "Na vypocet bola pouzita metoda TDMA" << endl;
	cout << "centralna schema - vypoctovy cas: " << chrono::duration_cast<chrono::nanoseconds>(koniec_central - zaciatok_central).count()*10E-10 << "s" << std::endl;
	cout << "upwind/upstream schema - vypoctovy cas: " << chrono::duration_cast<chrono::nanoseconds>(koniec_upw_ups - zaciatok_upw_ups).count()*10E-10 << "s" << std::endl;

	ofstream outfile1; //zápis do suboru csv na odovzdanie ulohy
	outfile1.open("priloha2.csv");

	//nastavenie presnosti vypisu
	outfile1.precision(11);
	outfile1.setf(std::ios::fixed, std::ios::floatfield);

	for (i = 0; i < 4; i++)
	{
		for (j = 0; j < n; j++)
		{
			if (j<n-1) outfile1 << riesenie[i][j] << ",";
			else outfile1 << riesenie[i][j];
		}
		outfile1 << "\n";
	}
	outfile1.close();

	system("PAUSE");

    return 0;
}

