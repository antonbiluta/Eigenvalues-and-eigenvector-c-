#include <iostream>
#include <Windows.h>
#include "Methods.h"
using namespace std;

#define N 3
#define	TwoPi  6.28318530717958648
const double eps = 1e-14;

static double _root3(double x)
{
	double s = 1.;
	while (x < 1.)
	{
		x *= 8.;
		s *= 0.5;
	}
	while (x > 8.)
	{
		x *= 0.125;
		s *= 2.;
	}
	double r = 1.5;
	r -= 1. / 3. * (r - x / (r * r));
	r -= 1. / 3. * (r - x / (r * r));
	r -= 1. / 3. * (r - x / (r * r));
	r -= 1. / 3. * (r - x / (r * r));
	r -= 1. / 3. * (r - x / (r * r));
	r -= 1. / 3. * (r - x / (r * r));
	return r * s;
}


double root3(double x)
{
	if (x > 0) return _root3(x); else
		if (x < 0) return-_root3(-x); else
			return 0.;
}


int SolveP3(double* x, double a, double b, double c) {
	//---------------------------------------------------------------------------
	// In case 3 действительные корни: => x[0], x[1], x[2], return 3
	//         2 действиткльных: x[0], x[1],          return 2
	//         1 действительный: x[0], x[1] � i*x[2], return 1

	// solve cubic equation x^3 + a*x^2 + b*x + c = 0
	double a2 = a * a;
	double q = (a2 - 3 * b) / 9;
	double r = (a * (2 * a2 - 9 * b) + 27 * c) / 54;
	// equation x^3 + q*x + r = 0
	double r2 = r * r;
	double q3 = q * q * q;
	double A, B;
	if (r2 <= (q3 + eps)) {//<<-- FIXED!
		if (q3 == 0) q3 = 1;
		double t = r / sqrt(q3);
		if (t < -1) t = -1;
		if (t > 1) t = 1;
		t = acos(t);
		a /= 3; q = -2 * sqrt(q);
		x[0] = q * cos(t / 3) - a;
		x[1] = q * cos((t + TwoPi) / 3) - a;
		x[2] = q * cos((t - TwoPi) / 3) - a;
		return(3);
	}
	else {

		//A =-pow(fabs(r)+sqrt(r2-q3),1./3); 
		A = -root3(fabs(r) + sqrt(r2 - q3));
		if (r < 0) A = -A;
		B = (A == 0 ? 0 : B = q / A);

		a /= 3;
		x[0] = (A + B) - a;
		x[1] = -0.5 * (A + B) - a;
		x[2] = 0.5 * sqrt(3.) * (A - B);
		if (fabs(x[2]) < eps) { x[2] = x[1]; return(2); }
		return(1);
	}
}


int solve(double a, double b, double c, double x) {
	if ((-1) * x * x * x + a * x * x + b * x + c == 0) {

		if ((-3) * x * x + 2 * a * x + b == 0) {
			if ((-6) * x + 2 * a == 0) return 3;
			else return 2;
		}
		else return 1;
	}
	else return 0;
}

// a-l b c
// d e-l f
// j h g-l

// (a-l)(e-l)(g-l) + d*h*c + j*b*f - c*j*(e-l) - f*h*(a-l) - (g-l)*b*d
// a*e*g - a*l*g - l*e*g + l*l*g - a*e*l +a*l*l + l*e*l - l*l*l + d*h*c + j*b*f - c*j*e + c*j*l - f*h*a + f*h*l - g*b*d + l*b*d

// -l*l*l + (a+g+e)*l*l + (- a*g - e*g - a*e + c*j + f*h + b*d) * l + (a*e*g + d*h*c + j*b*f - c*j*e - f*h*a - g*b*d)
double* getKoef(double** matr) {
	// a-l b c
	// d e-l f
	// j h g-l

	// -l*l*l + (a+g+e)*l*l + (- a*g - e*g - a*e + c*j + f*h + b*d) * l + (a*e*g + d*h*c + j*b*f - c*j*e - f*h*a - g*b*d)

	double b = matr[0][0] + matr[1][1] + matr[2][2];
	double c = matr[0][2] * matr[2][0] + matr[1][2] * matr[2][1] + matr[0][1] * matr[1][0] - matr[0][0] * matr[2][2] - matr[1][1] * matr[2][2] - matr[0][0] * matr[1][1];
	double d = matr[0][0] * matr[1][1] * matr[2][2] + matr[1][0] * matr[2][1] * matr[0][2] + matr[2][0] * matr[0][1] * matr[1][2] - matr[0][2] * matr[2][0] * matr[1][1] - matr[1][2] * matr[2][1] * matr[0][0] - matr[2][2] * matr[0][1] * matr[1][0];

	double x[3] = { b, c, d };
	return x;
}


void input(double** matr) {
	std::cout << "Введите значения матрицы 3х3";
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			std::cout << "a[" << i << "]" << j << "] = ";
			cin >> matr[i][j];
		}
	}
}


void outMatr(double** matr, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << "\t" << matr[i][j];
		}
		cout << endl;
	}
}


double** generateMatrix(double(*test)[3]) {

	double** matr = new double* [N];
	for (int i = 0; i < N; i++) {
		matr[i] = new double[N];
	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			matr[i][j] = test[i][j];
		}
	}

	return matr;

}


double** genTempMatr(double** matr, double lam) {

	double** copy_matr = new double* [N];
	for (int i = 0; i < N; i++) {
		copy_matr[i] = new double[N];
	}

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			if (i == j) copy_matr[i][j] = matr[i][j] - lam;
			else copy_matr[i][j] = matr[i][j];
		}
	}

	return copy_matr;
}

void outEquals(double b, double c, double d) {
	SetConsoleOutputCP(1251);
	cout << "Характеристическое уравнение:\n";
	SetConsoleOutputCP(65001);
	cout << u8"-\u03BB\u00B3 + " << b << u8"\u03BB\u00B2 + " << c << u8"\u03BB + " << d << endl << endl;
}

void outKorni(double* x) {
	SetConsoleOutputCP(1251);
	cout << "Корни:\n";
	SetConsoleOutputCP(65001);
	for (int i = 0; i < 3; i++) {
		cout << u8"\u03BB[" << (i + 1) << "] = ";
		if (i != 2) cout << x[i] << ";\t";
		else cout << x[i] << endl;
	}
}





int gauss(double** a, double* b, double* x, const int n) {
	double aim, r;

	for (int m = 0; m < n - 1; m++) {
		double amm = a[m][m];
		int im = m;

		for (int i = m; i < n; i++) {
			if (fabs(a[i][m]) > fabs(amm)) {
				amm = a[i][m];
				im = i;
			}
		}

		if (im != m) {
			r = b[im];
			b[im] = b[m];
			b[m] = r;

			for (int k = m; k < n; k++) {
				r = a[im][k];
				a[im][k] = a[m][k];
				a[m][k] = r;
			}
		}

		for (int k = m; k <= n - 1; k++)
			a[m][k] = a[m][k] / amm;

		b[m] = b[m] / amm;

		for (int i = m + 1; i < n; i++) {
			aim = a[i][m];
			for (int k = m; k < n; k++)
				a[i][k] = a[i][k] - a[m][k] * aim;
			b[i] = b[i] - b[m] * aim;
		}
	}

	a[2][2] = 0;
	double del = a[0][1];
	for (int j = 1; j < 3; j++) {
		a[0][j] -= a[1][j] * del;
	}

	x[0] = -(a[0][2]);
	x[1] = -(a[1][2]);
	x[2] = 1;

	return 0;
}