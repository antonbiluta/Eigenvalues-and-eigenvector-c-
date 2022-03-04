#include <iostream>
#include <math.h>
#include "Methods.h"
#include <Windows.h>
using namespace std;
#define N 3

//double test[3][3] = { {1,3,4},{13,15,19},{23,41,19} };
//double test[3][3] = { {5,1,2},{1,4, 1},{2,1,3} };
double test[3][3] = { {1,-3,3},{-2,-6, 13},{-1,-4,8} };
//double test[3][3] = { {4,6,-15},{1,3, -5},{1,2,-4} };

int main() {

	// Инициализируем нужные нам параметры
	double** matr = generateMatrix(test); // матрица (в методе тестовые варианты)
	double x[3]; // массив корней

	double b = getKoef(matr)[0], c = getKoef(matr)[1], d = getKoef(matr)[2]; // получаем коэф-ты характеристического уравнения
	int res = SolveP3(x, (-1) * b, (-1) * c, (-1) * d); // находим корни
	cout << endl;
	cout << "--------------------------------------------";
	cout << endl << endl;
	outMatr(matr, N);
	cout << endl;
	cout << "--------------------------------------------";
	cout << endl << endl;

	outEquals(b, c, d); // Вывод уравнения

	outKorni(x); // Вывод корней

	cout << endl;
	cout << "--------------------------------------------";
	cout << endl;

	SetConsoleOutputCP(1251);

	for (int k = 0; k < 3; k++) {

		double** copy_matr = genTempMatr(matr, x[k]);

		cout << endl << endl << "Матрица при "; SetConsoleOutputCP(65001);
		cout << u8"\u03BB = " << x[k] << ":\n"; SetConsoleOutputCP(1251);
		cout << endl;
		outMatr(copy_matr, N);
		cout << endl;

		double* v;
		double y[3] = { 0,0,0 };
		double va[3];
		gauss(copy_matr, y, va, 3);

		cout << "v = ( " << va[0] << "; " << va[1] << "; " << va[2] << " )";
		cout << endl;
	}

}