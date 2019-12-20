#include "NumberSolution.h"
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
using namespace std;

//Решение задачи Дирихле для уравнения Пуассона (2-мерный случай)
//Общий вид уравнения:
//U"xx + U"yy = -f(x, y)
//a <= x <= b, c <= y <= d
//Граничные условия:
//U(a, y) = M1(y), U(b, y) = M2(y)
//U(x, c) = M3(x), U(x, d) = M3(x)

int main(void)
{
	setlocale(LC_ALL, "rus");

	int		Nmax = 10000;		//Максимальное число итераций
	int		S = 0;				//Счетчик итераций
	double	eps = 0.000000001;	//Параметр требуемой точности
	double	epsMax = 0;			//Достигнутая точность
	int		n = 0, m = 0;		//Размерность сетки
	double	**V = NULL;			//Искомый вектор 
	double	**F = NULL;			//f(x, y) из дифференециального уравнения в узлах сетки
	double	a, b, c, d;			//Границы области определния уравнения
	double start, finish;
	int Exit = 1, Show = 0;

	//Начальная инициализация
	a = 0;
	b = 1;
	c = 0;
	d = 3;

	//while (Exit == 1)
	//{
		n = 0; 
		m = 0;
		Nmax = 0;
		
		system("cls");
		cout << "Введите размерность сетки:";
		while (n <= 0)
		{
			cout << endl << "n = ";
			cin >> n;
		}
		while (m <= 0)
		{
			cout << endl << "m = ";
			cin >> m;
		}
		while (Nmax <= 0)
		{
			cout << endl << "Введите количество итераций: ";
			cin >> Nmax;
		}
		system("cls");
		cout << "Идет расчет численного решения, пожалуйста подождите...";

		start = clock();
		V = MemoryAllocator(n + 1, m + 1);
		//Начальное решение(Заполнение правильное, проверено на тестовой задаче)
		FillStartSolution(V, n, m, a, b, c, d);
		//Применение метода Зейдаля
		ZeidelsMethod(V, n, m, a, b, c, d, eps, Nmax, epsMax, S);
		finish = clock();

		system("cls");
		//Справка
		cout << endl << "******************[СПРАВКА]*******************" << endl;
		cout << "Размерность сетки: (" << n << ", " << m << ")" << endl;
		cout << "X принадлежит интервалу [" << a << "..." << b << "]" << endl;
		cout << "Y принадлежит интервалу [" << c << "..." << d << "]" << endl;
		cout << "Шаги сетки: h = " << (b - a) / n << ", k = " << (d - c) / m << endl;
		cout << "Требуемая точность: " << eps << endl;
		cout << "Достигнутая точность: " << epsMax << endl;
		cout << "Проведено итераций: " << S << endl;
		cout << "Невязка решения: " << DiscrepancyOfSolution(V, n, m, a, b, c, d) << endl;
		cout << "Погрешность решения: " << CheckNumSolution(V, n, m, a, b, c, d) << endl;
		cout << "На расчет затрачено: " << finish - start << " (ms)" << endl;
		cout << "**********************************************" << endl;

		if ((n + 1) * (m + 1) < 47)
		{
			cout << "Для того, чтобы увидеть численное решение, нажмите 1: ";
			cin >> Show;

			if (Show == 1)
			{
				ShowSolution(V, n, m);
			}
		}
		//cout << endl << "Для продолжения работы нажмите 1: ";
		//cin >> Exit;

		MemoryCleaner(V, n);
	//}

	cout << endl;
	return 0;
}