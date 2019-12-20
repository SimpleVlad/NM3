#include "NumberSolution.h"
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
using namespace std;

//������� ������ ������� ��� ��������� �������� (2-������ ������)
//����� ��� ���������:
//U"xx + U"yy = -f(x, y)
//a <= x <= b, c <= y <= d
//��������� �������:
//U(a, y) = M1(y), U(b, y) = M2(y)
//U(x, c) = M3(x), U(x, d) = M3(x)

int main(void)
{
	setlocale(LC_ALL, "rus");

	int		Nmax = 10000;		//������������ ����� ��������
	int		S = 0;				//������� ��������
	double	eps = 0.000000001;	//�������� ��������� ��������
	double	epsMax = 0;			//����������� ��������
	int		n = 0, m = 0;		//����������� �����
	double	**V = NULL;			//������� ������ 
	double	**F = NULL;			//f(x, y) �� ������������������ ��������� � ����� �����
	double	a, b, c, d;			//������� ������� ���������� ���������
	double start, finish;
	int Exit = 1, Show = 0;

	//��������� �������������
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
		cout << "������� ����������� �����:";
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
			cout << endl << "������� ���������� ��������: ";
			cin >> Nmax;
		}
		system("cls");
		cout << "���� ������ ���������� �������, ���������� ���������...";

		start = clock();
		V = MemoryAllocator(n + 1, m + 1);
		//��������� �������(���������� ����������, ��������� �� �������� ������)
		FillStartSolution(V, n, m, a, b, c, d);
		//���������� ������ �������
		ZeidelsMethod(V, n, m, a, b, c, d, eps, Nmax, epsMax, S);
		finish = clock();

		system("cls");
		//�������
		cout << endl << "******************[�������]*******************" << endl;
		cout << "����������� �����: (" << n << ", " << m << ")" << endl;
		cout << "X ����������� ��������� [" << a << "..." << b << "]" << endl;
		cout << "Y ����������� ��������� [" << c << "..." << d << "]" << endl;
		cout << "���� �����: h = " << (b - a) / n << ", k = " << (d - c) / m << endl;
		cout << "��������� ��������: " << eps << endl;
		cout << "����������� ��������: " << epsMax << endl;
		cout << "��������� ��������: " << S << endl;
		cout << "������� �������: " << DiscrepancyOfSolution(V, n, m, a, b, c, d) << endl;
		cout << "����������� �������: " << CheckNumSolution(V, n, m, a, b, c, d) << endl;
		cout << "�� ������ ���������: " << finish - start << " (ms)" << endl;
		cout << "**********************************************" << endl;

		if ((n + 1) * (m + 1) < 47)
		{
			cout << "��� ����, ����� ������� ��������� �������, ������� 1: ";
			cin >> Show;

			if (Show == 1)
			{
				ShowSolution(V, n, m);
			}
		}
		//cout << endl << "��� ����������� ������ ������� 1: ";
		//cin >> Exit;

		MemoryCleaner(V, n);
	//}

	cout << endl;
	return 0;
}