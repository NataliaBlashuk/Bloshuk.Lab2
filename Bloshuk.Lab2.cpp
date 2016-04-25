#include <iostream>
#include <math.h>

double** AllocMemory(int);
void FreeMemory(double**, int);
void DisplayMatrix(double**, int);
double Exponent(double, double);
double MyFactorial(int);
double MyPower(long double, int);
void MatrixWithOurValues(double**, int, double);
void EllementsOfMatrix(double**, int);
void MatrixDifference(double**, double**, double**, int);
double MaxElementOfMatrixDifference(double**, int);

using namespace std;

int main()
{
	while (true)
	{
		system("cls");
		int n;
		double eps;
		while (true)
		{
			cout << "Enter the size of matrix n > 1 " << endl;
			cin >> n;
			cout << "Please, enter epsilon " << endl;
			cin >> eps;
			system("cls");
			if ((eps > 0 && eps < 1) || (n > 1)) break;
			cout << "Error! Please, enter the other epsilon " << endl;
		}
		double** a = AllocMemory(n);
		cout << "Matrix with our values" << endl;
		MatrixWithOurValues(a, n, eps);
		DisplayMatrix(a, n);
		double** b = AllocMemory(n);
		EllementsOfMatrix(b, n);
		cout << endl << "Matrix with values of standart form " << endl;
		DisplayMatrix(b, n);
		double** c = AllocMemory(n);
		MatrixDifference(a, b, c, n);
		cout << endl << "Matrix of difference" << endl;
		DisplayMatrix(c, n);
		double maxElement = MaxElementOfMatrixDifference(c, n);
		cout << "Epsilon is " << eps << endl;
		cout << "Max element in matrix of difference is " << endl << maxElement << endl;
		FreeMemory(a, n);
		FreeMemory(b, n);
		FreeMemory(c, n);

		char yes;
		cout << "If you would like to continue, please, press y or Y " << endl;
		cin >> yes;
		if (yes == 'y' || yes == 'Y')
			continue;
		break;
	}
	return 0;
}

double** AllocMemory(int n)
{
	double** a = new double*[n];
	for (int i = 0; i < n; i++)
		a[i] = new double[n];
	return a;
}
void FreeMemory(double** a, int n)
{
	for (int i = 0; i < n; i++)
		delete[]a[i];
	delete[]a;
}

void DisplayMatrix(double** a, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout.width(8);
			cout << a[i][j];
		}
		cout << endl;
	}
}

double Exponent(double x, double eps)
{
	double sum = 0.0, p = 1;
	if (x == 0)
		return 1;
	for (int i = 0; fabs(p) > eps || p == 0; i++)
	{
		p = MyPower(x, i) / MyFactorial(i);
		sum += p;
	}
	return sum;
}

double MyPower(long double x, int y)
{
	long double result = 1;
	if (y == 0)
	{
		return 1;
	}
	if (y > 0)
	{
		for (int i = 0; i < y; i++)
		{
			result = result * x;
		}
	}
	if (y < 0)
	{
		for (int i = 0; i > y; i--)
			result = result / x;
	}
	return result;
}

double MyFactorial(int x)
{
	long double res = 1;
	if (x < 0)
		return 1;
	if (x == 0)
		return 1;
	else
	{
		for (int i = 1; i <= x; i++)
		{
			res *= i;
		}
	}
	return res;
}


double Sinus(double x, double eps)
{
	long double n = 1, sum = 0.0;
	for (int i = 1; fabs(n) > eps; i++)
	{
		n = MyPower(-1, i + 1) * MyPower(x, 2 * i - 1) / MyFactorial(2 * i - 1);
		sum += n;
	}
	return sum;
}

void MatrixWithOurValues(double** a, int n, double eps)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j)
			{
				a[i][j] = ((2 * i - 1) + Exponent((2 * i), eps)) / (Sinus((4 * i), eps) + 1);
			}
			else
				a[i][j] = i - j;
		}
	}
}

void EllementsOfMatrix(double** b, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j)
			{
				double pi = 3.1415926535, k = 2 * i;
				while (k > pi)
					k -= 2 * pi;
				while (k < -pi)
					k += 2 * pi;
				b[i][j] = ((2 * i - 1) + exp(2 * i)) / (sin(4 * i) + 1);
			}
			else
				b[i][j] = i - j;
		}
	}
}

void MatrixDifference(double** a, double** b, double** c, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			c[i][j] = fabs(a[i][j] - b[i][j]);
		}
	}
}

double MaxElementOfMatrixDifference(double** c, int n)
{
	double max = c[0][0];
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (c[i][j] > max) max = c[i][j];
		}
	}
	return max;
}
