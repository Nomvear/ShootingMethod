// Shooting.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>


using namespace std;

double F(double x, double v)
{
	return exp(x)*(x*x+x-3) + v;
}

double FF(double x)
{
	return ((2 * pow(x, 3) - 3 * x * x - 12 * x + (2 * (7 * exp(1) + 6))/(2*exp(1) - 1) + 12) * exp(x)) / 6 - (19 * exp(1)) / (3 * (2 * exp(1) - 1));
}

double Runge_Kutta_4(double p)
{
	double n = 10000;
	double k0, k1, k2, k3;
	double h = 1. / n;
	double x_0, v_0, u_n, v_n, u_0;
	
	x_0 = 0;
	u_n = u_0 = 0;
	v_n = v_0 = p;

	k0 = v_0 + h * F(x_0, v_0);

	for (int i = 0; i <= n; i++)
	{
		x_0 = i * h;
		u_0 = u_n;
		v_0 = k0;

		k0 = v_0 + h * F(x_0, v_0);
		k1 = v_0 + h * F(x_0 + h / 2., v_0 + h * k0 / 2.);
		k2 = v_0 + h * F(x_0 + h / 2., v_0 + h * k1 / 2.);
		k3 = v_0 + h * F(x_0 + h, v_0 + h * k2);
		u_n = u_0 + (h / 6) * (k0 + 2 * k1 + 2 * k2 + k3);
		v_n = k0;
	}

	return (u_n + v_n);
	
}
double Find_Param()
{
	double p1 = 10;
	double p2 = -1;
	double p = (p2 - p1) / 2.;
	while (Runge_Kutta_4(p1) > 1e-10)
	{
		cout << Runge_Kutta_4(p1) << endl;
		p = (p2 + p1) / 2.;
		cout << p << endl;
		if (Runge_Kutta_4(p1) * Runge_Kutta_4(p) < 0)
			p2 = p;
		else if (Runge_Kutta_4(p) * Runge_Kutta_4(p2) < 0)
			p1 = p;
	}
	return p1;
}

int main()
{
	ofstream outsh;
	outsh.open("Shooting.txt");
	ofstream doutsh;
	doutsh.open("DiffShooting.txt");
	double p = Find_Param();
	cout << p;
	double k0, k1, k2, k3;
	double n = 1e4 +1 ;
	double h = 1. / n;
	double x_0, v_0, u_n, v_n, u_0;
	x_0 = 0; 
	u_n = u_0 = 0;
	v_n = v_0 = p;

	k0 = v_0 + h * F(x_0, v_0);

	for (int i = 0; i <= n; i++)
	{
		outsh << i * h << " " << u_n << endl;
		doutsh << i * h << " " << u_n - FF(i * h) << endl;
		cout << u_n - FF(i*h) << endl;

		x_0 = i * h;
		u_0 = u_n;
		v_0 = k0;

		k0 = v_0 + h * F(x_0, v_0);
		k1 = v_0 + h * F(x_0 + h / 2., v_0 + h * k0 / 2.);
		k2 = v_0 + h * F(x_0 + h / 2., v_0 + h * k1 / 2.);
		k3 = v_0 + h * F(x_0 + h, v_0 + h * k2);
		u_n = u_0 + (h / 6) * (k0 + 2 * k1 + 2 * k2 + k3);
		v_n = k0;
	}
	outsh.close();
	doutsh.close();
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
