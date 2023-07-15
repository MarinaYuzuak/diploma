#ifndef _FUNCTIONS_H
#define _FUNCTIONS_H

#include "Structures.h"

double f(Point3D p)
{
	return -2 * p.y * p.z;
}

double u(Point3D p)
{
	return p.x * p.x * p.y * p.z + p.x * p.y * p.z;
}

double theta(Point3D p, Bound facet)
{
	switch (facet)
	{
	case Bound::UPPER:
		return p.x * p.x * p.y + p.x * p.y;
	case Bound::FRONT:
		return -p.x * p.x * p.z - p.x * p.z;
	case Bound::BACK:
		return p.x  * p.x* p.z + p.x * p.z;
	case Bound::RIGHT:
		return 2 * p.x * p.y * p.z + p.y * p.z;
	case Bound::LOWER:
		return -p.x * p.x * p.y - p.x * p.y;
	default:
		break;
	}
}

double dotProduct(const std::vector<double> a, const std::vector<double> b)
{
	double res = 0.0;
	for (uint32_t i = 0; i < a.size(); i++)
	{
		res += a[i] * b[i];
	}
	return res;
}

vector<double> MultVectorByNum(vector<double> V, double num)
{
	vector<double> result = V;
	for (auto i = 0; i < V.size(); i++)
		result[i] *= num;

	return result;
}

void holessPrec(vector<int>& ig, vector<int>& jg, vector<double>& gg, vector<double>&
	di, vector<double>& gg_L, vector<double>& di_L, int n)
{
	di_L[0] = sqrt(di[0]);
	for (int i = 1; i < n; i++)
	{
		di_L[i] = di[i];
		for (int j = ig[i]; j < ig[i + 1]; j++)
		{
			int stolb = jg[j];
			int ind1 = ig[stolb];
			int ind2 = ig[i];
			int* ptr_j1 = &jg[ind1];
			int* ptr_j2 = &jg[ind2];
			int start = 0;
			int end = ig[stolb + 1] - ind1;
			gg_L[j] = gg[j];
			while (*ptr_j2 < stolb && start < end)
			{
				if (*ptr_j2 == *ptr_j1)
				{
					gg_L[j] -= gg_L[ind1] * gg_L[ind2];
					ind1++; ind2++; ptr_j1++; ptr_j2++; start++;
				}
				else
				{
					if (*ptr_j1 < *ptr_j2)
					{
						ptr_j1++;
						ind1++;
						start++;
					}
					else if (*ptr_j2 < *ptr_j1)
					{
						ptr_j2++;
						ind2++;
					}
				}
			}
			gg_L[j] /= di_L[stolb];
			di_L[i] -= gg_L[j] * gg_L[j];
		}
		di_L[i] = sqrt(di_L[i]);
	}
}

void matrixMult(vector<int>& ig, vector<int>& jg, vector<double>& ggl, vector<double>&
	ggu, vector<double>& di, vector<double>& x, vector<double>& res, int n)
{
	for (int i = 0; i < n; i++)
	{
		res[i] = di[i] * x[i];
	}
	for (int i = 0; i < n; i++)
	{
		for (int j = ig[i]; j <= ig[i + 1] - 1; j++)
		{
			res[i] += ggl[j] * x[jg[j]];
			res[jg[j]] += ggu[j] * x[i];
		}
	}
}

void matrixSub(vector<double>& m1, vector<double>& m2, vector<double>& res)
{
	for (size_t i = 0; i < m1.size(); i++)
	{
		res[i] = m1[i] - m2[i];
	}
}
void matrixMultByNum(double num, vector<double>& m, vector<double>& res)
{
	for (size_t i = 0; i < m.size(); i++)
	{
		res[i] = num * m[i];
	}
}
void slauHolessStright(vector<int>& ig, vector<int>& jg, vector<double>& gg,
	vector<double>& di, vector<double>& f, vector<double>& res, int n)
{
	for (int i = 0; i < n; i++)
	{
		res[i] = f[i];
		for (int j = ig[i]; j < ig[i + 1]; j++)
		{
			res[i] -= gg[j] * res[jg[j]];
		}
		res[i] /= di[i];
	}
}
void slauHolessReversal(vector<int>& ig, vector<int>& jg, vector<double>& gg,
	vector<double>& di, vector<double>& f, vector<double>& res, int n)
{
	res = f;
	for (int j = n - 1; j >= 0; j--)
	{
		res[j] /= di[j];
		for (int i = ig[j]; i < ig[j + 1]; i++)
		{
			res[jg[i]] -= gg[i] * res[j];
		}
	}
}


void slauHoless(vector<int>& ig, vector<int>& jg, vector<double>& gg, vector<double>&
	di, vector<double>& f, vector<double>& res, int n)
{
	vector<double> y(n);
	slauHolessStright(ig, jg, gg, di, f, y, n);
	slauHolessReversal(ig, jg, gg, di, y, res, n);
}

void matrixAdd(vector<double>& m1, vector<double>& m2, vector<double>& res)
{
	for (size_t i = 0; i < m1.size(); i++)
	{
		res[i] = m1[i] + m2[i];
	}
}


void cgmWithHolesPrec(vector<int>& ig, vector<int>& jg, vector<double>& gg,
	vector<double>& di, vector<double>& f, vector<double>& res, int n, int maxiter, double
	eps)
{
	for (auto i = 0; i < res.size(); i++)
	{
		res[i] = 0;
	}
	
	vector<double> gg_L(gg.size());
	vector<double> di_L(n);
	holessPrec(ig, jg, gg, di, gg_L, di_L, n);
	vector<double> xk1(n, 0);
	vector<double> rk(n), rk1(n), zk(n), Am(n), numA(n);

	double scalMult, f_norm = sqrt(dotProduct(f, f));
	vector<double> Azk1(n), M1rk(n);
	double alpha, beta, residual;
	int k = 1;

	matrixMult(ig, jg, gg, gg, di, xk1, Am, n);
	matrixSub(f, Am, rk1);
	slauHoless(ig, jg, gg_L, di_L, rk1, zk, n);

	while (true)
	{
		slauHoless(ig, jg, gg_L, di_L, rk1, Am, n);
		scalMult = dotProduct(Am, rk1); // ( M^-1*r(k-1), r(k-1) )
		matrixMult(ig, jg, gg, gg, di, zk, Azk1, n);
		alpha = scalMult / dotProduct(Azk1, zk);
		matrixMultByNum(alpha, zk, numA);
		matrixAdd(xk1, numA, res);
		matrixMultByNum(alpha, Azk1, numA);
		matrixSub(rk1, numA, rk);
		slauHoless(ig, jg, gg_L, di_L, rk, M1rk, n);
		beta = dotProduct(M1rk, rk) / scalMult;
		matrixMultByNum(beta, zk, numA);
		matrixAdd(M1rk, numA, zk);
		xk1 = res;
		rk1 = rk;
		k++;
		residual = sqrt(dotProduct(rk, rk)) / f_norm;

		if (k > maxiter)
		{
			std::cout << "The maximum number of iterations has been exceeded." <<
				std::endl;
			break;
		}
		if (residual < eps)
		{
			
			std::cout << "Residual exit." << std::endl;
			break;
		}
	}
}

#endif