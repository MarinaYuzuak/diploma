#ifndef _FUNCTIONS_H
#define _FUNCTIONS_H

#include "Structures.h"

double f(Point3D p)
{
	return 0; // пока что. потом, наверное, switch сделать для разных тестов. это тест xyz.
}

double u(Point3D p)
{
	return p.x * p.y * p.z;
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


#endif