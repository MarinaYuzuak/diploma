#ifndef _VECTOR_H
#define _VECTOR_H

#include <vector>
using std::vector;

class Vector
{
public:
	double operator[](const int index) { return V[index]; };

	Vector(vector<double> _v)
	{
		V = _v;
	};

	Vector(int _size)
	{
		Size = _size;
	}

	Vector()
	{

	}

	void setV(vector<double> v)
	{
		V = v;
	}

	vector<double> getV()
	{
		return V;
	}

private:
	vector<double> V;
	int Size;

};

#endif