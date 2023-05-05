#ifndef _VECTOR_H
#define _VECTOR_H

#include <vector>
using std::vector;

class Vector
{
public:
	double operator[](const int index) { return V[index]; };
	vector<double> MultVectorByNum(double num);

	Vector(vector<double> _v)
	{
		V = _v;
	};

	Vector(int _size)
	{
		Size = _size;
	}

	Vector(int size, int num)
	{
		V.resize(size);
		for (auto i = 0; i < size; i++)
			V[i] = num;
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

//vector<double> Vector::MultVectorByNum(double num)
//{
//	vector<double> result = V;
//	for (auto i = 0; i < V.size(); i++)
//		result[i] *= num;
//
//	return result;
//}

#endif