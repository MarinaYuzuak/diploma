#ifndef _STRUCTURES_H
#define _STRUCTURES_H

#include <iostream>
#include <vector>
using std::vector;

struct Point3D
{
	Point3D() :
		x(0.0), y(0.0), z(0.0)
	{}

	Point3D(double _x, double _y, double _z) :
		x(_x), y(_y), z(_z)
	{}

	double x, y, z;
};

struct Point2D
{
	Point2D() :
		x(0.0), y(0.0)
	{}

	Point2D(double _x, double _y) :
		x(_x), y(_y)
	{}

	double x, y;
};

struct Layer 
{
	Layer() :
		sigma(0.0), z0(0.0), z1(0.0)
	{}

	Layer(double _sigma, double _z0, double _z1):
		sigma(_sigma), z0(_z0), z1(_z1)
	{}

	double sigma;
	double z0, z1;     // Layer thickness.
};

struct Receiver
{
	Receiver() : 
		M({ 0, 0 }),
		N({ 0, 0 })
	{}

	Receiver(Point2D _M, Point2D _N) :
		M(_M), N(_N)
	{}

	Point2D M;
	Point2D N;
};

struct Source
{
	Source() :
		A({ 0, 0 }),
		B({ 0, 0 })
	{}

	Source(Point2D _A, Point2D _B) :
		A(_A), B(_B)
	{}


	Point2D A;
	Point2D B;
	int elA = 0; // номер конечного элемента, которому принадлежит точка А (нужно для вычисления глобального вектора b)
	int elB = 0; // номер конечного элемента, которому принадлежит точка B (нужно для вычисления глобального вектора b)
};

struct Anomaly
{
	Anomaly()
	{
		layer_number = 1;
		ref_plane.resize(4);
		sigma = 0.0; z0 = 0.0; z1 = 0.0;
	}

	std::size_t layer_number;
	vector<Point2D> ref_plane;     // Reference plane.
	double z0, z1;
	double sigma;
};

enum class Bound
{
	FRONT, BACK,
	LEFT, RIGHT,
	LOWER, UPPER,
	UNDEFINED
};

struct Node
{
	int num; // глобальный номер узла
	Bound bound;
	int p, s, r;

	Node() :
		num(0),
		p(0), s(0), r(0),
		bound(Bound::UNDEFINED)
	{}

	Node(int _num, int _p, int _s, int _r, Bound _bound) :
		num(_num),
		p(_p), s(_s), r(_r),
		bound(_bound)
	{}
};

struct Element
{
	vector<Node> nodes;
	double sigma;

	Element()
	{
		nodes.resize(8);
		sigma = 1;
	}
};

struct GridSettins
{
	int numInt;
	double koef;

	double x1;
	double x2;

	GridSettins()
	{
		x1 = 0;
		x2 = 0;
		numInt = 0;
		koef = 1;
	}

	GridSettins(double x1, double x2, int numInt, double koef)
	{
		this->x1 = x1;
		this->x2 = x2;
		this->numInt = numInt;
		this->koef = koef;
	}
};

#endif