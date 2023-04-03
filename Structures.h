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
	int pA = 0, rA = 0, sA = 0; // это нужно для простой сборки глобального b. УБРАТЬ ЭТО ОТСЮДА, СМ. ТЕТРАДЬ!!!
	int pB = 0, rB = 0, sB = 0; // это нужно для простой сборки глобального b. УБРАТЬ ЭТО ОТСЮДА, СМ. ТЕТРАДЬ!!!
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

struct Element
{
	vector<int> globalNumbers;
	double sigma;

	Element()
	{
		globalNumbers.resize(8);
		sigma = 1;
	}
};

enum class Bound
{
	FRONT, BACK, 
	LEFT, RIGHT,
	LOWER, UPPER
};

#endif