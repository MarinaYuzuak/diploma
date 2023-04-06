#ifndef _LIGHTWEIGHT_H
#define _LIGHTWEIGHT_H

#include <cmath>
#include "MatrixSparse.h"

bool compareByNum(Node& a, Node& b)
{
	return a.num < b.num;
}

bool isUnique(Node& a, Node& b)
{
	return a.num == b.num;
}

int getFirstIndex(int i)
{
	return i % 2;
}

int getSecondIndex(int i)
{
	return (int)floor(i / 2) % 2;
}

int getThirdIndex(int i)
{
	return (int)floor(i / 4);
}

Element GetNumbers(int p, int s, int r, int nx, int ny)
{
	Element e;
	e.nodes[0].num = r * nx * ny + s * nx + p;
	e.nodes[1].num = e.nodes[0].num + 1;
	e.nodes[2].num = r * nx * ny + (s + 1) * nx + p;
	e.nodes[3].num = e.nodes[2].num + 1;
	e.nodes[4].num = (r + 1) * nx * ny + s * nx + p;
	e.nodes[5].num = e.nodes[4].num + 1;
	e.nodes[6].num = (r + 1) * nx * ny + (s + 1) * nx + p;
	e.nodes[7].num = e.nodes[6].num + 1;

	return e;
}

void SetCoord(vector<double>& xGrid, vector<double>& yGrid, vector<double>& zGrid, int p, int s, int r, Element& e)
{
	e.nodes[0].p = { xGrid[p], yGrid[s], zGrid[r] };
	e.nodes[1].p = { xGrid[p + 1], yGrid[s], zGrid[r] };
	e.nodes[2].p = { xGrid[p], yGrid[s + 1], zGrid[r] };
	e.nodes[3].p = { xGrid[p + 1], yGrid[s + 1], zGrid[r] };
	e.nodes[4].p = { xGrid[p], yGrid[s], zGrid[r + 1] };
	e.nodes[5].p = { xGrid[p + 1], yGrid[s], zGrid[r + 1] };
	e.nodes[6].p = { xGrid[p], yGrid[s + 1], zGrid[r + 1] };
	e.nodes[7].p = { xGrid[p + 1], yGrid[s + 1], zGrid[r + 1] };
 }

bool isPointInside(double x_min, double x_max, double y_min, double y_max, Point2D p)
{
	if (p.x >= x_min && p.x <= x_max && p.y >= y_min && p.y <= y_max)
		return true;
	else
		return false;
}

double X_func(int num, double x, double xp,  double hx)
{
	switch (num)
	{
	case 0:
		return (xp - x) / hx;
	case 1:
		return (x - xp) / hx;

	default:
		break;
	}
}

double Y_func(int num, double y, double ys, double hy)
{
	switch (num)
	{
	case 0:
		return (ys - y) / hy;
	case 1:
		return (y - ys) / hy;

	default:
		break;
	}
}

// тут пока что не все!!! функция написана только для вычисления вектора правой части, где нужны билинейные б/ф
double basisBilinearFunc(int num, double x, double xp, double y, double ys, double hx, double hy)
{
	switch (num)
	{
	case 0:
		return X_func(0, x, xp, hx) * Y_func(0, y, ys, hy);
	case 1:
		return X_func(1, x, xp, hx) * Y_func(0, y, ys, hy);
	case 2:
		return X_func(0, x, xp, hx) * Y_func(1, y, ys, hy);
	case 3:
		return X_func(1, x, xp, hx) * Y_func(1, y, ys, hy);

	default:
		break;
	}
}

int binarySearch(const vector<int>& values, int value, int left, int right)
{
	while (left != right)
	{
		auto mid = (left + right) / 2 + 1;
		if (values[mid] > value)
			right = mid - 1;
		else
			left = mid;
	}
	if (values[left] == value)
		return static_cast<int>(left);
	return -1;
}

int binarySearch(const vector<Node>& values, int value, int left, int right)
{
	while (left != right)
	{
		auto mid = (left + right) / 2 + 1;
		if (values[mid].num > value)
			right = mid - 1;
		else
			left = mid;
	}
	if (values[left].num == value)
		return static_cast<int>(left);
	return -1;
}

void removeDuplicate(vector<Node>& vec)
{
	std::sort(vec.begin(), vec.end(), compareByNum);
	auto last = std::unique(vec.begin(), vec.end(), isUnique);
	vec.erase(last, vec.end());
}

template <typename T>
void removeDuplicate(vector<T>& vec)
{
	std::sort(vec.begin(), vec.end());
	auto last = std::unique(vec.begin(), vec.end());
	vec.erase(last, vec.end());
}

#endif