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

void SetCoord(int p, int s, int r, Element& e)
{
	e.nodes[0].p = p; e.nodes[0].s = s; e.nodes[0].r = r;
	e.nodes[1].p = p + 1; e.nodes[1].s = s; e.nodes[1].r = r;
	e.nodes[2].p = p; e.nodes[2].s = s + 1; e.nodes[2].r = r;
	e.nodes[3].p = p + 1; e.nodes[3].s = s + 1; e.nodes[3].r = r;
	e.nodes[4].p = p; e.nodes[4].s = s; e.nodes[4].r = r + 1;
	e.nodes[5].p = p + 1; e.nodes[5].s = s; e.nodes[5].r = r + 1;
	e.nodes[6].p = p; e.nodes[6].s = s + 1; e.nodes[6].r = r + 1;
	e.nodes[7].p = p + 1; e.nodes[7].s = s + 1; e.nodes[7].r = r + 1;
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

// тут пока что не все!!! функци€ написана только дл€ вычислени€ вектора правой части, где нужны билинейные б/ф
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

int isPointInsideAnomaly(Point3D p, Anomaly a)
{
	if ((p.z <= a.z1 && p.z >= a.z0)
		&& (p.y <= a.ref_plane[2].y && p.y >= a.ref_plane[0].y)
		&& (p.x <= a.ref_plane[1].x && p.x >= a.ref_plane[0].x))
		return 1;
	else
		return 0;
}

/// <summary>
/// ќпредел€ет, пересекаютс€ ли два паралеллепипеда.
/// </summary>
/// <param name="P1"></param>
/// ¬осемь координат точек первого паралеллепипеда в соответствии с локальной нумерацией.
/// <param name="P2"></param>
/// ¬осемь координат точек второго паралеллепипеда в соответствии с локальной нумерацией.
/// <returns></returns>
bool doParallelepipedsIntersect(vector<Point3D> P1, vector<Point3D> P2)
{
	auto left1 = P1[0].x; auto left2 = P2[0].x;
	auto right1 = P1[1].x; auto right2 = P2[1].x;

	auto front1 = P1[0].y; auto front2 = P2[0].y;
	auto back1 = P1[2].y; auto back2 = P2[2].y;

	auto upper1 = P1[4].z; auto upper2 = P2[4].z;
	auto lower1 = P1[0].z; auto lower2 = P2[0].z;

	if (left1 <= right2 && right1 >= left2 && front1 <= back2 && back1 >= front2 && upper1 >= lower2 && lower1 <= upper2)
		return true;
	else return false;

}

// пока что можно считать, что эти две функции вообще не работают))))) потому что координаты плюсовые
bool isElementIntersectLayer(double zLayer, double z0El, double z1El)
{
	if (z0El < zLayer && z1El > zLayer)
		return true;
	else
		return false;
}

bool isElementBelongsLayer(Layer l, double z0El, double z1El)
{
	if (z0El >= l.z0 && z1El <= l.z1)
		return true;
	else
		return false;
}

#endif