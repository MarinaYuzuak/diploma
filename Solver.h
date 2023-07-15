#ifndef _SOLVER_H
#define _SOLVER_H

#include <numbers>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include "functions.h"
#include "Matrix.h"
#include "MatrixSparse.h"

using std::filesystem::path;

#define _SIMPLE_TESTS
//#define _NORMAL_TESTS

class Solver
{
public:
	Solver(const path& _path)
	{
		InitDomainData(_path);

		auto&& [ig, jg] = GeneratePortrait();
		_stiffnessMatrix.setIg(ig); _stiffnessMatrix.setJg(jg);

		SetSigmaForElements();
		DefineBoundsNodes();
	};
	

private:
	vector<Point2D> _refPlane;
	vector<Anomaly> _anomals;
	vector<Layer> _layers;
	vector<Source> _sources;
	vector<Receiver> _receivers;

	vector<Point3D> _nodes;

	vector<GridSettins> _xSettings;
	vector<GridSettins> _ySettings;
	int _xNesting;
	int _yNesting;

	MatrixSparse _stiffnessMatrix;   
	vector<double> _solution;  
	vector<double> _b;

	vector<Node> _nodesDirichlet;  
	vector<vector<Node>> _nodesNeumann;  
	vector<int> _boundCond;

	vector<double> _zGrid;  
	vector<double> _xGrid;  
	vector<double> y_Grid;  
	vector<Element> _elements;  


private:
	void InitDomainData(const path& _path);  
	void FindElementForPoint();

	void GetSolutionOnSurface();
	void CalculateSolution();  

	std::pair<vector<int>, vector<int>> GeneratePortrait();
	void AssembleMatrix();  
	void AssembleVector();  
	Matrix BuildLocalG(double hx, double hy, double hz, double sigma);  
	Matrix BuildLocalC(double hx, double hy, double hz);  

	void ApplyDirichlet();  
	void ApplyNeumann();  

	void PerformGaussianReduction();  
	void ReduceRow(int row);  

	void DefineBoundsNodes();  
	void DefineGlobalNumbers();  
	void SetSigmaForElements();  
	void DefineFirstBoundNodes(int start, int elementWithWhichLastRowBegins, vector<int>& locNums, int addition, Bound facet, int dim);  
	void DefineSecondBoundNodes(int start, int elementWithWhichLastRowBegins, vector<int>& locNums, int addition, Bound facet, int dim);  
	void BuildZGrid(const path& _path);  
	void BuildXGrid();  
	void BuildYGrid();  
	void BuildDomain(const path& _path);

	vector<Point3D> GetElemPoints(Element e);

	int GetElemNumForPoint(Point3D p);
	double GetSolutionInPoint(Point3D p);
	double GetNumericalPotential(int recieverNum);
	double GetAnalyticalPotential(int sourceNum, int recieverNum);


	double GetVForIntersectionAnomaly(Element e, Anomaly a); 
	double GetVElement(Element e);  
	Point3D GetPoint(int p, int s, int r);  
	vector<Point3D> GetAnomalyPoints(Anomaly a);  

};

double Solver::GetNumericalPotential(int recieverNum)
{
	auto Vm = GetSolutionInPoint({ _receivers[recieverNum].M.x, _receivers[recieverNum].M.y, 0 });
	auto Vn = GetSolutionInPoint({ _receivers[recieverNum].N.x, _receivers[recieverNum].N.y, 0 });

	return Vm - Vn;
}

double Solver::GetAnalyticalPotential(int sN, int rN)
{
	return (1. / (2 * std::numbers::pi * _layers.back().sigma)) *
		(
			(1. / distance(_receivers[rN].M, _sources[sN].B) - 1. / distance(_receivers[rN].M, _sources[sN].A))

			- 

			(1. / distance(_receivers[rN].N, _sources[sN].B) - 1. / distance(_receivers[rN].N, _sources[sN].A))

		);

}

double Solver::GetSolutionInPoint(Point3D p)
{
	auto elemNum = GetElemNumForPoint(p);
	auto elemPoints = GetElemPoints(_elements[elemNum]);
	double res = 0;

	auto e = _elements[elemNum];
	for (auto i = 0; i < e.nodes.size(); i++)
	{
		auto a = getFirstIndex(i);
		auto b = getSecondIndex(i);
		auto c = getThirdIndex(i);

		Point3D pEl;

		a == 0 ? pEl.x = elemPoints[1].x : pEl.x = elemPoints[0].x;
		b == 0 ? pEl.y = elemPoints[2].y : pEl.y = elemPoints[0].y;
		c == 0 ? pEl.z = elemPoints[4].z : pEl.z = elemPoints[0].z;

		Point3D h{ elemPoints[1].x - elemPoints[0].x, elemPoints[2].y - elemPoints[0].y, elemPoints[4].z - elemPoints[0].z };

		res += _solution[e.nodes[i].num] * Psi(i, p, pEl, h);
	}
	return res;
}

int Solver::GetElemNumForPoint(Point3D p)
{
	int layerNum;

	for (auto i = 0; i < _layers.size(); i++)
	{
		if (p.z >= _layers[i].z0 && p.z <= _layers[i].z1)
		{
			layerNum = i;
			break;
		}
	}

	int Nz = 0;

	while (_zGrid[Nz] != _layers[layerNum].z1)
	{
		Nz++;
	}

	auto start = (_xGrid.size() - 1) * (y_Grid.size() - 1) * (Nz - 1);

	for (auto i = start; i < _elements.size(); i++)
	{
		auto points = GetElemPoints(_elements[i]);
		if (p.x >= points[0].x && p.x <= points[1].x
			&& p.y >= points[0].y && p.y <= points[2].y)
		{
			return i;
		}
	}
}

vector<Point3D> Solver::GetElemPoints(Element e)
{
	vector<Point3D> result;
	for (auto i = 0; i < e.nodes.size(); i++)
	{
		result.push_back(GetPoint(e.nodes[i].p, e.nodes[i].s, e.nodes[i].r));
	}

	return result;
}

void Solver::GetSolutionOnSurface()
{
	int nx = _xGrid.size();
	int ny = y_Grid.size();
	int nz = _zGrid.size();

	int elemWithWithLastLayerStarts = (nx - 1) * (ny - 1) * (nz - 2);
	auto start = _elements[elemWithWithLastLayerStarts].nodes[4].num;

	int increment = 1;

	int sumIntX = 0, sumIntY = 0;

	for (auto i = 0; i < _xSettings.size(); i++)
	{
		sumIntX += _xSettings[i].numInt;
	}

	for (auto i = 0; i < _ySettings.size(); i++)
	{
		sumIntY += _ySettings[i].numInt;
	}

	int xNest = _xNesting, yNest = _yNesting;
	xNest == 0 ? xNest++ : xNest *= 2;
	yNest == 0 ? yNest++ : yNest *= 2;

	std::ofstream fout;
	fout.open("surface_graf.txt");
	fout << "X\tY\tSolution" << std::endl;

	for (auto i = 0; i <= sumIntY; i++)
	{
		auto end = start + xNest * sumIntX + 1;

		for (auto j = start; j < end; j += xNest, start += xNest)
		{
			if (fabs(_solution[j]) > 1e-15)
			{
				fout << _nodes[j].x << "\t" << _nodes[j].y << "\t";
				fout << std::setprecision(15) << _solution[j] << std::endl;
			}
		}
		start -= xNest - 1;
		start += (yNest - 1) * (sumIntX * xNest + 1);
	}

	fout.close();
}

void Solver::FindElementForPoint()
{
	int nx = _xGrid.size();
	int ny = y_Grid.size();
	int nz = _zGrid.size();

	int elemWithWithLastLayerStarts = (nx - 1) * (ny - 1) * (nz - 2);

	for (auto i = 0; i < _sources.size(); i++)
	{
		bool is_found_A = false;
		bool is_found_B = false;
		int whereToSearchPoint = elemWithWithLastLayerStarts;

		for (auto s = 0; s < ny - 1; s++)
		{
			for (auto p = 0; p < nx - 1; p++)
			{
				if (!is_found_A)
				{
					is_found_A = isPointInside(_xGrid[p], _xGrid[p + 1], y_Grid[s], y_Grid[s + 1], _sources[i].A);
					if (is_found_A)
					{
						_sources[i].elA = whereToSearchPoint;
					}
				}

				if (!is_found_B)
				{
					is_found_B = isPointInside(_xGrid[p], _xGrid[p + 1], y_Grid[s], y_Grid[s + 1], _sources[i].B);
					if (is_found_B)
					{
						_sources[i].elB = whereToSearchPoint;
					}
				}

				if (is_found_A && is_found_B)
					break;

				whereToSearchPoint++;
			}
			if (is_found_A && is_found_B)
				break;
		}
	}
}

vector<Point3D> Solver::GetAnomalyPoints(Anomaly a)
{
	vector<Point3D> points(8);
	points[0] = Point3D{ a.ref_plane[0].x, a.ref_plane[0].y, a.z0 };
	points[1] = Point3D{ a.ref_plane[1].x, a.ref_plane[1].y, a.z0 };
	points[2] = Point3D{ a.ref_plane[2].x, a.ref_plane[2].y, a.z0 };
	points[3] = Point3D{ a.ref_plane[3].x, a.ref_plane[3].y, a.z0 };

	points[4] = Point3D{ a.ref_plane[0].x, a.ref_plane[0].y, a.z1 };
	points[5] = Point3D{ a.ref_plane[1].x, a.ref_plane[1].y, a.z1 };
	points[6] = Point3D{ a.ref_plane[1].x, a.ref_plane[2].y, a.z1 };
	points[7] = Point3D{ a.ref_plane[3].x, a.ref_plane[3].y, a.z1 };

	return points;
}

double Solver::GetVForIntersectionAnomaly(Element e, Anomaly a)
{
	vector<double> x, y, z;
	vector<Point3D> points;

	for (auto i = 0; i < e.nodes.size(); i++)
		points.push_back(GetPoint(e.nodes[i].p, e.nodes[i].s, e.nodes[i].r));

	x.push_back(a.ref_plane[0].x); x.push_back(a.ref_plane[1].x); 
	x.push_back(points[0].x); x.push_back(points[1].x);
	std::sort(x.begin(), x.end());

	y.push_back(a.ref_plane[0].y); y.push_back(a.ref_plane[2].y);
	y.push_back(points[0].y); y.push_back(points[2].y);
	std::sort(y.begin(), y.end());

	z.push_back(a.z0); z.push_back(a.z1);
	z.push_back(points[0].z); z.push_back(points[4].z);
	std::sort(z.begin(), z.end());

	double h = z[2] - z[1];
	double sOsn = (x[2] - x[1]) * (y[2] - y[1]);

	return (1. / 3.) * sOsn * h;
}

double Solver::GetVElement(Element e)
{
	double sOsn = (_xGrid[e.nodes[1].p] - _xGrid[e.nodes[0].p]) * (y_Grid[e.nodes[2].s] - y_Grid[e.nodes[0].s]);
	double h = _zGrid[e.nodes[4].r] - _zGrid[e.nodes[0].r];
	return (1. / 3.) * sOsn * h;
}

void Solver::SetSigmaForElements()
{
	vector<vector<Anomaly>> anomaliesWithWhichFiniteElementIntersects(_elements.size());

	for (auto i = 0; i < _elements.size(); i++)
	{
		bool intersectWithAnomal = false;
		for (auto j = 0; j < _anomals.size(); j++)
		{
			int indicator = 0;
			vector<Point3D> pointsEl(8);
			vector<Point3D> pointsAnomaly = GetAnomalyPoints(_anomals[j]);

			for (auto l = 0; l < _elements[i].nodes.size(); l++)
			{
				auto p = GetPoint(_elements[i].nodes[l].p, _elements[i].nodes[l].s, _elements[i].nodes[l].r);
				pointsEl[l] = p;
				indicator += isPointInsideAnomaly(p, _anomals[j]);
			}

			// КЭ полностью лежит в АО.
			if (indicator == 8)
			{
				_elements[i].sigma = _anomals[j].sigma;
				intersectWithAnomal = true;
				break; 
			}

			if (doParallelepipedsIntersect(pointsEl, pointsAnomaly))
			{
				intersectWithAnomal = true;
				anomaliesWithWhichFiniteElementIntersects[i].push_back(_anomals[j]);
			}
		}

		if (!intersectWithAnomal)
		{
			/*смотрим, какому слою (или слоям) принадлежит КЭ*/
			for (auto l = 0; l < _layers.size(); l++)
			{
				bool isBelong = isElementBelongsLayer(_layers[l], _zGrid[_elements[i].nodes[0].r], _zGrid[_elements[i].nodes[4].r]);
				if (isBelong)
				{
					_elements[i].sigma = _layers[l].sigma;
					isBelong = true;
					break;
				}
			}
		}
	}

  	for (auto i = 0; i < anomaliesWithWhichFiniteElementIntersects.size(); i++)
	{
		if (anomaliesWithWhichFiniteElementIntersects[i].size() != 0) // считаем долю объема от каждого АО и от каждого слоя, содержащего АО.
		{
			/*значит, что КЭ лежит в одном слое и пересекается с каким-то количеством АО.*/
			auto layer = anomaliesWithWhichFiniteElementIntersects[i][0].layer_number;
			double vEl = GetVElement(_elements[i]);
			double vFromLayer = vEl;

			_elements[i].sigma = 0;
			for (auto l = 0; l < anomaliesWithWhichFiniteElementIntersects[i].size(); l++)
			{
				double V = GetVForIntersectionAnomaly(_elements[i], anomaliesWithWhichFiniteElementIntersects[i][l]);
				vFromLayer -= V;
				_elements[i].sigma += (V / vEl) * anomaliesWithWhichFiniteElementIntersects[i][l].sigma;
			}
			vFromLayer = vFromLayer / vEl;
			_elements[i].sigma += vFromLayer * _layers[layer].sigma;	
		}
	}
}

Point3D Solver::GetPoint(int p, int s, int r)
{
	auto x = _xGrid[p];
	auto y = y_Grid[s];
	auto z = _zGrid[r];
	Point3D point{ x, y, z };

	return point;
}

void Solver::ApplyNeumann()
{
#ifdef _SIMPLE_TESTS
	Matrix M({ {4, 2, 2, 1}, {2, 4, 1, 2}, {2, 1, 4, 2}, {1, 2, 2, 4} });
	vector<double> t(4, 0.0);

	for (auto i = 0; i < _nodesNeumann.size(); i++)
	{
		for (auto j = 0; j < _nodesNeumann[i].size(); j++)
		{
			Point3D p = {_xGrid[_nodesNeumann[i][j].p], y_Grid[_nodesNeumann[i][j].s], _zGrid[_nodesNeumann[i][j].r] };
			t[j] = theta(p, _nodesNeumann[i][j].bound);
		}

		double hp, hr;
		switch (_nodesNeumann[i][0].bound)
		{
		case Bound::UPPER:
			hp = _xGrid[_nodesNeumann[i][1].p] - _xGrid[_nodesNeumann[i][0].p];
			hr = y_Grid[_nodesNeumann[i][2].s] - y_Grid[_nodesNeumann[i][0].s];
			break;
		case Bound::FRONT:
			hp = _xGrid[_nodesNeumann[i][1].p] - _xGrid[_nodesNeumann[i][0].p];
			hr = _zGrid[_nodesNeumann[i][2].r] - _zGrid[_nodesNeumann[i][0].r];
			break;
		case Bound::BACK:
			hp = _xGrid[_nodesNeumann[i][1].p] - _xGrid[_nodesNeumann[i][0].p];
			hr = _zGrid[_nodesNeumann[i][2].r] - _zGrid[_nodesNeumann[i][0].r];
			break;
		case Bound::LEFT:
			hp = y_Grid[_nodesNeumann[i][1].s] - y_Grid[_nodesNeumann[i][0].s];
			hr = _zGrid[_nodesNeumann[i][2].r] - _zGrid[_nodesNeumann[i][0].r];
			break;
		case Bound::RIGHT:
			hp = y_Grid[_nodesNeumann[i][1].s] - y_Grid[_nodesNeumann[i][0].s];
			hr = _zGrid[_nodesNeumann[i][2].r] - _zGrid[_nodesNeumann[i][0].r];
			break;
		case Bound::LOWER:
			hp = _xGrid[_nodesNeumann[i][1].p] - _xGrid[_nodesNeumann[i][0].p];
			hr = y_Grid[_nodesNeumann[i][2].s] - y_Grid[_nodesNeumann[i][0].s];
			break;
		default:
			break;
		}

		auto bLocal = M.MultMatrixByVector(t);
		bLocal = MultVectorByNum(bLocal, (hp * hr) / 36.);

		for (auto j = 0; j < _nodesNeumann[i].size(); j++)
		{
			_b[_nodesNeumann[i][j].num] += bLocal[j];
		}
	}
#endif

#ifdef _NORMAL_TESTS
	std::cout << std::endl;
#endif
}

void Solver::CalculateSolution()
{
	_solution.resize(_b.size(), 0.0);
	auto eps = 1e-21;
	auto maxIter = static_cast<uint32_t>(1e+7);
	auto loadVectorNorm = sqrt(dotProduct(_b, _b));
	auto z = _b;
	auto r = _b;

	auto p = _stiffnessMatrix.MultiplyByVector(z);

	auto rDotR = dotProduct(r, r);
	auto residual = sqrt(rDotR) / loadVectorNorm;

	for (uint32_t k = 1; k < maxIter && residual > eps; k++)
	{
		auto pDotP = dotProduct(p, p);
		auto a = dotProduct(p, r) / pDotP;

		for (uint32_t i = 0; i < _b.size(); i++)
		{
			_solution[i] += a * z[i];
			r[i] -= a * p[i];
		}
		auto Ar = _stiffnessMatrix.MultiplyByVector(r);
		auto b = -dotProduct(p, Ar) / pDotP;

		for (auto i = 0; i < _b.size(); i++)
		{
			z[i] = r[i] + b * z[i];
			p[i] = Ar[i] + b * p[i];
		}
		rDotR = dotProduct(r, r);
		residual = sqrt(rDotR) / loadVectorNorm;
	}
}

void Solver::DefineSecondBoundNodes(int start, int elementWithWhichLastRowBegins, vector<int>& locNums, int addition, Bound facet, int dim)
{
	int endRow = start + addition * dim;

	while (start < elementWithWhichLastRowBegins + addition)
	{
		for (auto i = start; i <= endRow; i += addition)
		{
			vector<Node> nodes(locNums.size());
			for (auto j = 0; j < locNums.size(); j++)
			{
				nodes[j] = _elements[i].nodes[locNums[j]];
				nodes[j].bound = facet;
			}

			_nodesNeumann.push_back(nodes);
		}

		switch (facet)
		{
		case Bound::FRONT:
			start += (_xGrid.size() - 1) * (y_Grid.size() - 1);
			break;
		case Bound::BACK:
			start += (_xGrid.size() - 1) * (y_Grid.size() - 1);
			break;
		default:
			start = endRow + addition;
			break;
		}
		endRow = start + addition * dim;
	}
}

void Solver::DefineFirstBoundNodes(int start, int elementWithWhichLastRowBegins, vector<int>& locNums, int addition, Bound facet, int dim)
{
	int endRow = start + addition * dim;

	while(start < elementWithWhichLastRowBegins + addition)
	{
		_nodesDirichlet.push_back(_elements[start].nodes[locNums[0]]);
		_nodesDirichlet.back().bound = facet;

		for (auto i = start; i <= endRow; i += addition)
		{
			_nodesDirichlet.push_back(_elements[i].nodes[locNums[1]]);
			_nodesDirichlet.back().bound = facet;
		}

		switch (facet)
		{
		case Bound::FRONT:
			start += (_xGrid.size() - 1) * (y_Grid.size() - 1);
			break;
		case Bound::BACK:
			start += (_xGrid.size() - 1) * (y_Grid.size() - 1);
			break;
		default:
			start = endRow + addition;
			break;
		}
		endRow = start + addition * dim;
	}

	start = elementWithWhichLastRowBegins;
	endRow = start + addition * dim;

	_nodesDirichlet.push_back(_elements[start].nodes[locNums[2]]);
	_nodesDirichlet.back().bound = facet;
	for (auto i = start; i <= endRow; i += addition)
	{
		_nodesDirichlet.push_back(_elements[i].nodes[locNums[3]]);
		_nodesDirichlet.back().bound = facet;
	}
}

void Solver::DefineBoundsNodes()
{
	Bound facet = Bound::FRONT;
	int start = 0;
	int elementWithWhichLastRowBegins = (_xGrid.size() - 1) * (y_Grid.size() - 1) * (_zGrid.size() - 2);
	vector<int> locNums = { 0, 1, 4, 5 };
	int addition = 1;
	int dim = _xGrid.size() - 2;
	int condition = _boundCond[0];
	if (condition == 1)
		DefineFirstBoundNodes(start, elementWithWhichLastRowBegins, locNums, addition, facet, dim);
	else if (condition == 2)
		DefineSecondBoundNodes(start, elementWithWhichLastRowBegins, locNums, addition, facet, dim);

	facet = Bound::BACK;
	start = (_xGrid.size() - 1) * (y_Grid.size() - 2);
	elementWithWhichLastRowBegins = (_xGrid.size() - 1) * ((y_Grid.size() - 1) * (_zGrid.size() - 1) - 1);
	locNums = {2, 3, 6, 7};
	condition = _boundCond[1];
	if (condition == 1)
		DefineFirstBoundNodes(start, elementWithWhichLastRowBegins, locNums, addition, facet, dim);
	else if (condition == 2)
		DefineSecondBoundNodes(start, elementWithWhichLastRowBegins, locNums, addition, facet, dim);

	facet = Bound::LEFT;
	start = 0;
	elementWithWhichLastRowBegins = (_xGrid.size() - 1) * (y_Grid.size() - 1) * (_zGrid.size() - 2);
	locNums = { 0, 2, 4, 6 };
	addition = _xGrid.size() - 1;
	dim = y_Grid.size() - 2;
	condition = _boundCond[2];
	if (condition == 1)
		DefineFirstBoundNodes(start, elementWithWhichLastRowBegins, locNums, addition, facet, dim);
	else if (condition == 2)
		DefineSecondBoundNodes(start, elementWithWhichLastRowBegins, locNums, addition, facet, dim);

	facet = Bound::RIGHT;
	start = _xGrid.size() - 2;
	elementWithWhichLastRowBegins = (_xGrid.size() - 1) * (y_Grid.size() - 1) * (_zGrid.size() - 2) + (_xGrid.size() - 2); // как я это вывела?))
	locNums = { 1, 3, 5, 7 };
	condition = _boundCond[3];
	if (condition == 1)
		DefineFirstBoundNodes(start, elementWithWhichLastRowBegins, locNums, addition, facet, dim);
	else if (condition == 2)
		DefineSecondBoundNodes(start, elementWithWhichLastRowBegins, locNums, addition, facet, dim);

	facet = Bound::LOWER;
	start = 0;
	elementWithWhichLastRowBegins = (_xGrid.size() - 1) * (y_Grid.size() - 2);
	locNums = { 0, 1, 2, 3 };
	addition = 1;
	dim = _xGrid.size() - 2;
	condition = _boundCond[4];
	if (condition == 1)
		DefineFirstBoundNodes(start, elementWithWhichLastRowBegins, locNums, addition, facet, dim);
	else if (condition == 2)
		DefineSecondBoundNodes(start, elementWithWhichLastRowBegins, locNums, addition, facet, dim);


	facet = Bound::UPPER;
	start = (_xGrid.size() - 1) * (y_Grid.size() - 1) * (_zGrid.size() - 2);
	elementWithWhichLastRowBegins = (_xGrid.size() - 1) * ((y_Grid.size() - 1) * (_zGrid.size() - 1) - 1);
	locNums = { 4, 5, 6, 7 };
	condition = _boundCond[5];
	if (condition == 1)
		DefineFirstBoundNodes(start, elementWithWhichLastRowBegins, locNums, addition, facet, dim);
	else if (condition == 2)
		DefineSecondBoundNodes(start, elementWithWhichLastRowBegins, locNums, addition, facet, dim);

	removeDuplicate(_nodesDirichlet);
}

void Solver::ApplyDirichlet()
{
	for (const auto& node : _nodesDirichlet)
	{
		auto point = GetPoint(node.p, node.s, node.r);

		_b[node.num] = u(point);
	}
}

void Solver::ReduceRow(int row)
{
	auto& di = _stiffnessMatrix.getDi();
	auto ggl = _stiffnessMatrix.getGgl();
	auto ig = _stiffnessMatrix.getIg();
	auto jg = _stiffnessMatrix.getJg();

	di[row] = 1.0;
	auto ibeg = ig[row];
	auto iend = ig[row + 1];

	for (int j = ibeg; j < iend; j++)
	{
		if (binarySearch(_nodesDirichlet, jg[j], 0, _nodesDirichlet.size() - 1) != -1) continue;
		auto multiplier = ggl[j];
		_b[jg[j]] -= _b[row] * multiplier;
	}
	for (int j = row + 1; j < _b.size(); j++)
	{
		if (binarySearch(_nodesDirichlet, j, 0, _nodesDirichlet.size() - 1) != -1) continue;
		auto jbeg = ig[j];
		auto jend = ig[j + 1];
		auto index = binarySearch(jg, row, jbeg, jend - 1);
		if (index != -1)
		{
			auto multiplier = ggl[index];
			_b[j] -= _b[row] * multiplier;
		}
	}
}


void Solver::PerformGaussianReduction()
{
	for (const auto& node : _nodesDirichlet)
	{
		auto row = node.num;
		ReduceRow(row);
	}
	for (const auto& node : _nodesDirichlet)
	{
		auto row = node;
		_stiffnessMatrix.zeroOutRowAndCol(row.num);
	}
}

void Solver::AssembleVector()
{
	vector<double> globB(_xGrid.size() * y_Grid.size() * _zGrid.size(), 0.0);

#ifdef _NORMAL_TESTS

	//// А будет с минусом
	for (auto i = 0; i < m_sources.size(); i++)
	{
		auto globNums = m_Elements[m_sources[i].elA].nodes;

		double xp = m_X[globNums[4].p]; double xp1 = m_X[globNums[5].p]; ;
		double hx = xp1 - xp;

		double ys = m_Y[globNums[4].s]; double ys1 = m_Y[globNums[6].s]; ;
		double hy = ys1 - ys;

		globB[globNums[4].num] -= X_func(0, m_sources[i].A.x, xp1, hx) * Y_func(0, m_sources[i].A.y, ys1, hy);
		globB[globNums[5].num] -= X_func(1, m_sources[i].A.x, xp, hx) * Y_func(0, m_sources[i].A.y, ys1, hy);
		globB[globNums[6].num] -= X_func(0, m_sources[i].A.x, xp1, hx) * Y_func(1, m_sources[i].A.y, ys, hy);
		globB[globNums[7].num] -= X_func(1, m_sources[i].A.x, xp, hx) * Y_func(1, m_sources[i].A.y, ys, hy);
	}

	// В будет с плюсом
	for (auto i = 0; i < m_sources.size(); i++)
	{
		auto globNums = m_Elements[m_sources[i].elB].nodes;

		double xp = m_X[globNums[4].p]; double xp1 = m_X[globNums[5].p]; ;
		double hx = xp1 - xp;

		double ys = m_Y[globNums[4].s]; double ys1 = m_Y[globNums[6].s]; ;
		double hy = ys1 - ys;

		globB[globNums[4].num] += X_func(0, m_sources[i].B.x, xp1, hx) * Y_func(0, m_sources[i].B.y, ys1, hy);
		globB[globNums[5].num] += X_func(1, m_sources[i].B.x, xp, hx) * Y_func(0, m_sources[i].B.y, ys1, hy);
		globB[globNums[6].num] += X_func(0, m_sources[i].B.x, xp1, hx) * Y_func(1, m_sources[i].B.y, ys, hy);
		globB[globNums[7].num] += X_func(1, m_sources[i].B.x, xp, hx) * Y_func(1, m_sources[i].B.y, ys, hy);
	}
#endif

#ifdef _SIMPLE_TESTS
	int nx = _xGrid.size();
	int ny = y_Grid.size();
	int nz = _zGrid.size();

	int nEl = 0;
	for (auto r = 0; r < nz - 1; r++)
	{
		for (auto s = 0; s < ny - 1; s++)
		{
			for (auto p = 0; p < nx - 1; p++)
			{
				double hx = _xGrid[p + 1] - _xGrid[p];
				double hy = y_Grid[s + 1] - y_Grid[s];
				double hz = _zGrid[r + 1] - _zGrid[r];

				auto localC = BuildLocalC(hx, hy, hz);
				auto glob_nums = _elements[nEl].nodes;

				vector<double> _f(glob_nums.size());

				for (auto i = 0; i < _f.size(); i++)
				{
					auto p = GetPoint(glob_nums[i].p, glob_nums[i].s, glob_nums[i].r);
					_f[i] = f(p);
				}

				auto localB = localC.MultMatrixByVector(_f);


				for (auto i = 0; i < _elements[nEl].nodes.size(); i++)
				{
					globB[glob_nums[i].num] += localB[i];
				}
				nEl++;
			}
		}
	}
#endif

	_b = globB;
}

void Solver::AssembleMatrix()
{
	vector<int> ig = _stiffnessMatrix.getIg();
	vector<int> jg = _stiffnessMatrix.getJg();

	vector<double> di(_xGrid.size() * y_Grid.size() * _zGrid.size(), 0.0);
	vector<double> ggl(jg.size(), 0.0);

	int nx = _xGrid.size();
	int ny = y_Grid.size();
	int nz = _zGrid.size();

	int nEl = 0;
	for (auto r = 0; r < nz - 1; r++)
	{
		for (auto s = 0; s < ny - 1; s++)
		{
			for (auto p = 0; p < nx - 1; p++)
			{
				double hx = _xGrid[p + 1] - _xGrid[p];
				double hy = y_Grid[s + 1] - y_Grid[s];
				double hz = _zGrid[r + 1] - _zGrid[r];

				Matrix local = BuildLocalG(hx, hy, hz, _elements[nEl].sigma);
				auto glob_nums = _elements[nEl].nodes;

				for (auto i = 0; i < _elements[0].nodes.size(); i++)
				{
					for (auto j = 0; j < _elements[0].nodes.size(); j++)
					{
						int glob_i = glob_nums[i].num; int glob_j = glob_nums[j].num;

						if (glob_i == glob_j)
						{
							di[glob_i] += local[i][i];
						}
						else if (glob_i > glob_j)
						{
							int ind;
							for (ind = ig[glob_i]; ig[glob_i + 1]; ind++)
							{
								if (jg[ind] == glob_j)
									break;
							}
							ggl[ind] += local[i][j];
						}
					}
				}
				nEl++;
			}
		}
	}
	_stiffnessMatrix.setGgl(ggl); _stiffnessMatrix.setDi(di);
}

std::pair<vector<int>, vector<int>> Solver::GeneratePortrait()
{
	int N = _xGrid.size() * y_Grid.size() * _zGrid.size(); // число глобальных б/ф (количество узлов)

	uint32_t bigNumberWhichIDontKnowButMoreIThinkAboutItMoreItSeemsToBeLessThanNTimes8 = N * 20;
	vector<vector<uint32_t>> list(2, std::vector<uint32_t>(bigNumberWhichIDontKnowButMoreIThinkAboutItMoreItSeemsToBeLessThanNTimes8, 0));
	std::vector<int> ig, jg;
	ig.reserve(N);
	jg.reserve(bigNumberWhichIDontKnowButMoreIThinkAboutItMoreItSeemsToBeLessThanNTimes8);
	//vector<uint32_t> listbeg(N * N, 0);
	vector<uint32_t> listbeg(N, 0);
	uint32_t listSize = 0;

	for (const auto elem : _elements)
	{
		for (uint32_t i = 0; i < elem.nodes.size(); i++)
		{
			auto k = elem.nodes[i].num;
			for (uint32_t j = i + 1; j < elem.nodes.size(); j++)
			{
				auto ind1 = k;
				auto ind2 = elem.nodes[j].num;
				if (ind2 < ind1)
				{
					ind1 = ind2;
					ind2 = k;
				}
				auto iaddr = listbeg[ind2];
				if (!iaddr)
				{
					listSize++;
					listbeg[ind2] = listSize;
					list[0][listSize] = ind1;
					list[1][listSize] = 0;
				}
				else
				{
					while (list[0][iaddr] < ind1 && list[1][iaddr] > 0)
					{
						iaddr = list[1][iaddr];
					}
					if (list[0][iaddr] > ind1)
					{
						listSize++;
						list[0][listSize] = list[0][iaddr];
						list[1][listSize] = list[1][iaddr];
						list[0][iaddr] = ind1;
						list[1][iaddr] = listSize;
					}
					else
					{
						if (list[0][iaddr] < ind1)
						{
							listSize++;
							list[1][iaddr] = listSize;
							list[0][listSize] = ind1;
							list[1][listSize] = 0;
						}
					}
				}
			}
		}
	}
	ig.push_back(0);

	for (uint32_t i = 0; i < N; i++)
	{
		ig.push_back(ig[i]);
		auto iaddr = listbeg[i];
		while (iaddr != 0)
		{
			jg.push_back(list[0][iaddr]);
			ig[i + 1]++;
			iaddr = list[1][iaddr];
		}
	}

	return { std::move(ig), std::move(jg) };
}

Matrix Solver::BuildLocalC(double hx, double hy, double hz)
{
	Matrix localC(8, 8);

	Matrix G1({ {1, -1}, {-1, 1 } });
	Matrix M1({ {1. / 3, 1. / 6}, {1. / 6, 1. / 3} });

	Matrix Gx(2, 2), Gy(2, 2), Gz(2, 2);
	Matrix Mx(2, 2), My(2, 2), Mz(2, 2);

	Gx.matrixByNum(G1, 1. / hx); Gy.matrixByNum(G1, 1. / hy); Gz.matrixByNum(G1, 1. / hz);
	Mx.matrixByNum(M1, hx); My.matrixByNum(M1, hy); Mz.matrixByNum(M1, hz);

	for (auto i = 0; i < 8; i++)
	{
		for (auto j = 0; j < 8; j++)
		{
			int ai = getFirstIndex(i); int aj = getFirstIndex(j);
			int bi = getSecondIndex(i); int bj = getSecondIndex(j);
			int ci = getThirdIndex(i); int cj = getThirdIndex(j);

			localC[i][j] = Mx[ai][aj] * My[bi][bj] * Mz[ci][cj];
		}
	}

	return localC;
}

Matrix Solver::BuildLocalG(double hx, double hy, double hz, double sigma)
{
	Matrix localG(8, 8);

	Matrix G1({ {1, -1}, {-1, 1 } });
	Matrix M1({ {1. / 3, 1. / 6}, {1. / 6, 1. / 3} });

	Matrix Gx(2, 2), Gy(2, 2), Gz(2, 2);
	Matrix Mx(2, 2), My(2, 2), Mz(2, 2);

	Gx.matrixByNum(G1, 1. / hx); Gy.matrixByNum(G1, 1. / hy); Gz.matrixByNum(G1, 1. / hz);
	Mx.matrixByNum(M1, hx); My.matrixByNum(M1, hy); Mz.matrixByNum(M1, hz);

	for (auto i = 0; i < 8; i++)
	{
		for (auto j = 0; j < 8; j++)
		{
			int ai = getFirstIndex(i); int aj = getFirstIndex(j);
			int bi = getSecondIndex(i); int bj = getSecondIndex(j);
			int ci = getThirdIndex(i); int cj = getThirdIndex(j);

			localG[i][j] = sigma * (Gx[ai][aj] * My[bi][bj] * Mz[ci][cj] +
				Mx[ai][aj] * Gy[bi][bj] * Mz[ci][cj] + Mx[ai][aj] * My[bi][bj] * Gz[ci][cj]);
		}
	}

	return localG;
}

/// <summary>
/// В соответствие локальной нумерациеи базисных 
/// функций на элементе ставит номер глобальной функции. 
/// Создает список конечных элементов. А ЕЩЕ находит, к какому 
/// конечному элементу приналежат точки источника.
/// </summary>
/// <returns></returns>
void Solver::DefineGlobalNumbers()
{
	int nx = _xGrid.size();
	int ny = y_Grid.size();
	int nz = _zGrid.size();

	_elements.resize((nx - 1) * (ny - 1) * (nz - 1));
	_nodes.resize(nx * ny * nz);

	int r = 0;
	int elN = 0;
	int numberOfElementFromWhichLastLayerStarts;

	for (r = 0; r < nz - 2; r++)
	{
		for (auto s = 0; s < ny - 1; s++)
		{
			for (auto p = 0; p < nx - 1; p++)
			{
				_elements[elN] = GetNumbers(p, s, r, nx, ny);
				SetCoord(p, s, r, _elements[elN]);

				for (auto i = 0; i < _elements[elN].nodes.size(); i++)
				{
					_nodes[_elements[elN].nodes[i].num] = GetPoint(_elements[elN].nodes[i].p, _elements[elN].nodes[i].s, _elements[elN].nodes[i].r);
				}
				elN++;
			}	
		}
	}

	numberOfElementFromWhichLastLayerStarts = elN;

	for (auto s = 0; s < ny - 1; s++)
	{
		for (auto p = 0; p < nx - 1; p++)
		{
			_elements[elN] = GetNumbers(p, s, r, nx, ny);
			SetCoord(p, s, r, _elements[elN]);


			for (auto i = 0; i < _elements[elN].nodes.size(); i++)
			{
				_nodes[_elements[elN].nodes[i].num] = GetPoint(_elements[elN].nodes[i].p, _elements[elN].nodes[i].s, _elements[elN].nodes[i].r);
			}
			elN++;
		}
	}
}

void Solver::BuildDomain(const path& _path)
{
#ifdef _NORMAL_TESTS
	buildZGrid(_path);
	buildXGrid();
	buildYGrid();
#endif

	DefineGlobalNumbers();
	FindElementForPoint();
}

void Solver::BuildYGrid()
{
	for (int i = 0; i < _ySettings.size(); i++)
	{
		double y1 = _ySettings[i].x1;
		double y2 = _ySettings[i].x2;

		int numInt = _ySettings[i].numInt;
		double koef = _ySettings[i].koef;

		y_Grid.push_back(y1);

		if (_yNesting == 1)
		{
			numInt *= 2;
			koef = sqrt(koef);
		}
		else if (_yNesting == 2)
		{
			numInt *= 4;
			koef = sqrt(sqrt(koef));
		}

		if (koef != 1.0)
		{
			for (int j = numInt; j > 1; j--)
			{
				double geomProgression = (pow(koef, j) - 1) / (koef - 1);
				double step = (y2 - y1) / geomProgression;
				y1 += step;
				y_Grid.push_back(y1);
			}
		}
		else
		{
			double step = (y2 - y1) / numInt;

			for (auto j = 0; j < numInt - 1; j++)
			{
				y1 += step;
				y_Grid.push_back(y1);
			}
		}
	}
	y_Grid.push_back(_ySettings.back().x2);

	std::ofstream fout;
	fout.open("y.txt");
	for (auto i = 0; i < y_Grid.size(); i++)
	{
		fout << y_Grid[i] << std::endl;
	}
	fout.close();
}

void Solver::BuildXGrid()
{
	for (int i = 0; i < _xSettings.size(); i++)
	{
		int numInt = _xSettings[i].numInt;
		double koef = _xSettings[i].koef;

		double x1 = _xSettings[i].x1;
		double x2 = _xSettings[i].x2;

		_xGrid.push_back(x1);

		if (_xNesting == 1)
		{
			numInt *= 2;
			koef = sqrt(koef);
		}
		else if (_xNesting == 2)
		{
			numInt *= 4;
			koef = sqrt(sqrt(koef));
		}

		if (koef != 1.0)
		{
			for (int j = numInt; j > 1; j--)
			{
				double geomProgression = (pow(koef, j) - 1) / (koef - 1);
				double step = (x2 - x1) / geomProgression;
				x1 += step;
				_xGrid.push_back(x1);
			}
		}
		else
		{
			double step = (x2 - x1) / numInt;

			for (auto j = 0; j < numInt - 1; j++)
			{
				x1 += step;
				_xGrid.push_back(x1);
			}
		}
	}
	_xGrid.push_back(_xSettings.back().x2);

	std::ofstream fout;
	fout.open("x.txt");
	for (auto i = 0; i < _xGrid.size(); i++)
	{
		fout << _xGrid[i] << std::endl;
	}
	fout.close();
}

void Solver::BuildZGrid(const path& _path)
{
	vector<double> ref_point;
	ref_point.push_back(_layers[0].z0);
	
	for (auto i = 0; i < _layers.size(); i++)
	{
		ref_point.push_back(_layers[i].z1);
	}

	std::ifstream fin;
	fin.open(_path / "z_settings.txt");

	if (fin.is_open())
	{
		int numInt;
		int nesting;
		double koef;
		double x1, x2;

		fin >> nesting;

		for (auto i = 0; i < _layers.size() - 1; i++)
		{
			fin >> numInt;
			x1 = ref_point[i]; 
			x2 = ref_point[i + 1];

			if (nesting == 1)
			{
				numInt *= 2;
			}
			else if (nesting == 2)
			{
				numInt *= 4;
			}

			double step = (x2 - x1) / numInt;
			_zGrid.push_back(x1);

			for (auto j = 0; j < numInt - 1; j++)
			{
				x1 += step;
				_zGrid.push_back(x1);
			}
		}

		x1 = *(std::prev(ref_point.end(), 2));
		x2 = ref_point.back();
		_zGrid.push_back(x1);

		fin >> numInt >> koef;

		if (nesting == 1)
		{
			numInt *= 2;
			koef = sqrt(koef);
		}
		else if (nesting == 2)
		{
			numInt *= 4;
			koef = sqrt(sqrt(koef));
		}

		for (int i = numInt; i > 1; i--)
		{
			double geomProgression = (pow(koef, i) - 1) / (koef - 1);
			double step = (x2 - x1) / geomProgression;
			x1 += step;
			_zGrid.push_back(x1);
		}
		_zGrid.push_back(x2);
		fin.close();
	}
	else
	{
		std::cerr << "File z_settings.txt was not open." << std::endl;
		std::exit(1);
	}

	std::ofstream fout;
	fout.open("z.txt");
	for (auto i = 0; i < _zGrid.size(); i++)
	{
		fout << _zGrid[i] << std::endl;
	}
	fout.close();
}

void Solver::InitDomainData(const path& _path)
{
#ifdef _SIMPLE_TESTS
	std::ifstream fin(_path / "domain.txt");
	if (fin.is_open())
	{
		auto build = [&fin](vector<double>& grid)
		{
			int size; fin >> size;
			for (auto i = 0; i < size; i++)
			{
				double num; fin >> num;
				grid.push_back(num);
			}
		};

		build(_xGrid); build(y_Grid); build(_zGrid);

		fin.close();
	}
	else
	{
		std::cerr << "File domain.txt was not open." << std::endl;
		std::exit(1);
	}
#endif

#ifdef _NORMAL_TESTS

	fin.open(_path / "x_settings.txt");
	if (fin.is_open())
	{
		int subAreas;
		int nesting;

		fin >> subAreas >> nesting;

		XNesting = nesting;
		XSettings.resize(subAreas);

		for (auto i = 0; i < XSettings.size(); i++)
		{
			double x1, x2;
			int numInt;
			double koef;

			fin >> x1 >> x2 >> numInt >> koef;
			XSettings[i] = GridSettins{ x1, x2, numInt, koef };
		}
		fin.close();
	}
	else
	{
		std::cerr << "File x_settings.txt was not open." << std::endl;
		std::exit(1);
	}

	fin.open(_path / "y_settings.txt");
	if (fin.is_open())
	{
		int subAreas;
		int nesting;

		fin >> subAreas >> nesting;
		YNesting = nesting;
		YSettings.resize(subAreas);

		for (auto i = 0; i < YSettings.size(); i++)
		{
			double y1, y2;
			int numInt;
			double koef;

			fin >> y1 >> y2 >> numInt >> koef;
			YSettings[i] = GridSettins{ y1, y2, numInt, koef };
		}
		fin.close();
	}
	else
	{
		std::cerr << "File y_settings.txt was not open." << std::endl;
		std::exit(1);
	}
#endif

	fin.open(_path / "layers.txt");
	if (fin.is_open())
	{
		std::size_t size; fin >> size;
		_layers.resize(size);

		for (auto i = 0; i < size; i++)
		{
			int indicator; fin >> indicator;
			if (indicator >= 0)
			{
				fin >> _layers[indicator].z0 >> _layers[indicator].z1 >> _layers[indicator].sigma;
			}
		}
		fin.close();
	}
	else
	{
		std::cerr << "File layers.txt was not open." << std::endl;
		std::exit(1);
	}

	fin.open(_path / "anomaly.txt");
	if (fin.is_open())
	{
		std::size_t size; fin >> size;
		_anomals.resize(size);

		for (auto i = 0; i < size; i++)
		{
			fin >> _anomals[i].layer_number >> _anomals[i].sigma;

			for (auto j = 0; j < _anomals[i].ref_plane.size(); j++)
			{
				fin >> _anomals[i].ref_plane[j].x >> _anomals[i].ref_plane[j].y;
			}
			fin >> _anomals[i].z0 >> _anomals[i].z1;
		}
		fin.close();
	}
	else
	{
		std::cerr << "File anomaly.txt was not open." << std::endl; 
		std::exit(1);
	}

	fin.open(_path / "receivers.txt");
	if (fin.is_open())
	{
		std::size_t size; fin >> size;
		_receivers.resize(size);

		for (auto i = 0; i < size; i++)
		{
			fin >> _receivers[i].M.x >> _receivers[i].M.y;
		}

		for (auto i = 0; i < size; i++)
		{
			fin >> _receivers[i].N.x >> _receivers[i].N.y;
		}
		fin.close();
	}
	else
	{
		std::cerr << "File receivers.txt was not open." << std::endl; 
		std::exit(1);
	}

	fin.open(_path / "sources.txt");
	if (fin.is_open())
	{
		std::size_t size; fin >> size;
		_sources.resize(size);

		for (auto i = 0; i < size; i++)
		{
			fin >> _sources[i].A.x >> _sources[i].A.y;
		}

		for (auto i = 0; i < size; i++)
		{
			fin >> _sources[i].B.x >> _sources[i].B.y;
		}
		fin.close();
	}
	else
	{
		std::cerr << "File sources.txt was not open." << std::endl; 
		std::exit(1);
	}

	fin.open(_path / "bounds.txt");
	if (fin.is_open())
	{
		int type;
		while (fin >> type)
		{
			_boundCond.push_back(type);
		}
		fin.close();
	}
	else
	{
		std::cerr << "File bounds.txt was not open." << std::endl;
		std::exit(1);
	}

	BuildDomain(_path);
}

#endif