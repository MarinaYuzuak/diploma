#ifndef _SOLVER_H
#define _SOLVER_H

#include <fstream>
#include <algorithm>
#include <iomanip>
//#include "Structures.h"
#include "functions.h"
//#include "lightweight.h"
#include "Matrix.h"
#include "MatrixSparse.h"
#include "Vector.h"
using std::filesystem::path;

//#define _SIMPLE_TESTS
#define _NORMAL_TESTS

class Solver
{
public:
	Solver(const path& _path)
	{
		initDomainData(_path);

		auto&& [ig, jg] = generatePortrait();
		m_A.setIg(ig); m_A.setJg(jg);

		SetSigmaForElements();

		/*int kolSloi0 = 0; int kolSloi1 = 0; int kolSloi2 = 0; int kolWithAnomal = 0;
		int kolAnomal1 = 0; int kolAnomal2 = 0;

		for (auto i = 0; i < m_Elements.size(); i++)
		{
			if (m_Elements[i].sigma == 0.03)
			{
				kolSloi0++;
			}
			else if (m_Elements[i].sigma == 0.1)
				kolSloi1++;
			else if (m_Elements[i].sigma == 0.2)
				kolSloi2++;
			else if (m_Elements[i].sigma == 0.04)
				kolAnomal1++;
			else if (m_Elements[i].sigma == 0.05)
				kolAnomal2++;
			else
			{
				std::cout << i << " " << m_Elements[i].sigma << std::endl;
				kolWithAnomal++;
			}
		}

		std::cout << "Kol-vo el v 1 sloe: " << kolSloi0 << std::endl;
		std::cout << "Kol-vo el v 2 sloe: " << kolSloi1 << std::endl;
		std::cout << "Kol-vo el v 3 sloe: " << kolSloi2 << std::endl;
		std::cout << "Kol-vo el v 1 AO: " << kolAnomal1 << std::endl;
		std::cout << "Kol-vo el v 2 AO: " << kolAnomal2 << std::endl;
		std::cout << "Kol-vo el, kot peresekautsa s AO: " << kolWithAnomal << std::endl;
		*/

		auto ak = m_Elements[7541];
		assembleMatrix();
		assembleVector();

		defineBoundsNodes();


		applyNeumann();
		applyDirichlet();

	    auto gg = m_A.getGgl(); auto di = m_A.getDi();
		auto b = m_b.getV();

		m_Solution.resize(b.size());
		CGM_with_holes_prec(ig, jg, gg, di, b, m_Solution, b.size(), 10000, 1e-15);

		/*CalculateSolution();

		std::ofstream fout;
		
		std::cout << "X Y Z" << std::endl;
		vector<double> q(m_Solution.size(), 0.0);

		vector<Point3D> points;

		for (auto r = 0; r < m_Z.size(); r++)
		{
			for (auto s = 0; s < m_Y.size(); s++)
			{
				for (auto p = 0; p < m_X.size(); p++)
				{
					double x = m_X[p]; double y = m_Y[s]; double z = m_Z[r];
					Point3D point{ x, y, z };
					points.push_back(point);
				}
			}
		}

		for (auto i = 0; i < points.size(); i++)
		{
			q[i] = u(points[i]);
		}


		for (auto i = 0; i < points.size(); i++)
		{
			std::cout << std::fixed << std::setprecision(3) << points[i].x << " " << points[i].y << " " << points[i].z << std::endl;
		}*/

		std::cout << std::endl << "solution" << std::endl;
		for (auto i = 0; i < m_Solution.size(); i++)
		{
			std::cout << std::setprecision(15) << std::scientific << m_Solution[i] << std::endl;
		}

	};
	

private:
	vector<Point2D> m_ref_plane; // класс Domain
	vector<Anomaly> m_anomal; // класс Domain
	vector<Layer> m_layers; // класс Domain
	vector<Source> m_sources; // класс Domain
	vector<Receiver> m_receivers; // класс Domain


	MatrixSparse m_A; // основной класс FemSolver (или что-то типа такого названия)
	Vector m_b; // основной класс FemSolver (или что-то типа такого названия)
	vector<double> m_Solution; // основной класс FemSolver (или что-то типа такого названия)

	vector<Node> m_nodesDirichlet; // что-то типа класса (название потом придумаю) FEMDomain
	vector<vector<Node>> m_nodesNeumann; // что-то типа класса (название потом придумаю) FEMDomain
	vector<int> m_boundCond; // что-то типа класса (название потом придумаю) FEMDomain
	vector<double> m_Z; // что-то типа класса (название потом придумаю) FEMDomain
	vector<double> m_X; // что-то типа класса (название потом придумаю) FEMDomain
	vector<double> m_Y; // что-то типа класса (название потом придумаю) FEMDomain
	vector<Element> m_Elements; // что-то типа класса (название потом придумаю) FEMDomain


private:
	void initDomainData(const path& _path); // класс Domain


	void CalculateSolution(); // основной класс FemSolver (или что-то типа такого названия)
	void assembleMatrix(); // основной класс FemSolver (или что-то типа такого названия)
	void assembleVector(); // основной класс FemSolver (или что-то типа такого названия)
	Matrix buildLocalG(double hx, double hy, double hz, double sigma); // основной класс FemSolver (или что-то типа такого названия)
	Matrix buildLocalC(double hx, double hy, double hz); // основной класс FemSolver (или что-то типа такого названия)
	void applyDirichlet(); // основной класс FemSolver (или что-то типа такого названия)
	void applyNeumann(); // основной класс FemSolver (или что-то типа такого названия)
	void performGaussianReduction(); // это тоже будет в FEMSolver, но приватной?
	void reduceRow(int row); // это тоже будет в FEMSolver, но приватной?
	std::pair<vector<int>, vector<int>> generatePortrait(); // это тоже будет в FEMSolver, но приватной? или просто в Functions ее запихать?


	void defineBoundsNodes(); // что-то типа класса (название потом придумаю) FEMDomain
	void defineGlobalNumbers(); // что-то типа класса (название потом придумаю) FEMDomain
	void SetSigmaForElements(); // что-то типа класса (название потом придумаю) FEMDomain
	void defineFirstBoundNodes(int start, int elementWithWhichLastRowBegins, vector<int>& locNums, int addition, Bound facet, int dim); // что-то типа класса (название потом придумаю) FEMDomain
	void defineSecondBoundNodes(int start, int elementWithWhichLastRowBegins, vector<int>& locNums, int addition, Bound facet, int dim); // что-то типа класса (название потом придумаю) FEMDomain
	void buildZGrid(const path& _path); // что-то типа класса (название потом придумаю) FEMDomain
	void buildXGrid(const path& _path); // что-то типа класса (название потом придумаю) FEMDomain
	void buildYGrid(const path& _path); // что-то типа класса (название потом придумаю) FEMDomain
	void buildDomain(const path& _path); // класс FEMDomain



	double GetVForIntersectionLayer(double h, Element e); // статический класс Functions
	double GetVForIntersectionAnomaly(Element e, Anomaly a); // статический класс Functions
	double GetVElement(Element e); // статический класс Functions
	Point3D GetPoint(int p, int s, int r); // статический класс Functions
	vector<Point3D> GetAnomalyPoints(Anomaly a); // статический класс Functions

};

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
	vector<Point3D> points; // теперь тут координаты всех узлов конечного элемента

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
	double sOsn = (m_X[e.nodes[1].p] - m_X[e.nodes[0].p]) * (m_Y[e.nodes[2].s] - m_Y[e.nodes[0].s]);
	double h = m_Z[e.nodes[4].r] - m_Z[e.nodes[0].r];
	return (1. / 3) * sOsn * h;
}

double Solver::GetVForIntersectionLayer(double h, Element e)
{
	double sOsn = (m_X[e.nodes[1].p] - m_X[e.nodes[0].p]) * (m_Y[e.nodes[2].s] - m_Y[e.nodes[0].s]);
	return (1. / 3.) * sOsn * h;
}

void Solver::SetSigmaForElements()
{
	vector<vector<Anomaly>> anomaliesWithWhichFiniteElementIntersects(m_Elements.size());
	//int inNulLayer = 0; int inFirstLayer = 0;

	for (auto i = 0; i < m_Elements.size(); i++)
	{
		bool intersectWithAnomal = false;
		for (auto j = 0; j < m_anomal.size(); j++)
		{
			int indicator = 0;
			vector<Point3D> pointsEl(8);
			vector<Point3D> pointsAnomaly = GetAnomalyPoints(m_anomal[j]);

			for (auto l = 0; l < m_Elements[i].nodes.size(); l++)
			{
				auto p = GetPoint(m_Elements[i].nodes[l].p, m_Elements[i].nodes[l].s, m_Elements[i].nodes[l].r);
				pointsEl[l] = p;
				indicator += isPointInsideAnomaly(p, m_anomal[j]);
			}

			// КЭ полностью лежит в АО.
			if (indicator == 8)
			{
				m_Elements[i].sigma = m_anomal[j].sigma;
				intersectWithAnomal = true;
				break; // зачем нам дальше проверять принадлежность КЭ другим АО, если он уже полностью содержится в одном из АО.
			}

			if (doParallelepipedsIntersect(pointsEl, pointsAnomaly))
			{
				intersectWithAnomal = true;
				anomaliesWithWhichFiniteElementIntersects[i].push_back(m_anomal[j]);

				//if (j == 0)
				//{
				//	std::cout << "intersect with AO in null-layer " << ++inNulLayer << std::endl;
				//}
				//else if (j == 1)
				//{
				//	std::cout << "intersect with AO in first-layer " << ++inFirstLayer << std::endl;
				//}
			}
		}

		if (!intersectWithAnomal)
		{
			/*смотрим, какому слою (или слоям) принадлежит КЭ*/
			for (auto l = 0; l < m_layers.size(); l++)
			{
				bool isBelong = isElementBelongsLayer(m_layers[l], m_Z[m_Elements[i].nodes[0].r], m_Z[m_Elements[i].nodes[4].r]);
				if (isBelong)
				{
					m_Elements[i].sigma = m_layers[l].sigma;
					isBelong = true;
					break;
				}
			}
		}
	}

	/*int num = 0;
	for (auto k = 0; k < anomaliesWithWhichFiniteElementIntersects.size(); k++)
	{
		if (anomaliesWithWhichFiniteElementIntersects[k].size() != 0)
			num++;
	}

	std::cout << std::endl << "kol-vo el kot intersect s ao: " << num << std::endl;*/

	for (auto i = 0; i < anomaliesWithWhichFiniteElementIntersects.size(); i++)
	{
		if (anomaliesWithWhichFiniteElementIntersects[i].size() != 0) // считаем долю объема от каждого АО и от каждого слоя, содержащего АО.
		{
			/*значит, что КЭ лежит в одном слое и пересекается с каким-то количеством АО.*/
			auto layer = anomaliesWithWhichFiniteElementIntersects[i][0].layer_number;
			double vEl = GetVElement(m_Elements[i]);
			double vFromLayer = vEl;

			m_Elements[i].sigma = 0;
			for (auto l = 0; l < anomaliesWithWhichFiniteElementIntersects[i].size(); l++)
			{
				double V = GetVForIntersectionAnomaly(m_Elements[i], anomaliesWithWhichFiniteElementIntersects[i][l]);
				vFromLayer -= V;
				m_Elements[i].sigma += (V / vEl) * anomaliesWithWhichFiniteElementIntersects[i][l].sigma;
			}
			vFromLayer = vFromLayer / vEl;
			m_Elements[i].sigma += vFromLayer * m_layers[layer].sigma;	
		}
	}
}

Point3D Solver::GetPoint(int p, int s, int r)
{
	auto x = m_X[p];
	auto y = m_Y[s];
	auto z = m_Z[r];
	Point3D point{ x, y, z };

	return point;
}

void Solver::applyNeumann()
{
#ifdef _SIMPLE_TESTS
	Matrix M({ {4, 2, 2, 1}, {2, 4, 1, 2}, {2, 1, 4, 2}, {1, 2, 2, 4} });
	vector<double> t(4, 0.0);
	auto b = m_b.getV();

	for (auto i = 0; i < m_nodesNeumann.size(); i++)
	{
		for (auto j = 0; j < m_nodesNeumann[i].size(); j++)
		{
			Point3D p = {m_X[m_nodesNeumann[i][j].p], m_Y[m_nodesNeumann[i][j].s], m_Z[m_nodesNeumann[i][j].r] };
			t[j] = theta(p, m_nodesNeumann[i][j].bound);
		}

		double hp, hr;
		switch (m_nodesNeumann[i][0].bound)
		{
		case Bound::UPPER:
			/// ну и хуйня, конечно, ПИЗДЕЦ. переделать!!!!!!!!
			hp = m_X[m_nodesNeumann[i][1].p] - m_X[m_nodesNeumann[i][0].p];
			hr = m_Y[m_nodesNeumann[i][2].s] - m_Y[m_nodesNeumann[i][0].s];
			break;
		case Bound::FRONT:
			hp = m_X[m_nodesNeumann[i][1].p] - m_X[m_nodesNeumann[i][0].p];
			hr = m_Z[m_nodesNeumann[i][2].r] - m_Z[m_nodesNeumann[i][0].r];
			break;
		case Bound::BACK:
			hp = m_X[m_nodesNeumann[i][1].p] - m_X[m_nodesNeumann[i][0].p];
			hr = m_Z[m_nodesNeumann[i][2].r] - m_Z[m_nodesNeumann[i][0].r];
			break;
		case Bound::LEFT:
			hp = m_Y[m_nodesNeumann[i][1].s] - m_Y[m_nodesNeumann[i][0].s];
			hr = m_Z[m_nodesNeumann[i][2].r] - m_Z[m_nodesNeumann[i][0].r];
			break;
		case Bound::RIGHT:
			hp = m_Y[m_nodesNeumann[i][1].s] - m_Y[m_nodesNeumann[i][0].s];
			hr = m_Z[m_nodesNeumann[i][2].r] - m_Z[m_nodesNeumann[i][0].r];
			break;
		case Bound::LOWER:
			hp = m_X[m_nodesNeumann[i][1].p] - m_X[m_nodesNeumann[i][0].p];
			hr = m_Y[m_nodesNeumann[i][2].s] - m_Y[m_nodesNeumann[i][0].s];
			break;
		default:
			break;
		}

		auto bLocal = M.MultMatrixByVector(t);
		bLocal = MultVectorByNum(bLocal, (hp * hr) / 36.);

		for (auto j = 0; j < m_nodesNeumann[i].size(); j++)
		{
			b[m_nodesNeumann[i][j].num] += bLocal[j];
		}
	}
	m_b.setV(b);
#endif

#ifdef _NORMAL_TESTS
	std::cout << std::endl;
#endif
}

void Solver::CalculateSolution()
{
	auto b = m_b.getV();
	auto size = b.size();

	m_Solution.resize(b.size(), 0.0);
	auto eps = 1e-21;
	auto maxIter = static_cast<uint32_t>(1e+7);
	auto loadVectorNorm = sqrt(dotProduct(b, b));
	auto z = b;
	auto r = b;

	auto p = m_A.MultiplyByVector(z);

	auto rDotR = dotProduct(r, r);
	auto residual = sqrt(rDotR) / loadVectorNorm;

	for (uint32_t k = 1; k < maxIter && residual > eps; k++)
	{
		auto pDotP = dotProduct(p, p);
		auto a = dotProduct(p, r) / pDotP;

		for (uint32_t i = 0; i < b.size(); i++)
		{
			m_Solution[i] += a * z[i];
			r[i] -= a * p[i];
		}
		auto Ar = m_A.MultiplyByVector(r);
		auto b = -dotProduct(p, Ar) / pDotP;

		for (auto i = 0; i < size; i++)
		{
			z[i] = r[i] + b * z[i];
			p[i] = Ar[i] + b * p[i];
		}
		rDotR = dotProduct(r, r);
		residual = sqrt(rDotR) / loadVectorNorm;
	}
}

void Solver::defineSecondBoundNodes(int start, int elementWithWhichLastRowBegins, vector<int>& locNums, int addition, Bound facet, int dim)
{
	int endRow = start + addition * dim;

	while (start < elementWithWhichLastRowBegins + addition)
	{
		for (auto i = start; i <= endRow; i += addition)
		{
			vector<Node> nodes(locNums.size());
			for (auto j = 0; j < locNums.size(); j++)
			{
				nodes[j] = m_Elements[i].nodes[locNums[j]];
				nodes[j].bound = facet;
			}

			m_nodesNeumann.push_back(nodes);
		}

		switch (facet)
		{
		case Bound::FRONT:
			start += (m_X.size() - 1) * (m_Y.size() - 1);
			break;
		case Bound::BACK:
			start += (m_X.size() - 1) * (m_Y.size() - 1);
			break;
		default:
			start = endRow + addition;
			break;
		}
		endRow = start + addition * dim;
	}
}

void Solver::defineFirstBoundNodes(int start, int elementWithWhichLastRowBegins, vector<int>& locNums, int addition, Bound facet, int dim)
{
	int endRow = start + addition * dim;

	//while (start < m_finitElements.size())
	while(start < elementWithWhichLastRowBegins + addition)
	{
		m_nodesDirichlet.push_back(m_Elements[start].nodes[locNums[0]]);
		m_nodesDirichlet.back().bound = facet;

		for (auto i = start; i <= endRow; i += addition)
		{
			m_nodesDirichlet.push_back(m_Elements[i].nodes[locNums[1]]);
			m_nodesDirichlet.back().bound = facet;
		}

		switch (facet)
		{
		case Bound::FRONT:
			start += (m_X.size() - 1) * (m_Y.size() - 1);
			break;
		case Bound::BACK:
			start += (m_X.size() - 1) * (m_Y.size() - 1);
			break;
		default:
			start = endRow + addition;
			break;
		}
		endRow = start + addition * dim;
	}

	start = elementWithWhichLastRowBegins;
	endRow = start + addition * dim;

	m_nodesDirichlet.push_back(m_Elements[start].nodes[locNums[2]]);
	m_nodesDirichlet.back().bound = facet;
	for (auto i = start; i <= endRow; i += addition)
	{
		m_nodesDirichlet.push_back(m_Elements[i].nodes[locNums[3]]);
		m_nodesDirichlet.back().bound = facet;
	}
}

// не функция а просто ужасс
void Solver::defineBoundsNodes()
{
	Bound facet = Bound::FRONT;
	int start = 0;
	int elementWithWhichLastRowBegins = (m_X.size() - 1) * (m_Y.size() - 1) * (m_Z.size() - 2);
	vector<int> locNums = { 0, 1, 4, 5 };
	int addition = 1;
	int dim = m_X.size() - 2;
	int condition = m_boundCond[0];
	if (condition == 1)
		defineFirstBoundNodes(start, elementWithWhichLastRowBegins, locNums, addition, facet, dim);
	else if (condition == 2)
		defineSecondBoundNodes(start, elementWithWhichLastRowBegins, locNums, addition, facet, dim);

	facet = Bound::BACK;
	start = (m_X.size() - 1) * (m_Y.size() - 2);
	elementWithWhichLastRowBegins = (m_X.size() - 1) * ((m_Y.size() - 1) * (m_Z.size() - 1) - 1);
	locNums = {2, 3, 6, 7};
	condition = m_boundCond[1];
	if (condition == 1)
		defineFirstBoundNodes(start, elementWithWhichLastRowBegins, locNums, addition, facet, dim);
	else if (condition == 2)
		defineSecondBoundNodes(start, elementWithWhichLastRowBegins, locNums, addition, facet, dim);

	facet = Bound::LEFT;
	start = 0;
	elementWithWhichLastRowBegins = (m_X.size() - 1) * (m_Y.size() - 1) * (m_Z.size() - 2);
	locNums = { 0, 2, 4, 6 };
	addition = m_X.size() - 1;
	dim = m_Y.size() - 2; // поч не эм икс??
	condition = m_boundCond[2];
	if (condition == 1)
		defineFirstBoundNodes(start, elementWithWhichLastRowBegins, locNums, addition, facet, dim);
	else if (condition == 2)
		defineSecondBoundNodes(start, elementWithWhichLastRowBegins, locNums, addition, facet, dim);

	facet = Bound::RIGHT;
	start = m_X.size() - 2;
	elementWithWhichLastRowBegins = (m_X.size() - 1) * (m_Y.size() - 1) * (m_Z.size() - 2) + (m_X.size() - 2); // как я это вывела?))
	locNums = { 1, 3, 5, 7 };
	condition = m_boundCond[3];
	if (condition == 1)
		defineFirstBoundNodes(start, elementWithWhichLastRowBegins, locNums, addition, facet, dim);
	else if (condition == 2)
		defineSecondBoundNodes(start, elementWithWhichLastRowBegins, locNums, addition, facet, dim);

	facet = Bound::LOWER;
	start = 0;
	elementWithWhichLastRowBegins = (m_X.size() - 1) * (m_Y.size() - 2);
	locNums = { 0, 1, 2, 3 };
	addition = 1;
	dim = m_X.size() - 2;
	condition = m_boundCond[4];
	if (condition == 1)
		defineFirstBoundNodes(start, elementWithWhichLastRowBegins, locNums, addition, facet, dim);
	else if (condition == 2)
		defineSecondBoundNodes(start, elementWithWhichLastRowBegins, locNums, addition, facet, dim);


	facet = Bound::UPPER;
	start = (m_X.size() - 1) * (m_Y.size() - 1) * (m_Z.size() - 2);
	elementWithWhichLastRowBegins = (m_X.size() - 1) * ((m_Y.size() - 1) * (m_Z.size() - 1) - 1);
	locNums = { 4, 5, 6, 7 };
	condition = m_boundCond[5];
	if (condition == 1)
		defineFirstBoundNodes(start, elementWithWhichLastRowBegins, locNums, addition, facet, dim);
	else if (condition == 2)
		defineSecondBoundNodes(start, elementWithWhichLastRowBegins, locNums, addition, facet, dim);

	removeDuplicate(m_nodesDirichlet);
}

void Solver::applyDirichlet()
{
	auto b = m_b.getV();

	for (auto node : m_nodesDirichlet)
	{
		auto point = GetPoint(node.p, node.s, node.r);

		b[node.num] = u(point);
	}

	m_b.setV(b);

	performGaussianReduction();
}

void Solver::reduceRow(int row)
{
	auto& di = m_A.getDi();
	auto ggl = m_A.getGgl();
	auto ig = m_A.getIg();
	auto jg = m_A.getJg();

	di[row] = 1.0;
	auto ibeg = ig[row];
	auto iend = ig[row + 1];

	auto b = m_b.getV();

	for (int j = ibeg; j < iend; j++)
	{
		if (binarySearch(m_nodesDirichlet, jg[j], 0, m_nodesDirichlet.size() - 1) != -1) continue;
		auto multiplier = ggl[j];
		b[jg[j]] -= b[row] * multiplier;
	}
	for (int j = row + 1; j < b.size(); j++)
	{
		if (binarySearch(m_nodesDirichlet, j, 0, m_nodesDirichlet.size() - 1) != -1) continue;
		auto jbeg = ig[j];
		auto jend = ig[j + 1];
		auto index = binarySearch(jg, row, jbeg, jend - 1);
		if (index != -1)
		{
			auto multiplier = ggl[index];
			b[j] -= b[row] * multiplier;
		}
	}

	m_b.setV(b);
}


void Solver::performGaussianReduction()
{
	for (const auto& node : m_nodesDirichlet)
	{
		auto row = node.num;
		reduceRow(row);
	}
	for (const auto& node : m_nodesDirichlet)
	{
		auto row = node;
		m_A.zeroOutRowAndCol(row.num);
	}
}

// пока что учет только точки А!!! СПРОСИТЬ
void Solver::assembleVector()
{
	vector<double> globB(m_X.size() * m_Y.size() * m_Z.size(), 0.0);

#ifdef _NORMAL_TESTS
	// А будет с плюсом
	for (auto i = 0; i < m_sources.size(); i++)
	{
		auto globNums = m_Elements[m_sources[i].elA].nodes;

		double xp = m_X[globNums[4].p]; double xp1 = m_X[globNums[5].p]; ;
		double hx = xp1 - xp;

		double ys = m_Y[globNums[4].s]; double ys1 = m_Y[globNums[6].s]; ;
		double hy = ys1 - ys;

		globB[globNums[4].num] += X_func(0, m_sources[i].A.x, xp1, hx) * Y_func(0, m_sources[i].A.y, ys1, hy);
		globB[globNums[5].num] += X_func(1, m_sources[i].A.x, xp, hx) * Y_func(0, m_sources[i].A.y, ys1, hy);
		globB[globNums[6].num] += X_func(0, m_sources[i].A.x, xp1, hx) * Y_func(1, m_sources[i].A.y, ys, hy);
		globB[globNums[7].num] += X_func(1, m_sources[i].A.x, xp, hx) * Y_func(1, m_sources[i].A.y, ys, hy);
	}

	// В будет с минусом
	for (auto i = 0; i < m_sources.size(); i++)
	{
		auto globNums = m_Elements[m_sources[i].elB].nodes;

		double xp = m_X[globNums[4].p]; double xp1 = m_X[globNums[5].p]; ;
		double hx = xp1 - xp;

		double ys = m_Y[globNums[4].s]; double ys1 = m_Y[globNums[6].s]; ;
		double hy = ys1 - ys;

		globB[globNums[4].num] -= X_func(0, m_sources[i].B.x, xp1, hx) * Y_func(0, m_sources[i].B.y, ys1, hy);
		globB[globNums[5].num] -= X_func(1, m_sources[i].B.x, xp, hx) * Y_func(0, m_sources[i].B.y, ys1, hy);
		globB[globNums[6].num] -= X_func(0, m_sources[i].B.x, xp1, hx) * Y_func(1, m_sources[i].B.y, ys, hy);
		globB[globNums[7].num] -= X_func(1, m_sources[i].B.x, xp, hx) * Y_func(1, m_sources[i].B.y, ys, hy);
	}
#endif

#ifdef _SIMPLE_TESTS
	int nx = m_X.size();
	int ny = m_Y.size();
	int nz = m_Z.size();

	int nEl = 0;
	for (auto r = 0; r < nz - 1; r++)
	{
		for (auto s = 0; s < ny - 1; s++)
		{
			for (auto p = 0; p < nx - 1; p++)
			{
				double hx = m_X[p + 1] - m_X[p];
				double hy = m_Y[s + 1] - m_Y[s];
				double hz = m_Z[r + 1] - m_Z[r];

				auto localC = buildLocalC(hx, hy, hz);
				auto glob_nums = m_Elements[nEl].nodes;

				vector<double> _f(glob_nums.size());

				for (auto i = 0; i < _f.size(); i++)
				{
					auto p = GetPoint(glob_nums[i].p, glob_nums[i].s, glob_nums[i].r);
					_f[i] = f(p);
				}

				auto localB = localC.MultMatrixByVector(_f);


				for (auto i = 0; i < m_Elements[nEl].nodes.size(); i++)
				{
					globB[glob_nums[i].num] += localB[i];
				}
				nEl++;
			}
		}
	}
#endif

	m_b.setV(globB); 
}

void Solver::assembleMatrix()
{
	vector<int> ig = m_A.getIg();
	vector<int> jg = m_A.getJg();

	vector<double> di(m_X.size() * m_Y.size() * m_Z.size(), 0.0);
	vector<double> ggl(jg.size(), 0.0);

	int nx = m_X.size();
	int ny = m_Y.size();
	int nz = m_Z.size();

	int nEl = 0;
	for (auto r = 0; r < nz - 1; r++)
	{
		for (auto s = 0; s < ny - 1; s++)
		{
			for (auto p = 0; p < nx - 1; p++)
			{
				double hx = m_X[p + 1] - m_X[p];
				double hy = m_Y[s + 1] - m_Y[s];
				double hz = m_Z[r + 1] - m_Z[r];

				Matrix local = buildLocalG(hx, hy, hz, m_Elements[nEl].sigma);
				auto glob_nums = m_Elements[nEl].nodes;

				for (auto i = 0; i < m_Elements[0].nodes.size(); i++)
				{
					for (auto j = 0; j < m_Elements[0].nodes.size(); j++)
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
	m_A.setGgl(ggl); m_A.setDi(di);
}

std::pair<vector<int>, vector<int>> Solver::generatePortrait()
{
	int N = m_X.size() * m_Y.size() * m_Z.size(); // число глобальных б/ф (количество узлов)

	uint32_t bigNumberWhichIDontKnowButMoreIThinkAboutItMoreItSeemsToBeLessThanNTimes8 = N * 20;
	vector<vector<uint32_t>> list(2, std::vector<uint32_t>(bigNumberWhichIDontKnowButMoreIThinkAboutItMoreItSeemsToBeLessThanNTimes8, 0));
	std::vector<int> ig, jg;
	ig.reserve(N);
	jg.reserve(bigNumberWhichIDontKnowButMoreIThinkAboutItMoreItSeemsToBeLessThanNTimes8);
	vector<uint32_t> listbeg(N * N, 0);
	uint32_t listSize = 0;

	for (const auto elem : m_Elements)
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

Matrix Solver::buildLocalC(double hx, double hy, double hz)
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

Matrix Solver::buildLocalG(double hx, double hy, double hz, double sigma)
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

	//for (auto i = 0; i < 8; i++)
	//{
	//	for (auto j = 0; j < 8; j++)
	//	{
	//		std::cout << localG[i][j] << " ";
	//	}
	//	std::cout << std::endl;
	//}
	//std::cout << std::endl;
	return localG;
}

/// <summary>
/// В соответствие локальной нумерациеи базисных 
/// функций на элементе ставит номер глобальной функции. 
/// Создает список конечных элементов. А ЕЩЕ находит, к какому 
/// конечному элементу приналежат точки источника.
/// </summary>
/// <returns></returns>
void Solver::defineGlobalNumbers()
{
	int nx = m_X.size();
	int ny = m_Y.size();
	int nz = m_Z.size();

	m_Elements.resize((nx - 1) * (ny - 1) * (nz - 1));

	int r = 0;
	int elN = 0;
	int numberOfElementFromWhichLastLayerStarts;

	for (r = 0; r < nz - 2; r++)
	{
		for (auto s = 0; s < ny - 1; s++)
		{
			for (auto p = 0; p < nx - 1; p++)
			{
				m_Elements[elN] = GetNumbers(p, s, r, nx, ny);
				SetCoord(p, s, r, m_Elements[elN++]);
			}	
		}
	}

	numberOfElementFromWhichLastLayerStarts = elN;

	for (auto s = 0; s < ny - 1; s++)
	{
		for (auto p = 0; p < nx - 1; p++)
		{
			m_Elements[elN] = GetNumbers(p, s, r, nx, ny);
			SetCoord(p, s, r, m_Elements[elN++]);
		}
	}

	for (auto i = 0; i < m_sources.size(); i++)
	{
		bool is_found_A = false;
		bool is_found_B = false;
		int whereToSearchPoint = numberOfElementFromWhichLastLayerStarts;

		for (auto s = 0; s < ny - 1; s++)
		{
			for (auto p = 0; p < nx - 1; p++)
			{
				if (!is_found_A)
				{
					is_found_A = isPointInside(m_X[p], m_X[p + 1], m_Y[s], m_Y[s + 1], m_sources[i].A);
					if (is_found_A)
					{
						m_sources[i].elA = whereToSearchPoint;
					}
				}

				if (!is_found_B)
				{
					is_found_B = isPointInside(m_X[p], m_X[p + 1], m_Y[s], m_Y[s + 1], m_sources[i].B);
					if (is_found_B)
					{
						m_sources[i].elB = whereToSearchPoint;
					}
				}

				if (is_found_A && is_found_B)
					break;

				whereToSearchPoint++; 
			}
			if (is_found_A && is_found_B) // пипец))))
				break;
		}
	}
}

void Solver::buildDomain(const path& _path)
{
#ifdef _NORMAL_TESTS
	buildZGrid(_path);
	buildXGrid(_path);
	buildYGrid(_path);
#endif

	defineGlobalNumbers();
}

void Solver::buildYGrid(const path& _path)
{
	std::ifstream fin;
	double y1, y2;
	int numInt;
	int subAreas;
	int nesting;
	double koef;

	fin.open(_path / "y_settings.txt");

	if (fin.is_open())
	{
		fin >> subAreas >> nesting;
		for (int i = 0; i < subAreas; i++)
		{
			fin >> y1 >> y2 >> numInt >> koef;
			m_Y.push_back(y1);

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

			for (int j = numInt; j > 1; j--)
			{
				double geomProgression = (pow(koef, j) - 1) / (koef - 1);
				double step = (y2 - y1) / geomProgression;
				y1 += step;
				m_Y.push_back(y1);
			}
		}
		m_Y.push_back(y2);
		fin.close();
	}
	else
	{
		std::cerr << "File y_settings.txt was not open." << std::endl;
		std::exit(1);
	}

	std::ofstream fout;
	fout.open("y.txt");
	for (auto i = 0; i < m_Y.size(); i++)
	{
		fout << m_Y[i] << std::endl;
	}
	fout.close();
}

void Solver::buildXGrid(const path& _path)
{
	std::ifstream fin;
	double x1, x2;
	int numInt;
	int subAreas;
	int nesting;
	double koef;

	fin.open(_path / "x_settings.txt");

	if (fin.is_open())
	{
		fin >> subAreas >> nesting;

		for (int i = 0; i < subAreas; i++)
		{
			fin >> x1 >> x2 >> numInt >> koef;
			m_X.push_back(x1);

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

			for (int j = numInt; j > 1; j--)
			{
				double geomProgression = (pow(koef, j) - 1) / (koef - 1);
				double step = (x2 - x1) / geomProgression;
				x1 += step;
				m_X.push_back(x1);
			}
		}
		m_X.push_back(x2);
		fin.close();
	}
	else
	{
		std::cerr << "File x_settings.txt was not open." << std::endl;
		std::exit(1);
	}

	std::ofstream fout;
	fout.open("x.txt");
	for (auto i = 0; i < m_X.size(); i++)
	{
		fout << m_X[i] << std::endl;
	}
	fout.close();
}

void Solver::buildZGrid(const path& _path)
{
	vector<double> ref_point;
	ref_point.push_back(m_layers[0].z0);
	
	for (auto i = 0; i < m_layers.size(); i++)
	{
		ref_point.push_back(m_layers[i].z1);
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

		for (auto i = 0; i < m_layers.size() - 1; i++)
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
			m_Z.push_back(x1);

			for (auto j = 0; j < numInt - 1; j++)
			{
				x1 += step;
				m_Z.push_back(x1);
			}
		}

		x1 = *(std::prev(ref_point.end(), 2));
		x2 = ref_point.back();
		m_Z.push_back(x1);

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
			m_Z.push_back(x1);
		}
		m_Z.push_back(x2);
		fin.close();
	}
	else
	{
		std::cerr << "File z_settings.txt was not open." << std::endl;
		std::exit(1);
	}

	std::ofstream fout;
	fout.open("z.txt");
	for (auto i = 0; i < m_Z.size(); i++)
	{
		fout << m_Z[i] << std::endl;
	}
	fout.close();
}

void Solver::initDomainData(const path& _path)
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

		build(m_X); build(m_Y); build(m_Z);

		fin.close();
	}
	else
	{
		std::cerr << "File domain.txt was not open." << std::endl;
		std::exit(1);
	}
#endif

#ifdef _NORMAL_TESTS
	std::ifstream fin(_path / "layers.txt");
	if (fin.is_open())
	{
		std::size_t size; fin >> size;
		m_layers.resize(size);

		for (auto i = 0; i < size; i++)
		{
			int indicator; fin >> indicator;
			if (indicator >= 0)
			{
				fin >> m_layers[indicator].z0 >> m_layers[indicator].z1 >> m_layers[indicator].sigma;
			}
		}
		fin.close();
	}
	else
	{
		std::cerr << "File layers.txt was not open." << std::endl; 
		std::exit(1);
	}
#endif

	fin.open(_path / "anomaly.txt");
	if (fin.is_open())
	{
		std::size_t size; fin >> size;
		m_anomal.resize(size);

		for (auto i = 0; i < size; i++)
		{
			fin >> m_anomal[i].layer_number >> m_anomal[i].sigma;

			for (auto j = 0; j < m_anomal[i].ref_plane.size(); j++)
			{
				fin >> m_anomal[i].ref_plane[j].x >> m_anomal[i].ref_plane[j].y;
			}
			fin >> m_anomal[i].z0 >> m_anomal[i].z1;
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
		m_receivers.resize(size);

		for (auto i = 0; i < size; i++)
		{
			fin >> m_receivers[i].M.x >> m_receivers[i].M.y;
		}

		for (auto i = 0; i < size; i++)
		{
			fin >> m_receivers[i].N.x >> m_receivers[i].N.y;
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
		m_sources.resize(size);

		for (auto i = 0; i < size; i++)
		{
			fin >> m_sources[i].A.x >> m_sources[i].A.y;
		}

		for (auto i = 0; i < size; i++)
		{
			fin >> m_sources[i].B.x >> m_sources[i].B.y;
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
			m_boundCond.push_back(type);
		}
	}
	else
	{
		std::cerr << "File bounds.txt was not open." << std::endl;
		std::exit(1);
	}
	buildDomain(_path);
}
#endif