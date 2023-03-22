#ifndef _SOLVER_H
#define _SOLVER_H

#include <fstream>
#include <algorithm>
#include "Structures.h"
#include "lightweight.h"
#include "Matrix.h"
#include "MatrixSparse.h"
#include "Vector.h"
using std::filesystem::path;

#define _SIMPLE_TESTS
//#define _NORMAL_TESTS

class Solver
{
public:
	Solver(const path& _path)
	{
		initDomainData(_path);
		buildDomain();

		auto&& [ig, jg] = generatePortrait();
		m_A.setIg(ig); m_A.setJg(jg);


		std::cout << std::endl << std::endl;
		assembleMatrix();
		assembleVector();

		m_nodesDirichlet = nodesWithFirstConditions();

		std::cout << "b:" << std::endl;
		auto _b = m_b.getV();
		for (auto i = 0; i < _b.size(); i++)
			std::cout << _b[i] << " ";

		applyDirichlet();

		std::cout << std::endl;
		std::cout << std::endl;

		std::cout << "b:" << std::endl;
		 _b = m_b.getV();
		for (auto i = 0; i < _b.size(); i++)
			std::cout << _b[i] << " ";

		std::cout << std::endl;
		std::cout << std::endl;

		std::cout << "ggl:" << std::endl;
		auto ggl = m_A.getGgl();
		for (auto i = 0; i < ggl.size(); i++)
			std::cout << ggl[i] << std::endl;
		std::cout << ggl.size();

		std::cout << "\n\ndi:" << std::endl;
		auto di = m_A.getDi();
		for (auto i = 0; i < di.size(); i++)
			std::cout << di[i] << std::endl;
	};
	

private:
	vector<Point3D> m_ref_plane; // Plane of reference.
	vector<Anomaly> m_anomal;
	vector<Layer> m_layers;
	vector<Source> m_sources;
	vector<Receiver> m_receivers;

	vector<Element> m_finitElements; // можно и лучше название, наверное))

	MatrixSparse m_A;
	Vector m_b;

	vector<int> m_nodesDirichlet;

	vector<double> m_zGrid;
	vector<double> m_xGrid;
	vector<double> m_yGrid;


private:
	void initDomainData(const path& _path);
	void buildZGrid();
	void defineGlobalNumbers();
	void buildDomain(); // название потом придумаю. чтобы в конструкторе не было всякой лишней фигни типа построения сетки

	Matrix buildLocalG(double hx, double hy, double hz, double sigma);
	Matrix buildLocalC(double hx, double hy, double hz);

	std::pair<vector<int>, vector<int>> generatePortrait(); // наверное, стоит ее отсюда убрать, ну потому что что-то не так...

	void assembleMatrix();
	void assembleVector();

	vector<int> nodesWithFirstConditions();

	void reduceRow(int row);

	void performGaussianReduction();


	void applyDirichlet();

	vector<double> buildXGrid();
	vector<double> buildYGrid();
};

void Solver::applyDirichlet()
{
	auto b = m_b.getV();

	for (auto node : m_nodesDirichlet)
	{
		b[node] = 0.0;
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
		auto row = node;
		reduceRow(row);
	}
	for (const auto& node : m_nodesDirichlet)
	{
		auto row = node;
		m_A.zeroOutRowAndCol(row);
	}
}


vector<int> Solver::nodesWithFirstConditions()
{
	vector<int> result;
	int kol_el_po_x = m_xGrid.size() - 1;
	int kol_el_po_y = m_yGrid.size() - 1;

	for (auto el_na_kot_ost = 0; el_na_kot_ost < m_finitElements.size(); el_na_kot_ost += kol_el_po_x * kol_el_po_y)
	{
		for (auto j = el_na_kot_ost; j <= kol_el_po_x - 1 + el_na_kot_ost; j++)
		{
			result.push_back(m_finitElements[j].globalNumbers[0]);
			result.push_back(m_finitElements[j].globalNumbers[1]);
			result.push_back(m_finitElements[j].globalNumbers[4]);
			result.push_back(m_finitElements[j].globalNumbers[5]);
		}

		int el_s_kot_nach_posl_ryad = kol_el_po_x * kol_el_po_y - kol_el_po_x + el_na_kot_ost;
		for (auto k = el_na_kot_ost; k <= el_s_kot_nach_posl_ryad;)
		{
			result.push_back(m_finitElements[k].globalNumbers[0]);
			result.push_back(m_finitElements[k].globalNumbers[2]);
			result.push_back(m_finitElements[k].globalNumbers[4]);
			result.push_back(m_finitElements[k].globalNumbers[6]);

			int posl_el_v_ryadu = k + kol_el_po_x - 1;
			result.push_back(m_finitElements[posl_el_v_ryadu].globalNumbers[1]);
			result.push_back(m_finitElements[posl_el_v_ryadu].globalNumbers[3]);
			result.push_back(m_finitElements[posl_el_v_ryadu].globalNumbers[5]);
			result.push_back(m_finitElements[posl_el_v_ryadu].globalNumbers[7]);

			k = posl_el_v_ryadu + 1;
		}

		for (auto l = kol_el_po_x * kol_el_po_y - kol_el_po_x + el_na_kot_ost; l <= kol_el_po_x * kol_el_po_y + el_na_kot_ost - 1; l++)
		{
			result.push_back(m_finitElements[l].globalNumbers[2]);
			result.push_back(m_finitElements[l].globalNumbers[3]);
			result.push_back(m_finitElements[l].globalNumbers[6]);
			result.push_back(m_finitElements[l].globalNumbers[7]);
		}
	}

	std::sort(result.begin(), result.end());
	auto last = std::unique(result.begin(), result.end());
	result.erase(last, result.end());

	return result;
}

// пока что учет только точки А!!! СПРОСИТЬ
void Solver::assembleVector()
{
	vector<double> globB(m_xGrid.size() * m_yGrid.size() * m_zGrid.size(), 0.0);

	for (auto i = 0; i < m_sources.size(); i++)
	{
		auto globNums = m_finitElements[m_sources[i].elA].globalNumbers;
		double xp = m_xGrid[m_sources[i].pA]; double xp1 = m_xGrid[m_sources[i].pA + 1]; ;
		double hx = xp1 - xp;

		double ys = m_yGrid[m_sources[i].sA]; double ys1 = m_yGrid[m_sources[i].sA + 1]; ;
		double hy = ys1 - ys;

		globB[globNums[4]] += X_func(0, m_sources[i].A.x, xp1, hx) * Y_func(0, m_sources[i].A.y, ys1, hy);
		globB[globNums[5]] += X_func(1, m_sources[i].A.x, xp, hx) * Y_func(0, m_sources[i].A.y, ys1, hy);
		globB[globNums[6]] += X_func(0, m_sources[i].A.x, xp1, hx) * Y_func(1, m_sources[i].A.y, ys, hy);
		globB[globNums[7]] += X_func(1, m_sources[i].A.x, xp, hx) * Y_func(1, m_sources[i].A.y, ys, hy);
	}

	m_b.setV(globB);
}

void Solver::assembleMatrix()
{
	vector<int> ig = m_A.getIg();
	vector<int> jg = m_A.getJg();

	vector<double> di(m_xGrid.size() * m_yGrid.size() * m_zGrid.size(), 0.0);
	vector<double> ggl(jg.size(), 0.0);

	int nx = m_xGrid.size();
	int ny = m_yGrid.size();
	int nz = m_zGrid.size();

	int nEl = 0;
	for (auto r = 0; r < nz - 1; r++)
	{
		for (auto s = 0; s < ny - 1; s++)
		{
			for (auto p = 0; p < nx - 1; p++)
			{
				double hx = m_xGrid[p + 1] - m_xGrid[p];
				double hy = m_yGrid[s + 1] - m_yGrid[s];
				double hz = m_zGrid[r + 1] - m_zGrid[r];

				Matrix local = buildLocalG(hx, hy, hz, m_finitElements[nEl].sigma);
				auto glob_nums = m_finitElements[nEl].globalNumbers;

				for (auto i = 0; i < m_finitElements[0].globalNumbers.size(); i++)
				{
					for (auto j = 0; j < m_finitElements[0].globalNumbers.size(); j++)
					{
						int glob_i = glob_nums[i]; int glob_j = glob_nums[j];

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
	int N = m_xGrid.size() * m_yGrid.size() * m_zGrid.size(); // число глобальных б/ф (количество узлов)

	uint32_t bigNumberWhichIDontKnowButMoreIThinkAboutItMoreItSeemsToBeLessThanNTimes8 = N * 20;
	vector<vector<uint32_t>> list(2, std::vector<uint32_t>(bigNumberWhichIDontKnowButMoreIThinkAboutItMoreItSeemsToBeLessThanNTimes8, 0));
	std::vector<int> ig, jg;
	ig.reserve(N);
	jg.reserve(bigNumberWhichIDontKnowButMoreIThinkAboutItMoreItSeemsToBeLessThanNTimes8);
	vector<uint32_t> listbeg(N * N, 0);
	uint32_t listSize = 0;

	for (const auto elem : m_finitElements)
	{
		for (uint32_t i = 0; i < elem.globalNumbers.size(); i++)
		{
			auto k = elem.globalNumbers[i];
			for (uint32_t j = i + 1; j < elem.globalNumbers.size(); j++)
			{
				auto ind1 = k;
				auto ind2 = elem.globalNumbers[j];
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
	int nx = m_xGrid.size();
	int ny = m_yGrid.size();
	int nz = m_zGrid.size();

	m_finitElements.resize((nx - 1) * (ny - 1) * (nz - 1));

	int r = 0;
	int elN = 0;
	int numberOfElementFromWhichLastLayerStarts;

	for (r = 0; r < nz - 2; r++)
	{
		for (auto s = 0; s < ny - 1; s++)
		{
			for (auto p = 0; p < nx - 1; p++)
				m_finitElements[elN++] = GetNumbers(p, s, r, nx, ny);
		}
	}

	numberOfElementFromWhichLastLayerStarts = elN;

	for (auto s = 0; s < ny - 1; s++)
	{
		for (auto p = 0; p < nx - 1; p++)
			m_finitElements[elN++] = GetNumbers(p, s, r, nx, ny);
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
					is_found_A = isPointInside(m_xGrid[p], m_xGrid[p + 1], m_yGrid[s], m_yGrid[s + 1], m_sources[i].A);
					if (is_found_A)
					{
						m_sources[i].elA = whereToSearchPoint;
						m_sources[i].pA = p; m_sources[i].sA = s; m_sources[i].rA = r;
					}
				}

				if (!is_found_B)
				{
					is_found_B = isPointInside(m_xGrid[p], m_xGrid[p + 1], m_yGrid[s], m_yGrid[s + 1], m_sources[i].B);
					if (is_found_B)
					{
						m_sources[i].elB = whereToSearchPoint;
						m_sources[i].pB = p; m_sources[i].sB = s; m_sources[i].rB = r;
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

void Solver::buildDomain()
{
#ifdef _NORMAL_TESTS
	buildZGrid();
	m_xGrid = buildXGrid();
	m_yGrid = buildYGrid();
#endif

	defineGlobalNumbers();
}

vector<double> Solver::buildYGrid()
{
	vector<double> source_points;
	vector<double> grid;

	for (auto i = 0; i < m_sources.size(); i++)
	{
		source_points.push_back(m_sources[i].A.y);
		source_points.push_back(m_sources[i].B.y);
	}

	std::sort(source_points.begin(), source_points.end());
	auto last = std::unique(source_points.begin(), source_points.end());
	source_points.erase(last, source_points.end());

	vector<double> condensation_points;
	condensation_points.push_back(m_ref_plane[0].y);

	double hs = 15.; double hb = 50.; // smoll and big step
	for (auto i = 0; i < source_points.size(); i++)
	{
		condensation_points.push_back(source_points[i] - hb);
		condensation_points.push_back(source_points[i] + hb);
	}
	condensation_points.push_back(m_ref_plane[1].y);

	auto start = condensation_points[0]; auto end = condensation_points[1];
	grid.push_back(condensation_points[0]);

	auto build = [hb, hs, &grid](double& start, double& end, std::size_t i)
	{
		if (start <= end)
		{
			double h;
			i % 2 == 0 ? h = hb : h = hs;

			for (auto j = 1; j <= static_cast<int>((end - start) / h); j++)
			{
				grid.push_back(start + j * h);
			}
			if (grid.back() != end)
				grid.push_back(end);

			start = end;
		}
	};

	for (auto i = 0; i != condensation_points.size() - 2; end = condensation_points[i + 2], i++)
	{
		build(start, end, i);
	}
	build(start, end, condensation_points.size() - 2);

	/*std::ofstream fout;
	fout.open("y.txt");
	for (auto i = 0; i < grid.size(); i++)
	{
		fout << grid[i] << std::endl;
	}*/

	return grid;
}

vector<double> Solver::buildXGrid()
{
	vector<double> source_points;
	vector<double> grid;

	for (auto i = 0; i < m_sources.size(); i++)
	{
		source_points.push_back(m_sources[i].A.x);
		source_points.push_back(m_sources[i].B.x);
	}

	std::sort(source_points.begin(), source_points.end());
	auto last = std::unique(source_points.begin(), source_points.end());
	source_points.erase(last, source_points.end());

	vector<double> condensation_points;
	condensation_points.push_back(m_ref_plane[2].x);

	double hs = 15.; double hb = 50.; // smoll and big step
	for (auto i = 0; i < source_points.size(); i++)
	{
		condensation_points.push_back(source_points[i] - hb);
		condensation_points.push_back(source_points[i] + hb);
	}
	condensation_points.push_back(m_ref_plane[0].x);

	auto start = condensation_points[0]; auto end = condensation_points[1];
	grid.push_back(condensation_points[0]);

	auto build = [hb, hs, &grid](double& start, double& end, std::size_t i)
	{
		if (start <= end)
		{
			double h;
			i % 2 == 0 ? h = hb : h = hs;

			for (auto j = 1; j <= static_cast<int>((end - start) / h); j++)
			{
				grid.push_back(start + j * h);
			}
			if (grid.back() != end)
				grid.push_back(end);

			start = end;
		}
	};

	for (auto i = 0; i != condensation_points.size() - 2; end = condensation_points[i + 2], i++)
	{
		build(start, end, i);
	}
	build(start, end, condensation_points.size() - 2);

	/*std::ofstream fout;
	fout.open("x.txt");
	for (auto i = 0; i < grid.size(); i++)
	{
		fout << grid[i] << std::endl;
	}*/

	return grid;
}

void Solver::buildZGrid()
{
	vector<double> ref_point;
	ref_point.push_back(m_layers[0].z0);

	for (auto i = 0; i < m_layers.size(); i++)
	{
		ref_point.push_back(m_layers[i].z1);
	}

	for (auto i = 0; i < m_anomal.size(); i++)
	{
		ref_point.push_back(m_anomal[i].z0);
		ref_point.push_back(m_anomal[i].z1);
	}
	m_zGrid = ref_point;

	std::sort(ref_point.begin(), ref_point.end());
	auto last = std::unique(ref_point.begin(), ref_point.end());
	ref_point.erase(last, ref_point.end());

	double k = 1.2;
	double k_increase = 0.6 / (ref_point.size() - 2); // см 5ю страницу рис.
	double h = (ref_point[1] - ref_point[0]) / 10.;    // Пусть пока так, потом надо как-то нормировать первый шаг в завис от толщины.
	std::size_t pos = 1;

	for (auto i = 0; i < ref_point.size() - 1; i++)
	{
		m_zGrid.push_back(ref_point[i]);
		m_zGrid.push_back(ref_point[i] + h);
		std::size_t n = (std::size_t)floor((log((ref_point[i + 1]) / m_zGrid[pos]) / log(k)));
		
		for (auto j = 0; j < n; j++)
		{
			m_zGrid.push_back(m_zGrid[pos++] * k);
		}
		h = ref_point[i + 1] - m_zGrid.back();
		k += k_increase;
		pos += 2;
	}	
	m_zGrid.push_back(ref_point.back());

	/*std::ofstream fout;
	fout.open("z.txt");
	for (auto i = 0; i < m_zGrid.size(); i++)
	{
		fout << m_zGrid[i] << std::endl;
	}*/
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

		build(m_xGrid); build(m_yGrid); build(m_zGrid);

		fin.close();
	}
	else
	{
		std::cerr << "File domain.txt was not open." << std::endl;
		std::exit(1);
	}
#endif

#ifdef _NORMAL_TESTS
	std::ifstream fin(_path / "domain.txt");
	if (fin.is_open())
	{
		m_ref_plane.resize(4); 

		for (auto i = 0; i < m_ref_plane.size(); i++)
		{
			fin >> m_ref_plane[i].x >> m_ref_plane[i].y >> m_ref_plane[i].z;
		}

		std::size_t size; fin >> size;
		m_layers.resize(size);

		for (auto i = 0; i < size; i++)
		{
			int indicator; fin >> indicator;
			if (indicator > 0)
			{
				fin >> m_layers[i].z0 >> m_layers[i].z1 >> m_layers[i].sigma;
			}
		}
		fin.close();
	}
	else
	{
		std::cerr << "File domain.txt was not open." << std::endl; 
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
}
#endif