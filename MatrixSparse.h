#ifndef _MATRIX_SPARSE_H
#define _MATRIX_SPARSE_H

using std::vector;
#include "lightweight.h"

class MatrixSparse
{
public:
	MatrixSparse()
	{

	};

	vector<double> MultiplyByVector(const std::vector<double>& vector);
	vector<vector<double>> getDenseMatrix();

	std::vector<double> operator *(const std::vector<double>& vector) { return MultiplyByVector(vector); }


	void setDi(vector<double> _di)
	{
		di = _di;
		m_N = di.size();
	}

	void setIg(vector<int> _ig)
	{
		ig = _ig;
	}

	void setJg(vector<int> _jg)
	{
		jg = _jg;
	}

	void setGgu(vector<double> _ggu)
	{
		ggu = _ggu;
	}

	void setGgl(vector<double> _ggl)
	{
		ggl = _ggl;
	}

	vector<int> getIg()
	{
		return ig;
	}

	vector<int> getJg()
	{
		return jg;
	}

	vector<double> getGgl()
	{
		return ggl;
	}

	vector<double> getGgu()
	{
		return ggu;
	}


	vector<double>& getDi()
	{
		return di;
	}

	void zeroOutRowAndCol(int row);

private:
	vector<int> ig, jg;
	vector<double> ggl, ggu;
	vector<double> di;

	int m_N;

};

void MatrixSparse::zeroOutRowAndCol(int row)
{
	auto ibeg = ig[row];
	auto iend = ig[row + 1];
	for (auto i = ibeg; i < iend; i++)
	{
		ggl[i] = 0.0;
	}
	for (auto j = row + 1; j < m_N; j++)
	{
		auto jbeg = ig[j];
		auto jend = ig[j + 1];
		auto index = binarySearch(jg, row, jbeg, jend - 1);
		if (index != -1)
		{
			ggl[index] = 0.0;
		}
	}
}

std::vector<double> MatrixSparse::MultiplyByVector(const std::vector<double>& vector)
{
	std::vector<double> result(vector.size(), 0.0);
	for (auto i = 0u; i < vector.size(); i++)
	{
		result[i] = di[i] * vector[i];
	}

	for (auto i = 0u; i < vector.size(); i++)
	{
		for (auto j = ig[i]; j < ig[i + 1]; j++)
		{
			result[i] += ggl[j] * vector[jg[j]];
			result[jg[j]] += ggl[j] * vector[i];
		}
	}
	return result;
}

vector<vector<double>> MatrixSparse::getDenseMatrix()
{
	vector<vector<double>> res(ig.size() - 1, vector<double>(ig.size() - 1));

	for (auto i = 0; i < ig.size() - 1; i++)
	{
		for (auto j = 0; j < ig.size() - 1; j++)
		{
			if (i > j)
			{
					int ind;
					int l;
					bool isFound = false;

					for (l = ig[i]; l < ig[i + 1]; l++)
					{
						if (jg[l] == j)
						{
							isFound = true;
							ind = l;
							break;
						}
					}
					if (isFound)
					{
						res[i][j] = ggl[l];
						res[j][i] = ggl[l];
					}
			}
			else if (i == j)
			{
				res[i][i] = di[i];
			}
		}
	}
	return res;
}

#endif