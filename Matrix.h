#ifndef _MATRIX_H
#define _MATRIX_H

#include <vector>
using std::vector;

class Matrix
{
public:
	vector<double>& operator[](const int index) { return M[index]; };

	Matrix(vector<vector<double>> _M)
	{
		M = _M;
		n = _M.size(); // ������
		m = _M[0].size(); // �������
	}


	Matrix(int _n, int _m)
	{
		n = _n; m = _m;
		M.resize(_n, vector<double>(_m));
	}

	void matrixByNum(Matrix R, double num)
	{
		for (int i = 0; i < R.getLine(); i++)
		{
			for (int j = 0; j < R.getColumn(); j++)
			{
				this->M[i][j] = R[i][j] * num;
			}
		}
	}

	vector<double> MultMatrixByVector(vector<double>& v)
	{
		vector<double> result(v.size(), 0.0);
		for (auto i = 0; i < M.size(); i++)
		{
			for (auto j = 0; j < M[0].size(); j++)
			{
				result[i] += M[i][j] * v[j];
			}
		}

		return result;
	}

	vector<vector<double>> getMatrix()
	{
		return M;
	}

	int getLine()
	{
		return n;
	}

	int getColumn()
	{
		return m;
	}

private:
	vector<vector<double>> M;
	int n; int m;
};

#endif