#ifndef _MATRIX_SPARSE_H
#define _MATRIX_SPARSE_H

using std::vector;

class MatrixSparse
{
public:
	MatrixSparse()
	{

	};

	void setDi(vector<double> _di)
	{
		di = _di;
		// ÌÛ Ë ıÂÌ¸... œ≈–≈ƒ≈À¿“‹
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

	//vector<double> getDi()
	//{
	//	return di;
	//}

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
	auto iend = jg[row + 1];
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


#endif