#include <filesystem>
#include "Solver.h"

int main()
{
	using ::std::filesystem::path;
	path dir_to_file = "domain";
	Solver s(dir_to_file);

	return 0;
}