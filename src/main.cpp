//#include <iostream>
//#include "First_EE.hpp"
#include "init.hpp"

#include <gtest/gtest.h>

//#include <boost/algorithm/algorithm.hpp>
//#include <nlohmann/json.hpp>


int main(int argc, char* argv[]) {


	Init init;

	int N1 = atoi(argv[1]);
	int nt = atoi(argv[2]);
	double dt = std::stod(argv[3]);

	init.initialize(N1, nt, dt);

	return 0;
}