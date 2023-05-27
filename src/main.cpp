//#include <iostream>
//#include "First_EE.hpp"
#include "init.hpp"

//#include <boost/algorithm/algorithm.hpp>
//#include <nlohmann/json.hpp>


int main(int argc, char* argv[]) {

	//inputFile myInput;

	//std::cout << myInput;

	//for (int i = 0; i < argc; ++i)
	//	std::cout << argv[i] << "\n";

	Init init;

	int N1 = atoi(argv[1]);
	int nt = atoi(argv[2]);


	init.initialize(N1);

	return 0;
}