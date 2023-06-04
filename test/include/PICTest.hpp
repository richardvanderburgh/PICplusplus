#include <gtest/gtest.h>
#include <init.hpp>

// Define a test case
//TEST(GTestExample, AdditionTest) {
//    // Test assertion
//    
//    EXPECT_EQ(2 + 2, 4);
//}

//TEST(FFTTest, fftTest) {
//
//	const int ng = 128;
//	std::vector<double> rhoRE(ng, 0.0);
//	std::vector<double> rhoCO(ng, 0.0);
//	std::vector<double> rhokRE(ng, 0.0);
//	std::vector<double> rhokCO(ng, 0.0);
//
//	std::vector<double> phiRE(ng, 0.0);
//	std::vector<double> phiCO(ng, 0.0);
//
//	std::vector<double> space(ng, 0.0);
//
//	space = linspace(0, 2*M_PI, ng);
//
//	double val = sin(55);
//
//	complex complexRho[ng];
//	complex complexRhok[ng];
//	complex complexPhik[ng];
//	complex complexPhi[ng];
//
//	for (int i = 0; i < ng; i++) {
//		complexRho[i] = sin(space[i]);
//		rhoRE[i] = complexRho[i].re();
//		rhoCO[i] = complexRho[i].im();
//
//		std::cout << "rhoRE[" << i << "] = " << rhoRE[i] << std::endl;
//		std::cout << "rhoCO[" << i << "] = " << rhoCO[i] << std::endl;
//	}
//
//	matplot::scatter(space, rhoRE, 1);
//	matplot::scatter(space, rhoCO, 1);
//
//	CFFT::Forward(complexRho, complexRhok, ng);
//
//	for (int i = 0; i < ng; i++) {
//		//complexRho[i] = sin(space[i]);
//		rhokRE[i] = complexRhok[i].re();
//		rhokCO[i] = complexRhok[i].im();
//
//		std::cout << "rhokRE[" << i << "] = " << rhokRE[i] << std::endl;
//		std::cout << "rhokCO[" << i << "] = " << rhokCO[i] << std::endl;
//	}
//
//	matplot::scatter(space, rhokRE, 1);
//	matplot::scatter(space, rhokCO, 1);
//
//	CFFT::Inverse(complexRhok, complexPhi, ng);
//
//	for (int i = 0; i < ng; i++) {
//		//complexRho[i] = sin(space[i]);
//		phiRE[i] = complexPhi[i].re();
//		phiCO[i] = complexPhi[i].im();
//
//		std::cout << "phiRE[" << i << "] = " << phiRE[i] << std::endl;
//		std::cout << "phiCO[" << i << "] = " << phiCO[i] << std::endl;
//	}
//
//	matplot::scatter(space, phiRE, 1);
//	matplot::scatter(space, phiCO, 1);
//}
//
TEST(PICTest, initTest)
{
	int N = 5;
	int nt = 5;
	double dt = 0.01;
	Init init;
	bool success = init.initialize(N, nt, dt);

	int ng = 32;

	std::vector<double> expectedElectricField{
		4.2697e-05,
	0.40062,
	0.6913,
	0.1527,
	-0.58858,
	-0.53469,
	-0.157,
	0.23765,
	0.58721,
	0.44277,
	-0.29788,
	-0.63965,
	-0.31838,
	0.076411,
	0.48233,
	0.73623,
	-0.0057763,
	-0.7427,
	-0.48293,
	-0.076372,
	0.31842,
	0.63924,
	0.2943,
	-0.44633,
	-0.58765,
	-0.23755,
	0.15686,
	0.53581,
	0.59793,
	-0.14336,
	-0.69018,
	-0.40077,
	4.2697e-05
	};

	bool fieldsMatch = true;
	for (int i = 0; i < ng; i++) {

		if (init.mPicData.frames[1].electricField[i] != expectedElectricField[i]) {
			fieldsMatch = false;
		}

	}


	EXPECT_TRUE(fieldsMatch);
}