#pragma once

#include <boost/algorithm/algorithm.hpp>
#include <fft.hpp>
#include <fftw3.h>
#include <matplot/matplot.h>
#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <corecrt_math_defines.h>
#include <iostream>
#include <thread>
#include <vector>

#include "Accel.hpp"
#include "Fields.hpp"
#include "SetRho.hpp"

class Init {
public:

	struct SpeciesData {

		int speciesIndex;
		std::vector<double> positions;
		std::vector<double> velocities;
	};

	//struct ParticalData {
	//	double position;
	//	double velocity;
	//	std::string species;
	//	int id;
	//};

	struct FrameData {
		//std::vector <ParticalData> speciesData;
		std::vector <SpeciesData> speciesData;

	};

	struct PicData {
		std::vector<FrameData> frameData;
	};

	void initialize(int N1, int nt) {
		//int ng, int nt, double L, double dt, std::vector<int> N, int nsp, std::vector<double> qm, std::vector<double> wp, std::vector<double> wc, int mplot
		// 
		// FIRST EE - Set initial input values
		std::string example = "Electron - Electron Stream";
		// Input Variables
		double L = 6.28318530717958; // Physical length of system in meters
		int nsp = 2; // Number of particle species
		//int nt = 5;//600; // Number of time steps
		double dt = .1; // Time step in seconds
		int epsi = 1; // 1 over epsilon naught(F / m) Epsilon normalized to 1
		const int ng = 32; // Number of spatial grid points - only change the power of 2
		int iw = 2; // 1 for nearest grid point, 2 for linear
		int a1 = 0; // Smoothing factor
		int a2 = 0;
		int e0 = 0; // Add uniform electric field e0* cos(w0 * time)
		int w0 = 0;

		// Plotting Intervals
		int irho = 0; // nt / 10 + 4;
		int iphi = 0; // nt / 10 + 4;
		int iE = 0; // nt / 10 + 4;
		int mplot = 3;
		int ixvx = 1; // 60
		int ivxvy = 0;
		int ifvx = 0;// nt / 10 + 4;
		int nplot = 30; // ? ?

		// //Species Input Variables
		std::vector<int> N = { N1, N1 }; // Number of simulation particles
		std::array<int, 2> wp = { 1, 1 }; // Plasma Frequency
		std::array<int, 2> wc = { 0, 0 }; // Cyclotron frequency
		std::array<int, 2> qm = { -1, -1 }; // q / m charge to mass ratio(C / kg)
		std::array<int, 2> vt1 = { 0, 0 }; // RMS thermal velocity for random velocities
		std::array<int, 2> vt2 = { 0, 0 }; // RMS thermal velocity for ordered velocities
		std::array<int, 2> nv2 = { 0, 0 };
		std::array<int, 2> nlg = { 1, 1 }; // Number of loading groups
		std::array<int, 2> v0 = { 1, -1 }; // Drift velocity
		std::array<int, 2> pch = { 0, 0 }; // species pitch angle
		int distribution = 1; // Distribution 0 = Cold 1 = Two - Stream

		// Perturbation
		std::array<int, 2> mode = { 1, 1 };
		std::array<double, 2> x1 = { 0.001, 0.001 };
		std::array<int, 2> v1 = { 0, 0 };
		std::array<int, 2> thetax = { 0, 0 };
		std::array<int, 2> thetav = { 0, 0 };

		//////////////////////////////////////////////////////////
		// INIT - Initial values for each species

		int cv = 8;
		int angle = 0; // Magnetic field angle

		double step_size = L / ng;
		std::vector<double> gridx;

		for (double x = 0.0; x <= L; x += step_size) {
			gridx.push_back(x);
		}

		double dx = gridx[1]; // spatial grid bin width
		double dtdx = dt / dx;
		int nxp1 = ng + 1;
		int nxp2 = ng + 2;
		int slx = ng;
		int npt = 0;

		for (int i = 0; i < nsp; i++) {
			npt += N[i];
		}

		std::vector<int> X1;
		std::vector<int> X2;
		std::vector<int> X3;

		for (int i = 1; i <= ng; i++) {
			X1.push_back(i);
		}

		for (int i = 2; i <= ng + 1; i++) {
			X2.push_back(i);
		}

		for (int i = 3; i <= ng + 2; i++) {
			X3.push_back(i);
		}

		int cs = cv * cv;

		std::vector<double> q(nsp, 0.0);
		std::vector<double> mass(nsp, 0.0);

		for (int i = 0; i < nsp; i++) {
			q[i] = ng / N[i] * (wp[i] * wp[i]) / qm[i];
		}

		for (int i = 0; i < nsp; i++) {
			mass[i] = q[i] / qm[i];
		}

		double theta = M_PI / 180 * angle;
		double costh = cos(theta);
		double sinth = sin(theta);
		double b0 = wc[1] / qm[1];
		double bx0 = b0 * costh;
		double by0 = b0 * sinth;

		std::vector<std::vector<double>> E(ng + 1, std::vector<double>(nt + 1, 0.0));

		// Initialize all elements of the array to zero
		for (int i = 0; i < ng + 1; i++) {
			for (int j = 0; j < nt + 1; j++) {
				E[i][j] = 0.0;
			}
		}
		// --Field Initialization --
		std::vector<double> ex(nxp2, 0.0);
		std::vector<double> ey(nxp2, 0.0);
		std::vector<double> ez(nxp2, 0.0);
		std::vector<double> by(nxp2, by0);
		std::vector<double> bz(nxp2, 0.0);
		std::vector<double> ajx(nxp2, 0.0);
		std::vector<double> ajy(nxp2, 0.0);
		std::vector<double> ajz(nxp2, 0.0);

		int maxN = *std::max_element(N.begin(), N.end());
		std::vector<std::vector<double>> x(nsp, std::vector<double>(maxN, 0.0));
		std::vector<std::vector<double>> vx(nsp, std::vector<double>(maxN, 0.0));
		std::vector<std::vector<double>> vy(nsp, std::vector<double>(maxN, 0.0));
		std::vector<std::vector<double>> vz(nsp, std::vector<double>(maxN, 0.0));

		//// Energies
		// Total time for plot
		std::vector<std::vector<double>> esem(nt + 1, std::vector<double>(mplot + 1, 0.0));
		std::vector<std::vector<double>> ke(nt + 1, std::vector<double>(nsp, 0.0));
		std::vector<std::vector<double>> p(nt + 1, std::vector<double>(nsp, 0.0));
		std::vector<std::vector<double>> de(nt + 1, std::vector<double>(nsp, 0.0));
		std::vector<std::vector<double>> therme(nt + 1, std::vector<double>(nsp, 0.0));
		std::vector<std::vector<double>> v_old(maxN, std::vector<double>(nsp));
		std::vector<std::vector<double>> vy_old(maxN, std::vector<double>(nsp));
		std::vector<double> ESE(nt + 1, 0.0);
		std::vector<double> esestot(nt + 1, 0.0);
		std::vector<double> te(nt + 1, 0.0);

		std::vector<double> ddx(nsp);
		std::vector<double> m(nsp);
		std::vector<double> T(nsp);

		for (int species = 0; species < nsp; species++)
		{
			q[species] = L * wp[species] * wp[species] / (epsi * N[species] * qm[species]);
			m[species] = q[species] / qm[species];
			T[species] = tan(-wc[species] * dt / 2.);
			ddx[species] = L / N[species];

			double ngr = N[species] / nlg[species];
			double nm = N[species] * m[species];
			double lg = L / nlg[species];

			// Set evenly spaced charge distribution distribution
			// remember that ddx is the width of the charge cloud
			for (int I = 1; I < N[species] + 1; I++) {
				x[species][I - 1] = (I - 0.5) * ddx[species];
			}

			for (int I = 0; I < N[species]; I++) {
				vx[species][I] = v0[species];
			}

			// Load ordered velocities in vx ("quiet start", or at least subdued).
			// Is set up for Maxwellian*v**nv2, but can do any smooth distribution.
			// Hereafter, ngr is preferably a power of 2.
			// First store indefinite integral of distribution function in x array.
			// Use midpoint rule -simple and quite accurate.

			int i1 = 0;

			if (vt2[species] != 0) {
				int vmax = 5 * vt2[species];
				double dv = 2 * vmax / (N[species] - 1);
				double vvnv2 = 1;
				x[species][1] = 0;

				for (int ith = 2; ith <= N[species]; ith++) {
					double vv = ((ith - 1.5) * dv - vmax) / vt2[species];

					if (nv2[species] != 0) {
						vvnv2 = pow(vv, nv2[species]);
					}

					double fv = vvnv2 * exp(-0.5 * pow(vv, 2));
					int i1 = ith - 1 + 1;
					x[species][i1] = x[species][i1 - 1] + std::max(fv, 0.0);
				}

				// For evenly spaced(half - integer multiples) values of the integral,
				// find corresponding velocity by inverse linear interpolation.

				double df = x[species][i1] / ngr;
				int i1 = 1;
				int j = 1;

				for (int ith = 0; ith < ngr; ith++) {
					double fv = (ith + 1 - 0.5) * df;

					while (fv >= x[species][j + 1]) {
						j++;
					}

					double vv = dv * (j - 1 - 1 + (fv - x[species][j - 1]) / (x[species][j] - x[species][j - 1])) - vmax;
					vx[species][i1] = vx[species][i1] + vv;
					i1++;
				}

				// For ordered velocities, scramble positions to reduce correlations.
				// Scrambling is done by bit - reversed counter - compare sorter in cpft.
				// xs = .000, .100, .010, .110, .001, .101, .011, .111, .0001.. (binary fractions)

				double xs = 0;

				for (int ith = 1; ith <= ngr; ith++) {
					int i1 = ith - 1 + 1;
					x[i1][species] = xs * lg + 0.5 * ddx[species];
					double xsi = 1.0;

					while (xs >= 0) {
						xsi = 0.5 * xsi;
						xs = xs - xsi;
					}

					xs = xs + 2.0 * xsi;
					i1 = ngr + 1 - 1;
				}
			}

			if (wc[species] != 0) {
				for (int ith = 1; ith <= ngr; ++ith) {
					int i1 = ith - 1 + 1;
					double vv = vx[species][i1];
					double theta = 2 * M_PI * x[species][i1] / lg;
					vx[species][i1] = vv * cos(theta);
					vy[species][i1] = vv * sin(theta);
				}
			}

			if (nlg[species] != 1) {
				int j = ngr + 1;
				double xs = 0;
				for (int ith = j; ith <= N[species]; ++ith) {
					xs = xs + lg;
					for (int jj = 1; jj <= ngr; ++jj) {
						int i1 = jj - 1 + 1;
						int i2 = i1 + ith - 1;
						x[species][i2] = x[species][i1] + xs;
						vx[species][i2] = vx[species][i1];
						if (wc[species] != 0) {
							vy[species][i2] = vy[species][i1];
						}
					}
				}
			}

			// Add random Maxwellian.
			if (vt1[species] != 0) {
				for (int I = 0; I < N[species]; ++I) {
					double rm = 0;
					for (int ith = 0; ith < 12; ++ith) {
						rm += (double)rand() / RAND_MAX;
					}
					rm -= 6;

					if (wc[species] != 0) {
						vy[species][I] += vt1[species] * rm;
					}
					vx[species][I] += vt1[species] * rm;
				}
			}

			for (int a = 0; a < N[species]; ++a) {
				double theta = 2 * M_PI * mode[species] * x[species][a] / L;
				x[species][a] = x[species][a] + x1[species] * cos(theta + thetax[species]);
				vx[species][a] = vx[species][a] + v1[species] * sin(theta + thetav[species]);

				if (x[species][a] >= L) {
					x[species][a] = x[species][a] - L;
				}
				if (x[species][a] < 0) {
					x[species][a] = x[species][a] + L;
				}
			}
		}
		std::vector<double> N_grid = linspace(0, L, 32);
		std::vector<double> gridt = linspace(0, (nt + 1) * dt, nt + 1);
		std::vector<double> xhist(nt, 0.0);
		std::vector<double> xvxt(nt, 0.0);
		std::vector<double> rhot(nt, 0.0);
		std::vector<double> phit(nt, 0.0);

		// End of INIT
		// 
		// SETRHO - c  Converts position to computer normalization and accumulates charge density.
		std::vector<double> qdx{ q[0] / dx, q[1] / dx };
		std::vector<double> rho(ng + 1, 0.0); // charge density at the spatial grid points
		std::vector<std::vector<double>> rhos(nsp, std::vector<double>(ng + 1, 0.0));
		std::vector<std::vector<double>> rho0(nsp, std::vector<double>(ng + 1, 0.0));
		std::vector<std::vector<double>> qjdata(nsp, std::vector<double>(ng + 1, 0.0));
		std::vector<std::vector<double>> qjp1data(nsp, std::vector<double>(ng + 1, 0.0));
		std::vector<std::vector<double>> drhodata(nsp, std::vector<double>(ng + 1, 0.0));

		for (int species = 0; species < nsp; species++) {

			retRho(species, ng, dx, N, qdx, rho, x, rho0, rhos, iw);
			// Velocities
			for (int K = 0; K < N[species]; ++K) {
				vx[species][K] *= dtdx;
			}
		}

		// Acceleration
		std::vector<double> a(ng + 1);

		int t = 0;
		double ael = 1;
		fields(rho, L, iw, dx, E, t, ng, a, ael);
		accel(nsp, dx, dt, t, q, m, ael, a, ng, N, x, vx);

		//BEGIN TIME LOOP //////////////////////////////////////////////////////////////////
		// Create a figure and set the x and y limits
		//matplot::figure("Animation");

		//PicData picData;
		nlohmann::json frames;

		for (int t = 1; t <= nt; t++) {

			accel(nsp, dx, dt, t, q, m, ael, a, ng, N, x, vx);
			move(nsp, rho, rho0, qdx, N, x, vx, ng);
			fields(rho, L, iw, dx, E, t, ng, a, ael);

			FrameData frameData;

			// Store the data for the current frame (for all species)
			for (int sp = 0; sp < nsp; sp++) {
				SpeciesData speciesData;
				speciesData.speciesIndex = sp;
				speciesData.positions = x[sp];
				speciesData.velocities = vx[sp];
				frameData.speciesData.push_back(speciesData);
			}

			nlohmann::json speciesArray;
			for (const auto& species : frameData.speciesData) {
				nlohmann::json speciesObject;
				speciesObject["speciesIndex"] = species.speciesIndex;
				speciesObject["positions"] = species.positions;
				speciesObject["velocities"] = species.velocities;
				speciesArray.push_back(speciesObject);
			}

			nlohmann::json frame;
			//nlohmann::json particle;
			for (const auto& species : frameData.speciesData) {
				nlohmann::json speciesObject;
				speciesObject["speciesIndex"] = species.speciesIndex;
				speciesObject["positions"] = species.positions;
				speciesObject["velocities"] = species.velocities;
				speciesArray.push_back(speciesObject);
			}

			frame["timeStep"] = t;  // Store the time step index
			frame["speciesData"] = speciesArray;

			//nlohmann::json frame;
			//frame["positions"] = x[0];
			//frame["velocities"] = vx[0];
			frames.push_back(frame);

			//picData.frameData.push_back(frameData);

			//// Create a scatter plot 
			//matplot::scatter(x[0], vx[0], 1);
			//matplot::hold(true);
			//matplot::scatter(x[1], vx[1], 1);

			//matplot::show();
			//matplot::hold(false);
			//std::cout << "t = " << t << std::endl;
		}
		std::cout << frames.dump() << std::endl;
	}
	
	std::vector<double> linspace(double start, double end, int num_points) {
		std::vector<double> result(num_points);
		double step = (end - start) / (num_points - 1);
		for (int i = 0; i < num_points; ++i) {
			result[i] = start + i * step;
		}
		return result;
	}
};