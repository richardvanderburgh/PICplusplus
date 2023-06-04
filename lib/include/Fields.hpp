void fields(std::vector<double>& rho, double L, int iw, double dx, std::vector<std::vector<double>>& E, int t, const int constNg, std::vector<double>& a, double ael)
{
	// FIELDS - E field solver
	const int ng = 32;
	std::vector<double> phi(ng + 1, 0.0);
	double l = L;

	// Transform charge density.
	rho[0] += rho[ng];
	rho[ng] = rho[0];

	complex complexRho[ng];
	complex complexRhok[ng];
	complex complexPhik[ng];

	for (int i = 0; i < ng; i++) { complexRho[i] = rho[i]; }
	CFFT::Forward(complexRho, complexRhok, ng);

	double tol = 1;
	for (int k = 0; k < ng; k++) {
		int ii;
		if (k == 0) {
			ii = tol;
		}
		else if (k <= ng / 2 && k >= 0) {
			ii = k;
		}
		else if (k > ng / 2) {
			ii = k - ng;
		}

		complexPhik[k] = complexRhok[k] / -std::pow(2.0 * M_PI * ii / L, 2.0);
		//esestot[t + 1] = (esestot[t + 1] - (phik[k] * rhok[k])) / 2;
	}

	complex complexPhi[ng];
	CFFT::Inverse(complexPhik, complexPhi, ng);

	// Perform FFT using FFTW
	//std::vector<std::complex<double>> output(ng);
	//fftw_plan plan = fftw_plan_dft_r2c_1d(ng, rho.data(), reinterpret_cast<fftw_complex*>(output.data()), FFTW_ESTIMATE);
	//fftw_execute(plan);
	//fftw_destroy_plan(plan);


	for (int i = 0; i < ng; i++) {
		phi[i] = -complexPhi[i].re();
	}

	phi[ng] = phi[0];

	//Update esem using phi and rhok
   //for (int k = 0; k < mplot + 1; k++) {
   //	esem(t + 1, k) -= phik(k) * rhok(k) / 2.0;
   //}
	if (iw == 1 || iw == 2)
	{
		for (int j = 1; j < ng; j++) {
			E[t][j] = (phi[j - 1] - phi[j + 1]) / (2.0 * dx);
		}
		E[t][0] = (phi[ng - 1] - phi[1]) / (2.0 * dx);
		E[t][ng] = E[t][0];
	}
	else if (iw == 3)
	{
		double dxi = 1.0 / dx;
		for (int j = 0; j < ng; j++) {
			E[t][j] = (phi[j] - phi[j + 1]) * dxi;
		}
		E[t][ng] = E[t][0];
	}

	ael = 1;
	for (int i = 0; i <= ng; i++) {
		a[i] = E[t][i];
	}
}