#include <stdio.h>
#include<math.h>

#define I(i,j,s) ((i)*(s)+(j))
double coeff[] = {
	 1.00E+03, 1.76E+05, 2.79E+07, 8.47E+10, 3.20E+11,
	 9.46E-03,-1.34E-01, 1.10E-01,-6.65E-02,-1.40E-01,
	 4.27E-02, 0.00E+00,-1.91E-01,-9.74E-02, 5.72E-02,
	 1.80E-03,-1.25E-02,-2.08E-04, 1.23E-03, 1.13E-03,
	-7.99E-06, 0.00E+00, 4.84E-06, 4.50E-06,-2.08E-06,
	 3.37E-03,-1.93E-02, 2.11E-03, 5.36E-03, 3.42E-03,
	2.26E-02 ,-6.00E-02, 2.21E-04, 2.56E-03,-1.16E-02,
	3.79E-02, -2.00E-02,-1.62E-01,-9.82E-02, 5.52E-02,
	-1.92E-02, 5.88E-02,-9.22E-02, 1.01E-01,-1.01E-02,
	-9.24E-03, 0.00E+00,-8.30E-03, 6.26E-03,-4.39E-02,
	-2.11E+02, 9.23E+01, 2.14E+02, 1.16E+01, 1.96E+02,
	 1.03E-02, 0.00E+00, 2.35E-01, 1.76E-01, 3.25E-02,
	 2.89E-02, 0.00E+00,-7.91E-02,-7.13E-02, 5.65E-02,
	-7.63E+01, 0.00E+00, 1.10E+02, 1.06E+02, 8.82E+01,
	-1.85E-01, 3.30E-01,-1.27E+00, 3.33E-01, 2.88E-01,
	-2.03E-02, 1.05E-01,-1.85E+00,-1.15E-01,-3.14E-02,
	 1.45E-02, 0.00E+00, 6.66E-02,-4.24E-03, 0.00E+00,
	-3.63E+00,-1.48E+01,-1.87E+02,-6.75E-01,-2.00E+02,
	-2.89E-02,-9.07E-02,-4.22E-02,-4.17E-02, 6.07E-02,
	-1.74E+02,-7.20E+01,-2.17E+02, 1.34E+02, 5.12E+01,
	-1.08E-01, 2.09E-01,-1.28E-01,-6.59E-02,-4.61E-02,
	-2.00E-03, 2.83E-02,-6.18E-03,-1.93E-02,-9.70E-03,
	 3.40E-03, 0.00E+00,-1.75E-02,-8.59E-03, 3.49E-03,
	-1.61E-02, 8.57E-02,-4.36E-02, 8.64E-02, 0.00E+00,
	-9.61E-03,-2.48E-02,-5.29E-02, 8.57E-02, 0.00E+00,
	-1.05E-01, 3.83E-01, 3.12E-01, 4.55E-02,-7.33E-02,
	 4.58E-03, 2.94E-02,-2.14E-02,-3.87E-02, 2.12E-02,
	 4.58E-03, 0.00E+00,-2.45E-02,-1.29E-02, 7.03E-03,
	 1.78E-02,-3.97E-03,-3.67E-03,-9.18E-02, 0.00E+00,
	-4.23E-03, 4.36E-02,-8.83E-02, 4.26E-02, 0.00E+00,
	-3.57E-03, 0.00E+00, 3.40E-02, 4.34E-02,-6.92E-03,
	-3.61E-04, 0.00E+00, 5.13E-03, 5.69E-03, 0.00E+00,
	 1.05E-02, 0.00E+00,-1.47E-02, 6.65E-03,-7.29E-03,
	 4.57E-03, 0.00E+00, 8.56E-03,-1.02E-02, 0.00E+00,
	-2.18E-05, 0.00E+00, 1.96E-03, 2.71E-03,-2.92E-03,
	 1.51E-03, 0.00E+00,-4.50E-03,-4.71E-03, 2.71E-03,
	 2.17E-03, 0.00E+00,-1.02E-01,-1.46E-02,-5.99E-03,
	-1.14E-06, 0.00E+00,-1.69E-05,-6.41E-06, 0.00E+00,
	 6.58E-05, 0.00E+00,-1.43E-03, 1.52E-03, 8.21E-03
};

double harmonic_expansion(int i, double *P, double f, double f_mean, double k, double t, double d)
{
	double omega, Omega, F1, M,Q;
	double G, AN1, AN2, SAN1, SAN2, D, SD, TD, beta;
	double pi = 3.1415926535897932;
		omega = 2 * pi / 24;
	Omega = 2 * pi / 365;

	F1 = coeff[I(3, i,5)] * (f - f_mean) + coeff[I(4, i, 5)] * pow((f - f_mean), 2) +coeff[I(5, i, 5)] * (f_mean - 150) + coeff[I(37, i, 5)] * pow((f_mean - 150) , 2);
	M = (coeff[I(6,i, 5)] + coeff[I(7,i, 5)]*P[I(0,2,6)])*k;

	if (i==0) { //isT
		M = M + coeff[I(38,i, 5)]*exp(k);
	}
	else {
		M = M + pow(coeff[I(38,i, 5)]*k ,2);
	}


	Q = coeff[I(1,i, 5)]*P[I(0, 2, 6)] + coeff[I(2,i, 5)]*P[I(0, 4, 6)] + coeff[I(36,i, 5)]*P[I(0, 1, 6)];

	AN1 = (coeff[I(8,i, 5)] + coeff[I(9,i, 5)]*P[I(0, 2, 6)])*cos(Omega*(d - coeff[I(10,i, 5)]));
	AN2 = (coeff[I(14,i, 5)]*P[I(0, 1, 6)] + coeff[I(15,i, 5)]* P[I(0,3, 6)] + coeff[I(16,i, 5)]*P[I(0, 5, 6)])*cos(Omega*(d - coeff[I(17,i, 5)]));

	SAN1 = (coeff[I(11,i, 5)] + coeff[I(12,i, 5)]* P[I(0, 2, 6)])*cos(2 * Omega*(d - coeff[I(13,i, 5)]));
	SAN2 = coeff[I(18,i, 5)]* P[I(0, 1, 6)] *cos(2 * Omega*(d - coeff[I(19,i, 5)]));

	D = (coeff[I(20,i, 5)]*P[I(1, 1, 6)] + coeff[I(21,i, 5)]*P[I(1, 3, 6)] + coeff[I(22,i, 5)]*P[I(1, 5, 6)] +
		(coeff[I(23,i, 5)]*P[I(1, 1, 6)] + coeff[I(24,i, 5)]*P[I(1, 2, 6)])*cos(Omega*(d - coeff[I(17,i, 5)])))*cos(omega*t) +
		(coeff[I(25,i, 5)]*P[I(1, 1, 6)] + coeff[I(26,i, 5)]*P[I(1, 3, 6)] + coeff[I(27,i, 5)]*P[I(1, 5, 6)] +
		(coeff[I(28,i, 5)]*P[I(1, 1, 6)] + coeff[I(29,i, 5)]*P[I(1, 2, 6)])*cos(Omega*(d - coeff[I(17,i, 5)])))*sin(omega*t);

	SD = (coeff[I(30, i, 5)] *P[I(2, 2, 6)] + coeff[I(31, i, 5)] *P[I(2, 3, 6)] *cos(Omega*(d - coeff[I(17, i, 5)])))*cos(2 * omega*t) +
		(coeff[I(32, i, 5)] *P[I(2, 2, 6)] + coeff[I(33, i, 5)] *P[I(2, 3, 6)] *cos(Omega*(d - coeff[I(17, i, 5)])))*sin(2 * omega*t);

	TD = coeff[I(34, i, 5)] *P[I(3, 3, 6)] *cos(3 * omega*t) + coeff[I(35, i, 5)] *P[I(3, 3, 6)] *sin(3 * omega*t);

	beta = 1; //~isH
	if (i<1||i>2) {beta = beta + F1;}

	G = 1 + F1 + M + Q + beta*(AN1 + AN2 + SAN1 + SAN2 + D + SD + TD);
	return G;
 }

double  calculate_t_inf_and_initial_concentrations(
	double *P, double f, double f_mean, double k, double t, double d,double *cH,double *cHe,
	double *cO,double *cN2)
{
	double t_inf;
	//double coeffsT, coeffsH, coeffsHe, coeffsO, coeffsN2;
	double GT, GH, GHe, GO, GN2;
	//coeffsT = coeffs(:, 1);
	//coeffsH = coeffs(:, 2);
	//coeffsHe = coeffs(:, 3);
	//coeffsO = coeffs(:, 4);
	//coeffsN2 = coeffs(:, 5);

	GT = harmonic_expansion(0, P, f, f_mean, k, t, d);
	GH = harmonic_expansion(1, P, f, f_mean, k, t, d);
	GHe = harmonic_expansion(2, P, f, f_mean, k, t, d);
	GO = harmonic_expansion(3, P, f, f_mean, k, t, d);
	GN2 = harmonic_expansion(4, P, f, f_mean, k, t, d);

	t_inf = coeff[I(0,0,5)]*GT;
	*cH = coeff[I(0,1, 5)]*exp(GH - 1);
	*cHe = coeff[I(0,2, 5)]*exp(GHe - 1);
	*cO = coeff[I(0,3, 5)]*exp(GO - 1);
	*cN2 = coeff[I(0,4, 5)]*exp(GN2 - 1);
	return t_inf;

}
//function[t_inf, T, c_H, c_He, c_O, c_N2] = dtm94(z, coeffs, P, f, f_mean, k, t, d)
//implements the DTM 94 model for high altiude density as desribed by
//Berger et al.in "Improvement of the empirical thermospheric model DTM:
//DTM94 - a comparative review of various temporal variations
//and prospects in space geodesy applications".

//Given:
//-an altiude, z
//-a set of coefficients, coeffs
//-a vector of Legendre polynomials evaluated at sin(latitude), P
//-the daily solar flux, f
//-the average solar flux, f_mean
//-the current geomagnetic index, k
//-the local solar time, t
//-the day number, d
//this model returns :
//-the thermopause temperature, t_inf
//-the temperature, T
//-the concentration of Hydrogen, c_H
//-the concentration of Helium, c_He
//-the concentration of Oxygen, c_O
//-the concentration of molecular Nitrogen, c_N2
//The altitude z must be a scalar given in km.
//The coefficients must be given in a 36x4 matrix in the order of
//temperature, He, O, N2.
//P must be a matrix of size 6x6, such that P(i, j) = the associated
//legendre polynomial of degree i - 1 and order j - 1 evaluated at
//sin(latitude).
//f and f_mean must be measured in 10 ^ -22 W m^-2 Hz^-1.
//k must be an integer between 0 and 9.
//t must be measured in hours, between 0 and 24.
//d is the day number in the year, between 0 and 365.
//The temperatures returned are given in degrees Kelvin, while the
//concentrations are given in particles / cm^3.
////dtm94 Drag Temperature Model for high altitude atmospheric density
extern double dtm94(double z,
	double *P, double f, double f_mean, double k, double t, double d,
	double *T,double *c_H,double *c_He,double *c_O, double *c_N2)
{
	double t_inf=0.0;
	double g, xi, cH_0, cHe_0, cO_0, cN2_0;
	double sigma, gammaH, gammaHe, gammaO, gammaN2;
	double alphaH, alphaHe, fH, fHe, fO, fN2;
	double T120 = 380;              //Temperature at 120 km[K]
	double R = 6.35677e3;           //Equitorial radius of Earth[km]
	double T120p = 14.348;          //Relative temperature vertical gradient[K / km]

	double MU = 398600.4418;        //Gravitational parameter of Earth[km ^ 3 s^-2]
	double BOLTZ = 1.38064852;      //Boltzmann constant[10 ^ 29 km ^ 2 kg s^-2 K - 1]

	double MASS_H = 1.6737236e2;    //Mass of Hydrogen[10 ^ 29 kg]
	double MASS_HE = 6.6464764e2;   //Mass of Helium[10 ^ 29 kg]
	double MASS_O = 2.6567626e3;    //Mass of Oxygen[10 ^ 29 kg]
	double MASS_N2 = 2 * 2.3258671e3; //Mass of di - Nitrogen[10 ^ 29 kg]

	g = MU / ((R + 120)*(R+120));
	xi = (z - 120)*(R + 120) / (z + R);
	t_inf= calculate_t_inf_and_initial_concentrations(P, f, 
		f_mean, k, t, d, &cH_0, &cHe_0, &cO_0, &cN2_0);

	sigma = T120p / (t_inf - T120);

	gammaH = MASS_H*g / (sigma*BOLTZ*t_inf);
	gammaHe = MASS_HE*g / (sigma*BOLTZ*t_inf);
	gammaO = MASS_O*g / (sigma*BOLTZ*t_inf);
	gammaN2 = MASS_N2*g / (sigma*BOLTZ*t_inf);

	alphaH = -0.38;
	alphaHe = -0.38;

	*T = t_inf - (t_inf - T120)*exp(-sigma*xi);

	fH = pow((T120 / *T) ,(1 + alphaH + gammaH)*exp(-sigma*gammaH*xi));
	fHe = pow((T120 / *T) , (1 + alphaHe + gammaHe)*exp(-sigma*gammaHe*xi));
	fO = pow((T120 / *T) , (1 + gammaO)*exp(-sigma*gammaO*xi));
	fN2 = pow((T120 / *T) , (1 + gammaN2)*exp(-sigma*gammaN2*xi));

	*c_H = cH_0*fH;
	*c_He = cHe_0*fHe;
	*c_O = cO_0*fO;
	*c_N2 = cN2_0*fN2;
	return t_inf;
}