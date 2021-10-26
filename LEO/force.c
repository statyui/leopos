
#include "leoorb.h"

#define R_EARTH		6378.137e3      // Radius Earth [m]; WGS-84
#define f_EARTH		(1.0 / 298.257223563) // Flattening; WGS-84   
#define R_SUN		696000.0e3       // Radius Sun [m]; Seidelmann 1992
#define R_MOON		1738.0e3       // Radius Moon [m]
#define  NHPCOEF   50

// Earth rotation (derivative of GMST at J2000; differs from inertial period by precession)
#define OMEGA_EARTH   7.2921158553e-5    // [rad/s]; Aoki 1982, NIMA 1997

// Gravitational coefficient

#define GM_SUN		1.32712438e+20    // [m^3/s^2]; IAU 1976
#define GM_MOON		((GM_EARTH) / 81.300587)// [m^3/s^2]; DE200
// Solar radiation pressure at 1 AU 
#define SRP0       4.560E-6  // [N/m^2] (~1367 W/m^2); IERS 96

	// Fractional part of a number (y=x-[x])
	//double Frac(double x) { return x - floor(x); };

// Earth gravity model (JGM3)

//GravModel Grav(N_JGM3,            // Degree and order
//	GM_Earth,          // Gravitational coefficient [m^3/s^2]
//	R_JGM3,            // Reference radius [m]
//	&CS_JGM3[0][0]);  // Unnormalized harmonic coefficients


//------------------------------------------------------------------------------
// 
// Illumination
//
// Purpose:
//
//   Computes the fractional illumination of a spacecraft in the 
//   vicinity of the Earth assuming a cylindrical shadow model 
// 
// Input/output:
// 
//   r               Spacecraft position vector [m]
//   r_Sun           Sun position vector [m]
//   <return>        Illumination factor:
//                     nu=0   Spacecraft in Earth shadow 
//                     nu=1   Spacecraft fully illuminated by the Sun
//
//------------------------------------------------------------------------------

int illumination(const double *r, const double *rsun)
{
	double nsun = 0.0;
	double esun[3] = { 0.0 };
	double s = 0.0;
	double es[3] = { 0.0 };
	int i;
	nsun = norm(rsun,3);
	for (i = 0; i < 3; i++) {
		esun[i] = rsun[i] / nsun;
	}
	s = dot(r, esun, 3);
	for (i = 0; i < 3; i++) {
		es[i] = r[i] - s*esun[i];
	}
	if (s > 0 || norm(es,3) > R_EARTH) {
		return 1;
	}
	else {
		return 0;
	}
	//Vector e_Sun = r_Sun / Norm(r_Sun);   // Sun direction unit vector
	//double s = Dot(r, e_Sun);      // Projection of s/c position 

	//return ((s>0 || Norm(r - s*e_Sun)>R_Earth) ? 1.0 : 0.0);
}




//------------------------------------------------------------------------------
//
// pointattract
//
// Purpose:
//
//   Computes the perturbational acceleration due to a point mass
//
// Input/Output:
//
//   r           Satellite position vector 
//   s           Point mass position vector
//   GM          Gravitational coefficient of point mass
//   <return>    Acceleration (a=d^2r/dt^2)
//
//------------------------------------------------------------------------------

int pointattract(const double *r, const double *s, double GM,double *a)
{
	double d[3] = { 0.0 };
	int i;
	//Vector d(3);
	////  Relative position vector of satellite w.r.t. point mass 
	//d = r - s;
	//// Acceleration 
	//return  (-GM) * (d / pow(Norm(d), 3) + s / pow(Norm(s), 3));
	for (i = 0; i < 3; i++) { d[i] = r[i] - s[i]; }
	for (i = 0; i < 3; i++) {
		a[i] = -GM*(d[i] / pow(norm(d, 3), 3) + s[i] / pow(norm(s, 3), 3));
	}
	return 1;

}


//------------------------------------------------------------------------------
//
// solarrad
//
// Purpose:
//
//   Computes the acceleration due to solar radiation pressure assuming 
//   the spacecraft surface normal to the Sun direction
//
// Input/Output:
//
//   r           Spacecraft position vector 
//   r_Sun       Sun position vector 
//   Area        Cross-section 
//   mass        Spacecraft mass
//   CR          Solar radiation pressure coefficient
//   P0          Solar radiation pressure at 1 AU 
//   AU          Length of one Astronomical Unit 
//   <return>    Acceleration (a=d^2r/dt^2)
//
// Notes:
//
//   r, r_sun, Area, mass, P0 and AU must be given in consistent units,
//   e.g. m, m^2, kg and N/m^2. 
//
//------------------------------------------------------------------------------
double solarrad(double *r, double *rsun, double area, double mass, double CR,double *a) {
	double d[3] = { 0.0 };
	int i;
	for (i = 0; i < 3; i++) {
		d[i] = r[i] - rsun[i];
	}
	for (i = 0; i < 3; i++) {
		a[i] = CR*(area / mass)*SRP0*(AU*AU)*d[i] / pow(norm(d, 3), 3);
	}
	return 1;
}
//Vector AccelSolrad(const Vector& r, const Vector& r_Sun,
//	double Area, double mass, double CR,
//	double P0, double AU)
//{
//
//	Vector d(3);
//
//	// Relative position vector of spacecraft w.r.t. Sun
//
//	d = r - r_Sun;
//
//	// Acceleration 
//
//	return  CR*(Area / mass)*P0*(AU*AU) * d / pow(Norm(d), 3);
//
//}

//------------------------------------------------------------------------------
//
// Density_HP
//
// Purpose:
//
//   Computes the atmospheric density for the modified Harris-Priester model.
//
// Input/Output:
//
//   Mjd_TT      Terrestrial Time (Modified Julian Date)
//   r_tod       Satellite position vector in the inertial system [m]
//   <return>    Density [kg/m^3]
//
//------------------------------------------------------------------------------

double Density_HP(double Mjd_TT, double *rsun, const double *rtod)
{
	// Constants

	const double upper_limit = 1000.0;           // Upper height limit [km]
	const double lower_limit = 100.0;           // Lower height limit [km]
	const double ra_lag = 0.523599;           // Right ascension lag [rad]
	const int    n_prm = 3;           // Harris-Priester parameter 
									  // 2(6) low(high) inclination

									  // Harris-Priester atmospheric density model parameters 
									  // Height [km], minimum density, maximum density [gm/km^3]


	const double h[NHPCOEF] = {
		100.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0,
		210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 270.0, 280.0, 290.0, 300.0,
		320.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0,
		520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 640.0, 660.0, 680.0, 700.0,
		720.0, 740.0, 760.0, 780.0, 800.0, 840.0, 880.0, 920.0, 960.0,1000.0 };
	const double c_min[NHPCOEF] = {
		4.974e+05, 2.490e+04, 8.377e+03, 3.899e+03, 2.122e+03, 1.263e+03,
		8.008e+02, 5.283e+02, 3.617e+02, 2.557e+02, 1.839e+02, 1.341e+02,
		9.949e+01, 7.488e+01, 5.709e+01, 4.403e+01, 3.430e+01, 2.697e+01,
		2.139e+01, 1.708e+01, 1.099e+01, 7.214e+00, 4.824e+00, 3.274e+00,
		2.249e+00, 1.558e+00, 1.091e+00, 7.701e-01, 5.474e-01, 3.916e-01,
		2.819e-01, 2.042e-01, 1.488e-01, 1.092e-01, 8.070e-02, 6.012e-02,
		4.519e-02, 3.430e-02, 2.632e-02, 2.043e-02, 1.607e-02, 1.281e-02,
		1.036e-02, 8.496e-03, 7.069e-03, 4.680e-03, 3.200e-03, 2.210e-03,
		1.560e-03, 1.150e-03 };
	const double c_max[NHPCOEF] = {
		4.974e+05, 2.490e+04, 8.710e+03, 4.059e+03, 2.215e+03, 1.344e+03,
		8.758e+02, 6.010e+02, 4.297e+02, 3.162e+02, 2.396e+02, 1.853e+02,
		1.455e+02, 1.157e+02, 9.308e+01, 7.555e+01, 6.182e+01, 5.095e+01,
		4.226e+01, 3.526e+01, 2.511e+01, 1.819e+01, 1.337e+01, 9.955e+00,
		7.492e+00, 5.684e+00, 4.355e+00, 3.362e+00, 2.612e+00, 2.042e+00,
		1.605e+00, 1.267e+00, 1.005e+00, 7.997e-01, 6.390e-01, 5.123e-01,
		4.121e-01, 3.325e-01, 2.691e-01, 2.185e-01, 1.779e-01, 1.452e-01,
		1.190e-01, 9.776e-02, 8.059e-02, 5.741e-02, 4.210e-02, 3.130e-02,
		2.360e-02, 1.810e-02 };

	//const Vector h(&Data_h[0], NHPCOEF);
	//const Vector c_min(&Data_c_min[0], NHPCOEF);
	//const Vector c_max(&Data_c_max[0], NHPCOEF);

	// Variables
	int    i, ih;                              // Height section variables        
	double height;                             // Earth flattening
	double dec_Sun, ra_Sun, c_dec;             // Sun declination, right asc.
	double c_psi2;                             // Harris-Priester modification
	double density, h_min, h_max, d_min, d_max;// Height, density parameters
											   //Vector r_Sun(3);                           // Sun position
											   //Vector u(3);                               // Apex of diurnal bulge
	double u[3] = { 0.0 };
	double pos[3] = { 0.0 };

	// Satellite height
	ecef2pos(rtod, pos);
	height = pos[2] / 1000.0;   //  [km]


	// Exit with zero density outside height model limits

	if (height >= upper_limit || height <= lower_limit)
	{
		return 0.0;
	}

	// Sun right ascension, declination

	//r_Sun = Sun(Mjd_TT);
	ra_Sun = atan2(rsun[1], rsun[0]);
	dec_Sun = atan2(rsun[2], sqrt(pow(rsun[0], 2) + pow(rsun[1], 2)));


	// Unit vector u towards the apex of the diurnal bulge
	// in inertial geocentric coordinates

	c_dec = cos(dec_Sun);
	u[0] = c_dec * cos(ra_Sun + ra_lag);
	u[1] = c_dec * sin(ra_Sun + ra_lag);
	u[2] = sin(dec_Sun);


	// Cosine of half angle between satellite position vector and
	// apex of diurnal bulge

	c_psi2 = 0.5 + 0.5 * dot(rtod, u, 3) / norm(rtod, 3);

	// Height index search and exponential density interpolation

	ih = 0;                           // section index reset
	for (i = 0; i<NHPCOEF - 1; i++)       // loop over N_Coef height regimes
	{
		if (height >= h[i] && height < h[i + 1])
		{
			ih = i;                       // ih identifies height section
			break;
		}
	}

	h_min = (h[ih] - h[ih + 1]) / log(c_min[ih + 1] / c_min[ih]);
	h_max = (h[ih] - h[ih + 1]) / log(c_max[ih + 1] / c_max[ih]);

	d_min = c_min[ih] * exp((h[ih] - height) / h_min);
	d_max = c_max[ih] * exp((h[ih] - height) / h_max);

	// Density computation

	density = d_min + (d_max - d_min)*pow(c_psi2, n_prm);


	return density * 1.0e-12;       // [kg/m^3]

}

//------------------------------------------------------------------------------
//
// atmosdrag
//
// Purpose:
//
//   Computes the acceleration due to the atmospheric drag.
//
// Input/Output:
//
//   Mjd_TT      Terrestrial Time (Modified Julian Date)
//   r           Satellite position vector in the inertial system [m]
//   v           Satellite velocity vector in the inertial system [m/s]
//   T           Transformation matrix to true-of-date inertial system
//   Area        Cross-section [m^2]
//   mass        Spacecraft mass [kg]
//   CD          Drag coefficient
//   <return>    Acceleration (a=d^2r/dt^2) [m/s^2]
//
//------------------------------------------------------------------------------
int atmosdrag(double Mjd_TT, double *r, double *v, double *rsun,double *T, double area, double mass, double CD, double *a)
//Vector AccelDrag(double Mjd_TT, const Vector& r, const Vector& v,
//	const Matrix& T,double Area, double mass, double CD)
{

	// Constants

	// Earth angular velocity vector [rad/s]
	const double omega[3] = { 0.0, 0.0, 7.29212e-5 };
	// Variables

	double v_abs, dens;
	double rtod[3] = { 0.0 };
	double vtod[3] = { 0.0 };
	double vrel[3] = { 0.0 };
	double atod[3] = { 0.0 };
	double ttrp[9] = { 0.0 };
	double crosv[3] = { 0.0 };
	int i;
	//Vector r_tod(3), v_tod(3);
	//Vector v_rel(3), a_tod(3);
	//Matrix T_trp(3, 3);


	//// Transformation matrix to ICRF/EME2000 system
	//T_trp = Transp(T);
	//// Position and velocity in true-of-date system
	//r_tod = T * r;
	//v_tod = T * v;
	matmul("NN", 3, 1, 3, 1.0, T, r, 0.0, rtod);
	matmul("NN", 3, 1, 3, 1.0, T, v, 0.0, vtod);
	cross3(omega, rtod,crosv);
	for (i = 0; i < 3; i++) { vrel[i] = vtod[i] - crosv[i]; }
	v_abs = norm(vrel, 3);
	// Velocity relative to the Earth's atmosphere
	//v_rel = v_tod - Cross(omega, r_tod);
	//v_abs = Norm(v_rel);

	// Atmospheric density due to modified Harris-Priester model
	dens = Density_HP(Mjd_TT,rsun, rtod);

	// Acceleration 
	//a_tod = -0.5*CD*(Area / mass)*dens*v_abs*v_rel;
	for (i = 0; i < 3; i++) { atod[i] = -0.5*CD*(area / mass)*dens*v_abs*vrel[i]; }
	matmul("TN", 3, 1, 3, 1.0, T, atod, 0.0, a);
	//return T_trp * a_tod;
	return 1;
}

//------------------------------------------------------------------------------
//
// AccelMain
//
// Purpose:
//
//   Computes the acceleration of an Earth orbiting satellite due to 
//    - the Earth's harmonic gravity field, 
//    - the gravitational perturbations of the Sun and Moon
//    - the solar radiation pressure and
//    - the atmospheric drag
// Input/Output:
//
//   tutc        time in UT1
//   r           Satellite position vector in the ICRF/EME2000 system
//   v           Satellite velocity vector in the ICRF/EME2000 system
//   area        Cross-section 
//   mass        Spacecraft mass
//   CR          Radiation pressure coefficient
//   CD          Drag coefficient
//   a>          Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system
//
//------------------------------------------------------------------------------
extern int forces(gtime_t tutc,const double *erpv, double *r,double *v, double area, double mass, double CR, double CD,double *a)
{

	double rsun[3] = { 0.0 };
	double rmoon[3] = { 0.0 };
	double T[9] = { 0.0 };
	double E[9] = { 0.0 };
	double U[9] = { 0.0 };
	double agf[3] = { 0.0 };
	double ast[3] = { 0.0 };
	double amt[3] = { 0.0 };
	double asr[3] = { 0.0 };
	int illu = 0;
	double aad[3] = { 0.0 };
	double gmst=0.0;
	int i;
	double mjd=0.0;
	mjd=time2mjd(tutc);
	eci2ecef(tutc, erpv, U, &gmst);
	sunmoonpos_eci(tutc, rsun, rmoon);
	//// Acceleration due to harmonic gravity field
	//Mjd_UT1 = Mjd_TT;
	//T = NutMatrix(Mjd_TT) * PrecMatrix(MJD_J2000, Mjd_TT);
	//E = GHAMatrix(Mjd_UT1) * T;

	//a = AccelHarmonic(r, E, Grav.GM, Grav.R_ref, Grav.CS, Grav.n_max, Grav.m_max);
	// Luni-solar perturbations 
	//r_Sun = Sun(Mjd_TT);
	//r_Moon = Moon(Mjd_TT);
	gravityJGM3(r, U, 20, 20, agf);
	pointattract(r, rsun, GM_SUN, ast);
	pointattract(r, rmoon, GM_MOON, amt);
	solarrad(r, rsun, area, mass, CR, asr);
	illu = illumination(r,rsun);
	atmosdrag(mjd, r, v,rsun, U, area, mass, CD, aad);
	for (i = 0; i < 3; i++) {
		a[i] = agf[i] + ast[i] + amt[i] + illu*asr[i] + aad[i];
	}
	//a += AccelPointMass(r, r_Sun, GM_Sun);
	//a += AccelPointMass(r, r_Moon, GM_Moon);
	//// Solar radiation pressure
	//a += Illumination(r, r_Sun)* AccelSolrad(r, r_Sun, Area, mass, CR, P_Sol, AU);
	//// Atmospheric drag
	//a += AccelDrag(Mjd_TT, r, v, T, Area, mass, CD);
	//// Acceleration
	return 1;

}

//------------------------------------------------------------------------------
//
// Accel
//
// Purpose:
//
//   Computes the acceleration of an Earth orbiting satellite due to 
//   the Earth's harmonic gravity field
//
// Input/Output:
//
//   Mjd_GPS     GPS Time (Modified Julian Date)
//   r           Satellite position vector in the pseudo-true-of-date system
//   n           Degree and order of gravity field
//   <return>    Acceleration (a=d^2r/dt^2) in the pseudo-true-of-date system
//
//------------------------------------------------------------------------------

int forcessimple(gtime_t tutc, const double *erpv, const double *r, int n,double *a)
{
	// Variables 
	double  U[9] = { 0.0 }, dU[9] = {0.0};
	double gmst=0.0;
	//Vector  a(3);
	//Matrix  U(3, 3);
	// Acceleration due to harmonic gravity field
	//eci2ecef(tutc, erpv, U, &gmst);
	ecef2qeci(tutc, erpv[2], U, dU);
	//Mjd_UT1 = Mjd_GPS + (IERS::UT1_UTC(Mjd_GPS) - IERS::GPS_UTC(Mjd_GPS)) / 86400.0;

	//U = R_z(GMST(Mjd_UT1));
	//a = AccelHarmonic(r, U, GM_Earth, Grav.R_ref, Grav.CS, n, n);
	gravityJGM3(r, U, n, n, a);
	// Acceleration
	return 1;

}

//------------------------------------------------------------------------------
//
// EccAnom
//
// Purpose:
//
//   Computes the eccentric anomaly for elliptic orbits
//
// Input/Output:
//
//   M         Mean anomaly in [rad]
//   e         Eccentricity of the orbit [0,1[
//   <return>  Eccentric anomaly in [rad]
//
//------------------------------------------------------------------------------

double EccAnom(double M, double e)
{

	// Constants
	const int maxit = 15;
	const double eps = 1e-10;

	// Variables
	int    i = 0;
	double E, f;

	// Starting value
	M = fmod(M, 2.0*PI);
	if (e<0.8) E = M; else E = PI;

	// Iteration
	do {
		f = E - e*sin(E) - M;
		E = E - f / (1.0 - e*cos(E));
		++i;
		if (i == maxit) {
			break;
		}
	} while (fabs(f) > eps);
	return E;

}

//------------------------------------------------------------------------------
//
// pv2kepler
//
// Purpose:
//
//   Computes the osculating Keplerian elements from the satellite state vector
//   for elliptic orbits
//
// Input/Output:
//
//   GM        Gravitational coefficient
//             (gravitational constant * mass of central body)
//   y         State vector (x,y,z,vx,vy,vz)
//   <return>  Keplerian elements (a,e,i,Omega,omega,M) with
//               a      Semimajor axis 
//               e      Eccentricity 
//               i      Inclination [rad]
//               Omega  Longitude of the ascending node [rad]
//               omega  Argument of pericenter  [rad]
//               M      Mean anomaly  [rad]
//
// Notes:
//
//   The state vector and GM must be given in consistent units, 
//   e.g. [m], [m/s] and [m^3/s^2]. The resulting unit of the semimajor
//   axis is implied by the unity of y, e.g. [m].
//
//   The function cannot be used with state vectors describing a circular
//   or non-inclined orbit.
//
//------------------------------------------------------------------------------

int pv2kepler(double GM, const double *y, double *kepler)
{

	// Variables
	//Vector  r(3), v(3), h(3);
	double r[3] = { 0.0 };
	double v[3] = { 0.0 };
	double h[3] = { 0.0 };
	double  H, u, R;
	double  eCosE, eSinE, e2, E, nu;
	double  a, e, i, Omega, omega, M;
	int k;
	trace(3, "pv2kepler:\n");
	for (k = 0; k < 3; k++) {
		r[k] = y[k];
		v[k] = y[k + 3];
	}
	cross3(r, v, h);
	R = norm(r, 3);                                // Distance           
	a = 1.0 / (2.0 / R - dot(v, v, 3) / GM);       // Semi-major axis   
	H = norm(h, 3); 
	if (H < 1e-12 || a < 0) {
		return 0; /*no rotating orbit*/
	}
	Omega = atan2(h[0], -h[1]);                     // Long. ascend. node 
	Omega = fmod(Omega, PI*2);
	if (Omega < 0) Omega = Omega + 2 * PI;
	i = atan2(sqrt(h[0]*h[0] + h[1]*h[1]), h[2]); // Inclination        
	u = atan2(r[2]*H, -r[0]*h[1] + r[1]*h[0]);    // Arg. of latitude   
	eCosE = 1.0 - R / a;                          // e*cos(E)           
	eSinE = dot(r, v,3) / sqrt(GM*a);             // e*sin(E)           
	e2 = eCosE*eCosE + eSinE*eSinE;
	e = sqrt(e2);                                 // Eccentricity 
	E = atan2(eSinE, eCosE);                      // Eccentric anomaly  
	M = fmod(E - eSinE, PI*2);     
	if (M < 0) M = M + 2 * PI;// Mean anomaly
	nu = atan2(sqrt(1.0 - e2)*eSinE, eCosE - e2);          // True anomaly
	omega = fmod(u - nu, PI*2);    
	if (omega < 0) omega = omega + 2 * PI;
	// Arg. of perihelion 
	// Keplerian elements vector
	kepler[0] = a; kepler[1] = e; kepler[2] = i;
	kepler[3] = Omega; kepler[4] = omega; kepler[5] = M;
	return 1;
}
//
// F : local function for use by FindEta()
// F = 1 - eta +(m/eta**2)*W(m/eta**2-l)
//

double F(double eta, double m, double l)
{
	// Constants
	const double eps = 100.0 * 1e-16;
	// Variables
	double  w, W, a, n, g;

	w = m / (eta*eta) - l;
	if (fabs(w)<0.1) { // Series expansion
		W = a = 4.0 / 3.0; n = 0.0;
		do {
			n += 1.0;  a *= w*(n + 2.0) / (n + 1.5);  W += a;
		} while (fabs(a) >= eps);
	}
	else {
		if (w > 0.0) {
			g = 2.0*asin(sqrt(w));
			W = (2.0*g - sin(2.0*g)) / pow(sin(g), 3);
		}
		else {
			g = 2.0*log(sqrt(-w) + sqrt(1.0 - w));  // =2.0*arsinh(sqrt(-w))
			W = (sinh(2.0*g) - 2.0*g) / pow(sinh(g), 3);
		}
	}
	return (1.0 - eta + (w + l)*W);
}   // End of function F


//------------------------------------------------------------------------------
//
// FindEta 
//
//   Computes the sector-triangle ratio from two position vectors and 
//   the intermediate time 
//
// Input/Output:
//
//   r_a        Position at time t_a
//   r_a        Position at time t_b
//   tau        Normalized time (sqrt(GM)*(t_a-t_b))
//   <return>   Sector-triangle ratio
//
//------------------------------------------------------------------------------

double FindEta(const double *r_a, const double *r_b, double tau)
{
	// Constants
	const int maxit = 30;
	const double delta = 100.0*1e-10;

	// Variables
	int    i;
	double kappa, m, l, s_a, s_b, eta_min, eta1, eta2, F1, F2, d_eta;

	// Auxiliary quantities

	s_a = norm(r_a,3);
	s_b = norm(r_b,3);

	kappa = sqrt(2.0*(s_a*s_b + dot(r_a, r_b,3)));

	m = tau*tau / pow(kappa, 3);
	l = (s_a + s_b) / (2.0*kappa) - 0.5;

	eta_min = sqrt(m / (l + 1.0));

	// Start with Hansen's approximation
	eta2 = (12.0 + 10.0*sqrt(1.0 + (44.0 / 9.0)*m / (l + 5.0 / 6.0))) / 22.0;
	eta1 = eta2 + 0.1;

	// Secant method
	F1 = F(eta1, m, l);
	F2 = F(eta2, m, l);

	i = 0;

	while (fabs(F2 - F1) > delta)
	{
		d_eta = -F2*(eta2 - eta1) / (F2 - F1);
		eta1 = eta2; F1 = F2;
		while (eta2 + d_eta <= eta_min)  d_eta *= 0.5;
		eta2 += d_eta;
		F2 = F(eta2, m, l); ++i;

		if (i == maxit) {
			assert(0);
			break;
		}
	}

	return eta2;
}

//------------------------------------------------------------------------------
//
// GetVo
//
//
//------------------------------------------------------------------------------

extern int GetVo(double GM, double dt, const double *r_a, const double *r_b, double *pv) {
	
	double  tau, eta, p;
	double  n, nu, E, u, R;
	double  s_a, s_b, s_0, fac, sinhH;
	double  cos_dnu, sin_dnu, ecos_nu, esin_nu;
	double  a, e, i, Omega, omega, M;
	double  e_a[3] = { 0.0 }, r_0[3] = { 0.0 }, r_c[3] = { 0.0 }, r_r[3] = { 0.0 },e_0[3] = { 0.0 }, W[3] = { 0.0 }, v[3] = { 0.0 };
	int k;


	for (k = 0; k < 3; k++) { r_0[k] = r_b[k] - r_a[k]; }
	
	
	
	//s_a = norm(r_a, 3);
	//for (k = 0; k < 3; k++) { e_a[k] = r_a[k] / s_a; }
	//s_b = norm(r_b, 3);
	//fac = dot(r_b, e_a, 3);
	//for (k = 0; k < 3; k++) { r_0[k] = r_b[k] - fac * e_a[k]; }
	//s_0 = norm(r_0, 3);
	//for (k = 0; k < 3; k++) { e_0[k] = r_0[k] / s_0; }

	//for (k = 0; k < 3; k++) { r_r[k] = r_c[k] * e_0[k]; }

	for (k = 0; k < 3; k++) {
		pv[k] = r_a[k];
		pv[k + 3] = r_0[k] / dt;
	}

	return 1;

}





//------------------------------------------------------------------------------
//
// Elements 
//
// Purpose:
//
//   Computes orbital elements from two given position vectors and 
//   associated times 
//
// Input/Output:
//
//   GM        Gravitational coefficient
//             (gravitational constant * mass of central body)
//   Mjd_a     Time t_a (Modified Julian Date)
//   Mjd_b     Time t_b (Modified Julian Date)
//   r_a       Position vector at time t_a
//   r_b       Position vector at time t_b
//   <return>  Keplerian elements (a,e,i,Omega,omega,M)
//               a      Semimajor axis 
//               e      Eccentricity 
//               i      Inclination [rad]
//               Omega  Longitude of the ascending node [rad]
//               omega  Argument of pericenter  [rad]
//               M      Mean anomaly  [rad]
//             at time t_a 
//
// Notes:
//
//   The function cannot be used with state vectors describing a circular
//   or non-inclined orbit.
//
//------------------------------------------------------------------------------

extern int p22kepler(double GM, double dt,const double *r_a, const double *r_b, double *kep)
{

	// Variables

	double  tau, eta, p;
	double  n, nu, E, u, R;
	double  s_a, s_b, s_0, fac, sinhH;
	double  cos_dnu, sin_dnu, ecos_nu, esin_nu;
	double  a, e, i, Omega, omega, M;
	double  e_a[3] = { 0.0 }, r_0[3] = { 0.0 }, e_0[3] = { 0.0 }, W[3] = { 0.0 }, v[3] = { 0.0 };
	int k;

	//Vector  e_a(3), r_0(3), e_0(3), W(3);
	trace(3, "p22kepler dt=%13.4f r_a=%13.4f %13.4f %13.4f r_b=%13.4f %13.4f %13.4f\n", dt, r_a[0], r_a[1], r_a[2], r_b[0], r_b[1], r_b[2]);
	// Calculate vector r_0 (fraction of r_b perpendicular to r_a) 
	// and the magnitudes of r_a,r_b and r_0

	s_a = norm(r_a,3);
	for (k = 0; k < 3; k++) { e_a[k] = r_a[k] / s_a; }
	s_b = norm(r_b,3);
	fac = dot(r_b, e_a,3);
	for (k = 0; k < 3; k++) {r_0[k] = r_b[k] - fac*e_a[k];}
	s_0 = norm(r_0,3);
	for (k = 0; k < 3; k++) { e_0[k] = r_0[k] / s_0; }

	// Inclination and ascending node 
	cross3(e_a, e_0,W);
	Omega = atan2(W[0], -W[1]);                     // Long. ascend. node 
	Omega = fmod(Omega, PI*2);
	if (Omega < 0) Omega = Omega + 2 * PI;
	i = atan2(sqrt(W[0]*W[0] + W[1]*W[1]), W[2]); // Inclination     
	if (i == 0.0)
		u = atan2(r_a[1], r_a[0]);
	else
		u = atan2(+e_a[2], -e_a[0]*W[1] + e_a[1]*W[0]);

	// Semilatus rectum
	//tau = sqrt(GM) * 86400.0*fabs(Mjd_b - Mjd_a);
	tau = sqrt(GM) * 86400.0*fabs(dt);
	eta = FindEta(r_a, r_b, tau);
	p = pow(s_a*s_0*eta / tau, 2);

	// Eccentricity, true anomaly and argument of perihelion
	cos_dnu = fac / s_b;
	sin_dnu = s_0 / s_b;

	ecos_nu = p / s_a - 1.0;
	esin_nu = (ecos_nu * cos_dnu - (p / s_b - 1.0)) / sin_dnu;

	e = sqrt(ecos_nu*ecos_nu + esin_nu * esin_nu);
	nu = atan2(esin_nu, ecos_nu);

	omega = fmod(u - nu, PI * 2);
	if (omega < 0) omega = omega + 2 * PI;
	// Perihelion distance, semimajor axis and mean motion
	a = p / (1.0 - e * e);
	n = sqrt(GM / fabs(a*a*a));

	// Mean anomaly and time of perihelion passage
	if (e < 1.0) {
		E = atan2(sqrt((1.0 - e)*(1.0 + e)) * esin_nu, ecos_nu + e * e);
		M = fmod(E - e * sin(E), 2 * PI);
		if (M < 0) M = M + 2 * PI;
	}
	else
	{
		sinhH = sqrt((e - 1.0)*(e + 1.0)) * esin_nu / (e + e * ecos_nu);
		M = e * sinhH - log(sinhH + sqrt(1.0 + sinhH * sinhH));
	}


	// Keplerian elements vector
	kep[0] = a;
	kep[1] = e;
	kep[2] = i;
	kep[3] = Omega;
	kep[4] = omega;
	kep[5] = M;
	return 1;

}
//------------------------------------------------------------------------------
//
// State
//
// Purpose:
//
//   Computes the satellite state vector from osculating Keplerian elements 
//   for elliptic orbits
//
// Input/Output:
//
//   GM        Gravitational coefficient
//             (gravitational constant * mass of central body)
//   Kep       Keplerian elements (a,e,i,Omega,omega,M) with
//               a      Semimajor axis 
//               e      Eccentricity 
//               i      Inclination [rad]
//               Omega  Longitude of the ascending node [rad]
//               omega  Argument of pericenter  [rad]
//               M      Mean anomaly at epoch [rad]
//   dt        Time since epoch
//   <return>  State vector (x,y,z,vx,vy,vz)
//
// Notes:
//
//   The semimajor axis a=Kep(0), dt and GM must be given in consistent units, 
//   e.g. [m], [s] and [m^3/s^2]. The resulting units of length and velocity  
//   are implied by the units of GM, e.g. [m] and [m/s].
//
//------------------------------------------------------------------------------
extern int kepler2pv(double GM, const double *kep, double dt,double *pv)
{

	// Variables

	double  a, e, i, Omega, omega, M, M0, n;
	double  E, cosE, sinE, fac, R, V;
	//Vector  r(3), v(3);
	//Matrix  PQW(3, 3);
	double r[3] = { 0.0 }, r1[3] = {0.0};
	double v[3] = { 0.0 }, v1[3] = {0.0};
	double PQW[9] = { 0.0 };
	double Rp[9] = { 0.0 }, Rq[9] = { 0.0 }, Rw[9] = {0.0}, Rpq[9] = { 0.0 };
	int k;
	trace(3, "kepler2pv:dt=%.3f\n", dt);
	// Keplerian elements at epoch
	a = kep[0];  Omega = kep[3];
	e = kep[1];  omega = kep[4];
	i = kep[2];  M0 = kep[5];

	// Mean anomaly 

	if (dt == 0.0) {
		M = M0;
	}
	else {
		n = sqrt(GM / (a*a*a));
		M = M0 + n*dt;
	};

	// Eccentric anomaly 
	E = EccAnom(M, e);

	cosE = cos(E);
	sinE = sin(E);

	// Perifocal coordinates
	fac = sqrt((1.0 - e)*(1.0 + e));
	R = a*(1.0 - e*cosE);  // Distance
	V = sqrt(GM*a) / R;    // Velocity
	
	//r = Vector(a*(cosE - e), a*fac*sinE, 0.0);
	//v = Vector(-V*sinE, +V*fac*cosE, 0.0);
	r[0] = a*(cosE - e);
	r[1] = a*fac*sinE;
	v[0] = -V*sinE;
	v[1] = +V*fac*cosE;

	// Transformation to reference system (Gaussian vectors)
	//PQW = R_z(-Omega) * R_x(-i) * R_z(-omega);
	//r = PQW*r;
	//v = PQW*v;
	Rz(-Omega, Rp);
	Rx(-i, Rq);
	Rz(-omega, Rw);

	matmul("NN", 3, 3, 3, 1.0, Rp, Rq, 0.0, Rpq);
	matmul("NN", 3, 3, 3, 1.0, Rpq, Rw, 0.0, PQW);
	matmul("NN", 3, 1, 3, 1.0, PQW, r, 0.0, r1);
	matmul("NN", 3, 1, 3, 1.0, PQW, v, 0.0, v1);
	// State vector 
	for (k = 0; k < 3; k++) {
		pv[k] = r1[k];
		pv[k + 3] = v1[k];
	}
	return 1;
	//return Stack(r, v);

}


//------------------------------------------------------------------------------
//
// StatePartials
//
// Purpose:
//
//   Computes the partial derivatives of the satellite state vector with respect
//   to the orbital elements for elliptic, Keplerian orbits
//
// Input/Output:
//
//   GM        Gravitational coefficient
//             (gravitational constant * mass of central body)
//   Kep       Keplerian elements (a,e,i,Omega,omega,M) at epoch with
//               a      Semimajor axis 
//               e      Eccentricity 
//               i      Inclination [rad]
//               Omega  Longitude of the ascending node [rad]
//               omega  Argument of pericenter  [rad]
//               M      Mean anomaly at epoch [rad]
//   dt        Time since epoch
//   <return>  Partials derivatives of the state vector (x,y,z,vx,vy,vz) at time
//             dt with respect to the epoch orbital elements
//
// Notes:
//
//   The semimajor axis a=Kep(0), dt and GM must be given in consistent units, 
//   e.g. [m], [s] and [m^3/s^2]. The resulting units of length and velocity  
//   are implied by the units of GM, e.g. [m] and [m/s].
//
//   The function cannot be used with circular or non-inclined orbit.
//
//------------------------------------------------------------------------------

int keplerpartials(double GM, const double *Kep, double dt,double *dYdA)
{

	// Variables

	int     k;
	double  a, e, i, Omega, omega, M, M0, n, dMda;
	double  E, cosE, sinE, fac, r, v, x, y, vx, vy;
	trace(3, "keplerpartials dt=%.3f\n", dt);
	double PQW[9] = { 0.0 };
	double P[3] = { 0.0 }, Q[3] = { 0.0 }, W[3] = { 0.0 };
	double e_z[3] = { 0.0 };
	double N[3] = { 0.0 };
	double dPdi[3] = { 0.0 }, dPdO[3] = { 0.0 }, dPdo[3] = { 0.0 };
	double dQdi[3] = { 0.0 }, dQdO[3] = { 0.0 }, dQdo[3] = { 0.0 };
	double dYda[6] = { 0.0 }, dYde[6] = { 0.0 }, dYdi[6] = { 0.0 };
	double dYdO[6] = { 0.0 }, dYdo[6] = { 0.0 }, dYdM[6] = { 0.0 };
	double Rp[9] = { 0.0 }, Rq[9] = { 0.0 }, Rw[9] = { 0.0 }, Rpq[9] = { 0.0 };
	double nN;
	// Keplerian elements at epoch

	a = Kep[0];  Omega = Kep[3];
	e = Kep[1];  omega = Kep[4];
	i = Kep[2];  M0 = Kep[5];

	// Mean and eccentric anomaly

	n = sqrt(GM / (a*a*a));
	M = M0 + n*dt;
	E = EccAnom(M, e);

	// Perifocal coordinates

	cosE = cos(E);
	sinE = sin(E);
	fac = sqrt((1.0 - e)*(1.0 + e));

	r = a*(1.0 - e*cosE);  // Distance
	v = sqrt(GM*a) / r;    // Velocity

	x = +a*(cosE - e); y = +a*fac*sinE;
	vx = -v*sinE;     vy = +v*fac*cosE;

	// Transformation to reference system (Gaussian vectors) and partials
	Rz(-Omega, Rp);
	Rx(-i, Rq);
	Rz(-omega, Rw);
	matmul("NN", 3, 3, 3, 1.0, Rp, Rq, 0.0, Rpq);
	matmul("NN", 3, 3, 3, 1.0, Rpq, Rw, 0.0, PQW);
	for (k = 0; k < 3; k++) { P[k] = PQW[k]; Q[k] = PQW[3 + k]; W[k] = PQW[6 + k]; }
	e_z[2] = 1;
	cross3(e_z, W, N);
	nN = norm(N, 3);
	for (k = 0; k < 3; k++) { N[k] = N[k] /nN ; }
	cross3(N, P, dPdi), cross3(e_z, P, dPdO);
	cross3(N, Q, dQdi), cross3(e_z, Q, dQdO);
	for (k = 0; k < 3; k++) { dPdo[k] = Q[k]; dQdo[k] = -P[k]; }

	// Partials w.r.t. semimajor axis, eccentricity and mean anomaly at time dt

	for (k = 0; k < 3; k++) {
		dYda[k] = (x / a)*P[k] + (y / a)*Q[k];
		dYda[k+3] = (-vx / (2*a))*P[k] + (-vy / (2*a))*Q[k];
		dYde[k] = (-a-pow(y/fac,2)/r)*P[k] + (x*y/(r*fac*fac))*Q[k];
		dYde[k + 3] = (vx*(2 * a*x + e*pow(y / fac, 2)) / (r*r))*P[k] 
			+ ((n / fac)*pow(a / r, 2)*(x*x / r - pow(y / fac, 2) / a))*Q[k];
		dYdM[k] = (vx*P[k] + vy*Q[k])/n;
		dYdM[k + 3] = (-n*pow(a / r, 3))*(x*P[k] + y*Q[k]);
		dYdi[k] = x*dPdi[k] + y*dQdi[k];
		dYdi[k + 3] = vx*dPdi[k] + vy*dQdi[k];
		dYdO[k] = x*dPdO[k] + y*dQdO[k];
		dYdO[k + 3] = vx*dPdO[k] + vy*dQdO[k];
		dYdo[k] = x*dPdo[k] + y*dQdo[k];
		dYdo[k + 3] = vx*dPdo[k] + vy*dQdo[k];
	}

	// Derivative of mean anomaly at time dt w.r.t. the semimajor axis at epoch
	trace(3, "dYda="); tracemat(3, dYda, 1, 6, 18, 12);
	trace(3, "dYde="); tracemat(3, dYde, 1, 6, 18, 12);
	trace(3, "dYdM="); tracemat(3, dYdM, 1, 6, 18, 12);
	trace(3, "dYdi="); tracemat(3, dYdi, 1, 6, 18, 12);
	trace(3, "dYdO="); tracemat(3, dYdO, 1, 6, 18, 12);
	trace(3, "dYdo="); tracemat(3, dYdo, 1, 6, 18, 12);
	dMda = -1.5*(n / a)*dt;

	// Combined partial derivative matrix of state with respect to epoch elements
	for (k = 0; k<6; k++) {
		dYdA[k+ 0*6] = dYda[k] + dYdM[k]*dMda;
		dYdA[k+ 1 * 6] = dYde[k];
		dYdA[k+ 2 * 6] = dYdi[k];
		dYdA[k+ 3 * 6] = dYdO[k];
		dYdA[k+ 4 * 6] = dYdo[k];
		dYdA[k+ 5 * 6] = dYdM[k];
	}

	return 1;

}

//------------------------------------------------------------------------------
//
// TwoBody
//
// Purpose:
//
//   Propagates a given state vector and computes the state transition matrix 
//   for elliptical Keplerian orbits
//
// Input/Output:
//
//   GM        Gravitational coefficient
//             (gravitational constant * mass of central body)
//   Y0        Epoch state vector (x,y,z,vx,vy,vz)_0
//   dt        Time since epoch
//   Y         State vector (x,y,z,vx,vy,vz)
//   dYdY0     State transition matrix d(x,y,z,vx,vy,vz)/d(x,y,z,vx,vy,vz)_0
//
// Notes:
//
//   The state vector, dt and GM must be given in consistent units, 
//   e.g. [m], [m/s] and [m^3/s^2]. The resulting units of length and velocity  
//   are implied by the units of GM, e.g. [m] and [m/s].
//
//   Due to the internal use of Keplerian elements, the function cannot be 
//   used with epoch state vectors describing a circular or non-inclined orbit.
//
//------------------------------------------------------------------------------

extern void pvpartial(const double * Y0, double dt,double *Y, double *dYdY0)
{

	int     k;
	double  a, e, i, n, sqe2, naa;
	double  P_aM, P_eM, P_eo, P_io, P_iO;
	double A0[6] = { 0.0 };
	double dY0dA0[36] = { 0.0 };
	double dYdA0[36] = { 0.0 };
	double dA0dY0[36] = { 0.0 };
	double GM = GM_EARTH;
	/*Vector  A0(6);
	Matrix  dY0dA0(6, 6), dYdA0(6, 6), dA0dY0(6, 6);*/
	trace(3, "pvpartial:\n");
	//Y0[0] = -2915758.85792398;
	//Y0[1] = 4622150.28654998;
	//Y0[2] = -4581931.3;
	//Y0[3] = -529.964910987115;
	//Y0[4] = -5411.83706208213;
	//Y0[5] = -5119.3057;
	     
	// Orbital elements at epoch
	pv2kepler(GM, Y0, A0);
	//A0 = Elements(GM, Y0);
	a = A0[0];  e = A0[1];  i = A0[2];
	n = sqrt(GM / (a*a*a));

	// Propagated state 
	kepler2pv(GM, A0, dt, Y);
	//Y = State(GM, A0, dt);
	trace(3, "Y0=\t"); tracemat(3, Y0, 1, 6, 13, 4);
	trace(3, "A0=\t"); tracemat(3, A0, 1, 6, 18, 12);
	// State vector partials w.r.t epoch elements
	//dY0dA0 = StatePartials(GM, A0, 0.0);
	//dYdA0 = StatePartials(GM, A0, dt);
	keplerpartials(GM, A0, 0.0, dY0dA0);
	keplerpartials(GM, A0, dt, dYdA0);

	trace(3, "dY0dA0=\n"); tracemat(3, dY0dA0, 6, 6, 18, 10);
	trace(3, "dYdA0=\n"); tracemat(3, dYdA0, 6, 6, 18, 10);
	// Poisson brackets
	sqe2 = sqrt((1.0 - e)*(1.0 + e));
	naa = n*a*a;

	P_aM = -2.0 / (n*a);                   // P(a,M)     = -P(M,a)
	P_eM = -(1.0 - e)*(1.0 + e) / (naa*e);     // P(e,M)     = -P(M,e)
	P_eo = +sqe2 / (naa*e);                // P(e,omega) = -P(omega,e)
	P_io = -1.0 / (naa*sqe2*tan(i));       // P(i,omega) = -P(omega,i)
	P_iO = +1.0 / (naa*sqe2*sin(i));       // P(i,Omega) = -P(Omega,i)

										   // Partials of epoch elements w.r.t. epoch state

	for (k = 0; k < 3; k++) {

		dA0dY0[0 + k * 6] = +P_aM*dY0dA0[k + 3+ 5*6];
		dA0dY0[0+(k + 3)*6] = -P_aM*dY0dA0[k+ 5*6];

		dA0dY0[1+ k*6] = +P_eo*dY0dA0[k + 3+ 4*6] + P_eM*dY0dA0[k + 3+ 5*6];
		dA0dY0[1+ (k + 3) * 6] = -P_eo*dY0dA0[k+4*6] - P_eM*dY0dA0[k+ 5*6];

		dA0dY0[2+ k * 6] = +P_iO*dY0dA0[k + 3+ 3*6] + P_io*dY0dA0[k + 3+ 4*6];
		dA0dY0[2+ (k + 3) * 6] = -P_iO*dY0dA0[k+ 3*6] - P_io*dY0dA0[k+ 4*6];

		dA0dY0[3+ k * 6] = -P_iO*dY0dA0[k + 3+ 2*6];
		dA0dY0[3+ (k + 3) * 6] = +P_iO*dY0dA0[k+ 2*6];

		dA0dY0[4+ k * 6] = -P_eo*dY0dA0[k + 3+ 1*6] - P_io*dY0dA0[k + 3+ 2*6];
		dA0dY0[4+ (k + 3) * 6] = +P_eo*dY0dA0[k+1*6] + P_io*dY0dA0[k+ 2*6];

		dA0dY0[5+ k * 6] = -P_aM*dY0dA0[k + 3+ 0*6] - P_eM*dY0dA0[k + 3+ 1*6];
		dA0dY0[5+ (k + 3) * 6] = +P_aM*dY0dA0[k+ 0*6] + P_eM*dY0dA0[k+ 1*6];

	};

	// State transition matrix
	trace(3, "dA0dY0=\n"); tracemat(3, dA0dY0, 6, 6, 18, 10);
	matmul("NN", 6, 6, 6, 1.0, dYdA0, dA0dY0, 0.0, dYdY0);
	trace(3, "dYdY0=\n"); tracemat(3, dYdY0, 6, 6, 18, 10);
	//dYdY0 = dYdA0 * dA0dY0;
	return;

}
//------------------------------------------------------------------------------
//
// Deriv
//
// Purpose:
// 
//   Computes the derivative of the state vector 
//
// Note:
//
//   pAux is expected to point to a variable of type AuxDataRecord, which is
//   used to communicate with the other program sections and to hold data 
//   between subsequent calls of this function
//
//------------------------------------------------------------------------------

void step(gtime_t tutc,const double *erpv, int n, double *y, double * yp)
{

	// Pointer to auxiliary data record
	// Time
	//double  Mjd_GPS = (*p).Mjd0_GPS + t / 86400.0;
	double r[3] = { 0.0 }, v[3] = { 0.0 }, a[3] = {0.0};
	int i;
	for (i = 0; i < 3; i++) { r[i] = y[i]; v[i] = y[i + 3];}
	forcessimple(tutc, erpv, r, n, a);
	//forces(tutc, erpv,r, v,5,1000,1.3,2.3,a);
	for (i = 0; i < 3; i++) { yp[i] = v[i]; yp[i + 3] = a[i]; }

	//for (i = 0; i < 3; i++) { yp[i] = v[i]; yp[i + 3] = a[i]; }

	//yp = Stack(v, Accel(Mjd_GPS, r, (*p).n_grav));
	 //yp = Stack(v, AccelMain(Mjd_GPS, r, v,5,1000,1.3,2.3));
};

/* numerical integration with rk4--------------------*/
// t:current time
// tutc: time
// erpv: earth rotation parameters
// n: 
extern void rk4(double t,gtime_t tutc, const double *erpv, int n, double *y)
{
	double k1[6], k2[6], k3[6], k4[6], w[6];
	double y1[6] = { 0.0 };
	memcpy(y1, y, 6 * sizeof(double));
	gtime_t tut1;
	int i;
	tut1 = tutc;
	step(tutc, erpv, n, y, k1);
	for (i = 0; i<6; i++) w[i] = y[i] + k1[i] * t / 2.0;
	tut1 = timeadd(tutc, t / 2.0);
	step(tut1, erpv, n, w, k2);
	for (i = 0; i<6; i++) w[i] = y[i] + k2[i] * t / 2.0;
	step(tut1, erpv, n, w, k3);
	for (i = 0; i<6; i++) w[i] = y[i] + k3[i] * t;
	tut1 = timeadd(tutc, t);
	step(tut1, erpv, n, w, k4);
	for (i = 0; i<6; i++) y[i] = y[i] + ((k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])*t / 6.0); 

	trace(3, "rk4: dy=%14.3f %14.3f %14.3f %14.3f %14.3f %14.3f\n", y[0] - y1[0], y[1] - y1[1], y[2] - y1[2], y[3] - y1[3], y[4] - y1[4], y[5] - y1[5]);
}