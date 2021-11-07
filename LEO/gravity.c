#include"leoorb.h"
//reference
//http://icgem.gfz-potsdam.de/tom_longtime
#define MAXGEOORDER 120
#define NJGM3		20
#define RJGM3		(6378.1363e3)  // Radius Earth [m]; JGM3

// Earth gravity field JGM3
// Gravitational coefficients C, S are efficiently stored in a single
// array CS. The lower triangle matrix CS holds the non-sectorial C
// coefficients C_n,m (n!=m). Sectorial C coefficients C_n,n are the 
// diagonal elements of CS and the upper triangular matrix stores
// the S_n,m (m!=0) coefficients in columns, for the same degree n.
// Mapping of CS to C, S is achieved through 
// C_n,m = CS(n,m), S_n,m = CS(m-1,n)
// Radius Earth [m]; JGM3
const double CS_JGM3[NJGM3 + 1][NJGM3 + 1] = {
	{ 1.000000e+00, 0.000000e+00, 1.543100e-09, 2.680119e-07,-4.494599e-07,
	-8.066346e-08, 2.116466e-08, 6.936989e-08, 4.019978e-08, 1.423657e-08,
	-8.128915e-08,-1.646546e-08,-2.378448e-08, 2.172109e-08, 1.443750e-08,
	4.154186e-09, 1.660440e-08,-1.427822e-08,-1.817656e-08, 7.160542e-11,
	2.759192e-09 },
	{ 0.000000e+00, 0.000000e+00,-9.038681e-07,-2.114024e-07, 1.481555e-07,
	-5.232672e-08,-4.650395e-08, 9.282314e-09, 5.381316e-09,-2.228679e-09,
	-3.057129e-09,-5.097360e-09, 1.416422e-09,-2.545587e-09,-1.089217e-10,
	-1.045474e-09, 7.856272e-10, 2.522818e-10, 3.427413e-10,-1.008909e-10,
	3.216826e-10 },
	{ -1.082627e-03,-2.414000e-10, 1.574536e-06, 1.972013e-07,-1.201129e-08,
	-7.100877e-09, 1.843134e-10,-3.061150e-09,-8.723520e-10,-5.633921e-10,
	-8.989333e-10,-6.863521e-10, 9.154575e-11, 3.005522e-10, 5.182512e-11,
	3.265044e-11,-4.271981e-11, 1.297841e-11,-4.278803e-12,-1.190759e-12,
	3.778260e-11 },
	{ 2.532435e-06, 2.192799e-06, 3.090160e-07, 1.005589e-07, 6.525606e-09,
	3.873005e-10,-1.784491e-09,-2.636182e-10, 9.117736e-11, 1.717309e-11,
	-4.622483e-11,-2.677798e-11, 9.170517e-13,-2.960682e-12,-3.750977e-12,
	1.116419e-12, 5.250141e-12, 2.159727e-12, 1.105860e-13,-3.556436e-13,
	-1.178441e-12 },
	{ 1.619331e-06,-5.087253e-07, 7.841223e-08, 5.921574e-08,-3.982396e-09,
	-1.648204e-09,-4.329182e-10, 6.397253e-12, 1.612521e-11,-5.550919e-12,
	-3.122269e-12, 1.982505e-12, 2.033249e-13, 1.214266e-12,-2.217440e-13,
	8.637823e-14,-1.205563e-14, 2.923804e-14, 1.040715e-13, 9.006136e-14,
	-1.823414e-14 },
	{ 2.277161e-07,-5.371651e-08, 1.055905e-07,-1.492615e-08,-2.297912e-09,
	4.304768e-10,-5.527712e-11, 1.053488e-11, 8.627743e-12, 2.940313e-12,
	-5.515591e-13, 1.346234e-13, 9.335408e-14,-9.061871e-15, 2.365713e-15,
	-2.505252e-14,-1.590014e-14,-9.295650e-15,-3.743268e-15, 3.176649e-15,
	-5.637288e-17 },
	{ -5.396485e-07,-5.987798e-08, 6.012099e-09, 1.182266e-09,-3.264139e-10,
	-2.155771e-10, 2.213693e-12, 4.475983e-13, 3.814766e-13,-1.846792e-13,
	-2.650681e-15,-3.728037e-14, 7.899913e-15,-9.747983e-16,-3.193839e-16,
	2.856094e-16,-2.590259e-16,-1.190467e-16, 8.666599e-17,-8.340023e-17,
	-8.899420e-19 },
	{ 3.513684e-07, 2.051487e-07, 3.284490e-08, 3.528541e-09,-5.851195e-10,
	5.818486e-13,-2.490718e-11, 2.559078e-14, 1.535338e-13,-9.856184e-16,
	-1.052843e-14, 1.170448e-15, 3.701523e-16,-1.095673e-16,-9.074974e-17,
	7.742869e-17, 1.086771e-17, 4.812890e-18, 2.015619e-18,-5.594661e-18,
	1.459810e-18 },
	{ 2.025187e-07, 1.603459e-08, 6.576542e-09,-1.946358e-10,-3.189358e-10,
	-4.615173e-12,-1.839364e-12, 3.429762e-13,-1.580332e-13, 7.441039e-15,
	-7.011948e-16, 2.585245e-16, 6.136644e-17, 4.870630e-17, 1.489060e-17,
	1.015964e-17,-5.700075e-18,-2.391386e-18, 1.794927e-18, 1.965726e-19,
	-1.128428e-19 },
	{ 1.193687e-07, 9.241927e-08, 1.566874e-09,-1.217275e-09,-7.018561e-12,
	-1.669737e-12, 8.296725e-13,-2.251973e-13, 6.144394e-14,-3.676763e-15,
	-9.892610e-17,-1.736649e-17, 9.242424e-18,-4.153238e-18,-6.937464e-20,
	3.275583e-19, 1.309613e-19, 1.026767e-19,-1.437566e-20,-1.268576e-20,
	-6.100911e-21 },
	{ 2.480569e-07, 5.175579e-08,-5.562846e-09,-4.195999e-11,-4.967025e-11,
	-3.074283e-12,-2.597232e-13, 6.909154e-15, 4.635314e-15, 2.330148e-15,
	4.170802e-16,-1.407856e-17,-2.790078e-19,-6.376262e-20,-1.849098e-19,
	3.595115e-20,-2.537013e-21, 4.480853e-21, 4.348241e-22, 1.197796e-21,
	-1.138734e-21 },
	{ -2.405652e-07, 9.508428e-09, 9.542030e-10,-1.409608e-10,-1.685257e-11,
	1.489441e-12,-5.754671e-15, 1.954262e-15,-2.924949e-16,-1.934320e-16,
	-4.946396e-17, 9.351706e-18,-9.838299e-20, 1.643922e-19,-1.658377e-20,
	2.905537e-21, 4.983891e-22, 6.393876e-22,-2.294907e-22, 6.437043e-23,
	6.435154e-23 },
	{ 1.819117e-07,-3.068001e-08, 6.380398e-10, 1.451918e-10,-2.123815e-11,
	8.279902e-13, 7.883091e-15,-4.131557e-15,-5.708254e-16, 1.012728e-16,
	-1.840173e-18, 4.978700e-19,-2.108949e-20, 2.503221e-20, 3.298844e-21,
	-8.660491e-23, 6.651727e-24, 5.110031e-23,-3.635064e-23,-1.311958e-23,
	1.534228e-24 },
	{ 2.075677e-07,-2.885131e-08, 2.275183e-09,-6.676768e-11,-3.452537e-13,
	1.074251e-12,-5.281862e-14, 3.421269e-16,-1.113494e-16, 2.658019e-17,
	4.577888e-18,-5.902637e-19,-5.860603e-20,-2.239852e-20,-6.914977e-23,
	-6.472496e-23,-2.741331e-23, 2.570941e-24,-1.074458e-24,-4.305386e-25,
	-2.046569e-25 },
	{ -1.174174e-07,-9.997710e-09,-1.347496e-09, 9.391106e-11, 3.104170e-13,
	3.932888e-13,-1.902110e-14, 2.787457e-15,-2.125248e-16, 1.679922e-17,
	1.839624e-18, 7.273780e-20, 4.561174e-21, 2.347631e-21,-7.142240e-22,
	-2.274403e-24,-2.929523e-24, 1.242605e-25,-1.447976e-25,-3.551992e-26,
	-7.473051e-28 },
	{ 1.762727e-08, 6.108862e-09,-7.164511e-10, 1.128627e-10,-6.013879e-12,
	1.293499e-13, 2.220625e-14, 2.825477e-15,-1.112172e-16, 3.494173e-18,
	2.258283e-19,-1.828153e-21,-6.049406e-21,-5.705023e-22, 1.404654e-23,
	-9.295855e-24, 5.687404e-26, 1.057368e-26, 4.931703e-27,-1.480665e-27,
	2.400400e-29 },
	{ -3.119431e-08, 1.356279e-08,-6.713707e-10,-6.451812e-11, 4.698674e-12,
	-9.690791e-14, 6.610666e-15,-2.378057e-16,-4.460480e-17,-3.335458e-18,
	-1.316568e-19, 1.643081e-20, 1.419788e-21, 9.260416e-23,-1.349210e-23,
	-1.295522e-24,-5.943715e-25,-9.608698e-27, 3.816913e-28,-3.102988e-28,
	-8.192994e-29 },
	{ 1.071306e-07,-1.262144e-08,-4.767231e-10, 1.175560e-11, 6.946241e-13,
	-9.316733e-14,-4.427290e-15, 4.858365e-16, 4.814810e-17, 2.752709e-19,
	-2.449926e-20,-6.393665e-21, 8.842755e-22, 4.178428e-23,-3.177778e-24,
	1.229862e-25,-8.535124e-26,-1.658684e-26,-1.524672e-28,-2.246909e-29,
	-5.508346e-31 },
	{ 4.421672e-08, 1.958333e-09, 3.236166e-10,-5.174199e-12, 4.022242e-12,
	3.088082e-14, 3.197551e-15, 9.009281e-17, 2.534982e-17,-9.526323e-19,
	1.741250e-20,-1.569624e-21,-4.195542e-22,-6.629972e-24,-6.574751e-25,
	-2.898577e-25, 7.555273e-27, 3.046776e-28, 3.696154e-29, 1.845778e-30,
	6.948820e-31 },
	{ -2.197334e-08,-3.156695e-09, 7.325272e-10,-1.192913e-11, 9.941288e-13,
	3.991921e-14,-4.220405e-16, 7.091584e-17, 1.660451e-17, 9.233532e-20,
	-5.971908e-20, 1.750987e-21,-2.066463e-23,-3.440194e-24,-1.487095e-25,
	-4.491878e-26,-4.558801e-27, 5.960375e-28, 8.263952e-29,-9.155723e-31,
	-1.237749e-31 },
	{ 1.203146e-07, 3.688524e-09, 4.328972e-10,-6.303973e-12, 2.869669e-13,
	-3.011115e-14, 1.539793e-15,-1.390222e-16, 1.766707e-18, 3.471731e-19,
	-3.447438e-20, 8.760347e-22,-2.271884e-23, 5.960951e-24, 1.682025e-25,
	-2.520877e-26,-8.774566e-28, 2.651434e-29, 8.352807e-30,-1.878413e-31,
	4.054696e-32 }
};
/* read gravity parameters ----------------------------------------------
* read gravity parameters
* args   : char   *file       I   IGS ERP file (IGS ERP ver.2)
*          erp_t  *erp        O   earth rotation parameters
* return : status (1:ok,0:file open error)
*-----------------------------------------------------------------------------*/
extern int readgravity(const char *file, int type, geo_t *geo)
{
	FILE *fp;
	geod_t *geo_data;
	double v[4] = { 0 };
	int id[2] = { 0 };
	char buff[256];

	trace(3, "readgravity: file=%s\n", file);

	if (!(fp = fopen(file, "r"))) {
		trace(2, "gravity model file open error: file=%s\n", file);
		return 0;
	}
	fgets(buff, sizeof(buff), fp);
		if (type == GEOMODEL_JGM3) {
			geo->gmconst = 398600.44150E+09;
			geo->radius = 6378136.30;
			geo->maxorder = 70;
			if (!strstr(buff, "JGM3")) {
				trace(2, "gravity model mispatch!\n");
				return 0;
			}
		}
		else if (type == GEOMODEL_EGM96) {
			geo->gmconst = 398600.44150E+09;
			geo->radius = 6378136.30;
			geo->maxorder = 120;
			if (!strstr(buff, "EGM96")) {
				trace(2, "gravity model mispatch!\n");
				return 0;
			}
		}
		else if (type == GEOMODEL_EGM2008) {
			geo->gmconst = 398600.44150E+09;
			geo->radius = 6378136.30;
			geo->maxorder = 120;
			if (!strstr(buff, "EGM2008")) {
				trace(2, "gravity model mispatch!\n");
				return 0;
			}
		}
		else if (type == GEOMODEL_EIGEN2) {
			geo->gmconst = 398600.44150E+09;
			geo->radius = 6378136.46;
			geo->maxorder = 120;
			if (!strstr(buff, "EIGEN2")) {
				trace(2, "gravity model mispatch!\n");
				return 0;
			}
		}
		else if (type == GEOMODEL_EIGEN_GL04C) {
			geo->gmconst = 398600.44150E+09;
			geo->radius = 6378136.46;
			geo->maxorder = 120;
			if (!strstr(buff, "EIGEN_GL04C")) {
				trace(2, "gravity model mispatch!\n");
				return 0;
			}
		}
	while (fgets(buff, sizeof(buff), fp)) {
		if (!strstr(buff, "gfc")) continue; //skip header
		if (sscanf(buff+4, "%2d %2d %lf %lf %lf %lf",
				id, id + 1, v, v + 1, v + 2, v + 3) < 6) {
				continue;
		}
		if (id[0] <= MAXGM&&id[1] <= MAXGM)
		{
			geo->C[id[0]][id[1]] = v[0];
			geo->S[id[0]][id[1]] = v[1];
		}
		//if (geo->n >= geo->nmax) {
		//	geo->nmax = geo->nmax <= 0 ? 128 : geo->nmax * 2;
		//	geo_data = (geod_t *)realloc(geo->data, sizeof(geod_t)*geo->nmax);
		//	if (!geo_data) {
		//		free(geo->data); geo->data = NULL; geo->n = geo->nmax = 0;
		//		fclose(fp);
		//		return 0;
		//	}
		//	geo->data = geo_data;
		//}

		//geo->data[geo->n].m = id[0];
		//geo->data[geo->n].n = id[1] ;
		//geo->data[geo->n].C = v[0] ;
		//geo->data[geo->n].S = v[1];
	}
	fclose(fp);
	return 1;
}
//------------------------------------------------------------------------------
//
// gravityJGM3
//
// Purpose:
//
//   Computes the acceleration due to the harmonic gravity field of the 
//   central body
//
// Input/Output:
//
//   r           Satellite position vector in the inertial system
//   E           Transformation matrix to body-fixed system
//   GM          Gravitational coefficient
//   R_ref       Reference radius 
//   CS          Spherical harmonics coefficients (un-normalized)
//   n_max       Maximum degree 
//   m_max       Maximum order (m_max<=n_max; m_max=0 for zonals, only)
//   <return>    Acceleration (a=d^2r/dt^2)
//
//------------------------------------------------------------------------------
extern double gravityJGM3(const double *r, double *U, double nmax, double mmax, double *a)
{

	// Local variables
	int     n, m;                           // Loop counters
	double  r_sqr, rho, Fac;               // Auxiliary quantities
	double  x0, y0, z0;                      // Normalized coordinates
	double  ax, ay, az;									 //double  ax, ay, az;                      // Acceleration vector 
	double  C, S;                           // Gravitational coefficients
	double rbf[3] = { 0.0 };  // Body-fixed position
	double abf[3] = { 0.0 };  // Body-fixed acceleration                      
	double V[(NJGM3 + 2)][(NJGM3 + 2)] = { 0.0 };// Harmonic functions
	double W[(NJGM3 + 2)][(NJGM3 + 2)] = { 0.0 };  // work array (0..n_max+1,0..n_max+1)
	int k = 0;

	// Body-fixed position 
	//r_bf = E * r;
	matmul("NN", 3, 1, 3, 1.0, U, r, 0.0, rbf);
	//for (k = 0; k < 3; k++) rbf[k] = r[k];

	// Auxiliary quantities

	r_sqr = dot(rbf, rbf, 3);               // Square of distance
	rho = RJGM3*RJGM3 / r_sqr;

	x0 = RJGM3 * rbf[0] / r_sqr;          // Normalized
	y0 = RJGM3 * rbf[1] / r_sqr;          // coordinates
	z0 = RJGM3 * rbf[2] / r_sqr;
	//
	// Evaluate harmonic functions 
	//   V_nm = (R_ref/r)^(n+1) * P_nm(sin(phi)) * cos(m*lambda)
	// and 
	//   W_nm = (R_ref/r)^(n+1) * P_nm(sin(phi)) * sin(m*lambda)
	// up to degree and order n_max+1
	//
	// Calculate zonal terms V(n,0); set W(n,0)=0.0
	V[0][0] = RJGM3 / sqrt(r_sqr);
	W[0][0] = 0.0;
	V[1][0] = z0 * V[0][0];
	W[1][0] = 0.0;
	for (n = 2; n <= nmax + 1; n++) {
		V[n][0] = ((2 * n - 1) * z0 * V[n - 1][0] - (n - 1) * rho * V[n - 2][0]) / n;
		W[n][0] = 0.0;
	};

	// Calculate tesseral and sectorial terms 
	for (m = 1; m <= mmax + 1; m++) {
		// Calculate V(m,m) .. V(n_max+1,m)
		V[m][m] = (2 * m - 1) * (x0*V[m - 1][m - 1] - y0*W[m - 1][m - 1]);
		W[m][m] = (2 * m - 1) * (x0*W[m - 1][m - 1] + y0*V[m - 1][m - 1]);
		if (m <= nmax) {
			V[m + 1][m] = (2 * m + 1) * z0 * V[m][m];
			W[m + 1][m] = (2 * m + 1) * z0 * W[m][m];
		};
		for (n = m + 2; n <= nmax + 1; n++) {
			V[n][m] = ((2 * n - 1)*z0*V[n - 1][m] - (n + m - 1)*rho*V[n - 2][m]) / (n - m);
			W[n][m] = ((2 * n - 1)*z0*W[n - 1][m] - (n + m - 1)*rho*W[n - 2][m]) / (n - m);
		};
	};

	//
	// Calculate accelerations ax,ay,az
	//

	ax = ay = az = 0.0;
	for (m = 0; m <= mmax; m++) {
		for (n = m; n <= nmax; n++) {
			if (m == 0) {
				C = CS_JGM3[n][0];   // = C_n,0
				ax -= C * V[n + 1][1];
				ay -= C * W[n + 1][1];
				az -= (n + 1)*C * V[n + 1][0];
			}
			else {
				C = CS_JGM3[n][m];   // = C_n,m 
				S = CS_JGM3[m - 1][n]; // = S_n,m 
				Fac = 0.5 * (n - m + 1) * (n - m + 2);
				ax += +0.5 * (-C * V[n + 1][m + 1] - S * W[n + 1][m + 1])
					+ Fac * (+C * V[n + 1][m - 1] + S * W[n + 1][m - 1]);
				ay += +0.5 * (-C * W[n + 1][m + 1] + S * V[n + 1][m + 1])
					+ Fac * (-C * W[n + 1][m - 1] + S * V[n + 1][m - 1]);
				az += (n - m + 1) * (-C * V[n + 1][m] - S * W[n + 1][m]);
			};
		}
	}

	
	//for (m = 0; m < 3; m++) {
		abf[0] = (GM_EARTH / (RJGM3*RJGM3))*ax;
		abf[1] = (GM_EARTH / (RJGM3*RJGM3))*ay;
		abf[2] = (GM_EARTH / (RJGM3*RJGM3))*az;
	//}
	matmul("TN", 3, 1, 3, 1.0, U, abf, 0.0, a);
	//a_bf = (GM_EARTH / (RJGM3*RJGM3)) * Vector(ax, ay, az);
	// Inertial acceleration 
	//return  Transp(E)*a_bf;
	return 1;

}
extern double gravityforce(const double *r, double *U, double nm,geo_t *geo, double *a)
{

	// Local variables
	int     i,n, m;                           // Loop counters
	double  r_sqr, rho, Fac;               // Auxiliary quantities
	double  x0, y0, z0;                      // Normalized coordinates
											 //double  ax, ay, az;                      // Acceleration vector 
	double  C, S;                           // Gravitational coefficients
	double rbf[3] = { 0.0 };  // Body-fixed position
	double abf[3] = { 0.0 };  // Body-fixed acceleration                      
	double V[(NJGM3 + 2)][(NJGM3 + 2)] = { 0.0 };// Harmonic functions
	double W[(NJGM3 + 2)][(NJGM3 + 2)] = { 0.0 };  // work array (0..n_max+1,0..n_max+1)

												   // Body-fixed position 
												   //r_bf = E * r;
	matmul("NN", 3, 1, 3, 1.0, U, r, 0.0, rbf);

	// Auxiliary quantities

	r_sqr = dot(rbf, rbf, 3);               // Square of distance
	rho = RJGM3*RJGM3 / r_sqr;

	x0 = RJGM3 * rbf[0] / r_sqr;          // Normalized
	y0 = RJGM3 * rbf[1] / r_sqr;          // coordinates
	z0 = RJGM3 * rbf[2] / r_sqr;
	//
	// Evaluate harmonic functions 
	//   V_nm = (R_ref/r)^(n+1) * P_nm(sin(phi)) * cos(m*lambda)
	// and 
	//   W_nm = (R_ref/r)^(n+1) * P_nm(sin(phi)) * sin(m*lambda)
	// up to degree and order n_max+1
	//
	// Calculate zonal terms V(n,0); set W(n,0)=0.0
	V[0][0] = RJGM3 / sqrt(r_sqr);
	W[0][0] = 0.0;
	V[1][0] = z0 * V[0][0];
	W[1][0] = 0.0;
	for (n = 2; n <= n + 1; n++) {
		V[n][0] = ((2 * n - 1) * z0 * V[n - 1][0] - (n - 1) * rho * V[n - 2][0]) / n;
		W[n][0] = 0.0;
	};

	// Calculate tesseral and sectorial terms 
	for (m = 1; m <= nm + 1; m++) {
		// Calculate V(m,m) .. V(n_max+1,m)
		V[m][m] = (2 * m - 1) * (x0*V[m - 1][m - 1] - y0*W[m - 1][m - 1]);
		W[m][m] = (2 * m - 1) * (x0*W[m - 1][m - 1] + y0*V[m - 1][m - 1]);
		if (m <= nm) {
			V[m + 1][m] = (2 * m + 1) * z0 * V[m][m];
			W[m + 1][m] = (2 * m + 1) * z0 * W[m][m];
		};
		for (n = m + 2; n <= nm + 1; n++) {
			V[n][m] = ((2 * n - 1)*z0*V[n - 1][m] - (n + m - 1)*rho*V[n - 2][m]) / (n - m);
			W[n][m] = ((2 * n - 1)*z0*W[n - 1][m] - (n + m - 1)*rho*W[n - 2][m]) / (n - m);
		};
	};

	//
	// Calculate accelerations ax,ay,az
	//
	//ax = ay = az = 0.0;
	for (i = 0; i <= geo->n; i++) {
		m = geo->n;
		n = geo->n;
		if (m > nm || n > nm) continue;
		C = geo->C[n][m];
			if (m == 0) {
				  // = C_n,0
				abf[0] -= C * V[n + 1][1];
				abf[1] -= C * W[n + 1][1];
				abf[2] -= (n + 1)*C * V[n + 1][0];
			}
			else {
				// = C_n,m 
				S = geo->S[n][m]; // = S_n,m 
				Fac = 0.5 * (n - m + 1) * (n - m + 2);
				abf[0] += +0.5 * (-C * V[n + 1][m + 1] - S * W[n + 1][m + 1])
					+ Fac * (+C * V[n + 1][m - 1] + S * W[n + 1][m - 1]);
				abf[1] += +0.5 * (-C * W[n + 1][m + 1] + S * V[n + 1][m + 1])
					+ Fac * (-C * W[n + 1][m - 1] + S * V[n + 1][m - 1]);
				abf[2] += (n - m + 1) * (-C * V[n + 1][m] - S * W[n + 1][m]);
			};
	}

	// Body-fixed acceleration
	for (m = 0; m < 3; m++) {
		abf[m] = (GM_EARTH / (RJGM3*RJGM3))*abf[m];
	}
	matmul("TN", 3, 1, 3, 1.0, U, abf, 0.0, a);
	//a_bf = (GM_EARTH / (RJGM3*RJGM3)) * Vector(ax, ay, az);
	// Inertial acceleration 
	//return  Transp(E)*a_bf;
	return 1;

}
