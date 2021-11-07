#include "leoorb.h"

typedef struct ObsType {
	gtime_t  time;    // GPS time
	double  r_obs[3];      // Measured position vector (WGS84, [m])
	double  v_obs[3];      // Measured velocity vector (WGS84, [m/s])
	double  r_ref[3];      // True position vector (WGS84, [m])
	double  v_ref[3];      // True velocity vector (WGS84, [m/s])
}obst;


#define Rz(t,X) do { \
    (X)[8]=1.0; (X)[2]=(X)[5]=(X)[6]=(X)[7]=0.0; \
    (X)[0]=(X)[4]=cos(t); (X)[3]=sin(t); (X)[1]=-(X)[3]; \
} while (0)

int readobs(FILE *fp, obst *obs) {
	char buf[200];
	double ep[6] = { 0.0 };
	gtime_t time;
	
	fgets(buf, 512, fp);
	if (feof(fp)) return 0;
	sscanf(buf, "%lf/%lf/%lf %lf:%lf:%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", ep+0, ep+1, ep+2, ep+3, ep+4, ep+5,
		obs->r_obs, obs->r_obs+1, obs->r_obs+2, obs->v_obs, obs->v_obs+1, obs->v_obs+2, obs->r_ref, obs->r_ref+1, obs->r_ref+2,
		obs->v_ref, obs->v_ref+1, obs->v_ref+2);
	time=epoch2time(ep);
	obs->time = time;
	return 1;

}
void measupd(double y, double g, double sigma, const double *H, double *x, double *P)
{
	double K[6] = { 0.0 }, PH[6] = {0.0};                  // Kalman gain
	double R = sigma*sigma;    // Inverse weight (measurement covariance)
	double HPH;
	double Kp[36] = {0.0};
	double xp[6] = { 0.0 }, Pp[36] = { 0.0 };
	int i;
	trace(3, "measupd:\n");
	trace(3, "H=\t"); tracemat(3, H, 1, 6, 13, 4);
	

	 // Kalman gain
	matmul("NN", 6, 1, 6, 1.0, P, H, 0.0, PH);//PH=PH^T
	trace(3, "PH=\t"); tracemat(3, PH, 1, 6, 13, 4);
	//K = P*G / (Inv_W + Dot(G, P*G));
	HPH = dot(H, PH, 6); //HPH=HPH^T
	trace(3, "HPH=%13.4lf\n",HPH);
	for (i = 0; i < 6; i++) { K[i] = PH[i] / (R + HPH); } //K=PH^T(HPH^T+R)
	trace(3, "K=\t"); tracemat(3, K, 1, 6, 13, 4);
	// State update
	//x = x + K*(z - g);
	for (i = 0; i < 6; i++) { xp[i] = x[i]+K[i]*(y-g); } //xp=x+K(y-Hx)
	trace(3, "y-Ax=%13.4lf  R=%13.4f \n", (y-g),R);
	trace(3, "xp=\t"); tracemat(3, xp, 1, 6, 13, 4);
	// Covariance update
	for (i = 0; i < 6; i++) Kp[i + i * 6] = 1;
	matmul("NN", 6, 6, 1, -1.0, K, H, 1.0, Kp); //Kp=I-KH
	matmul("NN", 6, 6, 6, 1.0, Kp, P, 0.0, Pp); //Pp=KpP
	trace(3, "I-KH=\t"); tracemat(3, Kp, 6, 6, 13, 4);
	trace(3, "Pp=\t"); tracemat(3, Pp, 6, 6, 13, 4);
	memcpy(x, xp, 6 * sizeof(double));
	memcpy(P, Pp, 36 * sizeof(double));
//	P = (Id(n) - Dyadic(K, G))*P;

}

/*int main(int argc, char* argv[]) {*/

int forcastmain(int argc, char* argv[]) {

	// Constants
	double dxdY[6] = { 0.0 };
	double dydY[6] = { 0.0 };
	double dzdY[6] = { 0.0 };
	dxdY[0] = 1; dydY[1] = 1; dzdY[2] = 1;

	// Variables

	int      Reject;                       // Flag for rejected data
	int       i,j;                            // Loop counter
	//configuration
	double    Step=30;                         // Nominal RK4 step size [s]
	double    sigma_xyz=100;                    // Measurement std. dev. [m]
	double    sigma_pos=1000;                    // A priori std. dev. [m]
	double    sigma_vel=10;                    // A priori std. dev. [m/s]
	double    w_pos=0.5;                        // State noise std. dev. [m]
	double    w_vel=0.0005;                        // State noise std. dev. [m/s]
	double    EditLevel=3.0;                    // Data edit level


	double    Sig_a;                        // Standard deviation
	gtime_t    tutc, t_old;
	double Y[6] = { 0.0 }, Y_old[6] = { 0.0 }; // State vector
	double r_est[3] = { 0.0 }, v_est[3] = { 0.0 }; // State vector
	double SigY[6] = { 0.0 }, dadY[6] = { 0.0 }; // Standard deviation,Partials of sma. wrt. state 
	double U[9] = { 0.0 }, dU[9] = { 0.0 }; //Sidereal Time matrix and derivative
	double Phi[36] = { 0.0 }, P[36] = { 0.0 }, Qt[36] = { 0.0 }, Q[36] = { 0.0 }; // State transition, covariance and
	double dr[3] = { 0.0 };
	double    h;                            // Step size 
	FILE *fpobs = NULL;
	FILE *fpsol = NULL;
	char buf[200];
	obst obs;
	double dt;
	double gps_utc = 0.0;
	double erpv[5] = { 0.0 };
	double gmst = 0.0;
	double r[3] = { 0.0 }, v[3] = { 0.0 }, dv[3] = { 0.0 };
	double rr[3] = { 0.0 }, vr[3] = { 0.0 }, dvr[3] = { 0.0 };
	double drinno[3] = { 0.0 };
	int epid = 0;
	erpv[2]=-0.05;

	traceopen("LEOorb.trace");
	tracelevel(4);

	// Initialize UT1-UTC and UTC-TAI time difference
	//IERS::Set(-0.05, -30.00, 0.0, 0.0);

	fpobs = fopen("RTOD.dat", "r");
	fgets(buf, 256, fpobs);
	readobs(fpobs, &obs);
	fpsol = fopen("RTOD.sol", "w");

	// Reference epoch (GPS time)
	//Mjd0 = Obs.Mjd_GPS;
	//Aux.Mjd0_GPS = Obs.Mjd_GPS;
	//t = 0.0;

	//// Transformation from WGS to pseudo-true-of-date system
	//Mjd_UT1 = (IERS::UT1_UTC(Mjd0) - IERS::GPS_UTC(Mjd0));
	//Mjd_UT1 = Mjd0 + (IERS::UT1_UTC(Mjd0)
	//	- IERS::GPS_UTC(Mjd0)) / 86400.0;

	//S(0, 1) = 1.0; S(1, 0) = -1.0;
	//U = R_z(GMST(Mjd_UT1));            // Earth rotation matrix
	//dU = omega_Earth*S*U;               // and derivative [1/s]
	gps_utc = getleaps(obs.time);
	tutc = timeadd(obs.time, gps_utc);
	ecef2qeci(tutc, erpv[2], U, dU);
	
	matmul("TN", 3, 1, 3, 1.0, U, obs.r_obs, 0,  r);
	matmul("TN", 3, 1, 3, 1.0, U, obs.v_obs, 0,  v);
	matmul("TN", 3, 1, 3, 1.0, dU, obs.r_obs, 0,  dv);
	for (i = 0; i < 3; i++) { v[i] = v[i] + dv[i]; }
	//r = Transp(U)*Obs.r_obs;
	//v = Transp(U)*Obs.v_obs;
	//v = Transp(U)*Obs.v_obs + Transp(dU)*Obs.r_obs;

	// A priori covariance
	for (i = 0; i < 3; i++) {
		P[i + i * 6] = sigma_pos*sigma_pos;
		P[i + 3 + (i + 3) * 6] = sigma_vel*sigma_vel;
		Q[i + i * 6] = w_pos*w_pos;
		Q[i + 3 + (i + 3) * 6] = w_vel*w_vel;
	}
	// Initialization
	for (i = 0; i < 3; i++) {
		Y[i] = r[i]; Y[i + 3] = v[i];
	}

	for (i = 1; i <= 2000; i++) {
		epid++;
		if (epid ==239) {
			h = 0;
		}
		
		// Previous step
		t_old = tutc;
		memcpy(Y_old, Y, 6 * sizeof(double));
		// Next observation
		if (!readobs(fpobs, &obs)) break;
		//for (i = 0; i < 3; i++) {
		//	Y[i] = obs.r_obs[i]; Y[i + 3] = obs.v_obs[i];
		//}
		trace(3, "processing epoch id=%d, time=%10.3f\n",epid, fmod(obs.time.time, 86400));
		tutc = timeadd(obs.time, gps_utc);
		dt = timediff(tutc, t_old);
		pvpartial(Y_old, dt, Y, Phi);
		// Propagation to measurement epoch

		// Integration to time of measurement  
		memcpy(Y, Y_old, 6 * sizeof(double));
		while (timediff(tutc,t_old)>0) {
			h = timediff(tutc,t_old);
			if (h>Step) h = Step;
			//rk4(h, t_old, erpv, 20, Y);
			t_old=timeadd(t_old, h);
		};
		tutc = timeadd(obs.time, gps_utc);
		ecef2qeci(tutc, erpv[2], U, dU);
		matmul("TN", 3, 1, 3, 1.0, U, obs.r_obs, 0,  r);
		matmul("TN", 3, 1, 3, 1.0, U, obs.v_obs, 0,  v);
		matmul("TN", 3, 1, 3, 1.0, dU, obs.r_obs, 0,  dv);
		for (i = 0; i < 3; i++) { v[i] = v[i] + dv[i]; }
		trace(2, "Phi=\n"); tracemat(3, Phi, 6, 6, 13, 4);
		matmul("NN", 6, 6, 6, 1.0, Phi, P, 0.0, Qt);
		trace(2, "Qt=\n"); tracemat(3, Qt, 6, 6, 13, 4);
		matmul("NT", 6, 6, 6, 1.0, Qt, Phi, 0.0, P);
		trace(2, "P=\n"); tracemat(3, P, 6, 6, 20, 4);
		for (i = 0; i < 6; i++) {
				P[i + i*6] = P[i + i * 6] + Q[i+i*6];
		}
		for (i = 0; i < 3; i++) {
			drinno[i] = r[i] - Y[i];
		}
		trace(3, "obsdif= %13.4f %13.4f %13.4f\n", drinno[0], drinno[1], drinno[2]);
		if (norm(drinno, 3) < 300) {
			// Estimated state vector in rotating, Earth-fixed system
			measupd(r[0], Y[0], sigma_xyz, dxdY, Y, P);
			measupd(r[1], Y[1], sigma_xyz, dydY, Y, P);
			measupd(r[2], Y[2], sigma_xyz, dzdY, Y, P);
		}
		else {
			trace(3, "outlier in obs detected! norm=%13.4f\n", norm(drinno, 3));
		}
	/*	r_est = U*Filter.State().slice(0, 2);
		v_est = U*Filter.State().slice(3, 5) + dU*Filter.State().slice(0, 2);*/
		matmul("TN", 3, 1, 3, 1.0, U, obs.r_ref, 0,  rr);
		matmul("TN", 3, 1, 3, 1.0, U, obs.v_ref, 0,  vr);
		matmul("TN", 3, 1, 3, 1.0, dU, obs.r_ref, 0,  dvr);
		for (i = 0; i < 3; i++) { vr[i] = vr[i] + dvr[i]; }
		fprintf(fpsol, "%14.3lf %14.3lf %14.3lf %14.3lf\n", fmod(tutc.time, 86400), Y[0] - rr[0], Y[1] - rr[1], Y[2] - rr[2]);

	}

	fclose(fpobs);
	return 0;

}
