#include "leoorb.h"
/* read attitude file (ascii)----------------------------------------------
* read attitude file
* args   : char   *file       I   IGS ERP file (IGS ERP ver.2)
*          erp_t  *erp        O   earth rotation parameters
* return : status (1:ok,0:file open error)
*-----------------------------------------------------------------------------*/
#define NMAX 4
#define MAXATTGAP 15
extern int readatta(const char *file, att_t *att)
{
	FILE *fp;
	attd_t *att_data;
	double v[5] = { 0 };
	char buff[256];
	int timei, flg, hflg = 1;
	char id, tag;
	double ep[6] = { 2000,1,1,12,0,0 };

	trace(3, "readatt: file=%s\n", file);

	if (!(fp = fopen(file, "r"))) {
		trace(2, "att file open error: file=%s\n", file);
		return 0;
	}

	while (fgets(buff, sizeof(buff), fp)) {
		if (sscanf(buff, "%d %c %c %lf %lf %lf %lf %lf %d",
			&timei, &id, &tag, v, v + 1, v + 2, v + 3, v + 4,&flg)<5) {
			continue;
		}
		if (att->n >= att->nmax) {
			att->nmax = att->nmax <= 0 ? 128 : att->nmax * 2;
			att_data = (attd_t *)realloc(att->data, sizeof(attd_t)*att->nmax);
			if (!att_data) {
				free(att->data); att->data = NULL; att->n = att->nmax = 0;
				fclose(fp);
				return 0;
			}
			att->data = att_data;
		}
		att->data[att->n].time = timeadd(epoch2time(ep),timei);
		memcpy(att->data[att->n].quat,v,4*sizeof(double));
		att->n++;
	}
	fclose(fp);
	return 1;
}

/* read attitude file (binary)----------------------------------------------
* read attitude file
* args   : char   *file       I   att file (IGS ERP ver.2)
*          att_t  *att        O   attitude file
* return : status (1:ok,0:file open error)
*-----------------------------------------------------------------------------*/
extern int readattb(const char *file, att_t *att)
{
	FILE *fp;
	attd_t *att_data;
	double v[14] = { 0 };
	char buff[256];
	double ep[6] = { 2000,1,1,12,0,0 };

	trace(3, "readatt: file=%s\n", file);

	if (!(fp = fopen(file, "r"))) {
		trace(2, "att file open error: file=%s\n", file);
		return 0;
	}
	while (fgets(buff, sizeof(buff), fp)) {
		if (sscanf(buff, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
			v, v + 1, v + 2, v + 3, v + 4, v + 5, v + 6, v + 7, v + 8, v + 9, v + 10, v + 11, v + 12, v + 13)<5) {
			continue;
		}
		if (att->n >= att->nmax) {
			att->nmax = att->nmax <= 0 ? 128 : att->nmax * 2;
			att_data = (attd_t *)realloc(att->data, sizeof(attd_t)*att->nmax);
			if (!att_data) {
				free(att->data); att->data = NULL; att->n = att->nmax = 0;
				fclose(fp);
				return 0;
			}
			att->data = att_data;
		}
		att->data[att->n].time = timeadd(epoch2time(ep),v[0]);
		memcpy(att->data[att->n].quat, v+1, 4 * sizeof(double));
		att->n++;
	}
	fclose(fp);
	return 1;
}

/* satellite position by precise ephemeris -----------------------------------*/
static int searchatt(gtime_t time, const att_t *att, double *qt)
{
	double t[NMAX + 1], p[4][NMAX + 1],  std = 0.0;
	int i, j, k, index;

	trace(4, "searchatt : time=%s\n", time_str(time, 3));

	qt[0] = qt[1] = qt[2] =qt[3] = 0.0;
    // 15 minutes calculation
	if (att->n<NMAX + 1 ||
		timediff(time, att->data[0].time)<-MAXATTGAP ||
		timediff(time, att->data[att->n- 1].time)>MAXATTGAP) {
		trace(3, "no valid attitude data %s \n", time_str(time, 0));
		return 0;
	}
	/* binary search */
	for (i = 0, j = att->n - 1; i<j;) {
		k = (i + j) / 2;
		if (timediff(att->data[k].time, time)<0.0) i = k + 1; else j = k;
	}
	index = i <= 0 ? 0 : i - 1;
	

	/* polynomial interpolation for orbit */
	i = index - (NMAX + 1) / 2;
	if (i<0) i = 0; else if (i + NMAX >= att->n) i = att->n - NMAX - 1;
	for (j = 0; j <= NMAX; j++) {
		t[j] = timediff(att->data[i + j].time, time);
	}
	for (j = 0; j <= NMAX; j++) {
		p[0][j] = att->data[i + j].quat[0];
		p[1][j] = att->data[i + j].quat[1];
		p[2][j] = att->data[i + j].quat[2];
		p[3][j] = att->data[i + j].quat[3];
	}
	for (i = 0; i<4; i++) {
		qt[i] = interppol(t, p[i], NMAX + 1);
	}

	return 1;
}
//static int qt2mat(double *quat, double *rm) {
//	//rm[0 + 0 * 3] = quat[0] * quat[0] + quat[1] * quat[1] - quat[2] * quat[2] - quat[3]* quat[3];
//	//rm[0 + 1 * 3] = 2 * (quat[1] * quat[2] + quat[0] * quat[3]);
//	//rm[0 + 2 * 3] = 2 * (quat[1] * quat[3] - quat[0] * quat[2]);
//
//	//rm[1 + 0 * 3] = 2 * (quat[1] * quat[2] - quat[0] * quat[3]);
//	//rm[1 + 1 * 3] = quat[0] * quat[0] - quat[1] * quat[1] + quat[2] * quat[2] - quat[3] * quat[3];
//	//rm[1 + 2 * 3] = 2 * (quat[2] * quat[3] + quat[0] * quat[1]);
//
//	//rm[2 + 0 * 3] = 2 * (quat[1] * quat[3] + quat[0] * quat[2]);
//	//rm[2 + 1 * 3] = 2 * (quat[2] * quat[3] - quat[0] * quat[1]);
//	//rm[2 + 2 * 3] = quat[0] * quat[0] - quat[1] * quat[1] - quat[2] * quat[2] + quat[3] * quat[3];
//
//	rm[0 + 0 * 3] = quat[0] * quat[0] - quat[1] * quat[1] - quat[2] * quat[2] + quat[3] * quat[3];
//	rm[0 + 1 * 3] = 2 * (quat[0] * quat[1] - quat[2] * quat[3]);
//	rm[0 + 2 * 3] = 2 * (quat[0] * quat[2] + quat[1] * quat[3]);
//
//	rm[1 + 0 * 3] = 2 * (quat[0] * quat[1] + quat[2] * quat[3]);
//	rm[1 + 1 * 3] = -quat[0] * quat[0] + quat[1] * quat[1] - quat[2] * quat[2] + quat[3] * quat[3];
//	rm[1 + 2 * 3] = 2 * (quat[1] * quat[2] - quat[0] * quat[3]);
//
//	rm[2 + 0 * 3] = 2 * (quat[0] * quat[2] - quat[1] * quat[3]);
//	rm[2 + 1 * 3] = 2 * (quat[1] * quat[2] + quat[0] * quat[3]);
//	rm[2 + 2 * 3] = -quat[0] * quat[0] - quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3];
//	return 1;
//}
static int qt2mat(double *quat1, double *rm) {
	int i;
	double absq = 0.0;
	double quat[4] = { 0.0 };
	quat[0] = quat1[3];
	for (i = 1; i < 4; i++)quat[i] = quat1[i - 1];
	absq = sqrt(quat[0] * quat[0] + quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]);
	rm[0 + 0 * 3] = quat[0] * quat[0] - quat[1] * quat[1] - quat[2] * quat[2] + quat[3] * quat[3];
	rm[0 + 1 * 3] = 2 * (quat[0] * quat[1] - quat[2] * quat[3]);
	rm[0 + 2 * 3] = 2 * (quat[0] * quat[2] + quat[1] * quat[3]);

	rm[1 + 0 * 3] = 2 * (quat[0] * quat[1] + quat[2] * quat[3]);
	rm[1 + 1 * 3] = -quat[0] * quat[0] + quat[1] * quat[1] - quat[2] * quat[2] + quat[3] * quat[3];
	rm[1 + 2 * 3] = 2 * (quat[1] * quat[2] - quat[0] * quat[3]);

	rm[2 + 0 * 3] = 2 * (quat[0] * quat[2] - quat[1] * quat[3]);
	rm[2 + 1 * 3] = 2 * (quat[1] * quat[2] + quat[0] * quat[3]);
	rm[2 + 2 * 3] = -quat[0] * quat[0] - quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3];
	for (i = 0; i < 9; i++) {
		rm[i] = rm[i] / absq;
	}
	return 1;
}

/* satellite antenna phase center offset ---------------------------------------
* compute satellite antenna phase center offset in ecef
* args   : gtime_t time       I   time (gpst)
*          double *rs         I   satellite position and velocity (ecef)
*                                 {x,y,z,vx,vy,vz} (m|m/s)
*          int    sat         I   satellite number
*          nav_t  *nav        I   navigation data
*          double *dant       I   satellite antenna phase center offset (ecef)
*                                 {dx,dy,dz} (m) (iono-free LC value)
* return : none
*-----------------------------------------------------------------------------*/
extern void leoantoff(gtime_t time, const double *rs,const nav_t *nav, double *dant)
{
	//const double *lam = nav->lam[sat - 1];
	//const pcv_t *pcv = nav->pcvs + sat - 1;
	//double antoff[3] = { -0.0004, -0.0004, -0.4516 };//GRACE
	double antoff[3] = { -0.0004, -0.0004, -0.4516 };//GRACE
	//double antoff[3] = { -0.0000, -0.0000, -0.0000 };
	double ex[3], ey[3], ez[3], es[3], r[3], rsun[3], gmst, erpv[5] = { 0 };
	double b2i[9] = { 0.0 }, U[9] = { 0.0 }, R[9] = {0.0},qt[4] = { 0.0 };
	gtime_t tutc = gpst2utc(time);
	double dant1[3] = { 0.0 };
	double ri[3] = { 0.0 }, vi[3] = {0.0};
	double dU[9] = { 0.0 };
	int i, j = 0, k = 1;
	//if (&nav->att == NULL) {
		/* sun position in ecef */
		eci2ecef(tutc,erpv,U,&gmst);
		matmul("TN", 3, 1, 3, 1.0, U, rs, 0.0, ri);
		//if (norm(rs + 3, 3) < 1) {
			sunmoonpos(tutc, erpv, rsun, NULL, &gmst);
			for (i = 0; i < 3; i++) r[i] = rsun[i] - rs[i];
			if (!normv3(r, es)) return;
		//}
		/*else {
			sint = U[0 + 1 * 3];
			cost = -U[0 + 0 * 3];
			dU[0 + 0 * 3] = -sint;
			dU[0 + 1 * 3] = cost;
			dU[1 + 0 * 3] = -cost;
			dU[1 + 1 * 3] = -sint;
			matmul("TN", 3, 1, 3, 1.0, U, rs + 3, 0.0, vi);
			matmul("NN", 3, 1, 3, OMGE, dU, ri, 1.0, vi);
			for (i = 0; i < 3; i++) r[i] = vi[i];
			if (!normv3(r, es)) return;
		}*/

		/* unit vectors of satellite fixed coordinates */
		for (i = 0; i < 3; i++) r[i] = -ri[i];
		if (!normv3(r, ez)) return;

		cross3(ez, es, r);
		if (!normv3(r, ey)) return;
		cross3(ey, ez, ex);
		for (i = 0; i < 3; i++) {
			dant1[i] = antoff[0] * ex[i] + antoff[1] * ey[i] + antoff[2] * ez[i];
		}
	//}
	//else {
		searchatt(time, &nav->att,qt);/* attitude of LEO satellite*/
		qt2mat(qt, b2i);
		eci2ecef(tutc, erpv, U, &gmst);
		matmul("NN", 3, 3, 3, 1.0,U, b2i, 0.0, R);
		matmul("NN", 3, 1, 3, -1.0, R, antoff, 0.0, dant);

	//}

	trace(4, "leoantoff: time=%s dant=%8.4f %8.4f %8.4f dantn=%8.4f %8.4f %8.4f\n", time_str(time, 3), dant[0], dant[1], dant[2], dant1[0], dant1[1], dant1[2]);
}