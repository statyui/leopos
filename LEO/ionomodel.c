#include "leoorb.h"
/* ionosphere model ------------------------------------------------------------
* compute ionospheric delay by broadcast ionosphere model (klobuchar model)
* args   : gtime_t t        I   time (gpst)
*          double *ion      I   iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3}
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
* return : ionospheric delay (L1) (m)
*-----------------------------------------------------------------------------*/
extern double ionmodel(gtime_t t, const double *ion, const double *pos,
	const double *azel)
{
	const double ion_default[] = { /* 2004/1/1 */
		0.1118E-07,-0.7451E-08,-0.5961E-07, 0.1192E-06,
		0.1167E+06,-0.2294E+06,-0.1311E+06, 0.1049E+07
	};
	double tt, f, psi, phi, lam, amp, per, x;
	int week;

	if (pos[2]<-1E3 || azel[1] <= 0) return 0.0;
	if (norm(ion, 8) <= 0.0) ion = ion_default;

	/* earth centered angle (semi-circle) */
	psi = 0.0137 / (azel[1] / PI + 0.11) - 0.022;

	/* subionospheric latitude/longitude (semi-circle) */
	phi = pos[0] / PI + psi*cos(azel[0]);
	if (phi> 0.416) phi = 0.416;
	else if (phi<-0.416) phi = -0.416;
	lam = pos[1] / PI + psi*sin(azel[0]) / cos(phi*PI);

	/* geomagnetic latitude (semi-circle) */
	phi += 0.064*cos((lam - 1.617)*PI);

	/* local time (s) */
	tt = 43200.0*lam + time2gpst(t, &week);
	tt -= floor(tt / 86400.0)*86400.0; /* 0<=tt<86400 */

									   /* slant factor */
	f = 1.0 + 16.0*pow(0.53 - azel[1] / PI, 3.0);

	/* ionospheric delay */
	amp = ion[0] + phi*(ion[1] + phi*(ion[2] + phi*ion[3]));
	per = ion[4] + phi*(ion[5] + phi*(ion[6] + phi*ion[7]));
	amp = amp<    0.0 ? 0.0 : amp;
	per = per<72000.0 ? 72000.0 : per;
	x = 2.0*PI*(tt - 50400.0) / per;

	return CLIGHT*f*(fabs(x)<1.57 ? 5E-9 + amp*(1.0 + x*x*(-0.5 + x*x / 24.0)) : 5E-9);
}
/* ionosphere mapping function -------------------------------------------------
* compute ionospheric delay mapping function by single layer model
* args   : double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
* return : ionospheric mapping function
*-----------------------------------------------------------------------------*/
extern double ionmapf(const double *pos, const double *azel)
{
	if (pos[2] >= HION) return 1.0;
	return 1.0 / cos(asin((RE_WGS84 + pos[2]) / (RE_WGS84 + HION)*sin(PI / 2.0 - azel[1])));
}
/* ionospheric pierce point position -------------------------------------------
* compute ionospheric pierce point (ipp) position and slant factor
* args   : double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          double re        I   earth radius (km)
*          double hion      I   altitude of ionosphere (km)
*          double *posp     O   pierce point position {lat,lon,h} (rad,m)
* return : slant factor
* notes  : see ref [2], only valid on the earth surface
*          fixing bug on ref [2] A.4.4.10.1 A-22,23
*-----------------------------------------------------------------------------*/
extern double ionppp(const double *pos, const double *azel, double re,
	double hion, double *posp)
{
	double cosaz, rp, ap, sinap, tanap;

	rp = re / (re + hion)*cos(azel[1]);
	ap = PI / 2.0 - azel[1] - asin(rp);
	sinap = sin(ap);
	tanap = tan(ap);
	cosaz = cos(azel[0]);
	posp[0] = asin(sin(pos[0])*cos(ap) + cos(pos[0])*sinap*cosaz);

	if ((pos[0]> 70.0*D2R&& tanap*cosaz>tan(PI / 2.0 - pos[0])) ||
		(pos[0]<-70.0*D2R&&-tanap*cosaz>tan(PI / 2.0 + pos[0]))) {
		posp[1] = pos[1] + PI - asin(sinap*sin(azel[0]) / cos(posp[0]));
	}
	else {
		posp[1] = pos[1] + asin(sinap*sin(azel[0]) / cos(posp[0]));
	}
	return 1.0 / sqrt(1.0 - rp*rp);
}