#include <stdarg.h>
#include "leoorb.h"

#define SQRT(x)     ((x)<=0.0?0.0:sqrt(x))

#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics?6:3)
#define NL(opt)     ((opt)->mode>PMODE_PPP_KINEMA?0:0) /*LEO CR CD parameter*/
#define NC(opt)     (NSYS)
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3))
#define NI(opt)     ((opt)->ionoopt==IONOOPT_EST?MAXSAT:0)
#define ND(opt)     ((opt)->nf>=3?1:0)
#define NR(opt)     (NP(opt)+NL(opt)+NC(opt)+NT(opt)+NI(opt)+ND(opt))
#define NB(opt)     (NF(opt)*MAXSAT)
#define NX(opt)     (NR(opt)+NB(opt))
/* state variable index */
#define IC(s,opt)   (NP(opt)+(s))
#define IT(opt)     (NP(opt)+NC(opt))
#define II(s,opt)   (NP(opt)+NC(opt)+NT(opt)+(s)-1)
#define ID(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt))
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1)

/* global variables ----------------------------------------------------------*/
static int statlevel = 0;          /* leo status output level (0:off) */
static FILE *fp_stat = NULL;       /* leo status file pointer */
static char file_stat[1024] = "";  /* leo status file original path */
static gtime_t time_stat = { 0 };    /* leo status file time */

									 /* standard deviation of state -----------------------------------------------*/
static double STD(leo_t *leo, int i)
{
	if (leo->sol.stat == SOLQ_FIX) return SQRT(leo->Pa[i + i*leo->nx]);
	return SQRT(leo->P[i + i*leo->nx]);
}

extern int openstat(const char *file, int level)
{
	gtime_t time = utc2gpst(timeget());
	char path[1024];

	trace(3, "openstat: file=%s level=%d\n", file, level);

	if (level <= 0) return 0;

	reppath(file, path, time, "", "");

	if (!(fp_stat = fopen(path, "w"))) {
		trace(1, "openstat: file open error path=%s\n", path);
		return 0;
	}
	strcpy(file_stat, file);
	time_stat = time;
	statlevel = level;
	return 1;
}
/* close solution status file --------------------------------------------------
* close solution status file
* args   : none
* return : none
*-----------------------------------------------------------------------------*/
extern void closestat(void)
{
	trace(3, "closestat:\n");

	if (fp_stat) fclose(fp_stat);
	fp_stat = NULL;
	file_stat[0] = '\0';
	statlevel = 0;
}


/* swap solution status file -------------------------------------------------*/
static void swapsolstat(void)
{
	gtime_t time = utc2gpst(timeget());
	char path[1024];

	if ((int)(time2gpst(time, NULL) / INT_SWAP_STAT) ==
		(int)(time2gpst(time_stat, NULL) / INT_SWAP_STAT)) {
		return;
	}
	time_stat = time;

	if (!reppath(file_stat, path, time, "", "")) {
		return;
	}
	if (fp_stat) fclose(fp_stat);

	if (!(fp_stat = fopen(path, "w"))) {
		trace(2, "swapsolstat: file open error path=%s\n", path);
		return;
	}
	trace(3, "swapsolstat: path=%s\n", path);
}
/* open solution status file ---------------------------------------------------
* open solution status file and set output level
* args   : char     *file   I   leo status file
*          int      level   I   leo status level (0: off)
* return : status (1:ok,0:error)
* notes  : file can constain time keywords (%Y,%y,%m...) defined in reppath().
*          The time to replace keywords is based on UTC of CPU time.
* output : solution status file record format
*
*   $POS,week,tow,stat,posx,posy,posz,posxf,posyf,poszf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          posx/posy/posz    : position x/y/z ecef (m) float
*          posxf/posyf/poszf : position x/y/z ecef (m) fixed
*
*   $VELACC,week,tow,stat,vele,veln,velu,acce,accn,accu,velef,velnf,veluf,accef,accnf,accuf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          vele/veln/velu    : velocity e/n/u (m/s) float
*          acce/accn/accu    : acceleration e/n/u (m/s^2) float
*          velef/velnf/veluf : velocity e/n/u (m/s) fixed
*          accef/accnf/accuf : acceleration e/n/u (m/s^2) fixed
*
*   $CLK,week,tow,stat,clk1,clk2,clk3,clk4
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          clk1     : receiver clock bias GPS (ns)
*          clk2     : receiver clock bias GLO-GPS (ns)
*          clk3     : receiver clock bias GAL-GPS (ns)
*          clk4     : receiver clock bias BDS-GPS (ns)
*
*   $ION,week,tow,stat,sat,az,el,ion,ion-fixed
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          sat      : satellite id
*          az/el    : azimuth/elevation angle(deg)
*          ion      : vertical ionospheric delay L1 (m) float
*          ion-fixed: vertical ionospheric delay L1 (m) fixed
*
*   $TROP,week,tow,stat,rcv,ztd,ztdf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          rcv      : receiver (1:rover,2:base station)
*          ztd      : zenith total delay (m) float
*          ztdf     : zenith total delay (m) fixed
*
*   $HWBIAS,week,tow,stat,frq,bias,biasf
*          week/tow : gps week no/time of week (s)
*          stat     : solution status
*          frq      : frequency (1:L1,2:L2,...)
*          bias     : h/w bias coefficient (m/MHz) float
*          biasf    : h/w bias coefficient (m/MHz) fixed
*
*   $SAT,week,tow,sat,frq,az,el,resp,resc,vsat,snr,fix,slip,lock,outc,slipc,rejc
*          week/tow : gps week no/time of week (s)
*          sat/frq  : satellite id/frequency (1:L1,2:L2,...)
*          az/el    : azimuth/elevation angle (deg)
*          resp     : pseudorange residual (m)
*          resc     : carrier-phase residual (m)
*          vsat     : valid data flag (0:invalid,1:valid)
*          snr      : signal strength (dbHz)
*          fix      : ambiguity flag  (0:no data,1:float,2:fixed,3:hold,4:ppp)
*          slip     : cycle-slip flag (bit1:slip,bit2:parity unknown)
*          lock     : carrier-lock count
*          outc     : data outage count
*          slipc    : cycle-slip count
*          rejc     : data reject (outlier) count
*
*-----------------------------------------------------------------------------*/
/* output solution status ----------------------------------------------------*/
extern void outsolstat(leo_t *leo)
{
	ssat_t *ssat;
	int  week;
	double tow;
	char buff[MAXSOLMSG + 1], id[32], *p = buff;
	double pos[3], vel[3], acc[3], *x;
	int i, j, k, n, nfreq, nf = NF(&leo->opt);
	double vara = 0.0;
	if (statlevel <= 0 || !fp_stat ) return;

	trace(3, "outsolstat:\n");

	/* swap solution status file */
	swapsolstat();

	/* write solution status */
	x = leo->sol.stat == SOLQ_FIX ? leo->xa : leo->x;
	tow = time2gpst(leo->sol.time, &week);
	/* receiver position */
	p += sprintf(p, "$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n", week, tow,
		leo->sol.stat, x[0], x[1], x[2], STD(leo, 0), STD(leo, 1), STD(leo, 2));

	/* receiver velocity and acceleration */
	if ((leo->opt.dynamics) && (norm(leo->sol.rp, 3) > 1)) {//Sheng Zhuang

		ecef2pos(leo->sol.rr, pos);
		ecef2enu(pos, leo->x + 3, vel);
		ecef2enu(pos, leo->x + 6, acc);
		p += sprintf(p, "$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,"
			"%.4f,%.5f,%.5f,%.5f\n", week, tow, leo->sol.stat, vel[0], vel[1],
			vel[2], acc[0], acc[1], acc[2], 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	}
	/* receiver clocks */
	i = IC(0, &leo->opt);
	p += sprintf(p, "$CLK,%d,%.3f,%d,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,\n",
		week, tow, leo->sol.stat, 1, x[i] * 1E9 / CLIGHT, x[i + 1] * 1E9 / CLIGHT, x[i + 2] * 1E9 / CLIGHT,
		x[i + 3] * 1E9 / CLIGHT,
		STD(leo, i)*1E9 / CLIGHT, STD(leo, i + 1)*1E9 / CLIGHT, STD(leo, i+2)*1E9 / CLIGHT, STD(leo, i + 3)*1E9 / CLIGHT);

	/* tropospheric parameters */
	if (leo->opt.tropopt == TROPOPT_EST || leo->opt.tropopt == TROPOPT_ESTG) {
		i = IT(&leo->opt);
		p += sprintf(p, "$TROP,%d,%.3f,%d,%d,%.4f,%.4f\n", week, tow, leo->sol.stat,
			1, x[i], STD(leo, i));
	}
	if (leo->opt.tropopt == TROPOPT_ESTG) {
		i = IT(&leo->opt);
		p += sprintf(p, "$TRPG,%d,%.3f,%d,%d,%.5f,%.5f,%.5f,%.5f\n", week, tow,
			leo->sol.stat, 1, x[i + 1], x[i + 2], STD(leo, i + 1), STD(leo, i + 2));
	}
	/* ionosphere parameters */
	if (leo->opt.ionoopt == IONOOPT_EST) {
		for (i = 0; i<MAXSAT; i++) {
			ssat = leo->ssat + i;
			if (!ssat->vs) continue;
			j = II(i + 1, &leo->opt);
			if (leo->x[j] == 0.0) continue;
			satno2id(i + 1, id);
			p += sprintf(p, "$ION,%d,%.3f,%d,%s,%.1f,%.1f,%.4f,%.4f\n", week, tow,
				leo->sol.stat, id, leo->ssat[i].azel[0] * R2D,
				leo->ssat[i].azel[1] * R2D, x[j], STD(leo, j));
		}
	}
	n = p - buff;
	buff[n] = '\0';

	fputs(buff, fp_stat);

	if (leo->sol.stat == SOLQ_NONE || statlevel <= 1) return;

	tow = time2gpst(leo->sol.time, &week);
	nfreq = leo->opt.mode >= PMODE_DGPS ? nf : 1;
	if (fabs(tow - 102660.0)<0.5) {
		i = 0;
	}
	/* write residuals and status */
	for (i = 0; i<MAXSAT; i++) {
		ssat = leo->ssat + i;
		if (!ssat->vs) continue;
		//if (!ssat->vsat[0]) continue;
		satno2id(i + 1, id);
		for (j = 0; j<nfreq; j++) {
			fprintf(fp_stat, "$SAT,%d,%.3f,%s,%d,%.1f,%.1f,%.4f,%.4f,%d,%.0f,%d,%d,%d,%d,%d,%d\n",
				week, tow, id, j + 1, ssat->azel[0] * R2D, ssat->azel[1] * R2D,
				ssat->resp[j], ssat->resc[j], ssat->vsat[j], ssat->snr[j] * 0.25,
				ssat->fix[j], ssat->slip[j] & 3, ssat->lock[j], ssat->outc[j],
				ssat->slipc[j], ssat->rejc[j]);
		}
	}
	/* ambiguity parameters */
	for (i = 0; i<MAXSAT; i++) for (j = 0; j<NF(&leo->opt); j++) {
		k = IB(i + 1, j, &leo->opt);
		if (leo->x[k] == 0.0) continue;
		ssat = leo->ssat + i;
		if (!ssat->vsat[0]) continue;
		satno2id(i + 1, id);
		vara = SQRT(leo->P[k + k*leo->nx]);
		if (leo->sol.stat == SOLQ_FIX) vara = SQRT(leo->Pa[k + k*leo->nx]);
		fprintf(fp_stat, "$AMB,%d,%.3f,%d,%s,%d,%.4f,%.4f\n", week, tow,
			leo->sol.stat, id, j + 1, leo->x[k], vara);
	}

}