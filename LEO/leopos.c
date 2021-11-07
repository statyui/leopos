/*------------------------------------------------------------------------------
* ppp.c : precise point positioning
*
*
* options : -DIERS_MODEL  use IERS tide model
*           -DOUTSTAT_AMB output ambiguity parameters to solution status
*
* references :
*    [1] D.D.McCarthy, IERS Technical Note 21, IERS Conventions 1996, July 1996
*    [2] D.D.McCarthy and G.Petit, IERS Technical Note 32, IERS Conventions
*        2003, November 2003
*    [3] D.A.Vallado, Fundamentals of Astrodynamics and Applications 2nd ed,
*        Space Technology Library, 2004
*    [4] J.Kouba, A Guide to using International GNSS Service (IGS) products,
*        May 2009
*    [5] RTCM Paper, April 12, 2010, Proposed SSR Messages for SV Orbit Clock,
*        Code Biases, URA
*    [6] MacMillan et al., Atmospheric gradients and the VLBI terrestrial and
*        celestial reference frames, Geophys. Res. Let., 1997
*    [7] G.Petit and B.Luzum (eds), IERS Technical Note No. 36, IERS
*         Conventions (2010), 2010
*    [8] J.Kouba, A simplified yaw-attitude model for eclipsing GPS satellites,
*        GPS Solutions, 13:1-12, 2009
*    [9] F.Dilssner, GPS IIF-1 satellite antenna phase center and attitude
*        modeling, InsideGNSS, September, 2010
*    [10] F.Dilssner, The GLONASS-M satellite yaw-attitude model, Advances in
*        Space Research, 2010
*    [11] IGS MGEX (http://igs.org/mgex)
*
* version : $Revision:$ $Date:$
* history : 2010/07/20 1.0  new
*                           added api:
*                               tidedisp()
*           2010/12/11 1.1  enable exclusion of eclipsing satellite
*           2012/02/01 1.2  add gps-glonass h/w bias correction
*                           move windupcorr() to leocmn.c
*           2013/03/11 1.3  add otl and pole tides corrections
*                           involve iers model with -DIERS_MODEL
*                           change initial variances
*                           suppress acos domain error
*           2013/09/01 1.4  pole tide model by iers 2010
*                           add mode of ionosphere model off
*           2014/05/23 1.5  add output of trop gradient in solution status
*           2014/10/13 1.6  fix bug on P0(a[3]) computation in tide_oload()
*                           fix bug on m2 computation in tide_pole()
*           2015/03/19 1.7  fix bug on ionosphere correction for GLO and BDS
*           2015/05/10 1.8  add function to detect slip by MW-LC jump
*                           fix ppp solutin problem with large clock variance
*           2015/06/08 1.9  add precise satellite yaw-models
*                           cope with day-boundary problem of satellite clock
*           2015/07/31 1.10 fix bug on nan-solution without glonass nav-data
*                           pppoutsolsat() -> pppoutstat()
*           2015/11/13 1.11 add L5-receiver-dcb estimation
*                           merge post-residual validation by rnx2leop_test
*                           support support option opt->pppopt=-GAP_RESION=nnnn
*           2016/01/22 1.12 delete support for yaw-model bug
*                           add support for ura of ephemeris
*-----------------------------------------------------------------------------*/
#include "leoorb.h"

#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0?0.0:sqrt(x))
#define MAX(x,y)    ((x)>(y)?(x):(y))
#define MIN(x,y)    ((x)<(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

#define MU_GPS   3.9860050E14     /* gravitational constant         ref [1] */
#define MU_GLO   3.9860044E14     /* gravitational constant         ref [2] */
#define MU_GAL   3.986004418E14   /* earth gravitational constant   ref [7] */
#define MU_CMP   3.986004418E14   /* earth gravitational constant   ref [9] */

#define MAX_ITER    8               /* max number of iterations */
#define MAX_STD_FIX 0.15            /* max std-dev (3d) to fix solution */
#define MIN_NSAT_SOL 4              /* min satellite number for solution */
#define THRES_REJECT 4.0            /* reject threshold of posfit-res (sigma) */

#define THRES_MW_JUMP 1

#define VAR_POS     SQR(30.0)       /* init variance receiver position (m^2) */
#define VAR_VEL     SQR(10.0)       /* init variance receiver position (m^2) */
#define VAR_CLK     SQR(60.0)       /* init variance receiver clock (m^2) */
#define VAR_ZTD     SQR( 0.6)       /* init variance ztd (m^2) */
#define VAR_GRA     SQR(0.01)       /* init variance gradient (m^2) */
#define VAR_DCB     SQR(30.0)       /* init variance dcb (m^2) */
#define VAR_BIAS    SQR(60.0)       /* init variance phase-bias (m^2) */
#define VAR_IONO    SQR(60.0)       /* init variance iono-delay */
#define VAR_GLO_IFB SQR( 0.6)       /* variance of glonass ifb */

#define ERR_SAAS    0.3             /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5             /* broadcast iono model error factor */
#define ERR_CBIAS   0.3             /* code bias error std (m) */
#define REL_HUMI    0.7             /* relative humidity for saastamoinen model */
#define GAP_RESION  120             /* default gap to reset ionos parameters (ep) */

#define EFACT_GPS_L5 10.0           /* error factor of GPS/QZS L5 */

#define MUDOT_GPS   (0.00836*D2R)   /* average angular velocity GPS (rad/s) */
#define MUDOT_GLO   (0.00888*D2R)   /* average angular velocity GLO (rad/s) */
#define EPS0_GPS    (13.5*D2R)      /* max shadow crossing angle GPS (rad) */
#define EPS0_GLO    (14.2*D2R)      /* max shadow crossing angle GLO (rad) */
#define T_POSTSHADOW 1800.0         /* post-shadow recovery time (s) */

/* number and index of states */
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
#define IC(s,opt)   (NP(opt)+(s))
#define IT(opt)     (NP(opt)+NC(opt))
#define II(s,opt)   (NP(opt)+NC(opt)+NT(opt)+(s)-1)
#define ID(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt))
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1)



static void errmsg(leo_t *leo, const char *format, ...)
{
	char buff[256], tstr[32];
	int n;
	va_list ap;
	time2str(leo->sol.time, tstr, 2);
	n = sprintf(buff, "%s: ", tstr + 11);
	va_start(ap, format);
	n += vsprintf(buff + n, format, ap);
	va_end(ap);
	n = n<MAXERRMSG - leo->neb ? n : MAXERRMSG - leo->neb;
	memcpy(leo->errbuf + leo->neb, buff, n);
	leo->neb += n;
	trace(2, "%s", buff);
}
/* exclude meas of eclipsing satellite (block IIA) ---------------------------*/
static void testeclipse(const obsd_t *obs, int n, const nav_t *nav, double *rs)
{
    double rsun[3],esun[3],r,ang,erpv[5]={0},cosa;
    int i,j;
    const char *type;
    
    trace(3,"testeclipse:\n");
    
    /* unit vector of sun direction (ecef) */
    sunmoonpos(gpst2utc(obs[0].time),erpv,rsun,NULL,NULL);
    normv3(rsun,esun);
    
    for (i=0;i<n;i++) {
        type=nav->pcvs[obs[i].sat-1].type;
        
        if ((r=norm(rs+i*6,3))<=0.0) continue;
        
        /* only block IIA */
        if (*type&&!strstr(type,"BLOCK IIA")) continue;
        
        /* sun-earth-satellite angle */
        cosa=dot(rs+i*6,esun,3)/r;
        cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
        ang=acos(cosa);
        
        /* test eclipse */
        if (ang<PI/2.0||r*sin(ang)>RE_WGS84) continue;
        
        trace(3,"eclipsing sat excluded %s sat=%2d\n",time_str(obs[0].time,0),
              obs[i].sat);
        
        for (j=0;j<3;j++) rs[j+i*6]=0.0;
    }
}
/* nominal yaw-angle ---------------------------------------------------------*/
static double yaw_nominal(double beta, double mu)
{
    if (fabs(beta)<1E-12&&fabs(mu)<1E-12) return PI;
    return atan2(-tan(beta),sin(mu))+PI;
}
/* yaw-angle of satellite ----------------------------------------------------*/
extern int yaw_angle(int sat, const char *type, int opt, double beta, double mu,
                     double *yaw)
{
    *yaw=yaw_nominal(beta,mu);
    return 1;
}
/* satellite attitude model --------------------------------------------------*/
static int sat_yaw(gtime_t time, int sat, const char *type, int opt,
                   const double *rs, double *exs, double *eys)
{
    double rsun[3],ri[6],es[3],esun[3],n[3],p[3],en[3],ep[3],ex[3],E,beta,mu;
    double yaw,cosy,siny,erpv[5]={0};
    int i;
    
    sunmoonpos(gpst2utc(time),erpv,rsun,NULL,NULL);
    
    /* beta and orbit angle */
    matcpy(ri,rs,6,1);
    ri[3]-=OMGE*ri[1];
    ri[4]+=OMGE*ri[0];
    cross3(ri,ri+3,n);
    cross3(rsun,n,p);
    if (!normv3(rs,es)||!normv3(rsun,esun)||!normv3(n,en)||
        !normv3(p,ep)) return 0;
    beta=PI/2.0-acos(dot(esun,en,3));
    E=acos(dot(es,ep,3));
    mu=PI/2.0+(dot(es,esun,3)<=0?-E:E);
    if      (mu<-PI/2.0) mu+=2.0*PI;
    else if (mu>=PI/2.0) mu-=2.0*PI;
    
    /* yaw-angle of satellite */
    if (!yaw_angle(sat,type,opt,beta,mu,&yaw)) return 0;
    
    /* satellite fixed x,y-vector */
    cross3(en,es,ex);
    cosy=cos(yaw);
    siny=sin(yaw);
    for (i=0;i<3;i++) {
        exs[i]=-siny*en[i]+cosy*ex[i];
        eys[i]=-cosy*en[i]-siny*ex[i];
    }
    return 1;
}
/* phase windup model --------------------------------------------------------*/  
static int model_phw(gtime_t time, int sat, const char *type, int opt,
                     const double *rs, const double *rr, double *phw)
{
	double exs[3],eys[3],ek[3],exr[3],eyr[3],eks[3],ekr[3],E[9];
	double dr[3],ds[3],drs[3],r[3],pos[3],cosp,ph;
	int i;

	if (opt<=0) return 1; /* no phase windup */

	if (norm(rr,3)<=0.0) return 0;

	/* satellite yaw attitude model */
	if (!sat_yaw(time,sat,type,opt,rs,exs,eys)) return 0;

	/* unit vector satellite to receiver */
	for (i=0;i<3;i++) r[i]=rr[i]-rs[i];
	if (!normv3(r,ek)) return 0;

	/* unit vectors of receiver antenna */
	ecef2pos(rr,pos);
	xyz2enu(pos,E);
	exr[0]= E[1]; exr[1]= E[4]; exr[2]= E[7]; /* x = north */
	eyr[0]=-E[0]; eyr[1]=-E[3]; eyr[2]=-E[6]; /* y = west  */

	/* phase windup effect */
	cross3(ek,eys,eks);
	cross3(ek,eyr,ekr);
	for (i=0;i<3;i++) {
		ds[i]=exs[i]-ek[i]*dot(ek,exs,3)-eks[i];
		dr[i]=exr[i]-ek[i]*dot(ek,exr,3)+ekr[i];
	}
	cosp=dot(ds,dr,3)/norm(ds,3)/norm(dr,3);
	if      (cosp<-1.0) cosp=-1.0;
	else if (cosp> 1.0) cosp= 1.0;
	if (fabs(fabs(cosp)-1.0)<1.0e-10) return 0;

	ph=acos(cosp)/2.0/PI;
	cross3(ds,dr,drs);
	if (dot(ek,drs,3)<0.0) ph=-ph;

	*phw=ph+floor(*phw-ph+0.5); /* in cycle */
	return 1;
}
/* measurement error variance ------------------------------------------------*/
static double varerr(int sat, int sys, double el, int freq, int type,
                     const prcopt_t *opt)
{
    double fact=1.0,sinel=sin(el);
    
    if (type==1) fact*=opt->eratio[freq==0?0:1];
    fact*=sys==SYS_GLO?EFACT_GLO:EFACT_GPS;
    
    if (sys==SYS_GPS) {
        if (freq==2) fact*=EFACT_GPS_L5; /* GPS/QZS L5 error factor */
    }
    if (opt->ionoopt==IONOOPT_IFLC) fact*=3.0;
    return SQR(fact*opt->err[1])+SQR(fact*opt->err[2]/sinel);
}
/* initialize state and covariance -------------------------------------------*/
static void initx(leo_t *leo, double xi, double var, int i)
{
    int j;
    leo->x[i]=xi;
    for (j=0;j<leo->nx;j++) {
        leo->P[i+j*leo->nx]=leo->P[j+i*leo->nx]=i==j?var:0.0;
    }
}
/* geometry-free phase measurement -------------------------------------------*/
static double gfmeas(const obsd_t *obs, const nav_t *nav)
{
    const double *lam=nav->lam[obs->sat-1];
    int i=(satsys(obs->sat,NULL)&(SYS_GAL|SYS_SBS))?2:1;
    
    if (lam[0]==0.0||lam[i]==0.0||obs->L[0]==0.0||obs->L[i]==0.0) return 0.0;
    return lam[0]*obs->L[0]-lam[i]*obs->L[i];
}
/* Melbourne-Wubbena linear combination --------------------------------------*/
static double mwmeas(const obsd_t *obs, const nav_t *nav)
{
    const double *lam=nav->lam[obs->sat-1];
    int i=(satsys(obs->sat,NULL)&(SYS_GAL|SYS_SBS))?2:1;
    
    if (lam[0]==0.0||lam[i]==0.0||obs->L[0]==0.0||obs->L[i]==0.0||
        obs->P[0]==0.0||obs->P[i]==0.0) return 0.0;
    return lam[0]*lam[i]*(obs->L[0]-obs->L[i])/(lam[i]-lam[0])-
           (lam[i]*obs->P[0]+lam[0]*obs->P[i])/(lam[i]+lam[0]);
}
/* antenna corrected measurements --------------------------------------------*/ 
static void corr_meas(const obsd_t *obs, const nav_t *nav, const double *azel,
                      const prcopt_t *opt, const double *dantr,
                      const double *dants, double phw, double *L, double *P,
                      double *Lc, double *Pc)
{
    const double *lam=nav->lam[obs->sat-1];
    double C1,C2;
    int i,sys;
    
    for (i=0;i<NFREQ;i++) {
        L[i]=P[i]=0.0;
        if (lam[i]==0.0||obs->L[i]==0.0||obs->P[i]==0.0) continue;
        if (testsnr(0,0,azel[1],obs->SNR[i]*0.25,&opt->snrmask)) continue; 
        
        /* antenna phase center and phase windup correction */
        L[i]=obs->L[i]*lam[i]-dants[i]-dantr[i]-phw*lam[i];
        P[i]=obs->P[i]       -dants[i]-dantr[i];
        
        /* P1-C1,P2-C2 dcb correction (C1->P1,C2->P2) */ 
        if (obs->code[i]==CODE_L1C) {
            P[i]+=nav->cbias[obs->sat-1][1];
        }
        else if (obs->code[i]==CODE_L2C||obs->code[i]==CODE_L2X||
                 obs->code[i]==CODE_L2L||obs->code[i]==CODE_L2S) {
            P[i]+=nav->cbias[obs->sat-1][2];
#if 0
            L[i]-=0.25*lam[i]; /* 1/4 cycle-shift */
#endif
        }
    }
    /* iono-free LC */
    *Lc=*Pc=0.0;
    sys=satsys(obs->sat,NULL);
    i=(sys&(SYS_GAL|SYS_SBS))?2:1; /* L1/L2 or L1/L5 */
    if (lam[0]==0.0||lam[i]==0.0) return;
    
    C1= SQR(lam[i])/(SQR(lam[i])-SQR(lam[0]));
    C2=-SQR(lam[0])/(SQR(lam[i])-SQR(lam[0]));
    
#if 0
    /* P1-P2 dcb correction (P1->Pc,P2->Pc) */
    if (sys&(SYS_GPS|SYS_GLO|SYS_QZS)) {
        if (P[0]!=0.0) P[0]-=C2*nav->cbias[obs->sat-1][0];
        if (P[1]!=0.0) P[1]+=C1*nav->cbias[obs->sat-1][0];
    }
#endif
    if (L[0]!=0.0&&L[i]!=0.0) *Lc=C1*L[0]+C2*L[i];
    if (P[0]!=0.0&&P[i]!=0.0) *Pc=C1*P[0]+C2*P[i];
}
/* detect cycle slip by LLI --------------------------------------------------*/
static void detslp_ll(leo_t *leo, const obsd_t *obs, int n)
{
    int i,j;
    
    trace(3,"detslp_ll: n=%d\n",n);
    
    for (i=0;i<n&&i<MAXOBS;i++) for (j=0;j<leo->opt.nf;j++) {
        if (obs[i].L[j]==0.0||!(obs[i].LLI[j]&3)) continue;
        
        trace(3,"detslp_ll: slip detected sat=%2d f=%d\n",obs[i].sat,j+1);
        
        leo->ssat[obs[i].sat-1].slip[j]=1;
    }
}
/* detect cycle slip by geometry free phase jump -----------------------------*/
static void detslp_gf(leo_t *leo, const obsd_t *obs, int n, const nav_t *nav)
{
    double g0,g1;
    int i,j,id;
    
    trace(3,"detslp_gf: n=%d\n",n);
    
    for (i=0;i<n&&i<MAXOBS;i++) {
        
        if ((g1=gfmeas(obs+i,nav))==0.0) continue;
        
        g0=leo->ssat[obs[i].sat-1].gf;/* geometry-free phase L1-L2 (m) */
        leo->ssat[obs[i].sat-1].gf=g1;/* geometry-free phase L1-L2 (m) */
		if (g0 == 0 || g1 == 0) id = 1;
		else id = 0;
        trace(6,"detslip_gf: sat=%2d gf0=%8.3f gf1=%8.3f %d dgf=%8.3f\n",obs[i].sat,g0,g1,id,g1-g0);
        
        if (g0!=0.0&&fabs(g1-g0)>leo->opt.thresslip) {//slip threshold of geometry-free phase (m)
            trace(3,"detslip_gf: slip detected sat=%2d gf=%8.3f-%8.3f=%8.3f\n",
                  obs[i].sat,g0,g1,g0-g1);
            
            for (j=0;j<leo->opt.nf;j++) leo->ssat[obs[i].sat-1].slip[j]|=1;
        }
    }
}
/* detect slip by Melbourne-Wubbena linear combination jump ------------------*/
static void detslp_mw(leo_t *leo, const obsd_t *obs, int n, const nav_t *nav)
{
    double w0,w1;
    int i,j,id;

    
    trace(3,"detslp_mw: n=%d\n",n);
    
    for (i=0;i<n&&i<MAXOBS;i++) {
        if ((w1=mwmeas(obs+i,nav))==0.0) continue;
        
        w0=leo->ssat[obs[i].sat-1].mw;
        leo->ssat[obs[i].sat-1].mw=w1;
		if (w0 == 0 || w1 == 0) id = 1;
		else id = 0;
        trace(6,"detslip_mw: sat=%2d mw0=%8.3f mw1=%8.3f %d dmw=%8.3f\n",obs[i].sat,w0,w1,id,w1-w0);
        
        if (w0!=0.0&&fabs(w1-w0)>THRES_MW_JUMP) {
            trace(3,"detslip_mw: slip detected sat=%2d mw=%8.3f-%8.3f= %8.3f \n",
                  obs[i].sat,w0,w1, w0-w1);
            
            for (j=0;j<leo->opt.nf;j++) leo->ssat[obs[i].sat-1].slip[j]|=1;
        }
    }
}
/* temporal update of position -----------------------------------------------*/
static void udpos_leokin(leo_t *leo, gtime_t time, const nav_t *nav)
{
	int i, j;
	double y[6] = { 0.0 }, y1[6] = { 0.0 };
	double dy[36] = { 0.0 };
	double Q[36] = { 0.0 }, Qt[36] = { 0.0 };
	double dant[3] = { 0.0 };
	double dt = 0.0;
	double erpv[5] = { 0.0 };
	double kep[6] = { 0.0 }, pv[6] = { 0.0 };


	/* initialize position for first epoch */ 
	if (norm(leo->x, 3) <= 0.0) {
		for (i = 0; i<3; i++) initx(leo, leo->sol.rr[i], VAR_POS, i);
	}
	if (norm(leo->sol.rp, 3) > 1) {
		p22kepler(GM_EARTH, leo->tt, leo->sol.rp, leo->sol.rr, kep);
		trace(3, "$KEP,%13.4f \t");
		tracemat(3, kep, 1, 6, 18, 10);
		kepler2pv(GM_EARTH, kep, leo->tt, pv);
	}
	/* kinmatic mode without dynamics */
	//leoantoff(leo->sol.time, pv, nav, dant);
	memcpy(y, leo->sol.rr, 6 * sizeof(double));
	trace(3, "udpos_ppp: y=\t");
	tracemat(3, y, 1, 6, 13, 4);
	if (!leo->opt.dynamics) {
		for (i = 0; i<3; i++) {
			initx(leo, leo->sol.rr[i], VAR_POS, i);
		}
	}
	else {
		if (norm(leo->x, 3) < 1e-8) {//init 
			for (i = 0; i < 3; i++) {
				initx(leo, leo->sol.rr[i], VAR_POS, i);
				initx(leo, leo->sol.rr[i + 3], VAR_VEL, i + 3);
			}
		}
	}

}
/* temporal update of position -----------------------------------------------*/
static void udpos_leordn(leo_t *leo, gtime_t time, const nav_t *nav)
{
	int i, j;
	double y[6] = { 0.0 }, y1[6] = { 0.0 };
	double dy[36] = { 0.0 };
	double Q[36] = { 0.0 }, Qt[36] = { 0.0 };
	double dant[3] = { 0.0 };
	double dt = 0.0;
	double erpv1[5] = { 0.0 }, erpv2[5] = { 0.0 };
	double kep[6] = { 0.0 }, pv[6] = { 0.0 };
	double  U1[9] = { 0.0 }, U2[9] = { 0.0 }, dU1[9] = { 0.0 }, dU2[9] = { 0.0 };
	double rr[3] = { 0.0 };  // Body-fixed position
	double rp[3] = { 0.0 };  // Body-fixed acceleration     
	gtime_t now;

	gtime_t tutc;
	
	now = gpst2utc(time);
	tutc = timeadd(time, -(leo->tt));
	tutc = gpst2utc(tutc);
	if (&nav->erp) geterp(&nav->erp, tutc, erpv1);
	if (&nav->erp) geterp(&nav->erp, now, erpv2);
	/* initialize position for first epoch */
	if (norm(leo->x, 3) <= 0.0) {
		for (i = 0; i<3; i++) initx(leo, leo->sol.rr[i], VAR_POS, i);
	}
	/* static ppp mode */
	if (leo->opt.mode == PMODE_PPP_STATIC) {
		for (i = 0; i<3; i++) {
			leo->P[i*(1 + leo->nx)] += SQR(leo->opt.prn[5])*fabs(leo->tt);
		}
		return;
	}
	if (norm(leo->sol.rp, 3) > 1) {


		//ecef2qeci(tutc, erpv1[2], U1, dU1);
		//matmul("TN", 3, 1, 3, 1.0, U1, leo->sol.rp, 0.0, rp);





		//Differential calculation
		ecef2qeci(now, erpv2[2], U2, dU2);
		ecef2qeci(tutc, erpv1[2], U1, dU1);
		matmul("TN", 3, 1, 3, 1.0, U1, leo->sol.rp, 0.0, rp);
		matmul("TN", 3, 1, 3, 1.0, U2, leo->sol.rr, 0.0, rr);


		//p22kepler(GM_EARTH, leo->tt, rp, rr, kep);//at time rp
		GetVo(GM_EARTH, leo->tt, rp, rr, pv);

		trace(3, "$KEP,%13.4f \t");
		tracemat(3, kep, 1, 6, 18, 10);
		//kepler2pv(GM_EARTH, kep, leo->tt, pv);//at time rp

		/* kinmatic mode without dynamics */
		//leoantoff(leo->sol.time, pv, nav, dant);

		memcpy(y, pv, 6 * sizeof(double));//Sheng Zhuang
	}

	trace(3, "udpos_ppp: y=\t");
	tracemat(3, y, 1, 6, 13, 4);
	if ((!leo->opt.dynamics) || (y[4] == 0)) {//Sheng Zhuang
		for (i = 0; i<3; i++) {
			initx(leo, leo->sol.rr[i], VAR_POS, i);
		}
	}
	else {
		if (norm(leo->x, 3) < 1e-8) {//init 
			for (i = 0; i < 3; i++) {
				initx(leo, leo->sol.rr[i], VAR_POS, i);
				initx(leo, leo->sol.rr[i + 3], VAR_VEL, i + 3);
			}
		}
		else {//propogate

			pvpartial(y, leo->tt, y1, dy);

			for (i = 0; i < 6; i++) {
				for (j = 0; j < 6; j++) {
					Qt[i + j * 6] = leo->P[i + j*leo->nx];
				}
			}
			matmul("NN", 6, 6, 6, 1.0, dy, Qt, 0.0, Qt);
			matmul("NT", 6, 6, 6, 1.0, Qt, dy, 0.0, Q);
			rk4(leo->tt, tutc, nav, 10, y,erpv2);
			
			memcpy(leo->x, y, 6 * sizeof(double));
			for (i = 0; i < 6; i++) {
				leo->x[i] = y[i];
				for (j = 0; j < 6; j++) {
					leo->P[i + j*leo->nx] = Q[i + j * 6] + 1.0;
				}
			}
		}
	}

}
/* temporal update of position -----------------------------------------------*/
static void udpos_ppp(leo_t *leo, gtime_t time, const nav_t *nav)
{
	int i;

	trace(3, "udpos_ppp:\n");

	/* fixed mode */ 
	if (leo->opt.mode == PMODE_PPP_FIXED) {
		for (i = 0; i<3; i++) initx(leo, leo->opt.ru[i], 1E-8, i);
		return;
	}
	/* initialize position for first epoch */
	if (norm(leo->x, 3) <= 0.0) {
		for (i = 0; i<3; i++) initx(leo, leo->sol.rr[i], VAR_POS, i);
	}
	/* static ppp mode */
	if (leo->opt.mode == PMODE_PPP_STATIC) {
		for (i = 0; i<3; i++) {
			leo->P[i*(1 + leo->nx)] += SQR(leo->opt.prn[5])*fabs(leo->tt);
		}
		return;
	}
	/* kinmatic mode without dynamics */
	for (i = 0; i<3; i++) {
		initx(leo, leo->sol.rr[i], VAR_POS, i);
	}
}
/* temporal update of clock --------------------------------------------------*/
static void udclk_ppp(leo_t *leo)
{
    double dtr;
    int i;
    
    trace(3,"udclk_ppp:\n");
    
    /* initialize every epoch for clock (white noise) */
    for (i=0;i<NSYS;i++) {
        //if (leo->opt.sateph==EPHOPT_PREC) {
        //    /* time of prec ephemeris is based gpst */
        //    /* negelect receiver inter-system bias  */
        //    dtr=leo->sol.dtr[0];
        //}
        //else {
            dtr=i==0?leo->sol.dtr[0]:leo->sol.dtr[0]+leo->sol.dtr[i];//leiwang
        //}

        initx(leo,CLIGHT*dtr,VAR_CLK,IC(i,&leo->opt));

    }
}
/* get meterological parameters ----------------------------------------------*/
static void getmet(double lat, double *met)
{
	static const double metprm[][10] = { /* lat=15,30,45,60,75 */
		{ 1013.25,299.65,26.31,6.30E-3,2.77,  0.00, 0.00,0.00,0.00E-3,0.00 },
		{ 1017.25,294.15,21.79,6.05E-3,3.15, -3.75, 7.00,8.85,0.25E-3,0.33 },
		{ 1015.75,283.15,11.66,5.58E-3,2.57, -2.25,11.00,7.24,0.32E-3,0.46 },
		{ 1011.75,272.15, 6.78,5.39E-3,1.81, -1.75,15.00,5.36,0.81E-3,0.74 },
		{ 1013.00,263.65, 4.11,4.53E-3,1.55, -0.50,14.50,3.39,0.62E-3,0.30 }
	};
	int i, j;
	double a;
	lat = fabs(lat);
	if (lat <= 15.0) for (i = 0; i<10; i++) met[i] = metprm[0][i];
	else if (lat >= 75.0) for (i = 0; i<10; i++) met[i] = metprm[4][i];
	else {
		j = (int)(lat / 15.0); a = (lat - j*15.0) / 15.0;
		for (i = 0; i<10; i++) met[i] = (1.0 - a)*metprm[j - 1][i] + a*metprm[j][i];
	}
}
/* tropospheric delay correction -----------------------------------------------
* compute sbas tropospheric delay correction (mops model)
* args   : gtime_t time     I   time
*          double   *pos    I   receiver position {lat,lon,height} (rad/m)
*          double   *azel   I   satellite azimuth/elavation (rad)
*          double   *var    O   variance of troposphric error (m^2)
* return : slant tropospheric delay (m) 
*-----------------------------------------------------------------------------*/
extern double sbstropcorr(gtime_t time, const double *pos, const double *azel,
	double *var)
{
	const double k1 = 77.604, k2 = 382000.0, rd = 287.054, gm = 9.784, g = 9.80665;
	static double pos_[3] = { 0 }, zh = 0.0, zw = 0.0;
	int i;
	double c, met[10], sinel = sin(azel[1]), h = pos[2], m;

	trace(4, "sbstropcorr: pos=%.3f %.3f azel=%.3f %.3f\n", pos[0] * R2D, pos[1] * R2D,
		azel[0] * R2D, azel[1] * R2D);

	if (pos[2]<-100.0 || 10000.0<pos[2] || azel[1] <= 0) {
		*var = 0.0;
		return 0.0;
	}
	if (zh == 0.0 || fabs(pos[0] - pos_[0])>1E-7 || fabs(pos[1] - pos_[1])>1E-7 ||
		fabs(pos[2] - pos_[2])>1.0) {
		getmet(pos[0] * R2D, met);
		c = cos(2.0*PI*(time2doy(time) - (pos[0] >= 0.0 ? 28.0 : 211.0)) / 365.25);
		for (i = 0; i<5; i++) met[i] -= met[i + 5] * c;
		zh = 1E-6*k1*rd*met[0] / gm;
		zw = 1E-6*k2*rd / (gm*(met[4] + 1.0) - met[3] * rd)*met[2] / met[1];
		zh *= pow(1.0 - met[3] * h / met[1], g / (rd*met[3]));
		zw *= pow(1.0 - met[3] * h / met[1], (met[4] + 1.0)*g / (rd*met[3]) - 1.0);
		for (i = 0; i<3; i++) pos_[i] = pos[i];
	}
	m = 1.001 / sqrt(0.002001 + sinel*sinel);
	*var = 0.12*0.12*m*m;
	return (zh + zw)*m;
}

/* temporal update of tropospheric parameters --------------------------------*/
static void udtrop_ppp(leo_t *leo)
{
    double pos[3],azel[]={0.0,PI/2.0},ztd,var;
    int i=IT(&leo->opt),j;
    
    trace(3,"udtrop_ppp:\n");
    
    if (leo->x[i]==0.0) {
        ecef2pos(leo->sol.rr,pos);
        ztd=sbstropcorr(leo->sol.time,pos,azel,&var);
        initx(leo,ztd,var,i);
        
        if (leo->opt.tropopt>=TROPOPT_ESTG) {
            for (j=i+1;j<i+3;j++) initx(leo,1E-6,VAR_GRA,j);
        }
    }
    else {
        leo->P[i+i*leo->nx]+=SQR(leo->opt.prn[2])*fabs(leo->tt);
        
        if (leo->opt.tropopt>=TROPOPT_ESTG) {
            for (j=i+1;j<i+3;j++) {
                leo->P[j+j*leo->nx]+=SQR(leo->opt.prn[2]*0.1)*fabs(leo->tt);
            }
        }
    }
}
/* temporal update of ionospheric parameters ---------------------------------*/
static void udiono_ppp(leo_t *leo, const obsd_t *obs, int n, const nav_t *nav)
{
    const double *lam;
    double ion,sinel,pos[3],*azel;
    char *p;
    int i,j,k,gap_resion=GAP_RESION;
    
    trace(3,"udiono_ppp:\n");
    
    if ((p=strstr(leo->opt.pppopt,"-GAP_RESION="))) {
        sscanf(p,"-GAP_RESION=%d",&gap_resion);
        
    }
    for (i=0;i<MAXSAT;i++) {
        j=II(i+1,&leo->opt);
        if (leo->x[j]!=0.0&&(int)leo->ssat[i].outc[0]>gap_resion) {
            leo->x[j]=0.0;
            
           
        }
    }
    for (i=0;i<n;i++) {
        j=II(obs[i].sat,&leo->opt);
        if (leo->x[j]==0.0) {
            k=satsys(obs[i].sat,NULL)==SYS_GAL?2:1;
            lam=nav->lam[obs[i].sat-1];
            if (obs[i].P[0]==0.0||obs[i].P[k]==0.0||lam[0]==0.0||lam[k]==0.0) {
                continue;
            }
            ion=(obs[i].P[0]-obs[i].P[k])/(1.0-SQR(lam[k]/lam[0]));
            ecef2pos(leo->sol.rr,pos);
            azel=leo->ssat[obs[i].sat-1].azel;
            ion/=ionmapf(pos,azel);
            initx(leo,ion,VAR_IONO,j);
        }
        else {
            sinel=sin(MAX(leo->ssat[obs[i].sat-1].azel[1],5.0*D2R));
            leo->P[j+j*leo->nx]+=SQR(leo->opt.prn[1]/sinel)*fabs(leo->tt);
        }
    }
}
/* temporal update of L5-receiver-dcb parameters -----------------------------*/
static void uddcb_ppp(leo_t *leo)
{
    int i=ID(&leo->opt);
    
    trace(3,"uddcb_ppp:\n");
    
    if (leo->x[i]==0.0) {
        initx(leo,1E-6,VAR_DCB,i);
    }
}
/* temporal update of phase biases -------------------------------------------*/
static void udbias_ppp(leo_t *leo, const obsd_t *obs, int n, const nav_t *nav)
{
    const double *lam;
    double L[NFREQ],P[NFREQ],Lc,Pc,bias[MAXOBS],offset=0.0,pos[3]={0};
    double ion,dantr[NFREQ]={0},dants[NFREQ]={0};
    int i,j,k,l,f,sat,slip[MAXOBS]={0},clk_jump=0;
    
    trace(3,"udbias  : n=%d\n",n);
    
    /* handle day-boundary clock jump */ 
    if (leo->opt.posopt[5]) {
        clk_jump=ROUND(time2gpst(obs[0].time,NULL)*10)%864000==0;                  
    }
    for (i=0;i<MAXSAT;i++) for (j=0;j<leo->opt.nf;j++) {
        leo->ssat[i].slip[j]=0;
    }
    /* detect cycle slip by LLI */
    detslp_ll(leo,obs,n);
    
    /* detect cycle slip by geometry-free phase jump */ 
    detslp_gf(leo,obs,n,nav);
    
    /* detect slip by Melbourne-Wubbena linear combination jump */ 
    detslp_mw(leo,obs,n,nav);
    
    ecef2pos(leo->sol.rr,pos);
    
    for (f=0;f<NF(&leo->opt);f++) {
        
        /* reset phase-bias if expire obs outage counter */ 
        for (i=0;i<MAXSAT;i++) {
            if (++leo->ssat[i].outc[f]>(unsigned int)leo->opt.maxout||
                leo->opt.modear==ARMODE_INST||clk_jump) {
                initx(leo,0.0,0.0,IB(i+1,f,&leo->opt));
            }
        }
        for (i=k=0;i<n&&i<MAXOBS;i++) {
            sat=obs[i].sat;
            j=IB(sat,f,&leo->opt);
            corr_meas(obs+i,nav,leo->ssat[sat-1].azel,&leo->opt,dantr,dants,
                      0.0,L,P,&Lc,&Pc);
            
            bias[i]=0.0;
            
            if (leo->opt.ionoopt==IONOOPT_IFLC) {
                bias[i]=Lc-Pc;
                slip[i]=leo->ssat[sat-1].slip[0]||leo->ssat[sat-1].slip[1];
            }
            else if (L[f]!=0.0&&P[f]!=0.0) {
                slip[i]=leo->ssat[sat-1].slip[f];
                l=satsys(sat,NULL)==SYS_GAL?2:1;
                lam=nav->lam[sat-1];
                if (obs[i].P[0]==0.0||obs[i].P[l]==0.0||
                    lam[0]==0.0||lam[l]==0.0||lam[f]==0.0) continue;
                ion=(obs[i].P[0]-obs[i].P[l])/(1.0-SQR(lam[l]/lam[0]));
                bias[i]=L[f]-P[f]+2.0*ion*SQR(lam[f]/lam[0]);
            }
            if (leo->x[j]==0.0||slip[i]||bias[i]==0.0) continue;
            
            offset+=bias[i]-leo->x[j];
            k++;
        }
        /* correct phase-code jump to ensure phase-code coherency */  
        if (k>=2&&fabs(offset/k)>0.0005*CLIGHT) {
            for (i=0;i<MAXSAT;i++) {
                j=IB(i+1,f,&leo->opt);
                if (leo->x[j]!=0.0) leo->x[j]+=offset/k;
            }
            trace(2,"phase-code jump corrected: %s n=%2d dt=%12.9fs\n",
                  time_str(leo->sol.time,0),k,offset/k/CLIGHT);
        }
        for (i=0;i<n&&i<MAXOBS;i++) {
            sat=obs[i].sat;
            j=IB(sat,f,&leo->opt);
            
            leo->P[j+j*leo->nx]+=SQR(leo->opt.prn[0])*fabs(leo->tt);
            
            if (bias[i]==0.0||(leo->x[j]!=0.0&&!slip[i])) continue;
            
            /* reinitialize phase-bias if detecting cycle slip */ 
            initx(leo,bias[i],VAR_BIAS,IB(sat,f,&leo->opt));
            
            /* reset fix flags */
            for (k=0;k<MAXSAT;k++) leo->ambc[sat-1].flags[k]=0;
            
            trace(5,"udbias_ppp: sat=%2d bias=%.3f\n",sat,bias[i]);
        }
    }
}

/* temporal update of states --------------------------------------------------*/
static void udstate_ppp(leo_t *leo, const obsd_t *obs, int n, const nav_t *nav)
{
	trace(3, "udstate_ppp: n=%d\n", n);
	gtime_t time = obs->time;

	/* temporal update of position */                          
	if (leo->opt.mode == PMODE_LEO_KIN) {
		udpos_leokin(leo, time, nav);
		//udpos_ppp(leo, time, nav);
	}
	else if (leo->opt.mode == PMODE_LEO_RDN) { 
		udpos_leordn(leo, time, nav);
	}
	else{
		udpos_ppp(leo, time, nav);
	}
    
    /* temporal update of clock */
    udclk_ppp(leo);
    
    /* temporal update of tropospheric parameters */ 
    if (leo->opt.tropopt==TROPOPT_EST||leo->opt.tropopt==TROPOPT_ESTG) {
        udtrop_ppp(leo);
    }
    /* temporal update of ionospheric parameters */                      
    if (leo->opt.ionoopt==IONOOPT_EST) {
        udiono_ppp(leo,obs,n,nav);
    }
    /* temporal update of L5-receiver-dcb parameters */
    if (leo->opt.nf>=3) {
        uddcb_ppp(leo);
    }
    /* temporal update of phase-bias */  
    udbias_ppp(leo,obs,n,nav);
}
/* satellite antenna phase center variation ----------------------------------*/
static void satantpcv(const double *rs, const double *rr, const pcv_t *pcv,
                      double *dant,double *nadir)
{
    double ru[3],rz[3],eu[3],ez[3],cosa;
    int i;
    
    for (i=0;i<3;i++) {
        ru[i]=rr[i]-rs[i];
        rz[i]=-rs[i];
    }
    if (!normv3(ru,eu)||!normv3(rz,ez)) return;


    cosa=dot(eu,ez,3);
    cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
    *nadir=acos(cosa);
    
    antmodel_s(pcv,*nadir,dant);
}
/* precise tropospheric model ------------------------------------------------*/
static double trop_model_prec(gtime_t time, const double *pos,
                              const double *azel, const double *x, double *dtdx,
                              double *var)
{
    const double zazel[]={0.0,PI/2.0};
    double zhd,m_h,m_w,cotz,grad_n,grad_e;
    
    /* zenith hydrostatic delay */
    zhd=tropmodel(time,pos,zazel,0.0);
    
    /* mapping function */
    m_h=tropmapf(time,pos,azel,&m_w);
    
    if (azel[1]>0.0) {
        
        /* m_w=m_0+m_0*cot(el)*(Gn*cos(az)+Ge*sin(az)): ref [6] */
        cotz=1.0/tan(azel[1]);
        grad_n=m_w*cotz*cos(azel[0]);
        grad_e=m_w*cotz*sin(azel[0]);
        m_w+=grad_n*x[1]+grad_e*x[2];
        dtdx[1]=grad_n*(x[0]-zhd);
        dtdx[2]=grad_e*(x[0]-zhd);
    }
    dtdx[0]=m_w;
    *var=SQR(0.01);
    return m_h*zhd+m_w*(x[0]-zhd);
}
/* tropospheric model ---------------------------------------------------------*/ 
static int model_trop(gtime_t time, const double *pos, const double *azel,
                      const prcopt_t *opt, const double *x, double *dtdx,
                      const nav_t *nav, double *dtrp, double *var)
{
    double trp[3]={0};
    
    if (opt->tropopt==TROPOPT_SAAS) {
        *dtrp=tropmodel(time,pos,azel,REL_HUMI);
        *var=SQR(ERR_SAAS);
        return 1;
    }
    if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
        matcpy(trp,x+IT(opt),opt->tropopt==TROPOPT_EST?1:3,1);
        *dtrp=trop_model_prec(time,pos,azel,trp,dtdx,var);
        return 1;
    }
    return 0;
}
extern double gravitycor(const int sys, const double *rr,const double *rs)
{
	double	rm;
	double	sm;
	double	d;

	rm = sqrt(rr[0] * rr[0] + rr[1] * rr[1] +
		rr[2] * rr[2]);
	sm = sqrt(rs[0] * rs[0] + rs[1] * rs[1] +
		rs[2] * rs[2]);
	d = sqrt((rs[0] - rr[0])*(rs[0] - rr[0]) +
		(rs[1] - rr[1])*(rs[1] - rr[1]) +
		(rs[2] - rr[2])*(rs[2] - rr[2]));


	return 2.0*MU_GPS / (CLIGHT*CLIGHT)*log((sm + rm + d) / (sm + rm - d));
}

/* ionospheric model ---------------------------------------------------------*/
static int model_iono(gtime_t time, const double *pos, const double *azel,
                      const prcopt_t *opt, int sat, const double *x,
                      const nav_t *nav, double *dion, double *var)
{
    static double iono_p[MAXSAT]={0},std_p[MAXSAT]={0};
    static gtime_t time_p;
    
    if (opt->ionoopt==IONOOPT_BRDC) {
        *dion=ionmodel(time,nav->ion_gps,pos,azel);
        *var=SQR(*dion*ERR_BRDCI);
        return 1;
    }
    if (opt->ionoopt==IONOOPT_EST) {
        *dion=x[II(sat,opt)];
        *var=0.0;
        return 1;
    }
    if (opt->ionoopt==IONOOPT_IFLC) {
        *dion=*var=0.0;
        return 1;
    }
    return 0;
}

/* phase and code residuals --------------------------------------------------*/ 
static int ppp_res(int post, const obsd_t *obs, int n, const double *rs,
                   const double *dts, const double *var_rs, const int *svh,
                   const double *dr, int *exc, const nav_t *nav,
                   const double *x, leo_t *leo, double *v, double *H, double *R,
                   double *azel)
{
    const double *lam;
    prcopt_t *opt=&leo->opt;
    double y,r,cdtr,bias,C,rr[3],pos[3],e[3],dtdx[3],L[NFREQ],P[NFREQ],Lc,Pc;
    double var[MAXOBS*2],dtrp=0.0,dion=0.0,vart=0.0,vari=0.0,dcb;
    double dantr[NFREQ]={0},dants[NFREQ]={0};
    double ve[MAXOBS*2*NFREQ]={0},vmax=0;
    char str[32];
	double nadir,dg;
	int ne = 0, obsi[MAXOBS * 2 * NFREQ] = { 0 }, frqi[MAXOBS * 2 * NFREQ], maxobs, maxfrq, rej;
	//int maxobs, maxfrq, rej;
    int i,j,k,sat,sys,nv=0,nx=leo->nx,stat=1;
    
    time2str(obs[0].time,str,2);
	trace(6, "ppp_res: %3d %s\n",post, str);
	for (i = 0; i < MAXSAT; i++){
		leo->ssat[i].tropmap = 0;
		for (j = 0; j < 3; j++) {
			leo->ssat[i].pcos[j] = 0;
		}
		for (j = 0; j < opt->nf; j++) {
		leo->ssat[i].vsat[j] = 0;
		leo->ssat[i].pcvs[j] = 0;
		leo->ssat[i].pcvr[j] = 0;

		}
    }

    for (i=0;i<3;i++) rr[i]=x[i]+dr[i];
    ecef2pos(rr,pos);
    
    for (i=0;i<n&&i<MAXOBS;i++) {
        sat=obs[i].sat;//satellite num
        lam=nav->lam[sat-1];//carrier wave lengths
        if (lam[j/2]==0.0||lam[0]==0.0) continue;
        
        if ((r=geodist(rs+i*6,rr,e))<=0.0||
            satazel(pos,e,azel+i*2)<opt->elmin) {
            exc[i]=1;
			trace(4, "%s sat=%2d low ele el=%4.1f\n", str, sat,azel[1 + i * 2] * R2D);
            continue;
        }
        if (!(sys=satsys(sat,NULL))||!leo->ssat[sat-1].vs||
            satexclude(obs[i].sat,svh[i],opt)||exc[i]) {
            exc[i]=1;
			trace(4, "%s sat=%2d sat excluded or invalid\n", str, sat);
            continue;
        }
		/* tropospheric model*/ 
		if (leo->opt.mode<=PMODE_PPP_FIXED) {
			model_trop(obs[i].time, pos, azel + i * 2, opt, x, dtdx, nav, &dtrp, &vart);
			leo->ssat[sat-1].tropmap = dtdx[0];
		}
        /* ionospheric model */ 
        if (!model_iono(obs[i].time,pos,azel+i*2,opt,sat,x,nav,&dion,&vari)) {
			trace(4, "%s sat=%2d invalid trop or iono model\n", str, sat);
            continue;
        }
        /* satellite and receiver antenna model */  
		if (opt->posopt[0]) {
			satantpcv(rs + i * 6, rr, nav->pcvs + sat - 1, dants, &nadir);
			leo->ssat[sat - 1].nadir = nadir;
			for (k = 0; k < NFREQ; k++) {
				leo->ssat[sat - 1].pcvs[k] = dants[k];
			}
		}
		if (opt->mode <= PMODE_PPP_FIXED) {
			antmodel(opt->pcvr, opt->antdel[0], azel + i * 2, opt->posopt[1], dantr);
			for (k = 0; k < NFREQ; k++) {
				leo->ssat[sat - 1].pcvr[k] = dantr[k];
			}
		}
        /* phase windup model */   
        if (!model_phw(leo->sol.time,sat,nav->pcvs[sat-1].type,
                       opt->posopt[2]?2:0,rs+i*6,rr,&leo->ssat[sat-1].phw)) {
            continue;
        }
		dg = gravitycor(0, rr, rs + i * 6);
		leo->ssat[sat - 1].dg = dg; //gravitation correction
        /* corrected phase and code measurements */  
        corr_meas(obs+i,nav,azel+i*2,&leo->opt,dantr,dants,
                  leo->ssat[sat-1].phw,L,P,&Lc,&Pc);
        
        /* stack phase and code residuals {L1,P1,L2,P2,...} */
        for (j=0;j<2*NF(opt);j++) {
            
            dcb=bias=0.0;
            
            if (opt->ionoopt==IONOOPT_IFLC) {
                if ((y=j%2==0?Lc:Pc)==0.0) continue;
            }
            else {
                if ((y=j%2==0?L[j/2]:P[j/2])==0.0) continue;
                
                /* receiver DCB correction for P2 */
                if (j/2==1) dcb=-nav->rbias[0][sys==SYS_GLO?1:0][0];
            }
            C=SQR(lam[j/2]/lam[0])*ionmapf(pos,azel+i*2)*(j%2==0?-1.0:1.0);
            
            for (k=0;k<nx;k++) H[k+nx*nv]=k<3?-e[k]:0.0;
            
            /* receiver clock */
            k=sys==SYS_GLO?1:(sys==SYS_GAL?2:(sys==SYS_CMP?3:0));
            cdtr=x[IC(k,opt)];
            H[IC(k,opt)+nx*nv]=1.0;
            
            if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) {
                for (k=0;k<(opt->tropopt>=TROPOPT_ESTG?3:1);k++) {
                    H[IT(opt)+k+nx*nv]=dtdx[k];
                }
            }
            if (opt->ionoopt==IONOOPT_EST) {
                if (leo->x[II(sat,opt)]==0.0) continue;
                H[II(sat,opt)+nx*nv]=C;
            }
            if (j/2==2&&j%2==1) { /* L5-receiver-dcb */
                dcb+=leo->x[ID(opt)];
                H[ID(opt)+nx*nv]=1.0;
            }
            if (j%2==0) { /* phase bias */
                if ((bias=x[IB(sat,j/2,opt)])==0.0) continue;
                H[IB(sat,j/2,opt)+nx*nv]=1.0;
            }
            /* residual */  
            v[nv]=y-(r+cdtr-CLIGHT*dts[i*2]+dtrp+C*dion+dcb+bias-dg);
            
            if (j%2==0) leo->ssat[sat-1].resc[j/2]=v[nv];
            else        leo->ssat[sat-1].resp[j/2]=v[nv];
            
            /* variance */
            var[nv]=varerr(obs[i].sat,sys,azel[1+i*2],j/2,j%2,opt)+
                    vart+SQR(C)*vari+var_rs[i];
            if (sys==SYS_GLO&&j%2==1) var[nv]+=VAR_GLO_IFB;
            
            trace(3,"%s sat=%2d %s%d res=%9.4f sig=%9.4f el=%5.1f dantr=%11.8f %11.8f  dants=%11.8f %11.8f\n",str,sat,
                  j%2?"P":"L",j/2+1,v[nv],sqrt(var[nv]),azel[1+i*2]*R2D,dantr[0],dantr[1],dants[0],dants[1]);
            
            /* reject satellite by pre-fit residuals */
            if (!post&&opt->maxinno>0.0&&fabs(v[nv])>opt->maxinno) {
                trace(6,"outlier (%d) rejected %s sat=%2d %s%d res=%9.4f el=%4.1f\n",
                      post,str,sat,j%2?"P":"L",j/2+1,v[nv],azel[1+i*2]*R2D);
                exc[i]=1; leo->ssat[sat-1].rejc[j%2]++;
                continue;
            }
            /* record large post-fit residuals */  
            // if (post&&fabs(v[nv])>sqrt(var[nv])*THRES_REJECT) {
            //    obsi[ne]=i; frqi[ne]=j; ve[ne]=v[nv]; ne++;
            // }
            if (j%2==0) leo->ssat[sat-1].vsat[j/2]=1;
            nv++;
        }
    }
    /* reject satellite with large and max post-fit residual */
    if (post&&ne>0) {
        vmax=ve[0]; maxobs=obsi[0]; maxfrq=frqi[0]; rej=0;
        for (j=1;j<ne;j++) {
            if (fabs(vmax)>=fabs(ve[j])) continue;
            vmax=ve[j]; maxobs=obsi[j]; maxfrq=frqi[j]; rej=j;
        }
        sat=obs[maxobs].sat;
        trace(2,"outlier (%d) rejected %s sat=%2d %s%d res=%9.4f el=%4.1f\n",
            post,str,sat,maxfrq%2?"P":"L",maxfrq/2+1,vmax,azel[1+maxobs*2]*R2D);
        exc[maxobs]=1; leo->ssat[sat-1].rejc[maxfrq%2]++; stat=0;
        ve[rej]=0;
    }
    /* constraint to local correction */
	for (i = 0; i < n; i++) {
		sat = obs[i].sat;
		if (!leo->ssat[sat - 1].vsat[0]) continue;
		trace(6, "ppp_res: %s %d %2d P=%9.4f L=%9.4f phw=%8.3f tmf=%8.3f nadir=%9.4f pcvs=%9.4f %9.4f pcvr=%9.4f %9.4f  dg=%9.4f \n", 
			str, post,sat, leo->ssat[sat - 1].resp[0], leo->ssat[sat - 1].resc[0],
			leo->ssat[sat - 1].phw, leo->ssat[sat - 1].tropmap, leo->ssat[sat - 1].nadir, leo->ssat[sat - 1].pcvs[0], 
			leo->ssat[sat - 1].pcvs[1], leo->ssat[sat - 1].pcvr[0], leo->ssat[sat - 1].pcvr[1], leo->ssat[sat - 1].dg);
	}
    for (i=0;i<nv;i++) for (j=0;j<nv;j++) {
        R[i+j*nv]=i==j?var[i]:0.0;
    }
    return post?stat:nv;
}
/* number of estimated states ------------------------------------------------*/
extern int pppnx(const prcopt_t *opt)
{
    return NX(opt);
}
/* update solution status ----------------------------------------------------*/
static void updatesol(leo_t *leo, const obsd_t *obs, nav_t *nav,int n, int stat)
{
    const prcopt_t *opt=&leo->opt;
    int i,j;
	double dant[3] = { 0 };
    
    /* test # of valid satellites */
    leo->sol.ns=0;
    for (i=0;i<n&&i<MAXOBS;i++) {
        for (j=0;j<opt->nf;j++) {
            if (!leo->ssat[obs[i].sat-1].vsat[j]) continue;
            leo->ssat[obs[i].sat-1].lock[j]++;
            leo->ssat[obs[i].sat-1].outc[j]=0;
            if (j==0) leo->sol.ns++;
        }
    }
    leo->sol.stat=leo->sol.ns<MIN_NSAT_SOL?SOLQ_NONE:stat;
	if (leo->opt.mode > PMODE_PPP_FIXED) {
		leoantoff(leo->sol.time, leo->x, nav, dant); 
	}
    if (leo->sol.stat==SOLQ_FIX) {
        for (i=0;i<3;i++) {
            leo->sol.rr[i]=leo->xa[i]-dant[i];
            leo->sol.qr[i]=(float)leo->Pa[i+i*leo->na];
        }
        leo->sol.qr[3]=(float)leo->Pa[1];
        leo->sol.qr[4]=(float)leo->Pa[1+2*leo->na];
        leo->sol.qr[5]=(float)leo->Pa[2];
    }
    else {
        for (i=0;i<3;i++) {
            leo->sol.rr[i]=leo->x[i] - dant[i];
            leo->sol.qr[i]=(float)leo->P[i+i*leo->nx];
        }
        leo->sol.qr[3]=(float)leo->P[1];
        leo->sol.qr[4]=(float)leo->P[2+leo->nx];
        leo->sol.qr[5]=(float)leo->P[2];
    }
    leo->sol.dtr[0]=leo->x[IC(0,opt)];
    leo->sol.dtr[1]=leo->x[IC(1,opt)]-leo->x[IC(0,opt)];
    
    for (i=0;i<n&&i<MAXOBS;i++) for (j=0;j<opt->nf;j++) {
        leo->ssat[obs[i].sat-1].snr[j]=obs[i].SNR[j];
    }
    for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) {
        if (leo->ssat[i].slip[j]&3) leo->ssat[i].slipc[j]++;
        if (leo->ssat[i].fix[j]==2&&stat!=SOLQ_FIX) leo->ssat[i].fix[j]=1;
    }
}

/* precise point positioning -------------------------------------------------*/
extern int leopos(leo_t *leo, const obsd_t *obs, int n, const nav_t *nav)
{
    const prcopt_t *opt=&leo->opt;
    double *rs,*dts,*var,*v,*H,*R,*azel,*xp,*Pp,dr[3]={0};
    char str[32];
    int i,j,nv,info,svh[MAXOBS],exc[MAXOBS]={0},stat=SOLQ_SINGLE;
	gtime_t time;
	double *dant, *rel;
	gtime_t tutc;
	int nu;
	char msg[128] = "";

    time2str(obs[0].time,str,2);
    trace(3,"pppos   : time=%s nx=%d n=%d\n",str,leo->nx,n);
	trace(4, "obs=\n"); traceobs(4, obs, n);
	/* count leo observations */  
	for (nu = 0; nu <n&&obs[nu].rcv == 1; nu++);
	time = leo->sol.time; /* previous epoch */
	if (!pntpos(obs, nu, nav, &leo->opt, &leo->sol, NULL, leo->ssat, msg)) {
		errmsg(leo, "point pos error (%s)\n", msg);
		if (!leo->opt.dynamics) {  
			outsolstat(leo);
			return 0;
		}
	}
	/*if the single solution has something wrong,it will equit*/

	if (time.time != 0) leo->tt = timediff(obs->time, time);
    rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel=zeros(2,n);
	dant = mat(3, n); rel = mat(1, n);
    for (i=0;i<MAXSAT;i++) for (j=0;j<opt->nf;j++) leo->ssat[i].fix[j]=0;


	
	tutc = gpst2utc(obs->time);
	if (&nav->erp) geterp(&nav->erp, tutc, leo->erpv);


    /* temporal update of ekf states */  
    udstate_ppp(leo,obs,n,nav);
    
    /* satellite positions and clocks */
    satposs(1,obs[0].time,obs,n,nav,leo->opt.sateph,rs,dts,var,svh,dant,rel);
    
    /* exclude measurements of eclipsing satellite (block IIA) */  
    if (leo->opt.posopt[3]) {
        testeclipse(obs,n,nav,rs);
    }
    /* earth tides correction */ 
    if (opt->tidecorr) {
        tidedisp(gpst2utc(obs[0].time),leo->x,opt->tidecorr==1?1:7,&nav->erp,
                 opt->odisp[0],dr);
    }
    nv=n*leo->opt.nf*2+MAXSAT+3;
    xp=mat(leo->nx,1); Pp=zeros(leo->nx,leo->nx);
    v=mat(nv,1); H=mat(leo->nx,nv); R=mat(nv,nv);
    
    for (i=0;i<MAX_ITER;i++) {
        
        matcpy(xp,leo->x,leo->nx,1);
        matcpy(Pp,leo->P,leo->nx,leo->nx);
        
        /* prefit residuals */  
        if (!(nv=ppp_res(0,obs,n,rs,dts,var,svh,dr,exc,nav,xp,leo,v,H,R,azel))) {
            trace(2,"%s ppp (%d) no valid obs data\n",str,i+1);
            break;
        }
        /* measurement update of ekf states */ 
        if ((info=filter(xp,Pp,H,v,R,leo->nx,nv))) {
            trace(2,"%s ppp (%d) filter error info=%d\n",str,i+1,info);
            break;
        }
		
        /* postfit residuals */  
        if (ppp_res(i+1,obs,n,rs,dts,var,svh,dr,exc,nav,xp,leo,v,H,R,azel)) {
            matcpy(leo->x,xp,leo->nx,1);
            matcpy(leo->P,Pp,leo->nx,leo->nx);
            stat=SOLQ_PPP;
            break;
        }
    }
    if (i>=MAX_ITER) {
        trace(2,"%s ppp (%d) iteration overflows\n",str,i);
    }
	/* update solution status */
	updatesol(leo, obs,nav, n, stat);
	outsolstat(leo);
	trace(3, "xa=\t"); tracemat(3, leo->x, 1, leo->nx, 14, 4);
    free(rs); free(dts); free(var); free(azel);
    free(xp); free(Pp); free(v); free(H); free(R);
	return 1;
}
