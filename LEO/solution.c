/*------------------------------------------------------------------------------
* solution.c : solution functions
*
*          Copyright (C) 2007-2015 by T.TAKASU, All rights reserved.
*
* reference :
*     [1] National Marine Electronic Association and International Marine
*         Electronics Association, NMEA 0183 version 4.10, August 1, 2012
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2007/11/03  1.0 new
*           2009/01/05  1.1  add function outsols(), outsolheads(),
*                            setsolformat(), outsolexs, outsolex
*           2009/04/02  1.2  add dummy fields in NMEA mesassage
*                            fix bug to format lat/lon as deg-min-sec
*           2009/04/14  1.3  add age and ratio field to solution
*           2009/11/25  1.4  add function readsolstat()
*           2010/02/14  1.5  fix bug on output of gpstime at week boundary
*           2010/07/05  1.6  added api:
*                                initsolbuf(),freesolbuf(),addsol(),getsol(),
*                                inputsol(),outprcopts(),outprcopt()
*                            modified api:
*                                readsol(),readsolt(),readsolstat(),
*                                readsolstatt(),outsolheads(),outsols(),
*                                outsolexs(),outsolhead(),outsol(),outsolex(),
*                            deleted api:
*                                setsolopt(),setsolformat()
*           2010/08/14  1.7  fix bug on initialize solution buffer (2.4.0_p2)
*                            suppress enu-solution if base pos not available
*                            (2.4.0_p3)
*           2010/08/16  1.8  suppress null record if solution is not avalilable
*                            (2.4.0_p4)
*           2011/01/23  1.9  fix bug on reading nmea solution data
*                            add api freesolstatbuf()
*           2016/07/30  1.15 suppress output if std is over opt->maxsolstd
*-----------------------------------------------------------------------------*/
#include <ctype.h>
#include "leoorb.h"


/* constants and macros ------------------------------------------------------*/

#define SQR(x)     ((x)<0.0?-(x)*(x):(x)*(x))
#define SQRT(x)    ((x)<0.0?0.0:sqrt(x))

#define MAXFIELD   64           /* max number of fields in a record */
#define MAXNMEA    256          /* max length of nmea sentence */

#define KNOT2M     0.514444444  /* m/knot */

static const int solq_nmea[]={  /* nmea quality flags to leolib sol quality */
    /* nmea 0183 v.2.3 quality flags: */
    /*  0=invalid, 1=gps fix (sps), 2=dgps fix, 3=pps fix, 4=leo, 5=float leo */
    /*  6=estimated (dead reckoning), 7=manual input, 8=simulation */
    
    SOLQ_NONE ,SOLQ_SINGLE, SOLQ_DGPS, SOLQ_PPP , SOLQ_FIX,
    SOLQ_FLOAT,SOLQ_DR    , SOLQ_NONE, SOLQ_NONE, SOLQ_NONE
};
/* solution option to field separator ----------------------------------------*/
static const char *opt2sep(const solopt_t *opt)
{
    if (!*opt->sep) return " ";
    else if (!strcmp(opt->sep,"\\t")) return "\t";
    return opt->sep;
}
/* separate fields -----------------------------------------------------------*/
static int tonum(char *buff, const char *sep, double *v)
{
    int n,len=(int)strlen(sep);
    char *p,*q;
    
    for (p=buff,n=0;n<MAXFIELD;p=q+len) {
        if ((q=strstr(p,sep))) *q='\0'; 
        if (*p) v[n++]=atof(p);
        if (!q) break;
    }
    return n;
}
/* sqrt of covariance --------------------------------------------------------*/
static double sqvar(double covar)
{
    return covar<0.0?-sqrt(-covar):sqrt(covar);
}
/* convert ddmm.mm in nmea format to deg -------------------------------------*/
static double dmm2deg(double dmm)
{
    return floor(dmm/100.0)+fmod(dmm,100.0)/60.0;
}
/* convert time in nmea format to time ---------------------------------------*/
static void septime(double t, double *t1, double *t2, double *t3)
{
    *t1=floor(t/10000.0);
    t-=*t1*10000.0;
    *t2=floor(t/100.0);
    *t3=t-*t2*100.0;
}
/* solution to covariance ----------------------------------------------------*/
static void soltocov(const sol_t *sol, double *P)
{
    P[0]     =sol->qr[0]; /* xx or ee */
    P[4]     =sol->qr[1]; /* yy or nn */
    P[8]     =sol->qr[2]; /* zz or uu */
    P[1]=P[3]=sol->qr[3]; /* xy or en */
    P[5]=P[7]=sol->qr[4]; /* yz or nu */
    P[2]=P[6]=sol->qr[5]; /* zx or ue */
}
/* covariance to solution ----------------------------------------------------*/
static void covtosol(const double *P, sol_t *sol)
{
    sol->qr[0]=(float)P[0]; /* xx or ee */
    sol->qr[1]=(float)P[4]; /* yy or nn */
    sol->qr[2]=(float)P[8]; /* zz or uu */
    sol->qr[3]=(float)P[1]; /* xy or en */
    sol->qr[4]=(float)P[5]; /* yz or nu */
    sol->qr[5]=(float)P[2]; /* zx or ue */
}

/* test solution status ------------------------------------------------------*/
static int test_solstat(const char *buff)
{
    if (strlen(buff)<7||buff[0]!='$') return 0;
    return !strncmp(buff+1,"POS" ,3)||!strncmp(buff+1,"VELACC",6)||
           !strncmp(buff+1,"CLK" ,3)||!strncmp(buff+1,"ION"   ,3)||
           !strncmp(buff+1,"TROP",4)||!strncmp(buff+1,"HWBIAS",6)||
           !strncmp(buff+1,"TRPG",4)||!strncmp(buff+1,"AMB"   ,3)||
           !strncmp(buff+1,"SAT" ,3);
}

/* decode solution time ------------------------------------------------------*/
static char *decode_soltime(char *buff, const solopt_t *opt, gtime_t *time)
{
    double v[MAXFIELD];
    char *p,*q,s[64]=" ";
    int n,len;
    
    trace(4,"decode_soltime:\n");
    
    if (!strcmp(opt->sep,"\\t")) strcpy(s,"\t");
    else if (*opt->sep) strcpy(s,opt->sep);
    len=(int)strlen(s);
    
    if (opt->posf==SOLF_STAT) {
        return buff;
    }
    /* yyyy/mm/dd hh:mm:ss or yyyy mm dd hh:mm:ss */
    if (sscanf(buff,"%lf/%lf/%lf %lf:%lf:%lf",v,v+1,v+2,v+3,v+4,v+5)>=6) {
        if (v[0]<100.0) {
            v[0]+=v[0]<80.0?2000.0:1900.0;
        }
        *time=epoch2time(v);
        if (opt->times==TIMES_UTC) {
            *time=utc2gpst(*time);
        }
        if (!(p=strchr(buff,':'))||!(p=strchr(p+1,':'))) return NULL;
        for (p++;isdigit((int)*p)||*p=='.';) p++;
        return p+len;
    }

    /* wwww ssss */
    for (p=buff,n=0;n<2;p=q+len) {
        if ((q=strstr(p,s))) *q='\0'; 
        if (*p) v[n++]=atof(p);
        if (!q) break;
    }
    if (n>=2&&0.0<=v[0]&&v[0]<=3000.0&&0.0<=v[1]&&v[1]<604800.0) {
        *time=gpst2time((int)v[0],v[1]);
        return p;
    }
    return NULL;
}
/* decode x/y/z-ecef ---------------------------------------------------------*/
static int decode_solxyz(char *buff, const solopt_t *opt, sol_t *sol)
{
    double val[MAXFIELD],P[9]={0};
    int i=0,j,n;
    const char *sep=opt2sep(opt);
    
    trace(4,"decode_solxyz:\n");
    
    if ((n=tonum(buff,sep,val))<3) return 0;
    
    for (j=0;j<3;j++) {
        sol->rr[j]=val[i++]; /* xyz */
    }
    if (i<n) sol->stat=(unsigned char)val[i++];
    if (i<n) sol->ns  =(unsigned char)val[i++];
    if (i+3<n) {
        P[0]=val[i]*val[i]; i++; /* sdx */
        P[4]=val[i]*val[i]; i++; /* sdy */
        P[8]=val[i]*val[i]; i++; /* sdz */
        if (i+3<n) {
            P[1]=P[3]=SQR(val[i]); i++; /* sdxy */
            P[5]=P[7]=SQR(val[i]); i++; /* sdyz */
            P[2]=P[6]=SQR(val[i]); i++; /* sdzx */
        }
        covtosol(P,sol);
    }
    if (i<n) sol->age  =(float)val[i++];
    if (i<n) sol->ratio=(float)val[i];
    
    sol->type=0; /* postion type = xyz */
    
    if (MAXSOLQ<sol->stat) sol->stat=SOLQ_NONE;
    return 1;
}
/* decode lat/lon/height -----------------------------------------------------*/
static int decode_solllh(char *buff, const solopt_t *opt, sol_t *sol)
{
    double val[MAXFIELD],pos[3],Q[9]={0},P[9];
    int i=0,n;
    const char *sep=opt2sep(opt);
    
    trace(4,"decode_solllh:\n");
    
    n=tonum(buff,sep,val);
    
    if (!opt->degf) {
        if (n<3) return 0;
        pos[0]=val[i++]*D2R; /* lat/lon/hgt (ddd.ddd) */
        pos[1]=val[i++]*D2R;
        pos[2]=val[i++];
    }
    else {
        if (n<7) return 0;
        pos[0]=dms2deg(val  )*D2R; /* lat/lon/hgt (ddd mm ss) */
        pos[1]=dms2deg(val+3)*D2R;
        pos[2]=val[6];
        i+=7;
    }
    pos2ecef(pos,sol->rr);
    if (i<n) sol->stat=(unsigned char)val[i++];
    if (i<n) sol->ns  =(unsigned char)val[i++];
    if (i+3<n) {
        Q[4]=val[i]*val[i]; i++; /* sdn */
        Q[0]=val[i]*val[i]; i++; /* sde */
        Q[8]=val[i]*val[i]; i++; /* sdu */
        if (i+3<n) {
            Q[1]=Q[3]=SQR(val[i]); i++; /* sdne */
            Q[2]=Q[6]=SQR(val[i]); i++; /* sdeu */
            Q[5]=Q[7]=SQR(val[i]); i++; /* sdun */
        }
        covecef(pos,Q,P);
        covtosol(P,sol);
    }
    if (i<n) sol->age  =(float)val[i++];
    if (i<n) sol->ratio=(float)val[i];
    
    sol->type=0; /* postion type = xyz */
    
    if (MAXSOLQ<sol->stat) sol->stat=SOLQ_NONE;
    return 1;
}
/* decode e/n/u-baseline -----------------------------------------------------*/
static int decode_solenu(char *buff, const solopt_t *opt, sol_t *sol)
{
    double val[MAXFIELD],Q[9]={0};
    int i=0,j,n;
    const char *sep=opt2sep(opt);
    
    trace(4,"decode_solenu:\n");
    
    if ((n=tonum(buff,sep,val))<3) return 0;
    
    for (j=0;j<3;j++) {
        sol->rr[j]=val[i++]; /* enu */
    }
    if (i<n) sol->stat=(unsigned char)val[i++];
    if (i<n) sol->ns  =(unsigned char)val[i++];
    if (i+3<n) {
        Q[0]=val[i]*val[i]; i++; /* sde */
        Q[4]=val[i]*val[i]; i++; /* sdn */
        Q[8]=val[i]*val[i]; i++; /* sdu */
        if (i+3<n) {
            Q[1]=Q[3]=SQR(val[i]); i++; /* sden */
            Q[5]=Q[7]=SQR(val[i]); i++; /* sdnu */
            Q[2]=Q[6]=SQR(val[i]); i++; /* sdue */
        }
        covtosol(Q,sol);
    }
    if (i<n) sol->age  =(float)val[i++];
    if (i<n) sol->ratio=(float)val[i];
    
    sol->type=1; /* postion type = enu */
    
    if (MAXSOLQ<sol->stat) sol->stat=SOLQ_NONE;
    return 1;
}
/* decode solution status ----------------------------------------------------*/
static int decode_solsss(char *buff, sol_t *sol)
{
    double tow,pos[3],std[3]={0};
    int i,week,solq;
    
    trace(4,"decode_solssss:\n");
    
    if (sscanf(buff,"$POS,%d,%lf,%d,%lf,%lf,%lf,%lf,%lf,%lf",&week,&tow,&solq,
               pos,pos+1,pos+2,std,std+1,std+2)<6) {
        return 0;
    }
    if (week<=0||norm(pos,3)<=0.0||solq==SOLQ_NONE) {
        return 0;
    }
    sol->time=gpst2time(week,tow);
    for (i=0;i<6;i++) {
        sol->rr[i]=i<3?pos[i]:0.0;
        sol->qr[i]=i<3?(float)SQR(std[i]):0.0f;
        sol->dtr[i]=0.0;
    }
    sol->ns=0;
    sol->age=sol->ratio=sol->thres=0.0f;
    sol->type=0; /* position type = xyz */
    sol->stat=solq;
    return 1;
}

/* decode solution position --------------------------------------------------*/
static int decode_solpos(char *buff, const solopt_t *opt, sol_t *sol)
{
    sol_t sol0={{0}};
    char *p=buff;
    
    trace(4,"decode_solpos: buff=%s\n",buff);
    
    *sol=sol0;
    
    /* decode solution time */
    if (!(p=decode_soltime(p,opt,&sol->time))) {
        return 0;
    }
    /* decode solution position */
    switch (opt->posf) {
        case SOLF_XYZ : return decode_solxyz(p,opt,sol);
        case SOLF_LLH : return decode_solllh(p,opt,sol);
        case SOLF_ENU : return decode_solenu(p,opt,sol);
    }
    return 0;
}
/* decode reference position -------------------------------------------------*/
static void decode_refpos(char *buff, const solopt_t *opt, double *rb)
{
    double val[MAXFIELD],pos[3];
    int i,n;
    const char *sep=opt2sep(opt);
    
    trace(3,"decode_refpos: buff=%s\n",buff);
    
    if ((n=tonum(buff,sep,val))<3) return;
    
    if (opt->posf==SOLF_XYZ) { /* xyz */
        for (i=0;i<3;i++) rb[i]=val[i];
    }
    else if (opt->degf==0) { /* lat/lon/hgt (ddd.ddd) */
        pos[0]=val[0]*D2R;
        pos[1]=val[1]*D2R;
        pos[2]=val[2];
        pos2ecef(pos,rb);
    }
    else if (opt->degf==1&&n>=7) { /* lat/lon/hgt (ddd mm ss) */
        pos[0]=dms2deg(val  )*D2R;
        pos[1]=dms2deg(val+3)*D2R;
        pos[2]=val[6];
        pos2ecef(pos,rb);
    }
}
/* decode solution -----------------------------------------------------------*/
static int decode_sol(char *buff, const solopt_t *opt, sol_t *sol, double *rb)
{
    char *p;
    
    trace(4,"decode_sol: buff=%s\n",buff);
    

    if (test_solstat(buff)) { /* decode solution status */
        return decode_solsss(buff,sol);
    }
    if (!strncmp(buff,COMMENTH,1)) { /* reference position */
        if (!strstr(buff,"ref pos")&&!strstr(buff,"slave pos")) return 0;
        if (!(p=strchr(buff,':'))) return 0;
        decode_refpos(p+1,opt,rb);
        return 0;
    }
    /* decode position record */
    return decode_solpos(buff,opt,sol);
}
/* decode solution options ---------------------------------------------------*/
static void decode_solopt(char *buff, solopt_t *opt)
{
    char *p;
    
    trace(4,"decode_solhead: buff=%s\n",buff);
    
    if (strncmp(buff,COMMENTH,1)&&strncmp(buff,"+",1)) return;
    
    if      (strstr(buff,"GPST")) opt->times=TIMES_GPST;
    else if (strstr(buff,"UTC" )) opt->times=TIMES_UTC;
    
    if ((p=strstr(buff,"x-ecef(m)"))) {
        opt->posf=SOLF_XYZ;
        opt->degf=0;
        strncpy(opt->sep,p+9,1);
        opt->sep[1]='\0';
    }
    else if ((p=strstr(buff,"latitude(d'\")"))) {
        opt->posf=SOLF_LLH;
        opt->degf=1;
        strncpy(opt->sep,p+14,1);
        opt->sep[1]='\0';
    }
    else if ((p=strstr(buff,"latitude(deg)"))) {
        opt->posf=SOLF_LLH;
        opt->degf=0;
        strncpy(opt->sep,p+13,1);
        opt->sep[1]='\0';
    }
    else if ((p=strstr(buff,"e-baseline(m)"))) {
        opt->posf=SOLF_ENU;
        opt->degf=0;
        strncpy(opt->sep,p+13,1);
        opt->sep[1]='\0';
    }

}
/* read solution option ------------------------------------------------------*/
static void readsolopt(FILE *fp, solopt_t *opt)
{
    char buff[MAXSOLMSG+1];
    int i;
    
    trace(3,"readsolopt:\n");
    
    for (i=0;fgets(buff,sizeof(buff),fp)&&i<100;i++) { /* only 100 lines */
        
        /* decode solution options */
        decode_solopt(buff,opt);
    }
}


/* output solution as the form of x/y/z-ecef ---------------------------------*/
static int outecef(unsigned char *buff, const char *s, const sol_t *sol,
                   const solopt_t *opt)
{
    const char *sep=opt2sep(opt);
    char *p=(char *)buff;
    
    trace(3,"outecef:\n");
    
    p+=sprintf(p,"%s%s%14.4f%s%14.4f%s%14.4f%s%3d%s%3d%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%6.2f%s%6.1f\n",
               s,sep,sol->rr[0],sep,sol->rr[1],sep,sol->rr[2],sep,sol->stat,sep,
               sol->ns,sep,SQRT(sol->qr[0]),sep,SQRT(sol->qr[1]),sep,SQRT(sol->qr[2]),
               sep,sqvar(sol->qr[3]),sep,sqvar(sol->qr[4]),sep,sqvar(sol->qr[5]),
               sep,sol->age,sep,sol->ratio);
    return p-(char *)buff;
}
/* output solution as the form of lat/lon/height -----------------------------*/
static int outpos(unsigned char *buff, const char *s, const sol_t *sol,
                  const solopt_t *opt)
{
    double pos[3],dms1[3],dms2[3],P[9],Q[9];
    const char *sep=opt2sep(opt);
    char *p=(char *)buff;
    
    trace(3,"outpos  :\n");
    
    ecef2pos(sol->rr,pos);
    soltocov(sol,P);
    covenu(pos,P,Q);
    if (opt->height==1) { /* geodetic height */
		pos[2] -= 0;
    }
    if (opt->degf) {
        deg2dms(pos[0]*R2D,dms1,5);
        deg2dms(pos[1]*R2D,dms2,5);
        p+=sprintf(p,"%s%s%4.0f%s%02.0f%s%08.5f%s%4.0f%s%02.0f%s%08.5f",s,sep,
                   dms1[0],sep,dms1[1],sep,dms1[2],sep,dms2[0],sep,dms2[1],sep,
                   dms2[2]);
    }
    else p+=sprintf(p,"%s%s%14.9f%s%14.9f",s,sep,pos[0]*R2D,sep,pos[1]*R2D);
    p+=sprintf(p,"%s%10.4f%s%3d%s%3d%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%6.2f%s%6.1f\n",
               sep,pos[2],sep,sol->stat,sep,sol->ns,sep,SQRT(Q[4]),sep,
               SQRT(Q[0]),sep,SQRT(Q[8]),sep,sqvar(Q[1]),sep,sqvar(Q[2]),
               sep,sqvar(Q[5]),sep,sol->age,sep,sol->ratio);
    return p-(char *)buff;
}
/* output solution as the form of e/n/u-baseline -----------------------------*/
static int outenu(unsigned char *buff, const char *s, const sol_t *sol,
                  const double *rb, const solopt_t *opt)
{
    double pos[3],rr[3],enu[3],P[9],Q[9];
    int i;
    const char *sep=opt2sep(opt);
    char *p=(char *)buff;
    
    trace(3,"outenu  :\n");
    
    for (i=0;i<3;i++) rr[i]=sol->rr[i]-rb[i];
    ecef2pos(rb,pos);
    soltocov(sol,P);
    covenu(pos,P,Q);
    ecef2enu(pos,rr,enu);
    p+=sprintf(p,"%s%s%14.4f%s%14.4f%s%14.4f%s%3d%s%3d%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%8.4f%s%6.2f%s%6.1f\n",
               s,sep,enu[0],sep,enu[1],sep,enu[2],sep,sol->stat,sep,sol->ns,sep,
               SQRT(Q[0]),sep,SQRT(Q[4]),sep,SQRT(Q[8]),sep,sqvar(Q[1]),
               sep,sqvar(Q[5]),sep,sqvar(Q[2]),sep,sol->age,sep,sol->ratio);
    return p-(char *)buff;
}
static int outsp3(unsigned char *buff, const sol_t *sol)
{
	char *p = (char *)buff;
	double ep[6] = { 0.0 };
	time2epoch(sol->time, ep);
	p += sprintf(p, "*  %4d %2d %2d %2d %2d %11.8f\nPL09%14.6f%14.6f%14.6f%14.6f %2d %2d %2d %3d\n",
		(int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], round(ep[5]), sol->rr[0] / 1000.0, sol->rr[1] / 1000.0, sol->rr[2] / 1000.0, sol->dtr[0]/CLIGHT*1e6, 0, 0, 0, 0);
	return p - (char *)buff;
}
/* output processing options ---------------------------------------------------
* output processing options to buffer
* args   : unsigned char *buff IO output buffer
*          prcopt_t *opt    I   processign options
* return : number of output bytes
*-----------------------------------------------------------------------------*/
extern int outprcopts(unsigned char *buff, const prcopt_t *opt)
{
    const int sys[]={SYS_GPS,SYS_GLO,SYS_GAL,0};
    const char *s1[]={"single","dgps","kinematic","static","moving-base","fixed",
                 "ppp-kinematic","ppp-static","ppp-fixed","leo-kin","leo-rdn"};
    const char *s2[]={"L1","L1+L2","L1+L2+L5","L1+L2+L5+L6","L1+L2+L5+L6+L7",
                      "L1+L2+L5+L6+L7+L8",""};
    const char *s3[]={"forward","backward","combined"};
    const char *s4[]={"off","broadcast","sbas","iono-free","estimation",
                      "ionex tec","qzs","lex","vtec_sf","vtec_ef","gtec",""};
    const char *s5[]={"off","saastamoinen","sbas","est ztd","est ztd+grad",""};
    const char *s6[]={"broadcast","precise","broadcast+sbas","broadcast+ssr apc",
                      "broadcast+ssr com","qzss lex",""};
    const char *s7[]={"gps","glonass","galileo","qzss","sbas",""};
    const char *s8[]={"off","continuous","instantaneous","fix and hold",""};
    const char *s9[]={"off","on","auto calib","external calib",""};
    int i;
    char *p=(char *)buff;
    
    trace(3,"outprcopts:\n");
    
    p+=sprintf(p,"%s pos mode  : %s\n",COMMENTH,s1[opt->mode]);
    
    if (PMODE_DGPS<=opt->mode&&opt->mode<=PMODE_FIXED) {
        p+=sprintf(p,"%s freqs     : %s\n",COMMENTH,s2[opt->nf-1]);
    }
    if (opt->mode>PMODE_SINGLE) {
        p+=sprintf(p,"%s solution  : %s\n",COMMENTH,s3[opt->soltype]);
    }
    p+=sprintf(p,"%s elev mask : %.1f deg\n",COMMENTH,opt->elmin*R2D);
    if (opt->mode>PMODE_SINGLE) {
        p+=sprintf(p,"%s dynamics  : %s\n",COMMENTH,opt->dynamics?"on":"off");
        p+=sprintf(p,"%s tidecorr  : %s\n",COMMENTH,opt->tidecorr?"on":"off");
    }
    if (opt->mode<=PMODE_FIXED) {
        p+=sprintf(p,"%s ionos opt : %s\n",COMMENTH,s4[opt->ionoopt]);
    }
    p+=sprintf(p,"%s tropo opt : %s\n",COMMENTH,s5[opt->tropopt]);
    p+=sprintf(p,"%s ephemeris : %s\n",COMMENTH,s6[opt->sateph]);
    if (opt->navsys!=SYS_GPS) {
        p+=sprintf(p,"%s navi sys  :",COMMENTH);
        for (i=0;sys[i];i++) {
            if (opt->navsys&sys[i]) p+=sprintf(p," %s",s7[i]);
        }
        p+=sprintf(p,"\n");
    }
    if (PMODE_KINEMA<=opt->mode&&opt->mode<=PMODE_FIXED) {
        p+=sprintf(p,"%s amb res   : %s\n",COMMENTH,s8[opt->modear]);
        if (opt->navsys&SYS_GLO) {
            p+=sprintf(p,"%s amb glo   : %s\n",COMMENTH,s9[opt->glomodear]);
        }
        if (opt->thresar[0]>0.0) {
            p+=sprintf(p,"%s val thres : %.1f\n",COMMENTH,opt->thresar[0]);
        }
    }
    if (opt->mode==PMODE_MOVEB&&opt->baseline[0]>0.0) {
        p+=sprintf(p,"%s baseline  : %.4f %.4f m\n",COMMENTH,
                   opt->baseline[0],opt->baseline[1]);
    }
    for (i=0;i<2;i++) {
        if (opt->mode==PMODE_SINGLE||(i>=1&&opt->mode>PMODE_FIXED)) continue;
        p+=sprintf(p,"%s antenna%d  : %-21s (%7.4f %7.4f %7.4f)\n",COMMENTH,
                   i+1,opt->anttype[i],opt->antdel[i][0],opt->antdel[i][1],
                   opt->antdel[i][2]);
    }
    return p-(char *)buff;
}
/* output solution header ------------------------------------------------------
* output solution header to buffer
* args   : unsigned char *buff IO output buffer
*          solopt_t *opt    I   solution options
* return : number of output bytes
*-----------------------------------------------------------------------------*/
extern int outsolheads(unsigned char *buff, const solopt_t *opt)
{
    const char *s2[]={"ellipsoidal","geodetic"};
    const char *s3[]={"GPST","UTC "},*sep=opt2sep(opt);
    char *p=(char *)buff;
    int timeu=opt->timeu<0?0:(opt->timeu>20?20:opt->timeu);
    
    trace(3,"outsolheads:\n");
    
    if (opt->posf==SOLF_STAT) {
        return 0;
    }
    if (opt->outhead&&opt->posf!=SOLF_SP3) {
        p+=sprintf(p,"%s (",COMMENTH);
        if      (opt->posf==SOLF_XYZ) p+=sprintf(p,"x/y/z-ecef=WGS84");
        else if (opt->posf==SOLF_ENU) p+=sprintf(p,"e/n/u-baseline=WGS84");
        else p+=sprintf(p,"lat/lon/height=WGS84/%s",s2[opt->height]);
        p+=sprintf(p,",Q=1:fix,2:float,3:sbas,4:dgps,5:single,6:ppp,ns=# of satellites)\n");
    }
    p+=sprintf(p,"%s  %-*s%s",COMMENTH,(opt->timef?16:8)+timeu+1,s3[opt->times],sep);
    
    if (opt->posf==SOLF_LLH) { /* lat/lon/hgt */
        if (opt->degf) {
            p+=sprintf(p,"%16s%s%16s%s%10s%s%3s%s%3s%s%8s%s%8s%s%8s%s%8s%s%8s%s%8s%s%6s%s%6s\n",
                       "latitude(d'\")",sep,"longitude(d'\")",sep,"height(m)",sep,
                       "Q",sep,"ns",sep,"sdn(m)",sep,"sde(m)",sep,"sdu(m)",sep,
                       "sdne(m)",sep,"sdeu(m)",sep,"sdue(m)",sep,"age(s)",sep,"ratio");
        }
        else {
            p+=sprintf(p,"%14s%s%14s%s%10s%s%3s%s%3s%s%8s%s%8s%s%8s%s%8s%s%8s%s%8s%s%6s%s%6s\n",
                       "latitude(deg)",sep,"longitude(deg)",sep,"height(m)",sep,
                       "Q",sep,"ns",sep,"sdn(m)",sep,"sde(m)",sep,"sdu(m)",sep,
                       "sdne(m)",sep,"sdeu(m)",sep,"sdun(m)",sep,"age(s)",sep,"ratio");
        }
    }
    else if (opt->posf==SOLF_XYZ) { /* x/y/z-ecef */
        p+=sprintf(p,"%14s%s%14s%s%14s%s%3s%s%3s%s%8s%s%8s%s%8s%s%8s%s%8s%s%8s%s%6s%s%6s\n",
                   "x-ecef(m)",sep,"y-ecef(m)",sep,"z-ecef(m)",sep,"Q",sep,"ns",sep,
                   "sdx(m)",sep,"sdy(m)",sep,"sdz(m)",sep,"sdxy(m)",sep,
                   "sdyz(m)",sep,"sdzx(m)",sep,"age(s)",sep,"ratio");
    }
    else if (opt->posf==SOLF_ENU) { /* e/n/u-baseline */
        p+=sprintf(p,"%14s%s%14s%s%14s%s%3s%s%3s%s%8s%s%8s%s%8s%s%8s%s%8s%s%8s%s%6s%s%6s\n",
                   "e-baseline(m)",sep,"n-baseline(m)",sep,"u-baseline(m)",sep,
                   "Q",sep,"ns",sep,"sde(m)",sep,"sdn(m)",sep,"sdu(m)",sep,
                   "sden(m)",sep,"sdnu(m)",sep,"sdue(m)",sep,"age(s)",sep,"ratio");
    }
	else if (opt->posf == SOLF_SP3) { /*SP3 */
		p += sprintf(p, "#cP%4d %2d %02d %2d %2d %11.8f      96 ORBIT IGS05 HLM  IGS\n", 2017, 1, 1, 0, 0, 0.0);
		p += sprintf(p, "## 1594      0.00000000   900.00000000 55402 0.0000000000000\n");
		p += sprintf(p, "+    1   L09  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n");
		p += sprintf(p, "+          0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n");
		p += sprintf(p, "+          0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n");
		p += sprintf(p, "+          0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n");
		p += sprintf(p, "+          0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n");
		p += sprintf(p, "++         0  2  2  2  2  2  2  2  2  2  2  2  2  2  3  2  2\n");
		p += sprintf(p, "++         2  2  2  2  2  2  2  4  3  2  2  2  2  2  2  0  0\n");
		p += sprintf(p, "++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n");
		p += sprintf(p, "++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n");
		p += sprintf(p, "++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n");
		p += sprintf(p, "%%c L  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n");
		p += sprintf(p, "%%c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n");
		p += sprintf(p, "%%f  1.2500000  1.025000000  0.00000000000  0.000000000000000\n");
		p += sprintf(p, "%%f  0.0000000  0.000000000  0.00000000000  0.000000000000000\n");
		p += sprintf(p, "%%i    0    0    0    0      0      0      0      0         0\n");
		p += sprintf(p, "%%i    0    0    0    0      0      0      0      0         0\n");
		p += sprintf(p, "/* FINAL ORBIT COMBINATION FROM WEIGHTED AVERAGE OF:         \n");
		p += sprintf(p, "/* cod emr esa gfz grg jpl mit ngs sio                       \n");
		p += sprintf(p, "/* REFERENCED TO IGS TIME (IGST) AND TO WEIGHTED MEAN POLE:  \n");
		p += sprintf(p, "/* PCV:IGS05_1585 OL/AL:FES2004  NONE     Y  ORB:CMB CLK:CMB \n");
	}
    return p-(char *)buff;
}
/* std-dev of soltuion -------------------------------------------------------*/
static double sol_std(const sol_t *sol)
{
    /* approximate as max std-dev of 3-axis std-devs */
    if (sol->qr[0]>sol->qr[1]&&sol->qr[0]>sol->qr[2]) return SQRT(sol->qr[0]);
    if (sol->qr[1]>sol->qr[2]) return SQRT(sol->qr[1]);
    return SQRT(sol->qr[2]);
}
/* output solution body --------------------------------------------------------
* output solution body to buffer
* args   : unsigned char *buff IO output buffer
*          sol_t  *sol      I   solution
*          double *rb       I   base station position {x,y,z} (ecef) (m)
*          solopt_t *opt    I   solution options
* return : number of output bytes
*-----------------------------------------------------------------------------*/
extern int outsols(unsigned char *buff, const sol_t *sol, const double *rb,
                   const solopt_t *opt)
{
    gtime_t time,ts={0};
    double gpst;
    int week,timeu;
    const char *sep=opt2sep(opt);
    char s[64];
    unsigned char *p=buff;
    
    trace(3,"outsols :\n");
    
    /* suppress output if std is over opt->maxsolstd */
    if (opt->maxsolstd>0.0&&sol_std(sol)>opt->maxsolstd) {
        return 0;
    }
    timeu=opt->timeu<0?0:(opt->timeu>20?20:opt->timeu);
    
    time=sol->time;
    if (opt->times>=TIMES_UTC) time=gpst2utc(time);
    
    if (opt->timef) time2str(time,s,timeu);
    else {
        gpst=time2gpst(time,&week);
        if (86400*7-gpst<0.5/pow(10.0,timeu)) {
            week++;
            gpst=0.0;
        }
        sprintf(s,"%4d%s%*.*f",week,sep,6+(timeu<=0?0:timeu+1),timeu,gpst);
    }
    switch (opt->posf) {
        case SOLF_LLH:  p+=outpos (p,s,sol,opt);   break;
        case SOLF_XYZ:  p+=outecef(p,s,sol,opt);   break;
        case SOLF_ENU:  p+=outenu(p,s,sol,rb,opt); break;
		case SOLF_SP3:  p += outsp3(p,sol);		   break;
    }
    return p-buff;
}

/* output processing option ----------------------------------------------------
* output processing option to file
* args   : FILE   *fp       I   output file pointer
*          prcopt_t *opt    I   processing options
* return : none
*-----------------------------------------------------------------------------*/
extern void outprcopt(FILE *fp, const prcopt_t *opt)
{
    unsigned char buff[MAXSOLMSG+1];
    int n;
    
    trace(3,"outprcopt:\n");
    
    if ((n=outprcopts(buff,opt))>0) {
        fwrite(buff,n,1,fp);
    }
}
/* output solution header ------------------------------------------------------
* output solution heade to file
* args   : FILE   *fp       I   output file pointer
*          solopt_t *opt    I   solution options
* return : none
*-----------------------------------------------------------------------------*/
extern void outsolhead(FILE *fp, const solopt_t *opt)
{
    unsigned char buff[MAXSOLMSG+1];
    int n;
    
    trace(3,"outsolhead:\n");
	if ((n = outsolheads(buff, opt)) > 0) {
		fwrite(buff, n, 1, fp);
	}
}
/* output solution body --------------------------------------------------------
* output solution body to file
* args   : FILE   *fp       I   output file pointer
*          sol_t  *sol      I   solution
*          double *rb       I   base station position {x,y,z} (ecef) (m)
*          solopt_t *opt    I   solution options
* return : none
*-----------------------------------------------------------------------------*/
extern void outsol(FILE *fp, const sol_t *sol, const double *rb,
                   const solopt_t *opt)
{
    unsigned char buff[MAXSOLMSG+1];
    int n;
    
    trace(3,"outsol  :\n");
    
    if ((n=outsols(buff,sol,rb,opt))>0) {
        fwrite(buff,n,1,fp);
    }
}
