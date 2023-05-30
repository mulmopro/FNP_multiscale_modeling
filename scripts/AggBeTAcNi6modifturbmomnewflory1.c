/*UDF for Fluent 6.3*/
/*PCL precipitation*/
/*version 1.1_test2 - Mixing, Nucleation and Growth*/

/*  User Defined Scalars : (23)*/

/*  0:weight of abscissa 1 */
/*  1:weighted abscissa 1 */
/*  2:weighted abscissa 2 */
/*  3:abscissa 1 ---------------------------------------NOT RESOLVED */
/*  4:abscissa 2 ---------------------------------------NOT RESOLVED */
/*  5:*/
/*  6:*/
/*  7:weighted order 0 moment in node 1 */
/*  8:weighted order 1 moment in node 1 */
/*  9:weighted order 2 moment in node 1 */
/*  10:weighted order 3 moment in node 1 */
/*  11:weighted order 0 moment in node 2 */
/*  12:weighted order 1 moment in node 2 */
/*  13:weighted order 2 moment in node 2 */
/*  14:weighted order 3 moment in node 2 */
/*  15:order 0 moment in node 1 --------------------------------NOT RESOLVED */
/*  16:order 1 moment in node 1 --------------------------------NOT RESOLVED */
/*  17:order 2 moment in node 1 --------------------------------NOT RESOLVED */
/*  18:order 3 moment in node 1 --------------------------------NOT RESOLVED */
/*  19:order 0 moment in node 2 --------------------------------NOT RESOLVED */
/*  20:order 1 moment in node 2 --------------------------------NOT RESOLVED */
/*  21:order 2 moment in node 2 --------------------------------NOT RESOLVED */
/*  22:order 3 moment in node 2 --------------------------------NOT RESOLVED */
/*  23:weighted order 4 moment in node 1 */
/*  24:weighted order 5 moment in node 1 */
/*  25:weighted order 4 moment in node 2 */
/*  26:weighted order 5 moment in node 2 */
/*  27:order 4 moment in node 1 --------------------------------NOT RESOLVED */
/*  28:order 5 moment in node 1 --------------------------------NOT RESOLVED */
/*  29:order 4 moment in node 2 --------------------------------NOT RESOLVED */
/*  30:order 5 moment in node 2 --------------------------------NOT RESOLVED */




/*  User Defined Memories : (50)*/

/*  0:mean mixture fraction */
/*  1:mean mixture fraction variance */
/*  2:viscosità cinematica*/
/*  3:source term for the weighted 1st & 2nd nodes */
/*  4:concentration of PCL in node 1, c_1, kmol m-3 */
/*  5:concentration of PCL in node 2, c_2, kmol m-3 */
/*  6:equilibrium concentration of PCL in node 1, c_e1, kmol m-3 */
/*  7:equilibrium concentration of PCL in node 2, c_e2, kmol m-3 */
/*  8:relative supersaturation [-] */
/*  9:nucleation rate [# m-3 s-1]*/
/*  10:mean concentration of PCL, c_pcl, kmol m-3 */
/*  11:abscissa 1 in node 1 [m] */
/*  12:abscissa 2 in node 1 [m] */
/*  13:weight 1 in node 1 - */
/*  14:weight 2 in node 1 - */
/*  15:abscissa 1 in node 2 [m] */
/*  16:abscissa 2 in node 2 [m] */
/*  17:weight 1 in node 2 - */
/*  18:weight 2 in node 2 - */
/*  19:source term of moment of order 0 in node 1 */
/*  20:source term of moment of order 1 in node 1 */
/*  21:source term of moment of order 2 in node 1 */
/*  22:source term of moment of order 3 in node 1 */
/*  23:source term of moment of order 0 in node 2 */
/*  24:source term of moment of order 1 in node 2 */
/*  25:source term of moment of order 2 in node 2 */
/*  26:source term of moment of order 3 in node 2 */
/*  27:moment of order 0 [m-3] */
/*  28:moment of order 1 [m-2] */
/*  29:moment of order 2 [m-1] */
/*  30:moment of order 3 - */
/*  31:volume-average diameter d43 [m] è in nm!!*/
/*  32:Pe number in environment 1 [-] */
/*  33:Pe number in environment 2[-] */
/*  34:Average Pe number [-]*/
/*  35:nucleation rate in node 2 [# m-3 s-1]*/
/*  36:local Reynolds number*/
/*  37:C_phi, function of Rel*/
/*  38:micromixing rate, gamma*/
/*  39:surface tension in node 1*/
/*  40:surface tension in node 2*/
/*  41:water molar fraction in node 1*/
/*  42:water molar fraction in node 2*/
/*  43:size of nucleated particles in node 1 [m]*/
/*  44:size of nucleated particles in node 2 [m]*/
/*  45:growth rate ambient 1- node 1*/
/*  46:growth rate ambient 1- node 2*/
/*  47:growth rate ambient 2- node 1*/
/*  48:growth rate ambient 2- node 2*/
/*  49:moment of order 4 */
/*  50:mean water fraction*/
/*  51:mean equilibrium concentration of PCL [kmol/m-3] errore nel'articolo Di pasquale 2012*/
/*  70 scala di Kolmogorov */
/*  71-78:turb/brown */
/*  79: initial supersaturation*/
/*  80: Stoke−Einstein diffusion coefficient [m2/s] */
/*  81: S-E diffusion coefficient per singola molecola PCL [m2/s] */
/*  82: kernel browniano di aggregazione [m3/s] */
/*  83: termine presente nel kernel turbolento [s-2] */
/*  84: kernel turbolento di aggregazione tra due singole molecole [m3/s] */ 


#include "udf.h"
/* #include "sg.h" */
/* #include "flow.h" */
/* #include "models.h" */
/* #include <math.h> */

#define node 2           /* number of nodes in quadrature approximation of the PBE */
#define NA 6.0221415e+23 /*Avogadro Number*/
#define u_max 1000
#define u_min 1

/* turbulent schmidt number for mixture fraction calculation */
real Sct   = 0.7; /* default value for mf */
/* Sct = RP_Get_Real("species/sct"); */

/* constant for turbulent viscosity (Boussinesq) */
/*real C_mu  = 0.09;*/

/* round to zero constant */
real eps   = 1.0e-5;/* <=1e-4 to have low residuals */
real eps2   = 1.0e-9;

/* initial concentration of PCL, mg/ml da cambiare nel caso!!!!!!!!!!!!!!!!!!!!!*/
real c_in = 6;

/*molecular volume, [m3]
real v_m = 2.06e-28;*/
/*molecular diameter, d_m = 2*pow((3/4/3.14159*v_m),(1/3)), [m]*/
/*real d_m = 7.490e-10;*/

/*Polycaprolactone properties*/
real M_w = 14000; /*kg/kmol*/
real rho_pcl = 1146; /*kg/m3*/

/*Boltzmann constant, J/K*/
real k_B = 1.3805e-23;

/*Molecular Weight Monomer*/
real MM_w = 114; /*kg/kmol*/

/*Temperature, K*/
real T = 303;

/* molar volume of acetonitrile [m^3/mol] !!!!!!!!!!!!!!!!!!!!!!*/
real VolA = 0.000053214;		/*dal Perry 8th ed (vedi file excel)*/

/* molar volume of PCL [m^3/mol]*/
real VolPCL = 0.012216405;

/* molar volume of water [m^3/mol]*/
real VolW = 0.000018115;

/* density of acetonitrile [kg/m^3!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
real rhoA = 771.45;		/*dal Perry 8th ed. (vedi file excel)*/

/* density of water [kg/m^3]*/
real rhoW = 993.68;

/* viscosity of acetonitrile [Pa*s]!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
real muA = 0.000326;		/*da Perry 8th ed. (vedi file excel)*/


/* viscosity of water [Pa*s]*/
real muW = 0.00085;

/* pi greco */
real pi = 3.1415926535897932384626433;

/* turbulent aggregation kernel constant */
real turb = 1.2944;


/* The following function uses the Product-Difference algorithm of Gordon */
/* On input mom[2*node] is a vector of "node" moments. */
/* On output mom[2*node] is rearranged as a unique vector that consists */
/* of diagonal elements d[node] and subdiagonal elements e[node] of the*/
/* Jacobi matrix. */
void PDalg(real mom[])
{
  real alpha[2*node+1], p[2*node+2][2*node+2];
  int i,j;
  
  for (j=1;j<=2*node+1;j++) {
    for (i=1;i<=2*node+1;i++) {
      p[i][j] = 0.;
    }
  }

  for (i=1;i<=2*node+1;i++) {
    if (i==1) {p[1][i] = 1.;}
    else {p[1][i] = 0.;
    }
  }

  for (i=1;i<=2*node;i++) {
    p[2][i] = pow(-1.,i-1)*mom[i];
  }

  for (j=3;j<=2*node+1;j++) {
    for (i=1;i<=2*node+2-j;i++) {
      p[j][i] = p[j-1][1]*p[j-2][i+1] - p[j-2][1]*p[j-1][i+1];
    }
  }


  alpha[1]=0.;
  for (i=2;i<=2*node;i++) {
    if ( p[i][1]*p[i-1][1] != 0.) {
      alpha[i] = p[i+1][1]/(p[i][1]*p[i-1][1]);
    }
    else {alpha[i] = 0.;}
  }


  for (i=1;i<=2*node;i++){
    mom[i]=0.;
  }

  for (i=1;i<=node;i++)
    {
      mom[i]=alpha[2*i]+alpha[2*i-1];
    }

  for (i=1;i<=node-1;i++)
    {
      mom[node+i]= -sqrt(fabs(alpha[2*i+1]*alpha[2*i]));
    }

  mom[2*node]=0.0;

}



#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
real sqrarg;
#define SQR2(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

real pythag(real a, real b)
{
  real absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQR2(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR2(absa/absb)));
}



/* The following routine function has been taken from Numerical recipes*/
/* and slightly modified */
/* Its aim is to diagonalize a tridiagonal matrix */
void tqli(real d[], real e[], int n, real z[node+1][node+1])
{
  real pythag(real a, real b);
  int m,l,iter,i,k;
  real s,r,p,g,f,dd,c,b;


  for (l=1;l<=n;l++) {
    iter=0;
    do {
      for (m=l;m<=n-1;m++) {
    dd=fabs(d[m])+fabs(d[m+1]);
    if ((fabs(e[m])+dd) == dd) break;
      }
      if (m != l) {
    if (iter == 30) break;
    iter++;
    g=(d[l+1]-d[l])/(2.0*e[l]);
    r=pythag(g,1.0);
    g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
    s=c=1.0;
    p=0.0;
    for (i=m-1;i>=l;i--) {
      f=s*e[i];
      b=c*e[i];
      e[i+1]=(r=pythag(f,g));
      if (r == 0.0) {
        d[i+1] -= p;
        e[m]=0.0;
        break;
      }
      s=f/r;
      c=g/r;
      g=d[i+1]-p;
      r=(d[i]-g)*s+2.0*c*b;
      d[i+1]=g+(p=s*r);
      g=c*r-b;
      for (k=1;k<=n;k++) {
        f=z[k][i+1];
        z[k][i+1]=s*z[k][i]+c*f;
        z[k][i]=c*z[k][i]-s*f;
      }
    }
    if (r == 0.0 && i >= l) continue;
    d[l] -= p;
    e[l]=g;
    e[m]=0.0;
      }
    } while (m != l);
  }
}

/* Funzione che restituisce il valore del coefficiente di Flory*/
/* I coefficienti costanti sono i valori ottenuti interpolando i risultati delle simulazioni MD*/
/* IMPORTANTE: la legge è stata ottenuta considerando la frazione molare di acetone
 dato che il programma calcola la frazione molare di acqua, devo usare (1-x) WARNING*/

real nuF(real x){

    real n;

   /* n= -0.15*(1-x)*(1-x) +0.45*(1-x)+0.30;  legge modificata secondo la pubblicazione di march_dipasq_altri_2014
                                                         vers_precedente : n= -0.19514*(1-x)*(1-x) +0.51256*(1-x)+0.30844; */
      n = -0.10*(1-x)*(1-x)+0.40*(1-x)+0.30; /*nuova legge per nu: vedi file excel marco*/
 		
    return n;
}

/* Funzione che restituisce il valore del coefficiente moltiplicativo della legge di Flory*/

real kF(real x){

    real k;

    k= 0.86*(0.0064*exp(-3.15*(1-x)));/* legge modificata secondo la pubblicazione di march_dipasq_altri_2014
                                      vers_precedente :k= 0.00612*exp(-4.0325*(1-x));
					Marco: riduzione da relazione teorica */
    return k;
}

/***********************************************************************************************/
/* DEFINITIONS OF PROPERTIES (Fluent macros)*/

/*definizione della densità in funzione della frazione di miscelazione*/
DEFINE_PROPERTY(density,c,t)
{
  real rho;

  rho = C_UDMI(c,t,0)*rhoA+(1- C_UDMI(c,t,0))*rhoW;/*in base a eq.8 articolo alessio utilizzando la frazione in volume di A*/

  return rho;
}

/*definizione della viscosità in funzione della frazione di miscelazione*/
DEFINE_PROPERTY(viscosity,c,t)
{
  real mu;
  real xA;

  xA= C_UDMI(c,t,0)/VolA/(C_UDMI(c,t,0)/VolA+(1-C_UDMI(c,t,0))/VolW);/* utilizzando la frazione in volume di A*/
  mu = exp(xA*log(muA)+(1-xA)*log(muW));/*eq.9 articolo alessio*/

  return mu;
}


DEFINE_DIFFUSIVITY(TURB_DIFF,c,t,i)
{
/* Note: the diffusivity in Fluent GUI requires not D, but rho*D 
as well as C_UDS_DIFF(c,t,i)*/

real d_turb, d_eff;

/* Laminar diff coeff */
d_eff = 1.0e-6;
d_eff *= C_R(c,t);

if (rp_turb)
{
/* Turbulent diffusivity (* rho) */
d_turb = C_MU_T(c,t) / Sct;

/* Effective diffusivity (* rho)*/
d_eff += d_turb;
}

return d_eff;
}


/* source term for the probability of 1st node, i.e. weight1 */

DEFINE_SOURCE(source_p1, c, t, dS, eqn)
{
  real s_p1;

  s_p1 = 0.;

  return s_p1;
}


/* source term for the weighted 1st node, i.e. abscissa1 of quadrature */

DEFINE_SOURCE(source_sxi_1, c, t, dS, eqn)
{
  real s_sxi1;
  
  s_sxi1 = C_R(c,t)*C_UDMI(c,t,3);

  return s_sxi1;
}


/* source term for the weighted 2nd node, i.e. abscissa2 of quadrature */

DEFINE_SOURCE(source_sxi_2, c, t, dS, eqn)
{
  real s_sxi2;

  s_sxi2 = -C_R(c,t)*C_UDMI(c,t,3);

  return s_sxi2;
}

/* source term of order 0 moment in node 1 */

DEFINE_SOURCE(source_m0_1, c, t, dS, eqn)
{
  real s_m0_1;

  s_m0_1 = C_R(c,t)*C_UDMI(c,t,19);

  return s_m0_1;
}


/* source term of order 1 moment in node 1 */

DEFINE_SOURCE(source_m1_1, c, t, dS, eqn)
{
  real s_m1_1;

  s_m1_1 = C_R(c,t)*C_UDMI(c,t,20);

  return s_m1_1;
}


/* source term of order 2 moment in node 1 */

DEFINE_SOURCE(source_m2_1, c, t, dS, eqn)
{
  real s_m2_1;

  s_m2_1 = C_R(c,t)*C_UDMI(c,t,21);

  return s_m2_1;
}


/* source term of order 3 moment in node 1 */

DEFINE_SOURCE(source_m3_1, c, t, dS, eqn)
{
  real s_m3_1;

  s_m3_1 = C_R(c,t)*C_UDMI(c,t,22);

  return s_m3_1;
}

/* source term of order 4 moment in node 1 */

DEFINE_SOURCE(source_m4_1, c, t, dS, eqn)
{
  real s_m4_1;

  s_m4_1 = C_R(c,t)*C_UDMI(c,t,52);

  return s_m4_1;
}

/* source term of order 5 moment in node 1 */

DEFINE_SOURCE(source_m5_1, c, t, dS, eqn)
{
  real s_m5_1;

  s_m5_1 = C_R(c,t)*C_UDMI(c,t,53);

  return s_m5_1;
}


/* source term of order 0 moment in node 2 */

DEFINE_SOURCE(source_m0_2, c, t, dS, eqn)
{
  real s_m0_2;

  s_m0_2 = C_R(c,t)*C_UDMI(c,t,23);

  return s_m0_2;
}


/* source term of order 1 moment in node 2 */

DEFINE_SOURCE(source_m1_2, c, t, dS, eqn)
{
  real s_m1_2;

  s_m1_2 = C_R(c,t)*C_UDMI(c,t,24);

  return s_m1_2;
}


/* source term of order 2 moment in node 2 */

DEFINE_SOURCE(source_m2_2, c, t, dS, eqn)
{
  real s_m2_2;

  s_m2_2 = C_R(c,t)*C_UDMI(c,t,25);

  return s_m2_2;
}


/* source term of order 3 moment in node 2*/

DEFINE_SOURCE(source_m3_2, c, t, dS, eqn)
{
  real s_m3_2;

  s_m3_2 = C_R(c,t)*C_UDMI(c,t,26);

  return s_m3_2;
}

/* source term of order 4 moment in node 2*/

DEFINE_SOURCE(source_m4_2, c, t, dS, eqn)
{
  real s_m4_2;

  s_m4_2 = C_R(c,t)*C_UDMI(c,t,54);

  return s_m4_2;
}

/* source term of order 5 moment in node 2*/

DEFINE_SOURCE(source_m5_2, c, t, dS, eqn)
{
  real s_m5_2;

  s_m5_2 = C_R(c,t)*C_UDMI(c,t,55);

  return s_m5_2;
}


DEFINE_ADJUST(adjust,domain)
{
  cell_t c;
  Thread *t;
  
  int i,k,j;
  real gamma;
  real p1,p2;
  real sxi_1, sxi_2, xi_1, xi_2, sxi_1m, sxi_2m, xi_1v, xi_2v;
  real Rel, C_phi;
  real c_1, c_2, c_e1, c_e2;
  real x_w1, x_w2;

  /*real dmy;*/

  real mom[2*node+1], mkept[2*node+1];
  real mom_p[2*node+3], momsec_p[2*node+3];
  real momsec[2*node+1], mseckept[2*node+1];

  real a_1[node+1],w_1[node+1];
  real a_2[node+1],w_2[node+1];
  real z_1[node+1][node+1], z_2[node+1][node+1];
  real d_1[node+1],d_2[node+1],subdiag_1[node+1],subdiag_2[node+1];

  real birth_1[2*node+3],death_1[2*node+3];
  real birth_2[2*node+3],death_2[2*node+3];
  real brown /*, cost_brown*/;
  real kturb=0;
  real m4_1,m4_2;
  real nu_1,nu_2, k_1, k_2;
  real mu_1,mu_2;
  real Rgn1,Rgn1_1,Rgn1_2;
/*
FILE * fp;
FILE * fp2;*/
  /* real v_m, d_m;
  real d_m2, v_m2, D_2, mu_1;
  real d_m1, v_m1, D_1, mu_2;*/
  /* real R1_1,R1_2; */ /*Dimension of a single chain in each environment*/
  /*real Sx1,Sy1,Sz1,Sx2,Sy2,Sz2,Rx,Ry,Rz,RgA,RgW; */

  void PDalg(real mom[]);
  void tqli(real d[], real e[], int n, real z[node+1][node+1]);

  real PBEexch[2*node+1];
  real kin_visc;
  
  int nb_uds = 14;

/*
fp = fopen ("Nucleation", "w");
fp2 = fopen ("Sources", "w");*/

  if (!Data_Valid_P())
    return;

    
  thread_loop_c(t,domain)
    {
      begin_c_loop(c,t)
    {
     
      for (i=0; i<=nb_uds; ++i)
        {
          if (C_UDSI(c,t,i) < 0.) {C_UDSI(c,t,i) = 0.;}
        }
    
      for (i=0; i<=6; ++i)
        {
          if (C_UDSI(c,t,i) > 1.-eps) {C_UDSI(c,t,i) = 1.;}
        }

      p1     = C_UDSI(c,t,0);

	if (p1<eps){p1=0.0;}

	p2 = 1. - p1;

	if (C_UDSI(c,t,1) < eps2){
		sxi_1m = 0.0;
		}
	else {
 		sxi_1m = C_UDSI(c,t,1);
		}

	if (C_UDSI(c,t,2) < eps2){
		sxi_2m = 0.0;
		}
	else {
		sxi_2m = C_UDSI(c,t,2);
		}
/*da articolo CFD modelling and scale-up of Confined Impinging Jet Reactors di Gavi (2007)*/
      kin_visc = C_MU_L(c,t)/C_R(c,t);
      C_UDMI(c,t,2) = kin_visc;

      Rel = C_K(c,t)/pow(kin_visc*C_D(c,t),0.5);
      C_UDMI(c,t,36) = Rel;

        
        if (Rel>0.2){
                       C_phi = 0.4093 + 0.6015*log10(Rel) +
                       0.5851*pow(log10(Rel),2) + 0.09472*pow(log10(Rel),3) -
                       0.3903*pow(log10(Rel),4) + 0.1461*pow(log10(Rel),5) -
                       0.01604*pow(log10(Rel),6);
                       }
        else
        { C_phi = 2.0;}

      C_UDMI(c,t,37) = C_phi;
          
        
      if (C_K(c,t) > 1.0e-3) {
        /*da articolo Di pasquale 2012*/
        gamma = C_phi/2.0*C_D(c,t)/C_K(c,t);}
      else {gamma=0.;}
     
              
      /* mixture fraction in node 1 */
      if (p1 > eps) {
        xi_1 = sxi_1m/p1;
      }
      else {
        xi_1 = 0.;
      }
      C_UDSI(c,t,3) = xi_1;
      
      /* mixture fraction in node 2 */
      if (p2 > eps) {
        xi_2 = sxi_2m/p2;
      }
      else {
        xi_2 = 0.;
      }
      C_UDSI(c,t,4) = xi_2;
  
      if (xi_1 < eps) {
		  xi_1v = 0.;
	  }
	  else {
      xi_1v = pow((1+((1/xi_1-1)*(rhoA/rhoW))),-1);/*frazione in volume di A nell'ambiente 1*/
	  }
	  if (xi_2 < eps) {
		  xi_2v = 0.;
	  }
	  else {
      xi_2v = pow((1+((1/xi_2-1)*(rhoA/rhoW))),-1);/*frazione in volume di A nell'ambiente 2*/
	  }
	  sxi_1 = p1 * xi_1v;
	  sxi_2 = p2 * xi_2v;      

      /* Source term for the weighted abcissas of the mixture fraction */

      /*C_UDMI(c,t,3) = gamma*(p1*sxi_2 - p2*sxi_1);*/
      C_UDMI(c,t,3) = gamma*(p1*sxi_2m - p2*sxi_1m);


      /* mean mixture fraction*/
      
      C_UDMI(c,t,0) = sxi_1+sxi_2; /* corrisponde alla frazione in volume di A media*/

      /* mixture fraction variance*/

      if (p1*p2 > eps*eps) {
         /*C_UDMI(c,t,1) = pow(p2*sxi_1-p1*sxi_2,2.)/(p1*p2);*/
         C_UDMI(c,t,1) = pow(p2*sxi_1m-p1*sxi_2m,2.)/(p1*p2);  /*eq.21 articolo Di Pasquale 2012*/
      }
      else {
           C_UDMI(c,t,1) = 0.;
      }

/*------------------------------------------------------------------------------*/
/*                  Definition of the operative variables                       */
/*------------------------------------------------------------------------------*/


/*water molar fraction*/

    x_w1= (1-xi_1v)/VolW/(xi_1v/VolA+(1-xi_1v)/VolW);
    x_w2= (1-xi_2v)/VolW/(xi_2v/VolA+(1-xi_2v)/VolW);/*come in definizione di viscosità utilizzando la frazione di volume di A nell'ambiente 1 e 2*/

    if (x_w1 < 0.) {x_w1 = 0.;}
    if (x_w2 < 0.) {x_w2 = 0.;}
    if (x_w1 > 1.-eps) {x_w1 = 1.;}
    if (x_w2 > 1.-eps) {x_w2 = 1.;}

    C_UDMI(c,t,41) = x_w1;
    C_UDMI(c,t,42) = x_w2;
    C_UDMI(c,t,50) = p1*x_w1 + p2*x_w2;

/*Viscosity*/
      mu_1 = exp((1.-x_w1)*log(muA)+x_w1*log(muW));
      mu_2 = exp((1.-x_w2)*log(muA)+x_w2*log(muW));

/* <<<-- Properties of the polymer  -->>>  */

/*equilibrium concentration--------------------------------------------------------------*/

    c_e1= 1200*exp(-14.533*x_w1)/M_w; /*kmol/m3 */
    c_e2= 1200*exp(-14.533*x_w2)/M_w; /*kmol/m3 errore nell'articolo Di pasquale 2012*/;

    C_UDMI(c,t,6) = c_e1;
    C_UDMI(c,t,7) = c_e2;
    C_UDMI(c,t,51) = p1*c_e1 +p2*c_e2;

/*PCL concentration*/ /*kmol/m3*/
/*eq.7 articolo alessio*/

    c_1= c_in/M_w*xi_1v /*- rho_pcl/M_w*pi/6.*C_UDSI(c,t,18)*/;  
    c_2= c_in/M_w*xi_2v /* - rho_pcl/M_w*pi/6.*C_UDSI(c,t,22)*/;  

    if (c_1 < 0.) {c_1= 0.;}
    if (c_2 < 0.) {c_2= 0.;}

    C_UDMI(c,t,4) = c_1;
    C_UDMI(c,t,5) = c_2;
    C_UDMI(c,t,10) = p1*c_1+p2*c_2;

/* Shape factors */

    /*dmy = (1.-x_w1);*/
    
    /*Sx1 = 0.514 + 1.450*dmy - 2.533*dmy*dmy + 1.366*dmy*dmy*dmy; 
    Sy1 = 0.293 - 0.710*dmy + 1.318*dmy*dmy - 0.748*dmy*dmy*dmy;
    Sz1 = 0.192 - 0.735*dmy + 1.206*dmy*dmy - 0.614*dmy*dmy*dmy;*/

    /*dmy = (1.-x_w2);*/

    /*Sx2 = 0.514 + 1.450*dmy - 2.533*dmy*dmy + 1.366*dmy*dmy*dmy; 
    Sy2 = 0.293 - 0.710*dmy + 1.318*dmy*dmy - 0.748*dmy*dmy*dmy;
    Sz2 = 0.192 - 0.735*dmy + 1.206*dmy*dmy - 0.614*dmy*dmy*dmy;*/

    nu_1= nuF(x_w1);  /*Flory's coefficient in this environment*/
    k_1 = kF(x_w1);   /*Constant of the Flory's law in this environment*/
    
    nu_2= nuF(x_w2);  /*Flory's coefficient in this environment*/
    k_2 = kF(x_w2);   /*Constant of the Flory's law in this environment*/

    /* v_m1 = 1.5836366 + 30.806046*x_w1; */  /*nm3*/
    /* v_m2 = 1.5836366 + 30.806046*x_w2;  */ /*nm3*/
    /* v_m1 = v_m1*1.E-27;    */    /* trasformazione da nm3 --> m3 */
    /* v_m2 = v_m2*1.E-27;    */    /* trasformazione da nm3 --> m3 */
    /* v_m  = p1*v_m1 + p2*v_m2; */
    /* d_m1 = pow(v_m1*3./4./pi,1./3.);*/
    /* d_m2 = pow(v_m2*3./4./pi,1./3.);*/
    /* d_m  = p1*d_m1 + p2*d_m2; */

      /*diffusion coefficient*/
      /*D_1 = k_B*T/(6.*pi*mu_1*d_m1);
      D_2 = k_B*T/(6.*pi*mu_2*d_m2);

      R1_1 = pow(k_1*pow(M_w,2.*nu_1),0.5); */ /* posso modificare direttamente questi */
      /* R1_2 = pow(k_2*pow(M_w,2.*nu_2),0.5); */ /* oppure riscrivermeli per fatti miei*/
      
    C_UDMI(c,t,32) = 0.0;
    C_UDMI(c,t,33) = 0.0;
    C_UDMI(c,t,8) = 0.0;

    C_UDMI(c,t,39) = 0.0;
    C_UDMI(c,t,40) = 0.0;
C_UDMI(c,t,52)=0.0; 
C_UDMI(c,t,53)=0.0;
C_UDMI(c,t,54)=0.0;
C_UDMI(c,t,55)=0.0;
C_UDMI(c,t,56)=0.0;
C_UDMI(c,t,57)=0.0;
C_UDMI(c,t,58)=0.0;

      C_UDMI(c,t,43) = 0.0;
      C_UDMI(c,t,44) = 0.0;
      C_UDMI(c,t,9)  = 0.0;
      C_UDMI(c,t,34) = 0.0;
      C_UDMI(c,t,35) = 0.0;


/* --------------- Quadrature method of moments ---------------*/

	for (i=1;i<=2*node;i++) {
		if (p1 > eps) {
		      C_UDSI(c,t,14+i)= C_UDSI(c,t,7-1+i)/p1;
		}
		else {
		C_UDSI(c,t,14+i)= 0.;
		}
		if (p2 > eps) {
		      C_UDSI(c,t,14+2*node+i)= C_UDSI(c,t,7-1+2*node+i)/p2;
		}
		else {
		      C_UDSI(c,t,14+2*node+i)= 0.;
		}
		mom_p[i] = C_UDSI(c,t,14+i);
		momsec_p[i] = C_UDSI(c,t,14+2*node+i);
	}

	for (i=1;i<2*node+1;i++){
            mom[i] = mkept[i] = mom_p[i];
      	momsec[i] = mseckept[i] = momsec_p[i];
      	}

	for (j=1;j<=node;j++) {
		for (i=1;i<=node;i++) {
			if (i==j) {z_1[i][j] = z_2[i][j] = 1.;}
			else {
			z_1[i][j] = z_2[i][j] = 0.;
			}
        	}
      	}  	
      	
	
	if (mom[1] != 0. && mom[2] > eps) /* condizione per stabilità */ {
		PDalg(mom);
		for (j=1;j<=node;j++){
        		d_1[j] = mom[j];
        		subdiag_1[j] = mom[node+j];
          	}
		tqli(d_1,subdiag_1,node,z_1);
        }
        else {
         for (j=1;j<=node;j++){
         	d_1[j] = 0.;
        	z_1[1][j] = 0.;
         }
        }


	if (momsec[1] != 0. && momsec[2] > eps) /* condizione per stabilità */  {
		PDalg(momsec);
		for (j=1;j<=node;j++){
        		d_2[j] = momsec[j];
        		subdiag_2[j] = momsec[node+j];
          	}
		tqli(d_2,subdiag_2,node,z_2);
        }
        else {
         for (j=1;j<=node;j++){
         	d_2[j] = 0.;
        	z_2[1][j] = 0.;
         }
        }

	
      w_1[1] = z_1[1][1]*z_1[1][1]*mkept[1];
      w_1[2] = z_1[1][2]*z_1[1][2]*mkept[1];

      w_2[1] = z_2[1][1]*z_2[1][1]*mseckept[1];
      w_2[2] = z_2[1][2]*z_2[1][2]*mseckept[1];
      
      
      if (fabs(w_1[2]) > 0.) {
        if (fabs(w_1[1])/fabs(w_1[2]) > 1.0e3) d_1[2] = d_1[1];
      }
      if (fabs(w_2[2]) > 0.) {
        if (fabs(w_2[1])/fabs(w_2[2]) > 1.0e3) d_2[2] = d_2[1];
      }

      C_UDMI(c,t,11) = a_1[1] = fabs(d_1[1]);
      C_UDMI(c,t,12) = a_1[2] = fabs(d_1[2]);
      C_UDMI(c,t,13) = w_1[1] = z_1[1][1]*z_1[1][1]*mkept[1];
      C_UDMI(c,t,14) = w_1[2] = z_1[1][2]*z_1[1][2]*mkept[1];

      
      C_UDMI(c,t,15) = a_2[1] = fabs(d_2[1]);
      C_UDMI(c,t,16) = a_2[2] = fabs(d_2[2]);
      C_UDMI(c,t,17) = w_2[1] = z_2[1][1]*z_2[1][1]*mseckept[1];
      C_UDMI(c,t,18) = w_2[2] = z_2[1][2]*z_2[1][2]*mseckept[1];

            
	C_UDMI(c,t,45) = 0.0;
	C_UDMI(c,t,46) = 0.0;
	C_UDMI(c,t,47) = 0.0;
	C_UDMI(c,t,48) = 0.0;

/*----------------------source term of moments-------------------------*/

    /* brown = NA*0.01*8.*k_B*T/3.;  */    
    /* cost_brown = brown/C_MU_L(c,t); */ /*2*k_B*T/3/mu*/

    for (k=1;k <= 2*node;k++) {
        PBEexch[k] = gamma*(p1*C_UDSI(c,t,6+2*node+k)-p2*C_UDSI(c,t,6+k));
      }

    for (k=1;k < 2*node+3;k++) {
        birth_1[k] = 0.;
        death_1[k] = 0.;
        birth_2[k] = 0.;
        death_2[k] = 0.;    
 

            for (i=1;i <= node;i++) {
                for (j=1;j <= node;j++) {
                    if (c_1>c_e1 && a_1[i]*a_1[j] != 0.) {
                            brown = NA*(2.*k_B*T/(mu_1*3.))*pow(sqrt(k_1*(pow(a_1[i]*M_w,2.*nu_1)))+sqrt(k_1*(pow(a_1[j]*M_w,2.*nu_1))),2)/(sqrt(k_1*(pow(a_1[i]*M_w,2.*nu_1)))*sqrt(k_1*(pow(a_1[j]*M_w,2.*nu_1))));
                            kturb=NA*1.2944*pow(C_D(c,t)/kin_visc,0.5)*pow(0.000000001*(sqrt(k_1*(pow(a_1[i]*M_w,2.*nu_1)))+sqrt(k_1*(pow(a_1[j]*M_w,2.*nu_1)))),3); /* coefficiente serve per passare da nm a m!! */	
                            C_UDMI(c,t,71+(j-1)*2+(i-1)) = kturb/brown; 
      		            birth_1[k] = birth_1[k] + 0.5*w_1[i]*w_1[j]*(brown+kturb)
                             *pow(a_1[i]+a_1[j],(k-1.));
/**pow(sqrt(k_1*(pow(a_1[i]*M_w,2.*nu_1)+sqrt(k_1*(pow(a_1[j]*M_w,2.*nu_1)),2.)/(sqrt(k_1*(pow(a_1[i]*M_w,2.*nu_1)*sqrt(k_1*(pow(a_1[j]*M_w,2.*nu_1));*/
                            death_1[k] = death_1[k] + w_1[i]*pow(a_1[i],k-1)*w_1[j]*(brown+kturb);
                       }
                   else {
                        birth_1[k] = birth_1[k];
                        death_1[k] = death_1[k];
                     }

                   if (c_2>c_e2 && a_2[i]*a_2[j] != 0.) {
                        brown = NA*(2.*k_B*T/(mu_2*3.))*pow(sqrt(k_2*(pow(a_2[i]*M_w,2.*nu_2)))+sqrt(k_2*(pow(a_2[j]*M_w,2.*nu_2))),2)/(sqrt(k_2*(pow(a_2[i]*M_w,2.*nu_2)))*sqrt(k_2*(pow(a_2[j]*M_w,2.*nu_2))));
                        kturb=NA*1.2944*pow(C_D(c,t)/kin_visc,0.5)*pow(0.000000001*(sqrt(k_2*(pow(a_2[i]*M_w,2.*nu_2)))+sqrt(k_2*(pow(a_2[j]*M_w,2.*nu_2)))),3);
                        C_UDMI(c,t,75+(j-1)*2+(i-1)) = kturb/brown; 
                        birth_2[k] = birth_2[k] + 0.5*w_2[i]*w_2[j]
                              *pow(a_2[i]+a_2[j],(k-1.))*(brown+kturb);

/* *pow(sqrt(k_2*(pow(a_2[i]*M_w,2.*nu_2)+sqrt(k_2*(pow(a_2[j]*M_w,2.*nu_2)),2.)/(sqrt(k_2*(pow(a_2[i]*M_w,2.*nu_2)*sqrt(k_2*(pow(a_2[j]*M_w,2.*nu_2));*/
                        death_2[k] = death_2[k] + w_2[i]*pow(a_2[i ],k-1)*w_2[j]*(brown+kturb);
                        }
                   else {
                        birth_2[k] = birth_2[k];
                        death_2[k] = death_2[k];
                    }

/*
            for (i=1;i <= node;i++) {
                for (j=1;j <= node;j++) {
                    if (a_1[i]*a_1[j] != 0.) {		
      		            birth_1[k] = birth_1[k] + 0.5*w_1[i]*w_1[j]
                            *pow(pow(a_1[i],3)+pow(a_1[j],3),(k-1.)/3.)*(cost_brown*(pow(a_1[i]+a_1[j], 2)/(a_1[i]*a_1[j])));
                        death_1[k] = death_1[k] + w_1[i]*pow(a_1[i],k-1)*w_1[j]*(cost_brown*(pow(a_1[i]+a_1[j], 2)/(a_1[i]*a_1[j])));
                       }
                   else {
                        birth_1[k] = birth_1[k];
                        death_1[k] = death_1[k];
                     }

                   if (a_2[i]*a_2[j] != 0.) {
                        birth_2[k] = birth_2[k] + 0.5*w_2[i]*w_2[j]
                              *pow(pow(a_2[i],3)+pow(a_2[j],3),(k-1.)/3.)*(cost_brown*(pow(a_2[i]+a_2[j], 2)/(a_2[i]*a_2[j])));
                        death_2[k] = death_2[k] + w_2[i]*pow(a_2[i],k-1)*w_2[j]*(cost_brown*(pow(a_2[i]+a_2[j], 2)/(a_2[i]*a_2[j])));
                        }
                   else {
                        birth_2[k] = birth_2[k];
                        death_2[k] = death_2[k];
                    }
*/
                  	      }
          		      }
       		     }
	        



      if (p1 > eps && p1 < 1.-eps ) {
        C_UDMI(c,t,19) = PBEexch[1] + p1*(birth_1[1] - death_1[1]);
        C_UDMI(c,t,20) = PBEexch[2] + p1*(birth_1[2] - death_1[2]);
        C_UDMI(c,t,21) = PBEexch[3] + p1*(birth_1[3] - death_1[3]);
        C_UDMI(c,t,22) = PBEexch[4] + p1*(birth_1[4] - death_1[4]);
      }
      else  {
        C_UDMI(c,t,19) = 0.;
        C_UDMI(c,t,20) = 0.;
        C_UDMI(c,t,21) = 0.;
        C_UDMI(c,t,22) = 0.;
      }

      if (p2 > eps && p2 < 1.-eps) {
        C_UDMI(c,t,23) = -PBEexch[1] + p2*(birth_2[1] - death_2[1]);
        C_UDMI(c,t,24) = -PBEexch[2] + p2*(birth_2[2] - death_2[2]);
        C_UDMI(c,t,25) = -PBEexch[3] + p2*(birth_2[3] - death_2[3]);
        C_UDMI(c,t,26) = -PBEexch[4] + p2*(birth_2[4] - death_2[4]);
      }
      else  {
        C_UDMI(c,t,23) = 0.;
        C_UDMI(c,t,24) = 0.;
        C_UDMI(c,t,25) = 0.;
        C_UDMI(c,t,26) = 0.;	
      }

      /* moment of order 0 */
      C_UDMI(c,t,27) = C_UDSI(c,t,7) + C_UDSI(c,t,11);
      /* moment of order 1 */
      C_UDMI(c,t,28) = C_UDSI(c,t,8) + C_UDSI(c,t,12);
      /* moment of order 2 */
      C_UDMI(c,t,29) = C_UDSI(c,t,9) + C_UDSI(c,t,13);
      /* moment of order 3 */
      C_UDMI(c,t,30) = C_UDSI(c,t,10) + C_UDSI(c,t,14);

    /*initial supersaturation articolo alessio*/
    C_UDMI(c,t,79) = (C_UDMI(c,t,28)*0.001)/C_UDMI(c,t,51);/*fattore serve per avere rapporto tra kmol/kmol*/

     /* moment of order 4 */

      m4_1 = w_1[1]*pow(a_1[1],4) + w_1[2]*pow(a_1[2],4);

      m4_2 = w_2[1]*pow(a_2[1],4) + w_2[2]*pow(a_2[2],4);

      C_UDMI(c,t,49)  = p1*m4_1 + p2*m4_2;

      C_UDMI(c,t,45) = sqrt(k_1*(pow(a_1[1]*M_w,2.*nu_1)));
      C_UDMI(c,t,46) = sqrt(k_1*(pow(a_1[2]*M_w,2.*nu_1)));
      C_UDMI(c,t,47) = sqrt(k_2*(pow(a_2[1]*M_w,2.*nu_2)));
      C_UDMI(c,t,48) = sqrt(k_2*(pow(a_2[2]*M_w,2.*nu_2)));

     /* Volume-averaged diameter eq.13 articolo alessio */
      if (C_UDMI(c,t,30) > 0.) {

        C_UDMI(c,t,31) = p1*(w_1[1]*C_UDMI(c,t,45)+w_1[2]*C_UDMI(c,t,46))/(w_1[1]+w_1[2])+p2*(w_2[1]*C_UDMI(c,t,47)+w_2[2]*C_UDMI(c,t,48))/(w_2[1]+w_2[2]);
      }
      else {C_UDMI(c,t,31) = 0.;};

 /* Stoke−Einstein diffusion coefficient [m2/s] */

     C_UDMI(c,t,80) = k_B*T/(6.*pi*C_MU_L(c,t)*1.0e-9*C_UDMI(c,t,31)); 

    /* Calcolo S-E diffusion coefficient per singola molecola PCL e kernel browniano di aggregazione tra due singole molecole */

    if (p1 > eps && p1 < 1.-eps ) {
         
	Rgn1_1 = sqrt(k_1*(pow(M_w,2.*nu_1))); }
     else {
	Rgn1_1 = 0.; }		
	        
     if (p2 > eps && p2 < 1.-eps) {	
	
	Rgn1_2 = sqrt(k_2*(pow(M_w,2.*nu_2))); }
     else {
        Rgn1_2 = 0.; }

     if (Rgn1_1*Rgn1_2!=0.) {
	 Rgn1 = p1*Rgn1_1+p2*Rgn1_2;
        C_UDMI(c,t,81) = k_B*T/(6.*pi*C_MU_L(c,t)*1.0e-9*Rgn1); /* S-E diffusion coefficient per singola molecola PCL [m2/s] */
        C_UDMI(c,t,82) = 4.*pi*(C_UDMI(c,t,81)+C_UDMI(c,t,81))*1.0e-9*(Rgn1+Rgn1); /* kernel browniano di aggregazione [m3/s] */
	C_UDMI(c,t,83) = C_D(c,t)/kin_visc; /* termine presente nel kernel turbolento [s-2] */
	C_UDMI(c,t,84) = turb*sqrt(C_D(c,t)/kin_visc)*pow(1.0e-9*(Rgn1+Rgn1),3); /* kernel turbolento di aggregazione tra due singole molecole [m3/s] */
	}
     
     else {
        Rgn1 = 0.;}

      /*Average Peclét Number, <Pe>, [-] vale nell'ipotesi che le dimensioni delle particelle siano inferiori alla scala di Kolmogorov.
                                          E siamo proprio in questa situazione (anche perchè le particelle hanno inerzia trascurabile: St<1)*/

      
      C_UDMI(c,t,32)=(3.*pi*mu_1*pow((w_1[1]*C_UDMI(c,t,45)+w_1[2]*C_UDMI(c,t,46))/(w_1[1]+w_1[2])*0.000000001,3.)*pow(C_D(c,t)/kin_visc,0.5))/(k_B*T);
      C_UDMI(c,t,33)=(3.*pi*mu_2*pow((w_2[1]*C_UDMI(c,t,47)+w_2[2]*C_UDMI(c,t,48))/(w_2[1]+w_2[2])*0.000000001,3.)*pow(C_D(c,t)/kin_visc,0.5))/(k_B*T);
      C_UDMI(c,t,34)=p1*C_UDMI(c,t,32)+p2*C_UDMI(c,t,33);
      

      C_UDMI(c,t,70)=pow(pow(kin_visc,3.)/C_D(c,t),0.25); /* scala di Kolmogorov*/


    }

    end_c_loop (c,t)
    }

	/*fclose(fp);*/

    }

/* THE END  */

