#include <stdio.h>
#include <math.h>
#include <ascii.h>
#include <cdisort.h>
#include <c_tzs.h>
#include "locate.h"
#if HAVE_LIBGSL 
  #include <gsl/gsl_sf_expint.h>
#endif

// include DISORT.MXD

/*  nlyr   = number of atmospheric layers                     */
/*  nlev_c = number of atmospheric levels on common z grid    */
/*  zd_c   = atmospheric levels on common z grid              */
/*  nzout  = number of output levels = number of black clouds */
/*  zout   = vertical level heights in km above surface       */

int c_tzs( int nlyr, float *dtauc, int nlev_c, float *zd_c,
	   int nzout, float *zout, 
	   float *ssalb, float *temper, 
	   float wvnmlo, float wvnmhi, 
	   int usrtau, int ntau, float *utau, 
	   int usrang, int numu, float *umu, int nphi, float *phi, 
	   float albedo, double btemp, float ttemp, 
	   float temis, int planck, 
	   int *prnt, char *header,                //  int prndis[7]; char header[127]
	   float *rfldir, float *rfldn, float *flup,
	   float *dfdt, float *uavg,
	   float ***uu, int quiet)
{
  /*    .. Internal variables  */
  int         lc, j, iu, lu;
  double      *tauc=NULL;
  double      tau1=0, tau2=0, etau1=0, etau2=0, etaumu=0, ei1=0, ei2=0, *rad=NULL;
  double      *radem=NULL, *pkag=NULL, *etau=NULL, rademtot=0.0;
  double      radrefl, radrefl1, radrefl2, radrefl11, radrefl22;
  double      *b0=NULL, *b1=NULL;
  double      *coeff1=NULL, *coeff2=NULL;
  double      sqtau1=0.0, sqtau2=0.0, mu=0.0, cutau1=0.0, tauetau1;
  int         *nlayers=NULL, maxnlyr=0;
  float       *alb=NULL, *bplanck=NULL;

  /* nlyr     = number of layers */
  /* nlyr + 1 = number of levels */

  /* initialise */

  if ( (rad     = (double *) calloc (ntau,    sizeof(double))) == NULL )
    return ASCII_NO_MEMORY;
  if ( (nlayers = (int *)    calloc (ntau,  sizeof(int))) == NULL )
    return ASCII_NO_MEMORY;

  if ( (alb      = (float *) calloc (ntau,  sizeof(float))) == NULL )
    return ASCII_NO_MEMORY;
  if ( (bplanck  = (float *) calloc (ntau,  sizeof(float))) == NULL )
    return ASCII_NO_MEMORY;

  if ( prnt[0] && strlen( header )!=0 ){
    fprintf(stderr,"**************\n TZS: %s\n*****************\n",header);
  }

  /* Initialise everything to 0. */
  for (iu=0;iu<numu;iu++){
    for (lu=0;lu<nzout;lu++){
      for (j=0;j<nphi;j++){
	uu[j][lu][iu]=0.;
      }
    }
  }
  for (lu=0;lu<nzout;lu++){
    dfdt[lu]   = 0.;
    flup[lu]   = 0.;
    rfldir[lu] = 0.;
    rfldn[lu]  = 0.;
    uavg[lu]   = 0.;
  }

  /*    Layers: 0...NLYR-1 (# layers = NLYR)                                                           */
  /*    Levels: 0...NLYR   (# levels = NLYR+1)                                                         */
  /*    LC = 0...NYR-1                                                                                 */
  /*    LC = 0 : TOA                                                                                   */
  /*    -------------------------------------------------- Level LC, TAUC(LC), TEMPER(LC)              */
  /*                                                                                                   */
  /*                  Layer LC                                                                         */
  /*                                                                                                   */
  /*    -------------------------------------------------- Level LC+1  , TAUC(LC+1)  , TEMPER(LC+1)    */
    
  /*    Planck radiation from surface with emissivity=1-ALBEDO and temperature BTEMP */
  /*     BPLANCK = (1-ALBEDO)*PLKAVG2( WVNMLO, WVNMHI, BTEMP )  */
  
  /*    Max number of layers counted from TOA used for all UTAU  */
  maxnlyr=-1;
        
  /*    Loop over user defined tau = black cloud top heights:  */
  /*    determine layer number corresponding to user tau       */
  /*    Rem.: NTAU = NZOUT checked in CHEKIN_TZS               */

  if (ntau != nzout ){
    errmsg_tzs("TZS--nzout != ntau",TZS_ERROR);
  }

  for (lu=0;lu<ntau;lu++){
    
    /* zd_c[0] = TOA                                                         */
    /* zd_c[nlyr+1] = SUR                                                    */
    /* locate gives back lc such that zd_c(lc) >= zout(lu) > zd_lc[lc+1]     */
    lc = flocate(zd_c,nlev_c,zout[lu]);

    //    fprintf(stderr,"lc = %d  zd_c[lc] = %f  zd_c[lc+1] = %f  zout[lu] = %f\n",lc,zd_c[lc],zd_c[lc+1],zout[lu]);
     
    if ( lc==-1 || lc==nlev_c ){
      errmsg_tzs("TZS--Zout outside atmospheric levels",TZS_ERROR);
    }
    
    /* Check real presence of this level in tauc ( max diff = 1mm )                         */
    /* Determine layer above given user level and store layer index w.r.t. zd_c in nlayers  */
    /* nlayers contains the layer index whose bottom is at zout[lu]                         */

    if ( fabs(zd_c[lc+1]-zout[lu])<=1.e-6 ) {

      /*          check lower boundary           */

      /*          nlevel = nlayers[lu] = LC+1-1  */
      /*          nlayer = nlevel - 1            */
      nlayers[lu]=lc;
      
    } else if ( fabs(zd_c[lc]-zout[lu])<=1.e-6 ) {
      
      /*          check upper boundary           */

      /*          nlevel = nlayers[lu] = LC      */
      /*          nlayer = nlevel - 1            */
      nlayers[lu]=(int) fmax(lc-1,0);
      /* nlayers could be -1 when black cloud is at TOA */
      
    } else {

      /* zout_interpolate */

      //      fprintf(stderr,"%d4 %4d %9.3f %9.3f %9.3f %9.3f %9.3f\n", 
      //	      lc,nlev_c,zd_c[lc+1],zd_c[lc],zout[lu],fabs(zd_c[lc+1]-zout[lu]),fabs(zd_c[lc]-zout[lu]));
      errmsg_tzs("TZS--User height is not contained in atmospheric levels: use zout_interpolate!",TZS_ERROR);                        

    }

//    fprintf(stderr,"nlev_c = %4d   zout[lu] = %9.3f  zd_c[nlayers[lu]] = %9.3f   zd_c[nlayers[lu]+1] = %9.3f   lu = %4d\n", 
//	    nlev_c,zout[lu],zd_c[nlayers[lu]],zd_c[nlayers[lu]+1],lu);
//
//    fprintf(stderr,"lu = %4d   maxnlyr = %4d   nlayers[lu] = %4d\n",lu,maxnlyr,nlayers[lu]);

    maxnlyr=(int) fmax(maxnlyr,nlayers[lu]);

  }

  /* Above we defined layer indices, now we add +1 because we want the number of layers */
  maxnlyr += 1;

  //  fprintf(stderr,"%s %4d\n","Max # layers = ",maxnlyr);

  if (maxnlyr<0) {
    errmsg_tzs("TZS--wrong zout",TZS_ERROR);
  }

  if ( (radem   = (double *) calloc (maxnlyr,    sizeof(double))) == NULL )
    return ASCII_NO_MEMORY;
  if ( (tauc    = (double *) calloc (maxnlyr+1,  sizeof(double))) == NULL )
    return ASCII_NO_MEMORY;
  if ( (pkag    = (double *) calloc (maxnlyr+1,  sizeof(double))) == NULL )
    return ASCII_NO_MEMORY;
  if ( (etau    = (double *) calloc (maxnlyr+1,  sizeof(double))) == NULL )
    return ASCII_NO_MEMORY;
	  
  if ( (b0      = (double *) calloc (maxnlyr,  sizeof(double))) == NULL )
    return ASCII_NO_MEMORY;
  if ( (b1      = (double *) calloc (maxnlyr,  sizeof(double))) == NULL )
    return ASCII_NO_MEMORY;
  if ( (coeff1  = (double *) calloc (maxnlyr,  sizeof(double))) == NULL )
    return ASCII_NO_MEMORY;
  if ( (coeff2  = (double *) calloc (maxnlyr,  sizeof(double))) == NULL )
    return ASCII_NO_MEMORY;

  /*    calculate cumulative optical depth tauc[lc] up to level lc and      */
  /*    Planck radiation for every level (not layer!) from 0 to MAXNLYR     */
  /*    Remark: layer 1 is the first layer from above=TOA??                   */
  /*            level 0 is the first level from above=TOA, level 0 = TOA    */

  /* lc = level count */
  lc       = 0;
  tauc[lc] = 0.0;
  pkag[lc] = (double) c_planck_func1( wvnmlo, wvnmhi, temper[lc] );

  //  fprintf(stderr,"lc = %d  maxnlyr = %d  pkag[lc] = %f  tauc[lc] = %13.6e dtauc[lc-1] = %13.6e temper[lc] = %9.3f\n",lc,maxnlyr,pkag[lc],tauc[lc],dtauc[lc-1],temper[lc]);

  for (lc=1;lc<=maxnlyr;lc++){
    /* We work here with atmospheric optical thickness as it is provided by the user. */
    /* If you want to use aborption optical thickness the you have to define it in    */
    /* the uvspec input file.                                                         */
    tauc[lc] = tauc[lc-1]+ ((double) dtauc[lc-1]);

    pkag[lc] = (double) c_planck_func1( wvnmlo, wvnmhi, temper[lc] );

    //    fprintf(stderr,"lc = %d  maxnlyr = %d  pkag[lc] = %f  tauc[lc] = %13.6e dtauc[lc-1] = %13.6e temper[lc] = %9.3f\n",lc,maxnlyr,pkag[lc],tauc[lc],dtauc[lc-1],temper[lc]);

  }
  
  /*    Print input information  */
  if ( prnt[0] ){ 
  
    prtinp_tzs(nlyr, maxnlyr, dtauc, ssalb, temper, 
	       wvnmlo, wvnmhi, ntau, utau, numu, umu,
	       nphi, phi, albedo, btemp, ttemp, temis,
	       tauc, nlev_c, zd_c, nzout, zout );
  }
  
  /*    check input dimensions and variables */
  chekin_tzs ( nlyr, maxnlyr, dtauc, nlev_c, zd_c,
	       nzout, zout,
	       ssalb, temper, wvnmlo,
	       wvnmhi, usrtau, ntau, utau, usrang,
	       numu, umu, nphi, phi, 
	       albedo, btemp, ttemp,
	       temis, planck, 
	       tauc, quiet );

  for (lu=0;lu<ntau;lu++){
    //    fprintf(stderr,"lu = %4d  nlayers[lu] = %4d   nlyr = %4d\n",lu,nlayers[lu],nlyr);
    /* check layer bottom level */
    if ( nlayers[lu]+1==nlyr ) {
      /*          real surface albedo = 1 - emissivity                                         */
      alb[lu]     = albedo;
      /*          Planck radiation from surface with emissivity=1-ALBEDO and temperature BTEMP */
      bplanck[lu] = (1-alb[lu])*c_planck_func1( wvnmlo, wvnmhi, btemp );            
      //      fprintf(stderr,"%4d  %4d  %12.10f  %12.9f  %12.8f\n",lu, nlayers[lu]+1,alb[lu],bplanck[lu],btemp);  
    } else {
      /*          Real surface albedo = 1 - emissivity ==> emissivity=1 (black cloud)          */
      alb[lu]     = 0.0;
      /*          Planck radiation of bottom level with emissivity=1 (black cloud)             */
      /*          and temperature temper( nlayers[lu]+1 ) of bottom level                      */
      bplanck[lu] = (float) pkag[nlayers[lu]+1];
    }
  }

  if ( prnt[0] ) {
    fprintf( stderr, "\n%4d %s", ntau, " User optical depths :");
    for (lu=0;lu<ntau;lu++)
      fprintf(stderr,"%14.8f",utau[lu]);
    fprintf(stderr,"\n");

    fprintf( stderr, "\n%4d %s" , ntau, " User level numbers :");
    for (lu=0;lu<ntau;lu++)
      fprintf(stderr,"%4d",nlayers[lu]+1);
    fprintf(stderr,"\n");

    fprintf( stderr, "\n%4d %s", ntau, " User levels / km :");
    for (lu=0;lu<ntau;lu++)
      fprintf(stderr,"%14.8f",zout[lu]);
    fprintf(stderr,"\n");

    fprintf( stderr, "\n%4d %s", ntau, " User levels / km :");
    for (lu=0;lu<ntau;lu++)
      fprintf(stderr,"%14.8f",zd_c[nlayers[lu]+1]);
    fprintf(stderr,"\n");
  }
      
  /* loop from TOA to lowest layer/level needed for calculations */
  for (lc=0;lc<maxnlyr;lc++){
    /*          assume linear variation of Planck radiation from one level                     */
    /*          to the next one: Planck=b0+tauc*b1                                             */
    b1[lc]=0.;

    /* CE: this threshold was commented out, but it seems that it is needed because otherwise we get 
     negative emissions in some layers */
    if (dtauc[lc]>1.e-4)
      b1[lc]=(pkag[lc+1]-pkag[lc])/dtauc[lc];
    /* In case of thin atmosphere but high temperature variation from one level to the next, dtauc is very small  */
    /* but Planck varies sensibly, so we have to take care. The threshold 1.e10 is not crucial because thin       */
    /* atmospheric layers emit few radiation.                                                                     */
    if (fabs(b1[lc])>1.e10){
      /* set b0 to average Planck emission, also not crucial */
      b0[lc]=0.5*(pkag[lc+1]+pkag[lc]);
      b1[lc]=0.0;	  
    } else {
      b0[lc]=pkag[lc]-tauc[lc]*b1[lc];
    }
    
    /* Compute coefficients like this even if coeff1[lc]=pkag[lc+1] and coeff2[lc]=pkag[lc]    */
    coeff1[lc] = b0[lc] + tauc[lc+1]*b1[lc];
    coeff2[lc] = b0[lc] + tauc[lc]*b1[lc];
    
    //    fprintf(stderr," LC %d Deltatau %e DIFFCOEFF1 %e COEFF1 %e pkag[lc+1] = %e\n",lc,dtauc[lc],coeff1[lc]-pkag[lc+1],coeff1[lc],pkag[lc+1]);
    //    fprintf(stderr," LC %d Deltatau %e DIFFCOEFF2 %e COEFF2 %e pkag[lc]   = %e\n",lc,dtauc[lc],coeff2[lc]-pkag[lc],coeff2[lc],pkag[lc]);
  }
  
  /*    Loop over UMUs    */
  for (iu=0;iu<numu;iu++){

    /*       Double precision UMU  */
    mu=(double) umu[iu];

//    /* loop from TOA to lowest layer/level needed for calculations */
//    for (lc=0;lc<=maxnlyr;lc++){
//      /*          Initialise ETAU = Attenuation factors from TOA to levels LC in direction MU  */
//      etau[lc]=exp(-tauc[lc]/mu);
//    }

    /*       Emission of atmosphere:   */
    /*       ***********************   */

    /*          Initialise ETAU = Attenuation factors from TOA to level LC in direction MU   */
    lc=0;
    etau[lc]=exp(-tauc[lc]/mu);

    /* loop from TOA to lowest layer/level needed for calculations */

    for (lc=0;lc<maxnlyr;lc++){

      /*          Initialise ETAU = Attenuation factors from TOA to level LC+1 in direction MU */
      etau[lc+1]=exp(-tauc[lc+1]/mu);

      /*          Attenuation factors from TOA to levels LC+1 resp. LC:   */
      /*           ETAU[lc+1]                                             */
      /*           ETAU[lc]                                               */

      /*          Radiation emitted from layer LC (between levels LC and LC+1) towards TOA   */
      radem[lc] = etau[lc]*coeff2[lc] - etau[lc+1]*coeff1[lc] + b1[lc]*mu*(etau[lc]-etau[lc+1]);

    }

    /*       Loop over user defined tau = black cloud top heights   */
    for (lu=0;lu<ntau;lu++){

      //      fprintf(stderr,"lu = %4d  %4d\n",lu,nlayers[lu]);

      /*          Initialise RAD to 0.              */
      rad[lu]   = 0.;
      /*          Initialise TAU1                   */
      /*          to optical thickness from TOA up to bottom of layer nlayers[lu] */
      /*           TAU1  = TAUC(NLYR)               */
      tau1      = tauc[nlayers[lu]+1]; 
      /*          Initialise ETAU1                  */
      etau1     = exp(-tau1);
      tauetau1  = tau1*etau1;
      /*          Initialise EI1                    */
      
#if HAVE_LIBGSL
      /* exponential integral does not converge for values larger about 600*/
      /* the integral for this case becomes 0*/
      if (tau1 > 500.0)
	ei1=0.0;
      else
	ei1       = gsl_sf_expint_Ei(-tau1);
      /* fprintf(stderr, "lu %d  tau1 %g ei1 %g\n",lu,  tau1, ei1);  */
#endif
#if !HAVE_LIBGSL 
  fprintf(stderr,"libRadtran was built without the gsl-library. Thus \n");
  fprintf(stderr,"tzs may not be used. Please install gsl on your system\n");
  fprintf(stderr,"and rebuild libRadtran.\n");
  exit(0);
#endif
  

      /*          Initialise SQTAU1                 */
      sqtau1    = tau1*tau1;
      /*          Initialise CUTAU1                 */
      cutau1    = tau1*sqtau1;
      /*          initialise RADREFL11              */
      radrefl11 = -2.0*etau1+tauetau1-tau1*tauetau1-cutau1*ei1;	 
            
      /*          Attenuation factor due to full atmospheric absorption in direction MU  */
      etaumu    = etau[nlayers[lu]+1];
      
      /* total atmospheric emission i.e. intergral from TOA to bottom */
      rademtot  = 0.0;

      /*          Contribution from surface + atmosphere              */
      for (lc=0;lc<=nlayers[lu];lc++){
	
	if (alb[lu]>0.0){

	  /*             Reflection from surface:   */
	  /*             ************************   */
	  tau2      = tau1;
	  etau2     = etau1;
	  ei2       = ei1;
	  sqtau2    = sqtau1;
	  radrefl22 = radrefl11;
	  radrefl2  = etau2*coeff2[lc]*(1.-tau2)-coeff2[lc]*sqtau2*ei2;
	  
	  if (lc==nlayers[lu]){
	    radrefl1  = coeff1[lc];
	    radrefl11 = -2.0;
	  } else {
	    tau1      = tauc[nlayers[lu]+1]-tauc[lc+1];
	    etau1     = exp(-tau1);	  
	    tauetau1  = tau1*etau1;
	    
#if HAVE_LIBGSL 
	    /* exponential integral does not converge for values larger about 600*/
	    /* the integral for this case becomes 0*/
	    if (tau1 > 500.0)
	      ei1=0.0;
	    else
	      ei1       = gsl_sf_expint_Ei(-tau1);
#endif
	    
	    sqtau1    = tau1*tau1;
	    cutau1    = tau1*sqtau1;
	    
	    radrefl11 = -2.0*etau1+tauetau1-tau1*tauetau1-cutau1*ei1;	 

	    radrefl1  = etau1*coeff1[lc]*(1.-tau1)-coeff1[lc]*sqtau1*ei1;
	  }
	  
	  radrefl   = radrefl1 - radrefl2;

	  /* fprintf(stderr,"%s %4d %4d %13.6e\n","LC NLAYERS RADREFL[lc]",lc,nlayers[lu],radrefl);  */
	  /* fprintf(stderr,"%s %4d %13.6e\n","LC RADREFL1",lc,radrefl1);                          */
	  /* fprintf(stderr,"%s %4d %13.6e\n","LC RADREFL2",lc,radrefl2);                          */
	  /* fprintf(stderr,"%s %4d %13.6e %13.6e\n","LC COEFF1 COEFF2",lc,coeff1[lc],coeff2[lc]); */
	  /* fprintf(stderr,"%s %4d %13.6e %13.6e %13.6e %13.6e\n","LC TAU1 TAU2 SQTAU1 SQTAU2",lc,tau1,tau2,sqtau1,sqtau2); */
	  /* fprintf(stderr,"%s %4d %13.6e %13.6e\n","LC EI1 EI2",lc,ei1,ei2); */
			  
	  radrefl   += b1[lc]/3.*( radrefl11 - radrefl22 );
	  
//	  fprintf(stderr,"%s %4d %13.6e\n","LC RADREFL[lc]",lc,radrefl);
//	  fprintf(stderr,"%s %4d %13.6e\n","LC B1[lc]",lc,b1[lc]);      
//	  fprintf(stderr,"%s %4d %13.6e\n","LC RADREFL11",lc,radrefl11);
//	  fprintf(stderr,"%s %4d %13.6e\n","LC RADREFL22",lc,radrefl22);

	  /*              RAD[lu]   = RAD[lu]+RADEM[lc]+ALBEDO*ETAUMU*RADREFL      */
	  rad[lu]   += radrefl;
	  
	}
	
	/* emission of atmosphere */
	rademtot  += radem[lc];

      }

      /* fprintf(stderr,"lu = %4d   rad[lu] refl = %13.6e\n",lu,rad[lu]); */
      
      /* after computation of radiance reflected at the surface we multiply by the albedo and the extinction from bottom to TOA */
      /*              RAD[lu]   = RAD[lu]+RADEM[lc]+ALBEDO*ETAUMU*RADREFL      */
      rad[lu]   *= alb[lu]*etaumu;

      /* fprintf(stderr,"lu = %4d   rad[lu] refl attenuated = %13.6e\n",lu,rad[lu]); */

      /*          Contribution from atmospheric emission of all levels           */
      rad[lu]   += rademtot;

      /* fprintf(stderr,"lu = %4d   rad[lu] refl attenuated + emission = %13.6e\n",lu,rad[lu]); */

      /*          contribution from surface emission (bottom):     */
      /*          ********************************************     */
      rad[lu]   += bplanck[lu]*etaumu;

      /* fprintf(stderr,"lu = %4d   rad[lu] refl attenuated + emission + surface = %13.6e\n",lu,rad[lu]); */

      /* isotropic radiance                               */
      for (j=0;j<nphi;j++){
	uu[j][lu][iu] = (float) rad[lu];
      }
      uavg[lu]      = (float) rad[lu];
      /*       End loop over user tau          */
    }
    
    /*    End loop over user umu    */
  }
  
  free(radem);
  free(tauc);
  free(pkag);
  free(etau);
  free(b0);
  free(b1);
  free(coeff1);
  free(coeff2);

  return 0;
}



/*============================= errmsg_tzs() ===============================*/

/*
 * Print out a warning or error message;  abort if type == TZS_ERROR
 */

#define MAX_WARNINGS 100

void errmsg_tzs(char *messag,
		int   type)
{
  static int
    warning_limit = FALSE,
    num_warnings  = 0;

  if (type == TZS_ERROR) {
    fprintf(stderr,"\n ******* ERROR >>>>>>  %s\n",messag);
    exit(1);
  }

  if (warning_limit) return;

  if (++num_warnings <= MAX_WARNINGS) {
    fprintf(stderr,"\n ******* WARNING >>>>>>  %s\n",messag);
  }
  else {
    fprintf(stderr,"\n\n >>>>>>  TOO MANY WARNING MESSAGES --  They will no longer be printed  <<<<<<<\n\n");
    warning_limit = TRUE;
  }

  return;
}

#undef MAX_WARNINGS



void prtinp_tzs( int nlyr, int maxnlyr, float *dtauc, float *ssalb, float *temper,
		 float wvnmlo, float wvnmhi, int ntau, float *utau, int numu, float *umu,
		 int nphi, float *phi, float albedo, double btemp, float ttemp, float temis,
		 double *tauc, int nlev_c, float *zd_c, int nzout, float *zout )
{      
  /*     print values of input variables                                 */
  /*                                                                     */
  /*   Called by- TZS                                                    */
  /* --------------------------------------------------------------------*/

  /*     .. Local Scalars ..  */
  int lc, lu;

  fprintf( stderr, "\n%s%4d\n"," No. computational layers        = ", nlyr );
  fprintf( stderr, "\n%s%4d\n"," No. computational layers needed = ", maxnlyr );
  
  fprintf( stderr, "\n%4d %s", ntau, " User optical depths :");
  for (lu=0;lu<ntau;lu++)
    fprintf(stderr,"%13.6e",utau[lu]);
  fprintf(stderr,"\n");
  
  fprintf( stderr, "\n%4d %s", numu, " User polar angle cosines :");
  for (lu=0;lu<numu;lu++)
    fprintf(stderr,"%9.5f",umu[lu]);
  fprintf(stderr,"\n");
  
  fprintf( stderr, "\n%4d %s", nphi, " User azimuthal angles :");
  for (lu=0;lu<nphi;lu++)
    fprintf(stderr,"%9.2f",phi[lu]);
  fprintf(stderr,"\n");
  
  fprintf( stderr, "%s%14.4f%14.4f\n%s%10.2f%s%8.4f\n%s%8.4f\n", "    Thermal emission in wavenumber interval :", wvnmlo,
	   wvnmhi,"    Bottom temperature = ", btemp,"    Top temperature = ", ttemp,"    Top emissivity = ", temis );
  
  fprintf( stderr, "%s%8.4f\n","    Bottom albedo (Lambertian) =", albedo );
  
  /* Print layer variables */
  
  fprintf( stderr, "\n%s\n%s\n%s\n",  
	   "            Layer        Total        Single",
	   "          Optical      Optical    Scattering",
	   "            Depth        Depth        Albedo   Temperature" );
  
  for (lc=0;lc<maxnlyr;lc++){
    fprintf( stderr,"%4d%13.6e%13.6e%14.3f%14.3f\n",  lc, dtauc[lc], tauc[lc], ssalb[lc],temper[lc]);        
  }
  
  fprintf( stderr,"%4d%13.6e%13.6e%14.3f%14.3f\n",  maxnlyr, 0.0, tauc[maxnlyr], 0.0,temper[maxnlyr]);        

  fprintf( stderr, "\n%s\n%s\n%s\n",  
	   "            Layer        Single",
	   "          Optical    Scattering",
	   "            Depth        Albedo   Temperature" );
  
  for (lc=0;lc<nlyr;lc++){
    fprintf( stderr,"%4d%13.6e%14.3f%14.3f\n",  lc, dtauc[lc], ssalb[lc],temper[lc]);        
  }
  
  fprintf( stderr,"%4d%13.6e%14.3f%14.3f\n",  nlyr, -1.0, -1.0, temper[nlyr]);        
  
  fprintf( stderr, "\n%s\n%s\n%s\n","      Atmospheric","            Level","             [km]");
  
  for (lc=0;lc<nlev_c;lc++){
    fprintf( stderr,"%4d     %8.3f\n", lc, zd_c[lc]);
  }
  
  fprintf( stderr, "\n%s\n%s\n%s\n","           Output","            Level","             [km]");
  
  for (lc=0;lc<nzout;lc++){
    fprintf( stderr,"%4d     %8.3f\n", lc, zout[lc]);
  }
  
}




void chekin_tzs( int nlyr, int maxnlyr, float *dtauc, int nlev_c, float *zd_c,
		 int nzout, float *zout, 
		 float *ssalb, float *temper, 
		 float wvnmlo, float wvnmhi, 
		 int usrtau, int ntau, float *utau, 
		 int usrang, int numu, float *umu, int nphi, float *phi, 
		 float albedo, double btemp, float ttemp, 
		 float temis, int planck,
		 double *tauc, int quiet )
{      
  /*     Checks the input dimensions and variables  */
  int          inperr = FALSE;
  int          lc, j, iu, lu;

  if ( planck == FALSE ) 
    inperr = c_write_bad_var(VERBOSE, "planck" );
	  
  if ( nlyr<1 ) 
    inperr = c_write_bad_var(VERBOSE, "NLYR" );

  if ( maxnlyr<1 ) 
    inperr = c_write_bad_var(VERBOSE, "MAXNLYR" );
			     
  if (nlyr+1!=nlev_c) 
    inperr = c_write_bad_var(VERBOSE, "NLYR/NLEV_C" );

  for (lc = 0; lc < nlyr; lc++) {
    if (dtauc[lc] < 0.) {
      inperr = c_write_bad_var(VERBOSE,"dtauc");
    }
    if (zd_c[lc] < 0.) {
      inperr = c_write_bad_var(VERBOSE,"zd");
    }
    if (ssalb[lc] < 0.0 || ssalb[lc] > 1.0) {
      inperr = c_write_bad_var(VERBOSE,"ssalb");
    }
    if (temper[lc] < 0.) {
      inperr = c_write_bad_var(VERBOSE,"temper");
    }
  }

  if (usrtau) {
    if (ntau < 1) {
      inperr = c_write_bad_var(VERBOSE,"ntau");
    }
    if ( nzout < ntau ) 
      inperr = c_write_bad_var(VERBOSE, "ntau/nzout" );
    for (lu = 0; lu < ntau; lu++) {
      /* Do a relative check to see if we are just beyond the bottom boundary */
      /* This might happen due to numerical rounding off problems.  ak20110224*/

      //      fprintf(stderr,"lu = %d  utau = %13.6e  tauc[maxnlyr] = %13.6e  zout = %f\n",lu,utau[lu],tauc[maxnlyr],zout[lu]);
      if (fabs(utau[lu]-tauc[maxnlyr]) <= 1.e-6*tauc[maxnlyr]) {
	utau[lu] = tauc[maxnlyr];
      }
      //      fprintf(stderr,"lu = %d  utau = %13.6e  tauc[maxnlyr] = %13.6e  zout = %f\n",lu,utau[lu],tauc[maxnlyr],zout[lu]);

      //      if(utau[lu] < 0. || utau[lu] > tauc[maxnlyr]) {
      if (utau[lu] < 0. || (utau[lu] - tauc[maxnlyr]) > 1.e-6*tauc[maxnlyr] ) {
	fprintf(stderr,"WARNING: lu = %d  utau = %13.6e  tauc[maxnlyr] = %13.6e  zout = %f\n",lu,utau[lu],tauc[maxnlyr],zout[lu]);
	inperr = c_write_bad_var(VERBOSE,"utau");
      }
      if ( zout[lu]<0.0 ) 
	inperr = c_write_bad_var(VERBOSE, "zout" );
    }
  } else {         
    ntau=1;
    utau = (float *) calloc(1,sizeof(float));
    utau[0]=0.;    
    nzout=1;
    zout = (float *) calloc(1,sizeof(float));
    zout[0]=zd_c[0];
  }      

  if (usrang){
    if (numu < 0) {
      inperr = c_write_bad_var(VERBOSE,"numu");
    }
    for (iu = 0; iu < numu; iu++) {
      if (umu[iu] < -1. || umu[iu] > 1. || umu[iu] == 0.) {
	inperr = c_write_bad_var(VERBOSE,"umu");
      }
      if (iu >= 1) {
	if (umu[iu] < umu[iu-1]) {
	  inperr = c_write_bad_var(VERBOSE,"umu");
	}
      }
    }
    if ( nphi<=0 ) 
      inperr = c_write_bad_var(VERBOSE,"nphi" );
    
    for (j=0; j <nphi; j++) {
      if (phi[j] < 0. || phi[j] > 360.) {
	inperr = c_write_bad_var(VERBOSE,"phi");
      }
    }
  } else {      
    numu=1;
    umu = (float *) calloc(1,sizeof(float));
    umu[0]=1.;
    nphi=1;
    phi = (float *) calloc(1,sizeof(float));
    phi[0]=0.;
  }

  if (albedo < 0. || albedo > 1.)
    inperr = c_write_bad_var(VERBOSE,"albedo");
  
  
  if ( wvnmlo<0.0 || wvnmhi<=wvnmlo )
    inperr = c_write_bad_var(VERBOSE, "wvnmlo,hi" );
			     
  if (temis < 0. || temis > 1.)
    inperr = c_write_bad_var(VERBOSE,"temis");
    
  if ( btemp < 0.0 ) 
    inperr = c_write_bad_var(VERBOSE, "btemp" );
      
  if ( ttemp < 0.0 ) 
    inperr = c_write_bad_var(VERBOSE, "ttemp" );

  if ( inperr )
    c_errmsg( "TZS--input and/or dimension errors", TZS_ERROR );

  if (!quiet){
    for (lc = 1; lc <= nlyr; lc++) {
      if (fabs(temper[lc]-temper[lc-1]) > 10.) {
	c_errmsg("check_inputs--vertical temperature step may be too large for good accuracy",TZS_WARNING);
      }
    }
  }

}
