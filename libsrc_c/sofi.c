/************************************************************************
 * $Id: sofi.c 3313 2017-12-11 14:27:51Z Claudia.Emde $
 *
 * MYSTIC - Monte Carlo code for the physically correct tracing of
 *          photons in cloudy atmospheres.
 *
 * Copyright (c) 2000-2012 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras
 *
 * Correspondence: bernhard.mayer@lmu.de
 *
 ************************************************************************/

/************************************************************************/
/*                                                                      */
/* Functions to calculate radiation in the umbral shadow during the     */
/* total eclipse at the 29.03.2006 using MYSTIC.                        */
/*                                                                      */   
/* Author: Claudia Emde                                                 */
/*----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <errno.h>

#include "ascii.h"
#include "sofi.h"
#include "uvspecrandom.h"
#include "mystic.h"


#ifndef PI
#define PI 3.14159265358979323846264338327
#endif


/* prototypes */
/*
static photon_struct *generate_photon_sofi(int source,
					   atmosphere_struct *atmos, 
					   albedo_struct *albedo,
					   sample_struct *sample, 
					   double *pd,
					   double sza, double phi,
					   int ipa);
*/



/***********************************************************************************/
/* Function: sample_photons_sofi_060329                                   @62_30i@ */
/* Description:                                                                    */
/*  Generate a probability distribution functions which corresponds to the         */
/*  irradiance distribution at the top pf the atmosphere during the total          */
/*  eclipse of the sun at 29-03-2006 in Turkey.                                    */
/*  The integrated probability is returned.                                        */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/
int sample_photons_sofi( float wvnmlo, float wvnmhi, double *pd )
{
  
  double *x,*r;
  double *Ic, *Sl, *alpha; 
  int  i=0, j=0;
  double I_int=0.;
  FILE *file=NULL;;
  
  /* Definition of constants:*/ 
  double planck=6.63e-34;
  double speed_of_light=2.998e8;
  double boltzmann=1.38e-23;
  double T_S=5740.0;
  
  /* Wavelength */
  double lambda=2.0/((wvnmlo+wvnmhi)*100.0);

  /* Sampling will be done in km steps, counted from X=0, corresponding to the 
     centre of the umbral shadow. */
  /* int Nd=5000; */
  int Nx=10000;
  
  /* Number of sun radius grid points */
  int Nr=1001;
  
  /*extraterrestrial solar irradiance */
  double E_0=1368.0;
  
  /* Intensity in centre of solar disk */
  double I_0;
  
  /* ratio between moon diameter and sun diameter */
  double ratio=1.0494; /*   29-03-2006 in Turkey.*/
  /*double ratio=0.965;  Ratio for solar eclipse over south pole */ 
  
  /* Sun disk radius normalized to x */
  double r_M=ratio;
  
  /* calculated from data from Nasa report, corresponds to approx 1km */
  double dx=0.0006; 
  /*double dx = 5.56e-4; */  
  
  /* Limb darkening coefficient*/
  double beta=3.0*planck*speed_of_light*sqrt(sqrt(2.0))/
    (8.0*boltzmann* lambda * T_S); 

  double epsilon=1.e-8;
  
  /* Distance between centre of moon and centre of sun*/ 
  x = (double*) calloc(Nx+1,sizeof(double));
  /* Irradiance as a function of x*/
  Ic = (double*) calloc(Nx+1,sizeof(double));
  
  /* Initialize x. This is *not* X in Koepke! Normalisation not included. 
     X=x/(r_M+r_S) */
  x[0]=0.;
  for (j=1; j < Nx+1; j++)
    x[j]= x[j-1] + dx;
  
  /* Grid for radius on solar disk, needed for integration of apparant irradiance. 
     Also in units of apparent moon radius.*/
  r = (double*) calloc(Nr+1,sizeof(double));
  /* Angle for integration of solar disk.*/
  alpha = (double*) calloc(Nr+1,sizeof(double));
  /* Solar limb darkening */
  Sl = (double*) calloc(Nr+1,sizeof(double)); 

  r[0]=0.;
   Sl[0]=1.;
  
  for (i=1; i<Nr+1; i++){
    
      r[i] = r[i-1]+1.0/(Nr-1);
      
      Sl[i] = (1+beta*sqrt(fabs(1.-r[i]*r[i])))/(1+beta);
    }     
  
  /* Print solar irradiance of uncovered solar disk*/ 
  /*for (i=0; i<Nr; i++)
    fprintf(stderr, " %g %g \n",  r[i], Sl[i]);*/
  
  /* Integrate over solar disk to get I_0, take numerical integration
     for consistency */
  for (i=1; i<Nr; i++)
    I_int+=2*PI*(Sl[i]*r[i]+Sl[i+1]*r[i+1])*0.5*(r[i+1]-r[i]);
  
  /* Intensity in centre of solar disk */
  I_0=E_0/I_int;
  
  /* Calculate irradiance as a function of X*/
  for (j=0; j<Nx; j++){
    
    /* Moon has passed the sun, so we get the total intensity */
    if( x[j] >= r_M + 1.0)
      Ic[j] = E_0; 
    
    /* Moon in front of sun, integration over partly covered disk*/
    else{
      
      for (i=0; i<Nr; i++){
        
        if ( x[j] < r_M + 1.0 ){
          
          if ( ( r[i] < (x[j] - r_M) && x[j] > r_M ) ||
               ( r[i] > ( r_M + x[j] ) && x[j] < 1.0 - r_M ) ) 
            /* second part applies for ring eclipse */
            
            alpha[i] = PI;
          
          else if ( r[i] < (r_M - x[j])+epsilon && x[j] < r_M+epsilon )
            alpha[i] = 0.0;
          
          else 
            alpha[i] = acos( ( r_M*r_M - r[i]*r[i] - x[j]*x[j] ) / ( 2*r[i]*x[j] ) );
          
        }
        else
          fprintf(stderr, "Forgotten case x %g r %g \n. This should not happen!!!",
                  x[j], r[i]);  
        
      }
      
      /* Integrate uncovered part */
      for (i=0; i<Nr; i++){
  
        /* Print alpha for x[j] */
        /* if( j==1 ) */
	/*   fprintf(stderr, " %g %g %g \n", x[j], r[i], alpha[i]); */
	
	Ic[j] += 2.0 * I_0 * ( ( alpha[i] * Sl[i] * r[i] ) + 
                               ( alpha[i+1] * Sl[i+1] * r[i+1] ) ) *
          0.5 * ( r[i+1] - r[i] );
        
      }
    }
  }
  
  /* Photon density function */
  for (j=0; j < Nx; j++)
    pd[j] = Ic[j]/E_0;
  
  /*  pd[i] = (Ic[i]+1.1e-7*I_0)/E_0;
      1e-7: Upper limit of corona intensity.*/
  
  /* write dummy profile file */
  if ( (file = fopen("solar_intensity_pd.dat", "w")) == NULL)  
    return ASCIIFILE_NOT_FOUND;
  
  /* Print intensity as a function of distance */
  for (j=0; j < Nx; j++)
    fprintf(file, "%g %g %g \n", x[j], Ic[j], pd[j]);
  fclose(file);
  
  /* Free variables.*/
  free(r);
  free(alpha);
  free(Sl);

  return 0;
}



/***********************************************************************************/
/* Function: generate_photon_sofi                                         @62_30i@ */
/* Description:                                                                    */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

void generate_photon_sofi ( atmosphere_struct* atmos, photon_struct* p,
                            double sza, double phi, double* pd )
{
  double sinphi=0, cosphi=0, cossza=0, sinsza=0;
  double alpha=0, x=0, y=0, r=0, dx=0, dy=0;
  int i_re=0;

  /* sine and cosine of solar zenith */
  cossza = cosd(sza);
  sinsza = sind(sza);

  /* Different definition of azimuth angle */
  phi = phi-180.;

  /* sine and cosine of solar azimuth */
  cosphi = cosd(phi);
  sinphi = sind(phi);

  /* Center of umbral shadow at TOA */
  r=tand(sza) * atmos->Z[atmos->Nz];
  dx=r*sinphi;
  dy=r*cosphi;

  /* x and y coordinates, (x,y)=0 corresponds to center of ellipse */
  x=p->x[0]-200e3+dx;
  y=p->x[1]-200e3+dy;

  /* Starting point in polar coordinates */
  r=sqrt(x*x+y*y);
  alpha = atand(x/y);

  /* Equivalent radius, if the ellipsoidal shadow is transformed back to a circle, 
     for which the probability density function is given.*/
  i_re=(int)fabs(sqrt(pow(r,2)/
                      (pow(cosd(alpha-phi),2)/pow(cossza,2)+
                       pow(sind(alpha-phi),2))))/1000.;

  /* photon weight */
  p->weight        *= pd[i_re];

  /* initialize direction using the specified solar zenith and azimuth angles */
  /* and copy direction to "initial direction" struct                         */
  init_direction (sinsza, cossza, sinphi, cosphi, &(p->dir));
  cp_direction (&(p->dir0), &(p->dir));

  /* vertical start position: TOA */
  /* Make always sure to call set_photon_z() AFTER the new photon     */
  /* direction has been assigned because the start index and position */
  /* might depend on direction                                        */
  set_photon_z (atmos->Z[atmos->Nz], atmos, p);
}
