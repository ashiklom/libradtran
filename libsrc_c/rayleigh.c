/*--------------------------------------------------------------------
 * $Id: rayleigh.c 2839 2013-01-31 15:13:42Z svn-kylling $
 * 
 * This file is part of libRadtran.
 * Copyright (c) 1997-2012 by Arve Kylling, Bernhard Mayer,
 *                            Claudia Emde, Robert Buras
 *
 * ######### Contact info: http://www.libradtran.org #########
 *
 * This program is free software; you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License   
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.        
 * 
 * This program is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of  
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the   
 * GNU General Public License for more details.                    
 * 
 * You should have received a copy of the GNU General Public License          
 * along with this program; if not, write to the Free Software                
 * Foundation, Inc., 59 Temple Place - Suite 330, 
 * Boston, MA 02111-1307, USA.
 *--------------------------------------------------------------------*/

#include <stdio.h>							    
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef PI
#define PI 3.14159265358979323846264338327
#endif


/***********************************************************************************/
/* Function: crs_rayleigh_bodhaine                                                 */
/* Description:                                                                    */
/*    Calculates the Rayleigh scattering cross section according to the formula    */
/*    by Bodhaine et al, `On Rayleigh optical depth calculations', J. Atm. Ocean   */
/*    Technol., 16, 1854-1861, 1999.                                               */
/*    If a depolarization factor is given as input (depol_user>=0), the Rayleigh   */
/*    cross sections are calculated accordingly, otherwise the depolarization      */
/*    is calculated.                                                               */
/*                                                                                 */
/* Parameters:                                                                     */
/*                                                                                 */
/*  Input:                                                                         */
/*  float* lambda:          Pointer to ary holding the wavelength in nanometers.   */
/*  int n_lambda:           Number of wavelength points.                           */
/*  float mixing_ratio_co2: CO2 mixing ratio                                       */
/*  float depol_user:       User defined depolarization factor                     */    
/*                                                                                 */
/* Output:                                                                         */
/*  float **depol:          Calculated depolarization factor                       */
/*  float **crs:            Rayleigh cross section                                 */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling ????                                                       */
/*         2010-08-10 C. Emde (small modifications)                                */
/*                                                                                 */
/***********************************************************************************/

int crs_rayleigh_bodhaine (/* Input */
                           float *lambda, 
                           int n_lambda, 
                           float mixing_ratio_co2, 
			   float depol_user,                                                             
                           float **depol, 
                           float **crs) 
{
  double from_nm_to_cm = 1.0e-7;  
  double from_nm_to_um = 1.0e-3;
  double N_s           = 2.546899e+19; 
  double ray_const     = 24*pow(PI,3)/pow(N_s,2);
  double F_N_2=0, F_O_2=0, F_air=0;
  double co2 =0;
  double ref_ratio = 0, n_300=0, n=0;  /* n = refractive index of air at desired */ 
                                       /*     co2 concentration                  */
  double lambda_cm=0.0, lambda_um=0.0;
  int iv=0;
  
  /* fprintf(stderr, " mixing_ratio_co2 %g depol %g \n" ,  mixing_ratio_co2, depol_user);  */
  
  co2 = mixing_ratio_co2*1.e-4; /* This is co2 in parts per volume by percent, top of page */
                                /* 1858 of Bodhaine et al paper.                           */

  for (iv=0; iv<n_lambda; iv++) {

    lambda_cm = lambda[iv] * from_nm_to_cm;   /* Funny to use two different units */
    lambda_um = lambda[iv] * from_nm_to_um;   /* for lambda in the same paper.....*/

    n_300     = (8060.51 + (2480990/(132.274-pow(lambda_um,-2)))
		 + (17455.7/(39.32957-pow(lambda_um,-2))))*1e-08 ;  /* Eq. (18) */
    n         = (1+0.54*(mixing_ratio_co2*1.e-6-0.0003))*n_300 + 1;  /* Eq. (19) */
    ref_ratio = pow(n*n-1,2)/pow(n*n+2,2);

    F_N_2 = 1.034 + 3.17e-4/pow(lambda_um,2);                                /* Eq. (5) */ 
    F_O_2 = 1.096 + 1.385e-3/pow(lambda_um,2) + 1.448e-4/pow(lambda_um,4);   /* Eq. (6) */
    F_air = (78.084*F_N_2 + 20.946*F_O_2 + 0.934*1.0 + co2*1.15)/(78.084+20.946+0.934+co2); /* Eq. (23) */
    
    if (depol_user >= 0){
      F_air = (6.0+3.0 * depol_user) / (6.0-7.0* depol_user);
      (*depol)[iv]=depol_user; 
    }
    else
      (*depol)[iv] = 6.0 * (F_air-1.0) / (3.0+7.0*F_air);
    
    (*crs)[iv] =  (ray_const / pow(lambda_cm,4)) * ref_ratio * F_air; 
  }

  return 0;
}

/***********************************************************************************/
/* Function: crs_rayleigh_bodhaine29                                               */
/* Description:                                                                    */
/*    Calculates the Rayleigh scattering cross section according to Eq. 29         */
/*    of Bodhaine et al, `On Rayleigh optical depth calculations', J. Atm. Ocean   */
/*    Technol., 16, 1854-1861, 1999.                                               */
/*                                                                                 */
/* Parameters:                                                                     */
/*                                                                                 */
/*  Input:                                                                         */
/*  float* lambda:          Pointer to ary holding the wavelength in nanometers.   */
/*  int n_lambda:           Number of wavelength points.                           */
/*                                                                                 */
/* Output:                                                                         */
/*  float **crs:            Rayleigh cross section                                 */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling ????                                                       */
/*         2010-08-10 C. Emde (small modifications)                                */
/*                                                                                 */
/***********************************************************************************/

int crs_rayleigh_bodhaine29 (/* Input */
                           float *lambda, 
                           int n_lambda, 
                           float **crs) 
{
  double from_nm_to_um = 1.0e-3;
  double lambda_um=0.0, lambda_sq=0.0;
  int iv=0;
  

  for (iv=0; iv<n_lambda; iv++) {

    lambda_um = lambda[iv] * from_nm_to_um; 

    lambda_sq = lambda_um*lambda_um;

    (*crs)[iv] = 1.e-28 * (1.0455996 - 341.29061/lambda_sq - 0.90230850*lambda_sq)/
      (1 + 0.0027059889/lambda_sq - 85.968563*lambda_sq);

  }

  return 0;
}


/***********************************************************************************/
/* Function: crs_rayleigh_nicolet                                                  */
/* Description:                                                                    */
/*    Calculates the Rayleigh scattering cross section according to the formula    */
/*    by Nicolet, M., `On the molecular scattering in the terrestrial atmosphere:  */
/*    an empirical formula for its calculation in the homosphere', Planet. Space   */
/*    Sci., 32, 1467-1468, 1984.                                                   */
/*                                                                                 */
/* Parameters:                                                                     */
/*  float* lambda:         Pointer to ary holding the wavelength in nanometers.    */
/*  int n_lambda:          Number of wavelength points.                            */
/*                                                                                 */
/* Return value:                                                                   */
/*  double pointer to ary holding the Rayleigh scattering cross section as a       */
/*  function of wavelength.                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                                 */
/***********************************************************************************/

int crs_rayleigh_nicolet (float *lambda, int n_lambda, float **crs) 
{
  float ray_const  = 4.02e-28;
  float from_nm_um = 1.0e-3;
  float lambda_um=0.0, x=0.0;
  int iv=0;
  
  for (iv=0; iv<n_lambda; iv++) {

    lambda_um = lambda[iv] * from_nm_um;

    if (lambda_um > 0.55)
      x = 0.04;
    else
      x = 0.389*lambda_um + 0.09426/lambda_um - 0.3228;
    
    (*crs)[iv] = ray_const / pow(lambda_um,4+x);
  }

  return 0;
}

/***********************************************************************************/
/* Function: crs_rayleigh_penndorf                                                 */
/* Description:                                                                    */
/*    Calculates the Rayleigh scattering cross section according to Penndorf.      */
/*    The crs formula is taken from the Appendix of Vountas et al., 'Ring effect:  */
/*    Impact of Rotational Raman scattering on radiative transfer in Earth's       */
/*    atmosphere', JQSRT, 60, 943-961, 1998. Depolarization and refractive index   */
/*    calculated according to Bodhaine et al, `On Rayleigh optical depth           */
/*    calculations', J. Atm. Ocean Technol., 16, 1854-1861, 1999.                  */
/*    If a depolarization factor is given as input (depol_user>=0), the Rayleigh   */
/*    cross sections are calculated accordingly, otherwise the depolarization      */
/*    is calculated.                                                               */
/*                                                                                 */
/* Parameters:                                                                     */
/*                                                                                 */
/*  Input:                                                                         */
/*  float* lambda:          Pointer to ary holding the wavelength in nanometers.   */
/*  int n_lambda:           Number of wavelength points.                           */
/*  float mixing_ratio_co2: CO2 mixing ratio                                       */
/*  float depol_user:       User defined depolarization factor                     */    
/*                                                                                 */
/* Output:                                                                         */
/*  float **depol:          Calculated depolarization factor                       */
/*  float **crs:            Rayleigh cross section                                 */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling ???                                                        */
/*         2010-08-02 C. Emde (small modifications)                                */
/*                                                                                 */
/***********************************************************************************/

int crs_rayleigh_penndorf (float *lambda, int n_lambda, float mixing_ratio_co2, 
			   float depol_user, float **depol, float **crs) 
{
  double from_nm_to_cm = 1.0e-7;  
  double from_nm_to_um = 1.0e-3;
  double N_s           = 2.546899e+19; 
  double ray_const     = 32*pow(PI,3)/(3*pow(N_s,2));
  double F_N_2=0, F_O_2=0, F_air=0;
  double co2 =0;
  double n_300=0, n=0;           /* n = refractive index of air at desired */ 
                                 /*     co2 concentration                  */
  double lambda_cm=0.0, lambda_um=0.0;
  int iv=0;

  co2 = mixing_ratio_co2*1.e-4; /* This is co2 in parts per volume by percent, top of page */
                                /* 1858 of Bodhaine et al paper.                           */

  for (iv=0; iv<n_lambda; iv++) {

    lambda_cm = lambda[iv] * from_nm_to_cm;
    lambda_um = lambda[iv] * from_nm_to_um;

    n_300     = (8060.51 + (2480990/(132.274-pow(lambda_um,-2)))
		 + (17455.7/(39.32957-pow(lambda_um,-2))))*1e-08 ;  /* Eq. (18) */
    n         = (1+0.54*(mixing_ratio_co2*1.e-6-0.0003))*n_300 + 1;  /* Eq. (19) */

    F_N_2 = 1.034 + 3.17e-4/pow(lambda_um,2);                                /* Eq. (5) */ 
    F_O_2 = 1.096 + 1.385e-3/pow(lambda_um,2) + 1.448e-4/pow(lambda_um,4);   /* Eq. (6) */
    F_air = (78.084*F_N_2 + 20.946*F_O_2 + 0.934*1.0 + co2*1.15)/(78.084+20.946+0.934+co2); /* Eq. (23) */
    
    if (depol_user >= 0){
      F_air = (6.0+3.0*depol_user) / (6.0-7.0*depol_user);
      (*depol)[iv]=depol_user; 
    }
    else
      (*depol)[iv] = 6.0 * (F_air-1.0) / (3.0+7.0*F_air);

    
    (*crs)[iv] =  ( ray_const *  pow(n-1,2)/pow(lambda_cm,4) ) *
      ((6.0+3.0*(*depol)[iv]) / (6.0-7.0*(*depol)[iv]));

  }

  return 0;
}


