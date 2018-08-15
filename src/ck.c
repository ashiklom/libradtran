/*--------------------------------------------------------------------
 * $Id: ck.c 3276 2017-07-04 14:16:36Z Claudia.Emde $
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

#include "uvspec.h"
#include "ckdfu.h"
#include "ascii.h"
#include "fortran_and_c.h"
#include "f77-uscore.h"
#include "avhrr_kratz.h"

#if HAVE_LIBNETCDF
#include <netcdf.h>
#endif



/* molecular constants */
#define M_AIR  28.977
#define M_O3   48.000
#define M_H2O  18.015

#define N_MOL_CKDFU 6.0235e23



/* prototypes of internal functions */
static int add_katopath (char *filename, int scheme);
static int x_section (int scheme, int ispcs, 
                      float ***press, float ***temp_avg, int nlev, int Nx, int Ny, 
                      int iband, char *path, int *ngauss, double ****xabs, 
		      double *weight, int quiet, int verbose);

static int x_section_h2o (int scheme, float ***press, float ***temp_avg, float ****dens,
			  float ****dens_avg,
                          int nlev, int Nx, int Ny, 
			  int iband, char *path, 
                          int *ngauss, double ****xabs, double *weight,
			  int quiet, int verbose);

int hunt(float *xx, int n, double *x, int *jlo);

static int kato_crs (float ***temp_avg, float ***press, float ****dens, 
		     float ****dens_avg, float *z, int nlev, int Nx, int Ny,  
                     ck_struct *crs, int iband, int ispcs, char *path, 
		     int quiet, int verbose, 
                     double ****extinction, double *weight, int *ngauss);

static int fu_profile (float ***temper, float ***press, float *z, 
		       float ****dens, float ****dens_avg, int nlev,
		       int Nx, int Ny,
                       int h2o_cont, float umco2, float umch4, float umn2o,
                       float umf11, float umf12, float umf22,
                       ck_profile *ck);

static int avhrr_kratz (int channel, int interval,
                        float *z0, float *p0, float *t0, float *u0, float *ux, 
                        int nlev, float co2, float n2o, float f11, float f12, float ch4, 
                        double **od, double *wght, int *ntau);
     

static int avhrr_kratz_profile (float ***temper, float ***press, 
				float *z, float ****dens, int nlev, 
                                float umco2, float umch4, float umn2o, 
				float umf11, float umf12,
                                ck_profile *ck);

static int generic_profile (char *filename, float *z, int nlev, int *nwght, ck_profile *ck);


/* complement filename for Kato et al. parameterization */
static int add_katopath (char *filename, int scheme)
{
  switch (scheme) {
  case CK_KATO:
    strcat(filename, "correlated_k/kato/");
    break;
  case CK_KATO2:
    strcat(filename, "correlated_k/kato2/");
    break;
  case CK_KATO2_96:
    strcat(filename, "correlated_k/kato2.hitran96/");
    break;
  case CK_KATO2ANDWANDJI:
    strcat(filename, "correlated_k/kato2andwandji/");
    break;
  default:
    fprintf (stderr, "Error, unsupported correlated_k scheme %d\n", scheme);
    return -1;
  }
  return 0;
}


/*****************************************************************/
/* Read the Kato et al. wavelength grid and cross section tables */
/* and write the data to ck_struct *ck                           */
/*****************************************************************/

int kato_readtables (int scheme, ck_struct *ck, char *path, char *photon_filename)
{
  int i=0, iband=0, status=0;
  int iv=0, iq1=0, iq2=0, iq3=0, iq4=0;
  int nq1=0, nq2=0, nq3=0, nq4=0;
  char filename[FILENAME_MAX]="";
  char solar_filename  [FILENAME_MAX]="";
  char thermal_filename[FILENAME_MAX]="";

  int n_tmp=0;
  double *tmp1=NULL, *tmp2=NULL, *tmp3=NULL, *tmp4=NULL, *tmp5=NULL, *tmp6=NULL;


  ck->scheme = scheme;

  /* read wavelength information */

  strcpy(filename, path);
  status = add_katopath (filename, scheme);
  if (status!=0)
    return status;
  strcat(filename, "wvl.dat");

  status = read_4c_file (filename, &tmp1, &tmp2, &tmp3, &tmp4, &(ck->n_wvl));
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }

  /* allocate memory for wavelengths */
  ck->wvlc  = calloc(ck->n_wvl+1, sizeof(double));
  ck->wvnlo = calloc(ck->n_wvl+1, sizeof(double)); 
  ck->wvnhi = calloc(ck->n_wvl+1, sizeof(double));


  /* 4 components (1,2,3,4) in Kato et al. [1999] */
  ck->comp = calloc (KATO_COMPONENTS+1, sizeof (double *));
  ck->wght = calloc (KATO_COMPONENTS+1, sizeof (int *));

  /* check structure (bands need to be sorted in ascending order  */
  /* and are assumed to start with 1) and copy data to wavelength */
  /* arrays.                                                      */

  for (i=0; i<ck->n_wvl; i++) {
    iband = (int) (tmp1[i]+0.5); 
    
    if (iband != i+1) {
      fprintf (stderr, " ... inconsistency in %s, band %d\n", filename, iband);
      return -1;
    }
    
    /* copy wavelength bands to final destination */
    (ck->wvlc)[iband]  = tmp2[i];
    (ck->wvnhi)[iband] = 1.0e7/tmp3[i];  /* convert nm to cm-1 */
    (ck->wvnlo)[iband] = 1.0e7/tmp4[i];  /* convert nm to cm-1 */
  }

  free(tmp1); free(tmp2); free(tmp3); free(tmp4);

  
  /* allocate memory for cross sections */
  ck->comp[KATO_CO2] = calloc (ck->n_wvl+1, sizeof(double));
  ck->comp[KATO_O3]  = calloc (ck->n_wvl+1, sizeof(double));
  ck->comp[KATO_O2]  = calloc (ck->n_wvl+1, sizeof(double));
  ck->comp[KATO_H2O] = calloc (ck->n_wvl+1, sizeof(double));
  
  ck->wght[KATO_CO2] = calloc (ck->n_wvl+1, sizeof(int));
  ck->wght[KATO_O3]  = calloc (ck->n_wvl+1, sizeof(int));
  ck->wght[KATO_O2]  = calloc (ck->n_wvl+1, sizeof(int));
  ck->wght[KATO_H2O] = calloc (ck->n_wvl+1, sizeof(int));
  
    
  /* read CO2 data */

  strcpy (filename, path);
  status = add_katopath (filename, scheme);
  if (status!=0)
    return status;
  strcat(filename, "co2.dat");

  status = read_3c_file (filename, &tmp1, &tmp2, &tmp3, &n_tmp);
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }
  
  for (i=0; i<n_tmp; i++) {
    iband = (int) (tmp1[i]+0.5); 
    if (iband > ck->n_wvl) {
      fprintf (stderr, " ... inconsistency in %s, band %d\n",
               filename, iband);
      return -1;
    }
    
    (ck->comp)[KATO_CO2][iband] = tmp3[i];
    (ck->wght)[KATO_CO2][iband] = (int) (tmp2[i]+0.5);
  }

  free(tmp1); free(tmp2); free(tmp3);



  /* read O3 data */
  strcpy (filename, path);
  status = add_katopath (filename, scheme);
  if (status!=0)
    return status;
  strcat(filename, "o3.dat");

  status = read_3c_file (filename, &tmp1, &tmp2, &tmp3, &n_tmp);
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }
  
  for (i=0; i<n_tmp; i++) {
    iband = (int) (tmp1[i]+0.5); 
    if (iband > ck->n_wvl) {
      fprintf (stderr, " ... inconsistency in %s, band %d\n",
               filename, iband);
      return -1;
    }
    
    (ck->comp)[KATO_O3][iband] = tmp3[i];
    (ck->wght)[KATO_O3][iband] = (int) (tmp2[i]+0.5);
  }

  free(tmp1); free(tmp2); free(tmp3);


  /* read O2 data */
  strcpy (filename, path);
  status = add_katopath (filename, scheme);
  if (status!=0)
    return status;
  strcat(filename, "o2.dat");

  status = read_3c_file (filename, &tmp1, &tmp2, &tmp3, &n_tmp);
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }
  
  for (i=0; i<n_tmp; i++) {
    iband = (int) (tmp1[i]+0.5); 
    if (iband > ck->n_wvl) {
      fprintf (stderr, " ... inconsistency in %s, band %d\n",
               filename, iband);
      return -1;
    }
    
    ck->comp[KATO_O2][iband] = tmp3[i];
    ck->wght[KATO_O2][iband] = (int) (tmp2[i]+0.5);
  }

  free(tmp1); free(tmp2); free(tmp3);


  /* read H2O data */
  strcpy (filename, path);
  status = add_katopath (filename, scheme);
  if (status!=0)
    return status;
  strcat(filename, "h2o.dat");

  status = read_3c_file (filename, &tmp1, &tmp2, &tmp3, &n_tmp);
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }
  
  for (i=0; i<n_tmp; i++) {
    iband = (int) (tmp1[i]+0.5); 
    if (iband > ck->n_wvl) {
      fprintf (stderr, " ... inconsistency in %s, band %d\n",
               filename, iband);
      return -1;
    }
    
    (ck->comp)[KATO_H2O][iband] = tmp3[i];
    (ck->wght)[KATO_H2O][iband] = (int) (tmp2[i]+0.5);
  }

  free(tmp1); free(tmp2); free(tmp3);


  /* read photon fractions */
  if (strlen(photon_filename)==0) {

    /* default solar photon file */
    strcpy (solar_filename, path);
    status = add_katopath (solar_filename, scheme);
    if (status!=0)
      return status;
    strcat(solar_filename, "x_solar.dat");

    /* default thermal photon file */
    strcpy(thermal_filename, path);
    status = add_katopath (thermal_filename, scheme);
    if (status!=0)
      return status;
    strcat(thermal_filename, "x_thermal.dat");
  }
  else {
    strcpy (solar_filename, photon_filename);
    strcpy (thermal_filename, photon_filename);
  }

  /* allocate memory for photon arrays */
  ck->x_solar   = calloc (ck->n_wvl+1, sizeof (double ****));
  ck->x_thermal = calloc (ck->n_wvl+1, sizeof (double ****));
  
  for (iv=1; iv<=ck->n_wvl; iv++) {
    
    nq1 = (ck->wght)[KATO_CO2][iv];
    nq2 = (ck->wght)[KATO_O3] [iv];
    nq3 = (ck->wght)[KATO_O2] [iv];
    nq4 = (ck->wght)[KATO_H2O][iv];

    ((double *****) ck->x_solar)  [iv] = calloc (nq1+1, sizeof(double ***));
    ((double *****) ck->x_thermal)[iv] = calloc (nq1+1, sizeof(double ***));
    
    for (iq1=1; iq1<=nq1; iq1++) {
      ((double *****) ck->x_solar)  [iv][iq1] = calloc (nq2+1, sizeof(double **));
      ((double *****) ck->x_thermal)[iv][iq1] = calloc (nq2+1, sizeof(double **));

      for (iq2=1; iq2<=nq2; iq2++) {
        ((double *****) ck->x_solar)  [iv][iq1][iq2] = calloc (nq3+1, sizeof(double *));
        ((double *****) ck->x_thermal)[iv][iq1][iq2] = calloc (nq3+1, sizeof(double *));
        
        for (iq3=1; iq3<=nq3; iq3++) {
          ((double *****) ck->x_solar)  [iv][iq1][iq2][iq3] = calloc (nq4+1, sizeof(double));
          ((double *****) ck->x_thermal)[iv][iq1][iq2][iq3] = calloc (nq4+1, sizeof(double));
        }
      }
    }
  }

  /* solar */
  status = read_6c_file (solar_filename, &tmp1, &tmp2, &tmp3, &tmp4, &tmp5, &tmp6, &n_tmp);
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, solar_filename);
    return status;
  }

  for (i=0; i<n_tmp; i++) {
    iv  = (int) (tmp1[i]+0.5); 
    iq1 = (int) (tmp2[i]+0.5); 
    iq2 = (int) (tmp3[i]+0.5); 
    iq3 = (int) (tmp4[i]+0.5); 
    iq4 = (int) (tmp5[i]+0.5); 

    ((double *****) ck->x_solar)[iv][iq1][iq2][iq3][iq4] = tmp6[i];
  }
  free(tmp1); free(tmp2); free(tmp3);
  free(tmp4); free(tmp5); free(tmp6);

  /* thermal */
  status = read_6c_file (thermal_filename, &tmp1, &tmp2, &tmp3, &tmp4, &tmp5, &tmp6, &n_tmp);
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, thermal_filename);
    return status;
  }

  for (i=0; i<n_tmp; i++) {
    iv  = (int) (tmp1[i]+0.5); 
    iq1 = (int) (tmp2[i]+0.5); 
    iq2 = (int) (tmp3[i]+0.5); 
    iq3 = (int) (tmp4[i]+0.5); 
    iq4 = (int) (tmp5[i]+0.5); 

    ((double *****) ck->x_thermal)[iv][iq1][iq2][iq3][iq4] = tmp6[i];
  }
  free(tmp1); free(tmp2); free(tmp3);
  free(tmp4); free(tmp5); free(tmp6);


  return 0;
}



/*****************************************************************/
/* Read the Fu et al. wavelength grid and cross section tables   */
/* and write the data to ck_struct *ck                           */
/*****************************************************************/

int fu_readtables (ck_struct *ck, char *path, char *photon_filename)
{
  int i=0, iband=0, iwght=0, status=0;
  int n_tmp=0;
  double *tmp1=NULL, *tmp2=NULL, *tmp3=NULL, *tmp4=NULL;
  char filename[FILENAME_MAX]="";

  char solar_filename  [FILENAME_MAX]="";
  char thermal_filename[FILENAME_MAX]="";


  ck->scheme = CK_FU;

  /* read wavelength  information */

  strcpy(filename, path);
  strcat(filename, "correlated_k/fu/wvl.dat");

  status = read_4c_file (filename, &tmp1, &tmp2, &tmp3, &tmp4, &(ck->n_wvl));
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }

  /* allocate memory for wavelengths */
  ck->wvlc  = calloc(ck->n_wvl+1, sizeof(double));
  ck->wvnlo = calloc(ck->n_wvl+1, sizeof(double));
  ck->wvnhi = calloc(ck->n_wvl+1, sizeof(double));

  /* Only one component is considered explicitely */
  ck->comp = calloc (FU_COMPONENTS+1, sizeof (double *));
  ck->wght = calloc (FU_COMPONENTS+1, sizeof (int *));

  /* check structure (bands need to be sorted in ascending order  */
  /* and are assumed to start with 1) and copy data to wavelength */
  /* arrays.                                                      */

  for (i=0; i<ck->n_wvl; i++) {
    iband = (int) (tmp1[i]+0.5); 
    
    if (iband != i+1) {
      fprintf (stderr, " ... inconsistency in %s, band %d\n", filename, iband);
      return -1;
    }
    
    /* copy and convert to nm */
    (ck->wvlc) [iband] = tmp2[i]*1000.0;
    (ck->wvnhi)[iband] = 1e4/tmp3[i];  /* convert micron to cm-1 */
    (ck->wvnlo)[iband] = 1e4/tmp4[i];  /* convert micron to cm-1 */
  }


  
  /* allocate memory for cross sections */
  ck->comp[FU_WGHT] = calloc (ck->n_wvl+1, sizeof(double));
  ck->wght[FU_WGHT] = calloc (ck->n_wvl+1, sizeof(int));

  /* free memory */
  free(tmp1); free(tmp2); free(tmp3); free(tmp4);
  
    
  /* read weights */
  strcpy(filename, path);
  strcat(filename, "correlated_k/fu/quad.dat");

  status = read_3c_file (filename, &tmp1, &tmp2, &tmp3, &n_tmp);
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }
  
  for (i=0; i<n_tmp; i++) {
    iband = (int) (tmp1[i]+0.5); 
    if (iband > ck->n_wvl) {
      fprintf (stderr, " ... inconsistency in %s, band %d\n",
               filename, iband);
      return -1;
    }
    
    (ck->comp)[FU_WGHT][iband] = tmp3[i];
    (ck->wght)[FU_WGHT][iband] = (int) (tmp2[i]+0.5);
  }

  /* free memory */
  free(tmp1); free(tmp2); free(tmp3);


  /* read photon fractions */
  if (strlen(photon_filename)==0) {
    /* default photon files */

    strcpy(solar_filename, path);
    strcat(solar_filename, "correlated_k/fu/x_solar.dat");

    strcpy(thermal_filename, path);
    strcat(thermal_filename, "correlated_k/fu/x_thermal.dat");
  }
  else {
    strcpy(solar_filename, photon_filename);
    strcpy(thermal_filename, photon_filename);
  }


  /* solar */
  status = read_3c_file (solar_filename, &tmp1, &tmp2, &tmp3, &n_tmp);
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, solar_filename);
    return status;
  }

  /* allocate memory */
  ck->x_solar = calloc (ck->n_wvl+1, sizeof (double *));
  for (i=0; i<=ck->n_wvl; i++)
    ((double **) ck->x_solar)[i] = calloc (FU_MAXINT+1, sizeof(double));

  for (i=0; i<n_tmp; i++) {
    iband = (int) (tmp1[i]+0.5); 
    iwght = (int) (tmp2[i]+0.5); 
    
    ((double **) ck->x_solar)[iband][iwght] = tmp3[i];
  }

  free(tmp1); free(tmp2); free(tmp3);

  /* thermal */
  status = read_3c_file (thermal_filename, &tmp1, &tmp2, &tmp3, &n_tmp);
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, thermal_filename);
    return status;
  }

  /* allocate memory */
  ck->x_thermal = calloc (ck->n_wvl+1, sizeof (double *));
  for (i=0; i<=ck->n_wvl; i++)
    ((double **) ck->x_thermal)[i] = calloc (FU_MAXINT+1, sizeof(double));

  for (i=0; i<n_tmp; i++) {
    iband = (int) (tmp1[i]+0.5); 
    iwght = (int) (tmp2[i]+0.5); 
    
    ((double **) ck->x_thermal)[iband][iwght] = tmp3[i];
  }

  free(tmp1); free(tmp2); free(tmp3);

  return 0;
}



/*****************************************************************/
/* Read the Kratz AVHRR wavelength grid and cross section tables */
/* and write the data to ck_struct *ck                           */
/*****************************************************************/

int avhrr_kratz_readtables (ck_struct *ck, char *path, char *photon_filename)
{
  int i=0, iband=0, iwght=0, status=0;
  char filename[FILENAME_MAX]="";
  int n_tmp=0;
  double *tmp1=NULL, *tmp2=NULL, *tmp3=NULL, *tmp4=NULL;

  char solar_filename  [FILENAME_MAX]="";
  char thermal_filename[FILENAME_MAX]="";


  ck->scheme = CK_AVHRR_KRATZ;

  /* read wavelength  information */

  strcpy(filename, path);
  strcat(filename, "correlated_k/kratz/avhrr_wvl.dat");

  status = read_4c_file (filename, &tmp1, &tmp2, &tmp3, &tmp4, &(ck->n_wvl));
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }

  /* allocate memory for wavelengths */
  ck->wvlc  = calloc(ck->n_wvl+1, sizeof(double));
  ck->wvnlo = calloc(ck->n_wvl+1, sizeof(double));
  ck->wvnhi = calloc(ck->n_wvl+1, sizeof(double));


  /* Only one component is considered explicitely */
  ck->comp = calloc (AVHRR_KRATZ_COMPONENTS+1, sizeof (double *));
  ck->wght = calloc (AVHRR_KRATZ_COMPONENTS+1, sizeof (int *));

  /* check structure (bands need to be sorted in ascending order  */
  /* and are assumed to start with 1) and copy data to wavelength */
  /* arrays.                                                      */

  for (i=0; i<ck->n_wvl; i++) {
    iband = (int) (tmp1[i]+0.5); 
    
    if (iband != i+1) {
      fprintf (stderr, " ... inconsistency in %s, band %d\n", filename, iband);
      return -1;
    }
    
    /* copy and convert to nm */
    (ck->wvlc)[iband] = tmp2[i]*1000.0;
    (ck->wvnhi)[iband] = 1.0e4/tmp3[i];  /* convert micron to cm-1 */
    (ck->wvnlo)[iband] = 1.0e4/tmp4[i];  /* convert micron to cm-1 */
  }

  free(tmp1); free(tmp2); free(tmp3); free(tmp4);

  
  /* allocate memory for cross sections */
  ck->comp[AVHRR_KRATZ_WGHT] = calloc (ck->n_wvl+1, sizeof(double));
  ck->wght[AVHRR_KRATZ_WGHT] = calloc (ck->n_wvl+1, sizeof(int));
  
    
  /* read weights */
  strcpy(filename, path);
  strcat(filename, "correlated_k/kratz/avhrr_quad.dat");

  status = read_3c_file (filename, &tmp1, &tmp2, &tmp3, &n_tmp);
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }
  
  for (i=0; i<n_tmp; i++) {
    iband = (int) (tmp1[i]+0.5); 
    if (iband > ck->n_wvl) {
      fprintf (stderr, " ... inconsistency in %s, band %d\n",
               filename, iband);
      return -1;
    }
    
    (ck->comp)[FU_WGHT][iband] = tmp3[i];
    (ck->wght)[FU_WGHT][iband] = (int) (tmp2[i]+0.5);
  }

  /* free memory */
  free(tmp1); free(tmp2); free(tmp3);


  /* read photon fractions */
  if (strlen(photon_filename)==0) {
    /* default photon files */

    strcpy(solar_filename, path);
    strcat(solar_filename, "correlated_k/kratz/x_solar.dat");

    strcpy(thermal_filename, path);
    strcat(thermal_filename, "correlated_k/kratz/x_thermal.dat");
  }
  else {
    strcpy(solar_filename, photon_filename);
    strcpy(thermal_filename, photon_filename);
  }


  /* solar */
  status = read_3c_file (solar_filename, &tmp1, &tmp2, &tmp3, &n_tmp);
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, solar_filename);
    return status;
  }

  /* allocate memory */
  ck->x_solar = calloc (ck->n_wvl+1, sizeof (double *));
  for (i=0; i<=ck->n_wvl; i++)
    ((double **) ck->x_solar)[i] = calloc (AVHRR_KRATZ_MAXINT+1, sizeof(double));

  for (i=0; i<n_tmp; i++) {
    iband = (int) (tmp1[i]+0.5); 
    iwght = (int) (tmp2[i]+0.5); 
    
    ((double **) ck->x_solar)[iband][iwght] = tmp3[i];
  }


  /* thermal */
  status = read_3c_file (thermal_filename, &tmp1, &tmp2, &tmp3, &n_tmp);
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, thermal_filename);
    return status;
  }

  /* allocate memory */
  ck->x_thermal = calloc (ck->n_wvl+1, sizeof (double *));
  for (i=0; i<=ck->n_wvl; i++)
    ((double **) ck->x_thermal)[i] = calloc (AVHRR_KRATZ_MAXINT+1, sizeof(double));

  for (i=0; i<n_tmp; i++) {
    iband = (int) (tmp1[i]+0.5); 
    iwght = (int) (tmp2[i]+0.5); 
    
    ((double **) ck->x_thermal)[iband][iwght] = tmp3[i];
  }

  return 0;
}



/*****************************************************************/
/* Read a generic CK table and write the data to ck_struct *ck   */
/*****************************************************************/

int ck_generic_readtables (ck_struct *ck, char *filename, char *photon_filename)
{

#if HAVE_LIBNETCDF

  int i=0, iq=0, iv=0, iband=0, iwght=0, status=0;
  int n_tmp=0;
  double *tmp1=NULL, *tmp2=NULL, *tmp3=NULL;

  int ncid=0;
  int idd_nwvn=0, idd_nlyr=0, idd_maxc=0;
  int id_wvl=0, id_wvnlo=0, id_wvnhi=0, id_nc=0;

  size_t dimlen=0;

  size_t start[1] = {0};
  size_t count[1] = {0};

  double *wvl=NULL, *wvnlo=NULL, *wvnhi=NULL;
  int *nc=NULL;

  ck->scheme = CK_FILE;

  /* open netcdf file */
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d opening netCDF file %s\n", status, filename);
    return status;
  }


  /* get dimension id for "nwvn" */
  status = nc_inq_dimid (ncid, "nwvn", &idd_nwvn);

  /* get dimension length for "nwvn" */
  status = nc_inq_dimlen (ncid, idd_nwvn, &dimlen);
  ck->n_wvl = dimlen;


  /* get dimension id for "nlyr" */
  status = nc_inq_dimid (ncid, "nlyr", &idd_nlyr);
  
  /* get dimension length for "nlyr" */
  status = nc_inq_dimlen (ncid, idd_nlyr, &dimlen);

  /* get dimension id for "maxc" */
  status = nc_inq_dimid (ncid, "maxc", &idd_maxc);
  
  /* get dimension length for "maxc" */
  status = nc_inq_dimlen (ncid, idd_maxc, &dimlen);

  count[0] = ck->n_wvl;
  
  /* allocate memory for fields */
  wvl   = calloc (ck->n_wvl, sizeof(double));
  wvnlo = calloc (ck->n_wvl, sizeof(double));
  wvnhi = calloc (ck->n_wvl, sizeof(double));
  nc    = calloc (ck->n_wvl, sizeof(int));

  ck->wvnlo = calloc (ck->n_wvl+1, sizeof(double));
  ck->wvnhi = calloc (ck->n_wvl+1, sizeof(double));
  ck->wvlc  = calloc (ck->n_wvl+1, sizeof(double));
  

  /* get variable id for "wvl" */
  status = nc_inq_varid (ncid, "wvl", &id_wvl);
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d reading %s: double wvl() missing!\n", status, filename);
    return status;
  }
  

  /* read "wvl" */
  status = nc_get_vara_double (ncid, id_wvl, start, count, wvl);
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }

  /* get variable id for "wvnlo" */
  status = nc_inq_varid (ncid, "wvnlo", &id_wvnlo);

  /* read "wvnlo" */
  status = nc_get_vara_double (ncid, id_wvnlo, start, count, wvnlo);
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }

  /* get variable id for "wvnhi" */
  status = nc_inq_varid (ncid, "wvnhi", &id_wvnhi);

  /* read "wvnhi" */
  status = nc_get_vara_double (ncid, id_wvnhi, start, count, wvnhi);
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }


  for (iv=1; iv<=ck->n_wvl; iv++) {
    ck->wvlc[iv]  = wvl  [iv-1];
    ck->wvnlo[iv] = wvnlo[iv-1];
    ck->wvnhi[iv] = wvnhi[iv-1];
  }


  ck->comp = calloc (GENERIC_COMPONENTS+1, sizeof(double *));
  ck->wght = calloc (GENERIC_COMPONENTS+1, sizeof(int *));

  ck->comp[GENERIC_WGHT] = calloc (ck->n_wvl+1, sizeof(double));
  ck->wght[GENERIC_WGHT] = calloc (ck->n_wvl+1, sizeof(int));


  /* get variable id for "nc" */
  status = nc_inq_varid (ncid, "nc", &id_nc);
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }

  /* read "nc" */
  status = nc_get_vara_int (ncid, id_nc, start, count, nc);
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }

  for (iv=1; iv<=ck->n_wvl; iv++) {
    ck->comp[GENERIC_WGHT][iv] = -1.0;
    ck->wght[GENERIC_WGHT][iv] = nc[iv-1];
  }

  free (wvl);
  free (wvnlo);
  free (wvnhi);
  free (nc);

  nc_close(ncid);

  /*************************/
  /* read photon fractions */
  /*************************/

  /* allocate memory */
  ck->x_solar   = calloc (ck->n_wvl+1, sizeof (double *));
  ck->x_thermal = calloc (ck->n_wvl+1, sizeof (double *));
  for (iv=1; iv<=ck->n_wvl; iv++) {
    ((double **) ck->x_solar)  [iv] = calloc (ck->wght[GENERIC_WGHT][iv]+1, sizeof(double));
    ((double **) ck->x_thermal)[iv] = calloc (ck->wght[GENERIC_WGHT][iv]+1, sizeof(double));

    /* default: distribute evenly */
    for (iq=1;iq<=ck->wght[GENERIC_WGHT][iv];iq++)
      ((double **) ck->x_solar)[iv][iq] = 1.0 / (double) ck->wght[GENERIC_WGHT][iv];
  }

  if (strlen(photon_filename)>0) {
    status = read_3c_file (photon_filename, &tmp1, &tmp2, &tmp3, &n_tmp);
    if (status!=0) {
      fprintf (stderr, "Error %d reading %s\n", status, photon_filename);
      return status;
    }

    for (i=0; i<n_tmp; i++) {
      iband = (int) (tmp1[i]+0.5); 
      iwght = (int) (tmp2[i]+0.5); 
      
      ((double **) ck->x_solar)  [iband][iwght] = tmp3[i];
      ((double **) ck->x_thermal)[iband][iwght] = tmp3[i];
    }
  }


#else
  fprintf (stderr, " ***********************************************************************\n");
  fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
  fprintf (stderr, " * use the generic correlated_k option. Please get netcdf and rebuild. *\n");
  fprintf (stderr, " ***********************************************************************\n");
  return -1;
#endif
  
  return 0;  /* ok */
}




int sbdart_readtables (ck_struct *ck, int nwvl, char *path, char *photon_filename, int quiet)
{
  int iv=0, iq=0;
  int i=0, iband=0, iwght=0, status=0;
  int n_tmp=0;
  double *tmp1=NULL, *tmp2=NULL, *tmp3=NULL;
  
  ck->scheme = CK_LOWTRAN;
  ck->n_wvl = nwvl;

  if (!quiet)
    fprintf (stderr, " ... using %d wavelengths internally\n", nwvl);


  /* read photon fractions */
  if (strlen(photon_filename)==0) {

    /* if no photon distribution is given,     */
    /* allocate memory and set all number to 1 */

    ck->x_solar = calloc (ck->n_wvl+1, sizeof (double *));
    for (iv=0; iv<=ck->n_wvl; iv++) {
      ((double **) ck->x_solar)[iv] = calloc (LOWTRAN_MAXINT+1, sizeof(double));
      
      for (iq=0; iq<=LOWTRAN_MAXINT; iq++)
        ((double **) ck->x_solar)[iv][iq] = 1;
    }
    
    ck->x_thermal = calloc (ck->n_wvl+1, sizeof (double *));
    for (iv=0; iv<=ck->n_wvl; iv++) {
      ((double **) ck->x_thermal)[iv] = calloc (LOWTRAN_MAXINT+1, sizeof(double));
      
      for (iq=0; iq<=LOWTRAN_MAXINT; iq++)
        ((double **) ck->x_thermal)[iv][iq] = 1;
    }
  }
  else {
    /* solar */
    fprintf (stderr, " ... reading solar photon numbers from %s\n",
             photon_filename);
    status = read_3c_file (photon_filename, &tmp1, &tmp2, &tmp3, &n_tmp);
    if (status!=0) {
      fprintf (stderr, "Error %d reading %s\n", status, photon_filename);
      return status;
    }
    
    /* allocate memory */
    ck->x_solar = calloc (ck->n_wvl+1, sizeof (double *));
    for (i=0; i<=ck->n_wvl; i++)
      ((double **) ck->x_solar)[i] = calloc (LOWTRAN_MAXINT+1, sizeof(double));
    
    for (i=0; i<n_tmp; i++) {
      iband = (int) (tmp1[i]+0.5); 
      iwght = (int) (tmp2[i]+0.5); 
      
      if (iband<=0 || iband >= ck->n_wvl+1) {
	fprintf (stderr, "Error, band %d out of range in %s\n", iband, photon_filename);
	return -1;
      }

      if (iwght<=0 || iwght >= LOWTRAN_MAXINT+1) {
	fprintf (stderr, "Error, quadrature point %d out of range in %s\n", iwght, photon_filename);
	return -1;
      }

      ((double **) ck->x_solar)[iband][iwght] = tmp3[i];
    }
    

    /* thermal */
    fprintf (stderr, " ... reading thermal photon numbers from %s\n",
             photon_filename);

    status = read_3c_file (photon_filename, &tmp1, &tmp2, &tmp3, &n_tmp);
    if (status!=0) {
      fprintf (stderr, "Error %d reading %s\n", status, photon_filename);
      return status;
    }
    
    /* allocate memory */
    ck->x_thermal = calloc (ck->n_wvl+1, sizeof (double *));
    for (i=0; i<=ck->n_wvl; i++)
      ((double **) ck->x_thermal)[i] = calloc (LOWTRAN_MAXINT+1, sizeof(double));
    
    for (i=0; i<n_tmp; i++) {
      iband = (int) (tmp1[i]+0.5); 
      iwght = (int) (tmp2[i]+0.5); 
      
      ((double **) ck->x_thermal)[iband][iwght] = tmp3[i];
    }

    free (tmp1); free (tmp2); free (tmp3);
  }

  return 0;
}
  
  
/*****************************************************************/
/* Calculate cross section for a given temperature, pressure,    */
/* and concentration according to Kato et al. [1999]             */
/*                                                               */
/* output:                                                       */
/*  double **extinction    extinktion coefficient                */
/*  double *weight         k-distr. quadrature weight            */
/*  int *ngauss            number of quadrature points           */
/*                                                               */
/*****************************************************************/

static int kato_crs (float ***temp_avg, float ***press, float ****dens, 
		     float ****dens_avg, float *z, int nlev, int Nx, int Ny,  
                     ck_struct *crs, int iband, int ispcs, char *path, 
		     int quiet, int verbose, 
                     double ****extinction, double *weight, int *ngauss)
{
  
  /* Local variables */
  int lc=0, ix=0, iy=0;
  int status=0;   

  if (crs->comp[ispcs][iband] < 0) {
    /* if correlated k in this band */
    *ngauss = crs->wght[ispcs][iband];
    
    if (ispcs == KATO_H2O) {
      if (iband==9 || iband==10 || iband==11 || iband==13 || iband == 15) {
        status = x_section (crs->scheme, ispcs, press, temp_avg, nlev, Nx, Ny, iband, 
                            path, ngauss, extinction, weight, quiet, verbose);
        if (status!=0) {
          fprintf (stderr, "Error %d returned by x_section\n", status);
          fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
          return status;
        }
      }
      else {
        status = x_section_h2o (crs->scheme, press, temp_avg, dens, dens_avg, nlev, Nx, Ny, 
				iband, path, ngauss, extinction, weight, quiet, verbose);
        if (status!=0) {
          fprintf (stderr, "Error %d returned by x_section_h2o\n", status);
          fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
          return status;
        }
      }
    }
    else {
      status = x_section (crs->scheme, ispcs, press, temp_avg, nlev, Nx, Ny, iband, 
                          path, ngauss, extinction, weight, quiet, verbose);
      if (status!=0) {
        fprintf (stderr, "Error %d returned by x_section\n", status);
        fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
        return status;
      }
    }
  }
  else {
    *ngauss = 1;
    weight[0]     = crs->wght[ispcs][iband];
    for (ix=0; ix<Nx; ix++){
      for (iy=0; iy<Ny; iy++){
	for (lc=0; lc<nlev-1; lc++){
	  extinction[ix][iy][lc][0] = crs->comp[ispcs][iband]; //3DAbs ??
	}
      }
    }
  }
  
  
  return status;
}




/*****************************************************************/
/* Calculate cross section for a given temperature, pressure,    */
/* and concentration                                             */
/*****************************************************************/

static int fu_profile (float ***temper, float ***press, float *z, 
		       float ****dens, float ****dens_avg, int nlev,
		       int Nx, int Ny,
                       int h2o_cont, float umco2, float umch4, float umn2o,
                       float umf11, float umf12, float umf22,
                       ck_profile *ck)
{

  /* float *rhoair=NULL, *rhoo3=NULL, *rhoh2o=NULL; */
  /* float *pres=NULL, *temp=NULL, *davg=NULL; */

  int  status=0;  //3DAbs taken out lc=0


  /* /\* allocate memory for density profiles *\/ */
  /* rhoair = calloc (nlev, sizeof(float)); */
  /* rhoh2o = calloc (nlev, sizeof(float)); */
  /* rhoo3  = calloc (nlev, sizeof(float)); */
  /* temp   = calloc (nlev, sizeof(float)); */
  /* pres   = calloc (nlev, sizeof(float)); */
  /* davg   = calloc (nlev, sizeof(float)); */

  
  /* for (ix=0; ix<Nx; ix++){ */
  /*   for (iy=0, iy<Ny, iy++){ */
  /*     for (lc=0; lc<nlev; lc++) { */
  /* 	rhoair[lc] = dens[MOL_AIR][ix][iy][lc] * M_AIR / N_MOL_CKDFU * 1.0e6; */
  /* 	rhoh2o[lc] = dens[MOL_H2O][ix][iy][lc] * M_H2O / N_MOL_CKDFU * 1.0e6; */
  /* 	rhoo3 [lc] = dens[MOL_O3][ix][iy][lc] * M_O3  / N_MOL_CKDFU * 1.0e6; */
  /* 	temp  [lc] = temper[ix][iy][lc]; */
  /* 	pres  [lc] = press[ix][iy][lc]; */
  /*     } */

  /* for (lc=0; lc<nlev-1; lc++) */
  /*   davg  [lc] = dens_avg[MOL_AIR][ix][iy][lc] * M_AIR  /  */
  /*     N_MOL_CKDFU * 1.0e6 * (z[lc]-z[lc+1]); */

  /* call Fu and Liou routines */

#if HAVE_FULIOU
  status = ckdfu (umco2, umch4, umn2o,
                  umf11, umf12, umf22,
                  //rhoair, rhoh2o, rhoo3, pres, temp, davg, z, nlev,
		  temper, press, z, dens, dens_avg, nlev,
		  Nx, Ny, h2o_cont, ck);
#else
  fprintf (stderr, "Error, Fu and Liou not supported!\n");
  return -1;
#endif

  if (status<0) {
    fprintf (stderr, "Error %d returned by ckdfu()\n", status);
    return status;
  }
 
  /* free (rhoair); */
  /* free (rhoh2o); */
  /* free (rhoo3); */
  /* free (temp); */
  /* free (pres); */
  /* free (davg); */

  return 0; /* ok */
}



/*****************************************************************/
/* C wrapper for David Kratz' AVHRR parameterization; should be  */
/* moved to a separate source file.                              */
/*****************************************************************/

static int avhrr_kratz (int channel, int interval,
                        float *z0, float *p0, float *t0, float *u0, float *ux, 
                        int nlev, float co2, float n2o, 
                        float f11, float f12, float ch4, 
                        double **od, double *wght, int *ntau)
{
  int i=0, j=0;

  /* allocate memory for output */
  float *tmp_od   = calloc (KRATZ_MAXNTAU*(nlev-1), sizeof(float));
  float *tmp_wght = calloc (KRATZ_MAXNTAU, sizeof(float));
  float **tmp_od2 = NULL;


#if HAVE_AVHRR
  switch (channel) {
  case 1:
    switch(interval) {
    case 1:
      F77_FUNC (avhrr11, AVHRR11) (z0, p0, t0, u0, ux, &nlev, &co2, &n2o, &f11, &f12, &ch4, 
               tmp_od, tmp_wght, ntau);
      break;
    case 2:
      F77_FUNC (avhrr12, AVHRR12) (z0, p0, t0, u0, ux, &nlev, &co2, &n2o, &f11, &f12, &ch4, 
                         tmp_od, tmp_wght, ntau);
      break;
    case 3:
      F77_FUNC (avhrr13, AVHRR13) (z0, p0, t0, u0, ux, &nlev, &co2, &n2o, &f11, &f12, &ch4, 
                         tmp_od, tmp_wght, ntau);
      break;
    case 4:
      F77_FUNC (avhrr14, AVHRR14) (z0, p0, t0, u0, ux, &nlev, &co2, &n2o, &f11, &f12, &ch4, 
                         tmp_od, tmp_wght, ntau);
      break;
    case 5:
      F77_FUNC (avhrr15, AVHRR15) (z0, p0, t0, u0, ux, &nlev, &co2, &n2o, &f11, &f12, &ch4, 
                         tmp_od, tmp_wght, ntau);
      break;
    default:
      fprintf (stderr, "error, AVHRR channel %d has no interval %d\n", 
               channel, interval);
      return -1;
    }
    break;
    
  case 2:
    switch(interval) {
    case 1:
      F77_FUNC (avhrr21, AVHRR21) (z0, p0, t0, u0, ux, &nlev, &co2, &n2o, &f11, &f12, &ch4, 
                         tmp_od, tmp_wght, ntau);
      break;
    case 2:
      F77_FUNC (avhrr22, AVHRR22) (z0, p0, t0, u0, ux, &nlev, &co2, &n2o, &f11, &f12, &ch4, 
                         tmp_od, tmp_wght, ntau);
      break;
    case 3:
      F77_FUNC (avhrr23, AVHRR23) (z0, p0, t0, u0, ux, &nlev, &co2, &n2o, &f11, &f12, &ch4, 
                         tmp_od, tmp_wght, ntau);
      break;
    case 4:
      F77_FUNC (avhrr24, AVHRR24) (z0, p0, t0, u0, ux, &nlev, &co2, &n2o, &f11, &f12, &ch4, 
                         tmp_od, tmp_wght, ntau);
      break;
    default:
      fprintf (stderr, "error, AVHRR channel %d has no interval %d\n", 
               channel, interval);
      return -1;
    }
    break;
    
  case 3:
    switch(interval) {
    case 1:
      F77_FUNC (avhrr31, AVHRR31) (z0, p0, t0, u0, ux, &nlev, &co2, &n2o, &f11, &f12, &ch4, 
                         tmp_od, tmp_wght, ntau);
      break;
    case 2:
      F77_FUNC (avhrr32, AVHRR32) (z0, p0, t0, u0, ux, &nlev, &co2, &n2o, &f11, &f12, &ch4, 
                         tmp_od, tmp_wght, ntau);
      break;
    case 3:
      F77_FUNC (avhrr33, AVHRR33) (z0, p0, t0, u0, ux, &nlev, &co2, &n2o, &f11, &f12, &ch4, 
                         tmp_od, tmp_wght, ntau);
      break;
    case 4:
      F77_FUNC (avhrr34, AVHRR34) (z0, p0, t0, u0, ux, &nlev, &co2, &n2o, &f11, &f12, &ch4, 
                         tmp_od, tmp_wght, ntau);
      break;
    case 5:
      F77_FUNC (avhrr35, AVHRR35) (z0, p0, t0, u0, ux, &nlev, &co2, &n2o, &f11, &f12, &ch4, 
                         tmp_od, tmp_wght, ntau);
      break;
    default:
      fprintf (stderr, "error, AVHRR channel %d has no interval %d\n", 
               channel, interval);
      return -1;
    }
    break;
    
  case 4:
    switch(interval) {
    case 1:
      F77_FUNC (avhrr41, AVHRR41) (z0, p0, t0, u0, ux, &nlev, &co2, &n2o, &f11, &f12, &ch4, 
                         tmp_od, tmp_wght, ntau);
      break;
    default:
      fprintf (stderr, "error, AVHRR channel %d has no interval %d\n", 
               channel, interval);
      return -1;
    }
    break;

  case 5:
    switch(interval) {
    case 1:
      F77_FUNC (avhrr51, AVHRR51) (z0, p0, t0, u0, ux, &nlev, &co2, &n2o, &f11, &f12, &ch4, 
                         tmp_od, tmp_wght, ntau);
      break;
    default:
      fprintf (stderr, "error, AVHRR channel %d has no interval %d\n", 
               channel, interval);
      return -1;
    }
    break;
  default:
    fprintf (stderr, "error, unkown AVHRR channel %d\n", channel);
    return -1;
  }
#else
  fprintf (stderr, "Error, AVHRR parameterization not supported!\n");
  return -1;
#endif

  /* convert fortran array to a 2D C array */
  /* tmp_od2 [0 ... nlyr-1][0 ... ntau-1]  */

  tmp_od2 = fortran2c_2D_float_ary(nlev-1, *ntau, tmp_od);
  
  /* copy arrays to their final destination */
    
  for (j=0; j<*ntau; j++) {
    wght[j] = tmp_wght[j];
    for (i=0; i<nlev-1; i++)   /* flip profile */
      od[nlev-2-i][j] = tmp_od2[i][j] / (z0[i+1]-z0[i]);
  }

  free(tmp_od);
  free(tmp_wght);

  for (i=0; i<nlev-1; i++)
    free(tmp_od2[i]);
  free(tmp_od2);

  return 0;
}




/*****************************************************************/
/* Calculate cross section for a given temperature, pressure,    */
/* and concentration; David Kratz AVHRR parameterization         */
/*****************************************************************/

static int avhrr_kratz_profile (float ***temper, float ***press, float *z, 
                                float ****dens, int nlev, 
                                float co2, float ch4, float n2o, 
                                float f11, float f12,
                                ck_profile *ck)
{

  float *z0=NULL, *p0=NULL, *t0=NULL, *u0=NULL, *ux=NULL;

  int lc=0;

  /* allocate memory for profiles */
  z0 = calloc (nlev, sizeof(float));
  p0 = calloc (nlev, sizeof(float));
  t0 = calloc (nlev, sizeof(float));
  u0 = calloc (nlev, sizeof(float));
  ux = calloc (nlev, sizeof(float));


  //3DAbs need this for 3D??, so far only 1D implemented
  /* calculate input profiles */
  for (lc=0;lc<nlev;lc++) {

    z0[nlev-lc-1] = z[lc];

    p0[nlev-lc-1] = press [0][0][lc]/1013.25;   /* convert pressure from mbar to atm */
    t0[nlev-lc-1] = temper[0][0][lc];

    u0[nlev-lc-1] = dens[MOL_H2O][0][0][lc] * H2O_MOLWGHT / NAVOGADRO * 1.0E5;
    ux[nlev-lc-1] = dens[MOL_O3] [0][0][lc] * O3_MOLWGHT  / NAVOGADRO * 1.0E5;

  } 

  /* call AVHRR function */

  /* channel 1 */

  avhrr_kratz(1, 5, z0, p0, t0, u0, ux, nlev, co2, n2o, f11, f12, ch4, 
              ck[0].crs[0][0], ck[0].weight, &(ck[0].ngauss));
  avhrr_kratz(1, 4, z0, p0, t0, u0, ux, nlev, co2, n2o, f11, f12, ch4, 
              ck[1].crs[0][0], ck[1].weight, &(ck[1].ngauss));
  avhrr_kratz(1, 3, z0, p0, t0, u0, ux, nlev, co2, n2o, f11, f12, ch4, 
              ck[2].crs[0][0], ck[2].weight, &(ck[2].ngauss));
  avhrr_kratz(1, 2, z0, p0, t0, u0, ux, nlev, co2, n2o, f11, f12, ch4, 
              ck[3].crs[0][0], ck[3].weight, &(ck[3].ngauss));
  avhrr_kratz(1, 1, z0, p0, t0, u0, ux, nlev, co2, n2o, f11, f12, ch4, 
              ck[4].crs[0][0], ck[4].weight, &(ck[4].ngauss));

  /* channel 2 */
  avhrr_kratz(2, 4, z0, p0, t0, u0, ux, nlev, co2, n2o, f11, f12, ch4,
              ck[5].crs[0][0], ck[5].weight, &(ck[5].ngauss));
  avhrr_kratz(2, 3, z0, p0, t0, u0, ux, nlev, co2, n2o, f11, f12, ch4,
              ck[6].crs[0][0], ck[6].weight, &(ck[6].ngauss));
  avhrr_kratz(2, 2, z0, p0, t0, u0, ux, nlev, co2, n2o, f11, f12, ch4,
              ck[7].crs[0][0], ck[7].weight, &(ck[7].ngauss));
  avhrr_kratz(2, 1, z0, p0, t0, u0, ux, nlev, co2, n2o, f11, f12, ch4,
              ck[8].crs[0][0], ck[8].weight, &(ck[8].ngauss));

  /* channel 3 */
  avhrr_kratz(3, 5, z0, p0, t0, u0, ux, nlev, co2, n2o, f11, f12, ch4,
              ck[9].crs[0][0], ck[9].weight, &(ck[9].ngauss));
  avhrr_kratz(3, 4, z0, p0, t0, u0, ux, nlev, co2, n2o, f11, f12, ch4,
              ck[10].crs[0][0], ck[10].weight, &(ck[10].ngauss));
  avhrr_kratz(3, 3, z0, p0, t0, u0, ux, nlev, co2, n2o, f11, f12, ch4,
              ck[11].crs[0][0], ck[11].weight, &(ck[11].ngauss));
  avhrr_kratz(3, 2, z0, p0, t0, u0, ux, nlev, co2, n2o, f11, f12, ch4,
              ck[12].crs[0][0], ck[12].weight, &(ck[12].ngauss));
  avhrr_kratz(3, 1, z0, p0, t0, u0, ux, nlev, co2, n2o, f11, f12, ch4,
              ck[13].crs[0][0], ck[13].weight, &(ck[13].ngauss));

  /* channel 4 */
  avhrr_kratz(4, 1, z0, p0, t0, u0, ux, nlev, co2, n2o, f11, f12, ch4,
              ck[14].crs[0][0], ck[14].weight, &(ck[14].ngauss));

  /* channel 5 */
  avhrr_kratz(5, 1, z0, p0, t0, u0, ux, nlev, co2, n2o, f11, f12, ch4,
              ck[15].crs[0][0], ck[15].weight, &(ck[15].ngauss));


  return 0; /* ok */
}




int sbdart_profile (float ***temper, float ***press, float *z, 
                    float ***h2o, float ***o3, 
                    float xn2, float xo2, float xco2,
                    float xch4, float xn2o, float xno2, 
		    int ck_o4abs, int ck_n2abs, int ck_coabs, int ck_so2abs, int ck_nh3abs, int ck_noabs, int ck_hno3abs,
                    float *dt, float *ssa, int nlev, 
                    float sza, float lambda, int iv, ck_profile *ck, 
                    int spectral_is)
{

  //3DAbs SBDART/LOWTRAN not yet included for 3D gas absorption

#if HAVE_LOWTRAN

  int lc=0, ig=0;
  float wl = lambda / 1000.0;   /* convert from nm to micron */
  
  /* sbtaugas parameters */
  #define SBTAUGAS_NK 3
  #include "sbtaugas.h"
  
  void F77_FUNC (sbtaugas, SBTAUGAS) (int *nlev, float *zi, float *pi, float *ti,
                           float *wh, float *wo, float *sc,
                           float *wl, float *sza,
                           float *xn2, float *xo2, float *xco2, 
                           float *xch4, float *xn2o, float *xno2, 
			   int *o4abs, int *n2abs, int *coabs, int *so2abs, int *nh3abs, int *noabs, int *hno3abs,
                           float *dtau, float *wght, int *npass);


  float *zi = calloc (nlev, sizeof (float));   /* inverted altitude profile    */
  float *ti = calloc (nlev, sizeof (float));   /* inverted temperature profile */
  float *pi = calloc (nlev, sizeof (float));   /* inverted pressure profile    */
            
  float *wh = calloc (nlev, sizeof (float));   /* water vapour concentration   */
  float *wo = calloc (nlev, sizeof (float));   /* ozone concentration          */

  float *sc = calloc (nlev, sizeof (float));   /* scattering optical depth     */
  
  float dtau[MXLY*SBTAUGAS_NK], wght[SBTAUGAS_NK];
  int npass=0;

  float **dtauc=NULL;

  if (nlev>MXLY) {
    fprintf (stderr, "Error, number of layer too large for sbdart. Increase MXLY\n");
    fprintf (stderr, "       in %s to at least %d and recompile!\n", "src/sbtaugas.param", nlev);
    return -1;
  }

  /* first check if we have enough memory allocated; */
  /* if not, reallocate ck_profile                   */
  if (nlev > ck[iv].nlev) {

    for (lc=0; lc<ck[iv].nlev-1; lc++)
      free (ck[iv].crs[0][0][lc]);
    free (ck[iv].crs[0][0]);
    
    ck[iv].crs[0][0] = calloc (nlev, sizeof(double *));
    for (lc=0; lc<nlev; lc++)
      ck[iv].crs[0][0][lc] = calloc (LOWTRAN_MAXINT, sizeof(double));

    ck[iv].nlev = nlev;
  }


  for (lc=0; lc<nlev; lc++) {
    zi[lc] = z      [nlev-lc-1];
    pi[lc] = press  [0][0][nlev-lc-1];
    ti[lc] = temper [0][0][nlev-lc-1];
    wh[lc] = h2o    [0][0][nlev-lc-1] * 1.0e6 * (M_H2O / N_MOL_CKDFU);
    wo[lc] = o3     [0][0][nlev-lc-1] * 1.0e6 * (M_O3  / N_MOL_CKDFU);
  }


  /* calculate scattering optical thickness; */
  /* optical thickness is defined per layer, */
  /* not per level as the above quantities   */
  for (lc=0; lc<nlev-1; lc++)
    sc[lc] = dt[nlev-lc-2] * ssa[nlev-lc-2];

  F77_FUNC (sbtaugas, SBTAUGAS) (&nlev, zi, pi, ti, wh, wo, sc, &wl, &sza, 
                      &xn2, &xo2, &xco2, 
                      &xch4, &xn2o, &xno2, 
		      &ck_o4abs, &ck_n2abs, &ck_coabs, &ck_so2abs, &ck_nh3abs, &ck_noabs, &ck_hno3abs,
                      dtau, wght, &npass);


  dtauc = fortran2c_2D_float_ary (SBTAUGAS_NK, MXLY, dtau);

  /* copy data to final destination */
  ck[iv].ngauss = npass;
  
  /* CE: For spectral importance sampling, each wavelength is*/
  /*   calculated 3 (LOWTRAN_MAXINT) */
  /*     times, if the wavelength requires only one band, each */
  /*     calculation gets the weight 1/3.  */
  if (spectral_is && npass==1)
    for (ig=0; ig<LOWTRAN_MAXINT; ig++){
      ck[iv].weight[ig] = 1.0/LOWTRAN_MAXINT;
      for (lc=0; lc<nlev-1; lc++)
        ck[iv].crs[0][0][lc][ig] = dtauc[0][lc+1] / (z[lc]-z[lc+1]);
    } 
  else{
    for (ig=0; ig<npass; ig++) {
      
      /* copy data to final destination */
      ck[iv].weight[ig] = wght[ig];
      
      /* ??? ignoring optical thickness in 1st (uppermost) level */
      /* ??? what does that mean anyway ?                        */
      for (lc=0; lc<nlev-1; lc++)
        ck[iv].crs[0][0][lc][ig] = dtauc[ig][lc+1] / (z[lc]-z[lc+1]);
    }
  }
  
  free (zi); free (pi); free (ti); free (wh); free (wo); free (sc);
  ASCII_free_float (dtauc, SBTAUGAS_NK);

  return 0;
#else
  fprintf (stderr, "Error, LOWTRAN/SBDART gas absorption not available.\n");
  return -1;
#endif
}
  



/*****************************************************************/
/* Read cross section profile                                    */
/*****************************************************************/

static int generic_profile (char *filename, float *z, int nlev, int *nwght,
                            ck_profile *ck)
{

#if HAVE_LIBNETCDF

  int status=0;
  int iv=0, ig=0, lu=0;
  double deltaz=0;

  int ncid=0;
  int idd_nwvn=0, idd_nlyr=0, idd_maxc=0;
  int id_nc=0, id_wght=0, id_k=0;

  size_t dimlen=0;

  size_t start[1] = {0};
  size_t count[1] = {0};

  size_t wghtstart[2] = {0,0};
  size_t wghtcount[2] = {0,0};

  size_t crsstart[3] = {0,0,0};
  size_t crscount[3] = {0,0,0};

  int nwvn=0, nlyr=0;

  int *nc=NULL;

  /* open netcdf file */
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d opening netCDF file %s\n", status, filename);
    return status;
  }


  /* get dimension id for "nwvn" */
  status = nc_inq_dimid (ncid, "nwvn", &idd_nwvn);
  
  /* get dimension length for "nwvn" */
  status = nc_inq_dimlen (ncid, idd_nwvn, &dimlen);
  nwvn = dimlen;



  /* get dimension id for "nlyr" */
  status = nc_inq_dimid (ncid, "nlyr", &idd_nlyr);
  
  /* get dimension length for "nlyr" */
  status = nc_inq_dimlen (ncid, idd_nlyr, &dimlen);
  nlyr = dimlen;

  if (nlyr != nlev-1) {
    fprintf (stderr, "Error, different number of layers in atmospheric profile (%d)\n", nlev-1);
    fprintf (stderr, "and absorption coefficient profile %s (%d)\n", filename, nlyr);
    return -1;
  }

  /* get dimension id for "maxc" */
  status = nc_inq_dimid (ncid, "maxc", &idd_maxc);
  
  /* get dimension length for "maxc" */
  status = nc_inq_dimlen (ncid, idd_maxc, &dimlen);

  count[0] = nwvn;
  
  /* allocate memory for fields */
  nc    = calloc (nwvn, sizeof(int));


  /* get variable id for "nc" */
  status = nc_inq_varid (ncid, "nc", &id_nc);
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }

  /* read "nc" */
  status = nc_get_vara_int (ncid, id_nc, start, count, nc);
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }


  /* get variable id for "weight" */
  status = nc_inq_varid (ncid, "weight", &id_wght);
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }

  /* get variable id for "k" */
  status = nc_inq_varid (ncid, "k", &id_k);
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }


  /* Loop over bands */
  for (iv=0; iv<nwvn; iv++) {
    
    /* copy data to final destination */
    ck[iv].ngauss = nc[iv];
    
    wghtstart[0] = iv;
    wghtstart[1] = 0;

    wghtcount[0] = 1;
    wghtcount[1] = nwght[iv+1];
    
    status = nc_get_vara_double (ncid, id_wght, wghtstart, wghtcount, ck[iv].weight);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading %s\n", status, filename);
      return status;
    }
    

    for (lu=0; lu<nlev-1; lu++) {

      deltaz = (z[lu] - z[lu+1]);

      crsstart[0]=iv;
      crsstart[1]=lu;
      crsstart[2]=0;
      
      crscount[0]=1;
      crscount[1]=1;
      crscount[2]=nwght[iv+1];
      
      status = nc_get_vara_double (ncid, id_k, crsstart, crscount, ck[iv].crs[0][0][lu]);
      if (status!=NC_NOERR) {
        fprintf (stderr, "Error %d reading %s\n", status, filename);
        return status;
      }
      
      for (ig=0; ig<nwght[iv+1]; ig++)
        ck[iv].crs[0][0][lu][ig] /= deltaz; //3DAbs ?? where needed
    }
  }



  free (nc);

  nc_close(ncid);

#else
  fprintf (stderr, " ***********************************************************************\n");
  fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
  fprintf (stderr, " * use the generic correlated_k option. Please get netcdf and rebuild. *\n");
  fprintf (stderr, " ***********************************************************************\n");
  return -1;
#endif
  
  return 0;  /* ok */
}




/*************************************************************************/
/* Setup correlated-k structures ck and crs_ck for given profiles.       */
/* ck holds the general information (wavelengths, number of weighting    */
/* points, ...) while crs_ck holds the detailed cross section            */
/* profiles for each component and wavelength.                           */
/*************************************************************************/

int crs_ck (ck_struct *ck,
            float *lambda, int n_lambda, 
            float ***temp, float ***temp_avg, float ***press, float *z, 
            float ****dens, float ****dens_avg, int nlev, int Nx, int Ny, 
            float *mixing_ratio, int h2ocont,
	    int ck_o4abs, int ck_n2abs, int ck_coabs, int ck_so2abs, 
	    int ck_nh3abs, int ck_noabs, int ck_hno3abs,
	    float *sza,
            char *path, char *filename, int first, int verbose,
            int nipa,
            crs_ck_out_struct *crs_ck,
            input_struct input, output_struct *output)
{
  
  int is=0, lc=0, iv=0, ipa=0, status=0, ix=0, iy=0;
  int n_components=-666, n_maxint=-666;
  
  float *dtauc_dummy = calloc (nlev, sizeof(float));
  float *ssalb_dummy = calloc (nlev, sizeof(float));

  /* Allocate memory for crs_ck */
  if (first) {
    switch(ck->scheme) {
    case CK_KATO:
    case CK_KATO2:
    case CK_KATO2_96:
    case CK_KATO2ANDWANDJI:
      /* absorption by O3, O2, CO2, and H2O */
      n_components = KATO_COMPONENTS;
      n_maxint     = KATO_MAXINT;
      break;
    case CK_FU:
      /* need only one absorption cross section; band overlap is handled internally */
      n_components = FU_COMPONENTS;
      n_maxint     = FU_MAXINT;
      break;
    case CK_AVHRR_KRATZ:
      /* need only one absorption cross section; band overlap is handled internally */
      n_components = AVHRR_KRATZ_COMPONENTS;
      n_maxint     = AVHRR_KRATZ_MAXINT;
      break;
    case CK_LOWTRAN:
      /* need only one absorption cross section; band overlap is handled internally */
      n_components = LOWTRAN_COMPONENTS;
      n_maxint     = LOWTRAN_MAXINT;
      break;
    case CK_FILE:
      /* need only one absorption cross section; band overlap is handled internally */
      n_components = GENERIC_COMPONENTS;
      break;
    default:
      fprintf (stderr, "Error: unsupported correlated-k scheme %d\n", ck->scheme);
      return -1;
      break;
    }

    crs_ck->profile = calloc(n_components, sizeof(ck_profile *));
    for (is=0; is<n_components; is++) 
      crs_ck->profile[is] = calloc (n_lambda, sizeof(ck_profile));

    switch(ck->scheme) {
    case CK_KATO:
    case CK_KATO2:
    case CK_KATO2_96:
    case CK_KATO2ANDWANDJI:
    case CK_FU:
    case CK_AVHRR_KRATZ:
    case CK_LOWTRAN:
      for (is=0; is<n_components; is++) {
        for (iv=0; iv<n_lambda; iv++) {
	  crs_ck->profile[is][iv].weight = calloc (n_maxint, sizeof(double));
	  crs_ck->profile[is][iv].crs    = calloc (Nx, sizeof(double ***));
	  for (ix=0; ix<Nx; ix++){
	    crs_ck->profile[is][iv].crs[ix]    = calloc (Ny, sizeof(double **));
	    for (iy=0; iy<Ny; iy++){
	      crs_ck->profile[is][iv].crs[ix][iy]    = calloc (nlev, sizeof(double *));
	      for (lc=0; lc<nlev; lc++)
		crs_ck->profile[is][iv].crs[ix][iy][lc] = calloc(n_maxint, sizeof(double));
	    }
	  }
	}
      }
      break;
    case CK_FILE:
      for (is=0; is<GENERIC_COMPONENTS; is++) {
        for (iv=0; iv<n_lambda; iv++) {
          crs_ck->profile[is][iv].crs    = calloc (Nx, sizeof(double *));
          crs_ck->profile[is][iv].weight = calloc (ck->wght[GENERIC_WGHT][iv+1], sizeof(double));
	  for (ix=0; ix<Nx; ix++){
	    crs_ck->profile[is][iv].crs[ix]    = calloc (Ny, sizeof(double **));
	    for (iy=0; iy<Ny; iy++){
	      crs_ck->profile[is][iv].crs[ix][iy]    = calloc (nlev, sizeof(double *));
	      for (lc=0; lc<nlev; lc++)
		crs_ck->profile[is][iv].crs[ix][iy][lc] = calloc (ck->wght[GENERIC_WGHT][iv+1], sizeof(double));
	    }
	  }
	}
      }
      break;
    default:
      fprintf (stderr, "Error: unsupported correlated-k scheme %d\n", ck->scheme);
      return -1;
      break;
    }

  } /* End of: if (first) {... */
   
  for (ipa=0; ipa<nipa; ipa++) {

    /* /\* additional verbose output of the atmosphere before calculation of abs.-cross-sections / abs-coefficients *\/ */
    /* if (verbose) { */
    /*   fprintf (stderr, "# lc |  z[km]  |  Pressure  | Temp.  |    Air      |   Ozone     |     O2      | Water vap.  |    CO2      |    NO2      |\n"); */
    /*   fprintf (stderr, "#    |         |   [hPa]    |  [K]   |   [cm-3]    |   [cm-3]    |   [cm-3]    |   [cm-3]    |   [cm-3]    |   [cm-3]    |\n"); */
    /*  */
    /*   for (lc=0; lc<nlev; lc++) { */
    /*     if (z[lc] >= output->alt.altitude) { */
    /*       fprintf (stderr, "%5d  %7.3f   %10.5f   %6.2f   %.5e   %.5e   %.5e   %.5e   %.5e   %.5e", */
    /* 	           lc, z[lc], press[lc], temp[lc],  */
    /*                dens[MOL_AIR][lc], dens[MOL_O3][lc],  dens[MOL_O2][lc], dens[MOL_H2O][lc], dens[MOL_CO2][lc], dens[MOL_NO2][lc]); */
    /*       fprintf (stderr, "\n"); */
    /*     } */
    /*   } */
    /* } */

    switch(ck->scheme) {
    case CK_KATO:
    case CK_KATO2:
    case CK_KATO2_96:
    case CK_KATO2ANDWANDJI:
      /* loop over bands */
      for (iv=0; iv<n_lambda; iv++) {
	
        /* CO2 */ 
	status = kato_crs (temp_avg, press, dens, dens_avg, z, nlev, Nx, Ny, 
			   ck, iv+1, KATO_CO2, path, input.quiet, verbose, 
			   crs_ck->profile[KATO_CO2-1][iv].crs,
			   crs_ck->profile[KATO_CO2-1][iv].weight,
			   &(crs_ck->profile[KATO_CO2-1][iv].ngauss));
	if (status!=0) {
          fprintf (stderr, "Error %d calculating absorption cross section of CO2\n", status);
          fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
          return status;
        }
	
        
        /* O3 */
        status = kato_crs (temp_avg, press, dens, dens_avg, z, nlev, Nx, Ny, 
                           ck, iv+1, KATO_O3, path, input.quiet, verbose,
                           crs_ck->profile[KATO_O3-1][iv].crs,
                           crs_ck->profile[KATO_O3-1][iv].weight,
                           &(crs_ck->profile[KATO_O3-1][iv].ngauss));
        if (status!=0) {
          fprintf (stderr, "Error %d calculating absorption cross section of O3\n", status);
          fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
          return status;
        }
        
        
        /* O2 */
        status = kato_crs (temp_avg, press, dens, dens_avg, z, nlev, Nx, Ny, 
                           ck, iv+1, KATO_O2, path, input.quiet, verbose,
                           crs_ck->profile[KATO_O2-1][iv].crs,
                           crs_ck->profile[KATO_O2-1][iv].weight,
                           &(crs_ck->profile[KATO_O2-1][iv].ngauss));
        if (status!=0) {
          fprintf (stderr, "Error %d calculating absorption cross section of O2\n", status);
          fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
          return status;
        }
        
        
        /* H2O */
        status = kato_crs (temp_avg, press, dens, dens_avg, z, nlev, Nx, Ny,
                           ck, iv+1, KATO_H2O, path, input.quiet, verbose,
                           crs_ck->profile[KATO_H2O-1][iv].crs,
                           crs_ck->profile[KATO_H2O-1][iv].weight,
                           &(crs_ck->profile[KATO_H2O-1][iv].ngauss));
        if (status!=0) {
          fprintf (stderr, "Error %d calculating absorption cross section of H2O\n", status);
          fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
          return status;
        }
      }
    break;

    case CK_FU:

      /* calculate absorption coeffitients */    
      status = fu_profile (temp, press, z, dens, dens_avg, nlev, Nx, Ny,
			   h2ocont,
                           mixing_ratio[MX_CO2],
                           mixing_ratio[MX_CH4], 
                           mixing_ratio[MX_N2O], 
                           mixing_ratio[MX_F11], 
                           mixing_ratio[MX_F12], 
                           mixing_ratio[MX_F22], 
                           crs_ck->profile[0]);
    
      if (status!=0) {
        fprintf (stderr, "Error %d returned by fu_profile()\n", status);
        fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
        return status;
      }
    
      break;
      
    case CK_AVHRR_KRATZ:
    
      /* calculate absorption coeffitients */ 
      status = avhrr_kratz_profile (temp, press, z, dens, nlev,
                                    mixing_ratio[MX_CO2],
                                    mixing_ratio[MX_CH4], 
                                    mixing_ratio[MX_N2O], 
                                    mixing_ratio[MX_F11],
                                    mixing_ratio[MX_F12],
                                    crs_ck->profile[0]);
    
      if (status!=0) {
        fprintf (stderr, "Error %d returned by avhrr_kratz_profile()\n", status);
        fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
        return status;
      }
    
      break;
      

    case CK_LOWTRAN:
  
      /* loop over bands */
      for (iv=0; iv<n_lambda; iv++) {
      
        if (first)
          for (is=0; is<LOWTRAN_COMPONENTS; is++)
            crs_ck->profile[is][iv].nlev = nlev;
      
        if (verbose) {
          fprintf (stderr, " ... mixing ratios passed to sbdart_profile():\n");
          fprintf (stderr, " ...  [O2]  = %f\n", mixing_ratio[MX_O2]);
          fprintf (stderr, " ...  [CO2] = %f\n", mixing_ratio[MX_CO2]);
          fprintf (stderr, " ...  [CH4] = %f\n", mixing_ratio[MX_CH4]);
          fprintf (stderr, " ...  [N2O] = %f\n", mixing_ratio[MX_N2O]);
          fprintf (stderr, " ...  [NO2] = %f\n", mixing_ratio[MX_NO2]);
        }

        /* calculate absorption coefficients */
        status = sbdart_profile (temp, press, z, dens[MOL_H2O], dens[MOL_O3], 
                                 -1.0,
                                 mixing_ratio[MX_O2], 
                                 mixing_ratio[MX_CO2], 
                                 mixing_ratio[MX_CH4], 
                                 mixing_ratio[MX_N2O], 
                                 mixing_ratio[MX_NO2], 
				 ck_o4abs, ck_n2abs, ck_coabs, 
				 ck_so2abs, ck_nh3abs, ck_noabs, ck_hno3abs,
                                 dtauc_dummy, ssalb_dummy, nlev, sza[iv], lambda[iv],
                                 iv, crs_ck->profile[0], input.rte.mc.spectral_is);



      
        if (status!=0) {
          fprintf (stderr, "Error %d returned by sbdart_profile()\n", status);
          fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
          return status;
        }
      }

      break;

    case CK_FILE:
    
      status = generic_profile (filename, z, nlev, ck->wght[GENERIC_WGHT], crs_ck->profile[0]);
    
      if (status!=0) {
        fprintf (stderr, "Error %d returned by generic_profile()\n", status);
        fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
        return status;
      }
    
      break;

    default:
      fprintf (stderr, "Error: unsupported correlated-k scheme %d\n", ck->scheme);
      fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
      return -1;
    
      break;
    }

    /* copy crs_ck->profile[is][iv].crs[lc][iq] to crs_ck->ipa_profile[ipa][is][iv].crs[lc][iq] */
  } /* end of the ipa-loop */
  
  free(dtauc_dummy);
  free(ssalb_dummy);

  return 0;
}


/* Rayleigh cross section for the 6 solar bands of the Fu and Liou [1992/93] */
/* parameterization; see Fu and Liou [1993], pg, 2017.                       */

double fu_rayleigh (double pressure, double temperature, int iv)
{
  /* These coefficients don't seem right, re-calculation gives  */
  /* 7.7333e-06 for the first                                   */

  double R[6] = {9.022E-6, 5.282E-7, 5.722E-8, 1.433E-8, 4.526E-9, 1.529E-9};

  
  if (iv>=6)
    return 0;
  else
    return R[iv] * pressure / temperature;
}



/* routines to read the Kato et al. tables; translated from */
/* Seiji Kato's fortran program using f2c                   */

static int x_section (int scheme, int ispcs, float ***press, float ***temp_avg, 
		      int nlev, int Nx, int Ny,
                      int iband, char *path, int *ngauss, double ****xabs, double *weight,
		      int quiet, int verbose)
{
  int status=0;
#if HAVE_LIBNETCDF
  

  int w_id=0, p_id=0, t_id=0, ncid=0, i=0, k=0, lc=0, ix=0, iy=0;
  int wdim_id=0, pdim_id=0, tdim_id=0;
  int T_out_of_range=0;
  int p_out_of_range=0;
  size_t nw=0, np=0, nt=0;

  float *p=NULL, *t=NULL;
  int ipress=0, itemp=0;

  char sband[5]="";
  char xfile[FILENAME_MAX]="";
  double frcnx=0, frcny=0, ofrcnx=0, ofrcny=0; 
  double presslog=0;
  double pressure=0, temperature=0;
  double xsect[80]      /* was [4][20] */;
  double *xsectt=NULL;
  int xabs_id=0;
  static int warn=1;

/* This subroutine is written by Seiji Kato   */
/*          The Pennsylvania State University */
/*          Meteorology Department            */
/*          e-mail kato@essc.psu.edu          */
/*          25 March 1996                     */

/* revised by Sina Lohmann                    */
/* 11 May 2005                                */

/*                                                                                 */
/* This subroutine returns absorption cross section in cm**2 per molecule          */
/* at n-gaussian points at given pressure, temperature.                            */
/* The unit of the cross section is cm^2 per molecule.                             */
/*                                                                                 */
/* Input                                                                           */
/*  ipress                                                                         */
/*  itemp                                                                          */
/*  presslog    log10 of pressure (mb), which you want to obtain cross sections.   */
/*  temp        log10 of temperature (K), which you want to obtain cross sections. */
/*  iband       Band number based on the attached table.                           */
/*  ngauss      Number of the Gaussian quadratures based on the attached table     */
/*                                                                                 */
/* Output                                                                          */
/*  xabs        Absorption cross-sections at the Gaussian quadrature points        */
/*              in cm**2 per molecule.                                             */
/* weight       Gaussian weights                                                   */
/*                                                                                 */
/* This subroutine returns absorption cross section in cm**2 per molecule          */
/* at n-gaussian points at given pressure, temperature, and number concentration   */
/* of species. The species has to be either 'CO2', 'O3', or 'O2'.                  */
/*  ispcs = 1  :  CO2                                                              */
/*  ispcs = 2  :  O3                                                               */
/*  ispcs = 3  :  O2                                                               */
/*  ispcs = 4  :  H2O                                                              */
/* The unit of the cross section is cm**2 per molecule.                            */

  
  /* Function Body */
  /* temperature grid */
  
  strcpy(xfile, path);
  status = add_katopath (xfile, scheme);
  if (status!=0)
    return status;
  strcat(xfile, "cross_section.table.");
  
  sprintf (sband, "%d", iband);
  
  switch (ispcs) {
  case KATO_CO2:
    strcat(xfile, "CO2.");
    break;
  case KATO_O3:
    strcat (xfile, "_O3.");
    break;
  case KATO_O2:
    strcat (xfile, "_O2.");
    break;
  case KATO_H2O:
    strcat (xfile, "H2O.");
    break;
  default:
    fprintf (stderr, "Error, wrong species\n");
    fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
    return -1;
  }
  
  if (iband <= 9)
    strcat(xfile, "0");
  
  strcat(xfile, sband);

  switch (scheme) {
  case CK_KATO:
  case CK_KATO2_96:
    strcat(xfile, ".cdf");
    break;
  case CK_KATO2:
  case CK_KATO2ANDWANDJI:
    strcat(xfile, ".hitran2k.cdf");
    break;
  default:
    fprintf (stderr, "Error, unsupported correlated_k scheme %d\n", scheme);
    return -1;
  }
  
  /* I want error message, but do not wish error to */
  /* be fatal.                                      */
  
  /* open netcdf file */
  status = nc_open (xfile, NC_NOWRITE, &ncid);
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d opening netCDF file %s\n", status, xfile);
    return status;
  }

  /* get dimension id's */
  status = nc_inq_dimid (ncid, "gauss",       &wdim_id);
  status = nc_inq_dimid (ncid, "pressure",    &pdim_id);
  status = nc_inq_dimid (ncid, "temperature", &tdim_id);

  /* get dimensions */
  status = nc_inq_dimlen (ncid, wdim_id, &nw);
  status = nc_inq_dimlen (ncid, pdim_id, &np);
  status = nc_inq_dimlen (ncid, tdim_id, &nt);

  /* get variable id's */
  status = nc_inq_varid (ncid, "weight",      &w_id);
  status = nc_inq_varid (ncid, "pressure",    &p_id);
  status = nc_inq_varid (ncid, "temperature", &t_id);
  status = nc_inq_varid (ncid, "cross_section", &xabs_id);

  p = calloc (np, sizeof(float));
  t = calloc (nt, sizeof(float));
  xsectt = calloc (np*nt*nw, sizeof(double));
  
  /* read pressure */
  status = nc_get_var_float (ncid, p_id, p);
  for (i=0; i<np; i++)
    p[i] = log10(p[i]);
  
  /* read temperature */
  status = nc_get_var_float (ncid, t_id, t);

  /* check if number of weights equals ngauss */
  if (nw!=*ngauss) {
    fprintf (stderr, "Fatal error! Inconsistent number of ck weights in\n");
    fprintf (stderr, "%s and %s/%s\n", xfile, path, "correlated_k/kato/xxx.dat");
    return -1;
  }

  /* read weights */
  status = nc_get_var_double (ncid, w_id, weight);

  /* read cross sections */
  status = nc_get_var_double (ncid, xabs_id, xsectt);  

  nc_close(ncid);

  
  /* loop over layers */
  for (ix=0; ix<Nx; ix++){
    for (iy=0; iy<Ny; iy++){
      for (lc=0; lc<nlev-1; lc++) {
    
	T_out_of_range = 0;
	p_out_of_range = 0;
	
	/* Parameter adjustments */
	--(xabs[ix][iy][lc]);
	
	/* average parameters over layers */
	pressure    = log_average(press[ix][iy][lc], press[ix][iy][lc+1]);         /* [mbar] */   /* double <- float */  
	/* pressure    = 0.5 * ( press[lc] + press[lc+1] ); */     /* [mbar] */   /* double <- float */  /* before 2006-03-06, UH */
	temperature = temp_avg[ix][iy][lc];                                /* [K]    */   /* double <- float */
	
	presslog = log10(pressure);
	
	/* get indices */
	hunt(p,(int) np, &presslog,    &ipress);  /* why not: out_of_range = */
	hunt(t,(int) nt, &temperature, &itemp);   /*          out_of_range = */
	/* UH, 2009 */
	
	/* get cross sections */
	
	for (i=1; i<=*ngauss; ++i) {
	  xsect[(i << 2)-4] = xsectt[(i-1)*np*nt + (itemp-1)*np + (ipress-1)];
	  xsect[(i << 2)-3] = xsectt[(i-1)*np*nt + (itemp-1)*np + ipress];
	  xsect[(i << 2)-1] = xsectt[(i-1)*np*nt + itemp*np + (ipress-1)];
	  xsect[(i << 2)-2] = xsectt[(i-1)*np*nt + itemp*np + ipress];
	}   
	
	/* Do two-dimensional interpolation on pressure and temperature plane. */
        
	frcnx = (presslog - p[ipress - 1]) / (p[ipress] - p[ipress - 1]);
	if (frcnx <= 1)
	  ofrcnx = 1.0 - frcnx;
	else {
	  ofrcnx = 0;
	  p_out_of_range=1;
	}
	
	frcny = (temperature - t[itemp - 1]) / (t[itemp] - t[itemp - 1]);
	if(frcny <= 1)
	  ofrcny = 1.0 - frcny;
	else {
	  ofrcny = 0;
	  T_out_of_range=1;
	}
	
	for (k = 1; k <= *ngauss; ++k) {
	  xabs[ix][iy][lc][k] = ofrcnx *ofrcny * xsect[(k << 2) - 4] + frcnx * 
	    ofrcny * xsect[(k << 2) - 3] + frcnx * frcny * 
	    xsect[(k<< 2) - 2] + ofrcnx * frcny * xsect[(k << 2) - 1];
	}
	
	if (T_out_of_range==1 || p_out_of_range==1) {
	  if (warn==1 && !quiet) { 
	    fprintf (stderr, "\n");
	    fprintf (stderr, "*** Warning, at least one of the model layers is outside the pressure or\n");
	    fprintf (stderr, "*** temperature range covered by the Kato et al. [1999] parameterization,\n");
	    fprintf (stderr, "*** see Fig. 2 and Table 2 of the publication for details. The respective\n");
	    fprintf (stderr, "*** cross sections will be set to zero. This warning typically occurs for\n");
	    fprintf (stderr, "*** profiles that extend high up into regions with very low pressure and/or\n");
	    fprintf (stderr, "*** high temperatures. If this is the case, this warning may be ignored as\n");
	    fprintf (stderr, "*** the molecular densities up there are so small that the absorption\n");
	    fprintf (stderr, "*** can be neglected for most purposes. To get rid of the warning\n");
	    fprintf (stderr, "*** you may cut the atmosphere at a reasonable altitude.\n");
	    fprintf (stderr, "*** data file: %s, x_section\n", xfile);
	    fprintf (stderr, "\n");
	    warn=0;
	  }
	  if (verbose) {
	    fprintf (stderr, "*** Warning, ");
	    if (T_out_of_range==1)
	      fprintf (stderr, "T[%d] = %f K ", lc, temperature);
	    if (T_out_of_range==1 && p_out_of_range==1)
	      fprintf (stderr, "and");
	    if (p_out_of_range==1)
	      fprintf (stderr, "p[%4d] = %f hPa ", lc, pow(10.0, presslog));
	    fprintf (stderr, "out of range, iband=%d, gas=",iband);
	    switch (ispcs) {
	    case KATO_CO2:
	      fprintf (stderr, "%s\n","CO2");
	      break;
	    case KATO_O3:
	      fprintf (stderr, "%s\n","O3");
	      break;
	    case KATO_O2:
	      fprintf (stderr, "%s\n","O2");
	      break;
	    case KATO_H2O:
	      fprintf (stderr, "%s\n","H2O");
	      break;
	    default:
	      fprintf (stderr, "Error, wrong species\n");
	      fprintf (stderr, "      (line %d, function %s in %s)\n", __LINE__, __func__, __FILE__);
	      return -1;
	    }
	  }
	}
	++(xabs[ix][iy][lc]);
      }
    }
  }
  
  /* free memory */
  free(p); free(t); free(xsectt);
  
#endif  
  return status;
}



static int x_section_h2o (int scheme, float ***press, float ***temp_avg, float ****dens, 
			  float ****dens_avg,
                          int nlev, int Nx, int Ny, int iband, char *path, 
                          int *ngauss, double ****xabs, double *weight, 
			  int quiet, int verbose)
{
  /* Initialized data */
  int status=0;
#if HAVE_LIBNETCDF
  
  int tmp=0;

  int w_id=0, p_id=0, t_id=0, r_id=0, ncid=0, i=0, k=0, lc=0, ix=0, iy=0;
  int wdim_id=0, pdim_id=0, tdim_id=0, rdim_id=0;
  int T_out_of_range=0;
  int p_out_of_range=0;
  size_t nw=0, np=0, nt=0, nr=0;

  float *p=NULL, *t=NULL, *r=NULL;
  int ipress=0, itemp=0, irho=0;

  double rholog=0, presslog=0;
  double pressure=0, temperature=0, density=0;

  char sband[5]="";
  char xfile[FILENAME_MAX]="";
  double frcnx=0, frcny=0, ofrcnx=0, ofrcny;
  double xsect[160]     /* was [8][20] */;
  double ansyrl=0, ansyru=0;
  double *xsectt=NULL;
  int xabs_id=0;
  static int warn=1;


  /* This subroutine is written by Seiji Kato   */
  /*          The Pannsylvania State University */
  /*          Meteorology Department            */
  /*          e-mail kato@essc.psu.edu          */
  /*          25 March 1996                     */
  
  /* revised by Sina Lohmann                    */
  /* 11 May 2005                                */


  /* Function Body */

  strcpy(xfile, path);
  status = add_katopath (xfile, scheme);
  if (status!=0)
    return status;
  strcat(xfile, "cross_section.table.H2O.");
  
  if (iband <= 9)
    strcat(xfile, "0");
  
  sprintf (sband, "%d", iband);

  strcat(xfile, sband);

  switch (scheme) {
  case CK_KATO:
  case CK_KATO2_96:
    strcat(xfile, ".cdf");
    break;
  case CK_KATO2:
  case CK_KATO2ANDWANDJI:
    strcat(xfile, ".hitran2k.cdf");
    break;
  default:
    fprintf (stderr, "Error, unsupported correlated_k scheme %d\n", scheme);
    return -1;
  }

  /* open netcdf file */
  status = nc_open (xfile, NC_NOWRITE, &ncid);
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d opening netCDF file %s\n", status, xfile);
    return status;
  }

  /* get dimension id's */
  status = nc_inq_dimid (ncid, "gauss", &wdim_id);
  status = nc_inq_dimid (ncid, "pressure", &pdim_id);
  status = nc_inq_dimid (ncid, "temperature", &tdim_id);
  status = nc_inq_dimid (ncid, "rho", &rdim_id);

  /* get dimensions */
  status = nc_inq_dimlen (ncid, wdim_id, &nw);
  status = nc_inq_dimlen (ncid, pdim_id, &np);
  status = nc_inq_dimlen (ncid, tdim_id, &nt);
  status = nc_inq_dimlen (ncid, rdim_id, &nr);

  /* get variable id's */
  status = nc_inq_varid (ncid, "weight",      &w_id);
  status = nc_inq_varid (ncid, "pressure",    &p_id);
  status = nc_inq_varid (ncid, "temperature", &t_id);
  status = nc_inq_varid (ncid, "rho",         &r_id);

  p = calloc (np, sizeof(float));
  t = calloc (nt, sizeof(float));
  r = calloc (nr, sizeof(float));
  xsectt = calloc (np*nt*nr*nw, sizeof(double));

  /* read pressure */
  status = nc_get_var_float (ncid, p_id, p);
  for (i=0; i<np; i++)
    p[i] = log10(p[i]);

  /* read temperature */
  status = nc_get_var_float (ncid, t_id, t);

  /* read rho */
  status = nc_get_var_float (ncid, r_id, r);
  for (i=0; i<nr; i++)
    r[i] = log10(r[i]);

  /* check if number of weights equals ngauss */
  if (nw!=*ngauss) {
    fprintf (stderr, "Fatal error! Inconsistent number of ck weights in\n");
    fprintf (stderr, "%s and %s/%s\n", xfile, path, "correlated_k/kato/xxx.dat");
    return -1;
  }

  /* read weights */
  status = nc_get_var_double (ncid, w_id, weight);
  
  /* get the id of cross sections */
  status = nc_inq_varid (ncid, "cross_section", &xabs_id);
  
  /* read cross sections */
  status = nc_get_var_double (ncid, xabs_id, xsectt);

  for (ix=0; ix<Nx; ix++){
    for (iy=0; iy<Ny; iy++){
      
      for (lc=0; lc<nlev-1; lc++) {
	
	T_out_of_range = 0;
	p_out_of_range = 0;
    
	/* Parameter adjustments */
	--(xabs[ix][iy][lc]); //3DAbs why is this needed ??
    
	/* average parameters over layers */
	
	pressure    = log_average(press[ix][iy][lc], press[ix][iy][lc+1]);     
	/* [mbar] */   /* double <- float */  
	/* pressure    = 0.5 * ( press[lc] + press[lc+1] ); */     /* [mbar] */ 
	/* double <- float */  /* before 2006-03-06, UH */
	temperature = temp_avg[ix][iy][lc];                                /* [K]    */  
	/* double <- float */
	density = dens_avg[MOL_H2O][ix][iy][lc];
	/* density = 0.5 * (dens[MOL_H2O][lc] + dens[MOL_H2O][lc+1]); */ 
	/* before 2006-03-06, UH */
	
	presslog = log10(pressure);
	rholog   = log10(density*1e6);   /* convert density from cm-3 to m-3 and calculate log */

	/* get indices */
	hunt(p, (int) np, &presslog, &ipress);     /* why not?: out_of_range = */
	hunt(t, (int) nt, &temperature, &itemp);   /*           out_of_range = */
	hunt(r, (int) nr, &rholog,   &irho);       /*           out_of_range = */
	tmp = irho;                                /* UH, 2009 */
	
	for (i=1; i<=*ngauss; i++)  {
	  while ((irho-1)<10 && i==1) {
	    if (xsectt[(i-1)*nt*np*nr+(irho-1)*np*nt+(itemp-1)*np+(ipress-1)] <= 0)
	      ++irho;
	    else
	      break;
	  }
	  xsect[(i << 3) - 8] = xsectt[(i-1)*nt*np*nr+(irho-1)*np*nt+(itemp-1)*np+(ipress-1)];
	}
	
	irho=tmp;
	for (i=1; i<=*ngauss; i++) {
	  while ((irho-1)<10 && i==1) {
	    if (xsectt[(i-1)*nt*np*nr+(irho-1)*np*nt+(itemp-1)*np+ipress] <= 0)
	      ++irho;
	    else
	      break;
	  }
	  xsect[(i << 3) - 7] = xsectt[(i-1)*nt*np*nr+(irho-1)*np*nt+(itemp-1)*np+ipress];
	}
	
	irho=tmp;
	for (i=1; i<=*ngauss; i++) {
	  while ((irho-1)<10 && i==1) {
	    if (xsectt[(i-1)*nt*np*nr+(irho-1)*np*nt+itemp*np+(ipress-1)] <= 0)
	      ++irho;
	    else
	      break;
	  }
	  xsect[(i << 3) - 5] = xsectt[(i-1)*nt*np*nr+(irho-1)*np*nt+itemp*np+(ipress-1)];
	}
	
	irho=tmp;
	for (i=1; i<=*ngauss; i++) {
	  while ((irho-1)<10 && i==1) {
	    if (xsectt[(i-1)*nt*np*nr+(irho-1)*np*nt+itemp*np+ipress] <= 0)
	      ++irho;
	    else
	      break;
	  }
	  xsect[(i << 3) - 6] = xsectt[(i-1)*nt*np*nr+(irho-1)*np*nt+itemp*np+ipress];
	}
	
	irho=tmp;
	if (irho <= 10) {
	  for (i=1; i<=*ngauss; i++) {
	    while ((irho)<10 && i==1) {
	      if (xsectt[(i-1)*nt*np*nr+irho*np*nt+(itemp-1)*np+(ipress-1)] <= 0)
		++irho;
	      else
		break;
	    }
	    xsect[(i << 3) - 4] = xsectt[(i-1)*nt*np*nr+irho*np*nt+(itemp-1)*np+(ipress-1)];
	    
	  }
	  
	  irho=tmp; 
	  for (i=1; i<=*ngauss; i++) {
	    while ((irho)<10 && i==1) {
	      if (xsectt[(i-1)*nt*np*nr+irho*np*nt+(itemp-1)*np+ipress] <= 0)
		++irho;
	      else
		break;
	    }
	    xsect[(i << 3) - 3] = xsectt[(i-1)*nt*np*nr+irho*np*nt+(itemp-1)*np+ipress];
	  }
	  
	  irho=tmp;
	  for (i=1; i<=*ngauss; i++) {
	    while ((irho)<10 && i==1) {
	      if (xsectt[(i-1)*nt*np*nr+irho*np*nt+itemp*np+(ipress-1)] <= 0)
		++irho;
	      else
		break;
	    }
	    xsect[(i << 3) - 1] = xsectt[(i-1)*nt*np*nr+irho*np*nt+itemp*np+(ipress-1)];
	  }
	  
	  irho=tmp;
	  for (i=1; i<=*ngauss; i++) {
	    while ((irho)<10 && i==1) {
	      if (xsectt[(i-1)*nt*np*nr+irho*np*nt+itemp*np+ipress] <= 0)
		++irho;
	      else
		break;
	    }
	    xsect[(i << 3) - 2] = xsectt[(i-1)*nt*np*nr+irho*np*nt+itemp*np+ipress];
	  }
	  
	  irho=tmp;
	  
	  /* Checking interpolation region */
	  frcnx = (presslog - p[ipress - 1]) / 
	    (p[ipress] - p[ipress - 1]);
	  if (frcnx <= 1)
	    ofrcnx = 1.0 - frcnx;
	  else {
	    ofrcnx = 0;
	    p_out_of_range=1; 
	  }
	  
	  frcny = (temperature - t[itemp - 1]) / (t[itemp] - t[itemp - 1]);
	  if(frcny <= 1)
	    ofrcny = 1.0 - frcny;
	  else {
	    ofrcny = 0;
	    T_out_of_range=1;
	  }
	  
	  for (k = 1; k <= *ngauss; ++k) {
	    ansyrl = ofrcnx * ofrcny * xsect[(k << 3) - 8] + 
	      frcnx * ofrcny * xsect[(k << 3) - 7] + frcnx * 
	      frcny * xsect[(k << 3) - 6] + ofrcnx * frcny * 
	      xsect[(k << 3) - 5];
	    
	    ansyru = ofrcnx * ofrcny * xsect[(k << 3) - 4] + 
	      frcnx * ofrcny * xsect[(k << 3) - 3] + frcnx * 
	      frcny * xsect[(k << 3) - 2] + ofrcnx * frcny * 
	      xsect[(k << 3) - 1];
	    
	    xabs[ix][iy][lc][k] = (ansyru - ansyrl) * (rholog - r[irho - 1]) / 
	      (r[irho] - r[irho - 1]) + ansyrl;
	  }
	} 
	else { /* Checking interpolation region */
	  frcnx = (presslog - p[ipress - 1]) / (p[ipress] - p[ipress - 1]);
	  if (frcnx <= 1)
	    ofrcnx = 1.0 - frcnx;
	  else {
	    ofrcnx = 0;
	    p_out_of_range=1; 
	  }
	  
	  frcny = (temperature - t[itemp - 1]) / (t[itemp] - t[itemp - 1]);
	  if (frcny <= 1)
	    ofrcny = 1.0 - frcny;
	  else {
	    ofrcny = 0;
	    T_out_of_range=1;
	  }
	  
	  for (k = 1; k <= *ngauss; ++k) {
	    xabs[ix][iy][lc][k] = ofrcnx * ofrcny * xsect[(k << 3) - 8] + 
	      frcnx * ofrcny * xsect[(k << 3) - 7] + frcnx * 
	      frcny * xsect[(k << 3) - 6] + ofrcnx * frcny * 
	      xsect[(k << 3) - 5];
	  }
	}
	
	if (T_out_of_range==1 || p_out_of_range==1) {
	  if (warn==1) {
	    fprintf (stderr, "\n");
	    fprintf (stderr, "*** Warning, at least one of the model layers is outside the pressure,\n");
	    fprintf (stderr, "*** temperature or concentration range covered by the Kato et al. [1999]\n");
	    fprintf (stderr, "*** parameterization,see Fig. 2 and Table 2 of the publication for details.\n");
	    fprintf (stderr, "*** The respective cross sections will be set to zero. This warning\n");
	    fprintf (stderr, "*** typically occurs for profiles that extend high up into regions with\n");
	    fprintf (stderr, "*** very low pressure and/or high temperatures. If this is the case,\n");
	    fprintf (stderr, "*** this warning may be ignored as the molecular densities up there are\n");
	    fprintf (stderr, "*** so small that the absorption can be neglected for most purposes.\n");
	    fprintf (stderr, "*** To get rid of the warning you may cut the atmosphere at a reasonable\n");
	    fprintf (stderr, "*** altitude.\n");
	    fprintf (stderr, "*** data file: %s, x_section_h2o\n", xfile);
	    fprintf (stderr, "\n");
	    warn=0;
	  }
	  if (verbose) {
	    fprintf (stderr, "*** Warning, ");
	    if (T_out_of_range==1)
	      fprintf (stderr, "T[%4d] = %f K ", lc, temperature);
	    if (T_out_of_range==1 && p_out_of_range==1)
	      fprintf (stderr, " and ");
	    if (p_out_of_range==1)
	      fprintf (stderr, "p[%4d] = %f hPa ", lc, pow(10.0, presslog));
	    fprintf (stderr, "out of range, iband=%d, gas=H2O\n",iband);
	  }
	}
        
	/* Parameter adjustments */
	++(xabs[ix][iy][lc]); //3DAbs ???
      }
    }
  }
  nc_close (ncid);
  
  /* free memory */
  free(p); free(t); free(r); free(xsectt);
#endif
  return status;
}



int hunt(float *xx, int n, double *x, int *jlo)
{
  static int ascnd=0;
  static int jm=0, jhi=0;
  
  /* Parameter adjustments */
  --xx;
  
  /* Function Body */
  ascnd = xx[n] > xx[1];
  
  /* set initial values */      
  *jlo = 0;                             
  jhi = n+1;
  
  /* if out ouf range set values to end of array and return */
  if(ascnd) {
    if (*x >= xx[n]) {
      *jlo=n-1;   /* temperature too high, warning will be given! */
      return 1;
    }
    if (*x <= xx[1]) {
      *jlo=1;  /* temperature too low, warning will be given! */
      return 1;
    }
  }
  else {
    if (*x <= xx[n]) {
      /* this case occurs for low pressures and low concentrations
	 warning will be given! */
      *jlo=n-1;
      return 1;
    }
    if (*x >= xx[1]) {	/* this case occurs for high pressures and high concentrations */
      *jlo=1;
      return 1;
    }	
  }
  
  /* Begin the final bisection phase */
  while (jhi - *jlo != 1) {
    jm = (jhi + *jlo) >> 1;
    
    if ((*x > xx[jm]) == ascnd)
      *jlo = jm;
    else
      jhi = jm;
  }
  return 0;
}

