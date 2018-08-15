/*--------------------------------------------------------------------
 * $Id: phasetable.c 3219 2016-05-20 10:27:37Z Claudia.Emde $
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

/* @45c@ */
/* @c45@ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "phasetable.h"
#include "numeric.h"
#include "miecalc.h"
#include "mystic.h"
#include "errors.h"

#ifndef PI
#define PI 3.14159265358979323846264338327
#endif

/* number of entries in the cumulative phase function lookup table */
#define NPHASE 20000


/* prototypes of internal functions */
static double HG (double g, double mu) ; /* BCA this is from mystic.c */
//static float Rayleigh_depol (double mu, float depol); /* BCA this is from mystic */

static inline double legendre_polynomial (double plm1, double plm2,
					  double mu, int l);

/***********************************************************************************/
/* Function: setup_Legendre_table                                         @45_30i@ */
/* Description:                                                                    */
/*  Calculate lookup-tables of phase function and cumulative phase function        */
/*  for a given wavelength index iv from the single scattering properties          */
/*  stored in ssprop                                                               */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i45_30@ */
/***********************************************************************************/

int setup_Legendre_table (phase_function_table *phase,
			  ssprop_struct *ssprop,
			  int iv,
			  double rmin, double rmax, 
			  double delta_scaling_mucut, double truncate, 
			  int quiet)
{
  int i=0, l=0, ip=0, status=0;
  int imin=0, imax=0;

  double r0=0, r0_old, dr=0;
  size_t nreff=0, nreffact=0, nphamat=0;

  double **mu=NULL;
  int *n=NULL;

  float **p=NULL;
  double **legen=NULL;
  
  if (rmin > rmax) {
    fprintf (stderr, "Error, found rmin > rmax; this message indicates that you told MYSTIC\n");
    fprintf (stderr, "to use a Mie phase function for a cloud of zero optical thickness;\n");
    fprintf (stderr, "Please change either of these two because this causes trouble internally!\n");
    return -1;
  }

  /* get radius intervals */
  r0     = ssprop->r0[iv];
  dr     = ssprop->dr[iv];
  nreff  = ssprop->nreff[iv];
  nphamat = ssprop->nphamat;


  /* check if the range of radii specified in the 3D input file is covered; */
  /* if only one radius is provided and its value is negative, the phase    */
  /* function is assumed to be constant with location.                      */
     
  if (nreff==1 && r0<0) { 
    fprintf (stderr, " ... found only one (negative) effective radius; assuming\n");
    fprintf (stderr, " ... phase function is constant with location.\n");
  }
  else {
    /* chris - numerical problem here with double precision check when an input */
    /* rmin and rmax are FLOAT precision and should not be treated as double!!! */
    if (rmin < r0*(1. - 1.E-7) || rmax > (r0 + (double) (nreff-1) * dr) * (1. + 1.E-7)) {
      fprintf (stderr, "Error, range of effective radii provided is not sufficient:\n");
      fprintf (stderr, " provided: %16.14f - %16.14f, required: %16.14f - %16.14f\n", 
	       r0, r0 + (double) (nreff-1) * dr, rmin, rmax);
      return -1;
    }
  }
  

  /* determine the required range of effective radii */
  imin=0;
  while (r0 + (double) imin * dr < rmin) {
    imin++;
    if (imin==nreff) {
      imin = nreff-1;
      break;
    }
  }

  imin = (imin>0?imin-1:0);
  
  imax=nreff-1;
  while (r0 + (double) imax * dr > rmax) {
    imax--;
    if (imax==-1) {
      imax = 0;
      break;
    }
  }
  
  imax = (imax<nreff-1?imax+1:nreff-1);

  r0_old = r0;
  r0 = r0_old + (double) imin * dr;

  nreffact = (imax-imin+1);

  if (!quiet)
    fprintf (stderr, " ... using %d radii between %g and %g\n", 
	     (int)nreffact, r0_old + (double) imin * dr, r0_old + (double) imax * dr);


  calloc_phase_function_table (phase, nreffact, r0, dr);
  
  for (i=0; i<nreffact; i++) {   /* calculate only those tables which are really needed */
    
    if (!quiet)
      fprintf (stderr, " ... calculating Legendre table for reff = %g\n", 
	       r0_old + (double) (imin+i) * dr);


    switch (ssprop->type) {
    case PHASE_MOMENTS:
     
      /* phase function is derived from legendre coefficients */
      /* first evaluate phase function at a predefined grid   */
      legen=calloc (nphamat, sizeof(double*));

      for (ip=0; ip<nphamat; ip++)
        legen[ip]=calloc(ssprop->nleg[iv][imin+i], sizeof(double) );
     
      for (ip=0; ip<nphamat; ip++)
        for (l=0; l<ssprop->nleg[iv][imin+i]; l++)
          legen[ip][l] = (double) ssprop->legen[iv][imin+i][ip][l];
      
      /* calculate table of phase function p(mu) */
      status = calc_Legendre_phase (legen, ssprop->nleg[iv][imin+i], nphamat, &mu, &p, &n);
      
      if (status!=0) {
	fprintf (stderr, "Error %d returned by calc_Legendre_phase()\n", status);
	return status;
      }
      
      status = calc_cumulative_table (mu, p, n, nphamat, truncate, phase->iphase[i],
                                      delta_scaling_mucut, quiet);
      if (status!=0) {
	fprintf (stderr, "Error %d calculating cumulative probability table\n", status);
	return status;
      }
      
      for(ip=0; ip<nphamat; ip++){
        free(legen[ip]);
        free(mu[ip]);
        free(p[ip]);
      }
      
      free(legen); free(mu); free(p); free(n);
      break;

    case PHASE_EXPLICIT:
    case PHASE_HYBRID:
      /* phase function has been given explicitly    */
      /* calculate table of phase function p(mu)     */

      status = calc_cumulative_table (ssprop->mu[iv][imin+i], ssprop->phase[iv][imin+i],
                                      ssprop->ntheta[iv][imin+i], ssprop->nphamat,
				      truncate, phase->iphase[i], delta_scaling_mucut,
				      quiet);
      if (status!=0) {
	fprintf (stderr, "Error %d calculating cumulative probability table\n", status);
	return status;
      }

      break;
      
    default:
      fprintf (stderr, "Error, unknown phase function specification %d\n", ssprop->type);
      return -1;
    }

    /* copy extinction coefficient to final destination */
    phase->ext[i]   = ssprop->extinc[iv][imin+i];

    /* copy scattering coefficient to final destination */
    phase->sca[i]   = ssprop->extinc[iv][imin+i] * ssprop->albedo[iv][imin+i];

    /* copy delta scaling factor to final destination */
    phase->f[i] = ssprop->f[iv][imin+i];

    /* copy delta scaling factor TZ to final destination */
    phase->dscale[i] = phase->iphase[i]->dscale;

    if (!quiet)
      fprintf (stderr, "DDSSCALE %f\n", phase->dscale[i]);

  }
  
  /* calculate extinction and scattering increments */
  for (i=0; i<nreffact-1; i++) {
    phase->dext   [i] = (phase->ext   [i+1] - phase->ext   [i]) / dr;
    phase->dsca   [i] = (phase->sca   [i+1] - phase->sca   [i]) / dr;
  }

  return 0;
}


/***********************************************************************************/
/* Function: create_iphase_from_HG                                        @45_30i@ */
/* Description:                                                                    */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i45_30@ */
/***********************************************************************************/

int create_iphase_from_HG (pft *iphase, double g1, int quiet)
{
  int it=0, status=0;
  double increment=0.0;

  float **theta=NULL;
  double **mu=NULL;
  int *n=NULL;

  float **p=NULL;

  p=calloc(1, sizeof(float *));
  theta=calloc(1, sizeof(float *));
  mu=calloc(1, sizeof(double *));
  n=calloc(1, sizeof(int));

  if (g1 < 0.99) {
    status = create_phase_from_HG (g1, n, &(theta[0]), &(mu[0]), &(p[0]), 0);
    if (status!=0) {
      fprintf(stderr,"Error in create_phase_from_HG\n");
      return -1;
    }
  }
  else {
    /* need to do highres phase function */
    n[0]=NPHASE;

    p[0]=calloc(n[0], sizeof(float));
    mu[0]=calloc(n[0], sizeof(double));


    increment = PI/(double) (n[0]-1);
  
    for (it=0; it<n[0]; it++) {
      mu[0][it] = -cos((double) it * increment);
      p[0][it] = HG (g1, mu[0][it]);
    }
  }

  status = calc_cumulative_table( mu, p, n, 1, -1.0, iphase, -1.0, quiet);

  if (status!=0) {
    fprintf(stderr,"Error in calc_cumulative_table\n");
    return -1;
  }

  free(p[0]);
  free(p);
  free(mu[0]);
  free(mu);
  free(theta[0]);
  free(theta);
  free(n);

  return 0;
}


static double HG (double g, double mu) 
{
  double temp = 1.0 + g*g - 2.0*g*mu;

  return (1.0-g*g) / (temp * sqrt(temp));
}

/***********************************************************************************/
/* Function: free_pft                                                     @45_30i@ */
/* Description:                                                                    */
/*  Free memory of struct pft; used by free_phase_function_table().                */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i45_30@ */
/***********************************************************************************/

void free_pft (pft *p)
{
  int j=0, ip=0;

  for (j=0; j<p->nscales; j++) {

    /* careful: free only what has actually been allocated; see comment in    */
    /* calloc_iphase(): "full phase matrix not needed for delta scaled stuff" */

    if (j==0) {
      for (ip=0; ip<p->nphamat; ip++) 
	free (p->p[j][ip]);
    }
    else
      free (p->p[j][0]);

    free (p->p[j]);
    free (p->F[j]);
    free (p->A[j]);
    free (p->B[j]);
    free (p->C[j]);
  }

  free (p->p);
  free (p->F);
  free (p->A);
  free (p->B);
  free (p->C);

  for (ip=0; ip<p->nphamat; ip++) 
    free (p->mu[ip]);
  free (p->mu);
  
}


/***********************************************************************************/
/* Function: calloc_phase_function_table                                  @45_30i@ */
/* Description:                                                                    */
/*  Allocate memory for struct phase_function_table                                */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i45_30@ */
/***********************************************************************************/

void calloc_phase_function_table (phase_function_table *p, int n, double r0, double dr)
{
  int i=0;

  p->n  = n;
  p->r0 = r0;
  p->dr = dr;

  p->iphase = calloc ( n, sizeof(pft *));

  p->ext  = calloc ( n, sizeof(double));
  p->dext = calloc ( n, sizeof(double));

  p->sca  = calloc ( n, sizeof(double));
  p->dsca = calloc ( n, sizeof(double));

  p->f   = calloc ( n, sizeof(double));
  p->dscale = calloc ( n, sizeof(double));

  for (i=0; i<n; i++)
    p->iphase[i] = calloc (1, sizeof(pft));
}


/***********************************************************************************/
/* Function: free_phase_function_table                                    @45_30i@ */
/* Description:                                                                    */
/*  Free memory of struct phase_function_table                                     */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i45_30@ */
/***********************************************************************************/

void free_phase_function_table (phase_function_table *p)
{
  int i=0;

  for (i=0; i<p->n; i++)
    free_pft (p->iphase[i]);

  free (p->iphase);

  free (p->ext);
  free (p->dext);

  free (p->sca);
  free (p->dsca);

  free (p->f);
  free (p->dscale);

}


/***********************************************************************************/
/* Function: calc_Legendre_phase                                          @45_30i@ */
/* Description:                                                                    */
/*  Calculate the phase function from Legendre coefficients                        */ 
/*  for a set of polar angles.                                                     */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i45_30@ */
/***********************************************************************************/

int calc_Legendre_phase (double **pmom, int nmom, int nphamat, double ***mu, float ***p, int **n)
{
  int l=0, it=0, ip=0;

  double pl=0, plm1=0, plm2=0;
  double factor=0, sum=0, increment=0;

  int L=nmom;
  double **f=NULL;

  f=calloc(nphamat, sizeof(double*));
  for (ip=0; ip<nphamat; ip++)
    f[ip]=calloc(nmom, sizeof(double) );

  *n=calloc (nphamat, sizeof(int));      
 
  /* multiply by (2l+1) to convert from 'moments of the phase function' */
  /* to Legendre coefficients;                                          */
  for (ip=0; ip<nphamat; ip++)
    for (l=0; l<L; l++)
      f[ip][l] = pmom[ip][l] * (double) (2*l+1);
  
  /* normalize P11 to 1 */
  factor = 1.0 / f[0][0];

  for (ip=0; ip<nphamat; ip++){
    for (l=0; l<L; l++)
      f[ip][l] *= factor;
  
    /* number of entries in the cumulative phase function lookup table */
    /* default: all phase function elements have the same number of entries*/
    (*n)[ip] = NPHASE;
  }
  *mu=calloc (nphamat, sizeof(double*)); 
  *p=calloc (nphamat, sizeof(double*));
  
  for (ip=0; ip<nphamat; ip++){
    (*mu)[ip] = calloc ((*n)[ip], sizeof(double));
    (*p)[ip] = calloc ((*n)[ip], sizeof(double));
  }

  
  for (ip=0; ip<nphamat; ip++){
    increment = PI/(double) ((*n)[ip]-1);
  
    for (it=0; it<(*n)[ip]; it++) {
      (*mu)[ip][it] = -cos((double) it * increment);
      
      plm2 = 1;
      plm1 = (*mu)[ip][it];
      
      sum = plm2*f[ip][0] + plm1*f[ip][1];
      
      for (l=2; l<L; l++) {
        pl = ( (double) (2*l-1) * (*mu)[ip][it] * plm1 - (double) (l-1) * plm2) / (double) l;
        sum += f[ip][l]*pl;
        
        plm2=plm1;
        plm1=pl;
      }
      
      (*p)[ip][it] = sum;
      /* fprintf(stderr, "YYY %d %g %g \n", ip, (*mu)[ip][it], (*p)[ip][it] ); */
      
      if ((*p)[ip][it]<0 && ip==0) {
        fprintf (stderr, "Warning, phase function %g negative at mu = %g!\n", 
                 (*p)[ip][it], (*mu)[ip][it]);
      }
      
    }
    
  }

  
  /* free memory */
  for (ip=0; ip<nphamat; ip++)
    free(f[ip]); 
  free(f);
 
  return 0;
}


/***********************************************************************************/
/* Function: calloc_iphase                                                @62_30i@ */
/* Description:                                                                    */
/*  Allocate iphase                                                                */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int calloc_iphase (pft *iphase, int nscales, int nphamat, int *n)
{
  int ip=0, j=0;

  iphase->nscales=nscales;
  iphase->nphamat=nphamat; 

  iphase->n = calloc ((size_t) nphamat, sizeof(int));
  for (ip=0; ip<nphamat; ip++)
    iphase->n[ip] = n[ip];

  iphase->mu = calloc ((size_t) nphamat, sizeof(double*));
  for (ip=0; ip<nphamat; ip++)
    iphase->mu[ip] = calloc ((size_t) iphase->n[ip], sizeof(double));
      
  iphase->p  = calloc ((size_t) iphase->nscales, sizeof(double**));
  iphase->F  = calloc ((size_t) iphase->nscales, sizeof(double*));
  iphase->A  = calloc ((size_t) iphase->nscales, sizeof(double*));
  iphase->B  = calloc ((size_t) iphase->nscales, sizeof(double*));
  iphase->C  = calloc ((size_t) iphase->nscales, sizeof(double*));

  for (j=0; j<iphase->nscales; j++) {
    /* full phase matrix not needed for delta scaled stuff */
    if (j==0) {
      iphase->p[j]  = calloc ((size_t) nphamat, sizeof(double*));
      for (ip=0; ip<nphamat; ip++) 
	iphase->p[j][ip]  = calloc ((size_t) iphase->n[ip], sizeof(double));
    }
    else {
      iphase->p[j]  = calloc ((size_t) 1, sizeof(double*));
      iphase->p[j][0]  = calloc ((size_t) iphase->n[0], sizeof(double));
    }

    iphase->F[j]  = calloc ((size_t) iphase->n[0], sizeof(double));
    iphase->A[j]  = calloc ((size_t) iphase->n[0], sizeof(double));
    iphase->B[j]  = calloc ((size_t) iphase->n[0], sizeof(double));
    iphase->C[j]  = calloc ((size_t) iphase->n[0], sizeof(double));
    
  }

  return 0;
}


/***********************************************************************************/
/* Function: free_iphase                                                  @62_30i@ */
/* Description:                                                                    */
/*  Free iphase                                                                    */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int free_iphase (pft *iphase)
{
  int ip=0, j=0;

  free(iphase->n);

  for (ip=0; ip<iphase->nphamat; ip++)
    free(iphase->mu[ip]);
  free(iphase->mu);
      
  for (j=0; j<iphase->nscales; j++) {
    /* full phase matrix not needed for delta scaled stuff */
    if (j==0) {
      for (ip=0; ip<iphase->nphamat; ip++) 
	free(iphase->p[j][ip]);
      free(iphase->p[j]);
    }
    else {
      free(iphase->p[j][0]);
      free(iphase->p[j]);
    }

    free(iphase->F[j]);
    free(iphase->A[j]);
    free(iphase->B[j]);
    free(iphase->C[j]);
  }

  free(iphase->p);
  free(iphase->F);
  free(iphase->A);
  free(iphase->B);
  free(iphase->C);

  return 0;
}


/***********************************************************************************/
/* Function: calc_iphase_coeffs                                           @62_30i@ */
/* Description:                                                                    */
/*  Calculated coefficients for iphase                                             */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

void calc_iphase_coeffs (pft *iphase, int SC_mode)
{
  double beta=0.0, alpha=0.0;
  int i=0;

  /* calculate coefficients for linear interpolation of the table */
  for (i=0; i<iphase->n[0]-1; i++) {
    beta  = (iphase->p[SC_mode][0][i+1]-iphase->p[SC_mode][0][i]) / (iphase->mu[0][i+1]-iphase->mu[0][i]);
    alpha = iphase->p[SC_mode][0][i] - beta * iphase->mu[0][i];

    iphase->C[SC_mode][i] = iphase->F[SC_mode][i] - alpha*iphase->mu[0][i] - 0.5*beta*iphase->mu[0][i]*iphase->mu[0][i];
    iphase->B[SC_mode][i] = alpha;
    iphase->A[SC_mode][i] = beta/2.0;
  }
}


/***********************************************************************************/
/* Function: calc_cumulative_table                                        @62_30i@ */
/* Description:                                                                    */
/*  Calculate a table of cumulative probabilities F1[], directions mu1[],          */
/*  and slopes D1[] for use with sc_mu(); the table is calculated from the         */
/*  phase function itself, not from its Legendre coefficients.                     */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int calc_cumulative_table (double **mu, float **p, int *n, int nphamat, double truncate, 
			   pft *iphase, double delta_scaling_mucut, int quiet)
{
  int status=0;  
  int i=0, it=0, ip=0;
  double *F=NULL;
  int nscales=0, itcut=0;
  double delta_scaling_factor=0;
  float scr=0;
  double dif=0;

  if (mu[0][0]!=-1.0) {
    fprintf (stderr, "Error, first mu in phase function must be -1.0!\n");
    fprintf (stderr, "Found %.16f\n", mu[0][0]);
    return -1;
  }

  if (mu[0][n[0]-1]!=1.0) {
    fprintf (stderr, "Error, last mu in phase function must be +1.0!\n");
    fprintf (stderr, "Found %.16f\n", mu[0][n[0]-1]);
    return -1;
  }

  if (delta_scaling_mucut>-1.0)
    nscales = 2;
  else
    nscales = 1;

  status = calloc_iphase (iphase, nscales, nphamat, n);
  if (status) {
    fprintf (stderr, "Error! calloc_iphase did not work! %d\n", status);
    return status;
  }

  /* calculate cumulative phase function by linear integration and normalize phase to 2.0 */
  F = calloc (n[0], sizeof(double));
  normalize_phase ( mu, p, F, n, nphamat, !quiet);

  /* copy to final destination */
  iphase->is_one_theta_grid = 1;
  for (ip=0; ip<nphamat; ip++) 
    for (i=0; i<iphase->n[ip]; i++) {
      iphase->mu[ip][i] = mu[ip][i];
      iphase->p[MCSC_MODE_NORMAL][ip][i]  = p[ip][i];
      /* test whether all phase matrix elements have the same angular grid */
      if (iphase->is_one_theta_grid)
	if (mu[ip][i] != mu[0][i])
	  iphase->is_one_theta_grid = 0;
    }
  
  for (i=0; i<iphase->n[0]; i++)    
    iphase->F[MCSC_MODE_NORMAL][i]  = F[i];
  
  /* calculate coefficients for linear interpolation of the table */
  calc_iphase_coeffs (iphase, MCSC_MODE_NORMAL);

  /* determine whether table is linear in theta */
  scr = (float) ( acos(iphase->mu[0][1]) - acos(iphase->mu[0][0]) );
  iphase->is_linear = 1;
  iphase->theta_0 = acos(iphase->mu[0][0]);
  for (i=2; i<iphase->n[0]; i++) {
    if (scr != (float) ( acos(iphase->mu[0][i]) - acos(iphase->mu[0][i-1]) ) ) {
      iphase->is_linear = 0;
      break;
    }
  }
  if (iphase->is_linear) {
    iphase->dthetainv = 1.0 / ( acos(iphase->mu[0][1]) - acos(iphase->mu[0][0]) );
    if (!quiet)
      fprintf(stderr,"phase function table is linear in theta, switching to fast table lookup!\n");
  }
  
  /**************************************************************************/
  /* von hier TZ ds ... */

  /* initialize delta-scaling factor to 0 = no delta-scaling */
  iphase->dscale = 0.0;

  if (delta_scaling_mucut>-1.0) {
    itcut=n[0];
    dif=99999.;
    for (i=0; i<iphase->n[0]; i++) {
      if (fabs(mu[0][i]-delta_scaling_mucut) < dif) { /*find closest mu in tabulated values*/
	dif=fabs(mu[0][i]-delta_scaling_mucut);
	itcut=i;
      }
    }

    if (!quiet)
      fprintf (stderr, " ... delta-scaling ... truncating phase function at mucut = %g\n", mu[0][itcut]);
        
    delta_scaling_factor = 2.0/F[itcut]; /* F, the integral of p, is 2 */
    iphase->dscale = 1.0-1.0/delta_scaling_factor;
    
    for (i=0; i<n[0]; i++) {
      F[i] =  F[i]*delta_scaling_factor; /*TZ ds: delta-scaled integrated phase function*/
      p[0][i] =  p[0][i]*delta_scaling_factor; /*TZ ds: delta-scaled phase function*/
      if (i>itcut) F[i] = 2.;
      if (i>itcut) p[0][i] = 0;
    }    

    for (i=0; i<iphase->n[0]; i++) {
      iphase->p[MCSC_MODE_DELTA_SCALE][0][i]  = p[0][i];
      iphase->F[MCSC_MODE_DELTA_SCALE][i]  = F[i];
    }  

    /* calculate coefficients for linear interpolation of the table  */
    /* for delta-scaled version of the phase function */
    calc_iphase_coeffs (iphase, MCSC_MODE_DELTA_SCALE);
  }
  /* ... bis hier TZ ds*/  
  /**************************************************************************/

  if (truncate>-1.0) {
    iphase->truncate = 1;
    iphase->mutrnc = truncate;
    
    it = locate (iphase->mu[0], iphase->n[0], iphase->mutrnc);

    if (it==iphase->n[0]-1 || iphase->mutrnc >= 1.0) {
      fprintf (stderr, " ... no need to truncate at %f\n", iphase->mutrnc);
      iphase->truncate = 0;
    }
    else {
      /* value of the phase function at the truncation angle */
      iphase->ptrnc = iphase->p[MCSC_MODE_NORMAL][0][it] + (iphase->mutrnc - iphase->mu[0][it]) / 
	(iphase->mu[0][it+1] - iphase->mu[0][it]) * (iphase->p[MCSC_MODE_NORMAL][0][it+1] - 
                                                     iphase->p[MCSC_MODE_NORMAL][0][it]);

      iphase->Ftrnc = iphase->F[MCSC_MODE_NORMAL][it] + (iphase->mutrnc - mu[0][it]) * 
        (iphase->ptrnc+p[0][it])/2.0;
      iphase->dtrnc = 2.0/(1.0-iphase->mutrnc) * ((2.0-iphase->Ftrnc)/(1.0-iphase->mutrnc) - 
                                                  iphase->ptrnc);
    }
  }
  
  /* free memory */
  free (F);

  return 0;
}


/***********************************************************************************/
/* Function: read_phase_function                                          @62_30i@ */
/* Description:                                                                    */
/*  Read phase function from a 2-column file.                                      */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int read_phase_function (char *filename, 
			 double *r0, double *dr, size_t *nreff, size_t *nphamat,
			 double **extinc, double **albedo, int **ntheta,
			 float ***theta, float ****phase, int *type, int quiet)
{
  int ir=0, itheta=0, status=0;

  *nreff  = 1;  /* number of effective radii       */ 
  *nphamat = 1;  /* number of phase matrix elements */
  
  /* allocate memory for fields */
  *extinc = calloc ((size_t) *nreff, sizeof (double));
  *albedo = calloc ((size_t) *nreff, sizeof (double));
  
  *r0 = -999;
  *dr = -999;
  
  *ntheta = calloc ((size_t) *nreff, sizeof (int));
  *theta  = calloc ((size_t) *nreff, sizeof (double *));
  *phase  = calloc ((size_t) *nreff, sizeof (double **));

  for (ir=0; ir<*nreff; ir++)
    *phase[ir] = calloc (*nphamat, sizeof (double *));
 
  if (!quiet) 
    fprintf (stderr, " ... reading two columns of data from %s\n", filename);

  status = read_2c_file_float (filename, &((*theta)[0]), &((*phase)[0][0]), &((*ntheta)[0]));
  
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }
  
  /* check if the user defined theta or theta */
  for (itheta=0; itheta<(*ntheta)[0]; itheta++) {
    if ((*theta)[0][itheta]<0 || (*theta)[0][itheta]>180.) {
      fprintf (stderr, "Error in %s: Found theta=%f, expecting theta, p(theta)\n", filename, (*theta)[0][itheta]);
      return -1;
    }
  }

  *type = PHASE_EXPLICIT;

  if (!quiet)
    fprintf (stderr, " ... read %d coefficients from %s\n", (*ntheta)[0], filename);

  return 0;
}


/***********************************************************************************/
/* Function: sort_theta_and_mu                                                     */
/* Description:                                                                    */
/*  Take several mu grids and return a common mu grid containing all values of the */
/*  input mu grids                                                                 */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value: (status < 0) => NOT OK, (status == 0) => OK                       */
/* Example:                                                                        */
/* Files: cloud.c                                                                  */
/* Known bugs:                                                                     */
/* Author: Robert Buras, 15.07.2009                                                */    
/*                                                                                 */
/***********************************************************************************/

int sort_theta_and_mu (int n_in, int *nmu_in,
		       float **theta_in, float *theta_out,
		       double **mu_in, double *mu_out)
{
  int i=0, imu=0, loop=0;
  float theta=0.0;
  double mu=0.0;
  int *imu_in=NULL;

  if (n_in==0)
    return 0;

  /* check that borders are identical and equal to -1 and +1 */
  for (i=0;i<n_in;i++) {
    if (mu_in[i][0]!=-1.0) {
      fprintf(stderr,"Error! mu[%d][0] is %e, not -1! Exiting...\n",
	      i,mu_in[i][0]);
      return -1;
    }
    if (mu_in[i][nmu_in[i]-1]!=1.0) {
      fprintf(stderr,"Error! mu[%d][(nmu=%d)-1] %e, is not 1! Exiting...\n",
	      i,nmu_in[i],mu_in[i][nmu_in[i]-1]);
      return -1;
    }
  }

  imu_in = calloc(n_in, sizeof(int));

  /* sort mu's */
  while (loop==0) {
    mu=1.0;
    theta=0.0;
    for (i=0;i<n_in;i++)
      if ( mu > mu_in[i][imu_in[i]] ) {
	mu = mu_in[i][imu_in[i]];
	theta = theta_in[i][imu_in[i]];
      }

    for (i=0;i<n_in;i++)
      if ( mu == mu_in[i][imu_in[i]] )
	imu_in[i]++;

    mu_out[imu] = mu;
    theta_out[imu] = theta;

    imu++;

    loop=1;
    for (i=0;i<n_in;i++)
      if (imu_in[i] < nmu_in[i])
	loop=0;
  }

  free(imu_in);

  return imu;

}


/***********************************************************************************/
/* Function: interpolate_phase_weighted                                            */
/* Description:                                                                    */
/*  Take several phase functions defined on different mu grids, and calculate the  */
/*  weighted mean phase function on a given common mu grid                         */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value: (status < 0) => NOT OK, (status == 0) => OK                       */
/* Example:                                                                        */
/* Files: cloud.c                                                                  */
/* Known bugs:                                                                     */
/* Author: Robert Buras, 15.07.2009                                                */    
/*                                                                                 */
/***********************************************************************************/

int interpolate_phase_weighted (int n_in, int *nmu_in, double **mu_in, float **phase_in,
				double *weight,
				int nmu_out, double *mu_out, float *phase_out,
				int returnmax )
{
  int *imu_in=NULL;
  int imu=0, i=0;
  float *value=NULL, interp_optprop_weight=0.0;

  /* this was a quick fix. It solved problems calling DISORT2 where
     some atmospheric layers had vacuum. These layers then had NANs in
     the phase function. Now 0s are returned. Need to check whether
     this can be done more elegantly; BCA */
  if (n_in == 1)
    if (nmu_in[0] <= 1) {
      for (imu=0; imu<nmu_out; imu++)
	phase_out[imu] = 1.;
      return 0;
    }

  imu_in = calloc(n_in, sizeof(int));
  for (i=0;i<n_in;i++)
    imu_in[i]=0;

  value = calloc(n_in, sizeof(float));

  for (imu=0; imu<nmu_out; imu++) {
    for (i=0;i<n_in;i++) {
      if ( mu_in[i][imu_in[i]] == mu_out[imu] )
	value[i] = phase_in[i][imu_in[i]++];
      else {
	/* interpolate in mu */
	interp_optprop_weight = (float) ( ( mu_out[imu] - mu_in[i][imu_in[i]-1] ) /
					  ( mu_in[i][imu_in[i]] - mu_in[i][imu_in[i]-1] ) );
	if ( interp_optprop_weight > 1. ) {
	  fprintf(stderr,"Error! something went wrong in sum_phase_weighted, imu %d i %d interp_optprop_weight %e\n",imu,i,interp_optprop_weight);
	  fprintf(stderr,"mu_out %e mu_in[i] %e, mu_in[i-1] %e imu_in %d\n",mu_out[imu],mu_in[i][imu_in[i]], mu_in[i][imu_in[i]-1],imu_in[i]);
	  return -1;
	}
	value[i] = interp_optprop_weight * phase_in[i][imu_in[i]]
	  + ( 1. - interp_optprop_weight ) * phase_in[i][imu_in[i]-1];
      }
    }

    if (returnmax) {
      /* derive maximum */
      phase_out[imu] = 0.0;
      for (i=0;i<n_in;i++)
	if (value[i] > phase_out[imu])
	  phase_out[imu] = value[i];
    }
    else {
      /* sum up using weight */
      phase_out[imu] = 0.0;
      for (i=0;i<n_in;i++)
	phase_out[imu] += weight[i] * value[i];
    }

  }

  free (value);
  free (imu_in);

  return 0;

}


/***********************************************************************************/
/* Function: normalize_phase                                                       */
/* Description:                                                                    */
/*  Take a phase function defined on a mu grid and normalize it to 2               */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value: (status < 0) => NOT OK, (status == 0) => OK                       */
/* Example:                                                                        */
/* Files: cloud.c                                                                  */
/* Known bugs:                                                                     */
/* Author: Robert Buras, 15.07.2009                                                */    
/*                                                                                 */
/***********************************************************************************/

void normalize_phase (double **mu, float **p, double *F, int *n, int nphamat, int verbose)
{
  double Ftot=0.0, factor=0.0;
  int i=0, ip=0;

  /* calculate cumulative phase function by linear integration */
  if (F!=NULL) {
    /* want to return cumulative phase function */
    F[0]=0;
    for (i=1; i<n[0]; i++)
      F[i] = F[i-1] + (mu[0][i]-mu[0][i-1])*(p[0][i]+p[0][i-1])/2.0;
    Ftot = F[n[0]-1];
  }
  else
    /* simply calculate Ftot */
    for (i=1; i<n[0]; i++)
      Ftot += (mu[0][i]-mu[0][i-1])*(p[0][i]+p[0][i-1])/2.0;

  /* normalize to 2.0 */
  if (Ftot!=2.0) {
    if (verbose)
      fprintf(stderr,"... renormalizing phase function by factor %f\n",2.0/Ftot);

    factor = 2.0/Ftot;

    for (ip=0; ip<nphamat; ip++)
      for (i=0; i<n[ip]; i++) 
        p[ip][i] *= factor;
    
    if (F!=NULL){
      for (i=0; i<n[0]; i++)
        F[i] *= factor;
    }
  }
}


/***********************************************************************************/
/* Function: sort_and_add_weighted_phase                                           */
/* Description:                                                                    */
/*  Take several phase functions defined on different mu grids, define a common    */
/*  mu grid for them, and calculate the weighted mean phase function               */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value: (status < 0) => NOT OK, (status == 0) => OK                       */
/* Example:                                                                        */
/* Files: cloud.c                                                                  */
/* Known bugs:                                                                     */
/* Author: Robert Buras, 15.07.2009                                                */    
/*                                                                                 */
/***********************************************************************************/

int sort_and_add_weighted_phase (int n_in, double *interp_optprop_weight,
				 int **ntheta, float ***theta, double ***mu, float ***phase,
				 int **ntheta_new, float ***theta_new, double ***mu_new,
				 float ***phase_new, int nphamat, int returnmax,
				 int quiet)
{
  int i=0, n_tot=0, ip=0, status=0;
  float *tmp_theta_new=NULL;
  double *tmp_mu_new=NULL;

  for (ip=0; ip<nphamat; ip++) {
    n_tot=0;
    for (i=0; i<n_in; i++)
      n_tot += ntheta[ip][i];

    tmp_theta_new = calloc(n_tot, sizeof(float));
    tmp_mu_new    = calloc(n_tot, sizeof(double));

    if (n_tot)
      (*ntheta_new)[ip] = sort_theta_and_mu (n_in, ntheta[ip],
					     theta[ip], tmp_theta_new,
					     mu[ip], tmp_mu_new );

    if ( (*ntheta_new)[ip] == -1 ) {
      fprintf(stderr,"Error, sort_theta_and_mu returned error!\n");
      return -1;
    }

    /* allocate theta dimension */
    (*theta_new)[ip]  = calloc ((*ntheta_new)[ip], sizeof(float));
    (*mu_new)   [ip]  = calloc ((*ntheta_new)[ip], sizeof(double));
    (*phase_new)[ip]  = calloc ((*ntheta_new)[ip], sizeof(float));

    /* copy new theta grid into target */
    for (i=0; i<(*ntheta_new)[ip]; i++) {
      (*theta_new) [ip][i] = tmp_theta_new [i];
      (*mu_new)    [ip][i] = tmp_mu_new    [i];
    }

    /* testing */
    /*    if (returnmax) {
      fprintf(stderr,"old grid n :");
      for (i=0;i<n_in;i++)
	fprintf(stderr,"%d ",ntheta[0][i]);
      fprintf(stderr,"\n new grid n: %d\n",(*ntheta_new)[0]);
      }*/

    /* interpolate phase */
    if (n_tot)
      status = interpolate_phase_weighted ( n_in, ntheta[ip], mu[ip], phase[ip],
					    interp_optprop_weight,
					    (*ntheta_new)[ip], (*mu_new)[ip], (*phase_new)[ip],
					    returnmax );

    if (status !=0) {
      fprintf(stderr,"interpolate_phase_weighted returned error %d\n",status);
      return status;
    }

    free(tmp_theta_new);
    free(tmp_mu_new);
  }

  normalize_phase(*mu_new, *phase_new, NULL, *ntheta_new, nphamat, !quiet);

  return 0;
}


/***********************************************************************************/
/* Function: create_phase_from_Rayleigh                                   @45_30i@ */
/* Description: tabellized rayleigh phase function                                 */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i45_30@ */
/***********************************************************************************/

int create_phase_from_Rayleigh (float depol, int *ntheta,
				float **theta, double **mu, float **phase,
				int only_backdir, int ip)
{
  int it=0;

  int nth=38;
  float th[] = {180, 162, 155, 150, 145, 138, 132, 127, 123, 119, 115, 111, 108, 105, 102, 99, 96, 
                93, 90, 87, 84, 81, 78, 75, 72, 69, 65, 61, 57, 53, 48, 45, 42, 35, 30, 25, 18, 0};
  double phase_matrix[6];

  *ntheta = nth;

  if (only_backdir)
    /* special case, for sslidar we only need backward direction */
    *ntheta=1;

  *theta = calloc(*ntheta, sizeof(float));
  *mu    = calloc(*ntheta, sizeof(double));
  *phase = calloc(*ntheta, sizeof(float));

  for (it=0; it<*ntheta; it++) {
    (*theta)[it] = th[it] /180. * PI;
    (*mu)[it] = cos((*theta)[it]);
    get_phase_matrix_rayleigh((*mu)[it], depol, NPHAMAT, phase_matrix);
    (*phase)[it] = phase_matrix[ip]; /*Rayleigh_depol ((*mu)[it], depol);*/
  }

  (*mu)[0] = -1.0;

  /* normalize */
  if (!only_backdir)
    normalize_phase(mu, phase, NULL, ntheta, 1, 0);

  return 0;
}


/***********************************************************************************/
/* Function: create_phase_from_HG                                         @45_30i@ */
/* Description: Create a phasetable for Henyey-Greenstein                          */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i45_30@ */
/***********************************************************************************/

int create_phase_from_HG (double g, int *ntheta, float **theta, double **mu, float **phase, int only_backdir)
{
  int it=0, status=0;

  /* BCA if g>0.99 should use 20.000 grid points equidistant?? */

  /*  if (g > 0.99) {
    fprintf(stderr,"Error, HG asymmetry parameters larger than 0.99 are not allowed! U want %e\n",g);
    return -1;
    } */

  status = choose_HG_grid (g, ntheta, theta);
  if (status) {
    fprintf(stderr,"Error, choose_HG_grid did not work, error %d!\n",status);
    return -1;
  }

  if (only_backdir)
    /* special case, for sslidar we only need backward direction */
    *ntheta=1;

  *mu    = calloc(*ntheta, sizeof(double));
  *phase =  calloc(*ntheta, sizeof(float));
  
  for (it=0; it<*ntheta; it++) {
    (*mu)[it] = cos((*theta)[it]);
    (*phase)[it] = HG (g, (*mu)[it]);
  }

  (*mu)[0] = -1.0;


  /* normalize */
  if (!only_backdir)
    normalize_phase(mu, phase, NULL, ntheta, 1, 0);

  return 0;
}


/***********************************************************************************/
/* Function: choose_HG_grid                                               @45_30i@ */
/* Description: Here the "perfect" grid is chosen for the HG grid according to g   */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i45_30@ */
/***********************************************************************************/

int choose_HG_grid ( double g, int *ntheta, float **theta)
{
  int n = 18;
  int n_values[18] = {15,27,39,58,71,91,109,138,192,198,203,216,231,238,242,272,286,337};
  float g_values[18] = {0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99};
  float *theta_values[18];

  float theta_values_010 [] = {180, 140, 122, 114, 107, 99, 92, 85, 78, 70, 62, 53, 43, 30, 0};
  float theta_values_020 [] = {180, 153, 141, 131, 123, 115, 108, 101, 98, 94, 91, 87, 84, 80, 76, 72, 68, 64, 60, 55, 50, 45, 39, 32, 23, 16, 0};
  float theta_values_030 [] = {180, 153, 141, 132, 124, 120, 116, 112, 109, 105, 102, 99, 95, 92, 88, 85, 81, 77, 73, 69, 65, 63, 61, 59, 56, 54, 51, 49, 46, 43, 40, 36, 32, 28, 23, 20, 16, 11, 0};
  float theta_values_040 [] = {180, 161, 153, 142, 137, 133, 129, 125, 121, 118, 114, 111, 107, 104, 100, 97, 94, 92, 90, 87, 85, 83, 81, 79, 77, 75, 73, 71, 69, 67, 65, 63, 61, 60, 58, 56, 55, 53, 51, 50, 48, 47, 45, 44, 42, 40, 38, 36, 34, 32, 29, 27, 24, 21, 17, 12, 8, 0};
  float theta_values_050 [] = {180, 154, 148, 143, 138, 134, 130, 126, 122, 119, 115, 112, 108, 105, 103, 101, 98, 95, 93, 91, 88, 86, 84, 83, 81, 79, 77, 75, 73, 71, 69, 67, 65, 63, 62, 60, 58, 57, 55, 53, 52, 50, 49, 47, 46, 44, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 31, 29, 28, 26, 25, 23, 21, 19, 17, 16, 14, 12, 10, 7, 0};
  float theta_values_060 [] = {180, 155, 149, 144, 139, 135, 131, 127, 123, 120, 116, 113, 109, 106, 103, 100, 97, 95, 93, 90, 88, 86, 85, 83, 81, 79, 77, 75, 73, 71, 69, 67, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 10, 9, 7, 5, 0};
  float theta_values_070 [] = {180, 162, 155, 149, 144, 139, 135, 131, 127, 123, 120, 116, 113, 109, 106, 103, 100, 97, 94, 93, 91, 89, 87, 86, 84, 82, 80, 78, 76, 74, 72, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25.5, 25, 24.5, 24, 23.5, 23, 22.5, 22, 21.5, 21, 20.5, 20, 19.5, 19, 18.5, 18, 17.5, 17, 16.5, 16, 15.5, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 4, 0};
  float theta_values_080 [] = {180, 162, 155, 149, 144, 139, 135, 131, 128, 124, 121, 117, 114, 112, 110, 107, 104, 101, 98, 95, 93, 92, 90, 88, 87, 85, 83, 81, 79, 77, 75, 73, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37.5, 37, 36.5, 36, 35.5, 35, 34.5, 34, 33.5, 33, 32.5, 32, 31.5, 31, 30.5, 30, 29.5, 29, 28.5, 28, 27.5, 27, 26.5, 26, 25.5, 25, 24.5, 24, 23.5, 23, 22.5, 22, 21.5, 21, 20.5, 20, 19.5, 19, 18.5, 18, 17.5, 17, 16.5, 16, 15.5, 15, 14.5, 14, 13.5, 13, 12.5, 12, 11.5, 11, 10.5, 10, 9.5, 9, 8.5, 8, 7.5, 7, 6.5, 6, 5.5, 5, 4.5, 4, 3, 2, 0};
  float theta_values_090 [] = {180, 162, 155, 149, 144, 140, 136, 132, 128, 124, 121, 117, 114, 112, 110, 107, 104, 101, 98, 96, 95, 93, 92, 90, 88, 87, 85, 83, 81, 79, 77, 75, 74, 73, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40.5, 40, 39.5, 39, 38.5, 38, 37.5, 37, 36.5, 36, 35.5, 35, 34.5, 34, 33.5, 33, 32.5, 32, 31.5, 31, 30.5, 30, 29.5, 29, 28.5, 28, 27.5, 27, 26.5, 26, 25.5, 25, 24.5, 24, 23.5, 23, 22.5, 22, 21.5, 21, 20.5, 20, 19.5, 19.2, 19, 18.7, 18.5, 18.2, 18, 17.7, 17.5, 17.2, 17, 16.7, 16.5, 16.2, 16, 15.7, 15.5, 15.2, 15, 14.7, 14.5, 14.2, 14, 13.7, 13.5, 13.2, 13, 12.7, 12.5, 12.2, 12, 11.7, 11.5, 11.2, 11, 10.7, 10.5, 10.2, 10, 9.7, 9.5, 9.2, 9, 8.8, 8.6, 8.5, 8.3, 8.1, 8, 7.8, 7.6, 7.5, 7.3, 7.1, 7, 6.8, 6.6, 6.5, 6.3, 6.1, 6, 5.8, 5.6, 5.5, 5.3, 5.1, 5, 4.8, 4.6, 4.5, 4.3, 4.1, 4, 3.8, 3.5, 3.3, 3, 2.8, 2.5, 2.3, 2, 1.6, 1, 0};
  float theta_values_091 [] = {180, 162, 155, 149, 144, 139, 135, 131, 128, 124, 121, 117, 114, 112, 110, 107, 104, 101, 98, 96, 95, 93, 92, 90, 88, 87, 85, 83, 81, 79, 77, 75, 74, 73, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40.5, 40, 39.5, 39, 38.5, 38, 37.5, 37, 36.5, 36, 35.5, 35, 34.5, 34, 33.5, 33, 32.5, 32, 31.5, 31, 30.5, 30, 29.5, 29, 28.5, 28, 27.5, 27, 26.5, 26, 25.5, 25, 24.5, 24, 23.5, 23, 22.5, 22, 21.5, 21, 20.5, 20, 19.7, 19.5, 19.2, 19, 18.7, 18.5, 18.2, 18, 17.7, 17.5, 17.2, 17, 16.7, 16.5, 16.2, 16, 15.7, 15.5, 15.2, 15, 14.7, 14.5, 14.2, 14, 13.7, 13.5, 13.2, 13, 12.7, 12.5, 12.2, 12, 11.7, 11.5, 11.2, 11, 10.7, 10.5, 10.2, 10, 9.8, 9.7, 9.5, 9.3, 9.1, 9, 8.8, 8.6, 8.5, 8.3, 8.1, 8, 7.8, 7.6, 7.5, 7.3, 7.1, 7, 6.8, 6.6, 6.5, 6.3, 6.1, 6, 5.8, 5.6, 5.5, 5.3, 5.1, 5, 4.8, 4.6, 4.5, 4.3, 4.1, 4, 3.8, 3.6, 3.5, 3.3, 3.1, 3, 2.8, 2.6, 2.3, 2, 1.6, 1.3, 1, 0};
  float theta_values_092 [] = {180, 162, 155, 149, 144, 140, 136, 132, 128, 124, 121, 117, 114, 112, 110, 107, 104, 101, 98, 96, 95, 93, 92, 90, 88, 87, 85, 83, 81, 79, 77, 75, 73, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40.5, 40, 39.5, 39, 38.5, 38, 37.5, 37, 36.5, 36, 35.5, 35, 34.5, 34, 33.5, 33, 32.5, 32, 31.5, 31, 30.5, 30, 29.5, 29, 28.5, 28, 27.5, 27, 26.5, 26, 25.5, 25, 24.5, 24, 23.5, 23, 22.5, 22, 21.5, 21, 20.5, 20.2, 20, 19.7, 19.5, 19.2, 19, 18.7, 18.5, 18.2, 18, 17.7, 17.5, 17.2, 17, 16.7, 16.5, 16.2, 16, 15.7, 15.5, 15.2, 15, 14.7, 14.5, 14.2, 14, 13.7, 13.5, 13.2, 13, 12.7, 12.5, 12.2, 12, 11.7, 11.5, 11.2, 11, 10.8, 10.7, 10.5, 10.3, 10.2, 10, 9.8, 9.7, 9.5, 9.3, 9.1, 9, 8.8, 8.6, 8.5, 8.3, 8.1, 8, 7.8, 7.6, 7.5, 7.3, 7.1, 7, 6.8, 6.6, 6.5, 6.3, 6.1, 6, 5.8, 5.6, 5.5, 5.3, 5.1, 5, 4.8, 4.6, 4.5, 4.3, 4.1, 4, 3.8, 3.6, 3.5, 3.3, 3.1, 3, 2.8, 2.6, 2.5, 2.3, 2.2, 2, 1.8, 1.6, 1.3, 1, 0};
  float theta_values_093 [] = {180, 162, 155, 149, 144, 140, 136, 132, 128, 124, 121, 117, 114, 112, 110, 107, 104, 101, 98, 96, 95, 93, 92, 90, 88, 87, 85, 83, 81, 79, 77, 75, 74, 73, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40.5, 40, 39.5, 39, 38.5, 38, 37.5, 37, 36.5, 36, 35.5, 35, 34.5, 34, 33.5, 33, 32.5, 32, 31.5, 31, 30.5, 30, 29.5, 29, 28.5, 28, 27.5, 27, 26.5, 26, 25.5, 25, 24.5, 24, 23.5, 23, 22.5, 22, 21.5, 21, 20.5, 20.2, 20, 19.7, 19.5, 19.2, 19, 18.7, 18.5, 18.2, 18, 17.7, 17.5, 17.2, 17, 16.7, 16.5, 16.2, 16, 15.7, 15.5, 15.2, 15, 14.7, 14.5, 14.2, 14, 13.7, 13.5, 13.2, 13, 12.7, 12.5, 12.2, 12, 11.7, 11.5, 11.3, 11.2, 11, 10.8, 10.7, 10.5, 10.3, 10.2, 10, 9.8, 9.6, 9.5, 9.3, 9.1, 9, 8.8, 8.6, 8.5, 8.3, 8.1, 8, 7.8, 7.6, 7.5, 7.3, 7.1, 7, 6.8, 6.6, 6.5, 6.3, 6.1, 6, 5.8, 5.6, 5.5, 5.4, 5.3, 5.2, 5.1, 5, 4.9, 4.8, 4.7, 4.6, 4.5, 4.4, 4.3, 4.2, 4.1, 4, 3.9, 3.8, 3.7, 3.6, 3.5, 3.4, 3.3, 3.2, 3.1, 3, 2.8, 2.6, 2.5, 2.3, 2.2, 2, 1.8, 1.6, 1.3, 1, 0.7, 0};
  float theta_values_094 [] = {180, 162, 155, 149, 144, 140, 136, 132, 128, 124, 121, 117, 114, 112, 110, 107, 104, 101, 98, 97, 95, 93, 92, 90, 88, 87, 85, 83, 81, 79, 77, 75, 74, 73, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40.5, 40, 39.5, 39, 38.5, 38, 37.5, 37, 36.5, 36, 35.5, 35, 34.5, 34, 33.5, 33, 32.5, 32, 31.5, 31, 30.5, 30, 29.5, 29, 28.5, 28, 27.5, 27, 26.5, 26, 25.5, 25, 24.5, 24, 23.5, 23, 22.5, 22, 21.5, 21, 20.7, 20.5, 20.2, 20, 19.7, 19.5, 19.2, 19, 18.7, 18.5, 18.2, 18, 17.7, 17.5, 17.2, 17, 16.7, 16.5, 16.2, 16, 15.7, 15.5, 15.2, 15, 14.7, 14.5, 14.2, 14, 13.7, 13.5, 13.2, 13, 12.7, 12.5, 12.2, 12, 11.8, 11.7, 11.5, 11.3, 11.2, 11, 10.8, 10.7, 10.5, 10.3, 10.2, 10, 9.8, 9.6, 9.5, 9.3, 9.1, 9, 8.8, 8.6, 8.5, 8.3, 8.1, 8, 7.8, 7.6, 7.5, 7.3, 7.1, 7, 6.9, 6.8, 6.7, 6.6, 6.5, 6.4, 6.3, 6.2, 6.1, 6, 5.9, 5.8, 5.7, 5.6, 5.5, 5.4, 5.3, 5.2, 5.1, 5, 4.9, 4.8, 4.7, 4.6, 4.5, 4.4, 4.3, 4.2, 4.1, 4, 3.9, 3.8, 3.7, 3.6, 3.5, 3.4, 3.3, 3.2, 3.1, 3, 2.9, 2.8, 2.7, 2.6, 2.5, 2.4, 2.3, 2.2, 2.1, 2, 1.9, 1.8, 1.6, 1.5, 1.3, 1.2, 1, 0.7, 0};
  float theta_values_095 [] = {180, 162, 155, 149, 144, 139, 135, 131, 128, 124, 121, 117, 114, 112, 110, 107, 104, 101, 98, 96, 95, 93, 92, 90, 88, 86, 85, 83, 81, 79, 77, 75, 74, 73, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41.5, 41, 40.5, 40, 39.5, 39, 38.5, 38, 37.5, 37, 36.5, 36, 35.5, 35, 34.5, 34, 33.5, 33, 32.5, 32, 31.5, 31, 30.5, 30, 29.5, 29, 28.5, 28, 27.5, 27, 26.5, 26, 25.5, 25, 24.5, 24, 23.5, 23, 22.5, 22, 21.5, 21, 20.7, 20.5, 20.2, 20, 19.7, 19.5, 19.2, 19, 18.7, 18.5, 18.2, 18, 17.7, 17.5, 17.2, 17, 16.7, 16.5, 16.2, 16, 15.7, 15.5, 15.2, 15, 14.7, 14.5, 14.2, 14, 13.7, 13.5, 13.2, 13, 12.7, 12.5, 12.2, 12, 11.8, 11.7, 11.5, 11.3, 11.2, 11, 10.8, 10.7, 10.5, 10.3, 10.2, 10, 9.8, 9.6, 9.5, 9.3, 9.1, 9, 8.8, 8.6, 8.5, 8.3, 8.1, 8, 7.8, 7.6, 7.5, 7.4, 7.3, 7.2, 7.1, 7, 6.9, 6.8, 6.7, 6.6, 6.5, 6.4, 6.3, 6.2, 6.1, 6, 5.9, 5.8, 5.7, 5.6, 5.5, 5.4, 5.3, 5.2, 5.1, 5, 4.9, 4.8, 4.7, 4.6, 4.5, 4.4, 4.3, 4.2, 4.1, 4, 3.9, 3.8, 3.7, 3.6, 3.5, 3.4, 3.3, 3.2, 3.1, 3, 2.9, 2.8, 2.7, 2.6, 2.5, 2.4, 2.3, 2.2, 2.1, 2, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1, 0.9, 0.7, 0.5, 0};
  float theta_values_096 [] = {180, 162, 155, 149, 144, 139, 135, 131, 128, 124, 121, 118, 115, 113, 111, 108, 105, 102, 99, 96, 95, 93, 91, 89, 87, 86, 84, 82, 80, 78, 76, 74, 73, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40.5, 40, 39.5, 39, 38.5, 38, 37.5, 37, 36.5, 36, 35.5, 35, 34.5, 34, 33.5, 33, 32.5, 32, 31.5, 31, 30.5, 30, 29.5, 29, 28.5, 28, 27.5, 27, 26.5, 26, 25.5, 25, 24.5, 24, 23.5, 23, 22.5, 22, 21.5, 21.2, 21, 20.7, 20.5, 20.2, 20, 19.7, 19.5, 19.2, 19, 18.7, 18.5, 18.2, 18, 17.7, 17.5, 17.2, 17, 16.7, 16.5, 16.2, 16, 15.7, 15.5, 15.2, 15, 14.7, 14.5, 14.2, 14, 13.7, 13.5, 13.2, 13, 12.7, 12.5, 12.3, 12.2, 12, 11.8, 11.7, 11.5, 11.3, 11.2, 11, 10.8, 10.7, 10.5, 10.3, 10.2, 10, 9.8, 9.6, 9.5, 9.3, 9.1, 9, 8.8, 8.6, 8.5, 8.3, 8.1, 8, 7.9, 7.8, 7.7, 7.6, 7.5, 7.4, 7.3, 7.2, 7.1, 7, 6.9, 6.8, 6.7, 6.6, 6.5, 6.4, 6.3, 6.2, 6.1, 6, 5.9, 5.8, 5.7, 5.6, 5.5, 5.4, 5.3, 5.2, 5.1, 5, 4.9, 4.8, 4.7, 4.6, 4.5, 4.4, 4.3, 4.2, 4.1, 4, 3.9, 3.8, 3.7, 3.6, 3.5, 3.4, 3.3, 3.2, 3.1, 3, 2.9, 2.8, 2.7, 2.6, 2.5, 2.4, 2.3, 2.2, 2.1, 2, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 1, 0.9, 0.8, 0.7, 0.5, 0};
  float theta_values_097 [] = {180, 162, 155, 149, 144, 140, 136, 132, 128, 124, 121, 118, 115, 113, 111, 108, 105, 102, 99, 96, 94, 93, 91, 89, 88, 86, 84, 82, 80, 78, 76, 74, 73, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40.5, 40, 39.5, 39, 38.5, 38, 37.5, 37, 36.5, 36, 35.5, 35, 34.5, 34, 33.5, 33, 32.5, 32, 31.5, 31, 30.5, 30, 29.5, 29, 28.5, 28, 27.5, 27, 26.5, 26, 25.5, 25, 24.5, 24, 23.5, 23, 22.5, 22, 21.5, 21.2, 21, 20.7, 20.5, 20.2, 20, 19.7, 19.5, 19.2, 19, 18.7, 18.5, 18.2, 18, 17.7, 17.5, 17.2, 17, 16.7, 16.5, 16.2, 16, 15.7, 15.5, 15.2, 15, 14.7, 14.5, 14.2, 14, 13.7, 13.5, 13.2, 13, 12.7, 12.5, 12.3, 12.2, 12, 11.8, 11.7, 11.5, 11.3, 11.2, 11, 10.8, 10.7, 10.5, 10.3, 10.2, 10, 9.8, 9.6, 9.5, 9.3, 9.1, 9, 8.8, 8.6, 8.5, 8.3, 8.2, 8.1, 8, 7.9, 7.8, 7.7, 7.6, 7.5, 7.4, 7.3, 7.2, 7.1, 7, 6.9, 6.8, 6.7, 6.6, 6.5, 6.4, 6.3, 6.2, 6.1, 6, 5.9, 5.8, 5.7, 5.6, 5.5, 5.4, 5.3, 5.2, 5.1, 5, 4.9, 4.8, 4.7, 4.6, 4.5, 4.4, 4.3, 4.2, 4.1, 4, 3.9, 3.8, 3.7, 3.6, 3.5, 3.45, 3.4, 3.35, 3.3, 3.25, 3.2, 3.15, 3.1, 3.05, 3, 2.95, 2.9, 2.85, 2.8, 2.75, 2.7, 2.65, 2.6, 2.55, 2.5, 2.45, 2.4, 2.35, 2.3, 2.25, 2.2, 2.15, 2.1, 2.05, 2, 1.95, 1.9, 1.85, 1.8, 1.75, 1.7, 1.65, 1.6, 1.55, 1.5, 1.45, 1.4, 1.35, 1.3, 1.25, 1.2, 1.15, 1.1, 1.05, 1, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0};
  float theta_values_098 [] = {180, 162, 155, 149, 144, 140, 136, 132, 128, 124, 121, 118, 115, 113, 111, 108, 105, 102, 99, 96, 95, 93, 91, 89, 88, 86, 84, 82, 80, 78, 76, 74, 73, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40.5, 40, 39.5, 39, 38.5, 38, 37.5, 37, 36.5, 36, 35.5, 35, 34.5, 34, 33.5, 33, 32.5, 32, 31.5, 31, 30.5, 30, 29.5, 29, 28.5, 28, 27.5, 27, 26.5, 26, 25.5, 25, 24.5, 24, 23.5, 23, 22.5, 22, 21.5, 21.2, 21, 20.7, 20.5, 20.2, 20, 19.7, 19.5, 19.2, 19, 18.7, 18.5, 18.2, 18, 17.7, 17.5, 17.2, 17, 16.7, 16.5, 16.2, 16, 15.7, 15.5, 15.2, 15, 14.7, 14.5, 14.2, 14, 13.7, 13.5, 13.2, 13, 12.8, 12.7, 12.5, 12.3, 12.2, 12, 11.8, 11.7, 11.5, 11.3, 11.2, 11, 10.8, 10.7, 10.5, 10.3, 10.1, 10, 9.8, 9.6, 9.5, 9.3, 9.1, 9, 8.8, 8.6, 8.5, 8.4, 8.3, 8.2, 8.1, 8, 7.9, 7.8, 7.7, 7.6, 7.5, 7.4, 7.3, 7.2, 7.1, 7, 6.9, 6.8, 6.7, 6.6, 6.5, 6.4, 6.3, 6.2, 6.1, 6, 5.9, 5.8, 5.7, 5.6, 5.5, 5.4, 5.3, 5.2, 5.1, 5, 4.9, 4.8, 4.7, 4.6, 4.5, 4.4, 4.3, 4.2, 4.1, 4, 3.95, 3.9, 3.85, 3.8, 3.75, 3.7, 3.65, 3.6, 3.55, 3.5, 3.45, 3.4, 3.35, 3.3, 3.25, 3.2, 3.15, 3.1, 3.05, 3, 2.95, 2.9, 2.85, 2.8, 2.75, 2.7, 2.65, 2.6, 2.55, 2.5, 2.45, 2.4, 2.35, 2.3, 2.25, 2.2, 2.15, 2.1, 2.05, 2, 1.95, 1.9, 1.85, 1.8, 1.75, 1.7, 1.65, 1.6, 1.55, 1.5, 1.45, 1.4, 1.35, 1.3, 1.25, 1.2, 1.15, 1.1, 1.05, 1, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.2, 0};
  float theta_values_099 [] = {180, 162, 155, 149, 144, 140, 136, 132, 128, 124, 121, 118, 115, 113, 111, 108, 105, 102, 99, 96, 94, 93, 91, 89, 88, 86, 84, 82, 80, 78, 76, 74, 73, 72, 71, 70, 69, 68, 67, 66, 65, 64, 63, 62, 61, 60, 59, 58, 57, 56, 55, 54, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41.5, 41, 40.5, 40, 39.5, 39, 38.5, 38, 37.5, 37, 36.5, 36, 35.5, 35, 34.5, 34, 33.5, 33, 32.5, 32, 31.5, 31, 30.5, 30, 29.5, 29, 28.5, 28, 27.5, 27, 26.5, 26, 25.5, 25, 24.5, 24, 23.5, 23, 22.5, 22, 21.5, 21.2, 21, 20.7, 20.5, 20.2, 20, 19.7, 19.5, 19.2, 19, 18.7, 18.5, 18.2, 18, 17.7, 17.5, 17.2, 17, 16.7, 16.5, 16.2, 16, 15.7, 15.5, 15.2, 15, 14.7, 14.5, 14.2, 14, 13.7, 13.5, 13.2, 13, 12.8, 12.7, 12.5, 12.3, 12.2, 12, 11.8, 11.7, 11.5, 11.3, 11.2, 11, 10.8, 10.7, 10.5, 10.3, 10.1, 10, 9.8, 9.6, 9.5, 9.3, 9.1, 9, 8.8, 8.6, 8.5, 8.4, 8.3, 8.2, 8.1, 8, 7.9, 7.8, 7.7, 7.6, 7.5, 7.4, 7.3, 7.2, 7.1, 7, 6.9, 6.8, 6.7, 6.6, 6.5, 6.4, 6.3, 6.2, 6.1, 6, 5.9, 5.8, 5.7, 5.6, 5.5, 5.4, 5.3, 5.2, 5.1, 5, 4.9, 4.8, 4.7, 4.6, 4.5, 4.4, 4.3, 4.25, 4.2, 4.15, 4.1, 4.05, 4, 3.95, 3.9, 3.85, 3.8, 3.75, 3.7, 3.65, 3.6, 3.55, 3.5, 3.45, 3.4, 3.35, 3.3, 3.25, 3.2, 3.15, 3.1, 3.05, 3, 2.95, 2.9, 2.85, 2.8, 2.75, 2.7, 2.65, 2.6, 2.55, 2.5, 2.45, 2.4, 2.35, 2.3, 2.25, 2.2, 2.15, 2.1, 2.05, 2, 1.95, 1.93, 1.9, 1.88, 1.85, 1.83, 1.8, 1.78, 1.75, 1.73, 1.7, 1.68, 1.65, 1.63, 1.6, 1.58, 1.55, 1.53, 1.5, 1.48, 1.45, 1.43, 1.4, 1.38, 1.35, 1.33, 1.3, 1.28, 1.25, 1.23, 1.2, 1.18, 1.15, 1.13, 1.1, 1.08, 1.05, 1.03, 1, 0.98, 0.95, 0.93, 0.9, 0.88, 0.85, 0.83, 0.82, 0.8, 0.78, 0.77, 0.75, 0.73, 0.72, 0.7, 0.68, 0.67, 0.65, 0.63, 0.62, 0.6, 0.58, 0.57, 0.55, 0.53, 0.52, 0.5, 0.48, 0.47, 0.45, 0.43, 0.42, 0.4, 0.38, 0.35, 0.33, 0.3, 0.28, 0.26, 0.23, 0.2, 0.16, 0.1, 0};

  int i=0,j=0;
  double gabs=0;

  theta_values[0] = theta_values_010;
  theta_values[1] = theta_values_020;
  theta_values[2] = theta_values_030;
  theta_values[3] = theta_values_040;
  theta_values[4] = theta_values_050;
  theta_values[5] = theta_values_060;
  theta_values[6] = theta_values_070;
  theta_values[7] = theta_values_080;
  theta_values[8] = theta_values_090;
  theta_values[9] = theta_values_091;
  theta_values[10] = theta_values_092;
  theta_values[11] = theta_values_093;
  theta_values[12] = theta_values_094;
  theta_values[13] = theta_values_095;
  theta_values[14] = theta_values_096;
  theta_values[15] = theta_values_097;
  theta_values[16] = theta_values_098;
  theta_values[17] = theta_values_099;

  gabs = (g>0?g:-g) - MC_EPSILON; /* in case of 0.99 could lie an epsilon above */

  for (i=0; i<n; i++) 
    if (gabs < g_values[i])
      break;

  if (i==n) {
    if (gabs > 1.0 - 2.*MC_EPSILON) {
      fprintf(stderr,"Error! your asymmetry parameter is %e, this is too close to 1 !!! This can not work! Exiting...\n",g);
      return -1;
    }

    i=n-1;
    fprintf(stderr,"================================================================\n");
    fprintf(stderr,"WARNING! You have chosen an asymmetry parameter g= %e > 0.99!\n",gabs);
    fprintf(stderr,"The gridding chosen for the phase function is not fine enough\n");
    fprintf(stderr,"to resolve such a spiky phase function! Your results might be\n");
    fprintf(stderr,"wrong! Please contact the developers if you really need such\n");
    fprintf(stderr,"large asymmetry parameters!\n");
    fprintf(stderr,"================================================================\n");
  }

  *ntheta = n_values[i];
  (*theta) = calloc(n_values[i], sizeof(float));

  if (g>0)
    for (j=0; j<*ntheta; j++)
      (*theta)[j] = theta_values[i][j] /180. * PI;
  else
    for (j=0; j<*ntheta; j++)
      (*theta)[j] = (180. - theta_values[i][*ntheta-1-j]) /180. * PI;

  return 0;
}


/***********************************************************************************/
/* Function: calc_pmom                                                    @62_30i@ */
/* Description:                                                                    */
/*  Calculates Legendre moments from the phase function using the super fast       */
/*  method from Buras, Dowling, Emde 201X.                                         */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

void calc_pmom (double *mu, double *phase, int ntheta, int nmom,
		double *pmom)
{
  int i=0, l=0;
  double p0=0.0, p1=0.0, p2=0.0;
  double C=0.0, D=0.0;
  double Cp=0.0, Dp=0.0;
  double mu2[ntheta];

  for (i=0; i<ntheta; i++)
    mu2[i] = mu[i] * mu[i];

  /* sum up for zeroth legendre coefficients with simple formula */
  for (i=0; i<ntheta-1; i++)
    pmom[0] += 0.5 * ( phase[i+1] + phase[i] ) * ( mu[i+1] - mu[i] );

  for (i=0; i<ntheta; i++) { /* ntheta is number of grid points */

    /* linear fit for phase function between grid points */
    if (i<ntheta-1) {
      Dp = ( phase[i+1] - phase[i] ) / ( mu[i+1] - mu[i] );
      Cp = phase[i] - mu[i] * Dp;
    }
    else {
      Dp = 0.0;
      Cp = 0.0;
    }

    /* sum up for first legendre coefficients with simple formula */
    pmom[1] += 0.5 * ( C - Cp ) * mu2[i]
      + 1./3. * ( D - Dp ) * mu2[i] * mu[i];

    p0 = 1.0;
    p1=mu[i];

    /* move to l+2 */
    for (l=0; l<2; l++) {
      p2=legendre_polynomial(p1, p0, mu[i], l+2);
      p0 = p1;
      p1 = p2;
    }
    
    /* calculate coefficients using the formula from the paper (CDE 201X) */    
    for (l=2; l<nmom; l++) {
      p2=legendre_polynomial(p1, p0, mu[i], l+2);
      
      pmom[l] += ( C - Cp ) * ( p1 - mu[i] * p0 ) / ((double)l)
	+ ( D - Dp ) * ( ( p1 * mu[i] * ((double)l+2) - p2 ) / ((double)l+1)
			 - p0 * mu2[i] ) / ((double)l-1);
      
      p0 = p1;
      p1 = p2;
    }

    D = Dp;
    C = Cp;

  } /* end for i */
}


/***********************************************************************************/
/* Function: legendre_polynomial                                          @62_30i@ */
/* Description:                                                                    */
/*  Calculates Legendre polynomial for given P_l-1, P_l-2 and cos(theta)           */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static inline double legendre_polynomial (double plm1, double plm2,
					  double mu, int l)
{
  return ( (2*(double)l-1) * mu * plm1 - ((double)l-1) * plm2 ) / ((double)l); 
}


/***********************************************************************************/
/* Function: get_phase_matrix_rayleigh                                    @62_30i@ */
/* Description:                                                                    */
/*  Calculates the scattering phase matrix for Rayleigh scattering.                */
/*  Depolarization is considered (cp. Hansen and Travis, 1974).                    */ 
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Claudia Emde                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

void get_phase_matrix_rayleigh ( double  mu,
				 double  depol,
				 int     np,
				 /* output */
				 double *phase_matrix )
{
  double delta=0.0, delta2=0.0; 

  delta=(1.0-depol)/(1.0+depol/2.0);

  /* Rayleigh phase matrix, indices as in Evans 1991, */
  /* formulas as in Hansen and Travis 1974, P11 is equivalent to phase function */
  /* defined in Rayleigh_depol */
  phase_matrix[0]=delta*(0.75*(1.0+mu*mu))+(1.0-delta);    /* P11 */

  if (np==NPHAMAT) {
    delta2=(1.0-2*depol)/(1.0-depol);
    phase_matrix[1]=delta*(0.75*(-1.0+mu*mu));             /* P12 = P21*/
    phase_matrix[2]=delta*1.5*mu;                          /* P33 */
    phase_matrix[3]=0.0;                                   /* P34 = -P43 */
    phase_matrix[4]=phase_matrix[0]-(1.0-delta);           /* P22 */
    phase_matrix[5]=phase_matrix[2]*delta2;                /* P44 */
  }
}


/***********************************************************************************/
/* Function: get_phase_matrix_pft                                         @62_30i@ */
/* Description:                                                                    */
/*  Extract the phase function or matrix for a given polar angle mu from a         */
/*  phase function table.                                                          */
/*  Merged version of get_phase and get_phase_pol                                  */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Claudia Emde, Robert Buras                                              */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int get_phase_matrix_pft ( pft    *iphase,
			   double  mu,
			   int     scaled,
			   int     np,
			   /* output */
			   double *p )
{
  int it=-1, ip=0;
  
  /* make sure that mu is physical, this should be checked directly at mu calculation site! */
  /* if (mu > 1.0) mu=1.0;
     if (mu < -1.0) mu=-1.0;
  */

  if (np>1)
    np=iphase->nphamat;

  for (ip=0; ip<np; ip++){
    /* if inside the truncation range ???? what */
    if (iphase->truncate && mu > iphase->mutrnc && ip==0)
      p[0] = iphase->ptrnc + iphase->dtrnc * (mu - iphase->mutrnc);
    else {
      if (!iphase->is_one_theta_grid || it==-1) { /* later on, it will only be calced once */
	/* This can work only if for all phase matrix elements the scattering */
	/* angle grid is the same */
	if (iphase->is_linear) {
	  /* this is a new version of getting it, it is much quicker,
	     but less precise, CHECK */
	  it = (int) ((acos(mu) - iphase->theta_0)*iphase->dthetainv);
	  /* careful with rounding errors at the array boundaries */
	  if (it>iphase->n[ip]-1 || it<0) {
	    fprintf(stderr, "ERROR! mu is outside phasetable! Probably only infinitesimally, fix this in get_phase\n");
	    fprintf(stderr, "it %d nit %d mu %e rit %e \n",it,iphase->n[ip]-1,mu,((acos(mu) - iphase->theta_0)*iphase->dthetainv));
	    return -1;
	  }
	}
	else
	  it = locate (iphase->mu[ip], iphase->n[ip], mu);
      }

      if (it==iphase->n[ip]-1)
        p[ip] = iphase->p[scaled][ip][iphase->n[ip]-1];
      else
        p[ip] = iphase->p[scaled][ip][it] + (mu - iphase->mu[ip][it]) / 
          (iphase->mu[ip][it+1] - iphase->mu[ip][it])
          * (iphase->p[scaled][ip][it+1] - iphase->p[scaled][ip][it]);
    }
    /* phase function is only scaled for ip=0, for ip>0 use original! */
    scaled=0;
  } 

  /* spherical particles */
  if (np==4 && NPHAMAT==6) {
    p[4] = p[0];
    p[5] = p[2];
  }

  return 0;
}


/***********************************************************************************/
/* Function: get_phase_matrix_pft_interpol_reff                           @62_30i@ */
/* Description:                                                                    */
/*  Calculate the probability/phase matrix for caoth from phase function table     */
/*  for scattering angle mu = cos(theta).                                          */
/*  Merged version of scattering_probability_cld and phase_matrix_cld              */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras, Claudia Emde                                              */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int get_phase_matrix_pft_interpol_reff ( double                mu,
					 double                reff,
					 phase_function_table *phase_cld,
					 int                   scaled,
					 int                   np,
					 /* Output */
					 double               *phase_matrix_cld )
{
  int i1=0, i2=0, ip=0;
  double r1=0.0, temp=0.0;
  double *phase_matrix_cld2=NULL;
  int status=0;

  if (phase_cld->n == 1){   /* only one phase function specified */
    status = get_phase_matrix_pft ( phase_cld->iphase[0], mu, scaled, np, phase_matrix_cld );
    if (status)
      return fct_err_out (status, "get_phase_matrix_pft", ERROR_POSITION);
  }
  else {
    /* determine lower and upper reff's */
    temp = (reff - phase_cld->r0) / phase_cld->dr;
    /* this replaces "floor" and "ceil", it is much quicker, but less precise */
    i1 = (int) (temp);
    i2 = i1;
    if ( (float) i1 != (float) temp ) i2++;

    /* careful with rounding errors at the array boundaries */
    i1 = (i1>phase_cld->n-1?phase_cld->n-1:i1);
    i2 = (i2>phase_cld->n-1?phase_cld->n-1:i2);

    i1 = (i1<0?0:i1);
    i2 = (i2<0?0:i2);

    status = get_phase_matrix_pft ( phase_cld->iphase[i1], mu, scaled, np, phase_matrix_cld );
    if (status)
      return fct_err_out (status, "get_phase_matrix_pft", ERROR_POSITION);

    if (i2>i1) {  /* need to do linear interpolation */

      phase_matrix_cld2=calloc(np, sizeof(double));

      status = get_phase_matrix_pft ( phase_cld->iphase[i2], mu, scaled, np, phase_matrix_cld2 );
      if (status)
	return fct_err_out (status, "get_phase_matrix_pft", ERROR_POSITION);

      r1 = phase_cld->r0 + (double) i1 * phase_cld->dr;
      for (ip=0; ip<np; ip++)
	/* bug fix, weighting with beta_ext * ssa was missing */
	phase_matrix_cld[ip] += ( reff - r1 ) * phase_cld->sca[i2] /
	  ( ( phase_cld->dr - reff + r1 ) * phase_cld->sca[i1]
	    + ( reff - r1 ) * phase_cld->sca[i2] )
	  * ( phase_matrix_cld2[ip] - phase_matrix_cld[ip] );

      free(phase_matrix_cld2);
    }
  }

  return 0;
}
