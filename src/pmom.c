/*--------------------------------------------------------------------
 * $Id: pmom.c 3030 2014-06-11 14:37:08Z Claudia.Emde $
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
#include <getopt.h>
#include <time.h>

#include "ascii.h"
#include "numeric.h"


#ifndef PI
#define PI 3.14159265358979323846264338327
#endif


#define PROGRAM "pmom"
#define VERSION "1.0"

void calc_pmom (double *mu, double *phase, int ntheta, int nmom,
		double *pmom);
static inline double legendre_polynomial (double p1, double p0,
					  double mu, int l);


/*
  Documentation in latex for the User's Guide. Extracted by sdoc2.awk 
  and input in the doc/tools.tex.
*/
/*
  <lpdoc>

\subsection{Perform Legendre decomposition of phase function - \codeidx{pmom}}

The \code{pmom} tool calculates the Legendre moments of a given phase function. 
The input must be provided as 2-column file, containing the scattering
angle grid in the first column and the phase function value in the second column.
The output of \code{pmom} are the Legendre moments. 

The \code{pmom} tool is invoked on the command for instance as 
\begin{Verbatim}[fontsize=\footnotesize]
  pmom [options] <filename>
\end{Verbatim}

The following optional arguments may be specified:
\begin{description}
\item[-h] Display help message.
\item[-l $<$number$>$] Number of Legendre moments to be computed. In order
  to obtain an accurate decomposition of the phase function, the last
  terms of the Legendre series should approach 0.
\item[-r $<$grid$>$] Specify scattering angle grid which is used
  internally (see below for more explanation).
\item[-c] Calculate coefficients instead of polynomials (these
  include the factor $(2l+1)$. \code{uvspec} requires Legendre
  coefficients.
\item[-n] Normalize the phase function before computing the Legendre
  moments.
\end{description}
 
You may specify the number of moments using the option \code{-l}.
Different scattering angle grid resolutions can be chosen using the
option \code{-r}. For moderate forward peaks, the standard grid
(\code{-r 1} - equidistant, 0.01 degrees step width) should be
sufficient. For phase functions with a very strong forward peak,
e.g. ice particle phase functions, the finest grid resolution
(\code{-r 2} - equidistant, 0.001 degrees step width) should be
specified. If the grid of the input file should be used for the
Legendre decomposition, please use \code{-r 3}; this option uses the
new speedy and exact method for Legendre decomposition (Buras,
Dowling, Emde 201X). Default. You may test \code{-r 4} and \code{-r 5}, in
this case non-equdistant grids with a finer resolution around the
forward peak are used.

You may test the accuracy of the Legendre decomposition by
using the tool \code{phase}:
\begin{Verbatim}[fontsize=\footnotesize]
  phase -c -d -s 1 pmom_outfile.dat 
\end{Verbatim}
 
</lpdoc>
*/

/*************************************************************/
/* print usage information                                   */
/*************************************************************/

static void print_usage (char *filename)
{
  fprintf (stderr, "%s %s - Calculate Legendre moments\n\n", PROGRAM, VERSION);
  fprintf (stderr, "written by Claudia Emde,\n");
  fprintf (stderr, "           DLR, e-Mail claudia.emde@dlr.de\n");
  fprintf (stderr, "Version %s finished February 29, 2010\n\n", VERSION);
  fprintf (stderr, "Optional arguments:\n");
  fprintf (stderr, "-h            display help message\n");
  fprintf (stderr, "-l <number>   number of legendre moments\n");
  fprintf (stderr, "-r <grid no>  scattering angle grid resolution (1-5)\n");
  fprintf (stderr, "-c            calculate coeffcients instead of polynomials,\n");
  fprintf (stderr, "              required to generate input files for uvspec.\n");
  fprintf (stderr, "-n            normalize phase function \n\n");
 
  fprintf (stderr, "USAGE: %s [options] <filename>\n\n", filename);
  fprintf (stderr, "%s calculates the Legendre moments of the phase function.\n\n", PROGRAM);
  fprintf (stderr, "The input must be provided as 2-column file, containing \n");
  fprintf (stderr, "the scattering angle grid in the first column and the \n");
  fprintf (stderr, "phase function value in the second column.\n");
  fprintf (stderr, "Output are the Legendre moments. \n\n");
  fprintf (stderr, "You may specify the number of moments using the option -l.\n\n");
  fprintf (stderr, "Different scattering angle grid resolutions can be \n");
  fprintf (stderr, "chosen. For moderate forward peaks, the standard grid \n");
  fprintf (stderr, "  (1 - equidistant, 0.01 degrees step width) \n");
  fprintf (stderr, "should be sufficient. \n");
  fprintf (stderr, "For phase functions with a very strong forward peak, \n");
  fprintf (stderr, "e.g. ice particle phase functions, the finest grid resolution \n");
  fprintf (stderr, "  (2 - equidistant, 0.001 degrees step width) \n");
  fprintf (stderr, "should be specified. \n\n");
  fprintf (stderr, "Grid 3 (default) uses the scattering angle grid provided by the \n");
  fprintf (stderr, "input file; this option uses the new speedy and exact \n");
  fprintf (stderr, "method for Legendre decomposition (Buras,Dowling,Emde 201X).\n");
  fprintf (stderr, "Grids 4 and 5 have a finer resolution around\n");
  fprintf (stderr, "the forward peak. \n \n");
  fprintf (stderr, "The accuracy of the Legendre decomposition may be tested \n");
  fprintf (stderr, "using the tool phase: \n");
  fprintf (stderr, "phase -c -d -s 1 pmom_outfile.dat \n"); 
  fprintf (stderr, "\n");
}

static int get_options (int argc, char **argv,
                        char *programname,
                        char *infilename,
                        int *n_legendre,
                        int *grid_res, 
                        int *coeff,
                        int *norm)
{
  int c=0;
  char *dummy=NULL;
  
  /* save name of program */
  strncpy (programname, argv[0], FILENAME_MAX);
  
  /* set defaults */
  strcpy (infilename, "");
  *n_legendre = 0;
  *grid_res = 3;
  *norm = 0;

  while ((c=getopt (argc, argv, "hl:r:cn")) != EOF)  {
    switch(c)  {
    case 'h':  /* help */
      print_usage (programname);
      return (-1);
      break;
    case 'l':  /* number of legendre polynomials */
      *n_legendre = strtol (optarg, &dummy, 10);
      break;
    case 'r':
      *grid_res = strtol (optarg, &dummy, 10);
      break;
    case 'c':  /*Legendre coefficients */
      *coeff=1;
      break;
    case 'n':
      *norm=1;
      break; 
    default:
      print_usage (programname);
      return (-1);
    }
  }

  /* check number of remaining command line arguments */
  if (argc - optind != 1)  {
    print_usage (programname);
    return -1;
  }
  
  strncpy (infilename, argv[optind+0], FILENAME_MAX);
  
  return 0;  /* if o.k. */
}


int main(int argc, char **argv)
{
  int status=0;
  char programname[FILENAME_MAX], infilename[FILENAME_MAX]="";
  
  int n_legendre=1000;
  int grid_res=0;
  int coeff=0; 
  int norm=0;
  double *u=NULL;
  double *Xi=NULL;
  double pint=0; 
  double *ph_int=NULL; 
  double p0_1=0, p0_2=0, p1_1=0, p1_2=0, p2_1=0, p2_2=0;
  double *a0=NULL, *a1=NULL, *a2=NULL, *a3=NULL;
  double *f=NULL;
  double *s=NULL;
  int L=0;
  int l=0, i=0, N_t=0;
  double *theta=NULL;
  int conservative=0;
  double factor=1.0;
  double clockstarttime=0.0;

  status = get_options (argc, argv, programname, infilename, &n_legendre,
                        &grid_res, &coeff, &norm);

  if (status!=0)
    return status;
  
    
  /* Read phase function data */
  read_2c_file (infilename, &f, &s, &L); 
  fprintf (stderr, " ... read %d data points from %s\n", L, infilename);
  
  if (status!=0) {
    fprintf (stderr, "error %d reading %s\n", status, infilename);
    return status;
  }
  
  fprintf (stderr, "... calculate %d legendre polynomials\n", n_legendre);
  
  switch(grid_res)
    {
    case 1:
      fprintf (stderr, "... use grid resolution %d\n", grid_res);
      N_t = 18000;
      theta = (double*) malloc((N_t+1) * sizeof(f[0])); 
      theta[0] = 0.0;
      for (i=1; i <= N_t; i++)
        theta[i] = theta[i-1] + 0.01;
      break;

    case 2:
      fprintf (stderr, "... use grid resolution %d.\n", grid_res);
      N_t = 180000;
      theta = (double*) malloc((N_t+1) * sizeof(f[0]));
      theta[0] = 0.0;
      for (i=1; i <= N_t; i++)
        theta[i] = theta[i-1] + 0.001;
      break; 
     
    case 3:
      conservative=1;
      fprintf (stderr, "... use grid resolution %d.\n", grid_res);
      theta = (double*) malloc(L*sizeof(f[0]));
      N_t = L-1;
      for (i=0; i<L; i++)
        theta[i] = f[i];
      factor=-0.5;
      break;
      
    case 4:
      fprintf (stderr, "... use grid resolution %d\n", grid_res);
      N_t = 27000; 
      theta = (double*) malloc((N_t+1) * sizeof(f[0]));
      theta[0] = 0.0;
      for (i=1; i <= 10000; i++)
        theta[i] = theta[i-1] + 0.001;
      for (i=10001; i <= N_t; i++)
        theta[i] = theta[i-1] + 0.01;
      break;
      
    case 5:
      fprintf (stderr, "... use grid resolution %d.\n", grid_res);
      fprintf (stderr, "warning: this is probably inaccurate.\n");
      N_t = 2700; 
      theta = (double*) malloc((N_t+1) * sizeof(f[0]));
      theta[0] = 0.0;
      for (i=1; i <= 1000; i++)
        theta[i] = theta[i-1] + 0.01;
      for (i=1001; i <= N_t; i++)
        theta[i] = theta[i-1] + 0.1;
      break;
      
    default: 
      fprintf (stderr, "Grid resolution not valid. \n");
      fprintf (stderr, "Use resolution 1 \n");
      N_t = 180000;
      theta = (double*) malloc((N_t+1) * sizeof(f[0]));
      theta[0] = 0.0;
      for (i=1; i <= N_t; i++)
        theta[i] = theta[i-1] + 0.001;
      break;
      
    }               

  u  = calloc (N_t+1, sizeof(double));
  Xi = calloc (n_legendre, sizeof(double));
  ph_int = calloc (N_t+1, sizeof(double));

  status = linear_coeffc (f, s, L, &a0, &a1, &a2, &a3);
  if (status!=0)  {
    fprintf (stderr, "Sorry cannot do linear interpolation\n");
    fprintf (stderr, "linear_coeffc() returned status %d\n", status);
    return status;
  }

  
  if(grid_res!=3)
    {
      for (i = 0; i<=N_t; i++)
        {
        status = calc_splined_value (theta[i], &(ph_int[i]), f, L, a0, a1, a2, a3);
        /*ph_int[i] = gsl_interp_eval(interp, f, s, theta[i], acc);
          fprintf(stdout, "%f %g \n", theta[i], ph_int[i]); */
        }
    }
  else
    ph_int=s;
  
  free(a0); free(a1); free(a2); free(a3); 


  for (l=0; l<=N_t; l++)
    u[l] = cos(theta[l]*PI/180.);

  /*Check if the phase function is normalized to 2.0. */
  for (l=0; l<N_t; l++)
    pint += 0.5*(ph_int[l]+ph_int[l+1])*fabs(u[l+1] - u[l]);
  
  if (fabs(2.0 - pint) > 1e-2) {
    fprintf (stderr, "*** Warning, phase function not normalized to 2 but to %.8f\n", pint);
  }
  
  if(norm){
    fprintf (stderr, "*** Normalizing the phase function to 2.\n");
    for (l=0; l<N_t; l++)
      ph_int[l] *= 2/pint;
  }
  
  pint=0.0;
  /* Check if the phase function is normalized to 2.0. */
  for (l=0; l<N_t; l++)
    pint += 0.5*(ph_int[l]+ph_int[l+1])*fabs(u[l+1] - u[l]);
  
  if (fabs(2.0 - pint) > 1e-2)
    fprintf (stderr, "*** Warning, phase function not normalized to 2 but to %.8f\n", pint);
  
  if (conservative) {
    clockstarttime=(double)clock();
    calc_pmom (u, ph_int, N_t+1, n_legendre, Xi);
    fprintf(stderr,"cputime calc %9.3f s\n", ((double) clock()-clockstarttime)/1e6);
  }
  else {
    for (i=0; i<N_t;i++) {
      p0_1 = 1.0;
      p0_2 = 1.0;
    
      Xi[0]+=0.5*0.5*(ph_int[i]+ph_int[i+1])*fabs(u[i+1]-u[i]); 
    
      p1_1=u[i];
      p1_2=u[i+1];
    
      Xi[1]+=0.5*0.5*(p1_1*ph_int[i]+p1_2*ph_int[i+1])*fabs(u[i+1]-u[i]);  
    
      for (l=2; l<n_legendre; l++) {
	p2_1=(2*(double)l-1)/(double)l*u[i]*p1_1-((double)l-1)/(double)l*p0_1; 
	p2_2=(2*(double)l-1)/(double)l*u[i+1]*p1_2-((double)l-1)/(double)l*p0_2;
      
	Xi[l]+=0.5*0.5*(p2_1*ph_int[i]+p2_2*ph_int[i+1])*fabs(u[i+1]-u[i]);
      
	p0_1 = p1_1;
	p0_2 = p1_2;
	p1_1 = p2_1;
	p1_2 = p2_2;
      }
    }
  }

  /* Print out the Legendre coefficients */
  for (l=0; l<n_legendre; l++) {
    if (coeff)
      /* Factor to obtain Legendre coefficients,*/ 
      fprintf(stdout,  "%g " , factor*Xi[l] * (2*(double)l + 1) );
    else
      fprintf(stdout,  "%g \n " , factor*Xi[l]);
  }
  

  free (u);
  free (Xi);
  free (ph_int);

  return 0;
}  

/* BCA: these functions are in phasetable.c, they should be used from there */

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

  for (i=0; i<ntheta-1; i++) { /* ntheta is number of grid points */

    /* sum up for zeroth and first legendre coefficients with simple
       formulas */
    pmom[0] += 0.5 * ( phase[i+1] + phase[i] )
      * ( mu[i+1] - mu[i] ); 

  }

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

    pmom[1] += 0.5 * ( C - Cp )* mu2[i]
      + 1./3. * ( D - Dp ) * mu2[i] * mu[i];

    p0 = 1.0;
    p1=mu[i];

    /* move to l+2 */
    for (l=0; l<2; l++) {
      p2=legendre_polynomial(p1, p0, mu[i], l+2);
      p0 = p1;
      p1 = p2;
    }

    /* calculate coefficients using the formula from the paper (CDE 2011) */
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


static inline double legendre_polynomial (double p1, double p0,
					  double mu, int l)
{
  return ( (2*(double)l-1) * mu * p1 - ((double)l-1) * p0 ) / ((double)l); 
}
