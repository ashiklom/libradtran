/*--------------------------------------------------------------------
 * $Id: rodents.c 2623 2011-12-23 10:52:38Z robert.buras $
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

/*--------------------------------------------------------------------
 * RODENTS: ROberts' Delta-EddingtoN Two-Stream
 *
 * by Robert Buras
 *
 * development in joint cooperation with Bernhard Mayer
 *
 * implementation partly by Ulrike Wissmeier
 *
 * Programmed using the book by Zdunkowski, Trautmann and Bott,
 * "Radiation in the Atmosphere", chapters 6.1-6.4
 *
 * Note that in the 2007 print there are misprints in Eqs. 6.50 and 6.88
 * Also, we implemented thermal emission differently than in chapter 6.5
 * of the book.
 *
 *--------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "solver.h"
#include "rodents.h"
#include "ascii.h"
#include "bandec.h"

#ifndef PI
#define PI 3.14159265358979323846264338327
#endif

static inline int err_out ( char *out, int status );
static void free_input(rodents_input_struct input);

/*************************************************************/
/* rodents, former tsm                                       */
/*************************************************************/

int rodents ( /* INPUT */
	      int     nlyr,
	      float  *dtau,
	      float  *omega_0,
	      float  *g, 
	      float  *f,
	      int     delta_method,
	      float  *temper,
	      float   wvnm_low,
	      float   wvnm_upp,
	      int     usrtau,
	      int     nzout,
	      float  *utau,
	      double  S_0,
	      double  mu_0,
	      double  albedo,
	      float   btemp,
	      int     useplanck,
	      /* INPUT NECESSARY FOR TIPA DIR */
	      int     tipa,
	      float  *tausol,
	      /* OUTPUT */
	      float  *fldn,
	      float  *flup,
	      float  *fldir,
	      float  *uavg)
{
  int status=0;
  rodents_input_struct input;
  rodents_output_struct output;

  if (delta_method!=RODENTS_DELTA_METHOD_OFF && tipa==2) {
    fprintf(stderr,"Error! Tipa DIR and delta scaling do not work together!\n");
    fprintf(stderr,"Exiting...\n");
    return -1;
  }

  output.e_minus       = calloc((size_t) nlyr+1, sizeof(float));
  output.e_plus        = calloc((size_t) nlyr+1, sizeof(float));
  output.s_direct      = calloc((size_t) nlyr+1, sizeof(float));
  output.s_direct_unsc = calloc((size_t) nlyr+1, sizeof(float));
    
  /* initialization */

  status = setup_factors ( nlyr, dtau, omega_0, g, f, delta_method,
			   temper, wvnm_low, wvnm_upp,
			   usrtau, nzout, utau,
			   S_0, mu_0, albedo, btemp, useplanck,
			   tipa, tausol,
			   &input);
  if (status)
    return err_out("Error %d returned from setup_factors!\n",status);

  /* solvers */
  status = multi_layer_solver (input, output);
  if (status)
    return err_out("Error %d returned from multi_layer_solver!\n",status);

  status = process_output (input, output, fldn, flup, fldir, uavg);
  if (status)
    return err_out("Error %d returned from process_output!\n",status);

  /* free memory */  

  free_input (input);

  free(output.e_minus);
  free(output.e_plus);
  free(output.s_direct);
  free(output.s_direct_unsc);
  
  return 0;
}


/*************************************************************/
/* setup_factors                                             */
/* number of flops per layer: (needs update)                 */
/* add : 14                                                  */
/* mult: 23                                                  */
/* div :  5                                                  */
/* exp :  2                                                  */
/* sqrt:  1                                                  */
/*************************************************************/

int setup_factors ( /* INPUT */
		    int     nlyr,
		    float  *dtau, 
		    float  *omega_0, 
		    float  *g, 
		    float  *f, 
		    int     delta_method,
		    float  *temper,
		    float   wvnm_low,
		    float   wvnm_upp,
		    int     usrtau,
		    int     nzout,
		    float  *utau,
		    double  S_0,
		    double  mu_0,
		    double  albedo,
		    float   btemp,
		    int     useplanck,
		    int     tipa,
		    float  *tausol,
		    /* OUTPUT */
		    rodents_input_struct *input)
{
  int i=0, iu=0;
  double bscr=0, b_mmu_0 =0;
  double dtau_d=0.0, B_1_d=0.0, g_d=0.0, omega_0_d=0.0;
  double B_1 = 0.0, B_0 = 0.0;
  float *B_atm=NULL;
  double alpha_1=0.0, alpha_2=0.0, alpha_3=0.0, alpha_4=0.0;
  double alpha_5=0.0, alpha_6=0.0;
  double c1=0.0, c2=0.0, c3=0.0;
  double A=0.0, lambda=0.0;
  double mu_0_inv=0.0;
  double den1=0.0, exp1 = 0.0, term1=0.0;
  double gamma_0=0.0, gamma_1=0.0, gamma_2=0.0;
  double pi_div_mu_bar=0.0;
  double *tau=NULL;
  double dusrtau1=0.0, dusrtau1_d=0.0, dusrtau2=0.0, dusrtau2_d=0.0;
  double nor=0.0, scr=0.0;

  double epsilon = 1e-6; /* RB-UW replace by general definition BCA */

  input->nlyr     = nlyr;
  input->S_0      = S_0;
  input->mu_0     = mu_0;
  input->albedo   = albedo;
  input->usrtau   = usrtau;
  input->nzout    = nzout;

  /* define whether solar is on or off */
  if ( S_0 > 0.0 && mu_0 > 0.0 )
    input->solar = 1;
  else
    input->solar = 0;

  /* prepare thermal emission */
  input->thermal  = 0;
  if (useplanck) {
    if (wvnm_low > 5000.0) { /* no thermal for lambda < 2mu */
      input->B_ground = 0.0;
    }    
    else {
      input->thermal=1;
      B_atm = calloc((size_t) nlyr+1, sizeof(float));
      for (i=0;i<nlyr+1;i++)
	F77_FUNC (cplkavg, CPLKAVG) (&(wvnm_low), &(wvnm_upp), &(temper[i]),
				     &(B_atm[i]));
      F77_FUNC (cplkavg, CPLKAVG) (&(wvnm_low), &(wvnm_upp), &(btemp),
				   &(input->B_ground));
    }
  }

  /* allocate memory */
  input->a11 = calloc( nlyr, sizeof(double) );
  input->a12 = calloc( nlyr, sizeof(double) );

  input->b1 = calloc( nlyr, sizeof(double) );
  input->b2 = calloc( nlyr, sizeof(double) );

  input->a13 = calloc( nlyr, sizeof(double) );
  input->a23 = calloc( nlyr, sizeof(double) );
  input->a33 = calloc( nlyr, sizeof(double) );

  /* unscaled extinction for direct radiation */
  input->a33_unsc = calloc( nlyr, sizeof(double) );

  /* output level specific solution parameters */
  if (input->usrtau) {
    /* allocate */
    input->ilyr_usr     = calloc( input->nzout, sizeof(int) );
    input->a11_usr      = calloc( input->nzout, sizeof(double) );
    input->a12_usr      = calloc( input->nzout, sizeof(double) );
    input->a13_usr      = calloc( input->nzout, sizeof(double) );
    input->a21_usr      = calloc( input->nzout, sizeof(double) );
    input->a22_usr      = calloc( input->nzout, sizeof(double) );
    input->a23_usr      = calloc( input->nzout, sizeof(double) );
    input->a33_usr      = calloc( input->nzout, sizeof(double) );
    input->a33_usr_unsc = calloc( input->nzout, sizeof(double) );
    input->b1_usr       = calloc( input->nzout, sizeof(double) );
    input->b2_usr       = calloc( input->nzout, sizeof(double) );

    tau    = calloc( input->nlyr+1, sizeof(double) );
    tau[0] = 0.0;
    for (i = 0; i < input->nlyr; i++)
      tau[i+1] = tau[i] + dtau[i];

    /* detect usrtau layers */
    for (iu = 0; iu < input->nzout; iu++) {
      for (i = 0; i < input->nlyr-1; i++)
	if (utau[iu] >= tau[i] && utau[iu] <= tau[i+1])
	  break;
      input->ilyr_usr[iu] = i;
    }

    iu=input->nzout-1;
  }

  /* choose delta scaling factor */
  if (delta_method==RODENTS_DELTA_METHOD_HG)
    for (i=0;i<nlyr;i++)
      f[i] = g[i] * g[i];
  if (delta_method==RODENTS_DELTA_METHOD_OFF)
    for (i=0;i<nlyr;i++)
      f[i] = 0.0;

  /* calculate parameters */

  mu_0_inv = 1./mu_0;

  pi_div_mu_bar = PI * 2.;

  for (i=0;i<nlyr;i++) {

    /* catch some singularities */ 
    if ( omega_0[i] >= 1.0 )
      omega_0[i] = 1.0 - epsilon;

    if ( omega_0[i] * g[i] == 1.0 )
      omega_0[i] *= 1.0 - epsilon;

    /* define linear fit to Planck function */
    if (useplanck) {
      if (dtau[i] > 0.01) {
	B_0 = B_atm[i];
	B_1 = ( B_atm[i+1] - B_atm[i] ) / dtau[i];
      }
      else {
	B_0 = 0.5 * ( B_atm[i] + B_atm[i+1] );
	B_1 = 0.0;
      }
    }
    else {
      B_0=0.0;
      B_1=0.0;
    }
    
    /* do delta scaling */
    dtau_d = dtau[i] * ( 1. - omega_0[i] * f[i] );
    B_1_d  = B_1 / ( 1. - omega_0[i] * f[i] );
    g_d    = ( g[i] - f[i] ) / ( 1. - f[i] );
    omega_0_d = omega_0[i] * ( 1. - f[i] )
      / ( 1. - f[i] * omega_0[i] );

    b_mmu_0 = 0.5 - 0.75 * g_d * mu_0;

    /* Eq. 6.64 */
    bscr = 0.5 - 0.375 * g_d;
    alpha_1 = 2. * ( 1. - omega_0_d * ( 1. - bscr ) ) - 0.25;
    alpha_2 = 2. * omega_0_d * bscr - 0.25;

    /* Eq. 6.69 */
    lambda = sqrt ( alpha_1 * alpha_1 - alpha_2 * alpha_2 );

    if ( lambda * dtau_d > 1e2 ) {
      input->a11[i] = 0.0;
      input->a12[i] = ( alpha_1 - lambda ) / alpha_2;
    }
    else { 
      /* exp1, term1 needed for A, a_matrix */
      exp1  = exp( lambda * dtau_d );
      term1 = alpha_2 / ( alpha_1 - lambda ) * exp1;

      /* Eq. 6.82 */
      A = 1.0 / ( term1 - 1. / term1 );

      /* Eq. 6.88 and own calculations */
      input->a11[i] = A * 2.0 * lambda / alpha_2;
      input->a12[i] = A * ( exp1 - 1. / exp1 );
    }

    /* thermal only */
    if (input->thermal) {
      gamma_0 = pi_div_mu_bar * ( 1. - omega_0_d );
      gamma_1 = gamma_0 / ( alpha_1 - alpha_2);
      gamma_2 = gamma_0 / ( lambda * lambda );
      c1 = gamma_1 * ( B_0 + B_1_d * dtau_d );
      c2 = gamma_2 * B_1_d;
      c3 = gamma_1 * B_0;
      input->b1[i] = - input->a11[i] * ( c1 + c2 ) - input->a12[i] * ( c3 - c2 ) + c2 + c3;
      input->b2[i] = - input->a12[i] * ( c1 + c2 ) - input->a11[i] * ( c3 - c2 ) + c1 - c2;
    }
    
    /* solar only */
    if (input->solar) {
      /* den1 needed for alpha_5 and alpha_6 */
      den1 = 1. / ( mu_0_inv * mu_0_inv - lambda * lambda );

      /* Eqs. 6.64 and 6.77 */
      alpha_3 = - omega_0_d * b_mmu_0;
      alpha_4 = omega_0_d + alpha_3;
      alpha_5 = ( ( alpha_1 - mu_0_inv ) * alpha_3 - alpha_2 * alpha_4 ) * den1;
      alpha_6 = ( alpha_2 * alpha_3 - ( alpha_1 + mu_0_inv ) * alpha_4 ) * den1;

      input->a33[i]      = exp ( - dtau_d  * mu_0_inv );   
      input->a33_unsc[i] = exp ( - dtau[i] * mu_0_inv );   

      input->a13[i] =
	+ alpha_5 * ( 1.0 - input->a33[i] * input->a11[i] )
	- alpha_6 * input->a12[i];

      input->a23[i] =
	- alpha_5 * input->a33[i] * input->a12[i]
	+ alpha_6 * ( input->a33[i] - input->a11[i] );

      if (tipa==2) { /* BCA the "2" should be TIPA_DIR, but TIPA_DIR is defined in src/uvspec.h and not in libsrc_c */
	if (i==0)
	  input->a33[i] = exp ( - tausol[i] * mu_0_inv );
	else
	  input->a33[i] = exp ( - ( tausol[i] - tausol[i-1] ) * mu_0_inv );

	/* for TIPA dir, no scaling possible */
	input->a33_unsc[i] = input->a33[i];
      }
    }

    if (input->usrtau && iu>=0) {
      while ( input->ilyr_usr[iu] == i ) {

	/* calculate parameters for sub-layer above utau */
	dusrtau2   = utau[iu] - tau[i];
	dusrtau2_d = dusrtau2 * ( 1. - omega_0[i] * f[i] );

	if ( lambda * dusrtau2_d > 1e2 ) {
	  input->a22_usr[iu] = 0.0;
	  input->a21_usr[iu] = ( alpha_1 - lambda ) / alpha_2;
	}
	else { 
	  /* exp1, term1 needed for A, a_matrix */
	  exp1  = exp( lambda * dusrtau2_d );
	  term1 = alpha_2 / ( alpha_1 - lambda ) * exp1;

	  /* Eq. 6.82 */
	  A = 1.0 / ( term1 - 1. / term1 );

	  /* Eq. 6.88 and own calculations */
	  input->a22_usr[iu] = A * 2.0 * lambda / alpha_2;
	  input->a21_usr[iu] = A * ( exp1 - 1. / exp1 );
	}

	/* thermal only */
	if (input->thermal) {
	  c1 = gamma_1 * ( B_0 + B_1_d * dusrtau2_d );
	  input->b2_usr[iu] = - input->a21_usr[iu] * ( c1 + c2 ) - input->a22_usr[iu] * ( c3 - c2 ) + c1 - c2;
	}

	/* solar only */
	if (input->solar) {
	  input->a33_usr[iu]      = exp ( - dusrtau2_d * mu_0_inv );   
	  input->a33_usr_unsc[iu] = exp ( - dusrtau2   * mu_0_inv );   

	  input->a23_usr[iu] =
	    - alpha_5 *   input->a33_usr[iu] * input->a21_usr[iu]
	    + alpha_6 * ( input->a33_usr[iu] - input->a22_usr[iu] );

	  if (tipa==2) { /* BCA the "2" should be TIPA_DIR, but TIPA_DIR is defined in src/uvspec.h and not in libsrc_c */
	    /* not yet implemented XXX */
	    return -1;
	  }
	}

	/* calculate parameters for sub-layer above utau */
	dusrtau1   = - utau[iu] + tau[i+1];
	dusrtau1_d = dusrtau1 * ( 1. - omega_0[i] * f[i] );

	if ( lambda * dusrtau1_d > 1e2 ) {
	  input->a11_usr[iu] = 0.0;
	  input->a12_usr[iu] = ( alpha_1 - lambda ) / alpha_2;
	}
	else { 
	  /* exp1, term1 needed for A, a_matrix */
	  exp1  = exp( lambda * dusrtau1_d );
	  term1 = alpha_2 / ( alpha_1 - lambda ) * exp1;

	  /* Eq. 6.82 */
	  A = 1.0 / ( term1 - 1. / term1 );

	  /* Eq. 6.88 and own calculations */
	  input->a11_usr[iu] = A * 2.0 * lambda / alpha_2;
	  input->a12_usr[iu] = A * ( exp1 - 1. / exp1 );
	}

	/* thermal only */
	if (input->thermal) {
	  c1 = gamma_1 * ( B_0 + B_1_d * ( dusrtau2_d + dusrtau1_d ) );
	  c3 = gamma_1 * ( B_0 + B_1_d * dusrtau2_d );
	  input->b1_usr[iu] = - input->a11_usr[iu] * ( c1 + c2 ) - input->a12_usr[iu] * ( c3 - c2 ) + c2 + c3;
	}

	/* solar only */
	if (input->solar) {
	  input->a13_usr[iu] =
	    + alpha_5 * ( 1.0 - exp ( - dusrtau1_d * mu_0_inv )
			  * input->a11_usr[iu] )
	    - alpha_6 * input->a12_usr[iu];

	  if (tipa==2) { /* BCA the "2" should be TIPA_DIR, but TIPA_DIR is defined in src/uvspec.h and not in libsrc_c */
	    /* not yet implemented XXX */
	    return -1;
	  }
	}

	/* now use these to calculate the "true" parameters */

	nor = 1. / ( 1. - input->a12_usr[iu] * input->a21_usr[iu] );

	scr               = ( input->b1_usr[iu] + input->a12_usr[iu] * input->b2_usr[iu] ) * nor;
	input->b2_usr[iu] = ( input->b2_usr[iu] + input->a21_usr[iu] * input->b1_usr[iu] ) * nor;
	input->b1_usr[iu] = scr;

	scr                = ( input->a13_usr[iu] * input->a33_usr[iu] + input->a12_usr[iu] * input->a23_usr[iu] ) * nor;
	input->a23_usr[iu] = ( input->a23_usr[iu] + input->a21_usr[iu] * input->a13_usr[iu] * input->a33_usr[iu] ) * nor;
	input->a13_usr[iu] = scr;

	input->a11_usr[iu] *= nor;
	input->a22_usr[iu] *= nor;

	input->a12_usr[iu] *= input->a22_usr[iu];
	input->a21_usr[iu] *= input->a11_usr[iu];

	iu--;
	if (iu==-1)
	  break;
      } /* end while */
    } /* end if usrtau */
  } /* end do ilyr */

  if (input->thermal)
    free(B_atm);
  if (input->usrtau)
    free(tau);

  return 0;
}


/*************************************************************/
/* multi_layer_solver                                        */
/*************************************************************/

int multi_layer_solver (rodents_input_struct   input,
			rodents_output_struct  output)
{
  double **matr, **mlu, *b;
  long *indx;
  int i=0, n=2*input.nlyr+2;

  /* calloc memory */
  matr = calloc(n,sizeof(double *));
  for (i=0;i<n;i++)
    matr[i] = calloc(5,sizeof(double));

  mlu = calloc(n,sizeof(double *));
  for (i=0;i<n;i++)
    mlu[i] = calloc(2,sizeof(double));

  indx = calloc(n,sizeof(long));
  b    = calloc(n,sizeof(double));

  /******** RBUW start of important stuff ********/

  /* solve direct radiation */
  if (input.S_0 != 0.0) {
    solve_direct(input.a33, input.a33_unsc, input.nlyr, input.S_0,
		 output.s_direct, output.s_direct_unsc);
  }

  /* setup matrix and vector for diffuse radiation */
  setup_equations (input, output.s_direct, matr, b);

  /* solve diffuse radiation */
  bandec(matr, n, mlu, indx);
  banbks(matr, n, mlu, indx, b);
  
  /* put result */
  cp_to_result(b, input.nlyr, output.e_minus, output.e_plus);

  /******** RBUW end of important stuff ********/

  for (i=0;i<n;i++) {
    free(mlu[i]);
    free(matr[i]);
  }
  free(b);
  free(indx);
  free(mlu);
  free(matr);

  return 0;
}


/*************************************************************/
/* solve_direct                                              */
/*************************************************************/

void solve_direct(double *a33, double *a33_unsc, int nlyr, double S_0,
		  float *S, float *S_unsc)
{
  int i=0;
  
  S[0]      = S_0;
  S_unsc[0] = S_0;
  
  for (i=0;i<nlyr;i++) {
    S[i+1]      = a33[i]      * S[i];
    S_unsc[i+1] = a33_unsc[i] * S_unsc[i];
  }
}


/*************************************************************/
/* setup_equations                                           */
/*************************************************************/

void setup_equations ( /* input */
		       rodents_input_struct  input,
		       float                *s_direct,
		       /* output */
		       double              **matr,
		       double               *rhs )
{
  int i=0, i2=0;

  /* the matrix: */

  matr[0][2] = -1.0;

  for (i=0;i<input.nlyr;i++) {
    i2=2*i+1;
    matr[i2  ][2] = -1.0;
    matr[i2  ][4] = input.a11[i];
    matr[i2  ][1] = input.a12[i];
    matr[i2+1][2] = -1.0;
    matr[i2+1][0] = input.a11[i];
    matr[i2+1][3] = input.a12[i];
  }

  i=input.nlyr;
  i2=2*i+1;
  matr[i2  ][2] = -1.0;
  matr[i2  ][1] = input.albedo;

  /* the RHS */

  if (input.thermal) {
    for (i=0;i<input.nlyr;i++) {
      i2=2*i+1;
      rhs[i2  ] = - input.b1[i];
      rhs[i2+1] = - input.b2[i];
    }

    i=input.nlyr;
    i2=2*i;
    rhs[i2+1] = - (1.0 - input.albedo) * PI * input.B_ground;
  }

  if (input.solar) {
    for (i=0;i<input.nlyr;i++) {
      i2=2*i+1;
      rhs[i2  ] += - input.a13[i] * s_direct[i];
      rhs[i2+1] += - input.a23[i] * s_direct[i];
    }

    i=input.nlyr;
    i2=2*i;
    rhs[i2+1] += - input.albedo * s_direct[i] * input.mu_0;
  }
}


/*************************************************************/
/* cp_to_result                                              */
/*************************************************************/

void cp_to_result(double *b,
		  int nlyr,
		  float *e_minus,
		  float *e_plus)
{
  int i=0;
  for (i=0;i<nlyr+1;i++) {
    e_minus[i] = b[2*i];
    e_plus [i] = b[2*i+1];
  }
}

/*************************************************************/
/* process_output                                            */
/*************************************************************/

int process_output ( /* input */
		     rodents_input_struct   input,
		     rodents_output_struct  output,
		     /* output */
		     float                 *fldn,
		     float                 *flup,
		     float                 *fldir,
		     float                 *uavg)
{
  int i=0, iu=0;
  double fldir_unsc=0.0;

  if (input.usrtau) {
    for (iu=0;iu<input.nzout;iu++) {
      i = input.ilyr_usr[iu];

      /* calculate fluxes at zouts */
      fldn [iu]  = input.a21_usr     [iu] * output.e_plus  [i+1]
                 + input.a22_usr     [iu] * output.e_minus [i  ]
                 + input.a23_usr     [iu] * output.s_direct[i  ]
	         + input.b2_usr      [iu];

      flup [iu]  = input.a11_usr     [iu] * output.e_plus  [i+1]
                 + input.a12_usr     [iu] * output.e_minus [i  ]
                 + input.a13_usr     [iu] * output.s_direct[i  ]
	         + input.b1_usr      [iu];

      fldir[iu]  = input.a33_usr     [iu] * output.s_direct[i  ];

      fldir_unsc = input.a33_usr_unsc[iu] * output.s_direct_unsc[i  ];

      /* unscale results */
      fldn [iu]  = fldn [iu] + ( fldir[iu] - fldir_unsc ) * input.mu_0;
      fldir[iu]  = fldir_unsc * input.mu_0;

      uavg [iu]  = NAN;
    }
  }
  else {
    /* output all levels */
    for (iu=0;iu<input.nlyr+1;iu++)  {
      fldn [iu] = output.e_minus [iu]
	+ ( output.s_direct[iu] - output.s_direct_unsc[iu] ) * input.mu_0;
      flup [iu] = output.e_plus  [iu];
      fldir[iu] = output.s_direct_unsc[iu] * input.mu_0;
      uavg [iu] = NAN;
    }
  }

  /* old code before correct interpolation. Might be needed for tipa==2 */
  /*
  double epsilon = 1e-6;
  for (i=0;i<nlyr+1;i++) 
    for (lu=0;lu<nzout;lu++)
      if (zout_sur[lu] <= z[i]+epsilon && zout_sur[lu] >= z[i]-epsilon)
      {
	fldn [lu] = output.e_minus [i]
	  + ( output.s_direct[i] - output.s_direct_unsc[i] ) * mu_0;
	flup [lu] = output.e_plus  [i];
	fldir[lu] = output.s_direct_unsc[i] * mu_0;
      }
  */

  return 0;
}

/*************************************************************/
/* free_input                                                */
/*************************************************************/

static void free_input(rodents_input_struct input)
{
  free(input.a11);
  free(input.a12);
  free(input.b1);
  free(input.b2);
  free(input.a13);
  free(input.a23);
  free(input.a33);
  free(input.a33_unsc);
  if (input.usrtau) {
    free(input.ilyr_usr);
    free(input.a11_usr);
    free(input.a12_usr);
    free(input.a13_usr);
    free(input.a21_usr);
    free(input.a22_usr);
    free(input.a23_usr);
    free(input.a33_usr);
    free(input.a33_usr_unsc);
    free(input.b1_usr);
    free(input.b2_usr);
  }
}

/*************************************************************/
/* err_out                                                   */
/*************************************************************/

static inline int err_out ( char *out, int status )
{
  fprintf(stderr,out,status);
  return status;
}
