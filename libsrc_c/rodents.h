/*--------------------------------------------------------------------
 * $Id: rodents.h 2623 2011-12-23 10:52:38Z robert.buras $
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

#ifndef __rodents_h
#define __rodents_h

#define RODENTS_DELTA_METHOD_OFF 0
#define RODENTS_DELTA_METHOD_HG  1
#define RODENTS_DELTA_METHOD_ON  2

#if defined (__cplusplus)
extern "C" {
#endif

typedef struct {
  int    thermal;
  int    solar;
  int    nlyr;
  double S_0;
  double mu_0;
  double albedo;
  int    usrtau;
  int    nzout;
  float  B_ground;
  double *a11, *a12, *a13, *a23, *a33, *a33_unsc;
  int    *ilyr_usr;
  double *a11_usr, *a12_usr, *a13_usr, *a21_usr, *a22_usr, *a23_usr, *a33_usr, *a33_usr_unsc;
  double *b1, *b2, *b1_usr, *b2_usr;
} rodents_input_struct;

typedef struct {
  float *e_minus, *e_plus, *s_direct, *s_direct_unsc;
  //  float *e_minus_usr, *e_plus_usr, *s_direct_usr, *s_direct_unsc_usr;
} rodents_output_struct;


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
	      float  *uavg);

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
		    rodents_input_struct *input);

int multi_layer_solver ( rodents_input_struct   input,
			 rodents_output_struct  output);

void solve_direct(double *a33, double *a33_u, int nlyr, double S_0,
		  float *S, float *S_u);

void setup_equations( /* input */
		      rodents_input_struct  input,
		      float                *s_direct,
		      /* output */
		      double              **matr,
		      double               *rhs );

void cp_to_result(double *b,
		  int nlyr,
		  float *e_minus,
		  float *e_plus);

int process_output ( /* input */
		     rodents_input_struct   input,
		     rodents_output_struct  output,
		     /* output */
		     float                 *fldn,
		     float                 *flup,
		     float                 *fldir,
		     float                 *uavg);

#endif
