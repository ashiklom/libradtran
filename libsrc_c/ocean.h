/*--------------------------------------------------------------------
 * $Id: ocean.h 2623 2011-12-23 10:52:38Z robert.buras $
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

#ifndef __ocean_h
#define __ocean_h

#if defined (__cplusplus)
extern "C" {
#endif
  
#include "f77-uscore.h"
  

  /* prototypes */

  float ocean_brdf (float wvnmlo, float wvnmhi, float mu, float mup, float dphi,
		    float u10, float pcl, float xsal, int firstcall);
  
  
  /* Fortran-type oceanc function which builds a search tree etc. */
  void F77_FCN  (oceanc) (float *wvnmlo, float *wvnmhi, float *mu, float *mup, float *dphi, 
			float *u10, float *pcl, float *xsal, int *firstcall, float *bdref);

  /* simple C wrapper for the standard oceabrdf() function, without search tree etc */
  double oceabrdfc (float wvnmlo, float wvnmhi, 
                    float mu, float mup, float dphi, 
                    float u10, float pcl, float xsal, float paw);
  

  /* simple C wrapper for the bpdf function by Mishchenko */
  int bpdf_tsang(double wind_speed, float wavelength, double mu_inc, double phi_inc,
                 double mu_esc, double phi_esc, double ***refl_mat);
    
  typedef struct {
    double re;
    double im;
  } ref_complex;
  
  
#if defined (__cplusplus)
}
#endif

#endif
