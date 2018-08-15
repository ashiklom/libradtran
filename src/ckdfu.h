/*--------------------------------------------------------------------
 * $Id: ckdfu.h 3276 2017-07-04 14:16:36Z Claudia.Emde $
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

#ifndef __ckdfu_h
#define __ckdfu_h

#include <stdio.h>
#include "uvspec.h"

/* int ckdfu (float umco2, float umch4, float umn2o,  */
/* 	   float umf11, float umf12, float umf22,  */
/* 	   float *rhoair, float *rhoh2o, float *rhoo3,  */
/* 	   float *pres, float *temp, float *davg, float *z, int nlev, int h2o_cont, */
/* 	   ck_profile *ck); */
int ckdfu (float umco2, float umch4, float umn2o, 
	   float umf11, float umf12, float umf22, 
	   float ***temper, float ***press, float *z, 
	   float ****dens, float ****dens_avg, int nlev,
	   int Nx, int Ny,
	   //float *rhoair, float *rhoh2o, float *rhoo3, 
	   //float *pres, float *temp, float *davg, float *z, int nlev,
	   int h2o_cont,
	   ck_profile *ck);


int ckdfucld (float *z, float *reff, float *lwc, int nlev,
	      float **tau, float **gg, float **ssa);

int ckdfuice (float *z, float *reff, float *lwc, int nlev,
	      float **tau, float **gg, float **ssa, int unscaled);

int ckdfuray (int ib, int ig, float u0, int nlev, float *tau);

#endif
