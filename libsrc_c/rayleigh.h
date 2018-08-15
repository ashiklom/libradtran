/*--------------------------------------------------------------------
 * $Id: rayleigh.h 2839 2013-01-31 15:13:42Z svn-kylling $
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

#ifndef __rayleigh_h
#define __rayleigh_h

int crs_rayleigh_nicolet(float *lambda, int n_lambda, float **crs );

int crs_rayleigh_bodhaine(float *lambda, int n_lambda, float mixing_ratio_co2, 
			  float depol_user, float **depol, float **crs );

int crs_rayleigh_bodhaine29(float *lambda, int n_lambda, float **crs );

int crs_rayleigh_penndorf(float *lambda, int n_lambda, float mixing_ratio_co2, 
			  float depol_user, float **depol, float **crs );

#endif
