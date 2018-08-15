/*--------------------------------------------------------------------
 * $Id: uvspecrandom.h 3246 2016-08-23 16:32:20Z bernhard.mayer $
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

#ifndef __uvspecrandom_h
#define __uvspecrandom_h

#if defined (__cplusplus)
extern "C" {
#endif

double uvspec_random ();
double uvspec_random_gauss (double sigma);
int init_uvspec_random (int *randomseed, int need_mt19937, int quiet);

#if defined (__cplusplus)
}
#endif

#endif
