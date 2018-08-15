/*--------------------------------------------------------------------
 * $Id: numeric.h 2623 2011-12-23 10:52:38Z robert.buras $
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

#ifndef __numeric_h
#define __numeric_h

#if defined (__cplusplus)
extern "C" {
#endif

#include "equation.h"    /* solve linear equation systems             */
#include "function.h"    /* some analytical and statistical functions */
#include "spl.h"         /* interpolating and approximating splines   */
#include "locate.h"      /* location of a value in an array           */
#include "linear.h"      /* linear interpolation                      */
#include "cnv.h"         /* convolution                               */
#include "regress.h"     /* some linear and nonlinear regression      */
#include "integrat.h"    /* numerical integration                     */
#include "table.h"       /* reading/interpolation of tables           */

#if defined (__cplusplus)
}
#endif

#endif









