/*--------------------------------------------------------------------
 * $Id: ipa.h 2623 2011-12-23 10:52:38Z robert.buras $
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

#ifndef __ipa_h
#define __ipa_h

#include <stdio.h>
#include "uvspec.h"

int setup_ipa ( input_struct   input,
		output_struct *output );

int read_caoth_ipa ( input_struct       input,
		     caoth_inp_struct   input_caoth,
		     float              altitude,
		     wl_out_struct      wl_out,
		     atm_out_struct     atm_out,
		     /* output */
		     int               *nipa,
		     float            **ipaweight,
		     caoth_out_struct  *output_caoth,
		     caoth_out_struct **output_caoth_ipa,
		     cf_out_struct     *cf_out );

int compare_levels_ipa_df ( float *z1,
			    int    n1,
			    float *z2,
			    int    n2 );
#endif



