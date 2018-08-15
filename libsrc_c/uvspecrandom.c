/*--------------------------------------------------------------------
 * $Id: uvspecrandom.c 3246 2016-08-23 16:32:20Z bernhard.mayer $
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

#include <unistd.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

#if HAVE_LIBGSL
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_randist.h>
#endif

/* global variable: random number generator */
#if HAVE_LIBGSL 
  gsl_rng *uvspecrng; 
#endif


/*****************************************************************/
/* uvspec random number generator; either the highly-recommended */
/* MT19937 if GSL is available, or random() otherwise.           */
/*****************************************************************/

double uvspec_random ()
{
  #if HAVE_LIBGSL 
    return gsl_rng_uniform (uvspecrng);
  #else
    return (double) random();
  #endif
}

double uvspec_random_gauss (double sigma)
{
  #if HAVE_LIBGSL 
    return gsl_ran_gaussian (uvspecrng, sigma);
  #else
    return (double) random();
  #endif
}

int init_uvspec_random (int *randomseed, int need_mt19937, int quiet)
{
  int rseed=0;

/* global variable, defined in uvspec_lex.l */
#if HAVE_LIBGSL 
  /*  extern gsl_rng *uvspecrng; */  /* global random number generator */
  const gsl_rng_type *T;
#endif


  /* always initialize standard rng */
  srandom(rseed);

  if (*randomseed>0)
    rseed = *randomseed;
  else
    rseed = (int) time(NULL) + (int) getpid();

#if HAVE_LIBGSL 
  T = gsl_rng_mt19937;
  gsl_rng_default_seed = rseed;
  uvspecrng = gsl_rng_alloc (T);
  if (!quiet && need_mt19937)
    fprintf (stderr, " ... using %s random number generator\n", gsl_rng_name (uvspecrng));
#else
  if (need_mt19937) {
    fprintf (stderr, "Error, GSL not found. The MT19937 random number generator is\n");
    fprintf (stderr, "required for MYSTIC simulations. Please install GSL and recompile!\n");
    return -1;
  }
#endif

  return rseed;
}
