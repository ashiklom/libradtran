/*--------------------------------------------------------------------
 * $Id: temperature2d.c 2931 2013-06-10 11:28:53Z robert.buras $
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
#include <float.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "uvspec.h"
#include "wcloud3d.h"
#include "ascii.h"
#include "mystic.h"


#ifndef PI
#define PI 3.14159265358979323846264338327
#endif


static int equal (double x, double y);


/***********************************************************************************/
/* Function: setup_temperature2D                                          @62_30i@ */
/* Description:                                                                    */
/*  Initialize the 2D temperature_struct with the temperature data read from       */
/*  filename.                                                                      */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int setup_temperature2D (char *filename,
			 int Nx, int Ny, double delX, double delY, 
			 surftemp_struct *surftemp,
			 int ixmin, int ixmax, int iymin, int iymax,
			 int quiet)
{
  int ia=0, status=0;
  
  if (!quiet)
    fprintf (stderr, " ... reading 2D surface temperature data from %s\n", filename);
  
  /* we use read_2D_albedo because file formats are identical for temperature and albedo  */
  status = read_2D_albedo (filename,
			   &(surftemp->Nx), &(surftemp->Ny),
			   &(surftemp->delX), &(surftemp->delY),
			   &(surftemp->temp2D),
			   ixmin, ixmax, iymin, iymax,
			   quiet);
  
  if (status!=0) {
    fprintf (stderr, "Error %d reading 2D surface temperature from %s\n", status, filename);
    return status;
  }
  
  if (!equal ((double) (surftemp->Nx) * surftemp->delX, (double) Nx * delX)) {
    fprintf (stderr, "Error, x-size of 2D surface temperature area (%g) does not equal x-size of cloud area (%g)\n",
	     (double) (surftemp->Nx) * surftemp->delX, (double) Nx * delX);
    return -1;
  }
  
  if (!equal ((double) (surftemp->Ny) * surftemp->delY, (double) Ny * delY)) {
    fprintf (stderr, "Error, y-size of surface temperature area (%g) does not equal y-size of cloud area (%g)\n",
	     (double) (surftemp->Ny) * surftemp->delY, (double) Ny * delY);
    return -1;
  }
  
  /* allocate ,memory for planck emission */
  surftemp->plkavg2D = calloc (surftemp->Nx, sizeof(double *));
  for (ia=0; ia<surftemp->Nx; ia++)
    surftemp->plkavg2D[ia] = calloc (surftemp->Ny, sizeof(double));

  return 0;
}


/***********************************************************************************/
/* Function: equal                                                        @62_30i@ */
/* Description:                                                                    */
/*  Two numbers equal?                                                             */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

static int equal (double x, double y)
{
  if ((float) x != (float) y)
    return 0;

  return 1;
}

