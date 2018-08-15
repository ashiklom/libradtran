/*--------------------------------------------------------------------
 * $Id: elevation2d.c 3216 2016-05-19 08:13:59Z bernhard.mayer $
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
#include "elevation2d.h"
#include "ascii.h"
#include "mystic.h"
#if HAVE_OPENGL
  #include "GLmystic.h"
#endif

#ifndef PI
#define PI 3.14159265358979323846264338327
#endif


static int equal (double x, double y);
static int read_2D_elevation (char *filename, 
			      int *Nx, int *Ny,
			      double *delX, double *delY,
			      double ***elevation, int quiet);


/***********************************************************************************/
/* Function: setup_elevation2D                                            @62_30i@ */
/* Description:                                                                    */
/*  Initialize the 2D elevation_struct with the elevation data read from           */
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

int setup_elevation2D (char *filename,
		       int Nx, int Ny, double delX, double delY, 
		       elevation_struct *elev,
		       int visualize, int quiet)
{
  int status=0;
  int ie=0, je=0;
  double z11=0, z12=0, z21=0, z22=0, temp=0;
  double **elevation2D=NULL;
  
  
  if (!quiet)
    fprintf (stderr, " ... reading 2D elevation data from %s\n", filename);

  status = read_2D_elevation (filename,
			      &(elev->Nx), &(elev->Ny),
			      &(elev->delX), &(elev->delY),
			      &elevation2D,
			      quiet);
  
  if (status!=0) {
    fprintf (stderr, "Error %d reading %s\n", status, filename);
    return status;
  }
    
  elev->Nx -= 1;
  elev->Ny -= 1;
  
  /* check if elevation domain size is equal to 3D grid */
  
  if (!equal ((double) (elev->Nx) * elev->delX, (double) Nx * delX)) {
    fprintf (stderr, "Error, x-size of elevation area (%g) does not equal x-size of sample area (%g)\n",
	     (double) (elev->Nx) * elev->delX, (double) Nx * delX);
    return -1;
  }
    
  if (!equal ((double) (elev->Ny) * elev->delY, (double) Ny * delY)) {
    fprintf (stderr, "Error, y-size of elevation area (%g) does not equal y-size of sample area (%g)\n",
	     (double) (elev->Ny) * elev->delY, (double) Ny * delY);
    return -1;
  }
  
  /* and to make sure that both grids are really equal, recalculate delX and delY in double precision                           */
  /* (even tiny discrepancies (smaller than single precision accuracy) cause errors as confirmed in May 2016 by Philipp Gregor) */
  elev->delX = (double) Nx * delX / (double) (elev->Nx);
  elev->delY = (double) Ny * delY / (double) (elev->Ny);


  /* check periodicity in x and y directions */

  for (ie=0; ie<=elev->Nx; ie++)
    if (fabs (elevation2D[ie][0] - elevation2D[ie][elev->Ny]) > 0) {
      fprintf (stderr, "Error, elevation grid not periodic!\n");
      fprintf (stderr, "elevation2D [%d][%d] = %g, elevation2D [%d][%d] = %g\n",
	       ie, 0, elevation2D[ie][0], ie, 
	       elev->Ny, elevation2D[ie][elev->Ny]);
      return -1;
    }

  for (je=0; je<=elev->Ny; je++)
    if (fabs (elevation2D[0][je] - elevation2D[elev->Nx][je]) > 0) {
      fprintf (stderr, "Error, elevation grid not periodic!\n");
      fprintf (stderr, "elevation2D [%d][%d] = %g, elevation2D [%d][%d] = %g\n",
	       0, je, elevation2D[0][je], 
               elev->Nx, je, elevation2D[elev->Nx][je]);
      return -1;
    }

  if (fabs (elevation2D[0][0] - elevation2D[elev->Nx][elev->Ny]) > 0) {
    fprintf (stderr, "Error, elevation grid not periodic!\n");
    fprintf (stderr, "elevation2D [%d][%d] = %g, elevation2D [%d][%d] = %g\n",
	     0, je, elevation2D[0][0], 
	     elev->Nx, je, elevation2D[elev->Nx][elev->Ny]);
    return -1;
  }
    
  #if HAVE_OPENGL
  if (visualize)  {
    GLmystic_calloc_shared (GLMYSTIC_ELEVATION, (elev->Nx+1) * (elev->Ny+1));
 
    GLmystic_setdomain (GLMYSTIC_ELEVATION,
			0, 0, 0, 
			(double) (elev->Nx) * elev->delX, 
			(double) (elev->Ny) * elev->delY, 
			0);

    GLmystic_write_elevation (elevation2D, elev->Nx, elev->Ny);
  }
  #endif


  /* calculate coefficients for bilinear interpolation of elevation */
  elev->surf = calloc ((size_t) elev->Nx, sizeof(surface *));
  
  /* maximum surface elevation */
  elev->surfmax = 0.0;
  
  for (ie=0; ie<elev->Nx; ie++) {
    (elev->surf)[ie] = calloc ((size_t) elev->Ny, sizeof(surface));
    
    for (je=0; je<elev->Ny; je++) {
      z11 = elevation2D[ie][je];
      z21 = elevation2D[ie+1][je];
      z12 = elevation2D[ie][je+1];
      z22 = elevation2D[ie+1][je+1];
      
      (elev->surf)[ie][je].a = (z21-z11) / elev->delX;
      (elev->surf)[ie][je].b = (z12-z11) / elev->delY;
      (elev->surf)[ie][je].c = (z22+z11-z12-z21) / elev->delY / elev->delX;
      (elev->surf)[ie][je].d = z11;
      
      /* calculate maximum surface elevation for this pixel */
      temp = z11;
      
      if (z21>temp)
	temp = z21;
      
      if (z12>temp)
	temp = z12;
      
      if (z22>temp)
	temp = z22;
      
      (elev->surf)[ie][je].upper = temp;
      if (temp > elev->surfmax)
	elev->surfmax = temp;
    }
  }
  

  for (ie=0; ie<=elev->Nx; ie++)
    free (elevation2D[ie]);
  
  free (elevation2D);

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


/*******************************************************/
/* Read data from a 2D elevation file.                 */
/*******************************************************/

static int read_2D_elevation (char *filename, 
			      int *Nx, int *Ny,
			      double *delX, double *delY,
			      double ***elevation, int quiet)
{
  int i=0, status=0;
  int ix=0, iy=0;

  double **value=NULL;
  int rows=0, min_columns=0, max_columns=0;


#if HAVE_LIBNETCDF
  struct stat buf;

  size_t start[1] = {0};
 
  size_t tstart[2] = {0,0};
  size_t tcount[2] = {0,0};
  
  size_t n=0;
 
  int    ncid   =0; 
  int    idd_nx =0, idd_ny =0;
  int    id_delX=0, id_delY=0, id_type=0;

  double *temp=NULL;
#endif



  /* try open as netCDF file */
  
#if HAVE_LIBNETCDF
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status==NC_NOERR) {

    if (!quiet)
      fprintf (stderr, " ... reading 2D elevation data from netCDF file %s\n", filename);

    /* determine file date and stop if the file was created before */
    /* September 7, 2007 where the format was changed              */
    stat (filename, &buf);
    
    gmtime(&buf.st_mtime);
    
    if (buf.st_mtime<1189187082) {
      fprintf (stderr, "\n");
      fprintf (stderr, "*** Error %s was last changed %s", filename, ctime(&buf.st_mtime));
      fprintf (stderr, "*** and the MYSTIC elevation2D convention has been changed on\n");
      fprintf (stderr, "*** September 7, 2007. It is likely that you use the \n");
      fprintf (stderr, "*** old convention! Please check the documentation!\n");
      fprintf (stderr, "\n");
      
      return -1;
    }
    
    fprintf (stderr, "\n");
    fprintf (stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf (stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf (stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf (stderr, "!!! Careful! The netcdf file format for 2D elevation  !!!\n");
    fprintf (stderr, "!!! was changed on September 7, 2007. In particular,  !!!\n");
    fprintf (stderr, "!!! x and y were interchanged so that an elevation    !!!\n");
    fprintf (stderr, "!!! file viewed with ncview should look like a map    !!!\n");
    fprintf (stderr, "!!! now, with (x,y)=(1,1) in the lower left corner.   !!!\n");
    fprintf (stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf (stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf (stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf (stderr, "\n");

    /* get dimension id for "nx" */
    status = nc_inq_dimid (ncid, "nx", &idd_nx);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading Nx from %s\n", status, filename);
      return status;
    }
    
    /* get dimension length for "nx" */
    status = nc_inq_dimlen (ncid, idd_nx, &n);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading Nx from %s\n", status, filename);
      return status;
    }

    *Nx=n;

    /* get dimension id for "ny" */
    status = nc_inq_dimid (ncid, "ny", &idd_ny);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading Ny from %s\n", status, filename);
      return status;
    }
    
    /* get dimension length for "ny" */
    status = nc_inq_dimlen (ncid, idd_ny, &n);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading Ny from %s\n", status, filename);
      return status;
    }

    *Ny=n;

    /* get variable id for "delX" */
    status = nc_inq_varid (ncid, "dx", &id_delX);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading delX from %s\n", 
	       status, filename);
      return status;
    }

    /* read "delX" */
    status = nc_get_var1_double (ncid, id_delX, start, delX);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading delX from %s\n", 
	       status, filename);
      return status;
    }


    /* get variable id for "delY" */
    status = nc_inq_varid (ncid, "dy", &id_delY);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading dy from %s\n", 
	       status, filename);
      return status;
    }

    /* read "delY" */
    status = nc_get_var1_double (ncid, id_delY, start, delY);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading dy from %s\n", 
	       status, filename);
      return status;
    }


    *delX *= 1000.0;
    *delY *= 1000.0;


    /* allocate memory for 2D elevation field */
    *elevation = calloc((size_t) *Nx, sizeof(double *));

    for (ix=0; ix<*Nx; ix++)
      (*elevation)[ix] = calloc((size_t) *Ny, sizeof(double));
    

    /* read data */

    /* get variable id for "elevation" */
    status = nc_inq_varid (ncid, "elevation", &id_type);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading elevation from %s\n", status, filename);
      return status;
    }

    temp = calloc (*Nx, sizeof(double));

    /* read type */
    tstart[1] = 0;   /* start with first row          */
    tcount[1] = *Nx; /* read *Nx elements for each iy */

    tcount[0] = 1;  /* read one line at a time */

    for (iy=0; iy<*Ny; iy++) {
      tstart[0] = iy;

      status = nc_get_vara_double (ncid, id_type, tstart, tcount, temp);
      if (status!=NC_NOERR) {
	fprintf (stderr, "Error %d reading elevation from %s\n", status, filename);
	return status;
      }

      for (ix=0; ix<*Nx; ix++)
	(*elevation)[ix][iy] = temp[ix] * 1000.0;
    }
    nc_close (ncid);

    free(temp);
    
    if (!quiet)
      fprintf (stderr, " ... read %d data points from %s\n", 
	       (*Nx)*(*Ny), filename);
  }
  else {
    if (!quiet)
      fprintf (stderr, " ... %s not in netCDF format, trying to open as ASCII\n", 
	       filename);
#endif


    /* read file; ASCII_file2double is a waste of memory in this */
    /* case, but it is very convenient because comments and      */
    /* everything are handled correctly.                         */
    
    status = ASCII_file2double (filename,
				&rows, &max_columns, &min_columns, 
				&value);
    
    if (status!=0) {
      fprintf (stderr, "Error %d reading file %s\n", status, filename);
      return status;
    }

    if (max_columns<4) {
      fprintf (stderr, "Error: found less than four columns in %s\n", 
	       filename);
      return -1;
    }
    
    if (rows<1) {
      if (!quiet)
	fprintf (stderr, " ... found only one row in %s\n", filename);
      return -1;
    }
    
    
    /* 1st line, number of grid points */
    *Nx = (int) (value[0][0] + 0.5);
    *Ny = (int) (value[0][1] + 0.5);
    
    
    /* 1st line, horizontal distances, convert from km to m */
    *delX = value[0][2] * 1000.0;
    *delY = value[0][3] * 1000.0;
    
    if (!quiet) {
      fprintf (stderr, " ... x-distance: %g\n", *delX);
      fprintf (stderr, " ... y-distance: %g\n", *delY);
    }
    
    /* allocate memory for elevation */
    *elevation = calloc((size_t) *Nx, sizeof(double *));
    
    for (ix=0; ix<*Nx; ix++)
      (*elevation)[ix] = calloc((size_t) *Ny, sizeof(double));
    
    
    for (i=1; i<rows; i++) {
      ix = (int) (value[i][0] + 0.5) - 1;
      iy = (int) (value[i][1] + 0.5) - 1;
      
      if (ix<0 || ix>=*Nx) {
	fprintf (stderr, "Error, ix = %d out of bounds in line %d\n", ix+1, i+1);
	return -1;
      }
      
      if (iy<0 || iy>=*Ny) {
	fprintf (stderr, "Error, iy = %d out of bounds in line %d\n", iy+1, i+1);
	return -1;
      }
      
      /* copy data to result arrays, convert from km to m */
      (*elevation)[ix][iy] = value[i][2] * 1000.0;
      
    }
    
    /* free memory */
    (void) ASCII_free_double (value, rows);

    if (!quiet)
      fprintf (stderr, " ... read %d data points from %s\n", 
	       rows-1, filename);

#if HAVE_LIBNETCDF
  }
#endif
    
  return 0;
}



/***********************************************************************************/
/* Function: free_elevation                                               @62_30i@ */
/* Description:                                                                    */
/*  Free memory of struct elevation_struct.                                        */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/

void free_elevation (elevation_struct *elev) 
{
  int ie=0;

  for (ie=0; ie<elev->Nx; ie++)
    free(elev->surf[ie]);
    
  free(elev->surf);

  free(elev);
}


