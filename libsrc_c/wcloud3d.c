/*--------------------------------------------------------------------
 * $Id: wcloud3d.c 3234 2016-07-08 15:04:27Z Claudia.Emde $
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
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <float.h>
#include "ascii.h"
#include "wcloud3d.h"

#if HAVE_LIBNETCDF
#include <netcdf.h>
#endif

#ifndef PI
#define PI 3.14159265358979323846264338327
#endif

/* Maximum number of characters per line */
#define MAX_LENGTH_OF_LINE     65536


/*******************************************************/
/* Read data from a 3D cloud description file.         */
/*******************************************************/
 
int read_3D_caoth ( char      *filename, 
		    int       *Nx,
		    int       *Ny,
		    int       *Nz,
		    double    *delX,
		    double    *delY,
		    float    **z, 
		    float  ****lwc,
		    float  ****reff,
		    float  ****ext,
		    float  ****g1, 
		    float  ****g2, 
		    float  ****ff, 
		    float  ****f, 
		    float  ****ssa, 
		    float  ****dscale, 
		    float     *rmin,
		    float     *rmax,
		    int       *cldproperties,
		    int      **threed,
		    int        ixmin,
		    int        ixmax,
		    int        iymin,
		    int        iymax,
		    int        quiet )
{
  int i=0, number=0, status=0;
  int ix=0, iy=0, iz=0;

  int rows=0, min_columns=0, max_columns=0, max_length=0;

  char line[MAX_LENGTH_OF_LINE+1]="";
  char *string=NULL;

  FILE *file=NULL;
  char **array=NULL;
  char *dummy=NULL;

#if HAVE_LIBNETCDF
  size_t tstart[3] = {0,0,0};
  size_t tcount[3] = {0,0,0};
  int ixstart=0, iystart=0;

  size_t n=0;
  int read_netcdf=0;

  int get_ext=0, get_g1=0, get_ssa=0, get_lwc=0, get_reff=0;

  int    ncid   =0, idd_nx =0, idd_ny =0, idd_nz=0;
  int    id_z=0, id_ext=0, id_lwc=0, id_reff=0, id_g1=0, id_ssa=0;

  float *temp=NULL;
#endif

#if HAVE_LIBNETCDF
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status==NC_NOERR) {
    read_netcdf=1;

    if (!quiet)
      fprintf (stderr, " ... reading Cloud data from netCDF file %s\n", filename);

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

    /* get dimension id for "nz" */
    status = nc_inq_dimid (ncid, "nz", &idd_nz);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading Nz from %s\n", status, filename);
      return status;
    }
    
    /* get dimension length for "nz" */
    status = nc_inq_dimlen (ncid, idd_nz, &n);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading Nz from %s\n", status, filename);
      return status;
    }

    *Nz=n;

    /* get dimension length for "cldproperties" */
    status = nc_get_att_int (ncid, NC_GLOBAL, "cldproperties", cldproperties);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading cldproperties from %s\n", status, filename);
      return status;
    }

    /* read "delX" */
    status = nc_get_att_double (ncid, NC_GLOBAL, "dx", delX);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading delX from %s\n", 
	       status, filename);
      return status;
    }

    /* read "delY" */
    status = nc_get_att_double (ncid, NC_GLOBAL, "dy", delY);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading dy from %s\n", 
	       status, filename);
      return status;
    }

    /* *delX *= 1000.0; */
    /* *delY *= 1000.0; */

    if (ixmin!=-1 && ixmax!=-1 && iymin!=-1 && iymax!=-1) {
      if ( ixmax>=*Nx || ixmin<0 || iymax>=*Ny || iymin<0 ) {
	fprintf(stderr,"Error! You have specified a region with `atmos_region` (%d,%d,%d,%d) which is larger than the file %s\n",ixmin,ixmax,iymin,iymax,filename);
	return -1;
      }

      *Nx = ixmax - ixmin + 1;
      *Ny = iymax - iymin + 1;
      ixstart = ixmin;
      iystart = iymin;
    }
    else {
      ixstart=0;
      iystart=0;
    }

  }
  else {
    if (!quiet)
      fprintf (stderr, " ... %s not in netCDF format, trying to open as ASCII\n", 
	       filename);
#endif
    
    /* check input file */
    status =  ASCII_checkfile (filename, 
			       &rows,
			       &min_columns,
			       &max_columns,
			       &max_length);
  
    if (status!=0) {
      fprintf (stderr, "Error %d reading %s\n",
	       status, filename);
      return status;
    }
    
    if (rows<2) {
      fprintf (stderr, "Error: found less than two rows in %s\n", filename);
      return -1;
    }

    if (min_columns<3) {
      fprintf (stderr, "Error: found less than three columns in %s\n", filename);
      return -1;
    }
  
    /* reset string to the beginning of the memory block */
    string = line;
  
    /* open file */
    if ( (file = fopen(filename, "r")) == NULL)  
      return ASCIIFILE_NOT_FOUND;
  

    /* first line */
    number=0;
    while (number==0) {
      dummy=fgets (string, MAX_LENGTH_OF_LINE, file);
      status = ASCII_parsestring (string, &array, &number);

      if (status!=0) {
	fprintf (stderr, "Error %d reading 1st line of %s\n", status, filename);
	return status;
      }
    }

    if (number<3) {
      fprintf (stderr, "Error: expected at least 3 values in 1st line, found only %d\n", number);
      return -1;
    }
  
    /* first 3 values are number of cells */
    *Nx = strtol (array[0], &dummy, 0);
    *Ny = strtol (array[1], NULL, 0);
    *Nz = strtol (array[2], NULL, 0);

    if (number>3)
      *cldproperties = strtol (array[3], NULL, 0);
    
    if (number>4) {
      fprintf (stderr, "*** WARNING: Found a 5th number in the first line of %s.\n", filename);
      fprintf (stderr, "*** Possibly you use the old format. Please be aware that wspec\n");
      fprintf (stderr, "*** doesn't exist anymore and the 5th number is ignored!\n");
    }
  
    free (array);

    /* 2nd line */
    number=0;
    while (number==0) {
      dummy=fgets (string, MAX_LENGTH_OF_LINE, file);
      status = ASCII_parsestring (string, &array, &number);

      if (status!=0) {
	fprintf (stderr, "Error %d reading 1st line of %s\n", status, filename);
	return status;
      }
    }

    if (number!=*Nz+1+2) {
      fprintf (stderr, "Error: found %d z levels, expected %d\n",
	       number-2, *Nz+1);
      return -1;
    }

    /* horizontal grid sizes, in km */
    *delX = strtod (array[0], NULL);
    *delY = strtod (array[1], NULL);

#if HAVE_LIBNETCDF
  }
#endif

  if (!quiet) {
    fprintf (stderr, " ... reading %d x %d x %d data points from %s\n",
	     *Nx, *Ny, *Nz, filename);
    fprintf (stderr, " ... cldproperties = %d\n", *cldproperties);
    fprintf (stderr, " ... cloud grid size in x direction: %g km\n", *delX);
    fprintf (stderr, " ... cloud grid size in y direction: %g km\n", *delY);
  }

  /* allocate memory for altitude levels */
  *z = calloc ((size_t) *Nz+1, sizeof(float));
  if (*z==NULL) {
    fprintf (stderr, "Error allocating memory for z\n");
    return -1;
  }
	     
  *threed = calloc ((size_t) *Nz, sizeof(int));
  if (*threed==NULL) {
    fprintf (stderr, "Error allocating memory for threed\n");
    return -1;
  }
	     
  /* allocate memory for optical depth, asymmetry factor, and single scattering albedo */
  *ext    = calloc((size_t) *Nz, sizeof(float **));
  if (*ext==NULL) {
    fprintf (stderr, "Error allocating memory for ext\n");
    return -1;
  }
  *g1     = calloc((size_t) *Nz, sizeof(float **));
  if (*g1==NULL) {
    fprintf (stderr, "Error allocating memory for g1\n");
    return -1;
  }
  *g2     = calloc((size_t) *Nz, sizeof(float **));
  if (*g2==NULL) {
    fprintf (stderr, "Error allocating memory for g2\n");
    return -1;
  }
  *ff     = calloc((size_t) *Nz, sizeof(float **));
  if (*ff==NULL) {
    fprintf (stderr, "Error allocating memory for ff\n");
    return -1;
  }
  *ssa    = calloc((size_t) *Nz, sizeof(float **));
  if (*ssa==NULL) {
    fprintf (stderr, "Error allocating memory for ssa\n");
    return -1;
  }
  *f      = calloc((size_t) *Nz, sizeof(float **));
  if (*f==NULL) {
    fprintf (stderr, "Error allocating memory for f\n");
    return -1;
  }
  *dscale = calloc((size_t) *Nz, sizeof(float **));
  if (*dscale==NULL) {
    fprintf (stderr, "Error allocating memory for dscale\n");
    return -1;
  }
  *lwc    = calloc((size_t) *Nz, sizeof(float **));
  if (*lwc==NULL) {
    fprintf (stderr, "Error allocating memory for lwc\n");
    return -1;
  }
  *reff   = calloc((size_t) *Nz, sizeof(float **));
  if (*reff==NULL) {
    fprintf (stderr, "Error allocating memory for reff\n");
    return -1;
  }
  
  for (iz=0; iz<*Nz; iz++) {
    (*ext)    [iz] = calloc((size_t) *Nx, sizeof(float *));
    if ((*ext)[iz]==NULL) {
      fprintf (stderr, "Error allocating memory for ext[%d]\n",iz);
      return -1;
    }
    (*g1)     [iz] = calloc((size_t) *Nx, sizeof(float *));
    if ((*g1)[iz]==NULL) {
      fprintf (stderr, "Error allocating memory for g1[%d]\n",iz);
      return -1;
    }
    (*g2)     [iz] = calloc((size_t) *Nx, sizeof(float *));
    if ((*g2)[iz]==NULL) {
      fprintf (stderr, "Error allocating memory for g2[%d]\n",iz);
      return -1;
    }
    (*ff)     [iz] = calloc((size_t) *Nx, sizeof(float *));
    if ((*ff)[iz]==NULL) {
      fprintf (stderr, "Error allocating memory for ff[%d]\n",iz);
      return -1;
    }
    (*ssa)    [iz] = calloc((size_t) *Nx, sizeof(float *));
    if ((*ssa)[iz]==NULL) {
      fprintf (stderr, "Error allocating memory for ssa[%d]\n",iz);
      return -1;
    }
    (*f)      [iz] = calloc((size_t) *Nx, sizeof(float *));
    if ((*f)[iz]==NULL) {
      fprintf (stderr, "Error allocating memory for g[%d]\n",iz);
      return -1;
    }
    (*dscale) [iz] = calloc((size_t) *Nx, sizeof(float *));
    if ((*dscale)[iz]==NULL) {
      fprintf (stderr, "Error allocating memory for dscale[%d]\n",iz);
      return -1;
    }
    (*lwc)    [iz] = calloc((size_t) *Nx, sizeof(float *));
    if ((*lwc)[iz]==NULL) {
      fprintf (stderr, "Error allocating memory for lwc[%d]\n",iz);
      return -1;
    }
    (*reff)   [iz] = calloc((size_t) *Nx, sizeof(float *));
    if ((*reff)[iz]==NULL) {
      fprintf (stderr, "Error allocating memory for reff[%d]\n",iz);
      return -1;
    }

    for (ix=0; ix<*Nx; ix++) {
      (*ext)    [iz][ix] = calloc((size_t) *Ny, sizeof(float)); 
      if ((*ext)[iz][ix]==NULL) {
	fprintf (stderr, "Error allocating memory for ext[%d][%d]\n",iz,ix);
	return -1;
      }
      (*g1)     [iz][ix] = calloc((size_t) *Ny, sizeof(float)); 
      if ((*g1)[iz][ix]==NULL) {
	fprintf (stderr, "Error allocating memory for g1[%d][%d]\n",iz,ix);
	return -1;
      }
      (*g2)     [iz][ix] = calloc((size_t) *Ny, sizeof(float)); 
      if ((*g2)[iz][ix]==NULL) {
	fprintf (stderr, "Error allocating memory for g2[%d][%d]\n",iz,ix);
	return -1;
      }
      (*ff)     [iz][ix] = calloc((size_t) *Ny, sizeof(float)); 
      if ((*ff)[iz][ix]==NULL) {
	fprintf (stderr, "Error allocating memory for ff[%d][%d]\n",iz,ix);
	return -1;
      }
      (*ssa)    [iz][ix] = calloc((size_t) *Ny, sizeof(float)); 
      if ((*ssa)[iz][ix]==NULL) {
	fprintf (stderr, "Error allocating memory for ssa[%d][%d]\n",iz,ix);
	return -1;
      }
      (*f)      [iz][ix] = calloc((size_t) *Ny, sizeof(float)); 
      if ((*f)[iz][ix]==NULL) {
	fprintf (stderr, "Error allocating memory for f[%d][%d]\n",iz,ix);
	return -1;
      }
      (*dscale) [iz][ix] = calloc((size_t) *Ny, sizeof(float)); 
      if ((*dscale)[iz][ix]==NULL) {
	fprintf (stderr, "Error allocating memory for dscale[%d][%d]\n",iz,ix);
	return -1;
      }
      (*lwc)    [iz][ix] = calloc((size_t) *Ny, sizeof(float));  
      if ((*lwc)[iz][ix]==NULL) {
	fprintf (stderr, "Error allocating memory for lwc[%d][%d]\n",iz,ix);
	return -1;
      }
      (*reff)   [iz][ix] = calloc((size_t) *Ny, sizeof(float));  
      if ((*reff)[iz][ix]==NULL) {
	fprintf (stderr, "Error allocating memory for reff[%d][%d]\n",iz,ix);
	return -1;
      }

	/* initialize forward HG fraction with 1 */
	for (iy=0; iy<*Ny; iy++)
	  (*ff)[iz][ix][iy] = 1.0; 
      }
  }


  /* initialize rmin, rmax */
  switch (*cldproperties) {
  case CLD_OPTPROP:
    *rmin=0; *rmax=0;  /* rmin, rmax not used */
    break;

  case CLD_EXTREFF:
  case CLD_LWCREFF:
    *rmin=+FLT_MAX;
    *rmax=-FLT_MAX;
    break;

  default:
    fprintf (stderr, "Error, unknown cloud properties %d in %s\n", *cldproperties, filename);
    return -1;
  }

#if HAVE_LIBNETCDF
  if (read_netcdf==1) {

    /* read data */

    /* set which parameters to read */
    switch (*cldproperties) {
    case CLD_OPTPROP:
      get_ext=1;
      get_g1=1;
      get_ssa=1;

      for (iy=0; iy<*Ny; iy++)
	for (ix=0; ix<*Nx; ix++)
	  for (iz=0; iz<*Nz; iz++)
	    (*ff)[iz][ix][iy]=1.0;

      break;
    case CLD_EXTREFF:
      get_ext=1;
      get_reff=1;
      break;
    case CLD_LWCREFF:
      get_lwc=1;
      get_reff=1;
      break;
    default:
      fprintf (stderr, "Error, unknown cloud properties %d in %s\n", *cldproperties, filename);
      return -1;
    }

    status = nc_inq_varid (ncid, "z", &id_z);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading z from %s\n", status, filename);
	return status;
    }

    tstart[0]=0;
    tcount[0]=*Nz+1;
    status = nc_get_vara_float (ncid, id_z, tstart, tcount, *z);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading z from %s\n", status, filename);
      return status;
    }


    /* get variable id for "ext" */
    if (get_ext) {
      status = nc_inq_varid (ncid, "ext", &id_ext);
      if (status!=NC_NOERR) {
	fprintf (stderr, "Error %d reading ext from %s\n", status, filename);
	return status;
      }
    }

    /* get variable id for "lwc" */
    if (get_lwc) {
      status = nc_inq_varid (ncid, "lwc", &id_lwc);
      if (status!=NC_NOERR) {
	fprintf (stderr, "Error %d reading lwc from %s\n", status, filename);
	return status;
      }
    }

    /* get variable id for "reff" */
    if (get_reff) {
      status = nc_inq_varid (ncid, "reff", &id_reff);
      if (status!=NC_NOERR) {
	fprintf (stderr, "Error %d reading reff from %s\n", status, filename);
	return status;
      }
    }

    /* get variable id for "g1" */
    if (get_g1) {
      status = nc_inq_varid (ncid, "g1", &id_g1);
      if (status!=NC_NOERR) {
	fprintf (stderr, "Error %d reading g1 from %s\n", status, filename);
	return status;
      }
    }

    /* get variable id for "ssa" */
    if (get_ssa) {
      status = nc_inq_varid (ncid, "ssa", &id_ssa);
      if (status!=NC_NOERR) {
	fprintf (stderr, "Error %d reading ssa from %s\n", status, filename);
	return status;
      }
    }

    temp = calloc (*Nz, sizeof(float));

    /* read type */
    tstart[2] = 0;   /* start with first z element       */
    tcount[2] = *Nz; /* read *Nz elements for each iy,ix */

    tcount[0] = 1;   /* read one y element at a time  */
    tcount[1] = 1;   /* read one x element at a time  */

    for (iy=0; iy<*Ny; iy++) {
      tstart[0] = iy + iystart;

      if (!quiet)
	fprintf (stderr, ".");

      for (ix=0; ix<*Nx; ix++) {
	tstart[1] = ix + ixstart;

	if (get_ext) {
	  status = nc_get_vara_float (ncid, id_ext, tstart, tcount, temp);
	  if (status!=NC_NOERR) {
	    fprintf (stderr, "Error %d reading iso from %s\n", status, filename);
	    return status;
	  }

	  for (iz=0; iz<*Nz; iz++) {
	    (*ext)[iz][ix][iy]=temp[iz] / 1000.0; /* convert from km-1 to m-1 */
	    if (temp[iz]!=0.0)
	      (*threed)[iz]=1;
	  }
	}

	if (get_lwc) {
	  status = nc_get_vara_float (ncid, id_lwc, tstart, tcount, temp);
	  if (status!=NC_NOERR) {
	    fprintf (stderr, "Error %d reading iso from %s\n", status, filename);
	    return status;
	  }

	  for (iz=0; iz<*Nz; iz++) {
	    (*lwc)[iz][ix][iy]=temp[iz];
	    if (temp[iz]!=0.0)
	      (*threed)[iz]=1;
	  }
	}

	if (get_reff) {
	  status = nc_get_vara_float (ncid, id_reff, tstart, tcount, temp);
	  if (status!=NC_NOERR) {
	    fprintf (stderr, "Error %d reading iso from %s\n", status, filename);
	    return status;
	  }

	  for (iz=0; iz<*Nz; iz++) {
	    (*reff)[iz][ix][iy]=temp[iz];
	    if ((*reff)[iz][ix][iy] < *rmin && (*reff)[iz][ix][iy] != 0.0)
	      *rmin = (*reff)[iz][ix][iy];
	    if ((*reff)[iz][ix][iy] > *rmax)
	      *rmax = (*reff)[iz][ix][iy];
	  }
	}

	if (get_g1) {
	  status = nc_get_vara_float (ncid, id_g1, tstart, tcount, temp);
	  if (status!=NC_NOERR) {
	    fprintf (stderr, "Error %d reading iso from %s\n", status, filename);
	    return status;
	  }

	  for (iz=0; iz<*Nz; iz++)
	    (*g1)[iz][ix][iy]=temp[iz];
	}

	if (get_ssa) {
	  status = nc_get_vara_float (ncid, id_ssa, tstart, tcount, temp);
	  if (status!=NC_NOERR) {
	    fprintf (stderr, "Error %d reading iso from %s\n", status, filename);
	    return status;
	  }

	  for (iz=0; iz<*Nz; iz++)
	    (*ssa)[iz][ix][iy]=temp[iz];
	}

      }
    }
    nc_close (ncid);

    free(temp);

    if (!quiet)
      fprintf (stderr, "\n");

  }
  else {
#endif

    /* rest of line 2 */
    /* set altitude levels */
    for (i=2; i<=*Nz+2; i++)
      (*z)[i-2] = (float) strtod(array[i], NULL);

    free(array);

    for (i=2; i<rows; i++) {
      number=0;
      while (number==0) {
	dummy=fgets (string, MAX_LENGTH_OF_LINE, file);
	status = ASCII_parsestring (string, &array, &number);
    
	if (status!=0) {
	  fprintf (stderr, "Error %d reading 1st line of %s\n", status, filename);
	  return status;
	}
      }

      switch (*cldproperties) {
      case CLD_OPTPROP:
	if (number<6) {
	  fprintf (stderr, "Error: found %d columns in row %d, expected 6\n", number, i);
	  return -1;
	}
	break;

      case CLD_EXTREFF:
	if (number<5) {
	  fprintf (stderr, "Error: found %d columns in row %d, expected 5\n", number, i);
	  return -1;
	}
	break;

      case CLD_LWCREFF:
	if (number<5) {
	  fprintf (stderr, "Error: found %d columns in row %d, expected 5\n", number, i);
	  return -1;
	}
	break;

      default:
	fprintf (stderr, "Error, unknown cloud properties %d in %s\n", *cldproperties, filename);
	return -1;
      }
      
      
      /* cell indices */
      ix = strtol(array[0], NULL, 0) - 1;
      iy = strtol(array[1], NULL, 0) - 1;
      iz = strtol(array[2], NULL, 0) - 1;

      if (ix<0 || ix>=*Nx) {
	fprintf (stderr, "Error, ix = %d out of bounds in line %d\n", ix+1, i+1);
	return -1;
      }

      if (iy<0 || iy>=*Ny) {
	fprintf (stderr, "Error, iy = %d out of bounds in line %d\n", iy+1, i+1);
	return -1;
      }

      if (iz<0 || iz>=*Nz) {
	fprintf (stderr, "Error, iz = %d out of bounds in line %d\n", iz+1, i+1);
	return -1;
      }


      /* copy data to result arrays */
      switch (*cldproperties) {
      case CLD_OPTPROP:
	(*ext)[iz][ix][iy] = strtod (array[3], NULL)   / 1000.0;  /* convert from km-1 to m-1 */
	(*g1) [iz][ix][iy] = strtod (array[4], NULL);
	(*g2) [iz][ix][iy] = 0;
	(*ff) [iz][ix][iy] = 1.0;
	(*ssa)[iz][ix][iy] = strtod (array[5], NULL);
	break;

      case CLD_EXTREFF:
	(*ext) [iz][ix][iy] = strtod (array[3], NULL) / 1000.0;  /* convert from km-1 to m-1 */
	(*reff)[iz][ix][iy] = strtod (array[4], NULL);
      
	if ((*reff)[iz][ix][iy] < *rmin)  *rmin = (*reff)[iz][ix][iy];
	if ((*reff)[iz][ix][iy] > *rmax)  *rmax = (*reff)[iz][ix][iy];
  
	break;

      case CLD_LWCREFF:
	(*lwc) [iz][ix][iy] = strtod (array[3], NULL);
	(*reff)[iz][ix][iy] = strtod (array[4], NULL);

	if ((*reff)[iz][ix][iy] < *rmin)  *rmin = (*reff)[iz][ix][iy];
	if ((*reff)[iz][ix][iy] > *rmax)  *rmax = (*reff)[iz][ix][iy];
	break;

      default:
	fprintf (stderr, "Error, unknown cloud properties %d in %s\n", *cldproperties, filename);
	return -1;
      }

      /* set flag indicating that this layer is 3D */
      (*threed)[iz] = 1;

      free(array);
    }


    /* close input file */
    fclose (file);

    if (!quiet)
      fprintf (stderr, " ... read %d data points from %s\n", 
	       rows-2, filename);

#if HAVE_LIBNETCDF
  }
#endif

  return 0;
}





/*******************************************************/
/* Read data from a 2D albedo description file.        */
/*******************************************************/

int read_2D_albedo ( char       *filename, 
		     int       *Nx,
		     int       *Ny,
		     double    *delX,
		     double    *delY,
		     double  ***albedo,
		     int        ixmin,
		     int        ixmax,
		     int        iymin,
		     int        iymax,
		     int        quiet )
{
  int i=0, status=0;
  int ix=0, iy=0;

  double **value=NULL;
  int rows=0, min_columns=0, max_columns=0;

#if HAVE_LIBNETCDF
  int ixstart=0, iystart=0;

  struct stat buf;
/* CP: Commented out since it's not needed */
/*  struct tm *tmstruct=NULL;*/

  size_t tstart[2] = {0,0};
  size_t tcount[2] = {0,0};
 
  size_t n=0;

  int    ncid   =0, idd_nx =0, idd_ny =0;
  int    id_type=0;

  float *temp=NULL;
#endif


  /* try to open as netCDF file */
  
#if HAVE_LIBNETCDF
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status==NC_NOERR) {

    if (!quiet)
      fprintf (stderr, " ... reading 2D albedo data from netCDF file %s\n", filename);

    /* determine file date and stop if the file was created before */
    /* September 7, 2007 where the format was changed              */
    stat (filename, &buf);
    
/* CP: Commented out since it's not needed */
    /*tmstruct=gmtime(&buf.st_mtime);*/
    
    if (buf.st_mtime<1189187082) {
      fprintf (stderr, "\n");
      fprintf (stderr, "*** Error %s was last changed %s", filename, ctime(&buf.st_mtime));
      fprintf (stderr, "*** and the MYSTIC albedo2D convention has been changed on\n");
      fprintf (stderr, "*** September 7, 2007. It is likely that you use the \n");
      fprintf (stderr, "*** old convention! Please check the documentation!\n");
      fprintf (stderr, "\n");
      
      return -1;
    }
    
    fprintf (stderr, "\n");
    fprintf (stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf (stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf (stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    fprintf (stderr, "!!! Careful! The netcdf file format for 2D albedo was !!!\n");
    fprintf (stderr, "!!! changed on September 7, 2007. In particular,      !!!\n");
    fprintf (stderr, "!!! x and y were interchanged so that an albedo file  !!!\n");
    fprintf (stderr, "!!! viewed with ncview should look like a map now,    !!!\n");
    fprintf (stderr, "!!! with (x,y)=(1,1) in the lower left corner.        !!!\n");
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

    *Ny = n;

    /* read "delX" */
    status = nc_get_att_double (ncid, NC_GLOBAL, "dx", delX);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading delX from %s\n", 
	       status, filename);
      return status;
    }

    /* read "delY" */
    status = nc_get_att_double (ncid, NC_GLOBAL, "dy", delY);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading dy from %s\n", 
	       status, filename);
      return status;
    }

    *delX *= 1000.0;
    *delY *= 1000.0;


    if (ixmin!=-1 && ixmax!=-1 && iymin!=-1 && iymax!=-1) {
      if ( ixmax>=*Nx || ixmin<0 || iymax>=*Ny || iymin<0 ) {
	fprintf(stderr,"Error! You have specified a region with `atmos_region` (%d,%d,%d,%d) which is larger than the file %s\n",ixmin,ixmax,iymin,iymax,filename);
	return -1;
      }

      *Nx = ixmax - ixmin + 1;
      *Ny = iymax - iymin + 1;
      ixstart = ixmin;
      iystart = iymin;
    }
    else {
      ixstart=0;
      iystart=0;
    }

    /* allocate memory for 2D albedo field */
    *albedo = calloc((size_t) *Nx, sizeof(double *));

    for (ix=0; ix<*Nx; ix++)
      (*albedo)[ix] = calloc((size_t) *Ny, sizeof(double));
    

    /* read data */

    /* get variable id for "albedo" */
    status = nc_inq_varid (ncid, "albedo", &id_type);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading albedo from %s\n", status, filename);
      return status;
    }

    temp = calloc (*Nx, sizeof(float));
    
    /* read type */
    tstart[1] = ixstart;   /* start with first row          */
    tcount[1] = *Nx; /* read *Nx elements for each iy */

    tcount[0] = 1;  /* read one line at a time */

    for (iy=0; iy<*Ny; iy++) {
      tstart[0] = iy + iystart;

      status = nc_get_vara_float (ncid, id_type, tstart, tcount, temp);
      if (status!=NC_NOERR) {
	fprintf (stderr, "Error %d reading albedo from %s\n", status, filename);
	return status;
      }

      for (ix=0; ix<*Nx; ix++)
	(*albedo)[ix][iy] = temp[ix];
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
      fprintf (stderr, "Error, found less than four columns in %s\n", 
	       filename);
      return -1;
    }
    
    if (rows<1) {
      if (!quiet)
	fprintf (stderr, "Error, no header in %s\n", filename);
      return -1;
  }
    
    
    /* 1st line, number of cells */
    *Nx = (int) (value[0][0] + 0.5);
    *Ny = (int) (value[0][1] + 0.5);
    
    
    /* 1st line, horizontal distances, convert from km to m */
    *delX = value[0][2] * 1000.0;
    *delY = value[0][3] * 1000.0;
    
    if (!quiet) {
      fprintf (stderr, " ... x-distance: %g\n", *delX);
      fprintf (stderr, " ... y-distance: %g\n", *delY);
    }
    
    /* allocate memory for albedo */
    *albedo = calloc((size_t) *Nx, sizeof(double *));
    
    for (ix=0; ix<*Nx; ix++)
      (*albedo)[ix] = calloc((size_t) *Ny, sizeof(double));
    
    
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
      
      /* copy data to result arrays */
      (*albedo)[ix][iy] = value[i][2];
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



/*******************************************************/
/* Read data from a 2D umu/phi description file.        */
/*******************************************************/

int read_2D_umu (char     *filename, 
		 int      *Nx,
		 int      *Ny,
		 double   *delX,
		 double   *delY,
		 double ***umu,
		 double ***phi,
		 int       ixmin,
		 int       ixmax,
		 int       iymin,
		 int       iymax,
		 int       quiet)
{
  int i=0, status=0;
  int ix=0, iy=0;

  double **value=NULL;
  int rows=0, min_columns=0, max_columns=0;

#if HAVE_LIBNETCDF
  int ixstart=0, iystart=0;

  size_t tstart[2] = {0,0};
  size_t tcount[2] = {0,0};
 
  size_t n=0;

  int    ncid   =0, idd_nx =0, idd_ny =0;
  int    id_type=0;

  float *temp=NULL;
#endif


  /* try to open as netCDF file */
  
#if HAVE_LIBNETCDF
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status==NC_NOERR) {

    if (!quiet)
      fprintf (stderr, " ... reading 2D umu data from netCDF file %s\n", filename);

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

    *Ny = n;

    /* read "delX" */
    status = nc_get_att_double (ncid, NC_GLOBAL, "dx", delX);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading delX from %s\n", 
	       status, filename);
      return status;
    }

    /* read "delY" */
    status = nc_get_att_double (ncid, NC_GLOBAL, "dy", delY);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading dy from %s\n", 
	       status, filename);
      return status;
    }

    *delX *= 1000.0;
    *delY *= 1000.0;

    if (ixmin!=-1 && ixmax!=-1 && iymin!=-1 && iymax!=-1) {
      if ( ixmax>=*Nx || ixmin<0 || iymax>=*Ny || iymin<0 ) {
	fprintf(stderr,"Error! You have specified a region with `atmos_region` (%d,%d,%d,%d) which is larger than the file %s\n",ixmin,ixmax,iymin,iymax,filename);
	return -1;
      }

      *Nx = ixmax - ixmin + 1;
      *Ny = iymax - iymin + 1;
      ixstart = ixmin;
      iystart = iymin;
    }
    else {
      ixstart=0;
      iystart=0;
    }

    /* allocate memory for 2D umu/phi fields */
    *umu = calloc((size_t) *Nx, sizeof(double *));
    *phi = calloc((size_t) *Nx, sizeof(double *));

    for (ix=0; ix<*Nx; ix++) {
      (*umu)[ix] = calloc((size_t) *Ny, sizeof(double));
      (*phi)[ix] = calloc((size_t) *Ny, sizeof(double));
    }

    /* read data */

    /* get variable id for "umu" */
    status = nc_inq_varid (ncid, "umu", &id_type);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading umu from %s\n", status, filename);
      return status;
    }

    temp = calloc (*Nx, sizeof(float));
    
    /* read type */
    tstart[1] = ixstart;   /* start with first row          */
    tcount[1] = *Nx; /* read *Nx elements for each iy */

    tcount[0] = 1;  /* read one line at a time */

    for (iy=0; iy<*Ny; iy++) {
      tstart[0] = iy + iystart;

      status = nc_get_vara_float (ncid, id_type, tstart, tcount, temp);
      if (status!=NC_NOERR) {
	fprintf (stderr, "Error %d reading umu from %s\n", status, filename);
	return status;
      }

      for (ix=0; ix<*Nx; ix++)
	(*umu)[ix][iy] = temp[ix];
    }

    /* get variable id for "phi" */
    status = nc_inq_varid (ncid, "phi", &id_type);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading phi from %s\n", status, filename);
      return status;
    }

    temp = calloc (*Nx, sizeof(float));
    
    /* read type */
    tstart[1] = ixstart;   /* start with first row          */
    tcount[1] = *Nx; /* read *Nx elements for each iy */

    tcount[0] = 1;  /* read one line at a time */

    for (iy=0; iy<*Ny; iy++) {
      tstart[0] = iy + iystart;

      status = nc_get_vara_float (ncid, id_type, tstart, tcount, temp);
      if (status!=NC_NOERR) {
	fprintf (stderr, "Error %d reading phi from %s\n", status, filename);
	return status;
      }

      for (ix=0; ix<*Nx; ix++)
        /* Switch from DISORT to MYSTIC convention */
	(*phi)[ix][iy] = temp[ix] + 180.0;
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
      fprintf (stderr, "Error, found less than four columns in %s\n", 
	       filename);
      return -1;
    }
    
    if (rows<1) {
      if (!quiet)
	fprintf (stderr, "Error, no header in %s\n", filename);
      return -1;
  }
    
    
    /* 1st line, number of cells */
    *Nx = (int) (value[0][0] + 0.5);
    *Ny = (int) (value[0][1] + 0.5);
    
    
    /* 1st line, horizontal distances, convert from km to m */
    *delX = value[0][2] * 1000.0;
    *delY = value[0][3] * 1000.0;
    
    if (!quiet) {
      fprintf (stderr, " ... x-distance: %g\n", *delX);
      fprintf (stderr, " ... y-distance: %g\n", *delY);
    }
    
    /* allocate memory for umu/phi */
    *umu = calloc((size_t) *Nx, sizeof(double *));
    *phi = calloc((size_t) *Nx, sizeof(double *));
    
    for (ix=0; ix<*Nx; ix++) {
      (*umu)[ix] = calloc((size_t) *Ny, sizeof(double));
      (*phi)[ix] = calloc((size_t) *Ny, sizeof(double));
    }
    
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
      
      /* copy data to result arrays */
      (*umu)[ix][iy] = value[i][2];
      /* Switch from DISORT to MYSTIC convention */
      (*phi)[ix][iy] = value[i][3] + 180.0; 
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


/******************************************************************/
/* Read data from a 2D surface (albedo or RPV) description file.  */
/******************************************************************/

int read_2D_surface_labels (char *filename, 
			    char **label_lib, int nlabel_lib,
			    int *Nx, int *Ny,
			    double *delX, double *delY,
			    unsigned char ***label, int **nlabel,
			    int *therewereCaMs, int quiet)
{
  int i=0, il=0, status=0;
  int ix=0, iy=0;

  int rows=0, min_columns=0, max_columns=0, max_length=0;

  char ***string=NULL;

#if HAVE_LIBNETCDF
  size_t tstart[2] = {0,0};
  size_t tcount[2] = {0,0};
 
  size_t n=0;

  int    ncid=0; 
  int    idd_nx=0, idd_ny=0;
  int    id_type=0;

  unsigned char *temp=NULL;
  char tempstr[4]="";
#endif

  *therewereCaMs=0;

  /* try open as netCDF file */
  
#if HAVE_LIBNETCDF
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status==NC_NOERR) {

    if (!quiet)
      fprintf (stderr, " ... reading surface data from netCDF file %s\n", filename);

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

    /* read "delX" */
    status = nc_get_att_double (ncid, NC_GLOBAL, "dx", delX);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading delX from %s\n", 
	       status, filename);
      return status;
    }

    /* read "delY" */
    status = nc_get_att_double (ncid, NC_GLOBAL, "dy", delY);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading dy from %s\n", 
	       status, filename);
      return status;
    }

    *delX *= 1000.0;
    *delY *= 1000.0;

    /* allocate memory for 2D BRDF field */
    *label  = calloc((size_t) *Nx, sizeof(char *));
    *nlabel = calloc((size_t) nlabel_lib, sizeof(int));

    for (ix=0; ix<*Nx; ix++)
      (*label)[ix] = calloc((size_t) *Ny, sizeof(char));

    /* read data */

    /* get variable id for "type" */
    status = nc_inq_varid (ncid, "type", &id_type);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading type from %s\n", status, filename);
      return status;
    }

    temp = calloc (*Nx, sizeof(unsigned char));

    /* read type */
    tstart[1] = 0;   /* start with first x element    */
    tcount[1] = *Nx; /* read *Nx elements for each iy */

    tcount[0] = 1;   /* read one y element at a time  */

    for (iy=0; iy<*Ny; iy++) {
      tstart[0] = iy;

      if (!quiet)
	fprintf (stderr, ".");
      
      status = nc_get_vara_uchar (ncid, id_type, tstart, tcount, temp);
      if (status!=NC_NOERR) {
	fprintf (stderr, "Error %d reading type from %s\n", status, filename);
	return status;
      }

      for (ix=0; ix<*Nx; ix++) {

	sprintf (tempstr, "%d", temp[ix]);
	
	for (il=0; il<nlabel_lib; il++) {
	  if (!strcmp(tempstr, label_lib[il]))
	    break;
	}

	if (il==nlabel_lib) {
	  fprintf (stderr, "Error, did not find an entry for %s\n", tempstr);
	  return -1;
	}
	
	if (il>UCHAR_MAX) {
	  fprintf (stderr, "Error, index %d larger than %d\n", 
		   il, UCHAR_MAX);
	  return -1;
	}
	  
	(*label) [ix][iy] = il; 
	(*nlabel)[il]++;
      }
    }
    nc_close (ncid);

    free(temp);
    
    if (!quiet)
      fprintf (stderr, "\n");
  }
  else {
    if (!quiet)
      fprintf (stderr, " ... %s not in netCDF format, trying to open as ASCII\n", 
	       filename);
#endif
    
    /* read file; ASCII_file2double is a waste of memory in this */
    /* case, but it is very convenient because comments and      */
    /* everything are handled correctly.                         */
    
    status = ASCII_checkfile (filename, 
			      &rows,
			      &min_columns,
			      &max_columns,
			      &max_length);
    
    if (status!=0) {
      fprintf (stderr, "Error %d reading %s\n", status, filename);
      return status;
    }
    
    if (max_columns<4) {
      fprintf (stderr, "Error, found less than four columns in %s\n", 
	       filename);
      return -1;
    }
    
    if (rows<1) {
      if (!quiet)
      fprintf (stderr, "Error, no header in %s\n", filename);
      return -1;
    }
    
    
    /* allocate memory */
    status = ASCII_calloc_string (&string,
				  rows,
				  max_columns,
				  max_length);
    
    if (status!=0) {
      fprintf (stderr, "Error %d allocating memory\n", status);
      return status;
    }
    
    
    /* read file to string array */
    status = ASCII_readfile (filename, string);
    
    if (status!=0) {
      fprintf (stderr, "Error %d reading %s\n", status, filename);
      return status;
    }
    
  
    /* 1st line, number of cells */
    *Nx = strtol(string[0][0], NULL, 0);
    *Ny = strtol(string[0][1], NULL, 0);
    
    
    /* 1st line, horizontal distances, convert from km to m */
    *delX = strtod(string[0][2], NULL) * 1000.0;
    *delY = strtod(string[0][3], NULL) * 1000.0;
    
    if (!quiet) {
      fprintf (stderr, " ... x-distance: %g\n", *delX);
      fprintf (stderr, " ... y-distance: %g\n", *delY);
    }
    
    /* allocate memory for albedo */
    *label  = calloc((size_t) *Nx, sizeof(char *));
    *nlabel = calloc((size_t) nlabel_lib, sizeof(int));

    for (ix=0; ix<*Nx; ix++)
      (*label)[ix] = calloc((size_t) *Ny, sizeof(char));
    
    
    for (i=1; i<rows; i++) {
      ix = strtol(string[i][0], NULL, 0) - 1;
      iy = strtol(string[i][1], NULL, 0) - 1;
      
      if (ix<0 || ix>=*Nx) {
	fprintf (stderr, "Error, ix = %d out of bounds in line %d\n", ix+1, i+1);
	return -1;
      }
      
      if (iy<0 || iy>=*Ny) {
	fprintf (stderr, "Error, iy = %d out of bounds in line %d\n", iy+1, i+1);
	return -1;
      }
      
      /* copy data to result arrays */
      for (il=0; il<nlabel_lib; il++) {
	if (!strcmp(string[i][2], label_lib[il]))
	  break;
      }
      
      if (il==nlabel_lib) {
	fprintf (stderr, "Error, did not find an entry for %s\n", string[i][2]);
	return -1;
      }
      
      /* if (il>UCHAR_MAX) { */
      /* 	fprintf (stderr, "Error, index %d larger than %d\n",  */
      /* 		 il, UCHAR_MAX); */
      /* 	return -1; */
      /* } */
	  
      (*label) [ix][iy] = il; 
      (*nlabel)[il]++;
    }

    /* free memory */
    (void) ASCII_free_string (string, rows, max_columns);

    if (!quiet)
      fprintf (stderr, " ... read %d data points from %s\n", 
	       rows-1, filename);


    if ((*nlabel)[0]>0)
      *therewereCaMs=1;

#if HAVE_LIBNETCDF
  }    
#endif

  return 0;
}




/*******************************************************/
/* Read data from a 2D Ross Li BRDF description file.  */
/*******************************************************/

int read_2D_rossli (char *filename, 
		    int *Nx, int *Ny,
		    double *delX, double *delY,
		    rossli_brdf_spec ***rossli,
		    int *therewereCaMs,
		    int hotspot, int isAmbralsFile,
		    int ixmin, int ixmax, int iymin, int iymax,
		    int quiet)
{
  int i=0, status=0;
  int ix=0, iy=0;
  int contains_CaM=0, isCaM=0;

#if HAVE_LIBNETCDF
  size_t tstart[2] = {0,0};
  size_t tcount[2] = {0,0};
  int ixstart=0, iystart=0;

  size_t n=0;

  int    ncid   =0, idd_nx =0, idd_ny =0;
  int    id_iscam=0;
  int    id_iso=0, id_geo=0, id_vol=0;

  int *isCaMv=NULL;
  float *temp=NULL;
#endif

  double **value=NULL;
  int rows=0, min_columns=0, max_columns=0;

  double vol_fac=1.0, iso_fac=1.0, geo_fac=1.0;

  *therewereCaMs=0;

  if (isAmbralsFile){
    vol_fac=3./4.;
    iso_fac=1./PI;
    geo_fac=1./PI;
  }

  /* try open as netCDF file */
  
#if HAVE_LIBNETCDF
  status = nc_open (filename, NC_NOWRITE, &ncid);
  if (status==NC_NOERR) {

    if (!quiet)
      fprintf (stderr, " ... reading Ross-Li data from netCDF file %s\n", filename);

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

    /* read "delX" */
    status = nc_get_att_double (ncid, NC_GLOBAL, "dx", delX);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading delX from %s\n", 
	       status, filename);
      return status;
    }

    /* read "delY" */
    status = nc_get_att_double (ncid, NC_GLOBAL, "dy", delY);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading dy from %s\n", 
	       status, filename);
      return status;
    }

    *delX *= 1000.0;
    *delY *= 1000.0;

    if (ixmin!=-1 && ixmax!=-1 && iymin!=-1 && iymax!=-1) {
      if ( ixmax>=*Nx || ixmin<0 || iymax>=*Ny || iymin<0 ) {
	fprintf(stderr,"Error! You have specified a region with `atmos_region` (%d,%d,%d,%d) which is larger than the file %s\n",ixmin,ixmax,iymin,iymax,filename);
	return -1;
      }

      *Nx = ixmax - ixmin + 1;
      *Ny = iymax - iymin + 1;
      ixstart = ixmin;
      iystart = iymin;
    }
    else {
      ixstart=0;
      iystart=0;
    }

    /* allocate memory for albedo */
    *rossli = calloc((size_t) *Nx, sizeof(rossli_brdf_spec *));

    for (ix=0; ix<*Nx; ix++)
      (*rossli)[ix] = calloc((size_t) *Ny, sizeof(rossli_brdf_spec));

    /* read data */

    /* get variable id for "iso" */
    status = nc_inq_varid (ncid, "iso", &id_iso);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading iso from %s\n", status, filename);
      return status;
    }

    /* get variable id for "geo" */
    status = nc_inq_varid (ncid, "geo", &id_geo);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading geo from %s\n", status, filename);
      return status;
    }

    /* get variable id for "vol" */
    status = nc_inq_varid (ncid, "vol", &id_vol);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading vol from %s\n", status, filename);
      return status;
    }

    /* get variable id for "iscam" */
    status = nc_inq_varid (ncid, "isCaM", &id_iscam);
    if (status!=NC_NOERR) {
      fprintf (stderr, "Error %d reading iscam from %s\n", status, filename);
      return status;
    }

    temp = calloc (*Nx, sizeof(float));
    isCaMv = calloc (*Nx, sizeof(int));

    /* read type */
    tstart[1] = ixstart;   /* start with first x element    */
    tcount[1] = *Nx; /* read *Nx elements for each iy */

    tcount[0] = 1;   /* read one y element at a time  */

    for (iy=0; iy<*Ny; iy++) {
      tstart[0] = iy + iystart;

      if (!quiet)
	fprintf (stderr, ".");
      
      status = nc_get_vara_int (ncid, id_iscam, tstart, tcount, isCaMv);
      if (status!=NC_NOERR) {
	fprintf (stderr, "Error %d reading isCaM from %s\n", status, filename);
	return status;
      }

      for (ix=0; ix<*Nx; ix++)
	if (isCaMv[ix]==1) {
	  (*rossli)[ix][iy].isCaM=1;
	  *therewereCaMs=1;
	}
	else
	  (*rossli)[ix][iy].isCaM=0;

      status = nc_get_vara_float (ncid, id_iso, tstart, tcount, temp);
      if (status!=NC_NOERR) {
	fprintf (stderr, "Error %d reading iso from %s\n", status, filename);
	return status;
      }

      for (ix=0; ix<*Nx; ix++)
	if (isCaMv[ix]==1)
	  (*rossli)[ix][iy].iso=0.0;
	else
	  (*rossli)[ix][iy].iso=iso_fac*temp[ix];

      status = nc_get_vara_float (ncid, id_geo, tstart, tcount, temp);
      if (status!=NC_NOERR) {
	fprintf (stderr, "Error %d reading geo from %s\n", status, filename);
	return status;
      }

      for (ix=0; ix<*Nx; ix++)
	if (isCaMv[ix]==1)
	  (*rossli)[ix][iy].geo=0.0;
	else
	  (*rossli)[ix][iy].geo=geo_fac*temp[ix];

      status = nc_get_vara_float (ncid, id_vol, tstart, tcount, temp);
      if (status!=NC_NOERR) {
	fprintf (stderr, "Error %d reading vol from %s\n", status, filename);
	return status;
      }

      for (ix=0; ix<*Nx; ix++)
	if (isCaMv[ix]==1)
	  (*rossli)[ix][iy].vol=0.0;
	else
	  (*rossli)[ix][iy].vol=vol_fac*temp[ix];

      for (ix=0; ix<*Nx; ix++)
	if (isCaMv[ix]==1)
	  (*rossli)[ix][iy].hotspot=0;
	else
	  (*rossli)[ix][iy].hotspot=hotspot;

    }
    nc_close (ncid);

    free(temp);
    free(isCaMv);

    if (!quiet)
      fprintf (stderr, "\n");
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

    if (max_columns<5) {
      fprintf (stderr, "Error, found less than five columns in %s\n", 
	       filename);
      return -1;
    }

    /* if more than 5 columns, the 6th column defines whether CaM or not */
    if (max_columns > 5)
      contains_CaM=1;
  
    if (rows<1) {
      if (!quiet)
	fprintf (stderr, "Error, no header in %s\n", filename);
      return -1;
    }


    /* 1st line, number of cells */
    *Nx = (int) (value[0][0] + 0.5);
    *Ny = (int) (value[0][1] + 0.5);


    /* 1st line, horizontal distances, convert from km to m */
    *delX = value[0][2] * 1000.0;
    *delY = value[0][3] * 1000.0;

    if (!quiet) {
      fprintf (stderr, " ... x-distance: %g\n", *delX);
      fprintf (stderr, " ... y-distance: %g\n", *delY);
    }

    /* allocate memory for albedo */
    *rossli = calloc((size_t) *Nx, sizeof(rossli_brdf_spec *));

    for (ix=0; ix<*Nx; ix++)
      (*rossli)[ix] = calloc((size_t) *Ny, sizeof(rossli_brdf_spec));


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

      if (contains_CaM)
	isCaM = (value[i][5] == 1);

      if (isCaM) {
	/* CaM */
	(*rossli)[ix][iy].iso = 0.0;
	(*rossli)[ix][iy].vol = 0.0;
	(*rossli)[ix][iy].geo = 0.0;
	(*rossli)[ix][iy].hotspot = 0;
	(*rossli)[ix][iy].isCaM = 1;
	*therewereCaMs=1;
      }
      else {
	/* copy data to result arrays */
	(*rossli)[ix][iy].iso = iso_fac*value[i][2];
	(*rossli)[ix][iy].vol = geo_fac*value[i][3];
	(*rossli)[ix][iy].geo = vol_fac*value[i][4];
	(*rossli)[ix][iy].hotspot = hotspot;
	(*rossli)[ix][iy].isCaM = 0;
      }
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
