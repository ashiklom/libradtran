/*--------------------------------------------------------------------
 * $Id: netCDF_functions.c 3183 2015-09-09 10:11:29Z Claudia.Emde $
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
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#if HAVE_LIBNETCDF
#include <netcdf.h>
#endif

#include "netCDF_functions.h"
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}


/***********************************************************************************/
/* Function: write_netCDF_1D_float                                                 */
/* Description:                                                                    */
/*  Write a 1D float (like latitude array), according attributes and dimension     */
/*  (number of array entries) to an opened netCDF result file.                     */
/*  The function checks, if this variables already                                 */
/*  exists in the netCDF header. If not, it is created.                            */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the netCDF input file                       */
/*     char *filename     filename of the result file                              */
/*     char *data_name    short name of the data set, e.g. "lat"                   */
/*     int *id_data       id number of the data (is created, of data is new)       */
/*     char *long_name    long name of the data set. e.g."latitude"                */
/*     char *units        units of the data set, e.g. "degrees east"               */
/*     float *data        the data set itself                                      */
/*     char *dim_name     name of the according dimension                          */
/*                        (in general the same as data_name), e.g. "lat"           */
/*     int id_dim         id number of the dimension                               */
/*     int n              number of array entries (dimension size)                 */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int write_netCDF_1D_float (int ncid, char *filename, 
                           char *data_name, int *id_data, char *long_name, char *units, float *data, 
                           char *dim_name,  int *id_dim,  int n)
{

#if HAVE_LIBNETCDF

  int status=NC_NOERR;
  size_t iv=0;

  char function_name[]="write_netCDF_1D_float";
  char file_name[]="netCDF_functions.c";

  /* get variable id for data_name */

  status = nc_inq_varid (ncid, data_name, id_data);
  if (status!=NC_NOERR) {

    /* variable does not exist, create it! */
    status = nc_redef(ncid); /* put in define mode */
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d putting netCDF file %s into define mode, in %s (%s)\n", status, filename, function_name, file_name);
      return status;
    }

    /* create dimension */
    status = nc_def_dim(ncid, dim_name, n, id_dim);
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d defining data '%s' in netCDF file %s, in %s (%s)\n", status, data_name, filename, function_name, file_name);
      return status;
    }

    status = nc_def_var (ncid, data_name, NC_FLOAT, 1, id_dim, id_data);
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d defining data %s in netCDF file %s, in %s (%s)\n", status, data_name, filename, function_name, file_name);
      return status;
    }

    status = nc_enddef(ncid); /* end define mode */
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d leaving the define mode for %s, in %s (%s) \n", status, filename, function_name, file_name);
      return status;    
    }
        
    /* write data into the variable */
    for ( iv=0; iv< n; iv++ ) {
      /* WRITE ONE VARIABLE */
      status = nc_put_var1_float(ncid, (*id_data), &iv, &(data[iv]));
      if (status != NC_NOERR) {
        fprintf (stderr, "Error %d writing %s in netCDF file %s, in %s (%s)\n", status, data_name, filename, function_name, file_name);
        return status;
      }
    }
  }
  else {
    /* yes, check if they are identical */
  }

  return status;

#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif

}

/***********************************************************************************/
/* Function: write_netCDF_3D_float                                                 */
/* Description:                                                                    */
/*  Write one result (float) (time, lat, lon) into a 3D array                      */
/*  The function checks, if this data set already exits  in  the netCDF header.    */
/*  If not it will be created with according attributes, units ... and             */
/*  initialiesed to 0.                                                             */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the netCDF input file                       */
/*     char *filename     filename of the result file                              */
/*     char *data_name    short name of the data set, e.g. "lat"                   */
/*     char *long_name    long name of the data set. e.g."latitude"                */
/*     char *units        units of the data set, e.g. "degrees east"               */
/*     float data         data to write                                            */
/*     int *id_dim_array  array of dimension id numbers                            */
/*     size_t *index_3D   where to write the data                                  */
/*     int ndim           number of dimenision (here always 3)                     */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Jan 2008   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int write_netCDF_3D_float (int ncid, char *filename,
                           char *data_name, char *long_name, char *units, float data, 
                           int *id_dim_array, size_t *index, int ndim)
{

#if HAVE_LIBNETCDF

  int status  = NC_NOERR;
  int id_data = NOT_DEFINED_INTEGER;

  char function_name[]="write_netCDF_3D_float";
  char file_name[]="netCDF_functions.c";

  status = create_float_variable_or_get_id (ncid, filename, 
                                            data_name, &id_data, ndim, id_dim_array, 
                                            long_name, units);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error %d creating float variable %s in netCDF file %s, in %s (%s)\n", status, long_name, filename, function_name, file_name);
    return status;
  }

  /* WRITE ONE VARIABLE */
  status = nc_put_var1_float(ncid, id_data, index, &(data));
  if (status != NC_NOERR) {
    fprintf (stderr, "Error %d writing %s in netCDF file %s, in %s (%s)\n", status, data_name, filename, function_name, file_name);
    return status;
  }

  return status;

#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif

}

/***********************************************************************************/
/* Function: write_netCDF_4D_float                                                 */
/* Description:                                                                    */
/*  Write one column of results (mostly z-column) into a 4D array                  */
/*  (time, z, lat, lon). The function checks, if this data set already exits in    */
/*  the netCDF header. If not it will be created with according attributes,        */
/*  units ...                                                                      */ 
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the netCDF input file                       */
/*     char *filename     filename of the result file                              */
/*     char *data_name    short name of the data set, e.g. "lat"                   */
/*     int *id_data       id number of the data (is created, of data is new)       */
/*     char *long_name    long name of the data set. e.g."latitude"                */
/*     char *units        units of the data set, e.g. "degrees east"               */
/*     float *data        the data set itself                                      */
/*     int *id_dim_array  array of dimension id numbers                            */
/*     size_t *index      specifies the position, where to write the entry         */
/*                        into the 4D file                                         */
/*     int ndim           number of dimenision (here always 4)                     */
/*     size_t *n_end      size of the dimensions                                   */
/*     int *zout_index    index of the zout levels                                 */
/*                        (nessesary, as some zout-levels, which are below the     */ 
/*                        ground altitude are obmitted)                            */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int write_netCDF_4D_float (int ncid, char *filename, 
                           char *data_name, char *long_name, char *units, float *data, 
                           int *id_dim_array, size_t *index, int ndim, size_t *n_end, int *zout_index, int *simulation_done)
{

#if HAVE_LIBNETCDF

  int status  = NC_NOERR;
  int id_data = NOT_DEFINED_INTEGER;
  size_t lev  = 0;

  char function_name[]="write_netCDF_4D_float";
  char file_name[]="netCDF_functions.c";

  /* check if this result is valid */
  for (lev=0; lev<n_end[POSITION_ZOUT_4D]; lev++) {
    /* replace infinite numbers with NaN (Not a number) */
    if ( isinf(data[lev]) )
      data[lev]=0.0/0.0;
    /* check for NaN (Not a Number) */
    if ( isnan(data[lev]) ) {
      (*simulation_done)=FALSE;
      fprintf (stderr, " *** Warning, unresonable result for '%s' (lev=%d)\n", long_name, (int) lev);
    }
  }

  status = create_float_variable_or_get_id (ncid, filename, 
                                            data_name, &id_data, ndim, id_dim_array, 
                                            long_name, units);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error %d creating float variable %s in netCDF file %s, in %s (%s)\n", status, long_name, filename, function_name, file_name);
    return status;
  }

  /* write data into the variable */
  for ( lev = 0; lev < n_end[POSITION_ZOUT_4D]; lev++ ) {
    index[POSITION_ZOUT_4D] = zout_index[lev];
  
    /* WRITE ONE VARIABLE */
    status = nc_put_var1_float(ncid, id_data, index, &(data[lev]));
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d writing %s in netCDF file %s, in %s (%s)\n", status, data_name, filename, function_name, file_name);
      return status;
    }
  }

  return status;

#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif

}

/***********************************************************************************/
/* Function: write_netCDF_5D_float                                                 */
/* Description:                                                                    */
/*  Write one column of results (mostly z-column) into a 5D array                  */
/*  (time, z, sza, lat, lon). The function checks, if this data set already exits  */
/*   in  the netCDF header. If not it will be created with according attributes,   */
/*  units ...                                                                      */
/*  If result is not valid (infinity or Not-a-number) simulation done will be      */ 
/*  set to FALSE == not done. infinity is replaced by Not-a-number.                */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the netCDF input file                       */
/*     char *filename     filename of the result file                              */
/*     char *data_name    short name of the data set, e.g. "lat"                   */
/*     int *id_data       id number of the data (is created, of data is new)       */
/*     char *long_name    long name of the data set. e.g."latitude"                */
/*     char *units        units of the data set, e.g. "degrees east"               */
/*     float *data        the data set itself                                      */
/*     int *id_dim_array  array of dimension id numbers                            */
/*     size_t *index      specifies the position, where to write the entry         */
/*                        into the 4D file                                         */
/*     int ndim           number of dimenision (here always 5)                     */
/*     size_t *n_end      size of the dimensions                                   */
/*     int *zout_index    index of the zout levels                                 */
/*                        (nessesary, as some zout-levels, which are below the     */ 
/*                        ground altitude are obmitted)                            */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/


int write_netCDF_5D_float (int ncid, char *filename, 
                           char *data_name, char *long_name, char *units, float **data, 
                           int *id_dim_array, size_t *index, int ndim, size_t *n_end, int *zout_index, int *simulation_done)

{

#if HAVE_LIBNETCDF

  int status  = NC_NOERR;
  int id_data = NOT_DEFINED_INTEGER;
  size_t lev  = 0;
  size_t iv  = 0;

  char function_name[]="write_netCDF_5D_float";
  char file_name[]="netCDF_functions.c";

  /* check if this result is valid */
  for (lev=0; lev<n_end[POSITION_ZOUT_5D]; lev++)
    for (iv=0; iv<n_end[POSITION_LAMBDA_5D]; iv++) {
      /* replace infinite numbers with NaN (Not a number) */
      if ( isinf(data[lev][iv]) )
        data[lev][iv]=0.0/0.0;
      /* check for NaN (Not a Number) */
      if ( isnan(data[lev][iv]) ) {
        (*simulation_done)=FALSE;
        fprintf (stderr, " *** Warning, unresonable result for '%s' (iv=%d, lev=%d)\n", long_name, (int) iv, (int) lev);
      }
    }


  status = create_float_variable_or_get_id (ncid, filename, 
                                            data_name, &id_data, ndim, id_dim_array, 
                                            long_name, units);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error %d creating float variable %s in netCDF file %s, in %s (%s)\n", status, long_name, filename, function_name, file_name);
    return status;
  }


  for ( lev = 0; lev < n_end[POSITION_ZOUT_5D]; lev++ ) {
    index[POSITION_ZOUT_5D] = zout_index[lev];

    for ( iv = 0; iv < n_end[POSITION_LAMBDA_5D]; iv++ ) {
      index[POSITION_LAMBDA_5D] = iv;

      /* WRITE ONE VARIABLE */
      status = nc_put_var1_float(ncid, id_data, index, &(data[lev][iv]));
      if (status != NC_NOERR) {
        fprintf (stderr, "Error %d writing %s in netCDF file %s, in %s (%s)\n", status, data_name, filename, function_name, file_name);
        return status;
      }
    }
  }

  return status;

#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif

}

/***********************************************************************************/
/* Function: write_netCDF_6D_float                                                 */
/* Description:                                                                    */
/*  Write one column of results (mostly z-column) into a 6D array                  */
/*  (time, z, sza, lat, lon). The function checks, if this data set already exits  */
/*   in  the netCDF header. If not it will be created with according attributes,   */
/*  units ...                                                                      */
/*  If result is not valid (infinity or Not-a-number) simulation done will be      */ 
/*  set to FALSE == not done. infinity is replaced by Not-a-number.                */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the netCDF input file                       */
/*     char *filename     filename of the result file                              */
/*     char *data_name    short name of the data set, e.g. "lat"                   */
/*     int *id_data       id number of the data (is created, of data is new)       */
/*     char *long_name    long name of the data set. e.g."latitude"                */
/*     char *units        units of the data set, e.g. "degrees east"               */
/*     float *data        the data set itself                                      */
/*     int *id_dim_array  array of dimension id numbers                            */
/*     size_t *index      specifies the position, where to write the entry         */
/*                        into the 4D file                                         */
/*     int ndim           number of dimenision (here always 6)                     */
/*     size_t *n_end      size of the dimensions                                   */
/*     int *zout_index    index of the zout levels                                 */
/*                        (nessesary, as some zout-levels, which are below the     */ 
/*                        ground altitude are obmitted)                            */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int write_netCDF_6D_float (int ncid, char *filename, 
                        char *data_name, char *long_name, char *units, float ***data, 
                        int *id_dim_array, size_t *index, int ndim, size_t *n_end, int *zout_index, int *simulation_done)
{

#if HAVE_LIBNETCDF

  int status  = NC_NOERR;
  int id_data = NOT_DEFINED_INTEGER;
  size_t lev  = 0;
  size_t iv  = 0;
  size_t iu  = 0;

  char function_name[]="write_netCDF_6D_float";
  char file_name[]="netCDF_functions.c";

  /* check if this result is valid */
  for (lev=0; lev<n_end[POSITION_ZOUT_6D]; lev++)  
    for (iu=0; iu<n_end[POSITION_UMU_6D]; iu++)
      for (iv=0; iv<n_end[POSITION_LAMBDA_6D]; iv++) {
        /* replace infinite numbers with NaN (Not a number) */
        if ( isinf(data[lev][iu][iv]) )
          data[lev][iu][iv]=0.0/0.0;
        /* check for NaN (Not a Number) */
        if ( isnan(data[lev][iu][iv]) ) {
          (*simulation_done)=FALSE;
          fprintf (stderr, " *** Warning, unresonable result for '%s' (iv=%d, lev=%d, iu=%d)\n", long_name, (int) iv, (int) lev, (int) iu);
        }
      }

  status = create_float_variable_or_get_id (ncid, filename, 
                                            data_name, &id_data, ndim, id_dim_array, 
                                            long_name, units);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error %d creating float variable %s in netCDF file %s, in %s (%s)\n", status, long_name, filename, function_name, file_name);
    return status;
  }

  for ( lev = 0; lev < n_end[POSITION_ZOUT_6D]; lev++ ) {
    index[POSITION_ZOUT_6D] = zout_index[lev];

    for ( iv = 0; iv < n_end[POSITION_LAMBDA_6D]; iv++ ) {
      index[POSITION_LAMBDA_6D] = iv;

      for ( iu = 0; iu < n_end[POSITION_UMU_6D]; iu++ ) {
        index[POSITION_UMU_6D] = iu;

        /* WRITE ONE VARIABLE */
        status = nc_put_var1_float(ncid, id_data, index, &(data[lev][iu][iv]));
        if (status != NC_NOERR) {
          fprintf (stderr, "Error %d writing %s in netCDF file %s, in %s (%s)\n", status, data_name, filename, function_name, file_name);
          return status;
        }
      }
    }
  }

  return status;

#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif

}

/***********************************************************************************/
/* Function: write_netCDF_7D_float                                                 */
/* Description:                                                                    */
/*  Write one column of results (mostly z-column) into a 7D array                  */
/*  (time, z, sza, lat, lon). The function checks, if this data set already exits  */
/*   in  the netCDF header. If not it will be created with according attributes,   */
/*  units ...                                                                      */
/*  If result is not valid (infinity or Not-a-number) simulation done will be      */ 
/*  set to FALSE == not done. infinity is replaced by Not-a-number.                */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the netCDF input file                       */
/*     char *filename     filename of the result file                              */
/*     char *data_name    short name of the data set, e.g. "lat"                   */
/*     int *id_data       id number of the data (is created, of data is new)       */
/*     char *long_name    long name of the data set. e.g."latitude"                */
/*     char *units        units of the data set, e.g. "degrees east"               */
/*     float *data        the data set itself                                      */
/*     int *id_dim_array  array of dimension id numbers                            */
/*     size_t *index      specifies the position, where to write the entry         */
/*                        into the 4D file                                         */
/*     int ndim           number of dimenision (here always 7)                     */
/*     size_t *n_end      size of the dimensions                                   */
/*     int *zout_index    index of the zout levels                                 */
/*                        (nessesary, as some zout-levels, which are below the     */ 
/*                        ground altitude are obmitted)                            */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int write_netCDF_7D_float (int ncid, char *filename, 
                        char *data_name, char *long_name, char *units, float ****data, 
                        int *id_dim_array, size_t *index, int ndim, size_t *n_end, int *zout_index, int *simulation_done)
{

#if HAVE_LIBNETCDF

  int status  = NC_NOERR;
  int id_data = NOT_DEFINED_INTEGER;
  size_t lev=0, iv=0, iu=0, ip=0 ;

  char function_name[]="write_netCDF_7D_float";
  char file_name[]="netCDF_functions.c";

  /* check if this result is valid */
  for (lev=0; lev<n_end[POSITION_ZOUT_6D]; lev++)  
    for ( ip = 0; ip < n_end[POSITION_PHI_7D]; ip++ )
      for (iu=0; iu<n_end[POSITION_UMU_6D]; iu++)
        for (iv=0; iv<n_end[POSITION_LAMBDA_6D]; iv++) {
          /* replace infinite numbers with NaN (Not a number) */
          if ( isinf(data[lev][ip][iu][iv]) )
            data[lev][ip][iu][iv]=0.0/0.0;
          /* check for NaN (Not a Number) */
          if ( isnan(data[lev][ip][iu][iv]) ) {
            (*simulation_done)=FALSE;
            fprintf (stderr, " *** Warning, unresonable result for '%s' (iv=%d, lev=%d, iu=%d, ip=%d)\n", long_name, (int)iv, (int)lev, (int)iu, (int)ip);
          }
        }

  status = create_float_variable_or_get_id (ncid, filename, 
                                            data_name, &id_data, ndim, id_dim_array, 
                                            long_name, units);
  if (status != NC_NOERR) {
    fprintf (stderr, "Error %d creating float variable %s in netCDF file %s, in %s (%s)\n", status, long_name, filename, function_name, file_name);
    return status;
  }

  for ( lev = 0; lev < n_end[POSITION_ZOUT_7D]; lev++ ) {
    index[POSITION_ZOUT_7D] = zout_index[lev];

    for ( iv = 0; iv < n_end[POSITION_LAMBDA_7D]; iv++ ) {
      index[POSITION_LAMBDA_7D] = iv;

      for ( iu = 0; iu < n_end[POSITION_UMU_7D]; iu++ ) {
        index[POSITION_UMU_7D] = iu;

        for ( ip = 0; ip < n_end[POSITION_PHI_7D]; ip++ ) {
          index[POSITION_PHI_7D] = ip;

          /* WRITE ONE VARIABLE */
          status = nc_put_var1_float(ncid, id_data, index, &(data[lev][ip][iu][iv]));
          if (status != NC_NOERR) {
            fprintf (stderr, "Error %d writing %s in netCDF file %s, in %s (%s)\n", status, data_name, filename, function_name, file_name);
            return status;
          }
        }
      }
    }
  }
  return status;

#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif

}

/***********************************************************************************/
/* Function: create_float_variable_or_get_id                                       */
/* Description:                                                                    */
/*  The function checks, if this data set already exits                            */
/*  in  the netCDF header. If not it will be created with according attributes,    */
/*  units ...                                                                      */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the netCDF input file                       */
/*     char *filename     filename of the result file                              */
/*     char *data_name    short name of the data set, e.g. "lat"                   */
/*     int *id_data       id number of the data (is created, of data is new)       */
/*     char *long_name    long name of the data set. e.g."latitude"                */
/*     char *units        units of the data set, e.g. "degrees east"               */
/*     int *id_dim_array  array of dimension id numbers                            */
/*     int ndim           number of dimenision                                     */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int create_float_variable_or_get_id ( int ncid, char *filename, 
                                      char *data_name, int *id_data, int ndim, int *id_dim_array, 
                                      char *long_name, char *units)
{

#if HAVE_LIBNETCDF

  int status = NC_NOERR;
  float missing_value = -1.E10;

  char function_name[]="create_float_variable_or_get_id";
  char file_name[]="netCDF_functions.c";

  /* get variable id for data_name */
  status = nc_inq_varid (ncid, data_name, id_data);

  if (status != NC_NOERR) {

    /* variable does not exist, create it! */
    status = nc_redef(ncid); /* put in define mode */
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d putting netCDF file %s into define mode, in %s (%s)\n", status, filename, function_name, file_name);
      return status;
    }

    /* create variable */
    status = nc_def_var (ncid, data_name, NC_FLOAT, ndim, id_dim_array, id_data);
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d defining data '%s' in netCDF file %s, in %s (%s)\n", status, data_name, filename, function_name, file_name);
      return status;
    }

    status = nc_put_att_text (ncid, *id_data, "long_name", strlen(long_name), long_name);
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d put attribute long_name to '%s' in netCDF file %s, in %s (%s)\n", status, data_name, filename, function_name, file_name);
      return status;
    }

    status = nc_put_att_text (ncid, *id_data, "units", strlen(units), units);
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d put attribute units to '%s' in netCDF file %s, in %s (%s)\n", status, data_name, filename, function_name, file_name);
      return status;
    }

    status = nc_put_att_float (ncid, *id_data, "missing_value", NC_FLOAT, 1, &(missing_value));
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d defining missing argument in netCDF file %s, in %s (%s)\n", status, filename, function_name, file_name);
      return status;
    }

    status = nc_put_att_float (ncid, *id_data, "_FillValue", NC_FLOAT, 1, &(missing_value));
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d defining missing argument in netCDF file %s, in %s (%s)\n", status, filename, function_name, file_name);
      return status;
    }

    status = nc_enddef(ncid); /* end define mode */
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d leaving the define mode for %s, in %s (%s) \n", status, filename, function_name, file_name);
      return status;    
    }
  }

  return NC_NOERR;

#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif

}

/***********************************************************************************/
/* Function: write_netCDF_3D_byte                                                  */
/* Description:                                                                    */
/*  Write one number of results (TRUE or FLASE) into a 3D array (simulation done)  */
/*  The function checks, if this data set already exits  in  the netCDF header.    */
/*  If not it will be created with according attributes, units ... and             */
/*  initialiesed to 0.                                                             */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the netCDF input file                       */
/*     char *filename     filename of the result file                              */
/*     char *data_name    short name of the data set, e.g. "lat"                   */
/*     char *long_name    long name of the data set. e.g."latitude"                */
/*     char *units        units of the data set, e.g. "degrees east"               */
/*     short data         data to write                                            */
/*     int *id_dim_array  array of dimension id numbers                            */
/*     size_t *index_3D   where to write the data                                  */
/*     int ndim           number of dimenision (here always 3)                     */
/* Return value:                                                                   */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int write_netCDF_3D_byte (int ncid, char *filename,
                          char *data_name, char *long_name, char *units, short data, 
                          int *id_dim_array, size_t *index_3D, int ndim)
{

#if HAVE_LIBNETCDF

  int status  = NC_NOERR;
  int id_data = NOT_DEFINED_INTEGER;

  short missing_value = 0;
  size_t ntime=0, nlat=0, nlon=0;
  int n_all = 0;
  int t=-999;
  signed char *tmp_char = NULL;

  char function_name[]="write_netCDF_3D_byte";
  char file_name[]="netCDF_functions.c";

  /* get variable id for data_name */
  status = nc_inq_varid (ncid, data_name, &id_data);

  if (status != NC_NOERR) {

    /* variable does not exist, create it! */
    status = nc_redef(ncid); /* put in define mode */
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d putting netCDF file %s into define mode, in %s (%s)\n", status, filename, function_name, file_name);
      return status;
    }

    /* create variable */
    status = nc_def_var (ncid, data_name, NC_BYTE, ndim, id_dim_array, &id_data);
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d defining data '%s' in netCDF file %s, in %s (%s)\n", status, data_name, filename, function_name, file_name);
      return status;
    }

    status = nc_put_att_text (ncid, id_data, "long_name", strlen(long_name), long_name);
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d put attribute long_name to '%s' in netCDF file %s, in %s (%s)\n", 
                        status, data_name, filename, function_name, file_name);
      return status;
    }

    status = nc_put_att_text (ncid, id_data, "units", strlen(units), units);
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d put attribute units to '%s' in netCDF file %s, in %s (%s)\n", 
                       status, data_name, filename, function_name, file_name);
      return status;
    }

    status = nc_put_att_short (ncid, id_data, "missing_value", NC_BYTE, 1, &(missing_value));
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d defining missing argument in netCDF file %s, in %s (%s)\n", status, filename, function_name, file_name);
      return status;
    }

    status = nc_enddef(ncid); /* end define mode */
    if (status != NC_NOERR) {
      fprintf (stderr, "Error %d leaving the define mode for %s, in %s (%s) \n", status, filename, function_name, file_name);
      return status;    
    }

    /* get dimension length for "ntime" */
    status = nc_inq_dimlen (ncid, id_dim_array[POSITION_TIME_3D], &ntime); 
    if (status!=NC_NOERR) { 
      fprintf (stderr, "Error %d reading ntime from %s, in %s (%s)\n", status, filename, function_name, file_name); 
      return status; 
    } 

    /* get dimension length for "nlat" */
    status = nc_inq_dimlen (ncid, id_dim_array[POSITION_LAT_3D], &nlat); 
    if (status!=NC_NOERR) { 
      fprintf (stderr, "Error %d reading nlat from %s, in %s (%s)\n", status, filename, function_name, file_name); 
      return status; 
    } 

    /* get dimension length for "nlon" */
    status = nc_inq_dimlen (ncid, id_dim_array[POSITION_LON_3D], &nlon); 
    if (status!=NC_NOERR) { 
      fprintf (stderr, "Error %d reading nlon from %s, in %s (%s)\n", status, filename, function_name, file_name); 
      return status; 
    } 

    n_all = ntime * nlat * nlon;

    if ((tmp_char  = calloc (n_all, sizeof (char))) == NULL) {
      fprintf (stderr,"Error, allocating memory for %s, in %s (%s)\n", data_name, function_name, file_name);
      return -1;
    }

    for (t=0; t<n_all; t++)
      tmp_char[t]=0;

    /* WRITE ONE VARIABLE */
    /* does not work on the cluster with unsigned char, that why I use signed chars here !?! */
    status = nc_put_var_schar(ncid, id_data, tmp_char); 
    if (status != NC_NOERR) { 
      fprintf (stderr, "Error %d writing %s in netCDF file %s, in %s (%s)\n", status, data_name, filename, function_name, file_name); 
      return status; 
    } 

    free(tmp_char);
  }

  /* WRITE ONE VARIABLE */
  status = nc_put_var1_short(ncid, id_data, index_3D, &(data));
  if (status != NC_NOERR) {
    fprintf (stderr, "Error %d writing %s to netCDF file %s, in %s (%s)\n", status, data_name, filename, function_name, file_name);
    return status;
  }

  return status;

#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif

}


/***********************************************************************************/
/* Function: read_3d_float                                                         */
/* Description:                                                                    */
/*  Reads a 3D data set from an netCDF file and stores it in an 3D float array,    */
/*  which will be automatically be allocated                                       */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the read netCDF file                        */
/*     char *file         filename of the read netCDF file                         */
/*     char *var_name     short name of the data set, e.g. "lat"                   */
/*     long n1            number of entries in the 1st dimension                   */
/*     long n2            number of entries in the 2nd dimension                   */
/*     long n3            number of entries in the 3rd dimension                   */
/* Return value:                                                                   */
/*     float ****float3d  read data                                                */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Jun 2008   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int read_3d_float ( int ncid, char *file, char *var_name, float ****float3d, 
                    long n1, long n2, long n3, char ***done, int verbose )
{

#if HAVE_LIBNETCDF

  int status=0;
  int t=0, j=0, i=0;
  int id_float3d=0;
  float *tmp_float3d = NULL;

  char function_name[]="read_3d_float";
  char file_name[]="netCDF_functions.c";

  if (verbose)
    fprintf(stderr," ... read_float_3d: var = %-10s, data_type = %-8s, ndim = %3d \n", 
                         var_name, "float", 3);


  /* allocate space for "float3d" */
  if (((*float3d)  = calloc (n1, sizeof (float ***)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'float3d' in %s (%s)\n", function_name, file_name);
    return -1;
  }

  for (t=0;t<n1;t++) {
    if (((*float3d)[t]  = calloc (n2, sizeof (float **)))   == NULL) {
      fprintf (stderr,"Error, allocating memory for 'float3d' in %s (%s)\n", function_name, file_name);
      return -2;
    }

    for (j=0;j<n2;j++) {
      if (((*float3d)[t][j] = calloc (n3, sizeof (float)))   == NULL) {
        fprintf (stderr,"Error, allocating memory for 'float3d' in %s (%s)\n", function_name, file_name);
        return -3;
      }
    }
  }

  /* allocate space for "tmp_float3d" */
  if ( ( tmp_float3d  = calloc (n1*n2*n3, sizeof (float)) ) == NULL ) {
    fprintf (stderr,"Error, allocating memory for 'tmp_float3d' in %s (%s)\n", function_name, file_name);
    return -1;
  }

  /* get id for float3d */
  status = nc_inq_varid (ncid, var_name, &id_float3d);
  if (status != NC_NOERR) {
    fprintf (stderr,"Error, get id for '%s' field in the netCDF result file %s, in %s (%s) \n", var_name, file, function_name, file_name);
    return -1;
  }

  /* read values from netCDF variable */
  status = nc_get_var_float(ncid, id_float3d, tmp_float3d);
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d reading '%s' from %s, in %s (%s)\n", status, var_name, file, function_name, file_name);
    return status;
  }


  for (t=0;t<n1;t++)
    for (j=0;j<n2;j++)
      for (i=0;i<n3;i++) { 
        (*float3d)[t][j][i] = tmp_float3d[ t*n2*n3 + j*n3 + i ];
          /* fprintf (stderr, "t=%4d, z=%4d, y=%4d, x=%4d, float3d=%8.2f \n",t,c2,j,i,(*float3d)[t][c2][j][i]); */
          /* fflush(stderr); */
        }

  free(tmp_float3d);

  return 0;

#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif

}


/***********************************************************************************/
/* Function: read_4d_float                                                         */
/* Description:                                                                    */
/*  Reads a 4D data set from an netCDF file and stores it in an 4D float array     */
/*  will will be accordingly allocated                                             */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the read netCDF file                        */
/*     char *file         filename of the read netCDF file                         */
/*     char *var_name     short name of the data set, e.g. "lat"                   */
/*     long n1            number of entries in the 1st dimension                   */
/*     long n2            number of entries in the 2nd dimension                   */
/*     long n3            number of entries in the 3rd dimension                   */
/*     long n4            number of entries in the 4th dimension                   */
/* Return value:                                                                   */
/*     float *****float4d read data                                                */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int read_4d_float ( int ncid, char *file, char *var_name, float *****float4d, 
                    long n1, long n2, long n3, long n4, char ***done, int verbose )
{

#if HAVE_LIBNETCDF

  int status=0;
  int t=0, c2=0, j=0, i=0;
  int id_float4d=0;
  float *tmp_float4d = NULL;

  char function_name[]="read_4d_float";
  char file_name[]="netCDF_functions.c";

  if (verbose)
    fprintf(stderr," ... read_float_4d: var = %-10s, data_type = %-8s, ndim = %3d \n", 
                         var_name, "float", 4);


  /* allocate space for "float4d" */
  if (((*float4d)  = calloc (n1, sizeof (float ***)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'float4d' in %s (%s)\n", function_name, file_name);
    return -1;
  }

  for (t=0;t<n1;t++) {
    if (((*float4d)[t]  = calloc (n2, sizeof (float **)))   == NULL) {
      fprintf (stderr,"Error, allocating memory for 'float4d' in %s (%s)\n", function_name, file_name);
      return -2;
    }

    for (c2=0;c2<n2;c2++) {
      if (((*float4d)[t][c2] = calloc (n3, sizeof (float *)))   == NULL) {
        fprintf (stderr,"Error, allocating memory for 'float4d' in %s (%s)\n", function_name, file_name);
        return -3;
      }

      for (j=0;j<n3;j++) {
        if (((*float4d)[t][c2][j] = calloc (n4, sizeof (float)))   == NULL) {
          fprintf (stderr,"Error, allocating memory for 'float4d' in %s (%s)\n", function_name, file_name);
          return -4;
        }
      }
    }
  }

  /* allocate space for "tmp_float4d" */
  if ( ( tmp_float4d  = calloc (n1*n2*n3*n4, sizeof (float)) ) == NULL ) {
    fprintf (stderr,"Error, allocating memory for 'tmp_float4d' in %s (%s)\n", function_name, file_name);
    return -1;
  }

  /* get id for float4d */
  status = nc_inq_varid (ncid, var_name, &id_float4d);
  if (status != NC_NOERR) {
    fprintf (stderr,"Error, get id for '%s' field in the netCDF result file %s, in %s (%s) \n", var_name, file, function_name, file_name);
    return -1;
  }

  /* read values from netCDF variable */
  status = nc_get_var_float(ncid, id_float4d, tmp_float4d);
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d reading '%s' from %s, in %s (%s)\n", status, var_name, file, function_name, file_name);
    return status;
  }


  for (t=0;t<n1;t++)
    for (c2=0;c2<n2;c2++)
      for (j=0;j<n3;j++)
        for (i=0;i<n4;i++) { 
          (*float4d)[t][c2][j][i] = tmp_float4d[ t*n2*n3*n4 + c2*n3*n4 + j*n4 + i ];
            /* fprintf (stderr, "t=%4d, z=%4d, y=%4d, x=%4d, float4d=%8.2f \n",t,c2,j,i,(*float4d)[t][c2][j][i]); */
            /* fflush(stderr); */
          }

  free(tmp_float4d);

  return 0;

#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif

}

/***********************************************************************************/
/* Function: read_5d_float                                                         */
/* Description:                                                                    */
/*  Reads a 5D data set from an netCDF file and stores it in an 5D float array     */
/*  will will be accordingly allocated                                             */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the read netCDF file                        */
/*     char *file         filename of the read netCDF file                         */
/*     char *var_name     short name of the data set, e.g. "lat"                   */
/*     long n1             number of entries in the 1st dimension                  */
/*     long n2             number of entries in the 2nd dimension                  */
/*     long n3             number of entries in the 3rd dimension                  */
/*     long n4             number of entries in the 4th dimension                  */
/*     long n5             number of entries in the 5th dimension                  */
/* Return value:                                                                   */
/*     float ******float5d read data                                               */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/*                                                                                 */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/


int read_5d_float ( int ncid, char *file, char *var_name, float ******float5d, 
                    long n1, long n2, long n3, long n4, long n5, char ***done, int verbose )
{

#if HAVE_LIBNETCDF

  int status=0;
  int t=0, c2=0, c3=0, j=0, i=0;
  int id_float5d=0;
  float *tmp_float5d = NULL;

  char function_name[]="read_5d_float";
  char file_name[]="netCDF_functions.c";

  if (verbose)
    fprintf(stderr," ... read_float_5d: var = %-10s, data_type = %-8s, ndim = %3d \n", 
                         var_name, "float", 5);


  /* allocate space for "float5d" */
  if (((*float5d)  = calloc (n1, sizeof (float ****)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'float5d' in %s (%s)\n", function_name, file_name);
    return -1;
  }

  for (t=0;t<n1;t++) {
    if (((*float5d)[t]  = calloc (n2, sizeof (float ***)))   == NULL) {
      fprintf (stderr,"Error, allocating memory for 'float5d' in %s (%s)\n", function_name, file_name);
      return -2;
    }

    for (c2=0;c2<n2;c2++) {
      if (((*float5d)[t][c2] = calloc (n3, sizeof (float **)))   == NULL) {
        fprintf (stderr,"Error, allocating memory for 'float5d' in %s (%s)\n", function_name, file_name);
        return -3;
      }

      for (c3=0;c3<n3;c3++) {
        if (((*float5d)[t][c2][c3] = calloc (n4, sizeof (float *)))   == NULL) {
          fprintf (stderr,"Error, allocating memory for 'float5d' in %s (%s)\n", function_name, file_name);
          return -4;
        }

        for (j=0;j<n4;j++) {
          if (((*float5d)[t][c2][c3][j] = calloc (n5, sizeof (float)))   == NULL) {
            fprintf (stderr,"Error, allocating memory for 'float5d' in %s (%s)\n", function_name, file_name);
            return -5;
          }
        }
      }
    }
  }

  /* allocate space for "tmp_float5d" */
  if ( ( tmp_float5d  = calloc (n1*n2*n3*n4*n5, sizeof (float)) ) == NULL ) {
    fprintf (stderr,"Error, allocating memory for 'tmp_float5d' in %s (%s)\n", function_name, file_name);
    return -1;
  }

  /* get id for float5d */
  status = nc_inq_varid (ncid, var_name, &id_float5d);
  if (status != NC_NOERR) {
    fprintf (stderr,"Error, get id for '%s' field in the netCDF result file %s, in %s (%s) \n", var_name, file, function_name, file_name);
    return -1;
  }

  /* read values from netCDF variable */
  status = nc_get_var_float(ncid, id_float5d, tmp_float5d);
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d reading '%s' from %s\n", status, var_name, file);
    return status;
  }


  for (t=0;t<n1;t++)
    for (c2=0;c2<n2;c2++)
      for (c3=0;c3<n3;c3++)
        for (j=0;j<n4;j++)
          for (i=0;i<n5;i++) { 
            (*float5d)[t][c2][c3][j][i] = tmp_float5d[ t*n2*n3*n4*n5 + c2*n3*n4*n5 + c3*n4*n5 + j*n5 + i ];
            /* fprintf (stderr, "t=%4d, lambda=%4d, z=%4d, y=%4d, x=%4d, float5d=%8.2f \n",t,c2,c3,j,i,(*float5d)[t][c2][c3][j][i]); */
            /* fflush(stderr); */
          }

  free(tmp_float5d);

  return 0;

#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif

}

/***********************************************************************************/
/* Function: read_6d_float                                                         */
/* Description:                                                                    */
/*  Reads a 6D data set from an netCDF file and stores it in an 6D float array     */
/*  will will be accordingly allocated                                             */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the read netCDF file                        */
/*     char *file         filename of the read netCDF file                         */
/*     char *var_name     short name of the data set, e.g. "lat"                   */
/*     long n1              number of entries in the 1st dimension                 */
/*     long n2              number of entries in the 2nd dimension                 */
/*     long n3              number of entries in the 3rd dimension                 */
/*     long n4              number of entries in the 4th dimension                 */
/*     long n5              number of entries in the 5th dimension                 */
/*     long n6              number of entries in the 6th dimension                 */
/* Return value:                                                                   */
/*     float *******float6d read data                                              */
/*     int status         == 0, if everthing is OK                                 */
/*                                                                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int read_6d_float ( int ncid, char *file, char *var_name, float *******float6d, 
                    long n1, long n2, long n3, long n4, long n5, long n6, char ***done, int verbose )
{

#if HAVE_LIBNETCDF

  int status=0;
  int t=0, c2=0, c3=0, c4=0, j=0, i=0;
  int id_float6d=0;
  float *tmp_float6d = NULL;

  char function_name[]="read_6d_float";
  char file_name[]="netCDF_functions.c";

  if (verbose)
    fprintf(stderr," ... read_float_6d: var = %-10s, data_type = %-8s, ndim = %3d \n", 
                         var_name, "float", 6);


  /* allocate space for "float6d" */
  if (((*float6d)  = calloc (n1, sizeof (float *****)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'float6d' in %s (%s)\n", function_name, file_name);
    return -1;
  }

  for (t=0;t<n1;t++) {
    if (((*float6d)[t]  = calloc (n2, sizeof (float ****)))   == NULL) {
      fprintf (stderr,"Error, allocating memory for 'float6d' in %s (%s)\n", function_name, file_name);
      return -2;
    }

    for (c2=0;c2<n2;c2++) {
      if (((*float6d)[t][c2] = calloc (n3, sizeof (float ***)))   == NULL) {
        fprintf (stderr,"Error, allocating memory for 'float6d' in %s (%s)\n", function_name, file_name);
        return -3;
      }

      for (c3=0;c3<n3;c3++) {
        if (((*float6d)[t][c2][c3] = calloc (n4, sizeof (float **)))   == NULL) {
          fprintf (stderr,"Error, allocating memory for 'float6d' in %s (%s)\n", function_name, file_name);
          return -4;
        }

        for (c4=0;c4<n4;c4++) {
          if (((*float6d)[t][c2][c3][c4] = calloc (n5, sizeof (float *)))   == NULL) {
            fprintf (stderr,"Error, allocating memory for 'float6d' in %s (%s)\n", function_name, file_name);
            return -5;
          }
          for (j=0;j<n5;j++) {
            if (((*float6d)[t][c2][c3][c4][j] = calloc (n6, sizeof (float)))   == NULL) {
              fprintf (stderr,"Error, allocating memory for 'float6d' in %s (%s)\n", function_name, file_name);
              return -6;
            }
          }
        }
      }
    }
  }

  /* allocate space for "tmp_float6d" */
  if ( ( tmp_float6d  = calloc (n1*n2*n3*n4*n5*n6, sizeof (float)) ) == NULL ) {
    fprintf (stderr,"Error, allocating memory for 'tmp_float6d' in %s (%s)\n", function_name, file_name);
    return -1;
  }

  /* get id for float6d */
  status = nc_inq_varid (ncid, var_name, &id_float6d);
  if (status != NC_NOERR) {
    fprintf (stderr,"Error, get id for '%s' field in the netCDF result file %s, in %s (%s)\n", var_name, file, function_name, file_name);
    return -1;
  }

  /* read values from netCDF variable */
  status = nc_get_var_float(ncid, id_float6d, tmp_float6d);
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d reading '%s' from %s\n", status, var_name, file);
    return status;
  }


  for (t=0;t<n1;t++)
    for (c2=0;c2<n2;c2++)
      for (c3=0;c3<n3;c3++)
        for (c4=0;c4<n4;c4++)
          for (j=0;j<n5;j++)
            for (i=0;i<n6;i++) {
              /* if (done[t][j][i]) { */
              (*float6d)[t][c2][c3][c4][j][i] = tmp_float6d[ t*n2*n3*n4*n5*n6 + c2*n3*n4*n5*n6 + c3*n4*n5*n6 + c4*n5*n6 + j*n6 +i ];
              /* fprintf (stderr, "t=%4d, lambda=%4d, z=%4d, umu=%4d, y=%4d, x=%4d, float6d=%8.2f \n",t,c2,c3,c4,j,i,(*float6d)[t][c2][c3][c4][j][i]); */
              /* fflush(stderr); */
              /* } */
            }

  free(tmp_float6d);

  return 0;

#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif

}

/***********************************************************************************/
/* Function: read_7d_float                                                         */
/* Description:                                                                    */
/*  Reads a 7D data set from an netCDF file and stores it in an 7D float array     */
/*  will will be accordingly allocated                                             */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the read netCDF file                        */
/*     char *file         filename of the read netCDF file                         */
/*     char *var_name     short name of the data set, e.g. "lat"                   */
/*     long n1              number of entries in the 1st dimension                 */
/*     long n2              number of entries in the 2nd dimension                 */
/*     long n3              number of entries in the 3rd dimension                 */
/*     long n4              number of entries in the 4th dimension                 */
/*     long n5              number of entries in the 5th dimension                 */
/*     long n6              number of entries in the 6th dimension                 */
/*     long n7              number of entries in the 7th dimension                 */
/* Return value:                                                                   */
/*     float *******float6d read data                                              */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/*                                                                                 */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int read_7d_float ( int ncid, char *file, char *var_name, float ********float7d, 
                    long n1, long n2, long n3, long n4, long n5, long n6, long n7, 
                    char ***done, int verbose )
{

#if HAVE_LIBNETCDF

  int status=0;
  int t=0, c2=0, c3=0, c4=0, c5=0, j=0, i=0;
  int id_float7d=0;
  float *tmp_float7d = NULL;

  char function_name[]="read_7d_float";
  char file_name[]="netCDF_functions.c";

  if (verbose)
    fprintf(stderr," ... read_float_7d: var = %-10s, data_type = %-8s, ndim = %3d \n", 
                         var_name, "float", 7);


  /* allocate space for "float7d" */
  if (((*float7d)  = calloc (n1, sizeof (float ******)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'float7d' in %s (%s)\n", function_name, file_name);
    return -1;
  }

  for (t=0;t<n1;t++) {
    if (((*float7d)[t]  = calloc (n2, sizeof (float *****)))   == NULL) {
      fprintf (stderr,"Error, allocating memory for 'float7d' in %s (%s)\n", function_name, file_name);
      return -2;
    }

    for (c2=0;c2<n2;c2++) {
      if (((*float7d)[t][c2] = calloc (n3, sizeof (float ****)))   == NULL) {
        fprintf (stderr,"Error, allocating memory for 'float7d' in %s (%s)\n", function_name, file_name);
        return -3;
      }

      for (c3=0;c3<n3;c3++) {
        if (((*float7d)[t][c2][c3] = calloc (n4, sizeof (float ***)))   == NULL) {
          fprintf (stderr,"Error, allocating memory for 'float7d' in %s (%s)\n", function_name, file_name);
          return -4;
        }

        for (c4=0;c4<n4;c4++) {
          if (((*float7d)[t][c2][c3][c4] = calloc (n5, sizeof (float **)))   == NULL) {
            fprintf (stderr,"Error, allocating memory for 'float7d' in %s (%s)\n", function_name, file_name);
            return -5;
          }
          for (c5=0;c5<n5;c5++) {
            if (((*float7d)[t][c2][c3][c4][c5] = calloc (n6, sizeof (float *)))   == NULL) {
              fprintf (stderr,"Error, allocating memory for 'float7d' in %s (%s)\n", function_name, file_name);
              return -6;
            }
            for (j=0;j<n6;j++) {
              if (((*float7d)[t][c2][c3][c4][c5][j] = calloc (n7, sizeof (float)))   == NULL) {
                fprintf (stderr,"Error, allocating memory for 'float7d' in %s (%s)\n", function_name, file_name);
                return -7;
              }
            }
          }
        }
      }
    }
  }

  /* allocate space for "tmp_float7d" */
  if ( ( tmp_float7d  = calloc (n1*n2*n3*n4*n5*n6*n7, sizeof (float)) ) == NULL ) {
    fprintf (stderr,"Error, allocating memory for 'tmp_float7d' in %s (%s)\n", function_name, file_name);
    return -1;
  }

  /* get id for float7d */
  status = nc_inq_varid (ncid, var_name, &id_float7d);
  if (status != NC_NOERR) {
    fprintf (stderr,"Error, get id for '%s' field in the netCDF result file %s, in %s (%s) \n", var_name, file, function_name, file_name);
    return -1;
  }

  /* read values from netCDF variable */
  status = nc_get_var_float(ncid, id_float7d, tmp_float7d);
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error %d reading '%s' from %s\n", status, var_name, file);
    return status;
  }


  for (t=0;t<n1;t++)
    for (c2=0;c2<n2;c2++)
      for (c3=0;c3<n3;c3++)
        for (c4=0;c4<n4;c4++)
          for (c5=0;c5<n5;c5++)
            for (j=0;j<n6;j++)
              for (i=0;i<n7;i++) {
                /* if (done[t][j][i]) { */
                (*float7d)[t][c2][c3][c4][c5][j][i] = tmp_float7d[ (((((t*n2 + c2)*n3 + c3)*n4 + c4)*n5 + c5)*n6 + j)*n7 + i ];
                /* fprintf (stderr, "t=%4d, lambda=%4d, z=%4d, umu=%4d, y=%4d, x=%4d, float7d=%8.2f \n",t,c2,c3,c4,c5,j,(*float7d)[t][c2][c3][c4][c5][j][i]); */
                /* fflush(stderr); */
                /* } */
              }

  free(tmp_float7d);

  return 0;

#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif

}

/***********************************************************************************/
/* Function: read_3d_char                                                          */
/* Description:                                                                    */
/*  Reads a 3D data set from an netCDF file and stores it in an 3D char array      */
/*  will will be accordingly allocated.                                            */
/*                                                                                 */
/* Parameters:                                                                     */
/*     int ncid           id number of the read netCDF file                        */
/*     char *file         filename of the read netCDF file                         */
/*     char *var_name     short name of the data set, e.g. "lat"                   */
/*     long ntime         number of entries in the 1st dimension                   */
/*     long nlat          number of entries in the 2nd dimension                   */
/*     long nlon          number of entries in the 3rd dimension                   */
/*                                                                                 */
/* Return value:                                                                   */
/*     char ****done      read data                                                */
/*     int status         == 0, if everthing is OK                                 */
/*                        < 0, if there was an error                               */
/* Example:                                                                        */
/* Files:    netCDF_functions.c                                                    */
/* Known bugs: -                                                                   */
/* Author:                                                                         */
/*    Oct 2007   U. Hamann     Created                                             */
/*                                                                                 */
/***********************************************************************************/

int read_3d_char(int ncid, char *file, char *var_name, char ****done, long ntime, long nlat, long nlon, int verbose )
{

#if HAVE_LIBNETCDF

  int status=0;
  long t=0, j=0, i=0;
  int id_done=0;
  signed char *tmp_done = NULL;

  char function_name[]="read_3d_char";
  char file_name[]="netCDF_functions.c";

  if (verbose)
    fprintf(stderr," ... read_char_3d: var = %-10s, data_type = %-8s, ndim = %3d \n", 
                         var_name, "char", 3);

  /* allocate space for "done" */
  if (((*done)  = calloc ((ntime), sizeof (char **)))   == NULL) {
    fprintf (stderr,"Error, allocating memory for 'done' in %s (%s)\n", function_name, file_name);
    return -1;
  }

  for (t=0;t<(ntime);t++) {
    if (((*done)[t]  = calloc ((nlat), sizeof (char *)))   == NULL) {
      fprintf (stderr,"Error, allocating memory for 'done' in %s (%s)\n", function_name, file_name);
      return -1;
    }

    for (j=0;j<(nlat);j++)
      if (((*done)[t][j] = calloc ((nlon), sizeof (char)))   == NULL) {
        fprintf (stderr,"Error, allocating memory for 'done' in %s (%s)\n", function_name, file_name);
        return -1;
      }
  }

  /* allocate space for "tmp_done" */
  if ( ( tmp_done  = calloc ((ntime)*(nlon)*(nlat), sizeof (char)) ) == NULL ) {
    fprintf (stderr,"Error, allocating memory for 'tmp_done' in %s (%s)\n", function_name, file_name);
    return -1;
  }

  /* get id for done */
  status = nc_inq_varid (ncid, "done", &id_done);
  if (status != NC_NOERR) {
    /* this is in most of the cases OK, that's the start of an worldloop job */
    /* fprintf (stderr,"Error, get id for 'done' field in the netCDF result file %s, in %s (%s)\n", file, function_name, file_name); */
    return status;
  }

  /* read values from netCDF variable */
  status = nc_get_var_schar(ncid, id_done, tmp_done);
  if (status!=NC_NOERR) {
    fprintf (stderr, "Error '%s' while reading '%s' from %s, in %s (%s)\n", nc_strerror(status), var_name, file, function_name, file_name);
    return status;
  }

  for (t=0;t<(ntime);t++)
    for (j=0;j<(nlat);j++)
      for (i=0;i<(nlon);i++) {
        (*done)[t][j][i] = tmp_done[ t*(nlon)*(nlat) + j*(nlon) + i];
        /* if ((*done)[t][j][i]) { */
        /*   fprintf (stderr, "read_3d_char: t=%3ld j=%3ld i=%3ld, done=%d \n",t,j,i,(*done)[t][j][i]); */
        /*   fflush(stderr); } */
      }

  free(tmp_done);

  return 0;

#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif

}

/***********************************************************************************/
/* Function: write_netcdf_type                                                     */
/* Description: Write a variable to netCDF file; variable might be either          */
/*        int, float or double with specific dimension.                            */
/*        rows is the number of (not empty) rows of the file, columns is the       */
/*        number of columns of the file and columns is the number of columns       */
/*        of all (not empty) lines. If variable is 3D with regular grid size,      */
/*        length is the number of entries of each column. If variable is 3d        */
/*        with irregular grid size, length is a twodimensional array of int        */
/* Parameters:                                                                     */
/*  int  ncid:          Id of netCDF file give by 				   */
/*			nc_create(char *filename, int cmode, int &ncid)            */
/*  value:       	Pointer to a variable                                      */
/*  int  rows:          Number of rows                                             */
/*  int  columns:       Number of columns                                          */
/*  int  length:        Length of column                                           */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                         Bettina Richter */
/***********************************************************************************/

int write_netcdf_int(int ncid, int value, char *var_name) {
#if HAVE_LIBNETCDF
  if (ncid == 0) return -1;
  int varid, retval;

  if ((retval = nc_redef(ncid)))      ERR(retval);
   /********** Define the variable. The type of the variable in this case is NC_INT (4-byte integer). *********/
  if ((retval = nc_def_var(ncid, var_name, NC_INT,0, 0, &varid)))	ERR(retval);
   /********** End define mode. This tells netCDF we are done defining metadata. **********/
  if ((retval = nc_enddef(ncid)))      ERR(retval);
   /********** Write the pretend data to the file. Although netCDF supports
    * reading and writing subsets of data, in this case we write all
    * the data in one operation. **********/
  if ((retval = nc_put_var_int(ncid, varid, &value)))	ERR(retval);

  return 0;
#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif
}

int write_netcdf_float(int ncid, float value, char *var_name) {
#if HAVE_LIBNETCDF
  if (ncid == 0) return -1;
  int varid, retval;

  if ((retval = nc_redef(ncid)))      ERR(retval);
  if ((retval = nc_def_var(ncid, var_name, NC_FLOAT,0, 0, &varid)))	ERR(retval);
  if ((retval = nc_enddef(ncid)))      ERR(retval);
  if ((retval = nc_put_var_float(ncid, varid, &value)))	ERR(retval);

  return 0;
#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif
}

int write_netcdf_1Dint(int ncid, int *value, int rows, char *var_name, char *row_name) {
#if HAVE_LIBNETCDF
  if (ncid == 0) return -1;
  int retval, varid, row_dimid;

  if ((retval = nc_redef(ncid)))      ERR(retval);
  if (( nc_inq_dimid(ncid, row_name, &row_dimid) != NC_NOERR )) 	{ if ((retval = nc_def_dim( ncid, row_name, rows, &row_dimid)))			ERR(retval); }
  int dimid[1]={row_dimid};
  if ((retval = nc_def_var(ncid, var_name, NC_INT,1,dimid, &varid)))	ERR(retval);
  if ((retval = nc_enddef(ncid)))      ERR(retval);
  if ((retval = nc_put_var_int(ncid, varid, &value[0])))	ERR(retval);

  return 0;
#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif
}

int write_netcdf_1Dfloat(int ncid, float *value, int rows, char *var_name, char *row_name) {
#if HAVE_LIBNETCDF
  if (ncid == 0) return -1;
  int retval, varid, row_dimid;

  if ((retval = nc_redef(ncid)))      ERR(retval);
  if (( nc_inq_dimid(ncid, row_name, &row_dimid) != NC_NOERR )) 	{ if ((retval = nc_def_dim( ncid, row_name, rows, &row_dimid)))			ERR(retval); }
   /********** The dimids array is used to pass the IDs of the dimensions of the variable. *********/
  int dimid[1]={row_dimid};
  if ((retval = nc_def_var(ncid, var_name, NC_FLOAT,1,dimid, &varid)))	ERR(retval);
  if ((retval = nc_enddef(ncid)))      ERR(retval);
  if ((retval = nc_put_var_float(ncid, varid, &value[0])))	ERR(retval);

  return 0;
#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif
}

int write_netcdf_2Dint(int ncid, int **value, int rows, int columns, char *var_name, char *row_name, char *column_name) {
#if HAVE_LIBNETCDF
  if (ncid == 0) return -1;

  int value_data[rows][columns];
  int row, column;
  for (row=0; row<rows; row++) {
    for (column=0; column<columns; column++)
	value_data[row][column] = value[row][column];
  }

  int retval, varid, row_dimid, column_dimid;
  int NDIMS=2;

  if ((retval = nc_redef(ncid)))      ERR(retval);
  if (( nc_inq_dimid(ncid, row_name, &row_dimid) != NC_NOERR )) 	{ if ((retval = nc_def_dim( ncid, row_name, rows, &row_dimid)))			ERR(retval); }
  if (( nc_inq_dimid(ncid, column_name, &column_dimid) != NC_NOERR )) 	{ if ((retval = nc_def_dim( ncid, column_name, columns, &column_dimid)))	ERR(retval); }
  int dimids[2] = {column_dimid, row_dimid};
  if ((retval = nc_def_var(ncid, var_name, NC_INT, NDIMS, dimids, &varid)))   ERR(retval);
  if ((retval = nc_enddef(ncid)))      ERR(retval);
  if ((retval = nc_put_var_int(ncid, varid, &value_data[0][0])))      ERR(retval);
 
  return 0;
#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif
}

int write_netcdf_2Dfloat(int ncid, float **value, int rows, int columns, char *var_name, char *row_name, char *column_name) {
#if HAVE_LIBNETCDF
  if (ncid == 0) return -1;

  float value_data[rows][columns];
  int row, column;
  for (row=0; row<rows; row++) {
    for (column=0; column<columns; column++)
	value_data[row][column] = value[row][column];
  }

  int retval, varid, row_dimid, column_dimid;
  int NDIMS=2;

  if ((retval = nc_redef(ncid)))      ERR(retval);
  if (( nc_inq_dimid(ncid, row_name, &row_dimid) != NC_NOERR )) 	{ if ((retval = nc_def_dim( ncid, row_name, rows, &row_dimid)))			ERR(retval); }
  if (( nc_inq_dimid(ncid, column_name, &column_dimid) != NC_NOERR )) 	{ if ((retval = nc_def_dim( ncid, column_name, columns, &column_dimid)))	ERR(retval); }
  int dimids[2] = {column_dimid, row_dimid};
  if ((retval = nc_def_var(ncid, var_name, NC_FLOAT, NDIMS, dimids, &varid)))   ERR(retval);
  if ((retval = nc_enddef(ncid)))      ERR(retval);
  if ((retval = nc_put_var_float(ncid, varid, &value_data[0][0])))      ERR(retval);
 
  return 0;
#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif
}

int write_netcdf_2Ddouble(int ncid, double **value, int rows, int columns, char *var_name, char *row_name, char *column_name) {
#if HAVE_LIBNETCDF
  if (ncid == 0) return -1;

  double value_data[rows][columns];
  int row, column;
  for (row=0; row<rows; row++) {
    for (column=0; column<columns; column++)
	value_data[row][column] = value[row][column];
  }

  int retval, varid, row_dimid, column_dimid;
  int NDIMS=2;

  if ((retval = nc_redef(ncid)))      ERR(retval);
  if (( nc_inq_dimid(ncid, row_name, &row_dimid) != NC_NOERR )) 	{ if ((retval = nc_def_dim( ncid, row_name, rows, &row_dimid)))			ERR(retval); }
  if (( nc_inq_dimid(ncid, column_name, &column_dimid) != NC_NOERR )) 	{ if ((retval = nc_def_dim( ncid, column_name, columns, &column_dimid)))	ERR(retval); }
  int dimids[2] = {column_dimid, row_dimid};
  if ((retval = nc_def_var(ncid, var_name, NC_DOUBLE, NDIMS, dimids, &varid)))   ERR(retval);
  if ((retval = nc_enddef(ncid)))      ERR(retval);
  if ((retval = nc_put_var_double(ncid, varid, &value_data[0][0])))      ERR(retval);
 
  return 0;
#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif
}


int write_netcdf_3Dfloat(int ncid, float ***value, int rows, int columns, int length, char *var_name, char *row_name, char *column_name, char *length_name) {
#if HAVE_LIBNETCDF

  if (ncid == 0) return -1;

  float value_data[rows][columns][length];
  int row, column, il;
  for (row=0; row<rows; row++) {
    for (column=0; column<columns; column++) {
      for (il=0; il<length; il++)
	value_data[row][column][il] = value[row][column][il];
    }
  }

  int retval, varid, row_dimid, column_dimid, length_dimid;
  int NDIMS=3;
  if ((retval = nc_redef(ncid)))      ERR(retval);
  if (( nc_inq_dimid(ncid, row_name, &row_dimid) != NC_NOERR )) 	{ if ((retval = nc_def_dim( ncid, row_name, rows, &row_dimid)))			ERR(retval); }
  if (( nc_inq_dimid(ncid, column_name, &column_dimid) != NC_NOERR )) 	{ if ((retval = nc_def_dim( ncid, column_name, columns, &column_dimid)))	ERR(retval); }
  if (( nc_inq_dimid(ncid, length_name, &length_dimid) != NC_NOERR )) 	{ if ((retval = nc_def_dim( ncid, length_name, length, &length_dimid)))	 	ERR(retval); }

  int dimids[3] = {length_dimid, column_dimid, row_dimid};
  if ((retval = nc_def_var(ncid, var_name, NC_FLOAT, NDIMS, dimids, &varid)))   ERR(retval);
  if ((retval = nc_enddef(ncid)))      ERR(retval);
  if ((retval = nc_put_var_float(ncid, varid, &value_data[0][0][0])))      ERR(retval);
 
  return 0;
#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif
}

int write_netcdf_3Dirr_row_float(int ncid, float ***value, int rows, int columns, int *length, char *var_name, char *row_name, char *column_name, char *length_name) {
#if HAVE_LIBNETCDF

  if (ncid == 0) return -1;

  int row, column, il;
  /**** Maximum length of value; in case of regular grid, length_max == length ****/
  int length_max;
  length_max=0;
  for (row=0; row<rows; row++) {
    if ( length_max < length[row] ) length_max = length[row];
  }

  /**** Copy ***value into a block of contiguous data values for netcdf ****/ 
  float value_data[rows][columns][length_max];
  for (row=0; row<rows; row++) {
    for (column=0; column<columns; column++) {
      for (il=0; il<length_max; il++)
	if ( il < length[row] ) value_data[row][column][il] = value[row][column][il];
	else value_data[row][column][il] = NOT_DEFINED_FLOAT;
    }
  }

  int retval, varid, row_dimid, column_dimid, length_dimid;
  int NDIMS=3;
  if ((retval = nc_redef(ncid)))      ERR(retval);

  if (( nc_inq_dimid(ncid, row_name, &row_dimid) != NC_NOERR )) 	{ if ((retval = nc_def_dim( ncid, row_name, rows, &row_dimid)))			ERR(retval); }
  if (( nc_inq_dimid(ncid, column_name, &column_dimid) != NC_NOERR )) 	{ if ((retval = nc_def_dim( ncid, column_name, columns, &column_dimid)))	ERR(retval); }
  if (( nc_inq_dimid(ncid, length_name, &length_dimid) != NC_NOERR )) 	{ if ((retval = nc_def_dim( ncid, length_name, length_max, &length_dimid)))	ERR(retval); }
  int dimids[3] = {length_dimid, column_dimid, row_dimid};
  if ((retval = nc_def_var(ncid, var_name, NC_FLOAT, NDIMS, dimids, &varid)))   ERR(retval);
  if ((retval = nc_enddef(ncid)))      ERR(retval);
  if ((retval = nc_put_var_float(ncid, varid, &value_data[0][0][0])))      ERR(retval);
 
  return 0;
#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif
}

int write_netcdf_3Dirrfloat(int ncid, float ***value, int rows, int columns, int **length, char *var_name, char *row_name, char *column_name, char *length_name) {
#if HAVE_LIBNETCDF

  if (ncid == 0) return -1;

  int row, column, il;
  /**** Maximum length of value; in case of regular grid, length_max == length ****/
  int length_max;
  length_max=0;
  for (row=0; row<rows; row++) {
    for (column=0; column<columns; column++) {
      if ( length_max < length[row][column] ) length_max = length[row][column];
    }
  }

  if ( length_max<=0 ) { fprintf(stderr, "ERROR: Could write %s to netcdf_file\n", var_name); return 0; }
  /**** Copy ***value into a block of contiguous data values for netcdf ****/ 
  float value_data[rows][columns][length_max];
  for (row=0; row<rows; row++) {
    for (column=0; column<columns; column++) {
      for (il=0; il<length_max; il++)
	if ( il < length[row][column] ) value_data[row][column][il] = value[row][column][il];
	else value_data[row][column][il] = NOT_DEFINED_FLOAT;
    }
  }

  int retval, varid, row_dimid, column_dimid, length_dimid;
  int NDIMS=3;
  if ((retval = nc_redef(ncid)))      ERR(retval);

  if (( nc_inq_dimid(ncid, row_name, &row_dimid) != NC_NOERR )) 	{ if ((retval = nc_def_dim( ncid, row_name, rows, &row_dimid)))			ERR(retval); }
  if (( nc_inq_dimid(ncid, column_name, &column_dimid) != NC_NOERR )) 	{ if ((retval = nc_def_dim( ncid, column_name, columns, &column_dimid)))	ERR(retval); }
  if (( nc_inq_dimid(ncid, length_name, &length_dimid) != NC_NOERR )) 	{ if ((retval = nc_def_dim( ncid, length_name, length_max, &length_dimid)))	ERR(retval); }
  int dimids[3] = {length_dimid, column_dimid, row_dimid};
  if ((retval = nc_def_var(ncid, var_name, NC_FLOAT, NDIMS, dimids, &varid)))   ERR(retval);
  if ((retval = nc_enddef(ncid)))      ERR(retval);
  if ((retval = nc_put_var_float(ncid, varid, &value_data[0][0][0])))      ERR(retval);
 
  return 0;
#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif
}

int write_netcdf_3Dirrdouble(int ncid, double ***value, int rows, int columns, int **length, char *var_name, char *row_name, char *column_name, char *length_name) {
#if HAVE_LIBNETCDF

  if (ncid == 0) return -1;

  int row, column, il;
  /**** Maximum length of value; in case of regular grid, length_max == length ****/
  int length_max;
  length_max=0;
  for (row=0; row<rows; row++) {
    for (column=0; column<columns; column++) {
      if ( length_max < length[row][column] ) length_max = length[row][column];
    }
  }
  if ( length_max<=0 ) { fprintf(stderr, "ERROR: Could write %s to netcdf_file\n", var_name); return 0; }

  /**** Copy ***value into a block of contiguous data values for netcdf ****/ 
  double value_data[rows][columns][length_max];
  for (row=0; row<rows; row++) {
    for (column=0; column<columns; column++) {
      for (il=0; il<length_max; il++)
	if ( il < length[row][column] ) value_data[row][column][il] = value[row][column][il];
	else value_data[row][column][il] = NOT_DEFINED_FLOAT;
    }
  }

  int retval, varid, row_dimid, column_dimid, length_dimid;
  int NDIMS=3;
  if ((retval = nc_redef(ncid)))      ERR(retval);

  if (( nc_inq_dimid(ncid, row_name, &row_dimid) != NC_NOERR )) 	{ if ((retval = nc_def_dim( ncid, row_name, rows, &row_dimid)))			ERR(retval); }
  if (( nc_inq_dimid(ncid, column_name, &column_dimid) != NC_NOERR )) 	{ if ((retval = nc_def_dim( ncid, column_name, columns, &column_dimid)))	ERR(retval); }
  if (( nc_inq_dimid(ncid, length_name, &length_dimid) != NC_NOERR )) 	{ if ((retval = nc_def_dim( ncid, length_name, length_max, &length_dimid)))	ERR(retval); }
  int dimids[3] = {length_dimid, column_dimid, row_dimid};
  if ((retval = nc_def_var(ncid, var_name, NC_DOUBLE, NDIMS, dimids, &varid)))   ERR(retval);
  if ((retval = nc_enddef(ncid)))      ERR(retval);
  if ((retval = nc_put_var_double(ncid, varid, &value_data[0][0][0])))      ERR(retval);
 
  return 0;
#else
    fprintf (stderr, " ***********************************************************************\n");
    fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
    fprintf (stderr, " * use any netCDF option. Please get netcdf and rebuild.               *\n");
    fprintf (stderr, " ***********************************************************************\n");
    return -1;
#endif
}




