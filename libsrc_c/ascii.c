/*--------------------------------------------------------------------
 * $Id: ascii.c 3279 2017-07-07 20:35:28Z Claudia.Emde $
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


/* @30c@ */
/* @c30@ */
/****************************************************************************/
/* In order to use the functions provided by the ascii library,   @30_10c@  */
/* #include <ascii.h> in your source code and link with libRadtran_c.a.     */
/*                                                                          */
/* @strong{Example:}                                                        */
/* Example for a source file:                                               */
/* @example                                                                 */
/*                                                                          */
/*   ...                                                                    */
/*   #include "../src_c/ascii.h"                                            */
/*   ...                                                                    */
/*                                                                          */
/* @end example                                                             */
/*                                                                          */
/* Linking of the executable, using the GNU compiler gcc:                   */
/* @example                                                                 */
/*                                                                          */
/*   gcc -o test test.c -lRadtran_c -L../lib                                */
/*                                                                          */
/* @end example                                                   @c30_10@  */
/****************************************************************************/


/****************************************************************************/
/* The ASCII library provides functions for parsing ASCII files   @30_20c@  */
/* containing arrays of data.                                               */
/* An ASCII file is read line by line. Each line is split into fields;      */
/* a field is an arbitrary combination of characters which does             */
/* neither contain the line separator ('CARRIAGE RETURN') nor the field     */
/* separator ('SPACE').                                                     */
/*                                                                          */
/* In detail:                                                               */
/* @itemize @asis                                                           */
/*   @item Lines are separated by 'CARRIAGE RETURN'.                        */
/*   @item Tokens (or columns) are separated by one or more 'SPACE's.       */
/*   @item Empty lines are simply ignored.                                  */
/*   @item \% and # are comment symbols; text between a comment symbol      */
/*     and the next line separator is ignored                               */
/*   @item A comment symbol which is not at the beginning of a line is      */
/*     only recognized after a field separator, but not within a field      */ 
/*   @item The number of fields may differ from line to line.               */   
/* @end itemize                                                             */
/*                                                                          */
/* A simple example for an ASCII file, which would be recognized as a       */
/* valid one-column or two-column ASCII file:                               */
/* @example                                                                 */
/*                                                                          */
/*   % This is an example for the input ASCII file for sdose,               */
/*   % the time integration program                                         */
/*   11.0     13.0   % the two hours around noon                            */
/*   10.0     14.0                                                          */
/*    9.0     15.0                                                          */
/*                                                                          */
/*   # total dose                                                           */
/*   -1.0     24.0    % integrate over maximum available time interval      */
/*                                                                          */
/*   # the following line shows, that an extra column does not matter       */
/*   2  3.4  17                                                             */
/*                                                                          */
/* @end example                                                             */
/*                                                                          */ 
/* For most purposes, ASCII_file2double and ASCII_free_double provide a     */
/* convenient way for parsing files.                                        */
/* @strong{Example:}                                                        */
/* @example                                                                 */
/*                                                                          */
/*  #include <stdio.h>                                                      */
/*  #include <ascii.h>                                                      */
/*                                                                          */ 
/*  int main(int argc, char ** argv)                                        */
/*  @{                                                                      */
/*    int rows=0, max_columns=0, min_columns=0;                             */
/*    int i=0, status=0;                                                    */
/*    double **value=NULL;                                                  */
/*                                                                          */
/*    status = ASCII_file2double ("test.dat",                               */
/*  	                          &rows,                                    */
/*  			          &max_columns,                             */
/*  			          &min_columns,                             */
/*  			          &value);                                  */
/*		                                                            */
/*    for (i=0; i<rows; i++) @{                                             */
/*      ... do something for each row of the matrix                         */
/*   @}                                                                     */
/*                                                                          */
/*    ASCII_free_double (value, rows);                                      */
/*                                                                          */
/*    return 0;                                                             */
/* @}                                                                       */
/*                                                                          */
/* @end example                                                             */
/*                                                                          */
/* For special purposes (ASCII files with 1,2,3, or 5 columns) there are    */
/* additionally functions read_1c_file, ... which facilitate the            */
/* access even more.                                                        */
/*                                                                 @c30_20@ */
/****************************************************************************/
									    
#include <ctype.h>							    
#include <stdio.h>							    
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <ascii.h>
#include <errno.h>


/**********************************/
/* Error codes                    */
/**********************************/

#define ASCIIFILE_NOT_FOUND      -1 
#define ASCII_NO_MEMORY          -2 
#define LESS_THAN_TWO_COLUMNS    -3 
#define LESS_THAN_THREE_COLUMNS  -4 
#define LESS_THAN_FOUR_COLUMNS   -5 
#define LESS_THAN_FIVE_COLUMNS   -6 



/********************************************/
/* Limits                                   */
/********************************************/

/* Maximum number of characters per line */
#define MAX_LENGTH_OF_LINE     1048576



/***********************************/
/* Character to indicate a comment */
/***********************************/

#define ASCII_COMMENT_1          '%'
#define ASCII_COMMENT_2          '#'


/*******************************************/
/* Structure needed by the ASCII_sortarray */
/* function.                               */
/*******************************************/

typedef struct {
  double *x;      /* a row of a matrix */ 
  int index;      /* sort index        */
}  sortarray_struct;


/*************************************************/
/* Structure needed by the ASCII_sortarray_float */
/* function.                                     */
/*************************************************/

typedef struct {
  float *x;       /* a row of a matrix */ 
  int index;      /* sort index        */
}  sortarray_float_struct;



/************************************/
/* Prototypes of internal functions */
/************************************/

static int ASCII_comment (char t);

static int sortarray_fnc_up (const void *x1, const void *x2);
static int sortarray_fnc_dn (const void *x1, const void *x2);


/***********************************************************************************/
/* Function: ASCII_checkfile_nrows                                        @30_30i@ */
/* Description:                                                                    */
/*  Check an ASCII file: count rows, minimum number of columns,                    */
/*  maximum number of columns and the maximum length of a string;                  */
/*  empty rows and characters after one of the comment symbols                     */
/*  (either ASCII_COMMENT_1 or ASCII_COMMENT_2) are ignored.                       */ 
/*  If parameter nrows is larger or equal 0, only the first non-empty,             */
/*  non-comment nrows are checked and the rest of the file is ignored.             */
/*                                                                                 */
/* Parameters:                                                                     */
/*  char *filename:    name of the file which should be checked                    */
/*  int  *rows:        number of rows found                                        */
/*  int  *min_columns: minimum number of columns, set by function                  */
/*  int  *max_columns: maximum number of columns, set by function                  */
/*  int  *max_length:  maximum length of a field, set by function                  */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i30_30@ */
/***********************************************************************************/

int ASCII_checkfile_nrows (char *filename, int nrows, 
			   int *rows,
			   int *min_columns,
			   int *max_columns,
			   int *max_length)
   
{
  FILE *f=NULL;
  char line[MAX_LENGTH_OF_LINE+1]="";
  char *token=NULL;
  char *string=NULL;
  int temp1=0, temp2=0;
  int min_col=INT_MAX, max_col=0, max_len=0, r=0;
  

  /* reset parameters */
  *min_columns=0;
  *max_columns=0;
  *max_length=0;
  *rows=0;

  string = line;

  
  if ( (f = fopen(filename, "r")) == NULL) { 
    fprintf (stderr, "Error, file '%s' not found!\n", filename);
    return ASCIIFILE_NOT_FOUND;
  }

  /* count rows and columns */
  while ( fgets (string, MAX_LENGTH_OF_LINE, f) != NULL )  {

    if ( (token = strtok (string, " \t\n")) != NULL ) { /* if not an empty line */ 
      /* if not a comment     */
      if (!ASCII_comment(token[0]))  {

	r++;
	temp1++;

	/* check for maximal string length */
	temp2 = strlen(token);
	max_len = (temp2>max_len ? temp2 : max_len);

	while ( (token = strtok (NULL, " \t\n")) != NULL)  {

	  /* check for maximal string length */
	  temp2 = strlen(token);
	  max_len = (temp2>max_len ? temp2 : max_len);

	  if (ASCII_comment(token[0]))  /* if comment */
	    break;
	  temp1++;
	}

	min_col = (temp1<min_col ? temp1 : min_col);
	max_col = (temp1>max_col ? temp1 : max_col);
	temp1 = 0;
      }
    }

    string = line;

    /* stop after reading nrows rows */
    if (nrows>=0 && r>=nrows)
      break;
  }

  if (r<=0)
    min_col=0; 

  *min_columns = min_col;
  *max_columns = max_col;
  *max_length  = max_len;
  *rows        = r; 
  

  fclose (f);
  return 0;
}




/***********************************************************************************/
/* Function: ASCII_checkfile                                              @30_30i@ */
/* Description:                                                                    */
/*  Check an ASCII file: count rows, minimum number of columns,                    */
/*  maximum number of columns and the maximum length of a string;                  */
/*  empty rows and characters after one of the comment symbols                     */
/*  (either ASCII_COMMENT_1 or ASCII_COMMENT_2) are ignored.                       */ 
/*                                                                                 */
/* Parameters:                                                                     */
/*  char *filename:    name of the file which should be checked                    */
/*  int  *rows:        number of rows found                                        */
/*  int  *min_columns: minimum number of columns, set by function                  */
/*  int  *max_columns: maximum number of columns, set by function                  */
/*  int  *max_length:  maximum length of a field, set by function                  */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                        @i30_30@ */
/***********************************************************************************/

int ASCII_checkfile (char *filename, 
		     int *rows,
		     int *min_columns,
		     int *max_columns,
		     int *max_length)
   
{
  int status=0;

  status = ASCII_checkfile_nrows (filename, -1, 
				  rows,
				  min_columns,
				  max_columns,
				  max_length);
  return status;
}




/***********************************************************************************/
/* Function: ASCII_calloc_string                                          @30_30i@ */
/* Description: Allocate memory for a two-dimensional array of strings.            */
/*                                                                                 */
/* Parameters:                                                                     */
/*  char ****string:  Pointer to a two-dimensional array of strings; memory        */
/*                    for string is allocated automatically                        */
/*  int rows:         Number of rows, specified by the caller                      */
/*  int columns:      Number of columns, specified by the caller                   */
/*  int length:       Maximum length of a string, specified by the caller          */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_calloc_string (char ****string,
			 int rows,
			 int columns,
			 int length)
{
  int i=0, j=0;

  if ((*string = (char ***) calloc (rows, sizeof (char **))) == NULL)
    return ASCII_NO_MEMORY;

  for (i=0; i<rows; i++)  {
    if (((*string)[i] = (char **) calloc (columns, sizeof (char *))) == NULL)
      return ASCII_NO_MEMORY;

    for (j=0; j<columns; j++)
      if (((*string)[i][j] = (char *) calloc (length+1, sizeof (char))) == NULL) 
	return ASCII_NO_MEMORY;
    
  }
  
  return 0;
}    




/***********************************************************************************/
/* Function: ASCII_free_string                                            @30_30i@ */
/* Description: Free memory, which has been allocated with ASCII_calloc_string.    */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_free_string (char ***string, 
		       int rows,
		       int columns)

{
  int i=0, j=0;

  for (i=0; i<rows; i++)  {
    for (j=0; j<columns; j++)
      free (string[i][j]);

    free (string[i]);
  }

  free (string);

  return 0;
}    


/***********************************************************************************/
/* Function: ASCII_calloc_char                                            @30_30i@ */
/* Description: Allocate memory for a two-dimensional array of chars  .            */
/*                                                                                 */
/* Parameters:                                                                     */
/*  char ***value:    Pointer to a two-dimensional array of chars; memory          */
/*                    for char is allocated automatically                        */
/*  int rows:         Number of rows, specified by the caller                      */
/*  int columns:      Number of columns, specified by the caller                   */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_calloc_char (char ***value,
       int rows,
       int columns)
{
  int i=0;

  if ((*value = (char **) calloc (rows, sizeof (char *))) == NULL)
    return ASCII_NO_MEMORY;

  for (i=0; i<rows; i++)  
    if ( ((*value)[i] = (char *) calloc (columns, sizeof (char))) == NULL)
      return ASCII_NO_MEMORY;

  return 0;
}


/***********************************************************************************/
/* Function: ASCII_free_char                                              @30_30i@ */
/* Description: Free memory, which has been allocated with ASCII_calloc_char.      */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_free_char (char **value, 
		       int rows)
{
  int i=0;

  for (i=0; i<rows; i++)  
    free (value[i]);
  free (value);

  return 0;
}


/***********************************************************************************/
/* Function: ASCII_calloc_int                                             @30_30i@ */
/* Description: Allocate memory for a two-dimensional array of int.                */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_calloc_int (int ***value,
		      int rows,
		      int columns)
{
  int i=0;
  
  if ( (*value = (int **) calloc ((size_t) rows, sizeof (int *))) == NULL )
    return ASCII_NO_MEMORY;

  for (i=0; i<rows; i++)  
    if ( ((*value)[i] = (int *) calloc ((size_t) columns, 
					   sizeof (int))) == NULL )
      return ASCII_NO_MEMORY;

  return 0;
}    


/***********************************************************************************/
/* Function: ASCII_calloc_int_3D                                          @30_30i@ */
/* Description: Allocate memory for a three-dimensional array of int.              */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_calloc_int_3D (int ****value,
                         int rows,
                         int columns,
                         int length)
{
  int i=0, j=0;
  
  if ( (*value = (int ***) calloc (rows, sizeof (int **))) == NULL )  
    return ASCII_NO_MEMORY;
  
    
  for (i=0; i<rows; i++) {
    if ( ((*value)[i] = (int **) calloc (columns, sizeof (int *))) == NULL )
      return ASCII_NO_MEMORY;
    for (j=0; j<columns; j++)
      if (((*value)[i][j] = (int *) calloc (length, sizeof (int))) == NULL) 
	return ASCII_NO_MEMORY;
    
  }

  return 0;
  
}    


/***********************************************************************************/
/* Function: ASCII_calloc_int_4D                                          @30_30i@ */
/* Description: Allocate memory for a four-dimensional array of integer.           */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                   @i30_30@ */
/***********************************************************************************/

int ASCII_calloc_int_4D  (int *****value,
			  int rows,
			  int columns,
			  int length,
			  int fourth_dimension)
{
  int i=0, j=0, k=0;

  if ( (*value = (int ****) calloc (rows, sizeof (int ***))) == NULL )  
    return ASCII_NO_MEMORY;

    
  for (i=0; i<rows; i++) {
    if ( ((*value)[i] = (int ***) calloc (columns, sizeof (int **))) == NULL )
      return ASCII_NO_MEMORY;
    for (j=0; j<columns; j++) {
      if (((*value)[i][j] = (int **) calloc (length, sizeof (int *))) == NULL) 
	return ASCII_NO_MEMORY;
      for (k=0; k<length; k++)
	if (((*value)[i][j][k] = (int *) calloc (fourth_dimension, sizeof (int))) == NULL) 
	  return ASCII_NO_MEMORY;
    }    
  }

  return 0;
}    


/***********************************************************************************/
/* Function: ASCII_calloc_double                                          @30_30i@ */
/* Description: Allocate memory for a two-dimensional array of double.             */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_calloc_double (double ***value,
			 int rows,
			 int columns)
{
  int i=0;
  
  if (rows>0) {
    if ( (*value = (double **) calloc ((size_t) rows, sizeof (double *))) == NULL )
      return ASCII_NO_MEMORY;

    if (columns>0)
      for (i=0; i<rows; i++)  
	if ( ((*value)[i] = (double *) calloc ((size_t) columns, 
					       sizeof (double))) == NULL )
	  return ASCII_NO_MEMORY;
  }

  return 0;
}    


/***********************************************************************************/
/* Function: ASCII_calloc_double_3D                                       @30_30i@ */
/* Description: Allocate memory for a three-dimensional array of double.           */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                   @i30_30@ */
/***********************************************************************************/

int ASCII_calloc_double_3D  (double ****value,
			 int rows,
			 int columns,
			 int length)
{
  int i=0, j=0;

  if ( (*value = (double ***) calloc (rows, sizeof (double **))) == NULL )  
    return ASCII_NO_MEMORY;

    
  for (i=0; i<rows; i++) {
    if ( ((*value)[i] = (double **) calloc (columns, sizeof (double *))) == NULL )
      return ASCII_NO_MEMORY;
    for (j=0; j<columns; j++)
      if (((*value)[i][j] = (double *) calloc (length, sizeof (double))) == NULL) 
	return ASCII_NO_MEMORY;
    
  }

  return 0;
}    


/***********************************************************************************/
/* Function: ASCII_calloc_double_3D_arylen                                @30_30i@ */
/* Description: Allocate memory for a three-dimensional array of double.           */
/*              index range of the last index is variable                          */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                   @i30_30@ */
/***********************************************************************************/

int ASCII_calloc_double_3D_arylen (double ****value,
				   int rows,
				   int columns,
				   int *length)
{
  int i=0, j=0;

  if ( (*value = (double ***) calloc (rows, sizeof (double **))) == NULL )  
    return ASCII_NO_MEMORY;

    
  for (i=0; i<rows; i++) {
    if ( ((*value)[i] = (double **) calloc (columns, sizeof (double *))) == NULL )
      return ASCII_NO_MEMORY;

    for (j=0; j<columns; j++)
      if (length[j]>0)  /* need to check because some OSs don't like calloc (0) */
	if (((*value)[i][j] = (double *) calloc (length[j], sizeof (double))) == NULL) 
	  return ASCII_NO_MEMORY;
  }

  return 0;
}


/***********************************************************************************/
/* Function: ASCII_calloc_double_4D_arylen                                @30_30i@ */
/* Description: Allocate memory for a four-dimensional array of double.            */
/*              index range of the last index is variable                          */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Ulrich Hamann                                                  @i30_30@ */
/***********************************************************************************/

int ASCII_calloc_double_4D_arylen (double *****value,
				   int rows,
				   int columns,
                                   int levels,
				   int *length)
{
  int i=0, j=0, k=0;

  if ( (*value = (double ****) calloc (rows, sizeof (double ***))) == NULL )  
    return ASCII_NO_MEMORY;

  for (i=0;i<rows;i++) {
    if ( ((*value)[i] = (double ***) calloc (columns, sizeof (double **))) == NULL )
        return ASCII_NO_MEMORY;
    
    for (j=0; j<columns; j++) {
      if ( ((*value)[i][j] = (double **) calloc (levels, sizeof (double *))) == NULL )
        return ASCII_NO_MEMORY;

      for (k=0; k<levels; k++)
        if (length[k]>0)  /* need to check because some OSs don't like calloc (0) */
	  if (((*value)[i][j][k] = (double *) calloc (length[k], sizeof (double))) == NULL) 
	    return ASCII_NO_MEMORY;
    }
  }

  return 0;
}    

/***********************************************************************************/
/* Function: ASCII_calloc_double_5D_arylen                                         */
/* Description: Allocate memory for a five-dimensional array of double.            */
/*              index range of the last index is variable                          */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Claudia Emde (adapted from ASCII_calloc_double_3D_arylen)               */
/***********************************************************************************/

int ASCII_calloc_double_5D_arylen (double ******value,
				   int dim1,
				   int dim2,
				   int dim3,
                                   int dim4,
				   int *dim5)
{
  int i=0, j=0, k=0, l=0;

  if ( (*value = (double *****) calloc (dim1, sizeof (double ****))) == NULL )  
    return ASCII_NO_MEMORY;

  for (i=0;i<dim1;i++) {
    if ( ((*value)[i] = (double ****) calloc (dim2, sizeof (double ***))) == NULL )
        return ASCII_NO_MEMORY;
    
    for (j=0; j<dim2; j++) {
      if ( ((*value)[i][j] = (double ***) calloc (dim3, sizeof (double **))) == NULL )
        return ASCII_NO_MEMORY;

      for (k=0; k<dim3; k++) {
	if ( ((*value)[i][j][k] = (double **) calloc (dim4, sizeof (double *))) == NULL )
	  return ASCII_NO_MEMORY;

	for (l=0; l<dim4; l++)
	  if (dim5[l]>0)  /* need to check because some OSs don't like calloc (0) */
	    if (((*value)[i][j][k][l] = (double *) calloc (dim5[l],
							   sizeof (double))) == NULL) 
	      return ASCII_NO_MEMORY;
      }
    }
  }
  return 0;
}    


/***********************************************************************************/
/* Function: ASCII_calloc_double_3D_arylen_restricted                     @30_30i@ */
/* Description: Allocate memory for a three-dimensional array of double. Only      */
/*   a restricted range (marked by columns_lower and columns_upper) is actually    */
/*   allocated; the function is very special - maybe not of much general use.      */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                   @i30_30@ */
/***********************************************************************************/

int ASCII_calloc_double_3D_arylen_restricted (double ****value,
					      int rows,
					      int columns,
					      int columns_lower, int columns_upper,
					      int *length)
{
  int i=0, j=0;

  if ( (*value = (double ***) calloc (rows, sizeof (double **))) == NULL )  
    return ASCII_NO_MEMORY;

    
  for (i=0; i<rows; i++) {
    if ( ((*value)[i] = (double **) calloc (columns, sizeof (double *))) == NULL )
      return ASCII_NO_MEMORY;

    for (j=columns_lower; j<=columns_upper; j++)
      if (length[j]>0)  /* need to check because some OSs don't like calloc (0) */
	if (((*value)[i][j] = (double *) calloc (length[j], sizeof (double))) == NULL) 
	  return ASCII_NO_MEMORY;
  }

  return 0;
}

/***********************************************************************************/
/* Function: ASCII_calloc_double_5D_arylen_restricted                     @30_30i@ */
/* Description: Allocate memory for a five-dimensional array of double. Only      */
/*   a restricted range (marked by columns_lower and columns_upper) is actually    */
/*   allocated; the function is very special - maybe not of much general use.      */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Claudia Emde (extended 3D version)                             @i30_30@ */
/***********************************************************************************/

int ASCII_calloc_double_5D_arylen_restricted (double ******value,
					      int dim1,
					      int dim2,
					      int dim3, 
					      int dim4,
					      int dim4_lower, int dim4_upper,
					      int *dim5)
{
  int i=0, j=0, k=0, l=0;

  if ( (*value = (double *****) calloc (dim1, sizeof (double ****))) == NULL )  
    return ASCII_NO_MEMORY;
    
  for (i=0; i<dim1; i++) {
    if ( ((*value)[i] = (double ****) calloc (dim2, sizeof (double ***))) == NULL )
      return ASCII_NO_MEMORY;
    
    for (j=0; j<dim2; j++) {
      if ( ((*value)[i][j] = (double ***) calloc (dim3, sizeof (double **))) == NULL )
	return ASCII_NO_MEMORY;
      
      for (k=0; k<dim3; k++) {
	if ( ((*value)[i][j][k] = (double **) calloc (dim4, sizeof (double *))) == NULL )
	  return ASCII_NO_MEMORY;
	
	for (l=dim4_lower; l<=dim4_upper; l++)
	  if (dim5[l]>0)  /* need to check because some OSs don't like calloc (0) */
	    if (((*value)[i][j][k][l] = (double *) calloc (dim5[l], sizeof (double))) == NULL) 
	      return ASCII_NO_MEMORY;
      }
    }
  }
  
  return 0;
}



/***********************************************************************************/
/* Function: ASCII_calloc_float_3D                                        @30_30i@ */
/* Description: Allocate memory for a three-dimensional array of float.            */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                   @i30_30@ */
/***********************************************************************************/

int ASCII_calloc_float_3D  (float ****value,
			 int rows,
			 int columns,
			 int length)
{
  int i=0, j=0;

  if ( (*value = (float ***) calloc (rows, sizeof (float **))) == NULL )  
    return ASCII_NO_MEMORY;

    
  for (i=0; i<rows; i++) {
    if ( ((*value)[i] = (float **) calloc (columns, sizeof (float *))) == NULL )
      return ASCII_NO_MEMORY;
    for (j=0; j<columns; j++)
      if (((*value)[i][j] = (float *) calloc (length, sizeof (float))) == NULL) 
	return ASCII_NO_MEMORY;
    
  }

  return 0;
}


/***********************************************************************************/
/* Function: ASCII_calloc_float_4D                                        @30_30i@ */
/* Description: Allocate memory for a four-dimensional array of float.             */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                   @i30_30@ */
/***********************************************************************************/

int ASCII_calloc_float_4D  (float *****value,
			 int rows,
			 int columns,
			 int length,
			 int fourth_dimension)
{
  int i=0, j=0, k=0;

  if ( (*value = (float ****) calloc (rows, sizeof (float ***))) == NULL )  
    return ASCII_NO_MEMORY;

    
  for (i=0; i<rows; i++) {
    if ( ((*value)[i] = (float ***) calloc (columns, sizeof (float **))) == NULL )
      return ASCII_NO_MEMORY;
    for (j=0; j<columns; j++) {
      if (((*value)[i][j] = (float **) calloc (length, sizeof (float *))) == NULL) 
	return ASCII_NO_MEMORY;
      for (k=0; k<length; k++)
	if (((*value)[i][j][k] = (float *) calloc (fourth_dimension, sizeof (float))) == NULL) 
	  return ASCII_NO_MEMORY;
    }    
  }

  return 0;
}


/***********************************************************************************/
/* Function: ASCII_calloc_double_4D                                       @30_30i@ */
/* Description: Allocate memory for a four-dimensional array of double.            */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                   @i30_30@ */
/***********************************************************************************/

int ASCII_calloc_double_4D  (double *****value,
			     int rows,
			     int columns,
			     int length,
			     int fourth_dimension)
{
  int i=0, j=0, k=0;

  if ( (*value = (double ****) calloc (rows, sizeof (double ***))) == NULL )  
    return ASCII_NO_MEMORY;

    
  for (i=0; i<rows; i++) {
    if ( ((*value)[i] = (double ***) calloc (columns, sizeof (double **))) == NULL )
      return ASCII_NO_MEMORY;
    for (j=0; j<columns; j++) {
      if (((*value)[i][j] = (double **) calloc (length, sizeof (double *))) == NULL) 
	return ASCII_NO_MEMORY;
      for (k=0; k<length; k++)
	if (((*value)[i][j][k] = (double *) calloc (fourth_dimension, sizeof (double))) == NULL) 
	  return ASCII_NO_MEMORY;
    }    
  }

  return 0;
}


/***********************************************************************************/
/* Function: ASCII_calloc_float_5D                                        @30_30i@ */
/* Description: Allocate memory for a five-dimensional array of float.             */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                   @i30_30@ */
/***********************************************************************************/

int ASCII_calloc_float_5D  (float ******value,
			    int rows,
			    int columns,
			    int length,
			    int fourth_dimension,
			    int fifth_dimension)
{
  int i=0, j=0, k=0, l=0;

  if ( (*value = (float *****) calloc (rows, sizeof (float ****))) == NULL )
    if (errno==ENOMEM)
      return ASCII_NO_MEMORY;

  for (i=0; i<rows; i++) {
    if ( ((*value)[i] = (float ****) calloc (columns, sizeof (float ***))) == NULL )
      if (errno==ENOMEM)
	return ASCII_NO_MEMORY;
    for (j=0; j<columns; j++) {
      if (((*value)[i][j] = (float ***) calloc (length, sizeof (float **))) == NULL)
	if (errno==ENOMEM)
	  return ASCII_NO_MEMORY;
      for (k=0; k<length; k++) {
	if (((*value)[i][j][k] = (float **) calloc (fourth_dimension, sizeof (float *))) == NULL)
	  if (errno==ENOMEM)
	    return ASCII_NO_MEMORY;
	for (l=0; l<fourth_dimension; l++) 
	  if (((*value)[i][j][k][l] = (float *) calloc (fifth_dimension, sizeof (float))) == NULL)
	    if (errno==ENOMEM)
	      return ASCII_NO_MEMORY;
      }
    }    
  }
  
  return 0;
}


/***********************************************************************************/
/* Function: ASCII_calloc_double_5D                                        @30_30i@ */
/* Description: Allocate memory for a five-dimensional array of float.             */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Claudia Emde , adapted from Arve Kylling                                @i30_30@ */
/***********************************************************************************/
int ASCII_calloc_double_5D  (double ******value,
			    int rows,
			    int columns,
			    int length,
			    int fourth_dimension,
			    int fifth_dimension)
{
  int i=0, j=0, k=0, l=0;

  if ( (*value = (double *****) calloc (rows, sizeof (double ****))) == NULL )
    if (errno==ENOMEM)
      return ASCII_NO_MEMORY;

  for (i=0; i<rows; i++) {
    if ( ((*value)[i] = (double ****) calloc (columns, sizeof (double ***))) == NULL )
      if (errno==ENOMEM)
	return ASCII_NO_MEMORY;
    for (j=0; j<columns; j++) {
      if (((*value)[i][j] = (double ***) calloc (length, sizeof (double **))) == NULL)
	if (errno==ENOMEM)
	  return ASCII_NO_MEMORY;
      for (k=0; k<length; k++) {
	if (((*value)[i][j][k] = (double **) calloc (fourth_dimension, sizeof (double *))) == NULL)
	  if (errno==ENOMEM)
	    return ASCII_NO_MEMORY;
	for (l=0; l<fourth_dimension; l++) 
	  if (((*value)[i][j][k][l] = (double *) calloc (fifth_dimension, sizeof (double))) == NULL)
	    if (errno==ENOMEM)
	      return ASCII_NO_MEMORY;
      }
    }    
  }
  
  return 0;
}    


/***********************************************************************************/
/* Function: ASCII_calloc_float_6D                                        @30_30i@ */
/* Description: Allocate memory for a six-dimensional array of float.              */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling (extended to 6D by Claudia Emde                   @i30_30@ */
/***********************************************************************************/

int ASCII_calloc_float_6D  (float *******value,
			    int rows,
			    int columns,
			    int length,
			    int fourth_dimension,
			    int fifth_dimension, 
                            int sixth_dimension)
{
  int i=0, j=0, k=0, l=0, m=0;

  if ( (*value = (float ******) calloc (rows, sizeof (float *****))) == NULL )
    if (errno==ENOMEM)
      return ASCII_NO_MEMORY;

  for (i=0; i<rows; i++) {
    if ( ((*value)[i] = (float *****) calloc (columns, sizeof (float ****))) == NULL )
      if (errno==ENOMEM)
	return ASCII_NO_MEMORY;
    for (j=0; j<columns; j++) {
      if (((*value)[i][j] = (float ****) calloc (length, sizeof (float ***))) == NULL)
	if (errno==ENOMEM)
	  return ASCII_NO_MEMORY;
      for (k=0; k<length; k++) {
	if (((*value)[i][j][k] = (float ***) calloc (fourth_dimension, sizeof (float **))) == NULL)
	  if (errno==ENOMEM)
	    return ASCII_NO_MEMORY;
	for (l=0; l<fourth_dimension; l++){
	  if (((*value)[i][j][k][l] = (float **) calloc (fifth_dimension, sizeof (float *))) == NULL)
	    if (errno==ENOMEM)
	      return ASCII_NO_MEMORY;
          for (m=0; m<fifth_dimension; m++){
            if (((*value)[i][j][k][l][m] = (float *) calloc (sixth_dimension, sizeof (float))) == NULL)
              if (errno==ENOMEM)
                return ASCII_NO_MEMORY;  
          }
        }
      }    
    }
  }
  
  return 0;
}    


/***********************************************************************************/
/* Function: ASCII_calloc_float                                           @30_30i@ */
/* Description: Allocate memory for a two-dimensional array of float.              */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_calloc_float  (float ***value,
			 int rows,
			 int columns)
{
  int i=0;

  if (rows>0) {
    if ( (*value = (float **) calloc (rows, sizeof (float *))) == NULL )  
      return ASCII_NO_MEMORY;
    
    if (columns>0) 
      for (i=0; i<rows; i++) 
	if ( ((*value)[i] = (float *) calloc (columns, sizeof (float))) == NULL )
	  return ASCII_NO_MEMORY;
  }

  return 0;
}  

/***********************************************************************************/
/* Function: ASCII_calloc_short                                           @30_30i@ */
/* Description: Allocate memory for a two-dimensional array of short.              */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_calloc_short  (short ***value,
			 int rows,
			 int columns)
{
  int i=0;

  if (rows>0) {
    if ( (*value = (short **) calloc (rows, sizeof (short *))) == NULL )  
      return ASCII_NO_MEMORY;
    
    if (columns>0) 
      for (i=0; i<rows; i++) 
	if ( ((*value)[i] = (short *) calloc (columns, sizeof (short))) == NULL )
	  return ASCII_NO_MEMORY;
  }

  return 0;
}    
  



/***********************************************************************************/
/* Function: ASCII_free_int                                               @30_30i@ */
/* Description: Free memory, which has been allocated with ASCII_calloc_int.       */
/* Parameters:                                                                     */
/*  int **value:  Two-dimensional array of int                                     */ 
/*  int rows:        Number of rows, specified by caller                           */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_free_int (int **value, 
		       int rows)
{
  int i=0;

  for (i=0; i<rows; i++) 
    free (value[i]);
  free (value);

  return 0;
}    


/***********************************************************************************/
/* Function: ASCII_free_int_3D                                            @30_30i@ */
/* Description: Free memory, which has been allocated with ASCII_calloc_int_3D     */
/* Parameters:                                                                     */
/*  int ***value:  Three-dimensional array of int                                  */ 
/*  int rows:        Number of rows, specified by caller                           */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_free_int_3D (int ***value, 
                       int rows,
                       int columns)
{
  int i=0, j=0;

  for (i=0; i<rows; i++){
    for (j=0; j<columns; j++)
      free (value[i][j]);
    free (value[i]);
  }
  free (value);
  return 0;
}


/***********************************************************************************/
/* Function: ASCII_free_int_4D                                            @30_30i@ */
/* Description: Free memory, which has been allocated with ASCII_calloc_double_4D. */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                   @i30_30@ */
/***********************************************************************************/

int ASCII_free_int_4D (int ****value, 
		       int rows,
		       int columns,
		       int length)
{
  int i=0, j=0, k=0;

  for (i=0; i<rows; i++) {
    for (j=0; j<columns; j++) {
      for (k=0; k<length; k++)
	free (value[i][j][k]);
      free (value[i][j]);
    }
    free (value[i]);
  }

  free (value);
  return 0;
}


/***********************************************************************************/
/* Function: ASCII_free_double                                            @30_30i@ */
/* Description: Free memory, which has been allocated with ASCII_calloc_double.    */
/* Parameters:                                                                     */
/*  double **value:  Two-dimensional array of double                               */ 
/*  int rows:        Number of rows, specified by caller                           */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_free_double (double **value, 
		       int rows)
{
  int i=0;

  for (i=0; i<rows; i++) 
    free (value[i]);
  free (value);

  return 0;
}    



/***********************************************************************************/
/* Function: ASCII_free_float                                             @30_30i@ */
/* Description: Free memory, which has been allocated with ASCII_calloc_float.     */
/* Parameters:                                                                     */
/*  float **value:   Two-dimensional array of float                                */ 
/*  int rows:        Number of rows, specified by caller                           */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_free_float (float **value, 
		      int rows)
{
  int i=0;

  for (i=0; i<rows; i++) 
    free (value[i]);

  free (value);

  return 0;
}    

/***********************************************************************************/
/* Function: ASCII_free_short                                             @30_30i@ */
/* Description: Free memory, which has been allocated with ASCII_calloc_short.     */
/* Parameters:                                                                     */
/*  float **value:   Two-dimensional array of short                                */ 
/*  int rows:        Number of rows, specified by caller                           */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_free_short (short **value, 
		      int rows)
{
  int i=0;

  for (i=0; i<rows; i++) 
    free (value[i]);

  free (value);

  return 0;
}


/***********************************************************************************/
/* Function: ASCII_free_double_3D                                         @30_30i@ */
/* Description: Free memory, which has been allocated with ASCII_calloc_double_3D. */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                   @i30_30@ */
/***********************************************************************************/

int ASCII_free_double_3D (double ***value, 
			  int rows,
			  int columns)
{
  int i=0, j=0;

  for (i=0; i<rows; i++) {
    for (j=0; j<columns; j++)
      free (value[i][j]);
    free (value[i]);
  }

  free (value);
  return 0;
}    

/***********************************************************************************/
/* Function: ASCII_free_float_3D                                          @30_30i@ */
/* Description: Free memory, which has been allocated with ASCII_calloc_float_3D.  */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                   @i30_30@ */
/***********************************************************************************/

int ASCII_free_float_3D (float ***value, 
			 int rows,
			 int columns)
{
  int i=0, j=0;

  for (i=0; i<rows; i++) {
    for (j=0; j<columns; j++)
      free (value[i][j]);
    free (value[i]);
  }

  free (value);
  return 0;
}    


/***********************************************************************************/
/* Function: ASCII_free_float_4D                                          @30_30i@ */
/* Description: Free memory, which has been allocated with ASCII_calloc_float_4D.  */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                   @i30_30@ */
/***********************************************************************************/

int ASCII_free_float_4D (float ****value, 
			 int rows,
			 int columns,
			 int length)
{
  int i=0, j=0, k=0;

  for (i=0; i<rows; i++) {
    for (j=0; j<columns; j++) {
      for (k=0; k<length; k++)
	free (value[i][j][k]);
      free (value[i][j]);
    }
    free (value[i]);
  }

  free (value);
  return 0;
}    

/***********************************************************************************/
/* Function: ASCII_free_double_4D                                         @30_30i@ */
/* Description: Free memory, which has been allocated with ASCII_calloc_double_4D. */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                   @i30_30@ */
/***********************************************************************************/

int ASCII_free_double_4D (double ****value, 
			  int rows,
			  int columns,
			  int length)
{
  int i=0, j=0, k=0;

  for (i=0; i<rows; i++) {
    for (j=0; j<columns; j++) {
      for (k=0; k<length; k++)
	free (value[i][j][k]);
      free (value[i][j]);
    }
    free (value[i]);
  }

  free (value);
  return 0;
}    

/***********************************************************************************/
/* Function: ASCII_free_float_5D                                          @30_30i@ */
/* Description: Free memory, which has been allocated with ASCII_calloc_float_5D.  */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Arve Kylling                                                   @i30_30@ */
/***********************************************************************************/

int ASCII_free_float_5D (float *****value, 
			 int rows,
			 int columns,
			 int length,
			 int fourth_dimension)
{
  int i=0, j=0, k=0, l=0;

  for (i=0; i<rows; i++) {
    for (j=0; j<columns; j++) {
      for (k=0; k<length; k++) {
	for (l=0; l<fourth_dimension; l++) 
	  free (value[i][j][k][l]);
	free (value[i][j][k]);
      }
      free (value[i][j]);
    }
    free (value[i]);
  }

  free (value);
  return 0;
}    

/***********************************************************************************/
/* Function: ASCII_free_double_5D                                          @30_30i@ */
/* Description: Free memory, which has been allocated with ASCII_calloc_double_5D.  */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Claudia Emde                                                   @i30_30@ */
/***********************************************************************************/

int ASCII_free_double_5D (double *****value, 
			 int rows,
			 int columns,
			 int length,
			 int fourth_dimension)
{
  int i=0, j=0, k=0, l=0;

  for (i=0; i<rows; i++) {
    for (j=0; j<columns; j++) {
      for (k=0; k<length; k++) {
	for (l=0; l<fourth_dimension; l++) 
	  free (value[i][j][k][l]);
	free (value[i][j][k]);
      }
      free (value[i][j]);
    }
    free (value[i]);
  }

  free (value);
  return 0;
}    

/***********************************************************************************/
/* Function: ASCII_free_float_6D                                          @30_30i@ */
/* Description: Free memory, which has been allocated with ASCII_calloc_float_6D.  */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Ulrich Hamann (modified ASCII_free_float_5D, Arve Kylling)     @i30_30@ */
/***********************************************************************************/

int ASCII_free_float_6D (float ******value, 
			 int rows,
			 int columns,
			 int length,
			 int fourth_dimension,
                         int fifth_dimension)
{
  int i=0, j=0, k=0, l=0, m=0;

  for (i=0; i<rows; i++) {
    for (j=0; j<columns; j++) {
      for (k=0; k<length; k++) {
	for (l=0; l<fourth_dimension; l++) {
 	  for (m=0; m<fifth_dimension; m++)
	    free (value[i][j][k][l][m]);
	  free (value[i][j][k][l]);
	}
	free (value[i][j][k]);
      }
      free (value[i][j]);
    }
    free (value[i]);
  }

  free (value);
  return 0;
}    

/***********************************************************************************/
/* Function: ASCII_free_float_7D                                          @30_30i@ */
/* Description: Free memory, which has been allocated with ASCII_calloc_float_7D.  */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Ulrich Hamann (modified ASCII_free_float_5D, Arve Kylling)     @i30_30@ */
/***********************************************************************************/

int ASCII_free_float_7D (float *******value, 
			 int rows,
			 int columns,
			 int length,
			 int fourth_dimension,
                         int fifth_dimension,
                         int sixth_dimension)
{
  int i=0, j=0, k=0, l=0, m=0, n=0;

  for (i=0; i<rows; i++) {
    for (j=0; j<columns; j++) {
      for (k=0; k<length; k++) {
	for (l=0; l<fourth_dimension; l++) {
	  for (m=0; m<fifth_dimension; m++) {
	    for (n=0; n<sixth_dimension; n++)
		free (value[i][j][k][l][m][n]);
            free (value[i][j][k][l][m]);
	  }
	  free (value[i][j][k][l]);
	}
	free (value[i][j][k]);
      }
      free (value[i][j]);
    }
    free (value[i]);
  }

  free (value);
  return 0;
}    





/***********************************************************************************/
/* Function: ASCII_readfile_nrows                                         @30_30i@ */
/* Description: Read an ASCII file into a two-dimensional array of strings;        */
/*		before calling ASCII_readfile, the file must be parsed with        */
/*		ASCII_checkfile, and  memory must be allocated with                */
/*              ASCII_calloc_string. If parameter nrows is larger or equal 0,      */
/*              only the first non-empty, non-comment nrows are read and the       */
/*              rest of the file is ignored.                                       */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_readfile_nrows (char *filename, int nrows,
			  char ***array)
{
  FILE *f=NULL;
  char line[MAX_LENGTH_OF_LINE+1]="";
  char *string=NULL;
  char *t=NULL;
  int row=0, column=0;

  
  string = line;

  
  if ( (f = fopen(filename, "r")) == NULL)  
    return ASCIIFILE_NOT_FOUND;
  

  while ( fgets (string, MAX_LENGTH_OF_LINE, f) != NULL )  {

    column=0;

    if ( (t = strtok (string, " \t\n") ) != NULL)  {  /* if not an empty line */ 
      if (!ASCII_comment(t[0]))  { /* if not a comment */
	strcpy (array[row][column++], t);
	
	while ( (t = strtok (NULL, " \t\n") ) != NULL)  {
	  if (ASCII_comment(t[0]))     /* if comment */
	    break;
	  strcpy (array[row][column++], t);
	}

	row++;
      }
    }

    string = line;

    /* stop after reading nrows rows */
    if (nrows>=0 && row>=nrows)
      break;
  }


  fclose (f);
  return 0;
}


/***********************************************************************************/
/* Function: ASCII_readfile                                               @30_30i@ */
/* Description: Read an ASCII file into a two-dimensional array of strings;        */
/*		before calling ASCII_readfile, the file must be parsed with        */
/*		ASCII_checkfile, and  memory must be allocated with                */
/*              ASCII_calloc_string.                                               */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_readfile (char *filename, 
		    char ***array)
{
  int status=0;

  status = ASCII_readfile_nrows (filename, -1, array);

  return status;
}



/***********************************************************************************/
/* Function: ASCII_string2double                                          @30_30i@ */
/* Description: Convert a two-dimensional array of strings to a two-dimensional    */
/*		array of double.                                                   */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_string2double (double **value, 
			 char *** string,
			 int rows,
			 int columns)
{
  int i=0, j=0;
  char *dummy=NULL;

  for (i=0; i<rows; i++) 
    for (j=0; j<columns; j++)  {

      if (string[i][j][0] == 0)  
	value[i][j] = NAN;
      else  
	value[i][j] = strtod (string[i][j], &dummy);

    }

  return 0;
}    


/***********************************************************************************/
/* Function: ASCII_string2int                                             @30_30i@ */
/* Description: Convert a two-dimensional array of strings to a two-dimensional    */
/*		array of int.                                                      */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_string2int (int **value, 
		      char *** string,
		      int rows,
		      int columns)
{
  int i=0, j=0;

  for (i=0; i<rows; i++) 
    for (j=0; j<columns; j++)  {

      if (string[i][j][0] == 0)  
	value[i][j] = -999;
      else  
	value[i][j] = atoi (string[i][j]);

    }

  return 0;
}    


/***********************************************************************************/
/* Function: ASCII_string2float                                           @30_30i@ */
/* Description: Convert a two-dimensional array of strings to a two-dimensional    */
/*		array of float.                                                    */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_string2float (float **value, 
			char *** string,
			int rows,
			int columns)
{
  int i=0, j=0;
  char *dummy=NULL;

  for (i=0; i<rows; i++)  { 
    for (j=0; j<columns; j++)  {
      if (string[i][j][0] == 0)   {
	value[i][j] = (float) NAN;
      }
      else  {
	value[i][j] = (float) strtod (string[i][j], &dummy);
      }
    }
  }
  return 0;
}    



/************************************************************************************/
/* Function: ASCII_file2string                                             @30_30i@ */
/* Description: Parse an ASCII file and store data in a twodimensional array        */
/*        value[row][column]; memory allocation for value is done automatically.    */
/*        rows is the number of (not empty) rows of the file, max_columns is the    */
/*        maximal number of columns of the file and min_columns is the minimal      */
/*        number of columns of all (not empty) lines; the dimension of the array    */
/*        value is rows * max_columns; rows with less than                          */
/*        max_columns columns are filled up with ""; the allocated memory           */
/*        can be freed with ASCII_free_string (string, rows, max_columns)           */
/* Parameters:                                                                      */
/*  char *filename:    Name of the file which should be parsed                      */
/*  int  *rows:        Number of rows, set by function                              */
/*  int  *min_columns: Minimum number of columns, set by function                   */
/*  int  *max_columns: Maximum number of columns, set by function                   */
/*  char ***string:    Pointer to a two-dimensional array of strings,               */
/*                     string [0 ... rows-1][0 ... max_columns-1][0...n_character-1]*/
/*                     Memory is allocated automatically.                           */
/*                                                                                  */
/* Return value:                                                                    */
/*  0  if o.k., <0 if error                                                         */
/*                                                                                  */
/* Example:                                                                         */
/* Files:                                                                           */
/* Known bugs:                                                                      */
/* Author: Ulrich Hamann, 17.12.2007                                       @i30_30@ */
/************************************************************************************/

int ASCII_file2string (char *filename,   
		       int *rows,        
		       int *max_columns, 
		       int *min_columns, 
		       char ****string)
{
  int status=0;
  int max_length=0;

  /* count rows and columns of ASCII file <filename> */
  if ( (status = ASCII_checkfile (filename, 
				  rows, 
				  min_columns, 
				  max_columns, 
				  &max_length)) != 0)
    return status;
  


  /* allocate memory for string array */
  if ( (status = ASCII_calloc_string (string, 
				      *rows, 
				      *max_columns, 
				      max_length)) != 0 )
    return status;
  
  /* read ASCII file to string array */
  if ( (status = ASCII_readfile (filename, (*string))) != 0)
    return status;

  return 0;  /* everything ok */
} 



/***********************************************************************************/
/* Function: ASCII_file2double_nrows                                      @30_30i@ */
/* Description: Parse an ASCII file and store data in a twodimensional array       */
/*        value[row][column]; memory allocation for value is done automatically.   */
/*        rows is the number of (not empty) rows of the file, max_columns is the   */
/*        maximal number of columns of the file and min_columns is the minimal     */
/*        number of columns of all (not empty) lines; the dimension of the array   */
/*        value is rows * max_columns; strings that cannot be interpreted as       */
/*        floating point number are converted to 0; rows with less than            */
/*        max_columns columns are filled up with NAN; the allocated memory         */
/*        can be freed with ASCII_free_double (value, rows).                       */
/*        If parameter nrows is larger or equal 0, only                            */
/*        the first non-empty, non-comment nrows are read and the                  */
/*        rest of the file is ignored.                                             */
/* Parameters:                                                                     */
/*  char *filename:     Name of the file which should be parsed                    */
/*  int  *rows:         Number of rows, set by function                            */
/*  int  *min_columns:  Minimum number of columns, set by function                 */
/*  int  *max_columns:  Maximum number of columns, set by function                 */
/*  double ***value:    Pointer to a two-dimensional array of double,              */
/*                      value [0 ... rows-1][0 ... max_columns-1].                 */
/*                      Memory is allocated automatically.                         */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_file2double_nrows (char *filename, int nrows,
			     int *rows,        
			     int *max_columns, 
			     int *min_columns, 
			     double ***value)  
{
  char ***string = NULL;
  int status=0;
  int max_length=0;

  /* count rows and columns of ASCII file <filename> */
  if ( (status = ASCII_checkfile_nrows (filename, nrows,
					rows, 
					min_columns, 
					max_columns, 
					&max_length)) != 0)
    return status;
  


  /* allocate memory for string array */
  if ( (status = ASCII_calloc_string (&string, 
				      *rows, 
				      *max_columns, 
				      max_length)) != 0 )
    return status;
  
  /* read ASCII file to string array */
  if ( (status = ASCII_readfile_nrows (filename, nrows, string)) != 0)
    return status;

  /* allocate memory for double array */
  if ( (status = ASCII_calloc_double (value, *rows, *max_columns)) != 0 )
    return status;

  /* convert string array to double array */
  if ( (status = ASCII_string2double (*value, string, *rows, *max_columns)) != 0)
    return status; 



  /* free memory of string array */
  ASCII_free_string (string, *rows, *max_columns);

  return 0;  /* everything ok */
} 


/***********************************************************************************/
/* Function: ASCII_file2double                                            @30_30i@ */
/* Description: Parse an ASCII file and store data in a twodimensional array       */
/*        value[row][column]; memory allocation for value is done automatically.   */
/*        rows is the number of (not empty) rows of the file, max_columns is the   */
/*        maximal number of columns of the file and min_columns is the minimal     */
/*        number of columns of all (not empty) lines; the dimension of the array   */
/*        value is rows * max_columns; strings that cannot be interpreted as       */
/*        floating point number are converted to 0; rows with less than            */
/*        max_columns columns are filled up with NAN; the allocated memory         */
/*        can be freed with ASCII_free_double (value, rows).                       */
/* Parameters:                                                                     */
/*  char *filename:     Name of the file which should be parsed                    */
/*  int  *rows:         Number of rows, set by function                            */
/*  int  *min_columns:  Minimum number of columns, set by function                 */
/*  int  *max_columns:  Maximum number of columns, set by function                 */
/*  double ***value:    Pointer to a two-dimensional array of double,              */
/*                      value [0 ... rows-1][0 ... max_columns-1].                 */
/*                      Memory is allocated automatically.                         */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_file2double (char *filename,   
		       int *rows,        
		       int *max_columns, 
		       int *min_columns, 
		       double ***value)  
{
  int status=0;
  
  status = ASCII_file2double_nrows (filename, -1,
				    rows,        
				    max_columns, 
				    min_columns, 
				    value);
  
  return status;
} 



/***********************************************************************************/
/* Function: ASCII_file2int                                               @30_30i@ */
/* Description: Parse an ASCII file and store data in a twodimensional array       */
/*        value[row][column]; memory allocation for value is done automatically.   */
/*        rows is the number of (not empty) rows of the file, max_columns is the   */
/*        maximal number of columns of the file and min_columns is the minimal     */
/*        number of columns of all (not empty) lines; the dimension of the array   */
/*        value is rows * max_columns; strings that cannot be interpreted as       */
/*        floating point number are converted to 0; rows with less than            */
/*        max_columns columns are filled up with NAN; the allocated memory         */
/*        can be freed with ASCII_free_int (value, rows).                          */
/* Parameters:                                                                     */
/*  char *filename:     Name of the file which should be parsed                    */
/*  int  *rows:         Number of rows, set by function                            */
/*  int  *min_columns:  Minimum number of columns, set by function                 */
/*  int  *max_columns:  Maximum number of columns, set by function                 */
/*  int  ***value:      Pointer to a two-dimensional array of int,                 */
/*                      value [0 ... rows-1][0 ... max_columns-1].                 */
/*                      Memory is allocated automatically.                         */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_file2int (char *filename,
		    int *rows,
		    int *max_columns,
		    int *min_columns,
		    int ***value)
{
  char ***string = NULL;
  int status=0;
  int max_length=0;

  /* count rows and columns of ASCII file <filename> */
  if ( (status = ASCII_checkfile (filename, 
				  rows, 
				  min_columns, 
				  max_columns, 
				  &max_length)) != 0)
    return status;
  


  /* allocate memory for string array */
  if ( (status = ASCII_calloc_string (&string, 
				      *rows, 
				      *max_columns, 
				      max_length)) != 0 )
    return status;
  
  /* read ASCII file to string array */
  if ( (status = ASCII_readfile (filename, string)) != 0)
    return status;

  /* allocate memory for double array */
  if ( (status = ASCII_calloc_int (value, *rows, *max_columns)) != 0 )
    return status;

  /* convert string array to double array */
  if ( (status = ASCII_string2int (*value, string, *rows, *max_columns)) != 0)
    return status; 



  /* free memory of string array */
  ASCII_free_string (string, *rows, *max_columns);

  return 0;  /* everything ok */
} 



/***********************************************************************************/
/* Function: ASCII_file2float                                             @30_30i@ */
/* Description: Read an ASCII file and store data in a twodimensional array        */
/*        value[row][column]; memory allocation for value is done automatically.   */
/*        rows is the number of (not empty) rows of the file, max_columns is the   */
/*        maximal number of columns of the file and min_columns is the minimal     */
/*        number of columns of all (not empty) lines; the dimension of the array   */
/*        value is rows * max_columns; strings that cannot be interpreted as       */
/*        floating point number are converted to 0; rows with less than            */
/*        max_columns columns are filled up with NAN; the allocated memory         */
/*        can be freed with ASCII_free_float (value, rows).                        */
/*                                                                                 */
/* Parameters:                                                                     */
/*  char *filename:     Name of the file which should be parsed                    */
/*  int  *rows:         Number of rows, set by function                            */
/*  int  *min_columns:  Minimum number of columns, set by function                 */
/*  int  *max_columns:  Maximum number of columns, set by function                 */
/*  float ***value:     Pointer to a two-dimensional array of float,               */
/*                      value [0 ... rows-1][0 ... max_columns-1].                 */
/*                      Memory is allocated automatically.                         */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_file2float (char *filename,
		      int *rows,
		      int *max_columns,
		      int *min_columns,
		      float ***value)
{
  char ***string = NULL;
  int status=0;
  int max_length=0;

  /* count rows and columns of ASCII file <filename> */
  if ( (status = ASCII_checkfile (filename, 
				  rows, 
				  min_columns, 
				  max_columns, 
				  &max_length)) != 0)
    return status;
  

  /* allocate memory for string array */
  if ( (status = ASCII_calloc_string (&string, 
				      *rows, 
				      *max_columns, 
				      max_length)) != 0 )
    return status;
  

  
  /* read ASCII file to string array */
  if ( (status = ASCII_readfile (filename, string)) != 0)
     return status;


  /* allocate memory for float array */
  if ( (status = ASCII_calloc_float (value, *rows, *max_columns)) != 0 )
    return status;

  /* convert string array to float array */
  if ( (status = ASCII_string2float (*value, string, *rows, *max_columns)) != 0)
    return status; 

  /* free memory of string array */
  ASCII_free_string (string, *rows, *max_columns);

  
  return 0;  /* everything ok */
} 



/***********************************************************************************/
/* Function: ASCII_column                                                 @30_30i@ */
/* Description: Extract a specified column from a two-dimensional array of double. */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  Pointer to the column.                                                         */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

double *ASCII_column (double **value, int rows, int column)  
{
  int i=0;
  double *col = (double *) calloc (rows, sizeof(double));

  for (i=0; i<rows; i++)
    col[i] = value[i][column];

  return col;
}
  


/***********************************************************************************/
/* Function: ASCII_column_float                                           @30_30i@ */
/* Description: Extract a specified column from a two-dimensional array of float.  */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  Pointer to the column.                                                         */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

float *ASCII_column_float (float **value, int rows, int column)  
{
  int i=0;
  float *col = (float *) calloc (rows, sizeof(float));

  for (i=0; i<rows; i++)
    col[i] = value[i][column];

  return col;
}
  

/***********************************************************************************/
/* Function: ASCII_row                                                    @30_30i@ */
/* Description: Extract a specified row from a two-dimensional array of double.    */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  Pointer to the row.                                                            */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

double *ASCII_row (double **value, int columns, int row)  
{
  int i=0;
  double *col = (double *) calloc (columns, sizeof(double));

  for (i=0; i<columns; i++)
    col[i] = value[row][i];

  return col;
}
  

/***********************************************************************************/
/* Function: read_1c_file                                                 @30_30i@ */
/* Description: Read an ASCII file with (at least) 1 column.                       */
/*              Only the first column is returned in array first;                  */
/*              n is the number of values returned.                                */
/*              Memory allocation for first is done automatically;                 */
/*              field can be freed with a simple free().                           */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int read_1c_file (char *filename, 
		  double **first, int *n)
{
  int max_columns=0, min_columns=0;
  double **data=NULL;
  int status=0;

  /* read file */
  if ( (status = ASCII_file2double (filename, n, 
				    &max_columns, &min_columns, &data)) != 0)
    return status;


  /* check, if at least two columns */
  if (min_columns < 1)   {
    ASCII_free_double (data, *n);
    return LESS_THAN_TWO_COLUMNS;
  }
  
  *first  = ASCII_column (data, *n, 0);
  
  ASCII_free_double (data, *n);

  return 0;
}




/***********************************************************************************/
/* Function: read_1c_file_float                                           @30_30i@ */
/* Description: Read an ASCII file with (at least) 1 column.                       */
/*              Only the first column is returned in array first;                  */
/*              n is the number of values returned.                                */
/*              Memory allocation for first is done automatically;                 */
/*              field can be freed with a simple free().                           */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int read_1c_file_float (char *filename, 
			float **first, int *n)
{
  int max_columns=0, min_columns=0;
  float **data=NULL;
  int status=0;

  /* read file */
  if ( (status = ASCII_file2float (filename, n, 
				   &max_columns, &min_columns, &data)) != 0)
    return status;


  /* check, if at least two columns */
  if (min_columns < 1)   {
    ASCII_free_float (data, *n);
    return LESS_THAN_TWO_COLUMNS;
  }
  
  *first  = ASCII_column_float (data, *n, 0);
  
  ASCII_free_float (data, *n);

  return 0;
}




/***********************************************************************************/
/* Function: read_2c_file                                                 @30_30i@ */
/* Description: Read an ASCII file with (at least) 2 columns.                      */
/*              Only the first two column are returned in arrays first and second. */
/*              n is the number of values returned.                                */
/*              Memory allocation for first and second is done automatically;      */
/*              fields can be freed with a simple free().                          */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int read_2c_file (char *filename, 
		  double **first, double **second, int *n)
{
  int max_columns=0, min_columns=0;
  double **data=NULL;
  int status=0;

  /* read file */
  if ( (status = ASCII_file2double (filename, n, 
				    &max_columns, &min_columns, &data)) != 0)
    return status;


  /* check, if at least two columns */
  if (min_columns < 2)   {
    ASCII_free_double (data, *n);
    return LESS_THAN_TWO_COLUMNS;
  }
  
  *first  = ASCII_column (data, *n, 0);
  *second = ASCII_column (data, *n, 1);
  
  ASCII_free_double (data, *n);

  return 0;
}



/***********************************************************************************/
/* Function: read_2c_file_float                                           @30_30i@ */
/* Description: Read an ASCII file with (at least) 2 columns to a float array.     */
/*              Only the first two column are returned in arrays first and second. */
/*              n is the number of values returned.                                */
/*              Memory allocation for first and second is done automatically;      */
/*              fields can be freed with a simple free().                          */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int read_2c_file_float (char *filename, 
			float **first, float **second, int *n)
{
  int max_columns=0, min_columns=0;
  float **data=NULL;
  int status=0;

  /* read file */
  if ( (status = ASCII_file2float (filename, n, 
				   &max_columns, &min_columns, &data)) != 0)
    return status;


  /* check, if at least two columns */
  if (min_columns < 2)   {
    ASCII_free_float (data, *n);
    return LESS_THAN_TWO_COLUMNS;
  }
  
  *first  = ASCII_column_float (data, *n, 0);
  *second = ASCII_column_float (data, *n, 1);
  
  ASCII_free_float (data, *n);

  return 0;
}



/***********************************************************************************/
/* Function: read_3c_file                                                 @30_30i@ */
/* Description: Read an ASCII file with (at least) 3 columns.                      */
/*              Only the first three columns are returned in arrays first, second, */
/*              and third. n is the number of values returned.                     */
/*              Memory allocation for first, second, and third is done             */
/*              automatically; fields can be freed with a simple free().           */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int read_3c_file (char *filename, 
		  double **first, double **second, double **third, int *n)
{
  int max_columns=0, min_columns=0;
  double **data=NULL;
  int status=0;

  /* read file */
  if ( (status = ASCII_file2double (filename, n, 
				    &max_columns, &min_columns, &data)) != 0)
    return status;


  /* check, if at least three columns */
  if (min_columns < 3)   {
    ASCII_free_double (data, *n);
    return LESS_THAN_THREE_COLUMNS;
  }
  
  *first  = ASCII_column (data, *n, 0);
  *second = ASCII_column (data, *n, 1);
  *third  = ASCII_column (data, *n, 2);
  
  ASCII_free_double (data, *n);

  return 0;
}


/***********************************************************************************/
/* Function: read_3c_file_float                                           @30_30i@ */
/* Description: Read an ASCII file with (at least) 3 columns to a float array.     */
/*              Only the first three columns are returned in arrays first, second, */
/*              and third. n is the number of values returned.                     */
/*              Memory allocation for first, second, and third is done             */
/*              automatically; fields can be freed with a simple free().           */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int read_3c_file_float (char *filename, 
			float **first, float **second, float **third, int *n)
{
  int max_columns=0, min_columns=0;
  float **data=NULL;
  int status=0;

  /* read file */
  if ( (status = ASCII_file2float (filename, n, 
				    &max_columns, &min_columns, &data)) != 0)
    return status;


  /* check, if at least three columns */
  if (min_columns < 3)   {
    ASCII_free_float (data, *n);
    return LESS_THAN_THREE_COLUMNS;
  }
  
  *first  = ASCII_column_float (data, *n, 0);
  *second = ASCII_column_float (data, *n, 1);
  *third  = ASCII_column_float (data, *n, 2);
  
  ASCII_free_float (data, *n);

  return 0;
}


/***********************************************************************************/
/* Function: read_4c_file                                                 @30_30i@ */
/* Description: Read an ASCII file with (at least) 4 columns.                      */
/*              Only the first four columns are returned in arrays first, second,  */
/*              third, and fourth. n is the number of values returned.             */
/*              Memory allocation for first, second, third, and fourth is done     */
/*              automatically; fields can be freed with a simple free().           */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int read_4c_file (char *filename, 
		  double **first, double **second, double **third, double **fourth, int *n)
{
  int max_columns=0, min_columns=0;
  double **data=NULL;
  int status=0;

  /* read file */
  if ( (status = ASCII_file2double (filename, n, 
				    &max_columns, &min_columns, &data)) != 0)
    return status;


  /* check, if at least four columns */
  if (min_columns < 4)   {
    ASCII_free_double (data, *n);
    return LESS_THAN_FOUR_COLUMNS;
  }
  
  *first  = ASCII_column (data, *n, 0);
  *second = ASCII_column (data, *n, 1);
  *third  = ASCII_column (data, *n, 2);
  *fourth = ASCII_column (data, *n, 3);
  
  ASCII_free_double (data, *n);

  return 0;
}



/***********************************************************************************/
/* Function: read_4c_file_float                                           @30_30i@ */
/* Description: Read an ASCII file with (at least) 4 columns to a float array.     */
/*              Only the first four columns are returned in arrays first, second,  */
/*              third, and fourth. n is the number of values returned.             */
/*              Memory allocation for first, second, third, and fourth is done     */
/*              automatically; fields can be freed with a simple free().           */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int read_4c_file_float (char *filename, 
			float **first, float **second, float **third, float **fourth, int *n)
{
  int max_columns=0, min_columns=0;
  float **data=NULL;
  int status=0;

  /* read file */
  if ( (status = ASCII_file2float (filename, n, 
				   &max_columns, &min_columns, &data)) != 0)
    return status;


  /* check, if at least four columns */
  if (min_columns < 4)   {
    ASCII_free_float (data, *n);
    return LESS_THAN_FOUR_COLUMNS;
  }
  
  *first  = ASCII_column_float (data, *n, 0);
  *second = ASCII_column_float (data, *n, 1);
  *third  = ASCII_column_float (data, *n, 2);
  *fourth = ASCII_column_float (data, *n, 3);
  
  ASCII_free_float (data, *n);

  return 0;
}





/***********************************************************************************/
/* Function: read_5c_file                                                 @30_30i@ */
/* Description: Read an ASCII file with (at least) 5 columns.                      */
/*              Only the first five columns are returned in arrays first, second,  */
/*              third, fourth and fifth. n is the number of values returned.       */
/*              Memory allocation for first, second, third, fourth, and fifth      */
/*              is done automatically; fields can be freed with a simple free().   */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int read_5c_file (char *filename, 
		  double **first, double **second, double **third, 
		  double **fourth, double **fifth, int *n)
{
  int max_columns=0, min_columns=0;
  double **data=NULL;
  int status=0;

  /* read file */
  if ( (status = ASCII_file2double (filename, n, 
				    &max_columns, &min_columns, &data)) != 0)
    return status;


  /* check, if at least five columns */
  if (min_columns < 5)   {
    ASCII_free_double (data, *n);
    return LESS_THAN_THREE_COLUMNS;
  }
  
  *first  = ASCII_column (data, *n, 0);
  *second = ASCII_column (data, *n, 1);
  *third  = ASCII_column (data, *n, 2);
  *fourth = ASCII_column (data, *n, 3);
  *fifth  = ASCII_column (data, *n, 4);
  
  ASCII_free_double (data, *n);

  return 0;
}


/***********************************************************************************/
/* Function: read_5c_file_float                                           @30_30i@ */
/* Description: Read an ASCII file with (at least) 5 columns to a float array.     */
/*              Only the first five columns are returned in arrays first, second,  */
/*              third, fourth, and fifth. n is the number of values returned.      */
/*              Memory allocation for first, second, third, and fourth is done     */
/*              automatically; fields can be freed with a simple free().           */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int read_5c_file_float (char *filename, 
			float **first, float **second, float **third, float **fourth, float **fifth, int *n)
{
  int max_columns=0, min_columns=0;
  float **data=NULL;
  int status=0;

  /* read file */
  if ( (status = ASCII_file2float (filename, n, 
				   &max_columns, &min_columns, &data)) != 0)
    return status;


  /* check, if at least five columns */
  if (min_columns < 5)   {
    ASCII_free_float (data, *n);
    return LESS_THAN_FIVE_COLUMNS;
  }
  
  *first  = ASCII_column_float (data, *n, 0);
  *second = ASCII_column_float (data, *n, 1);
  *third  = ASCII_column_float (data, *n, 2);
  *fourth = ASCII_column_float (data, *n, 3);
  *fifth  = ASCII_column_float (data, *n, 4);
  
  ASCII_free_float (data, *n);

  return 0;
}


/***********************************************************************************/
/* Function: read_6c_file                                                 @30_30i@ */
/* Description: Read an ASCII file with (at least) 6 columns.                      */
/*              Only the first six columns are returned in arrays first, second,   */
/*              third, fourth, fifth, and sixth. n is the number of values         */
/*              returned. Memory allocation for the result arrays                  */
/*              is done automatically; fields can be freed with a simple free().   */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int read_6c_file (char *filename, 
		  double **first, double **second, double **third, 
		  double **fourth, double **fifth, double **sixth, int *n)
{
  int max_columns=0, min_columns=0;
  double **data=NULL;
  int status=0;

  /* read file */
  if ( (status = ASCII_file2double (filename, n, 
				    &max_columns, &min_columns, &data)) != 0)
    return status;


  /* check, if at least six columns */
  if (min_columns < 6)   {
    ASCII_free_double (data, *n);
    return LESS_THAN_THREE_COLUMNS;
  }
  
  *first  = ASCII_column (data, *n, 0);
  *second = ASCII_column (data, *n, 1);
  *third  = ASCII_column (data, *n, 2);
  *fourth = ASCII_column (data, *n, 3);
  *fifth  = ASCII_column (data, *n, 4);
  *sixth  = ASCII_column (data, *n, 5);
  
  ASCII_free_double (data, *n);

  return 0;
}

/***********************************************************************************/
/* Function: read_7c_file                                                 @30_30i@ */
/* Description: Read an ASCII file with (at least) 6 columns.                      */
/*              Only the first six columns are returned in arrays first, second,   */
/*              third, fourth, fifth, and sixth. n is the number of values         */
/*              returned. Memory allocation for the result arrays                  */
/*              is done automatically; fields can be freed with a simple free().   */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int read_7c_file (char *filename, 
		  double **first, double **second, double **third, 
		  double **fourth, double **fifth, double **sixth, 
                  double **seventh, int *n)
{
  int max_columns=0, min_columns=0;
  double **data=NULL;
  int status=0;

  /* read file */
  if ( (status = ASCII_file2double (filename, n, 
				    &max_columns, &min_columns, &data)) != 0)
    return status;


  /* check, if at least seven columns */
  if (min_columns < 7)   {
    ASCII_free_double (data, *n);
    return LESS_THAN_THREE_COLUMNS;
  }
  
  *first  = ASCII_column (data, *n, 0);
  *second = ASCII_column (data, *n, 1);
  *third  = ASCII_column (data, *n, 2);
  *fourth = ASCII_column (data, *n, 3);
  *fifth  = ASCII_column (data, *n, 4);
  *sixth  = ASCII_column (data, *n, 5);
  *seventh = ASCII_column (data, *n, 6);

  ASCII_free_double (data, *n);

  return 0;
}

/***********************************************************************************/
/* Function: read_7c_file_float                                           @30_30i@ */
/* Description: Read an ASCII file with (at least) 6 columns.                      */
/*              Only the first six columns are returned in arrays first, second,   */
/*              third, fourth, fifth, and sixth. n is the number of values         */
/*              returned. Memory allocation for the result arrays                  */
/*              is done automatically; fields can be freed with a simple free().   */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int read_7c_file_float (char *filename, 
			float **first, float **second, float **third, 
			float **fourth, float **fifth, float **sixth, 
			float **seventh, int *n)
{
  int max_columns=0, min_columns=0;
  float **data=NULL;
  int status=0;

  /* read file */
  if ( (status = ASCII_file2float (filename, n, 
				    &max_columns, &min_columns, &data)) != 0)
    return status;


  /* check, if at least seven columns */
  if (min_columns < 7)   {
    ASCII_free_float (data, *n);
    return LESS_THAN_THREE_COLUMNS;
  }
  
  *first  = ASCII_column_float (data, *n, 0);
  *second = ASCII_column_float (data, *n, 1);
  *third  = ASCII_column_float (data, *n, 2);
  *fourth = ASCII_column_float (data, *n, 3);
  *fifth  = ASCII_column_float (data, *n, 4);
  *sixth  = ASCII_column_float (data, *n, 5);
  *seventh = ASCII_column_float (data, *n, 6);

  ASCII_free_float (data, *n);

  return 0;
}

/***********************************************************************************/
/* Function: read_8c_file                                                 @30_30i@ */
/* Description: Read an ASCII file with (at least) 8 columns.                      */
/*              Only the first eight columns are returned in arrays first, second, */
/*              third, fourth, fifth, sixth, seventh, and eigths. n is the number  */
/*              of values returned. Memory allocation for the result arrays        */
/*              is done automatically; fields can be freed with a simple free().   */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int read_8c_file (char *filename, 
		  double **first, double **second, double **third, 
		  double **fourth, double **fifth, double **sixth, 
		  double **seventh, double **eighth, int *n)
{
  int max_columns=0, min_columns=0;
  double **data=NULL;
  int status=0;

  /* read file */
  if ( (status = ASCII_file2double (filename, n, 
				    &max_columns, &min_columns, &data)) != 0)
    return status;


  /* check, if at least eight columns */
  if (min_columns < 8)   {
    ASCII_free_double (data, *n);
    return LESS_THAN_THREE_COLUMNS;
  }
  
  *first   = ASCII_column (data, *n, 0);
  *second  = ASCII_column (data, *n, 1);
  *third   = ASCII_column (data, *n, 2);
  *fourth  = ASCII_column (data, *n, 3);
  *fifth   = ASCII_column (data, *n, 4);
  *sixth   = ASCII_column (data, *n, 5);
  *seventh = ASCII_column (data, *n, 6);
  *eighth  = ASCII_column (data, *n, 7);
  
  ASCII_free_double (data, *n);

  return 0;
}

/***********************************************************************************/
/* Function: read_8c_file_float                                           @30_30i@ */
/* Description: Read an ASCII file with (at least) 8 columns.                      */
/*              Only the first eight columns are returned in arrays first, second, */
/*              third, fourth, fifth, sixth, seventh, and eigths. n is the number  */
/*              of values returned. Memory allocation for the result arrays        */
/*              is done automatically; fields can be freed with a simple free().   */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int read_8c_file_float (char *filename, 
			float **first, float **second, float **third, 
			float **fourth, float **fifth, float **sixth, 
			float **seventh, float **eighth, int *n)
{
  int max_columns=0, min_columns=0;
  float **data=NULL;
  int status=0;

  /* read file */
  if ( (status = ASCII_file2float (filename, n, 
				    &max_columns, &min_columns, &data)) != 0)
    return status;


  /* check, if at least eight columns */
  if (min_columns < 8)   {
    ASCII_free_float (data, *n);
    return LESS_THAN_THREE_COLUMNS;
  }
  
  *first   = ASCII_column_float (data, *n, 0);
  *second  = ASCII_column_float (data, *n, 1);
  *third   = ASCII_column_float (data, *n, 2);
  *fourth  = ASCII_column_float (data, *n, 3);
  *fifth   = ASCII_column_float (data, *n, 4);
  *sixth   = ASCII_column_float (data, *n, 5);
  *seventh = ASCII_column_float (data, *n, 6);
  *eighth  = ASCII_column_float (data, *n, 7);
  
  ASCII_free_float (data, *n);

  return 0;
}



/***********************************************************************************/
/* Function: read_9c_file                                                 @30_30i@ */
/* Description: Read an ASCII file with (at least) 9 columns.                      */
/*              Only the first eight columns are returned in arrays first, second, */
/*              third, fourth, fifth, sixth, seventh, eigths, and ninth. n is the  */
/*              number of values returned. Memory allocation for the result arrays */
/*              is done automatically; fields can be freed with a simple free().   */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int read_9c_file (char *filename, 
		  double **first, double **second, double **third, 
		  double **fourth, double **fifth, double **sixth, 
		  double **seventh, double **eighth, double **ninth, int *n)
{
  int max_columns=0, min_columns=0;
  double **data=NULL;
  int status=0;

  /* read file */
  if ( (status = ASCII_file2double (filename, n, 
				    &max_columns, &min_columns, &data)) != 0)
    return status;


  /* check, if at least nine columns */
  if (min_columns < 9)   {
    ASCII_free_double (data, *n);
    return LESS_THAN_THREE_COLUMNS;
  }
  
  *first   = ASCII_column (data, *n, 0);
  *second  = ASCII_column (data, *n, 1);
  *third   = ASCII_column (data, *n, 2);
  *fourth  = ASCII_column (data, *n, 3);
  *fifth   = ASCII_column (data, *n, 4);
  *sixth   = ASCII_column (data, *n, 5);
  *seventh = ASCII_column (data, *n, 6);
  *eighth  = ASCII_column (data, *n, 7);
  *ninth   = ASCII_column (data, *n, 8);
  
  ASCII_free_double (data, *n);

  return 0;
}



/***********************************************************************************/
/* Function: read_2r_file                                                 @30_30i@ */
/* Description: Read the first 2 rows of an ASCII file.                            */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int read_2r_file (char *filename,   
		  int *rows,        
		  int *max_columns, 
		  int *min_columns, 
		  double ***value)  
{
  char ***string = NULL;
  int status=0;
  int max_length=0;

  /* count rows and columns of ASCII file <filename> */
  if ( (status = ASCII_checkfile (filename, 
				  rows, 
				  min_columns, 
				  max_columns, 
				  &max_length)) != 0)
    return status;
  


  /* allocate memory for string array */
  if ( (status = ASCII_calloc_string (&string, 
				      *rows, 
				      *max_columns, 
				      max_length)) != 0 )
    return status;
  
  /* read ASCII file to string array */
  if ( (status = ASCII_readfile (filename, string)) != 0)
    return status;

  /* allocate memory for double array */
  if ( (status = ASCII_calloc_double (value, *rows, *max_columns)) != 0 )
    return status;

  /* convert string array to double array */
  if ( (status = ASCII_string2double (*value, string, *rows, *max_columns)) != 0)
    return status; 



  /* free memory of string array */
  ASCII_free_string (string, *rows, *max_columns);

  return 0;  /* everything ok */
} 






/***********************************************************************************/
/* Function: substr                                                       @30_30i@ */
/* Description: Create substring starting at position start with length            */
/*              length. Result is written to buffer (which MUST be                 */
/*              allocated before).                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  Pointer to the substring.                                                      */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

char *substr (char *buffer, char *string, int start, int length)
{
  strncpy (buffer, string+start, length);                    
  buffer[length]=0;  /* end of string */

  return buffer;
}





/***********************************************************************************/
/* Function: ASCII_parse                                                  @30_30i@ */
/* Description: Parse string to an array of single words. Memory for an array      */
/*              of string pointers is allocated automatically. array[i]            */ 
/*              points to the address of word #i in string!                        */
/*		Word separator is specified in separator.                          */
/*		Number of words is returned in number.                             */
/*		Characters following the comment character are ignored.            */
/* Parameters:                                                                      */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_parse (char *string, char *separator, char ***array, int *number)
{

  char *start=NULL;
  char *t=NULL;
  char *save=NULL;
  char **temp=NULL;

  /* save start address of string */
  start = string;

  /* save string */
  save = (char *) calloc (strlen(string)+1, sizeof (char));
  strcpy (save, string);


  /* reset number */
  *number=0;
  
  /* count words */
  if ( (t = strtok (string, separator) ) != NULL)  {  /* if not an empty line */ 
    if (!ASCII_comment(t[0]))  {  /* if not a comment     */
      (*number)++;
      while ( (t = strtok (NULL, separator) ) != NULL)  {
	if (ASCII_comment(t[0]))     /* if comment */
	  break;

	(*number)++;
      }
    }
  }

  if (*number==0)   {
    free(save);
    return 0;     /* no words found. but o.k. */ 
  }

  /* restore string */
  string = start;

  /* allocate memory for an array of *number character pointers */
  temp = (char **) calloc (*number, sizeof(char *));

  /* restore string */
  strcpy (string, save);

  /* reset *number */
  *number=0;  

  /* now set array pointers */
  if ( (t = strtok (string, separator) ) != NULL)  {
    if (!ASCII_comment(t[0]))  {             
      temp[(*number)++] = t;
      
      while ( (t = strtok (NULL, separator) ) != NULL)  {
	if (ASCII_comment(t[0]))    
	  break;

	temp[(*number)++] = t;
      }
    }
  }


  /* copy temp to array */
  *array = temp;

  free(save);
  return 0;   /* if o.k. */
}



/***********************************************************************************/
/* Function: ASCII_parsestring                                            @30_30i@ */
/* Description: For compatibility reasons: ASCII_parsestring is just               */
/*              a call to ASCII_parse() with field separator set to " \t\v\f\r\n"  */ 
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                @i30_30@ */
/***********************************************************************************/

int ASCII_parsestring (char *string, char ***array, int *number)
{
  int status=0;

  status = ASCII_parse (string, " \t\v\f\r\n", array, number);
  return status;
}



/***********************************************************************************/
/* Function: ASCII_comment                                                         */
/* Description: Check if an ASCII character is one of the comment characters.      */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/*  1, if comment, and 0, if not a comment.                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/***********************************************************************************/

static int ASCII_comment (char t) {

  if (t == ASCII_COMMENT_1 || t == ASCII_COMMENT_2)
    return 1;
  else 
    return 0;
}




/***********************************************************************************/
/* Function: ASCII_sortarray                                                       */
/* Description:                                                                    */
/*  Sort a 2D double array by a specified column.                                  */
/*                                                                                 */
/* Parameters:                                                                     */
/*  double **value:   2D array to be sorted; value[row][column]                    */
/*  int  rows:        number of rows                                               */
/*  int  columns:     number of columns                                            */
/*  int  key:         key column number (starting at 0!)                           */
/*  char down:        0: ascending order, >0: descending order                     */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                                 */
/***********************************************************************************/

int ASCII_sortarray (double **value, int rows, int columns, int index, char down)
{
  int i=0, j=0;
  sortarray_struct *temp = calloc(rows, sizeof(sortarray_struct));
  
  /* check if index columns exists */
  if (index > columns) {
    fprintf (stderr, "error, sort index = %d is larger than columns = %d\n",
	     index, columns);
    return -1;
  }

  /* copy 2D array to an array of structures which can be */
  /* passed to qsort(); each array element contains a     */
  /* row of the matrix and the sort index.                */

  for (i=0; i<rows; i++) {
    temp[i].x     = calloc (columns, sizeof(double));
    temp[i].index = index;
    
    for (j=0; j<columns; j++)
      temp[i].x[j] = value[i][j];
  } 

  /* the array is now stored in a structure which can be used as */
  /* input to qsort().                                           */
  
  if (down)
    qsort (temp, rows, sizeof(sortarray_struct), sortarray_fnc_dn);
  else 
    qsort (temp, rows, sizeof(sortarray_struct), sortarray_fnc_up);
  
  for (i=0; i<rows; i++)
    for (j=0; j<columns; j++)
      value[i][j] = temp[i].x[j];


  return 0;   /* if o.k. */
}



/***********************************************************************************/
/* Function: ASCII_sortarray_float                                                 */
/* Description:                                                                    */
/*  Sort a 2D float array by a specified column.                                   */
/*                                                                                 */
/* Parameters:                                                                     */
/*  float **value:   2D array to be sorted; value[row][column]                     */
/*  int  rows:        number of rows                                               */
/*  int  columns:     number of columns                                            */
/*  int  key:         key column number (starting at 0!)                           */
/*  char down:        0: ascending order, >0: descending order                     */
/*                                                                                 */
/* Return value:                                                                   */
/*  0  if o.k., <0 if error                                                        */
/*                                                                                 */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/*                                                                                 */
/***********************************************************************************/

int ASCII_sortarray_float (float **value, int rows, int columns, int index, char down)
{
  int i=0, j=0;
  sortarray_float_struct *temp = calloc(rows, sizeof(sortarray_struct));
  
  /* check if index columns exists */
  if (index > columns) {
    fprintf (stderr, "error, sort index = %d is larger than columns = %d\n",
	     index, columns);
    return -1;
  }

  /* copy 2D array to an array of structures which can be */
  /* passed to qsort(); each array element contains a     */
  /* row of the matrix and the sort index.                */

  for (i=0; i<rows; i++) {
    temp[i].x     = calloc (columns, sizeof(float));
    temp[i].index = index;
    
    for (j=0; j<columns; j++)
      temp[i].x[j] = value[i][j];
  } 

  /* the array is now stored in a structure which can be used as */
  /* input to qsort().                                           */
  
  if (down)
    qsort (temp, rows, sizeof(sortarray_struct), sortarray_fnc_dn);
  else 
    qsort (temp, rows, sizeof(sortarray_struct), sortarray_fnc_up);
  
  for (i=0; i<rows; i++)
    for (j=0; j<columns; j++)
      value[i][j] = temp[i].x[j];


  return 0;   /* if o.k. */
}



/***********************************************************************************/
/* Function: sortarray_fnc_up                                                      */
/* Description: Internal function, required by ASCII_sortarray();                  */
/*              see description of standard C library function qsort().            */
/*                                                                                 */
/* Parameters:                                                                     */    
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/***********************************************************************************/

static int sortarray_fnc_up (const void *x1, const void *x2)
{
  sortarray_struct *c1 = (sortarray_struct *) x1;
  sortarray_struct *c2 = (sortarray_struct *) x2;

  if (c1->x[c1->index] > c2->x[c2->index])
    return 1;

  if (c1->x[c1->index] < c2->x[c2->index])
    return -1;

  return 0;  
}



/***********************************************************************************/
/* Function: sortarray_fnc_dn                                                      */
/* Description: Internal function, required by ASCII_sortarray();                  */
/*              see description of standard C library function qsort().            */
/*                                                                                 */
/* Parameters:                                                                     */    
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author:                                                                         */
/***********************************************************************************/

static int sortarray_fnc_dn (const void *x1, const void *x2)
{
  sortarray_struct *c1 = (sortarray_struct *) x1;
  sortarray_struct *c2 = (sortarray_struct *) x2;

  if (c1->x[c1->index] < c2->x[c2->index])
    return 1;

  if (c1->x[c1->index] > c2->x[c2->index])
    return -1;

  return 0;  /* if equal */
}


