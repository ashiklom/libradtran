/*--------------------------------------------------------------------
 * $Id: table.c 2623 2011-12-23 10:52:38Z robert.buras $
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
#include "table.h"

/**************************************************************/
/* Read Stamnes table from file filename.                     */
/* zenith [0..n_zenith-1] is an array of solar zenith angles, */ 
/* yy [0..n_yy-1] an array of yys found in the       */
/* file. The two-dimensional integer array table holds the    */
/* corresponding table columns.                               */
/* Memory for arrays is allocated automatically.              */
/**************************************************************/

int read_table (char *filename, TABLE **table)
{
  int i=0, j=0, status=0;
  int rows=0, columns=0, max_columns=0;
  double **table_r=NULL;


  /* read file to two-dimensional double array table */
  status = ASCII_file2double (filename,
			      &rows, &max_columns, &columns,
			      &table_r);

  if (status != 0) 
    return status;

  /* check if a rectangular matrix */  
  if (columns != max_columns)  { 
    ASCII_free_double (table_r, rows);
    fprintf (stderr, " ... read_table(): %s is not a rectangular matrix!\n", filename);
    fprintf (stderr, " ... minimum number of columns: %d, maximum number of columns: %d\n", 
	     columns, max_columns);

    return NOT_A_VALID_TABLE;
  }
    

  /* allocate memory for table */
  *table = calloc (1, sizeof(TABLE));
  
  /* first row/column contains yy/xx */
  (*table)->n_xx  = columns-1;
  (*table)->n_yy  = rows-1;
 
  /* allocate memory for arrays */
  (*table)->xx     = (double *)  calloc ((*table)->n_xx, sizeof (double));
  (*table)->yy     = (double *)  calloc ((*table)->n_yy, sizeof (double));
  (*table)->table  = (double **) calloc ((*table)->n_yy, sizeof (double *));

  for (i=0; i<(*table)->n_yy; i++)
    (*table)->table[i] = (double *) calloc ((*table)->n_xx,  sizeof (double));
  
  /* copy first row to array xx */
  for (i=1; i<=(*table)->n_xx; i++)  
    (*table)->xx[i-1] = table_r[0][i];

  /* copy first column to array yy */
  for (j=1; j<=(*table)->n_yy; j++)  
    (*table)->yy[j-1] = table_r[j][0];



  /* copy the remainder to 2dim double array table */
  for (i=1; i<=(*table)->n_yy; i++) 
    for (j=1; j<=(*table)->n_xx; j++)
      (*table)->table[i-1][j-1] = table_r[i][j];
    

  /* free memory of double array */
  ASCII_free_double (table_r, rows);
 
  return 0;  /* if o.k. */
}

/**************************************************************/
/* free memory of data read from table.               */
/**************************************************************/

void free_table (TABLE *table)
{
  int i=0;

  for (i=0; i<table->n_yy; i++) 
    free (table->table[i]);

  free(table->table);

  free(table->xx);
  free(table->yy);

  free(table);
}

/***************************************************************/
/* table_calculate uses the data read from a table             */
/* to calculate the table value for given xx and yy.           */
/*                                                             */
/* xx[0..n_xx-1]  xx values found in table.                    */
/* yy[0..n_yy-1]  yy's found in  table.                        */
/* table[0..n_yy-1][0..n_xx-1]   table columns found in table. */
/* calc_xx   xx to be used for calculation                     */
/* calc_yy   yy to be used for calculation                     */
/* column                                                      */
/***************************************************************/

int table_calculate ( TABLE *table,
		       double calc_xx, double calc_yy,
		       double *column)
{
  int counter=0, status=0;
  double *a0=NULL, *a1=NULL, *a2=NULL, *a3=NULL;
  double *b0=NULL, *b1=NULL, *b2=NULL, *b3=NULL;
  double lower=0, upper=0;

  /* reset table column */
  *column = -1;

  /* look for neighbour points in  array of yys */
  while (counter < table->n_yy && table->yy[counter] < calc_yy)
    counter++;

  /* check if calc_yy within table limits */ 
  if (counter<=0 || counter >= table->n_yy)
    return REQUESTED_VALUE_OUT_OF_RANGE;


  /* calculate spline coefficients for table as function of xx */
  status = spline_coeffc (table->xx, table->table[counter-1], table->n_xx, 
			  &a0, &a1, &a2, &a3);

  if (status!=0)  {
    free(a0);
    free(a1);
    free(a2);
    free(a3);
    return ERROR_SPLINING;
  }
    

  status = spline_coeffc (table->xx, table->table[counter], table->n_xx, 
			  &b0, &b1, &b2, &b3);
  
  if (status!=0)  {
    free(a0);
    free(a1);
    free(a2);
    free(a3);
    free(b0);
    free(b1);
    free(b2);
    free(b3);
    return ERROR_SPLINING;
  }

  /* calculate interpolated data */
  status =  calc_splined_value (calc_xx, &lower, 
				table->xx, table->n_xx, 
				a0, a1, a2, a3);

  free(a0);
  free(a1);
  free(a2);
  free(a3);
  

  if (status!=0)  {
    free(b0);
    free(b1);
    free(b2);
    free(b3);
    return ERROR_SPLINING;
  }

  status =  calc_splined_value (calc_xx, &upper, 
				table->xx, table->n_xx, 
				b0, b1, b2, b3);

  free(b0);
  free(b1);
  free(b2);
  free(b3);


  if (status!=0) 
    return ERROR_SPLINING;


  /* now interpolate linearly for yy */
  *column = lower + (upper - lower) / (table->yy[counter] - table->yy[counter-1]) *
    (calc_yy - table->yy[counter-1]);

  return 0;  /* if o.k. */
}
