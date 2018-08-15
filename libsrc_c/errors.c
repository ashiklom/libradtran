/*--------------------------------------------------------------------
 * $Id: errors.c 2993 2013-11-29 15:50:11Z robert.buras $
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
#include "errors.h"

/***********************************************************************************/
/* Functions: err_out                                                     @62_30i@ */
/* Description:                                                                    */
/*   short cuts error output and return status                                     */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int err_out (char *output, int status)
{
  fprintf (stderr,output,status);
  return status;
}


/***********************************************************************************/
/* Functions: mem_err_out                                                 @62_30i@ */
/* Description:                                                                    */
/*   short cuts memory allocation error output and return status                   */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int mem_err_out (char *output, int line, char *fff, char *file)
{
  fprintf (stderr,"Error when trying to allocate %s in (line %d, function '%s' in '%s')\n",
	   output, line, fff, file);
  return -1;
}


/***********************************************************************************/
/* Functions: fct_err_out                                                 @62_30i@ */
/* Description:                                                                    */
/*   short cuts error output and return status.                                    */
/*                                                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Robert Buras                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int fct_err_out (int status, char *function_name, int line, char *func, char *filename )
{
  fprintf (stderr, "Error %d returned by function %s in (line %d, function '%s', file '%s')\n",
	   status, function_name, line, func, filename);
  return status;
}


/***********************************************************************************/
/* Function: err_file_out                                                 @62_30i@ */
/* Description:                                                                    */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: RPB                                                                     */
/*                                                                        @i62_30@ */
/***********************************************************************************/

int err_file_out ( int status, char *filename, int line, char *func, char *cfilename )
{
  fprintf (stderr,"Error %d returned when reading file %s in (line %d, function '%s', file '%s')\n", status, filename, line, func, cfilename);
  return status;
}
