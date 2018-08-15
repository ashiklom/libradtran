/*--------------------------------------------------------------------
 * $Id: errors.h 2993 2013-11-29 15:50:11Z robert.buras $
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

#ifndef __errors_h
#define __errors_h

#define ERROR_POSITION __LINE__, (char*) __func__, __FILE__

int err_out (char *output, int status);

int mem_err_out (char *output, int line, char *fff, char *file);

int fct_err_out (int status, char *function_name, int line, char *func, char *filename );

int err_file_out ( int status, char *filename, int line, char *func, char *cfilename );

#endif /* _ERRORS_H */

