/************************************************************************/
/* avhrr_kratz.h                                                        */
/*                                                                      */
/* C header file for David Kratz' AVHRR parameterization                */
/*                                                                      */
/* Author: Bernhard Mayer,                                              */
/*         NCAR, bmayer@ucar.edu                                        */
/*                                                                      */
/*----------------------------------------------------------------------*/
/* Copyright (C) 1999 Bernhard Mayer                                    */
/*                                                                      */
/* This program is free software; you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation; either version 1, or (at your option)  */
/* any later version.                                                   */
/*                                                                      */
/* This program is distributed in the hope that it will be useful,      */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE. See the         */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* To obtain a copy of the GNU General Public License write to the      */
/* Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139,   */
/* USA.                                                                 */
/*----------------------------------------------------------------------*/
/************************************************************************/

#ifndef __avhrr_kratz_h
#define __avhrr_kratz_h

#if defined (__cplusplus)
extern "C" {
#endif

#include "f77-uscore.h"


/* prototypes */
/* molecular weights, required for David Kratz' AVHRR parameterization */
/* the numbers are taken from Kratz' routines for consitency reasons   */

#define KRATZ_MAXNTAU 20

#define H2O_MOLWGHT  18.01534
#define  O3_MOLWGHT  47.9982
#define CO2_MOLWGHT  44.00995
#define CH4_MOLWGHT  16.04303
#define N2O_MOLWGHT  44.0128
#define  O2_MOLWGHT  31.9988
#define F11_MOLWGHT  137.3685    /* CFCl3  */
#define F12_MOLWGHT  120.9139    /* CF2Cl2 */
#define NAVOGADRO    6.023e+23   /* as defined by David Kratz! */


void F77_FUNC  (avhrr11, AVHRR11) (float *z0, float *p0, float *t0, float *u0, float *ux, 
			 int *nlev, float *co2, float *n2o, 
			 float *f11, float *f12, float *ch4, 
			 float *tmp_od, float *tmp_wght, int *ntau);
void F77_FUNC  (avhrr12, AVHRR12) (float *z0, float *p0, float *t0, float *u0, float *ux, 
			 int *nlev, float *co2, float *n2o, 
			 float *f11, float *f12, float *ch4, 
			 float *tmp_od, float *tmp_wght, int *ntau);
void F77_FUNC  (avhrr13, AVHRR13) (float *z0, float *p0, float *t0, float *u0, float *ux, 
			 int *nlev, float *co2, float *n2o, 
			 float *f11, float *f12, float *ch4, 
			 float *tmp_od, float *tmp_wght, int *ntau);
void F77_FUNC  (avhrr14, AVHRR14) (float *z0, float *p0, float *t0, float *u0, float *ux, 
			 int *nlev, float *co2, float *n2o, 
			 float *f11, float *f12, float *ch4, 
			 float *tmp_od, float *tmp_wght, int *ntau);
void F77_FUNC  (avhrr15, AVHRR15) (float *z0, float *p0, float *t0, float *u0, float *ux, 
			 int *nlev, float *co2, float *n2o, 
			 float *f11, float *f12, float *ch4, 
			 float *tmp_od, float *tmp_wght, int *ntau);

void F77_FUNC  (avhrr21, AVHRR21) (float *z0, float *p0, float *t0, float *u0, float *ux, 
			 int *nlev, float *co2, float *n2o, 
			 float *f11, float *f12, float *ch4, 
			 float *tmp_od, float *tmp_wght, int *ntau);
void F77_FUNC  (avhrr22, AVHRR22) (float *z0, float *p0, float *t0, float *u0, float *ux, 
			 int *nlev, float *co2, float *n2o, 
			 float *f11, float *f12, float *ch4, 
			 float *tmp_od, float *tmp_wght, int *ntau);
void F77_FUNC  (avhrr23, AVHRR23) (float *z0, float *p0, float *t0, float *u0, float *ux, 
			 int *nlev, float *co2, float *n2o, 
			 float *f11, float *f12, float *ch4, 
			 float *tmp_od, float *tmp_wght, int *ntau);
void F77_FUNC  (avhrr24, AVHRR24) (float *z0, float *p0, float *t0, float *u0, float *ux, 
			 int *nlev, float *co2, float *n2o, 
			 float *f11, float *f12, float *ch4, 
			 float *tmp_od, float *tmp_wght, int *ntau);

void F77_FUNC  (avhrr31, AVHRR31) (float *z0, float *p0, float *t0, float *u0, float *ux, 
			 int *nlev, float *co2, float *n2o, 
			 float *f11, float *f12, float *ch4, 
			 float *tmp_od, float *tmp_wght, int *ntau);
void F77_FUNC  (avhrr32, AVHRR32) (float *z0, float *p0, float *t0, float *u0, float *ux, 
			 int *nlev, float *co2, float *n2o, 
			 float *f11, float *f12, float *ch4, 
			 float *tmp_od, float *tmp_wght, int *ntau);
void F77_FUNC  (avhrr33, AVHRR33) (float *z0, float *p0, float *t0, float *u0, float *ux, 
			 int *nlev, float *co2, float *n2o, 
			 float *f11, float *f12, float *ch4, 
			 float *tmp_od, float *tmp_wght, int *ntau);
void F77_FUNC  (avhrr34, AVHRR34) (float *z0, float *p0, float *t0, float *u0, float *ux, 
			 int *nlev, float *co2, float *n2o, 
			 float *f11, float *f12, float *ch4, 
			 float *tmp_od, float *tmp_wght, int *ntau);
void F77_FUNC  (avhrr35, AVHRR35) (float *z0, float *p0, float *t0, float *u0, float *ux, 
			 int *nlev, float *co2, float *n2o, 
			 float *f11, float *f12, float *ch4, 
			 float *tmp_od, float *tmp_wght, int *ntau);

void F77_FUNC  (avhrr41, AVHRR41) (float *z0, float *p0, float *t0, float *u0, float *ux, 
			 int *nlev, float *co2, float *n2o, 
			 float *f11, float *f12, float *ch4, 
			 float *tmp_od, float *tmp_wght, int *ntau);

void F77_FUNC  (avhrr51, AVHRR51) (float *z0, float *p0, float *t0, float *u0, float *ux, 
			 int *nlev, float *co2, float *n2o, 
			 float *f11, float *f12, float *ch4, 
			 float *tmp_od, float *tmp_wght, int *ntau);



#if defined (__cplusplus)
}
#endif

#endif






