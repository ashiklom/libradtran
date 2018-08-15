/*--------------------------------------------------------------------
 * $Id: ocean.c 3252 2016-11-30 23:25:50Z Bernhard.Mayer $
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

/****************************************************************************/
/* This file is an interface for the Fortran ocean BRDF routine oceabrdf.   */
/* It turns out that DISORT2 requires many calls to this routine, something */
/* like 140,000 calls for an ordinary 16 streams disort run. As the         */
/* oceabrdf routine is quite complicated, the calculation of the BRDF       */
/* causes the major part of the computational cost of the disort call.      */
/*                                                                          */
/* In many applications, disort is called several times for the same        */
/* wavelength during one uvspec run (all correlated_k applications).        */
/* In this case, the same 140,000 oceabrdf calls are repeated for each      */
/* individual wavelength.                                                   */
/*                                                                          */
/* ocean stores all quadruplets (mu, mup, dphi, bdref) that are calculated  */
/* in a search tree (code courtesy of Kearnighan and Ritchie).              */
/* A search tree is an efficient way of storing and extracting large        */
/* amounts of data the amount of which is not known a priori. Data are      */
/* added with addtree() and may be extracted again with searchtree().       */
/* Hence, for each now input triple (mu, mup, dphi) it is checked if the    */
/* value is already available in the tree. If yes, it is got from there;    */
/* if not, it is calculated and stored in the tree.                         */
/*                                                                          */
/* The bad news is that such trees are most efficient if the data to be     */
/* stored in the tree come in completely random (unsorted) order. This is   */
/* unfortunately not the case here because DISORT2 calls oceabrdf for       */
/* well sorted mu, mup, and dphi values. Hence, in a first step, the data   */
/* are gathered in a different tree (rtree) during the first DISORT2 call   */
/* for each wavelength (DISORT2 has been modified to pass the number of     */
/* the DISORT2 run since since the start of uvspec to oceanc). A function   */
/* rtree2tree then copies these data to their final destination in fully    */
/* randomized order.                                                        */
/*                                                                          */
/* A first test showed a speed improvement of a factor 3-4 for a            */
/* 16 streams DISORT2 run for the 32 KATO wavelength bands.                 */
/*                                                                          */
/* Removed this after nearly 15 years. It doesn't help at all with REPTRAN  */
/* where each wavelength is calculated only once and BRDFs may not be       */
/* recycled.                                                                */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "ocean.h"
#include "uvspecrandom.h"
#include "f77-uscore.h"
#include "fortran_and_c.h"

/* search tree node; this will become a 3D tree: "value" stores the  */
/* polar angle mu while "next" points to another tree where "value"  */
/* stores mup and "next" points to a third tree where "value" stores */
/* dphi and bdref finally holds the BRDF.                            */

struct node {
  float value;
  float bdref;
  
  struct node *left;
  struct node *right;
  struct node *next;
};


/* randomizing tree node, similar to the node above but one-dimensional */
struct rnode {
  float data[4];
  
  struct rnode *left;
  struct rnode *right;
};


/* prototypes */

static struct node *addtree  (struct node *p, float x, float y, float z, float bdref);
static struct node *addtree2 (struct node *p, float y, float z, float bdref);
static struct node *addtree3 (struct node *p, float z, float bdref);

static float searchtree  (struct node *p, float x, float y, float z);
static float searchtree2 (struct node *p, float y, float z);
static float searchtree3 (struct node *p, float z);

static void treefree  (struct node *p);
static void treefree2 (struct node *p);
static void treefree3 (struct node *p);

struct rnode *addrtree (struct rnode *p, float mu, float mup, float dphi, float bdref);
static void rtree2tree (struct rnode *source, struct node **dest);
static void rtreefree  (struct rnode *p);

/* prototype of the Fortran oceabrdf function */
void F77_FUNC  (oceabrdf, OCEABRDF) (float *wvnmlo, float *wvnmhi, float *mu, float *mup, float *dphi, 
			  float *u10, float *pcl, float *xsal, float *paw, float *bdref);


void F77_FUNC (rmatr, RMATR) (ref_complex *cn1, ref_complex *cn2,
                              double *sigma2, double *mu_inc, double *phi_inc,
                              double *mu_esc, double *phi_esc,
                              double *refl_mat);


void F77_FUNC (indwat, INDWAT) (float *lambda, float *dummy_sal, 
                                float *refre, float *refim);

/* simple C wrapper function for the standard oceabrdf call */ 
double oceabrdfc (float wvnmlo, float wvnmhi, 
		  float mu, float mup, float dphi, 
		  float u10, float pcl, float xsal,
                  float paw) 
{
  float bdref=0;

  F77_FUNC  (oceabrdf, OCEABRDF) (&wvnmlo, &wvnmhi, &mu, &mup, &dphi, &u10, &pcl, &xsal, &paw, &bdref);
  return (double) bdref;
}


/* add a node to the search tree, see Kernighan and Ritchie */
static struct node *addtree (struct node *p, float x, float y, float z, float bdref)
{
  if (p==NULL) {
    p = calloc (1, sizeof (struct node));
    p->value = x;
    p->next = addtree2 (p->next, y, z, bdref);
  } 
  else {
    if (x==p->value) {
      p->next = addtree2 (p->next, y, z, bdref);
    }
    else {
      if (x<p->value)
	p->left = addtree (p->left, x, y, z, bdref);
      else 
	p->right = addtree (p->right, x, y, z, bdref);
    }
  }
  return p;
}


/* add a node to the first subtree; required by addtree */
static struct node *addtree2 (struct node *p, float y, float z, float bdref)
{
  if (p==NULL) {
    p = calloc (1, sizeof (struct node));
    p->value = y;
    p->bdref=bdref;
    p->next = addtree3 (p->next, z, bdref);
  } 
  else {
    if (y==p->value) {
      p->next = addtree3 (p->next, z, bdref);
    }
    else {
      if (y<p->value)
	p->left = addtree2 (p->left, y, z, bdref);
      else 
	p->right = addtree2 (p->right, y, z, bdref);
    }
  }
  return p;
}


/* add a node to the second subtree; required by addtree */
static struct node *addtree3 (struct node *p, float z, float bdref)
{
  if (p==NULL) {
    p = calloc (1, sizeof (struct node));
    p->value = z;
    p->bdref=bdref;
  } 
  else {
    if (z<p->value)
      p->left = addtree3 (p->left, z, bdref);

    if (z>p->value)
      p->right = addtree3 (p->right, z, bdref);
  }
  return p;
}



/* search the search tree, see Kernighan and Ritchie */
static float searchtree (struct node *p, float x, float y, float z)
{
  if (p==NULL) {
    return -999;
  } 
  else {
    if (x==p->value) {
      return searchtree2 (p->next, y, z);
    }
    else {
      if (x<p->value)
	return searchtree (p->left, x, y, z);
      else 
	return searchtree (p->right, x, y, z);
    }
  }
}


/* search the first subtree, required by searchtree */
static float searchtree2 (struct node *p, float y, float z)
{
  if (p==NULL) {
    return -999;
  } 
  else {
    if (y==p->value) {
      return searchtree3 (p->next, z);
    }
    else {
      if (y<p->value)
	return searchtree2 (p->left, y, z);
      else 
	return searchtree2 (p->right, y, z);
    }
  }
}


/* search the second subtree, required by searchtree */
static float searchtree3 (struct node *p, float z)
{
  if (p==NULL) {
    return -999;
  } 
  else {
    if (z==p->value) {
      return p->bdref;
    }
    else {
      if (z<p->value)
	return searchtree3 (p->left, z);
      else 
	return searchtree3 (p->right, z);
    }
  }
}




/* free randomization tree */
static void rtreefree (struct rnode *p)
{
  if (p != NULL) {
    rtreefree (p->left);
    rtreefree (p->right);
    free (p);
  }
}


/* free search tree */
static void treefree (struct node *p)
{
  if (p != NULL) {
    treefree (p->left);
    treefree (p->right);

    treefree2 (p->next);
    free (p);
  }
}

/* free first subtree, required by treefree */
static void treefree2 (struct node *p)
{
  if (p != NULL) {
    treefree2 (p->left);
    treefree2 (p->right);

    treefree3 (p->next);
    free (p);
  }
}

/* free second subtree, required by treefree */
static void treefree3 (struct node *p)
{
  if (p != NULL) {
    treefree3 (p->left);
    treefree3 (p->right);
    free (p);
  }
}


/* copy randomization tree to search tree */
static void rtree2tree (struct rnode *source, struct node **dest)
{
  if (source!=NULL) {
    rtree2tree (source->left, dest);
    *dest = addtree (*dest, source->data[0], source->data[1], source->data[2], source->data[3]);
    rtree2tree (source->right, dest);
  }
}



/* add node to the radomization tree; new nodes are added at random branches; */
/* that is, the tree is completely unsorted which is what we want             */
struct rnode *addrtree (struct rnode *p, float mu, float mup, float dphi, float bdref)
{
  if (p==NULL) {
    p = calloc (1, sizeof (struct rnode));
    p->data[0] = mu;
    p->data[1] = mup;
    p->data[2] = dphi;
    p->data[3] = bdref;
  } 
  else {
    if (uvspec_random() < 0.5)
      p->left  = addrtree (p->left, mu, mup, dphi, bdref);
    else 
      p->right = addrtree (p->right, mu, mup, dphi, bdref);
  }
  return p;
}



/* prototype of the Fortran oceabrdf function */
void F77_FUNC  (oceabrdf, OCEABRDF) (float *wvnmlo, float *wvnmhi, float *mu, float *mup, float *dphi, 
				     float *u10, float *pcl, float *xsal, float *paw, float *bdref);

/* new version without search tree, after November 30, 2016 */
float ocean_brdf (float wvnmlo, float wvnmhi, float mu, float mup, float dphi, 
		  float u10, float pcl, float xsal, int firstcall)
{
  float paw;
  float bdref=0.0;

  /* NOTE! paw = 0.0 means that wind comes from sun! This cannot be fixed! BCA */
  paw = 0.0;

  /* if anything (wavelength, salinity etc has changed, */
  /* we need to create a completely new search tree     */

  F77_FUNC  (oceabrdf, OCEABRDF) ( &wvnmlo, &wvnmhi, &mu, &mup, &dphi, &u10, &pcl, &xsal, &paw, &bdref );
  return bdref;
}

/* old version with search tree, before November 30, 2016 */
float ocean_brdf_old (float wvnmlo, float wvnmhi, float mu, float mup, float dphi, 
		      float u10, float pcl, float xsal, int firstcall)
{
  static struct node  *root  = NULL;
  static struct rnode *rtree = NULL;
  static float swvnmlo=0, swvnmhi=0, su10=0, spcl=0, sxsal=0;
  static int randomize=1, randomized=0, randomrun=0;
  static float paw;
  float bdref=0.0;

  /* NOTE! paw = 0.0 means that wind comes from sun! This cannot be fixed! BCA */
  paw = 0.0;

  /* if anything (wavelength, salinity etc has changed, */
  /* we need to create a completely new search tree     */

  if ( swvnmlo != wvnmlo ||
       swvnmhi != wvnmhi ||
       su10    != u10    ||
       spcl    != pcl    ||
       sxsal   != xsal)  {
    
    swvnmlo = wvnmlo;
    swvnmhi = wvnmhi;
    su10    = u10;
    spcl    = pcl;  
    sxsal   = xsal;

    randomize  = 1;
    randomized = 0;
    randomrun = firstcall;

    /* free memory of the search and randomization trees */
    rtreefree (rtree);
    rtree = NULL;
    
    treefree (root);
    root = NULL;
  }
    
  if (firstcall != randomrun)
    randomize=0;

  /* add data to randomization tree during the first DISORT2 call */
  if (randomize) {
    F77_FUNC  (oceabrdf, OCEABRDF) ( &wvnmlo, &wvnmhi, &mu, &mup, &dphi, &u10, &pcl, &xsal, &paw, &bdref );

    rtree = addrtree (rtree, mu, mup, dphi, bdref);
    return bdref;
  }


  /* after creating the randomization tree, copy the data to the search tree */
  if (!randomized) {
    rtree2tree (rtree, &root);
	
    randomized = 1;
  }


  /* this is the main function of the program: */
  /* look if the data exist in the tree ...    */
  bdref = searchtree (root, mu, mup, dphi);
    

  /* ... if not, calculate them and store them in the search tree */
  if (bdref==-999) {
    
    F77_FUNC  (oceabrdf, OCEABRDF) ( &wvnmlo, &wvnmhi, &mu, &mup, &dphi, &u10, &pcl, &xsal, &paw, &bdref );
    
    root = addtree (root, mu, mup, dphi, bdref);
  }

  return bdref;
}





/***********************************************************************************/
/* Function: bpdf_tsang                                                            */
/* Description:                                                                    */
/*      Calculate reflection matrix (Tsang, 1985) using Mishchenko's code.         */
/*      The reflection matrix is for an ocean surface and depends on the wind      */
/*      speed. Shadowing of the waves is taken into account.                       */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Claudia Emde                                                            */
/*                                                                        @i62_30@ */
/***********************************************************************************/
int bpdf_tsang(double wind_speed, float wavelength, double mu_inc, double phi_inc,
               double mu_esc, double phi_esc, double ***refl_mat){
  
  int status = 0;
  
  int i1=0, i2=0, j=0;
  
  ref_complex cn1; 
  ref_complex cn2;

  double *tmp_refmat; 
  double sigma2=0.0;
  
  float dummy_sal=-1;
  float refre=0.0;
  float refim=0.0; 
  /*float T=288.200; surface temperature, only needed in Warren's */
  /* refractive index is used. */

  /* mean surface slope, should be calculated from wind speed 
     according to Eq. (18) Mishchenko and Travis (1997) */
  /* Numbers correspond to original Cox and Munk (1954) definition */
  sigma2=0.5*(0.003+0.00512*wind_speed);
  
  tmp_refmat = c2fortran_2D_double_ary(4,4,*refl_mat);

  /* To be very accurate, cn1 should be set according to refractive*/
  /* index profile */
  cn1.re=1.0;
  cn1.im=0.0; 

  
  /* wavelength dependent refractive index of water */
 
  /* this is used in Cox and Munk BRDF  */
  F77_FUNC (indwat, INDWAT) (&wavelength, &dummy_sal, &refre, &refim);
  
  /* This function uses the Warren database */
  /* F77_FUNC (indwat, WREFWAT) (&wavelength, &T, &refre, &refim); */
  cn2.re=refre;
  cn2.im=refim;

  /* Set values for refractive index, uncomment these lines if you want to do this */
  /* cn2.re=1.34; */
  /* cn2.im=0.0; */
   
  /*convention is clockwise in MYSTIC, anticlockwise in Mishchenkos program!*/
  /* phi must be replaced by -phi */
  phi_inc*=-1.;
  phi_esc*=-1.;
  F77_FUNC (rmatr, RMATR) (&cn1, &cn2, &sigma2, &mu_inc, &phi_inc, 
                           &mu_esc, &phi_esc, tmp_refmat);
  
  /* The indices are correct this way, the function *fortran2c_2D_double_ary_noalloc* transforms the matrix */
  for (i1=0; i1<4; i1++)
    for (i2=0; i2<4; i2++) 
      (*refl_mat)[i2][i1] = tmp_refmat[j++];

  free(tmp_refmat);
  
  return status; 
}

void F77_FUNC (oceanfort, OCEANFORT) (float *wvnmlo, float *wvnmhi, float *mu, float *mup, float *dphi, 
				      float *u10, float *pcl, float *xsal, int *firstcall, float *bdref)
{
  *bdref = (float) ocean_brdf (*wvnmlo, *wvnmhi, *mu, *mup, *dphi, 
			       *u10, *pcl, *xsal, *firstcall);
}
