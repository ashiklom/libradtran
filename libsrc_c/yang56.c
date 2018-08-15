/*--------------------------------------------------------------------
 * $Id: yang56.c 3146 2015-08-04 18:33:14Z robert.buras $
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


/* This module was implemented by Albano Gonzales, */
/* eMail aglezf@ull.es                             */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "yang56.h"
#include "ascii.h"

/*agf*/
#define MAXWL 80   /* Number of wavelengths in tables (max) */

#define EXT56NF -2  /* File ext56.table not found    */
#define SSA56NF -3  /* File ext56.table not found    */
#define G56NF   -4  /* File ext56.table not found    */
#define WLOUT   -5  /* Wavelength out of range       */
#define HABITNF -6  /* Habit not found               */
#define REFFOUT -7  /* Effective radius out of range */
/*agf*/

/* This is a fix needed at MIM, need to test whether this also works on DLR */
#if COMPILING_CONDOR_AT_MIM
#define fscanf __fscanf
#endif

/* prototypes of internal functions */
static double dm2area (double dm, char *habit);
static double dm2vol  (double dm, char *habit);


/*Global variables for 8-13um parametrization (only for hexagons)*/
static float wlp[12]={8.0,8.5,9.0,9.5,10.0,10.5,10.8,11.0,11.5,12.0,12.5,13.0};
static float qep[3][12]=
  {{1.989,1.979,1.982,1.992,1.984,2.032,2.036,2.029,2.012,2.006,2.000,1.998},
   {4.336,6.117,5.219,3.891,5.208,-3.643,-6.777,-5.179,-1.431,1.283,2.787,3.201},
   {-27.76,-48.77,-42.48,-37.43,-81.94,-40.75,5.749,3.789,-10.61,-28.0,-33.23,-26.64}};
static float ssp[3][12]=
  {{0.5294,0.5205,0.5198,0.5191,0.5116,0.523,0.5384,0.5475,0.5562,0.5633,0.5556,0.5678},
   {0.2536,1.421,1.178,1.04,1.77,-0.8619,-2.454,-2.552,-1.974,-1.919,-1.605,-2.051},
   {15.41,7.036,9.676,10.63,-1.811,1.271,12.25,14.17,9.011,9.232,7.246,12.37}};
static float gp[2][12]=
  {{0.8317,0.8217,0.8399,0.856,0.8917,0.9271,0.9092,0.908,0.8922,0.8519,0.8553,0.829},
   {0.03279,0.03565,0.0312,0.02779,0.02092,0.01387,0.01771,0.01501,0.01394,0.02355,0.01932,0.02682}};



/*agf*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++
Function yang56_old
 (original code by Albano Gonzales)
 Computes the parametrization of single scattering properties
 for 7 habits, for 0.2<wavelength<4.8 um.

Input:  
 wavel -> wavelength
 reff ->  effective radius (could be computed using de2re)
 habit -> name of habit: Solid-Column 0, Hollow-Column 1, Rough-Aggregate 2,
 		         Rosette-4 3, Rosette-6 4, Plate 5 or Dendrite?
			 (case insensitive)

Output:
 *exti -> pointer to interp. volume extinction coefficient per ice water content.
 *ssai -> pointer to interp. single scattering albedo.
 *gi ->   pointer to interp. asymmetry parameter.
 *g2i ->  pointer to interp. 2nd asymmetry parameter (double Henyey-Greenstaein).
 *fi ->   pointer to interp. factor f in double Henyey-Greenstaein function.

++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


int yang56_old (float wavel, float reff, int nhabit, char *path, float *exti, 
		float *ssai, float *gi, float *g2i, float *fi)
{
  static short filere=0;  /*Only read data files once*/
  static short flagwl=0;
  static float ext[IC_HABIT_NN][4][MAXWL];
  static float ssa[IC_HABIT_NN][4][MAXWL];
  static float g[IC_HABIT_NN][6][MAXWL];
  static float wl[MAXWL];   
  static char habit_class[IC_HABIT_NN][80];
  static int index;
  static float wlmin, wlmax;
  
  char habit_name[80];
  FILE *fp;

  char filename[256];
 
  int i,j;
  float auxd, auxd1, auxd2, maxreff, minreff;
  int ini,fin,med;

  /* check whether reff is in the right range */
  switch (nhabit)  {
  case IC_HABIT_SOLID_COLUMN:   
    maxreff=MAX_REFF_SC; 
    minreff=MIN_REFF_SC; 
    break;
  case IC_HABIT_HOLLOW_COLUMN:  
    maxreff=MAX_REFF_HC; 
    minreff=MIN_REFF_HC; 
    break;
  case IC_HABIT_ROUGH_AGGREGATE:     
    maxreff=MAX_REFF_RA; 
    minreff=MIN_REFF_RA; 
    break;
  case IC_HABIT_ROSETTE_4:      
    maxreff=MAX_REFF_R4; 
    minreff=MIN_REFF_R4; 
    break;
  case IC_HABIT_ROSETTE_6:      
    maxreff=MAX_REFF_R6; 
    minreff=MIN_REFF_R6; 
    break;
  case IC_HABIT_PLATE:         
    maxreff=MAX_REFF_PL; 
    minreff=MIN_REFF_PL; 
    break;
  default:
    fprintf (stderr, "Error, unknown habit %d\n", nhabit);
    return HABITNF;
  }

  if (reff>maxreff) {
    fprintf (stderr, "Error, effective radius %f larger than maximum\n",reff);
    fprintf (stderr, "allowed value %6.2f for this habit!\n", maxreff);
    return REFFOUT;
  }
  if (reff<minreff) {
    fprintf (stderr, "Error, effective radius %f smaller than minimum\n",reff);
    fprintf (stderr, "allowed value %6.2f for this habit!\n", minreff);
    return REFFOUT;
  }

  /*Read the data files*/
  if(filere==0)
  {
    filere++;

    /*extinction coefficient*/
    sprintf(filename,"%sext56.table",path);
    
    if ((fp=fopen(filename,"r"))==NULL) {
      fprintf(stderr,"Error, %s not found\n", filename);
      return EXT56NF;
    }
    
    index=0;

    if ( ! fscanf(fp,"%s",habit_class[0]) ) {
      fprintf (stderr, "Error reading file %s\n", filename);
      return -1;
    }
    fseek(fp,0,0);
    do{
      if ( ! fscanf(fp,"%s %f %f %f %f %f", habit_name, wl+index, 
		    &ext[0][0][index], &ext[0][1][index], &ext[0][2][index], 
		    &ext[0][3][index]) ) {
	fprintf (stderr, "Error reading file %s\n", filename);
	return -1;
      }
      index++;

    }while(strcasecmp(habit_name,habit_class[0])==0);

    index--;

    ext[1][0][0] = ext[0][0][index]; 
    ext[1][1][0] = ext[0][1][index]; 
    ext[1][2][0] = ext[0][2][index]; 
    ext[1][3][0] = ext[0][3][index]; 

    for(j=1;j<IC_HABIT_NN;j++)
      for(i=1;i<index;i++)
        if ( ! fscanf(fp,"%s %f %f %f %f %f",habit_name,&auxd, 
		      &ext[j][0][i], &ext[j][1][i], &ext[j][2][i], 
		      &ext[j][3][i]) ) {
	  fprintf (stderr, "Error reading file %s\n", filename);
	  return -1;
	}

    wlmin=wl[0];
    wlmax=wl[index-1];
   
    fclose(fp);

    /*------------------------------------*/
    /*single scattering albedo*/
    sprintf(filename,"%sssa56.table",path);
    
    if ((fp=fopen(filename,"r"))==NULL) {
      fprintf(stderr,"Error, %s not found\n", filename);
      return SSA56NF;
    }

    for(j=0;j<IC_HABIT_NN;j++)
      for(i=0;i<index;i++)
        if ( ! fscanf(fp,"%s %f %f %f %f %f",habit_name,&auxd, 
               &ssa[j][0][i], &ssa[j][1][i], &ssa[j][2][i], 
		      &ssa[j][3][i]) ) {
	  fprintf (stderr, "Error reading file %s\n", filename);
	  return -1;
	}

    fclose(fp);

    /*------------------------------------*/
    /*asymetry parameter*/
    sprintf(filename,"%sg56.table",path);
    
    if ((fp=fopen(filename,"r"))==0) {
      fprintf(stderr,"Error, %s not found\n", filename);
      return G56NF;
    }

    for(j=0;j<IC_HABIT_NN;j++)
      for(i=0;i<index;i++)
        if ( ! fscanf(fp,"%s %f %f %f %f %f %f %f",habit_name,&auxd, 
		      &g[j][0][i], &g[j][1][i], &g[j][2][i], 
		      &g[j][3][i], &g[j][4][i], &g[j][5][i]) ) {
	  fprintf (stderr, "Error reading file %s\n", filename);
	  return -1;
	}

    fclose(fp);

  }
 
  /*Linear interpolation of parameters*/
  if(((wavel<wlmin)||(wavel>wlmax))&&((wavel<wlp[0])||(wavel>wlp[11]+0.5)))
  {
    if(!flagwl){
      fprintf(stderr,"Error, wavelength %f micron out of yang parametrization range\n",
	      wavel);
      flagwl++;
    }
    return WLOUT;
    /*
    if(wavel<wlmin)
      wavel=wlmin;
    else
      wavel=wlmax;
      */
  }

  if(wavel<wlp[0]){
    ini=0;
    fin=index-1;
    med=(fin+ini)/2;
    
    while((med!=ini)&&(med!=fin))
      {
	if(wavel<=wl[med])
	  fin=med;
	else
	  ini=med;
	med=(fin+ini)/2;
      }
    
    auxd=1.0/reff;
    auxd1=ext[nhabit][0][ini]+ext[nhabit][1][ini]*auxd+
      ext[nhabit][2][ini]*auxd*auxd+ext[nhabit][3][ini]*auxd*auxd*auxd;  
    auxd2=ext[nhabit][0][fin]+ext[nhabit][1][fin]*auxd+
      ext[nhabit][2][fin]*auxd*auxd+ext[nhabit][3][fin]*auxd*auxd*auxd;
    *exti=auxd1+(wavel-wl[ini])*(auxd2-auxd1)/(wl[fin]-wl[ini]);

    auxd=reff;
    auxd1=ssa[nhabit][0][ini] + ssa[nhabit][1][ini]*auxd + ssa[nhabit][2][ini]*(auxd*auxd) + ssa[nhabit][3][ini]*(auxd*auxd*auxd);  
    auxd2=ssa[nhabit][0][fin] + ssa[nhabit][1][fin]*auxd + ssa[nhabit][2][fin]*auxd*auxd + ssa[nhabit][3][fin]*auxd*auxd*auxd;
    *ssai=auxd1+(wavel-wl[ini])*(auxd2-auxd1)/(wl[fin]-wl[ini]);
    
    auxd=reff;
    auxd1=g[nhabit][0][ini]+g[nhabit][1][ini]*auxd+
      g[nhabit][2][ini]*(auxd*auxd)+g[nhabit][3][ini]*(auxd*auxd*auxd);  
    auxd2=g[nhabit][0][fin]+g[nhabit][1][fin]*auxd+
      g[nhabit][2][fin]*auxd*auxd+g[nhabit][3][fin]*auxd*auxd*auxd;
    *gi=auxd1+(wavel-wl[ini])*(auxd2-auxd1)/(wl[fin]-wl[ini]);
    
    auxd1=g[nhabit][4][ini];
    auxd2=g[nhabit][4][fin];
    *g2i=auxd1+(wavel-wl[ini])*(auxd2-auxd1)/(wl[fin]-wl[ini]);
    
    auxd1=g[nhabit][5][ini];
    auxd2=g[nhabit][5][fin];
    *fi=auxd1+(wavel-wl[ini])*(auxd2-auxd1)/(wl[fin]-wl[ini]);
  }
  else{
    ini=0;
    fin=11;
    med=(fin+ini)/2;

    wavel= wavel > 13.0 ? 13.0 : wavel;
    
    while((med!=ini)&&(med!=fin))
      {
	if(wavel<=wlp[med])
	  fin=med;
	else
	  ini=med;
	med=(fin+ini)/2;
      }    

    auxd=1.0/(2.0*reff);
    auxd1=(3.0*1000.0/(0.9167*4.0*reff))*(qep[2][ini]*auxd*auxd+qep[1][ini]*auxd+qep[0][ini]);  
    auxd2=(3.0*1000.0/(0.9167*4.0*reff))*(qep[2][fin]*auxd*auxd+qep[1][fin]*auxd+qep[0][fin]);
    *exti=auxd1+(wavel-wlp[ini])*(auxd2-auxd1)/(wlp[fin]-wlp[ini]);

    auxd=1.0/(2.0*reff);
    auxd1=(ssp[2][ini]*auxd*auxd+ssp[1][ini]*auxd+ssp[0][ini]);  
    auxd2=(ssp[2][fin]*auxd*auxd+ssp[1][fin]*auxd+ssp[0][fin]);
    *ssai=auxd1+(wavel-wlp[ini])*(auxd2-auxd1)/(wlp[fin]-wlp[ini]);
    
    auxd=2.0*reff;
    auxd1=gp[0][ini]*pow(auxd,gp[1][ini]); 
    auxd2=gp[0][fin]*pow(auxd,gp[1][fin]);
    *gi=auxd1+(wavel-wlp[ini])*(auxd2-auxd1)/(wlp[fin]-wlp[ini]);
  
    *g2i=0.5;  /*Arbitrary*/
    *fi=1.0;   /*Only forward*/
  }

  return 0;
}



/*agf*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++
Function yang56.
 Computes the parametrization of single scattering properties
 for 7 habits, for 0.2<wavelength<4.8 um.

Input:  
 wavel -> wavelength
 reff ->  effective radius (could be computed using de2re)
 habit -> name of habit: Solid-Column 0, Hollow-Column 1, Rough-Aggregate 2,
 		         Rosette-4 3, Rosette-6 4, Plate 5 or Dendrite?
			 (case insensitive)

Output:
 *exti -> pointer to interp. volume extinction coefficient per ice water content.
 *ssai -> pointer to interp. single scattering albedo.
 *gi ->   pointer to interp. asymmetry parameter.
 *g2i ->  pointer to interp. 2nd asymmetry parameter (double Henyey-Greenstaein).
 *fi ->   pointer to interp. factor f in double Henyey-Greenstaein function.

++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int yang56 (float wavelength, float reff, int habit, char *path, float *exti, 
	    float *ssai, float *gi, float *g2i, float *fi, int newkey)
{
  static char fileread [IC_HABIT_NN];
  static int first=1;

  static double *a0[IC_HABIT_NN], *a1[IC_HABIT_NN], *a2[IC_HABIT_NN], *a3[IC_HABIT_NN];
  static double *b0[IC_HABIT_NN], *b1[IC_HABIT_NN], *b2[IC_HABIT_NN], *b3[IC_HABIT_NN];
  static double *c0[IC_HABIT_NN], *c1[IC_HABIT_NN], *c2[IC_HABIT_NN], *c3[IC_HABIT_NN];
  static double *g2[IC_HABIT_NN], *f[IC_HABIT_NN];

  static double **wlmin=NULL, **wlmax=NULL;
  static int nwvn=0;

  double minreff=0, maxreff=0;

  char filename[FILENAME_MAX]="";
  char extstr[FILENAME_MAX]="";
 
  int status=0, iwvn=0;


  /* minimum and maximum wavelength of each habit */
  if (first) {
    wlmin = calloc (IC_HABIT_NN, sizeof(double *));
    wlmax = calloc (IC_HABIT_NN, sizeof(double *));
    first=0;
  }

  /* check whether reff is in the right range */
  switch (habit)  {
  case IC_HABIT_SOLID_COLUMN:   
    maxreff=MAX_REFF_SC; 
    minreff=MIN_REFF_SC; 
    strcpy (extstr, "Solid-Column");
    break;

  case IC_HABIT_HOLLOW_COLUMN:  
    maxreff=MAX_REFF_HC; 
    minreff=MIN_REFF_HC; 
    strcpy (extstr, "Hollow-Column");
    break;

  case IC_HABIT_ROUGH_AGGREGATE:     
    maxreff=MAX_REFF_RA; 
    minreff=MIN_REFF_RA; 
    strcpy (extstr, "Rough-Aggregate");
    break;

  case IC_HABIT_ROSETTE_4:      
    maxreff=MAX_REFF_R4; 
    minreff=MIN_REFF_R4; 
    strcpy (extstr, "Rosette-4");
    break;

  case IC_HABIT_ROSETTE_6:      
    maxreff=MAX_REFF_R6; 
    minreff=MIN_REFF_R6; 
    strcpy (extstr, "Rosette-6");
    break;

  case IC_HABIT_PLATE:         
    maxreff=MAX_REFF_PL; 
    minreff=MIN_REFF_PL; 
    strcpy (extstr, "Plate");
    break;

  case IC_HABIT_DENDRITE:         
    maxreff=MAX_REFF_DE; 
    minreff=MIN_REFF_DE; 
    strcpy (extstr, "Dendrite");
    break;

  case IC_HABIT_DROXTAL:         
    maxreff=MAX_REFF_DR; 
    minreff=MIN_REFF_DR; 
    strcpy (extstr, "Droxtal");
    break;

  case IC_HABIT_SPHEROID:         
    maxreff=MAX_REFF_SP; 
    minreff=MIN_REFF_SP; 
    strcpy (extstr, "Spheroid");
    break;

  default:
    fprintf (stderr, "Error, unknown habit %d\n", habit);
    return HABITNF;
  }

  if (reff>maxreff) {
    fprintf (stderr, "Error, effective radius %f larger than maximum\n",reff);
    fprintf (stderr, "allowed value %6.2f for this habit!\n", maxreff);
    return REFFOUT;
  }

  if (reff<minreff) {
    fprintf (stderr, "Error, effective radius %f smaller than minimum\n",reff);
    fprintf (stderr, "allowed value %6.2f for this habit!\n", minreff);
    return REFFOUT;
  }

  /* Read the data files */
  if (!fileread[habit]) {
    fileread[habit]++;

    /*------------------------*/
    /* extinction coefficient */
    /*------------------------*/
    if (newkey)
      sprintf (filename, "%sext.new.%s", path, extstr);
    else 
      sprintf (filename, "%sext.%s", path, extstr);

    status = read_6c_file (filename, 
			   &wlmin[habit], &wlmax[habit], 
			   &a0[habit], &a1[habit], &a2[habit], &a3[habit], &nwvn);

    if (status!=0) {
      fprintf (stderr, "Error %d reading cirrus cloud extinction from %s\n", 
	       status, filename);
      return status;
    }


    /*--------------------------*/
    /* single scattering albedo */
    /*--------------------------*/
    if (newkey)
      sprintf (filename, "%sssa.new.%s", path, extstr);
    else
      sprintf (filename, "%sssa.%s", path, extstr);
      

    status = read_6c_file (filename, 
			   &wlmin[habit], &wlmax[habit], 
			   &b0[habit], &b1[habit], &b2[habit], &b3[habit], &nwvn);

    if (status!=0) {
      fprintf (stderr, "Error %d reading cirrus cloud single scattering albedo from %s\n", 
	       status, filename);
      return status;
    }


    /*---------------------*/
    /* asymmetry parameter */
    /*---------------------*/
    if (newkey)
      sprintf (filename, "%sg.new.%s", path, extstr);
    else
      sprintf (filename, "%sg.%s", path, extstr);
    
    status = read_8c_file (filename, 
			   &wlmin[habit], &wlmax[habit], 
			   &c0[habit], &c1[habit], &c2[habit], &c3[habit], &g2[habit], &f[habit], &nwvn);

    if (status!=0) {
      fprintf (stderr, "Error %d reading cirrus cloud single scattering albedo from %s\n", 
	       status, filename);
      return status;
    }
  }


  /* determine the tabulated wavelength interval containing the actual wavelength */
  if (wavelength < wlmin[habit][0] || wavelength > wlmax[habit][nwvn-1]) {
    fprintf (stderr, "Error, wavelength %f micron not covered by Yang/Key parameterization\n", wavelength);
    fprintf (stderr, "Allowed range is %.1f - %.1f micron\n", wlmin[habit][0], wlmax[habit][nwvn-1]);
    return -1;
  }

  /* determine wavelength index */
  for (iwvn=0; iwvn<nwvn-1; iwvn++)
    if (wavelength < wlmax[habit][iwvn])
      break;
  
  *exti = a0[habit][iwvn] + a1[habit][iwvn]/reff + a2[habit][iwvn]/reff/reff + a3[habit][iwvn]/reff/reff/reff;
  *ssai = b0[habit][iwvn] + b1[habit][iwvn]*reff + b2[habit][iwvn]*reff*reff + b3[habit][iwvn]*reff*reff*reff;
  *gi   = c0[habit][iwvn] + c1[habit][iwvn]*reff + c2[habit][iwvn]*reff*reff + c3[habit][iwvn]*reff*reff*reff;
  *g2i  = g2[habit][iwvn];
  *fi   = f [habit][iwvn];


  /* For safety reasons check if asymmetry parameter is within -0.999999 to +0.9999999. */
  /* Josef Gasteiger found that around 2.8 micron g might be larger than 1 for rough    */
  /* aggregates close to the upper limit of the size range (104.3 micron); 07.07.2012   */             

  if (*gi > 0.999999) {
    *gi = 0.999999;
    fprintf (stderr, "Warning, yang56() asymmetry parameter out of range; setting to reasonable value\n");
  }

  if (*g2i > 0.999999) {
    *g2i = 0.999999;
    fprintf (stderr, "Warning, yang56() asymmetry parameter out of range; setting to reasonable value\n");
  }

  if (*gi < -0.999999) {
    *gi = -0.999999;
    fprintf (stderr, "Warning, yang56() asymmetry parameter out of range; setting to reasonable value\n");
  }

  if (*g2i < -0.999999) {
    *g2i = -0.999999;
    fprintf (stderr, "Warning, yang56() asymmetry parameter out of range; setting to reasonable value\n");
  }



  return 0;
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++
Function dm2area.
Computes the projected area of an ice particle from its maximum
dimension and habit.

Input:  

dm -> maximum dimension
habit -> name of habit: Solid-Column, Hollow-Column, Rough-Aggregate,
 			Rosette-4, Rosette-6, Plate or Dendrite. 
			(case insensitive)

Output: projected area.


++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

static double dm2area (double dm, char *habit)
{
  double a[IC_HABIT_NN][5]=
	{{0.33401E+0, 0.36477E+0, 0.30855E+0, -0.55631E-1, 0.30162E-2},
	 {0.33401E+0, 0.36477E+0, 0.30855E+0, -0.55631E-1, 0.30162E-2},
	 {-0.47737E+0, 0.10026E+1, -0.10030E-2, 0.15166E-3, -0.78433E-5},
	 {0.15909E+0, 0.84308E+0, 0.70161E-2, -0.11003E-2, 0.45161E-4},
	 {0.14195E+0, 0.84394E+0, 0.72125E-2, -0.11219E-2, 0.45819E-4},
	 {0.43773E+0, 0.75497E+0, 0.19033E-1, 0.35191E-3, -0.70782E-4},
	 {0.43773E+0, 0.75497E+0, 0.19033E-1, 0.35191E-3, -0.70782E-4}/*!!!*/ 
	};
  double area=0, auxd=0;
  long i=0, nhabit=0;
  char habit_class[IC_HABIT_NN][80];

  sprintf(habit_class[IC_HABIT_SOLID_COLUMN ],"%s","Solid-Column");
  sprintf(habit_class[IC_HABIT_HOLLOW_COLUMN],"%s","Hollow-Column");
  sprintf(habit_class[IC_HABIT_ROUGH_AGGREGATE],   "%s","Rough-Aggregate");
  sprintf(habit_class[IC_HABIT_ROSETTE_4],    "%s","Rosette-4");
  sprintf(habit_class[IC_HABIT_ROSETTE_6],    "%s","Rosette-6");
  sprintf(habit_class[IC_HABIT_PLATE],       "%s","Plate");
  sprintf(habit_class[IC_HABIT_DENDRITE],    "%s","Dendrite");

  for(i=0;i<IC_HABIT_NN;i++)
    if(strcasecmp(habit,habit_class[i])==0)
    {
      nhabit=i;
      break;
    }
  if(i==IC_HABIT_NN)
  {
    printf("Habit '%s' not found\n",habit);
    return -1.0;
  }

  auxd=1.0;
  area=0.0;
  for(i=0;i<5;i++)
  {
    area+=a[nhabit][i]*auxd;
    auxd*=log(dm);
  }
  area=exp(area);  /*Diameter of an equivalent sphere*/

  area=(area*area/4.0)*3.1416;
  
  return area;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++
Function dm2vol.
Computes the volume of an ice particle from its maximum
dimension and habit.

Input:  

dm -> maximum dimension
habit -> name of habit: Solid-Column, Hollow-Column, Rough-Aggregate,
 			Rosette-4, Rosette-6, Plate or Dendrite. 
			(case insensitive)

Output: volume.


++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

static double dm2vol(double dm, char *habit)
{
  double a[IC_HABIT_NN][5]=
	{{0.30581E+0, 0.26252E+0, 0.35458E+0, -0.63202E-1, 0.33755E-2},
	 {0.24568E+0, 0.26202E+0, 0.35479E+0, -0.63236E-1, 0.33773E-2},
	 {-0.70160E+0, 0.99215E+0, 0.29322E-2, -0.40492E-3, 0.18841E-4},
	 {-0.97940E-1, 0.85683E+0, 0.29483E-2, -0.14341E-2, 0.74627E-4},
	 {-0.10318E+0, 0.86290E+0, 0.70665E-3, -0.11055E-2, 0.57906E-4},
	 {0.31228E+0, 0.80874E+0, 0.29287E-2, -0.44378E-3, 0.23109E-4},
	 {0.31228E+0, 0.80874E+0, 0.29287E-2, -0.44378E-3, 0.23109E-4} /*!!!*/
	};
  double vol=0, auxd=0;
  long i=0, nhabit=0;
  char habit_class[IC_HABIT_NN][80];

  sprintf(habit_class[IC_HABIT_SOLID_COLUMN ],"%s","Solid-Column");
  sprintf(habit_class[IC_HABIT_HOLLOW_COLUMN],"%s","Hollow-Column");
  sprintf(habit_class[IC_HABIT_ROUGH_AGGREGATE],   "%s","Rough-Aggregate");
  sprintf(habit_class[IC_HABIT_ROSETTE_4],    "%s","Rosette-4");
  sprintf(habit_class[IC_HABIT_ROSETTE_6],    "%s","Rosette-6");
  sprintf(habit_class[IC_HABIT_PLATE],       "%s","Plate");
  sprintf(habit_class[IC_HABIT_DENDRITE],    "%s","Dendrite");

  for(i=0;i<IC_HABIT_NN;i++)
    if(strcasecmp(habit,habit_class[i])==0)
    {
      nhabit=i;
      break;
    }
  if(i==IC_HABIT_NN)
  {
    printf("Habit '%s' not found\n",habit);
    return -1.0;
  }

  auxd=1.0;
  vol=0.0;
  for(i=0;i<5;i++)
  {
    vol+=a[nhabit][i]*auxd;
    auxd*=log(dm);
  }
  vol=exp(vol);  /*Diameter of an equivalent sphere*/

  vol=vol*vol*vol*3.1416/6.0;
  
  return vol;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++
Function dm2re.
Computes the "effective radius" of an ice particle from
its maximum dimension and habit.

Input:  

dm -> maximum dimension
habit -> name of habit: Solid-Column, Hollow-Column, Rough-Aggregate,
 			Rosette-4, Rosette-6, Plate or Dendrite. 
			(case insensitive)

Output: effective radius.

++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

double dm2re (double dm, char *habit)
{
  double re = (3.0/4.0) * dm2vol (dm, habit) / dm2area (dm, habit);
  return re;
}
/*agf*/
