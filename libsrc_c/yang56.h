/************************************************************************
 * $Id: yang56.h 3092 2015-03-31 09:30:57Z Claudia.Emde $
 ************************************************************************/

#ifndef __yang56_h
#define __yang56_h


#define IC_HABIT_NN           15  /* Number of habits                         */
#define IC_HABIT_SOLID_COLUMN  0 
#define IC_HABIT_HOLLOW_COLUMN 1
#define IC_HABIT_ROUGH_AGGREGATE    2 
#define IC_HABIT_ROSETTE_4     3
#define IC_HABIT_ROSETTE_6     4
#define IC_HABIT_PLATE        5
#define IC_HABIT_DENDRITE     6
#define IC_HABIT_DROXTAL      7
#define IC_HABIT_SPHEROID     8
#define IC_HABIT_GHM          9
/* further habits for Yang2013 */
#define IC_HABIT_COLUMN_8ELEMENTS      10
#define IC_HABIT_HOLLOW_BULLET_ROSETTE 11
#define IC_HABIT_PLATE_10ELEMENTS      12
#define IC_HABIT_PLATE_5ELEMENTS       13
#define IC_HABIT_SOLID_BULLET_ROSETTE  14

/* roughness values for Yang2013 */
#define IC_ROUGHNESS_NN       3
#define IC_ROUGHNESS_SMOOTH   0
#define IC_ROUGHNESS_MODERATE 1
#define IC_ROUGHNESS_SEVERE   2


#define MAX_REFF_SC  84.22   /* Maximum effective radius for solid column    */
#define MAX_REFF_HC  70.24   /* Maximum effective radius for hollow column   */
#define MAX_REFF_RA 108.10   /* Maximum effective radius for rough aggregate */
#define MAX_REFF_R4  45.30   /* Maximum effective radius for rosette 4       */
#define MAX_REFF_R6  46.01   /* Maximum effective radius for rosette 6       */
#define MAX_REFF_PL  48.18   /* Maximum effective radius for plate           */
#define MAX_REFF_DE   1.88   /* Maximum effective radius for dendrite        */
#define MAX_REFF_DR 293.32   /* Maximum effective radius for droxtal         */
#define MAX_REFF_SP 203.39   /* Maximum effective radius for spheroid        */

#define MIN_REFF_SC   5.96   /* Minimum effective radius for solid column    */
#define MIN_REFF_HC   4.97   /* Minimum effective radius for hollow column   */
#define MIN_REFF_RA   3.55   /* Minimum effective radius for rough aggregate */
#define MIN_REFF_R4   2.77   /* Minimum effective radius for rosette 4       */
#define MIN_REFF_R6   2.85   /* Minimum effective radius for rosette 6       */
#define MIN_REFF_PL   4.87   /* Minimum effective radius for plate           */
#define MIN_REFF_DE   0.45   /* Minimum effective radius for dendroid        */
#define MIN_REFF_DR   9.48   /* Minimum effective radius for droxtal         */
#define MIN_REFF_SP   6.58   /* Minimum effective radius for spheroid        */

int yang56 (float wavel, float reff, int nhabit, char *path, float *exti, 
	    float *ssai, float *gi, float *g2i, float *fi, int newkey);

double dm2re (double dm, char *habit);

#endif
