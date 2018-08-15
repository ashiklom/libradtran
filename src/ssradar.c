/*--------------------------------------------------------------------
 * $Id: ssradar.c 2623 2013-XX-XX XX:XX:XXZ christian.pause $
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

/* TODO: Add Attenuation by hydrometeors
         Check Possibility of Mie Correction
         Marshall-Palmer distributions
         Add an atmosphere */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errors.h>
#if HAVE_LIBNETCDF
#include <netcdf.h>
#endif
#include "ascii.h"
#include "integrat.h"

#ifndef PI
#define PI 3.1415926535897932
#endif

#define DIST_NONE       0
#define DIST_GAMMA      1
#define DIST_LOGNORMAL  2
#define DIST_FILE       3
#define DIST_AER        4

typedef struct {
  int      kindOfDistribution;
  char*    sizeDistributionFilename;
  double   height;
  double   thickness;
  double   effectiveRadius; 
  double   distributionParameterAlpha;
  double   distributionStep;
  double   distributionMaxRadiusFactor;
  double   reflectivityFactor;
  double*  dropletDistribution;
  double*  radiusGrid;
  int      numberOfGridPoints;
  double   temperature;
  int      waterPhase;
  double   LWC;
} distributionParameterStruct;

typedef struct {
  int                           numberOfReflectiveLayers;
  distributionParameterStruct** allDistributionParameters;
  double                        radarZPosition;
  double                        cosZenithAngleUmu;
  int                           numberOfRangeGates;
  double                        rangeGateWidth;
  double                        firstRangeGate;
  double                        wavelengthInMicrometer; 
  int                           output;                     
  double                        groundAltitude;
  double                        toaAltitude;
} ssRadarStruct;


int singleScatteringRadar ( ssRadarStruct* ssRadar, char* inputFilename, char* outputFilename );
int readSSRadarFile ( ssRadarStruct *ssRadar, char* inputFilename);
int findEmptyLayers(double height0, double height1, double thickness);
double floatZeros (float f, int acc);
double linearInterpolation(double a, double b, double c, double d, double e);
int writeSSRadarOutput(double* rangeGatesVertical,
                       double* rangeGates,
                       double* rangeGateReflecticityFactor,
                       double  radarZPosition,
                       double  rangeGateWidth,
                       double  wavelengthInMicrometer,
                       double  cosZenithAngleUmu,
                       char*   outputFilename,
                       int     numberOfRangeGates,
                       int     output);
void exitSSRadar(distributionParameterStruct*** allDistributionParameters,
                double** rangeGates,
                double** rangeGatesVertical,
                double** rangeGateReflectivityFactor,
                double** reflectiveLayersHeights, 
                double** reflectiveLayersThicknesses,
                double** reflectiveLayersReflectivities,
                double** allLayersHeights, 
                double** allLayersThicknesses,
                double** allLayersReflectivities,
                int** indexArray,
                int numberOfReflectiveLayers);
double getWaterDensity(double temperature, int phase);
double gamma_ln(float xx);
double double_max(double x, double y);
int calc_distribution_radius_grid(/* Input */
                                  int distribution,
                                  char* sd_filename,
                                  double r_eff_min,
                                  double r_eff_max,
                                  double lambda,      /*wavelength in micrometer*/
                                  double dx_max,
                                  double n_r_max,
                                  int verbose,
                                  /* Output */
                                  int* n_dens,
                                  double** radius,
                                  double** number_dens);
int calc_distribution_number_dens (/* Input */
                                   int distribution,
                                   double alpha,
                                   double r_eff,
                                   double rho_medium,
                                   double LWC,
                                   int n_dens,
                                   double* radius,
                                   int verbose,
                                   /* Output */
                                   double* number_dens);

void printArrayDouble(double *A, int n, char* C);
void printArrayInt(int *A, int n, char* C);
int locate2 (double *xx, int n, double x, double umu);

/***********************************************************************************/
/* Function: main                                                         @62_30i@ */
/* Description: The main function obviously.                                       */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Christian Pause                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/
int main( int argc, const char* argv[] ) {
  char* inputFilename = NULL;
  char* outputFilename = NULL;
  inputFilename = calloc(256,sizeof(char));
  outputFilename = calloc(256,sizeof(char));
  ssRadarStruct* ssRadar = NULL;
  ssRadar = calloc(1,sizeof(ssRadarStruct));
  if (argc > 1) {
    if (!strcmp(argv[1],"-h")) {
      fprintf(stderr,"===============================================\n");
      fprintf(stderr,"Usage: ./ssradar inputfilename [outputfilename]\n");
      fprintf(stderr,"Outputfilename is optional.\n");
      fprintf(stderr,"See libradtran.pdf for details.\n");
      fprintf(stderr,"===============================================\n");
    }
    else 
      strcpy(inputFilename,argv[1]);
  }
  if (argc > 2)
    strcpy(outputFilename,argv[2]);
  if (argc < 2 || argc > 3) {
    fprintf(stderr,"===============================================\n");
    fprintf(stderr,"Usage: ./ssradar inputfilename [outputfilename]\n");
    fprintf(stderr,"Outputfilename is optional.\n");
    fprintf(stderr,"See libradtran.pdf for details.\n");
    fprintf(stderr,"===============================================\n");
  }
  else if (singleScatteringRadar(ssRadar,inputFilename,outputFilename))
    fprintf(stderr,"Error in ssradar.\n");
  free(ssRadar);
  free(inputFilename);
  free(outputFilename);
}

/***********************************************************************************/
/* Function: singleScatteringRadar                                        @62_30i@ */
/* Description: Simple radar simulator, calculates reflectivity directly from      */
/*              droplet distribution.                                              */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Christian Pause                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/
int singleScatteringRadar ( ssRadarStruct* ssRadar, char* inputFilename, char* outputFilename )
{
  int i,j;
  int addIndex = 0;
  double* rangeGates = NULL;
  double* rangeGatesVertical = NULL;
  double* rangeGateReflectivityFactor = NULL;
  int verbose = 0;
  double waterDensity;
  double* reflectiveLayersHeights = NULL;
  double* reflectiveLayersThicknesses = NULL;
  double* reflectiveLayersReflectivities = NULL;
  double* allLayersHeights = NULL;
  double* allLayersThicknesses = NULL;
  double* allLayersReflectivities = NULL;
  int numberOfAllLayers = 0;
  double umuFactor = 0.;
  double lowerBorder = 0.;
  double upperBorder = 0.;
  int umuIndex = 0;
  int umuIndexR = 0;
  double temp = 0.;
  int counter = 0;
  int* indexArray = NULL;

  if (verbose) {  
    fprintf(stderr,"Welcome to Single Scattering Radar!\n");
    fprintf(stderr,"This program calculates the pure Rayleigh scattering reflectivity factor\nin an arbitrary plan-parallel atmospheric environment.\n");
  }

  /* read in Data */
  if (verbose)
    fprintf(stderr,"Reading data file.\n");
  if (readSSRadarFile ( ssRadar, inputFilename ) ) {
    fprintf(stderr,"Error in readSSRadarFile\n");
    exitSSRadar(&ssRadar->allDistributionParameters, &rangeGates, &rangeGatesVertical, &rangeGateReflectivityFactor,
                &reflectiveLayersHeights, &reflectiveLayersThicknesses, &reflectiveLayersReflectivities,
                &allLayersHeights, &allLayersThicknesses, &allLayersReflectivities,
                &indexArray, ssRadar->numberOfReflectiveLayers);
    return -1;
  }

  /* allocate stuff */
  rangeGates = calloc((ssRadar->numberOfRangeGates+1), sizeof(double));
  rangeGatesVertical = calloc((ssRadar->numberOfRangeGates+1), sizeof(double));
  rangeGateReflectivityFactor = calloc((ssRadar->numberOfRangeGates+1), sizeof(double));
  reflectiveLayersHeights = calloc(ssRadar->numberOfReflectiveLayers,sizeof(double));
  reflectiveLayersThicknesses = calloc(ssRadar->numberOfReflectiveLayers,sizeof(double));
  reflectiveLayersReflectivities = calloc(ssRadar->numberOfReflectiveLayers,sizeof(double));
  indexArray = calloc(ssRadar->numberOfRangeGates+1,sizeof(int));  

  /* set auxillary variables related to umu */
  if (ssRadar->cosZenithAngleUmu > 0.) {
    umuFactor = 1.;
    umuIndex = 0;
    umuIndexR = ssRadar->numberOfRangeGates;
  }
  else {
    umuFactor = -1.;
    umuIndex = ssRadar->numberOfRangeGates;
    umuIndexR = 0;
  }

  if (verbose)
    fprintf(stderr,"Calculating droplet distributions and reflectivities.\n");
  /* calculate droplet distribution and reflectivity factor for each reflective layer */
  for (i=0;i<ssRadar->numberOfReflectiveLayers;i++) {
    if ( calc_distribution_radius_grid( ssRadar->allDistributionParameters[i]->kindOfDistribution,
                                   ssRadar->allDistributionParameters[i]->sizeDistributionFilename,
                                   ssRadar->allDistributionParameters[i]->effectiveRadius, 
                                   ssRadar->allDistributionParameters[i]->effectiveRadius, 
                                   ssRadar->wavelengthInMicrometer,
                                   ssRadar->allDistributionParameters[i]->distributionStep,
                                   ssRadar->allDistributionParameters[i]->distributionMaxRadiusFactor,
                                   verbose,
                                   &ssRadar->allDistributionParameters[i]->numberOfGridPoints,
                                   &ssRadar->allDistributionParameters[i]->radiusGrid,
                                   &ssRadar->allDistributionParameters[i]->dropletDistribution) ) {
      exitSSRadar(&ssRadar->allDistributionParameters, &rangeGates, &rangeGatesVertical, &rangeGateReflectivityFactor,
                  &reflectiveLayersHeights, &reflectiveLayersThicknesses, &reflectiveLayersReflectivities,
                  &allLayersHeights, &allLayersThicknesses, &allLayersReflectivities,
                  &indexArray, ssRadar->numberOfReflectiveLayers);
      fprintf(stderr,"Error in calc_distribution_radius_grid\n");
      return -1;
    }
    if ((waterDensity = getWaterDensity(ssRadar->allDistributionParameters[i]->temperature,
                                        ssRadar->allDistributionParameters[i]->waterPhase)) == -9999.) {
      exitSSRadar(&ssRadar->allDistributionParameters, &rangeGates, &rangeGatesVertical, &rangeGateReflectivityFactor,
                  &reflectiveLayersHeights, &reflectiveLayersThicknesses, &reflectiveLayersReflectivities,
                  &allLayersHeights, &allLayersThicknesses, &allLayersReflectivities,
                  &indexArray, ssRadar->numberOfReflectiveLayers);
      fprintf(stderr,"Errorneous temperature in layer %d\n",i+1);
      return -1;
    }
    if ( calc_distribution_number_dens( ssRadar->allDistributionParameters[i]->kindOfDistribution,
                                   ssRadar->allDistributionParameters[i]->distributionParameterAlpha,
                                   ssRadar->allDistributionParameters[i]->effectiveRadius,
                                   waterDensity,
                                   ssRadar->allDistributionParameters[i]->LWC,
                                   ssRadar->allDistributionParameters[i]->numberOfGridPoints, 
                                   ssRadar->allDistributionParameters[i]->radiusGrid,
                                   verbose,
                                   ssRadar->allDistributionParameters[i]->dropletDistribution) ) {
      exitSSRadar(&ssRadar->allDistributionParameters, &rangeGates, &rangeGatesVertical, &rangeGateReflectivityFactor,
                  &reflectiveLayersHeights, &reflectiveLayersThicknesses, &reflectiveLayersReflectivities,
                  &allLayersHeights, &allLayersThicknesses, &allLayersReflectivities,
                  &indexArray, ssRadar->numberOfReflectiveLayers);
      fprintf(stderr,"Error in calc_distribution_number_dens\n");
      return -1;
    }
    for (j=0;j<ssRadar->allDistributionParameters[i]->numberOfGridPoints;j++) {
      ssRadar->allDistributionParameters[i]->reflectivityFactor += ssRadar->allDistributionParameters[i]->dropletDistribution[j]
        * pow(2.*ssRadar->allDistributionParameters[i]->radiusGrid[j], 6.);

    }
  }

  /* Create range gates. */
  rangeGates[0] = ssRadar->firstRangeGate + ssRadar->radarZPosition;
  rangeGatesVertical[0] = ssRadar->radarZPosition + umuFactor*(ssRadar->firstRangeGate - ssRadar->radarZPosition) * ssRadar->cosZenithAngleUmu;
  for(i=1;i<ssRadar->numberOfRangeGates;i++) {
    rangeGates[i] = rangeGates[i-1] + umuFactor * ssRadar->rangeGateWidth;
    rangeGatesVertical[i] = rangeGatesVertical[i-1] + ssRadar->rangeGateWidth * ssRadar->cosZenithAngleUmu;
  }
  rangeGates[ssRadar->numberOfRangeGates] = umuFactor * ssRadar->radarZPosition/ssRadar->cosZenithAngleUmu - rangeGates[ssRadar->numberOfRangeGates-1];
  rangeGatesVertical[ssRadar->numberOfRangeGates] = ssRadar->radarZPosition - rangeGatesVertical[ssRadar->numberOfRangeGates-1];

  /* find total number of layers, i.e. empty and reflective layers between ground 
     and the highest reflective layer. An additional layers is added on top and/or
     below the ground if needed to cover all range gates.  */

  if (verbose)
    fprintf(stderr,"Creating layers and range gates.\n");
    /* Create range gates. Last index is the end of the last range gate. */
  rangeGates[0] = ssRadar->firstRangeGate;
  rangeGatesVertical[0] = ssRadar->radarZPosition + ssRadar->firstRangeGate * ssRadar->cosZenithAngleUmu;
  for(i=1;i<ssRadar->numberOfRangeGates+1;i++) {
    rangeGates[i] = rangeGates[i-1] + ssRadar->rangeGateWidth;
    rangeGatesVertical[i] = rangeGatesVertical[i-1] + ssRadar->rangeGateWidth * ssRadar->cosZenithAngleUmu;
  }

  /* copy data from main struct into handier arrays */
  for(i=0;i<ssRadar->numberOfReflectiveLayers;i++) {
    reflectiveLayersHeights[i] = ssRadar->allDistributionParameters[i]->height;
    reflectiveLayersThicknesses[i] = ssRadar->allDistributionParameters[i]->thickness;
    reflectiveLayersReflectivities[i] = ssRadar->allDistributionParameters[i]->reflectivityFactor;
  }

  /* determine upper and lower border */
  temp = (rangeGatesVertical[umuIndex] < ssRadar->groundAltitude) ? rangeGatesVertical[umuIndex] : ssRadar->groundAltitude;
  lowerBorder = (temp < reflectiveLayersHeights[0]) ? temp : reflectiveLayersHeights[0];
  upperBorder = (rangeGatesVertical[umuIndexR] > reflectiveLayersHeights[ssRadar->numberOfReflectiveLayers-1] + reflectiveLayersThicknesses[ssRadar->numberOfReflectiveLayers-1]) ? rangeGatesVertical[umuIndexR] : reflectiveLayersHeights[ssRadar->numberOfReflectiveLayers-1] + reflectiveLayersThicknesses[ssRadar->numberOfReflectiveLayers-1];

  /* count number of all layers */
  if (lowerBorder < reflectiveLayersHeights[0]) 
    numberOfAllLayers++;
  if (lowerBorder < ssRadar->groundAltitude && ssRadar->groundAltitude < reflectiveLayersHeights[0])
    numberOfAllLayers++;
  numberOfAllLayers += ssRadar->numberOfReflectiveLayers;
  for(i=0;i<ssRadar->numberOfReflectiveLayers-1;i++)
    numberOfAllLayers += findEmptyLayers(reflectiveLayersHeights[i],reflectiveLayersHeights[i+1],reflectiveLayersThicknesses[i]);
  if (upperBorder > reflectiveLayersHeights[ssRadar->numberOfReflectiveLayers-1] + reflectiveLayersThicknesses[ssRadar->numberOfReflectiveLayers-1])
    numberOfAllLayers++;

  /* allocate arrays for all layers */
  allLayersHeights = calloc(numberOfAllLayers+1,sizeof(double));
  allLayersThicknesses = calloc(numberOfAllLayers,sizeof(double));
  /* CE: increased length of allLayersReflectivities by 1, because of memory error below. Not sure whether this is correct ... */
  allLayersReflectivities = calloc(numberOfAllLayers+1,sizeof(double));

  /* fill arrays for all layers */
  if (lowerBorder < reflectiveLayersHeights[0]) {
    if (lowerBorder < ssRadar->groundAltitude && ssRadar->groundAltitude < reflectiveLayersHeights[0]) {
      allLayersHeights[0] = lowerBorder;
      allLayersThicknesses[0] = ssRadar->groundAltitude - lowerBorder;
      allLayersReflectivities[0] = 1.e+28;
      allLayersHeights[1] = ssRadar->groundAltitude;
      allLayersThicknesses[1] = reflectiveLayersHeights[0]-ssRadar->groundAltitude;
      allLayersReflectivities[1] = 0.;
      counter = 2;
    }
    else {
      allLayersHeights[0] = lowerBorder;
      allLayersThicknesses[0] = reflectiveLayersHeights[0]-lowerBorder;
      allLayersReflectivities[0] = 0.;
      counter = 1;
    }
  }
  i=counter;
  for(j=0;j<ssRadar->numberOfReflectiveLayers;j++) {
    allLayersHeights[i] = reflectiveLayersHeights[j];
    allLayersThicknesses[i] = reflectiveLayersThicknesses[j];
    allLayersReflectivities[i] = reflectiveLayersReflectivities[j];
    i++;
    if(j < ssRadar->numberOfReflectiveLayers-1 && findEmptyLayers(reflectiveLayersHeights[j],reflectiveLayersHeights[j+1],reflectiveLayersThicknesses[j])) {
      allLayersHeights[i] = reflectiveLayersHeights[j] + reflectiveLayersThicknesses[j];
      allLayersThicknesses[i] = reflectiveLayersHeights[j+1] - allLayersHeights[i];
      allLayersReflectivities[i] = 0.;
      i++;
    }
    counter = i;
  }

  if (upperBorder > reflectiveLayersHeights[ssRadar->numberOfReflectiveLayers-1] + reflectiveLayersThicknesses[ssRadar->numberOfReflectiveLayers-1]) {
    allLayersHeights[counter] = reflectiveLayersHeights[ssRadar->numberOfReflectiveLayers-1] + reflectiveLayersThicknesses[ssRadar->numberOfReflectiveLayers-1];
    allLayersThicknesses[counter] = upperBorder - allLayersHeights[counter];
    allLayersReflectivities[counter] = 0.;
  }

  allLayersHeights[numberOfAllLayers] = allLayersHeights[numberOfAllLayers-1] + allLayersThicknesses[numberOfAllLayers-1]; 

  /* fill index array which contains the index of the respective layer each range gates begins in */
  for(i=0;i<ssRadar->numberOfRangeGates+1;i++) 
    indexArray[i] = locate2(allLayersHeights,numberOfAllLayers+1,rangeGatesVertical[i],umuFactor);

  /* calculate reflectivity factor in each range gate, interpolation is performed by dividing 
     the fraction of the range falling into one layer by the length of the range gate */
  if (ssRadar->cosZenithAngleUmu > 0.) {
    for(i=0;i<ssRadar->numberOfRangeGates;i++) {
      if (indexArray[i] < numberOfAllLayers) {
        if (indexArray[i+1] == indexArray[i]) {
          rangeGateReflectivityFactor[i] = allLayersReflectivities[indexArray[i]];
        }
        else if (indexArray[i+1] > indexArray[i]) {
          addIndex = 0;
          rangeGateReflectivityFactor[i] = ((allLayersHeights[indexArray[i]+1] - rangeGatesVertical[i])
                                         /   ssRadar->rangeGateWidth) * allLayersReflectivities[indexArray[i]];
          if (indexArray[i+1] > indexArray[i]+1) {
            for(j=1;j<indexArray[i+1]-indexArray[i];j++) {
              rangeGateReflectivityFactor[i] += ((allLayersThicknesses[indexArray[i]+j])
                                              /   ssRadar->rangeGateWidth) * allLayersReflectivities[indexArray[i]+j];
              addIndex = j;
            }
          }
	  rangeGateReflectivityFactor[i] += ((rangeGatesVertical[i] + ssRadar->rangeGateWidth - allLayersHeights[indexArray[i]+addIndex+1] )
                                          /   ssRadar->rangeGateWidth) * allLayersReflectivities[indexArray[i]+addIndex+1];
        }
      }
    } 
  }
  else {
    for(i=0;i<ssRadar->numberOfRangeGates;i++) {
      if (indexArray[i] < numberOfAllLayers) {
        if (indexArray[i+1] == indexArray[i])
          rangeGateReflectivityFactor[i] = allLayersReflectivities[indexArray[i]];
        else if (indexArray[i+1] < indexArray[i]) {
          addIndex = 0;
          rangeGateReflectivityFactor[i] = ((- allLayersHeights[indexArray[i]] + rangeGatesVertical[i])
                                         /   ssRadar->rangeGateWidth) * allLayersReflectivities[indexArray[i]];
          if (indexArray[i+1] < indexArray[i]+1) {
            for(j=1;j<indexArray[i+1]-indexArray[i];j++) {
              rangeGateReflectivityFactor[i] += ((allLayersThicknesses[indexArray[i]-j-1])
                                              /   ssRadar->rangeGateWidth) * allLayersReflectivities[indexArray[i]-j];
              addIndex = j;
            }
          }
          rangeGateReflectivityFactor[i] += ((- rangeGatesVertical[i+1] + allLayersHeights[indexArray[i]-addIndex] )
                                          /   ssRadar->rangeGateWidth) * allLayersReflectivities[indexArray[i]-1-addIndex];
        }
      }
    } 
  }

  if (verbose) {
    printArrayDouble(allLayersHeights,numberOfAllLayers+1,"All Layers Heights");
    printArrayDouble(allLayersThicknesses,numberOfAllLayers,"All Layers Thicknesses");
    printArrayDouble(allLayersReflectivities,numberOfAllLayers,"All Layers Reflectivities");
    printArrayInt(indexArray,ssRadar->numberOfRangeGates+1,"Index Array");
    printArrayDouble(rangeGates,ssRadar->numberOfRangeGates+1,"Range Gates");
    printArrayDouble(rangeGatesVertical,ssRadar->numberOfRangeGates+1,"Range Gates Vertical");
    printArrayDouble(rangeGateReflectivityFactor,ssRadar->numberOfRangeGates,"Range Gate Reflectivities");
  }

  /* Logarithmitize RRF, set signal from ground to 100 dBZ */
  for (i=0;i<ssRadar->numberOfRangeGates;i++) {
    if (rangeGateReflectivityFactor[i] > 0.) 
      rangeGateReflectivityFactor[i] =  10.*log10(rangeGateReflectivityFactor[i]/pow(1000.,6.));
    else
      rangeGateReflectivityFactor[i] = 0.;
  }

  if (verbose)
    fprintf(stderr,"Writing output.\n");
  /* Write output files. */
  if (writeSSRadarOutput(rangeGatesVertical,
                         rangeGates,
                         rangeGateReflectivityFactor,
                         ssRadar->radarZPosition,
                         ssRadar->rangeGateWidth,
                         ssRadar->wavelengthInMicrometer,
                         ssRadar->cosZenithAngleUmu,
                         outputFilename,
                         ssRadar->numberOfRangeGates,
                         ssRadar->output)) {
    exitSSRadar(&ssRadar->allDistributionParameters, &rangeGates, &rangeGatesVertical, &rangeGateReflectivityFactor,
                &reflectiveLayersHeights, &reflectiveLayersThicknesses, &reflectiveLayersReflectivities,
                &allLayersHeights, &allLayersThicknesses, &allLayersReflectivities,
                &indexArray, ssRadar->numberOfReflectiveLayers);
    fprintf(stderr,"Error in writeSSRadarOutput.\n");
    return -1;
  }
  /* Done. */
  exitSSRadar(&ssRadar->allDistributionParameters, &rangeGates, &rangeGatesVertical, &rangeGateReflectivityFactor,
              &reflectiveLayersHeights, &reflectiveLayersThicknesses, &reflectiveLayersReflectivities,
              &allLayersHeights, &allLayersThicknesses, &allLayersReflectivities,
              &indexArray, ssRadar->numberOfReflectiveLayers);
  if (verbose)
    fprintf(stderr,"All done.\n");
  return 0;
}

/***********************************************************************************/
/* Function: readSSRadarFile                                              @62_30i@ */
/* Description: Reads a SSRadar file                                               */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Christian Pause                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/
int readSSRadarFile ( ssRadarStruct *ssRadar, char* inputFilename) {

  int i;
  int minColumns=0, maxColumns=0, rows=0, maxLength=0;
  char*** string = NULL;

  if (ASCII_checkfile ( inputFilename, &rows, &minColumns, &maxColumns, &maxLength ) ) {
    fprintf(stderr,"Error in ASCII_checkfile.\n");
    return -1;
  }
  if (rows == 0) {
    fprintf(stderr,"Error, no reflective layers in ssradar input file!\n");
    return -1;
  }
  ssRadar->numberOfReflectiveLayers = rows-2;
  if (ASCII_calloc_string ( &string, rows, maxColumns, maxLength ) ) {
    fprintf(stderr,"Error in ASCII_calloc_string.\n");
    return -1;
  }
  if (ASCII_readfile ( inputFilename, string ) ) {
    fprintf(stderr,"Error in ASCII_readfile.\n");
    return -1;
  }

  ssRadar->allDistributionParameters = calloc(ssRadar->numberOfReflectiveLayers, sizeof(distributionParameterStruct*));
  for(i=0;i<rows-2;i++) {
    ssRadar->allDistributionParameters[i] = calloc(1, sizeof(distributionParameterStruct));
  }

  /* Format of strings: string[row][column][char] */
  /* Format of input file:
     first line:
     format: wvl in mm, z-pos of radar, zenith angle in degrees, altitude
             d          d               d                        d       
     second line (output: -1 =  text, 1 = mmclx, 0 = both):
     format: No of range gates, width of range gates, first range gate, output
             i                  d                     d                 i       
     following lines:
     format: height, thickness, reff, disttype, alpha, max. radius factor, temp (C), phase, LWC (g/km^3)
             d       d          d     str       d      d                   d         str    d            */
 
  ssRadar->wavelengthInMicrometer = 1000.*atof(string[0][0]);
  ssRadar->radarZPosition = atof(string[0][1]);
  ssRadar->cosZenithAngleUmu = cos(PI*atof(string[0][2])/180.);
  if (ssRadar->cosZenithAngleUmu == 0.) {
    fprintf(stderr,"Error in SSRadar, zenith angle must not be 90 degrees!\n");
    return -1;
  }
  ssRadar->groundAltitude = atof(string[0][3]);
  ssRadar->numberOfRangeGates = atoi(string[1][0]);
  ssRadar->rangeGateWidth = atof(string[1][1]);
  ssRadar->firstRangeGate = atof(string[1][2]);
  ssRadar->output = atoi(string[1][3]);
  for (i=2;i<rows;i++) {
    ssRadar->allDistributionParameters[i-2]->height = atof(string[i][0]);
    ssRadar->allDistributionParameters[i-2]->thickness = atof(string[i][1]);
    ssRadar->allDistributionParameters[i-2]->effectiveRadius = atof(string[i][2]);
    ssRadar->allDistributionParameters[i-2]->distributionParameterAlpha = atof(string[i][4]);
    if (!strcmp( string[i][3], "") || !strcmp( string[i][3],"mono") || !strcmp(string[i][3], "Mono") || 
        !strcmp(string[i][3],"monodisperse") || !strcmp(string[i][3], "Monodisperse"))
      ssRadar->allDistributionParameters[i-2]->kindOfDistribution = DIST_NONE;
    else if (!strcmp(string[i][3], "gamma") || !strcmp(string[i][3], "Gamma")) 
      ssRadar->allDistributionParameters[i-2]->kindOfDistribution = DIST_GAMMA;
    else if (!strcmp(string[i][3], "lognormal") || !strcmp(string[i][3], "Lognormal") || !strcmp(string[i][3], "log") || 
             !strcmp(string[i][3], "Log") || !strcmp(string[i][3], "logn") || !strcmp(string[i][3], "Logn"))
      ssRadar->allDistributionParameters[i-2]->kindOfDistribution = DIST_LOGNORMAL;
    else {
      ssRadar->allDistributionParameters[i-2]->kindOfDistribution = DIST_FILE;
      ssRadar->allDistributionParameters[i-2]->sizeDistributionFilename = calloc(256,sizeof(char));
      strcpy(ssRadar->allDistributionParameters[i-2]->sizeDistributionFilename, string[i][3]);
    }
    ssRadar->allDistributionParameters[i-2]->distributionMaxRadiusFactor = atof(string[i][5]);
    ssRadar->allDistributionParameters[i-2]->temperature = atof(string[i][6]);
    ssRadar->allDistributionParameters[i-2]->LWC = 0.001*atof(string[i][8]); /* kg per cubic meter */
    if (!strcmp(string[i][7], "water") || !strcmp(string[i][7], "Water")) {
      ssRadar->allDistributionParameters[i-2]->waterPhase = 1;
    }
    else if (!strcmp(string[i][7], "ice") || !strcmp(string[i][7], "Ice")) {
      ssRadar->allDistributionParameters[i-2]->waterPhase = -1;
    }
    else {
      fprintf(stderr,"Error in SSRadar, unknown phase!\n");
      return -1;
    }
    if ((ssRadar->allDistributionParameters[i-2]->waterPhase == 1 && 
        (ssRadar->allDistributionParameters[i-2]->temperature > 55. ||
         ssRadar->allDistributionParameters[i-2]->temperature < -34.)) ||
        (ssRadar->allDistributionParameters[i-2]->waterPhase == -1 && 
        (ssRadar->allDistributionParameters[i-2]->temperature > 0. ||
         ssRadar->allDistributionParameters[i-2]->temperature < -100.))) {
    fprintf(stderr,"Error in SSRadar, wrong combination of temperature and phase %d %f!\n",
            ssRadar->allDistributionParameters[i-2]->waterPhase,
            ssRadar->allDistributionParameters[i-2]->temperature);
    return -1;
    }
    ssRadar->allDistributionParameters[i-2]->distributionStep = 0.0005;
  }
  /*  if ((ssRadar->radarZPosition > ssRadar->firstRangeGate && ssRadar->cosZenithAngleUmu > 0.)
		  || (ssRadar->radarZPosition < ssRadar->firstRangeGate && ssRadar->cosZenithAngleUmu < 0.)) {
	  fprintf(stderr,"Error, First Range Gate and Radar Z-Position not compatible!\n");
	  return -1;
  } */
  for (i=0;i<ssRadar->numberOfReflectiveLayers-1;i++) {
    if (ssRadar->allDistributionParameters[i+1]->height < ssRadar->allDistributionParameters[i]->height) {
      fprintf(stderr,"Error, layers must be in ascending order.\n");
      return -1;
    }
  }
  if (ASCII_free_string ( string, rows, maxColumns ) ) {
    fprintf(stderr,"Error in ASCII_free_string.\n");
    return -1;
  }
  return 0;
}

/***********************************************************************************/
/* Function: writeSSRadarOutput                                           @62_30i@ */
/* Description: outputs ssradar into txt and mmclx, modified from mystic.c         */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: ...,Christian Pause                                                     */
/*                                                                        @i62_30@ */
/***********************************************************************************/
int writeSSRadarOutput(double* rangeGatesVertical,
                       double* rangeGates,
                       double* rangeGateReflectivityFactor,
                       double  radarZPosition,
                       double  rangeGateWidth,
                       double  wavelengthInMicrometer,
                       double  cosZenithAngleUmu,
                       char*   outputFilename,
                       int     numberOfRangeGates,
                       int     output) {


/* output: 1 = only mmclx, -1 = only text, 0 = both */

#if HAVE_LIBNETCDF
  int status=0;
  int ncid=0, iid=0;
  int id_Nt=0;
  int id_var[100];
  char* mmclxFilename;
  char name[FILENAME_MAX], unit[10];
  int dimids[10];
  size_t start[10];
  size_t count[10];
  double range[numberOfRangeGates];

  double *Zerange = NULL;
  Zerange = calloc(2,sizeof(double));
  Zerange[0] = -60.;
  Zerange[1] = 80.;

  int *dB = NULL;
  dB = calloc(1,sizeof(int));
  *dB = 1;

  double theta = acos(cosZenithAngleUmu);

  double za[1];
  double wvl[1];

  mmclxFilename = calloc(256,sizeof(char));
#endif
  int i=0;
  FILE *ssRadarOutputFile = NULL;
  char* ssRadarOutputFilename;
  ssRadarOutputFilename = calloc(256,sizeof(char));

  if (*outputFilename != *ssRadarOutputFilename) {
    if (output > -1) {
#if HAVE_LIBNETCDF

      for (i=0; i<10; i++) {
        start[i]=0;
        count[i]=1;
      }
  
      mmclxFilename = strtok(outputFilename,".");
      strcat(mmclxFilename,".mmclx");
  
      status = nc_create(mmclxFilename, NC_CLOBBER, &ncid);
      if (status)
        return err_out ("Error %d creating mmclx file \n", status);
    
      /* define global attributes and dimensions */
    
      za[0] = 180*theta/PI;
      wvl[0] = wavelengthInMicrometer/1000.;
      nc_put_att_double(ncid, NC_GLOBAL, "zenith_angle_in_degree", NC_DOUBLE, 1, za);
      nc_put_att_double(ncid, NC_GLOBAL, "wavelength_in_mm", NC_DOUBLE, 1, wvl);
      
      status = nc_def_dim(ncid, "range", numberOfRangeGates, &id_Nt);
      if (status)
        return err_out ("Error %d creating dimension range \n", status);
    
      /* Define arrays */
    
      iid=0;
      /* range */
      dimids[0] = id_Nt;
    
      status = nc_def_var(ncid, "range", NC_DOUBLE, 1, dimids, &id_var[iid]);
      if (status)
        return  err_out ("Error %d defining variable range \n", status);
    
      strcpy(name,"range gates");
      strcpy(unit,"m");
      nc_put_att_text(ncid, id_var[iid], "long_name", strlen(name), name);
      nc_put_att_text(ncid, id_var[iid], "units", strlen(unit), unit);
    
      /* Ze */
      iid++;
      dimids[0] = id_Nt;
    
      status = nc_def_var(ncid, "Ze", NC_DOUBLE, 1, dimids, &id_var[iid]);
      if (status)
        return  err_out ("Error %d defining variable Ze \n", status);
    
      strcpy(name,"equivalent radar reflectivity factor Ze of hydrometeors");
      strcpy(unit,"mm^6/m^3");
      nc_put_att_text(ncid, id_var[iid], "long_name", strlen(name), name);
      nc_put_att_text(ncid, id_var[iid], "units", strlen(unit), unit);
      nc_put_att_double(ncid, id_var[iid], "yrange", NC_DOUBLE, 2, Zerange);
      nc_put_att_int(ncid, id_var[iid], "db", NC_INT, 1, dB);
    
      /* done */
      status = nc_close(ncid);
      if (status)
        return  err_out ("Error %d closing netcdf file \n", status);
    
      /* Write arrays */
    
      for (i=0; i<10; i++) {
        start[i]=0;
        count[i]=1;
      }
    
      status = nc_open(mmclxFilename, NC_WRITE, &ncid);
      if (status)
        return err_out ("Error %d opening netcdf file \n", status);
    
      iid=0;
      /* range */
      for (i=0;i<numberOfRangeGates-1;i++)
        range[i] = (rangeGates[i] + 0.5*(rangeGates[i+1] - rangeGates[i]));
      range[numberOfRangeGates-1] = (rangeGates[numberOfRangeGates-1] + 0.5*(rangeGates[numberOfRangeGates-1] - rangeGates[numberOfRangeGates-2]));
      for (i=0;i<numberOfRangeGates;i++) {
        start[0] = i;
        status = nc_put_vara_double(ncid, id_var[iid],
                                    start, count, &range[i]);
      if (status)
        return  err_out ("Error %d writing variable range \n", status);
      }
    
      iid++;
      for (i=0;i<numberOfRangeGates;i++) {
        start[0] = i;
        status = nc_put_vara_double(ncid, id_var[iid], start, count,
                                    &rangeGateReflectivityFactor[i]);
        if (status)
          return  err_out ("Error %d writing variable Ze \n", status);
      }
  
      free(Zerange);
      free(dB);
      status = nc_close(ncid);
      if (status)
        return  err_out ("Error %d closing netcdf file \n", status);
    
#else
      fprintf (stderr, " ***********************************************************************\n");
      fprintf (stderr, " * You have built uvspec without libnetcdf and hence cannot            *\n");
      fprintf (stderr, " * use write_mmclx_locest_file. Please get netcdf and rebuild.         *\n");
      fprintf (stderr, " ***********************************************************************\n");
      return -1;
#endif
    }
    if (output < 1) {
  
      ssRadarOutputFilename = strtok(outputFilename,".");
      strcat(ssRadarOutputFilename,".out");
      if ((ssRadarOutputFile = fopen (ssRadarOutputFilename, "w")) == NULL) {
        fprintf (stderr, "Error opening %s for writing\n",ssRadarOutputFilename);
        return -1;
      }
  
      fprintf(ssRadarOutputFile,"Wavelength: %d um\n", (int)wavelengthInMicrometer );
      fprintf(ssRadarOutputFile,"Height of radar: %e m\n", radarZPosition );
      fprintf(ssRadarOutputFile,"Zenith Angle: %e\n", 180.*acos(cosZenithAngleUmu)/PI );
      fprintf(ssRadarOutputFile,"Number of Range Gates: %d\n", numberOfRangeGates );
      fprintf(ssRadarOutputFile,"Range Gate Width: %e m\n", rangeGateWidth );
      fprintf(ssRadarOutputFile,"First Range Gate Radial Distance: %e m\n", rangeGates[0] );
      fprintf(ssRadarOutputFile,"First Range Gate Z-Position: %e m\n", rangeGatesVertical[0] );
      fprintf(ssRadarOutputFile,"RG No \t RG Dist [m] \t RG Z-Pos [m] \t Z [dBZ]\n" );
      for (i=0;i<numberOfRangeGates;i++)
        fprintf(ssRadarOutputFile,"%d \t %e \t %e \t %f\n",i+1,rangeGates[i],rangeGatesVertical[i],rangeGateReflectivityFactor[i] );
  
      fclose(ssRadarOutputFile);
    }
  }

  fprintf(stdout,"Wavelength: %d um\n", (int)wavelengthInMicrometer );
  fprintf(stdout,"Height of radar: %e m\n", radarZPosition );
  fprintf(stdout,"Zenith Angle: %e\n", 180.*acos(cosZenithAngleUmu)/PI );
  fprintf(stdout,"Number of Range Gates: %d\n", numberOfRangeGates );
  fprintf(stdout,"Range Gate Width: %e m\n", rangeGateWidth );
  fprintf(stdout,"First Range Gate Radial Distance: %e m\n", rangeGates[0] );
  fprintf(stdout,"First Range Gate Z-Position: %e m\n", rangeGatesVertical[0] );
  fprintf(stdout,"RG No \t RG Dist [m] \t RG Z-Pos [m] \t Z [dBZ]\n" );
  for (i=0;i<numberOfRangeGates;i++)
    fprintf(stdout,"%d \t %e \t %e \t %f\n",i+1,rangeGates[i],rangeGatesVertical[i],rangeGateReflectivityFactor[i] );

  return 0;

}


/***********************************************************************************/
/* Function: findEmptyLayers                                             @62_30i@ */
/* Description: Finds empty layers                                                 */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Christian Pause                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/
int findEmptyLayers(double height0, double height1, double thickness) {
  if (height0 + thickness == height1)
    return 0; 
  else if (height0 + thickness < height1)
    return 1;
  else return -1;
}



/***********************************************************************************/
/* Function: linearInterpolation                                         @62_30i@ */
/* Description: linear interpolation                                               */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Christian Pause                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/
double linearInterpolation(double a, double b, double c, double d, double e) {
  return a + ((b-a)/(d-c))*(e-c);
}


/***********************************************************************************/
/* Function: exitSSRadar                                                  @62_30i@ */
/* Description: free allocated stuff on exit                                       */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Christian Pause                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/
void exitSSRadar(distributionParameterStruct*** allDistributionParameters,
                double** rangeGates,
                double** rangeGatesVertical,
                double** rangeGateReflectivityFactor,
                double** reflectiveLayersHeights, 
                double** reflectiveLayersThicknesses,
                double** reflectiveLayersReflectivities,
                double** allLayersHeights, 
                double** allLayersThicknesses,
                double** allLayersReflectivities,
                int** indexArray,
                int numberOfReflectiveLayers) {

  int i;
  if (numberOfReflectiveLayers != 0) {
    for(i=0;i<numberOfReflectiveLayers;i++) {
      free(((*allDistributionParameters)[i])->sizeDistributionFilename);
      free(((*allDistributionParameters)[i])->dropletDistribution);
      free(((*allDistributionParameters)[i])->radiusGrid);
      free((*allDistributionParameters)[i]);
    }
    free(*allDistributionParameters);
  }
  free(*rangeGates);
  free(*rangeGatesVertical);
  free(*rangeGateReflectivityFactor);
  free(*reflectiveLayersHeights);
  free(*reflectiveLayersThicknesses);
  free(*reflectiveLayersReflectivities); 
  free(*indexArray); 
  free(*allLayersHeights); 
  free(*allLayersThicknesses); 
  free(*allLayersReflectivities); 
  return;
}


/***********************************************************************************/
/* Function: getWaterDensity                                              @62_30i@ */
/* Description: return density of water for a given temperature                    */
/* water: 0 - 55 degree C, supercooled water: -34 - 0 degree C,                    */
/*   ice: -100 - 0 degree C                                                        */
/* Parameters:                                                                     */
/* Return value:                                                                   */
/* Example:                                                                        */
/* Files:                                                                          */
/* Known bugs:                                                                     */
/* Author: Christian Pause                                                         */
/*                                                                        @i62_30@ */
/***********************************************************************************/
double getWaterDensity(double temperature, int phase) {
  double iceDensityTable[6] = {929.,927.,935.,923.,920.,917.,};
  double waterDensityTable[12] = {999.84, 999.96, 999.7 , 999.1 , 998.2 , 997.04,
                               995.64, 994.03, 992.21, 990.22, 988.04, 985.7 };
  double supercooledWaterDensityTable[7] = {977.5,983.9,989.5,993.5,996.3,998.2,999.3};
  double iceTemperatureTable[6] = {-100.,-80.,-60.,-40.,-20.,0.};
  double waterTemperatureTable[12] = {0.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.};
  double supercooledWaterTemperatureTable[7] = {-34.,-30.,-25.,-20.,-15.,-10.,-5.};
  int index,i;
  /* supercooled water density data taken from Hare and Sorensen (1987) */
  if (temperature >= 0. && phase == 1) {
    for (i=0;i<12;i++) {
      if(waterTemperatureTable[i] == temperature) {
        return waterDensityTable[i];
      }
    }
    index = locate(waterTemperatureTable,12,temperature);
    return linearInterpolation(waterDensityTable[index],waterDensityTable[index+1],
                                waterTemperatureTable[index],waterTemperatureTable[index+1],temperature);
  }
  else if (temperature <= 0. && phase == 1) {
    for (i=0;i<12;i++) {
      if(supercooledWaterTemperatureTable[i] == temperature) {
        return supercooledWaterDensityTable[i];
      }
    }
    index = locate(supercooledWaterTemperatureTable,7,temperature);
    return linearInterpolation(supercooledWaterDensityTable[index],supercooledWaterDensityTable[index+1],
                                supercooledWaterTemperatureTable[index],supercooledWaterTemperatureTable[index+1],temperature);
  }
  else if (temperature <= 0. && phase == -1) {
    for (i=0;i<12;i++) {
      if(iceTemperatureTable[i] == temperature) {
        return iceDensityTable[i];
      }
    }
    index = locate(iceTemperatureTable,6,temperature);
    return linearInterpolation(iceDensityTable[index],iceDensityTable[index+1],
                                iceTemperatureTable[index],iceTemperatureTable[index+1],temperature);
  }
  else
    return -9999.;
}

/***********************************************************************************/
/* Prints an array of doubles with name */
/* Only for diagnostic purposes */
void printArrayDouble(double *A, int n, char* C) {
  int i;
  for(i=0;i<n;i++) {
    fprintf(stdout,"%s",C);
    fprintf(stdout," %e\n",A[i]);
  }
}

/***********************************************************************************/
/* Prints an array of ints with name */
/* Only for diagnostic purposes */
void printArrayInt(int *A, int n, char* C) {
  int i;
  for(i=0;i<n;i++) {
    fprintf(stdout,"%s",C);
    fprintf(stdout," %d\n",A[i]);
  }
}

/***********************************************************************************/
/* locates a value x in an array xx of length n */
/* returns n if value is over the highest bond */
int locate2 (double *xx, int n, double x, double umu) {
  int i;
  if (umu > 0.) {
    if (x > xx[n-1])
      return n;
    else {
      for(i=0;i<n-1;i++) {
        if (x >= xx[i] && x < xx[i+1])
          return i;
      if (x >= xx[n-1])
        return n-1;
      }
    }
  }
  else {
    if (x > xx[n-1])
      return n;
    else {
      for(i=0;i<n-1;i++) {
        if (x > xx[i] && x <= xx[i+1]) {
          return i;
      }
      if (x > xx[n-1])
        return n-1;
      }
    }
  }
  return 0;
}



/***********************************************************************************
 * calc_distribution_radius_grid
 * 
 * Calculate a radius grid for the size distribution over which the
 * optical properties are intgrated. This grid must be sufficiently
 * fine, otherwise "ripples" in the phase function are not averaged
 * out. Here the sampling is done in the same way as in cloudprp by
 * F. Evans.
 *
 * @author ???
 * @date   ???
 * taken and modified from mie.c
 */
int calc_distribution_radius_grid(/* Input */
                                  int distribution, 
                                  char* sd_filename,
                                  double r_eff_min,  /* radius given in um */
                                  double r_eff_max,  /* radius given in um */
                                  double lambda,      /*wavelength in micrometer*/
                                  double dx_max,
                                  double n_r_max,
                                  int verbose,
                                  /* Output */
                                  int* n_dens,
                                  double** radius,  /* radius given in um */
                                  double** number_dens 
                                  )
{
  double radmin=0.0, radmax=0.0, rad=0.0, delta_rad=0.0;
  double x=0.0, delta_x=0.0, dx=0.0;
  int is=0;
  int status=0;
 
  switch(distribution) {
  
  case DIST_NONE:
    /* single radius calculation, nothing to do */
    *n_dens = 1;
    if (((*radius) = calloc(*n_dens, sizeof(double)))==NULL) {
      fprintf (stderr, "Error: Allocation of output->dist.radius (make_distribution)\n");
      return -2;
    }
    if (((*number_dens) = calloc(*n_dens, sizeof(double)))==NULL) {
      fprintf (stderr, "Error: Allocation of output->dist.numb_dens (make_distribution)\n");
      return -3;
    } 
    break;
    
  case DIST_FILE:
    /* read size distribution from file */
    status = read_2c_file (sd_filename, &(*radius), &(*number_dens), &(*n_dens));
 
    if (status != 0) {
      fprintf (stderr, "error %d reading %s\n", status, sd_filename);
      return status;
    }

    if (verbose)
      fprintf (stderr, " ... read %d data points from %s\n", 
               *n_dens, sd_filename);
    
    break;
  case DIST_GAMMA:
  case DIST_LOGNORMAL:
    radmin = 0.02 * r_eff_min; /* cover the whole droplet spectrum   */
    radmax = n_r_max * r_eff_max; /* in order to save mie-calculations  */
 
    dx=dx_max;
    
    /* count number of radii */
    rad = radmin;
    is = 0;
    while (rad <= radmax) {
      x = 2*PI*rad / lambda; 
      delta_x = double_max(dx,dx*sqrt(x));
      delta_rad = delta_x * lambda / (2*PI);
      rad = rad + delta_rad;
      is = is + 1;
    }
    *n_dens = is;
   
    if (verbose)
      fprintf(stderr, "size distributions: r_min %g um, r_max %g um, number of grid points %d \n", radmin, radmax, *n_dens); 
    
    if (((*radius) = calloc(*n_dens, sizeof(double)))==NULL) {
      fprintf (stderr, "Error: Allocation of output->dist.radius (make_distribution)\n");
      return -2;
    }
    if (((*number_dens) = calloc(*n_dens, sizeof(double)))==NULL) {
      fprintf (stderr, "Error: Allocation of output->dist.numb_dens (make_distribution)\n");
      return -3;
    } 

    /* calculate radius grid */
    rad = radmin;
    is = 0;
    while (rad <= radmax) {
      (*radius)[is]      = rad;
      x = 2*PI*rad / lambda; 
      delta_x = double_max(dx,dx*sqrt(x));
      delta_rad = delta_x * lambda / (2*PI);
      rad = rad + delta_rad;
      is = is + 1;
    }
    break; 
    
  case DIST_AER: 
    /* Do nothing here. Radius grid and number densities are calculated in
       aerosol_size_distribution. */
    break;
  default:
    fprintf (stderr, "Error in calc_distribution_radius_grid: unrecognized distribution %d\n", distribution);
    return -1;
    break;
  } 
  
  /* check if data points are sorted in ascending order */
  switch(distribution) {
  case DIST_NONE:
    break;
  case DIST_FILE:
  case DIST_GAMMA:
  case DIST_LOGNORMAL:
    /* for gamma and lognormal this is not nessesary, but it costs almost no time */
    
    for (is=0; is<*n_dens-1; is++)
      if ( (*radius)[is+1]< (*radius)[is]) {
        if (distribution == DIST_FILE)
          fprintf (stderr, "Error: %s not sorted in ascending order\n", 
                   sd_filename);
        else
          fprintf (stderr, "Error: radius in distribution not sorted in ascending order, internal program bug!!!\n");
        return -5;
      }
    
    if ((*radius)[0]==0) {
      (*radius)[0] = 0.1 * (*radius)[1];
      fprintf (stderr, " ... cannot handle r=0; setting first data point to r=%g\n",
               (*radius)[0]);
    }
    break;
  default:
    fprintf (stderr, "Error in calc_distribution_radius_grid: unrecognized distribution %d\n", distribution);
    break;
  }
  return status;
}

/***********************************************************************************
 * calc_distribution_number_dens
 *
 * Calculate particle number densities for each radius in the size
 * distribution.
 *
 * @author ???
 * @date   ???
 * taken and modified from mie.c
 *
 * r_eff is in um!!!
 * LWC is in kg/m^3!!!
 * rho_medium is in kg/m^3!!!
 */
int calc_distribution_number_dens (/* Input */
                                   int distribution,
                                   double alpha,
                                   double r_eff,
                                   double rho_medium,
                                   double LWC,
                                   int n_dens,
                                   double* radius,
                                   int verbose, 
                                   /* Output */
                                   double* number_dens
                                   )
{
  double A=0.0, B=0.0;
  int is=0;
  int status=0;
  double integral = 0.;

  if (rho_medium <= 0.0) {
    fprintf (stderr, "Error, rho_medium = %f is 0 or negative!\n", rho_medium);
    return -1;
  }
  switch(distribution) {

  case DIST_NONE:
    number_dens[0] = (3.*LWC)/(4.*PI*rho_medium*pow(0.000001*r_eff,3.));
    radius[0] = r_eff;
    break;

  case DIST_FILE:
    for (is=0;is<n_dens-1;is++)
      radius[is] += 0.5*(radius[is+1]-radius[is]);
    radius[n_dens-1] += 0.5*(radius[n_dens-1]-radius[n_dens-2]);

    /* normalize distribution and get LWC into it */
    for (is=0;is<n_dens;is++) {
	    integral += number_dens[is]*pow(radius[is],3.);
    }
    for (is=0;is<n_dens;is++) {
      number_dens[is] *= (3.*LWC*pow(radius[is],3.))/(4*PI*rho_medium*pow(0.000001*radius[is],3.)*integral);
    }
    break;
    
  case DIST_GAMMA:
  
    B = (alpha+3.0)/r_eff;
    A = (0.75/PI)*pow(B,(alpha+4))/exp(gamma_ln(alpha+4.));
    /* calculate number density */
    for (is=0;is<n_dens;is++) {
      number_dens[is] = A * pow(radius[is], alpha) * exp(-B*radius[is]);
      integral += number_dens[is]*pow(radius[is],3.);
    }
    for (is=0;is<n_dens;is++) {
      number_dens[is] *= (3.*LWC*pow(radius[is],3.))/(4*PI*rho_medium*pow(0.000001*radius[is],3.)*integral);
    }

    break;
    
  case DIST_LOGNORMAL:
    
    if (verbose)
      fprintf(stderr, "log-normal distribution, alpha=%f, r_mod=%f \n", alpha, r_eff);
    
    /* specified r_eff is taken as r_mod for log normal distribution */
    for (is=0;is<n_dens;is++) {
      A=(log((radius)[is])-log(r_eff))/log(alpha); 
      number_dens[is] = 1.0 /
      (sqrt(2.0*PI)*log(alpha) * radius[is]) * exp(-0.5*A*A); 
      integral += number_dens[is]*pow(radius[is],3.);
    }
    for (is=0;is<n_dens;is++) {
      number_dens[is] *= (3.*LWC*pow(radius[is],3.))/(4*PI*rho_medium*pow(0.000001*radius[is],3.)*integral);
    }

    break;
  default:
    fprintf (stderr, "Error in make_distribution : unrecognized distribution %d\n", distribution);
    status=-1;
    break;
  } /* switch(distribution) */
  
  return status;
}

/* taken from mie.c */
double gamma_ln(float xx)
{
  int j=0;
  double *cof;
  double stp = 2.50662827465;
  double half = 0.5, one=1.0, fpf = 5.5, x, tmp, ser;
  double gammaln = 0.0;

  cof  = (double *) calloc (6, sizeof(double));
  cof[0] =  76.18009173;
  cof[1] = -86.50532033;
  cof[2] =  24.01409822;
  cof[3] =  -1.231739516;
  cof[4] =   0.120858003;
  cof[5] =  -0.53638200;

  x = xx - one;
  tmp = x + fpf;
  tmp = (x+half)*log(tmp)-tmp;
  ser=one;
  for (j=0; j<=5; j++) {
    x = x + one;
    ser = ser + cof[j]/x; 
  }
  gammaln=tmp+log(stp*ser);
 
  free(cof); 
  
  return gammaln;
}

double double_max(double x, double y)
{
  if (x >= y)
    return x;
  else
    return y;
}

