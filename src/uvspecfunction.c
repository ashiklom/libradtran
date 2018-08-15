/************************************************************************
 * $Id: uvspecfunction.c 3330 2017-12-20 09:27:56Z svn-kylling $
 ************************************************************************/

#include "ascii.h"
#include <stdio.h>
#include <math.h>

#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <errno.h>
#include <string.h>

int calluvspec(char *outfilename);


int main(int argc, char **argv)
{
  extern FILE *yyin;
  FILE *fin=NULL;
  int i=0, status=0;
  char infilename[FILENAME_MAX]="", outfilename[FILENAME_MAX]="";

  double *heat=NULL, *z=NULL;
  int filedes = -1;

  int PID=1;
  int rc=0;
  int n=0;

  int nprocesses=8;
  double *SZA = calloc (nprocesses, sizeof(double));

  for (i=0; i<nprocesses; i++)
    SZA[i] = 90.0 / (double) nprocesses * (double) i;



  /* start child process */
  for (i=0; i<nprocesses; i++) {

    if (PID != 0)   /* parent process  */
      PID = fork();

    if (PID == 0)  {  /* child process  */

      sprintf (infilename, "tmp%d.XXXXXX", (int) getpid());
      sprintf (outfilename, "tmp%d.XXXXXX", (int) getpid());
      filedes = mkstemp (infilename);
      if(filedes<1)	{
	fprintf(stderr, "\n Creation of temp file %s failed with error [%s]\n", infilename,strerror(errno));
	return 1;
      }
      filedes = mkstemp (outfilename);
      if(filedes<1)	{
	fprintf(stderr, "\n Creation of temp file %s failed with error [%s]\n", outfilename,strerror(errno));
	return 1;
      }
      
      fprintf (stderr, " ... writing to %s\n", infilename);

      fin = fopen (infilename,"w");

      fprintf (fin, "data_files_path ../data/\n");
      fprintf (fin, "atmosphere_file midlatitude_summer \n");
      fprintf (fin, "\n");
      fprintf (fin, "sza %f\n", SZA[i]);
      fprintf (fin, "\n");
      fprintf (fin, "\n");
      fprintf (fin, "\n");
      fprintf (fin, "rte_solver disort\n");
      fprintf (fin, "number_of_streams 32\n");
      fprintf (fin, "\n");
      fprintf (fin, "source solar\n");
      fprintf (fin, "\n");
      fprintf (fin, "mol_abs_param kato2\n");
      fprintf (fin, "output_process sum\n");
      fprintf (fin, "\n");
      fprintf (fin, "\n");
      fprintf (fin, "output_user zout heat\n");
      fprintf (fin, "\n");
      fprintf (fin, "zout  0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 16.0 17.0 18.0 19.0 20.0 21.0 22.0 23.0 24.0 25.0 27.50 30.0 32.50 35.0 37.50 40.0 42.50 45.0 47.50 50.0 \n");
      
      fclose(fin);

      yyin = fopen (infilename,"r");
      status = calluvspec (outfilename);

      if (status!=0)
	fprintf (stderr, "Error, calluvspec returned %d\n", status);
      
      read_2c_file (outfilename, &z, &heat, &n);
      for (i=0;i<n;i++)
	fprintf (stdout, "%f %f\n", z[i], heat[i]);
      
      fclose (yyin);
      
      free(z);
      free(heat);
      
      remove (infilename);
      remove (outfilename);

      exit(status);
    }
  }
  
  wait(&rc);
  fprintf (stderr, "done!\n");


  return status;
}
