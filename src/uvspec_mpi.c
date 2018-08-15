#ifdef HAVE_LIBTENSTREAM
#include <mpi.h>
#include "petscsys.h"
#include "tenstream.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_LIBTENSTREAM
static char help[] = "This is the C wrapper interface to the Tenstream solver environment.\n\n";
#endif

int calluvspec(char *outfilename);

int main(int argc, char **argv)
{
#ifdef HAVE_LIBTENSTREAM
  extern FILE *yyin;
  char outfilename[FILENAME_MAX]="";
  int numprocs,myid;
#endif
  char infilename[FILENAME_MAX]="";
  int ierr=0;

  int have_inputfile=0;

  if(argc>=1) have_inputfile=1;
  //TODO check if this is really a input file...

  if(!have_inputfile) {
    printf("We need a Libradtran Input file to work with, please specify one with -f <filename>\n");
    abort();
  }

  sprintf(infilename, argv[1]);
  printf("Using File %s as libradtran input\n",infilename);

#ifdef HAVE_LIBTENSTREAM
  MPI_Init(&argc,&argv);
  PetscInitialize(&argc,&argv,(char*)0,help);
  PetscInitializeFortran();

  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if(myid==0) {
//    sprintf(infilename, "run_box70.inp");
    sprintf(outfilename,"std.out");

    yyin = fopen (infilename,"r");
    ierr = calluvspec (outfilename);

    if (ierr!=0)
      fprintf (stderr, "Error, calluvspec returned %d\n", ierr);

    fclose (yyin);

    int func_index = MPI_FUNC_TENSTREAM_FINALIZE;
    MPI_Bcast(&func_index, 1, MPI_INT, 0, MPI_COMM_WORLD);

  }
  else
  {
    ierr = tenstream_slave_loop();
  }

  PetscFinalize();
  MPI_Finalize();
#else
  fprintf(stderr,"It seems your installation does currently not support the Tenstream solver... check your installation...\n");
  ierr = -1;
#endif
  return ierr;
}
