/************************************************************************
 * $Id: bandec.c 2623 2011-12-23 10:52:38Z robert.buras $
 ************************************************************************/

/*--------------------------------------------------------------------
 * This subroutine was implemented for the two-stream model of
 * Robert Buras (rodents)
 * ulrike 16.06.2010
 *--------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#define SWAP(a,b) {dum=(a);(a)=(b);(b)=dum;}
#define TINY 1.0e-20

#include "bandec.h"

/***********************************************************************/
/* Given an n x n band diagonal matrix A with m1 subdiagonal rows and  */
/* m2 superdiagonal rows compactly stored in the array		       */
/* a[0..n-1][0..4] as described in the comment for routine banmul,     */
/* this routine constructs an LU decomposition of a rowwise permuation */
/* of A. The upper triangular matrix replaces a, while the lower       */
/* triangular matrix is returned in al[0..n-1][0..1]. indx[0..n-1] is  */
/* an output vector which records the row permutation effected by the  */
/* partial pivoting; d is output as +-1 depending on whether the       */
/* number of row interchanges was even or odd, respectively. This      */
/* routine is used in combination with banbks to solve band-diagonal   */
/* sets of equations.                                                  */
/***********************************************************************/

void bandec(double **a, long n, double **al,
	    long indx[])
{
  long i,j,k,l;
  double dum;

  l=2;
  for (i=0;i<2;i++) { /* Rearrange the storage a bit. */
    for (j=2-i;j<5;j++) a[i][j-l]=a[i][j];
    l--;
    for (j=4-l;j<5;j++) a[i][j]=0.0;
  }

  l=2;
  for (k=0;k<n;k++) { /* For each row ... */
    dum=a[k][0];
    i=k;
    if (l < n) l++;
    for (j=k+1;j<l;j++) { /* Find the pivot element. */
      if (fabs(a[j][0]) > fabs(dum)) {
	dum=a[j][0];
	i=j;
      }
    }
    indx[k]=i;
    if (dum == 0.0)
      a[k][0]=TINY;
    /* Matrix is algorithmically singular, but proceed anyway
       with TINY pivot (desirable in some applications). */
    if (i != k) { /* Interchange rows. */
      for (j=0;j<5;j++)
	SWAP(a[k][j],a[i][j]);
    }
    for (i=k+1;i<l;i++) { /* Do the elimination. */
      dum=a[i][0]/a[k][0];
      al[k][i-k-1]=dum;
      for (j=1;j<5;j++)
	a[i][j-1] = a[i][j] - dum * a[k][j];
      a[i][4]=0.0;
    }
  }
}

/***********************************************************************/
/* Given the arrays a, al, and indx as returned from bandec, and given */
/* a right-hand side vector b[0..n-1], solves the band diagonal linear */
/* equations Ax=b. The solution vector x overwrites b[0..n-1]. The     */
/* other input arrays are not modified, and can be left in place for   */
/* successive calls with different right-hand sides.                   */
/***********************************************************************/

void banbks(double **a, long n, double **al,
	    long indx[], double b[])
{
  long i=0,k=0,l=0;
  double dum;

  l=2;
  for (k=0;k<n;k++) { /* Forward substitution, unscrambling the
			  permuted rows as we go. */
    if (indx[k] != k)
      SWAP(b[k],b[i]);
    if (l < n) l++;
    for (i=k+1;i<l;i++)
      b[i] -= al[k][i-k-1] * b[k];
  }

  l=1;
  for (i=n-1;i>=0;i--) { /* Backsubstitution. */
    for (k=1;k<l;k++)
      b[i] -= a[i][k]*b[k+i];
    b[i]/=a[i][0];
    if (l < 5) l++;
  }
}
