/* xytohex.c */
/*
 * This set of routines provide conversion between 
 * regular coordinate system (x,y) and coordinate
 * system of square cells covering all 2D space. 
 * The spacing between cells in x direction is equal 
 * to 1 which is also the side to side cell size.
 * The cell in origin of (x,y) plane is cell number
 * one. The center of cell number two is at (1,0).
 * Then cells of the first hexagonal ring are numbered 
 * in clockwise direction around (x,y) origin. The center 
 * of cell number ten is at (2,0) and cells of the 
 * second square ring are numbered in clockwise 
 * direction again. This numbering process is continued
 * until all 2D space is covered. There are two routines
 * available for a user. Routine  
 *
 * xy_to_ns(double *x, double *y, int *n);
 *
 * maps (x,y) coordinates onto cell number n and (x,y)
 * pair of coordinates in translated reference frame 
 * which has origin at the center of the cell number n.
 * Routine 
 *
 * nh_to_xy(int n, double *x, double *y);
 *
 * performs inverse transformation. Using given cell
 * number, n, it returns (x,y) coordinates of the
 * cell center. Because all cells are numbered by 
 * integers >0, n must be >0. If routine is called 
 * with negative or zero n, n value will be re-assigned 
 * to 1 (origin).
 *
 * March 25, 2008
 *
 * sjf
 */ 

#include <assert.h>
#include "xytosquare.h"

#define MAX(x,y) ((x)>(y)?(x):(y))

void xy_to_ns(double *x, double *y, int *ns)
/*
 * converts (x,y) into square cell number,
 * and returns (x,y) in the reference frame in which 
 * cell center is placed in origin. It is assumed that
 * spacing of cells in x direction is equal to 1. Side
 * to side cell size is therefore equal to 1 too.
 *
 * sjf - March 25, 2008
 */
{
  const double rx[8] = {  1.0,  1.0,  0.0, -1.0, -1.0, -1.0,  0.0,  1.0 };
  const double ry[8] = {  0.0, -1.0, -1.0, -1.0,  0.0,  1.0,  1.0,  1.0 };
  const double ex[8] = {  0.0, -1.0, -1.0,  0.0,  0.0,  1.0,  1.0,  0.0 };
  const double ey[8] = { -1.0,  0.0,  0.0,  1.0,  1.0,  0.0,  0.0, -1.0 };

  const double fabsx = fabs(*x);
  const double fabsy = fabs(*y);
  const int iring = lround(MAX(fabsx,fabsy));
  if(iring==0)*ns = 1;
  else
    {
      const double r = double(iring);
      for(unsigned iside=0;iside<8;iside++)
	{
	  const double sx = *x - r*rx[iside];
	  const double sy = *y - r*ry[iside];
	  const double esx = sx*ex[iside] + sy*ey[iside];
	  const double esy = sx*ey[iside] - sy*ex[iside];
	  const int iesx = lround(esx);
	  const int iesy = lround(esy);
	  if((iesy==0)&&(iesx>=0)&&(iesx<=iring))
	    {
	      *x = sx - iesx*ex[iside];
	      *y = sy - iesx*ey[iside];
	      *ns = (2*iring-1)*(2*iring-1)+iside*iring+iesx+1;
	      break;
	    }
	  assert(iside != 7);
	}
    }
  return;
}

/* by Jim Ulery */
static unsigned julery_isqrt(unsigned long val) {
  unsigned long temp, g=0, b = 0x8000, bshft = 15;
  do {
    if (val >= (temp = (((g << 1) + b)<<bshft--))) {
      g += b;
      val -= temp;
    }
  } while (b >>= 1);
  return g;
}

#include<stdio.h>

void ns_to_xy(int ns, double *x, double *y)
/*
 * Converts square cell number into (x,y) coordinates.
 * Cell number must always be integer > 0. If it is <= 0,
 * it is assigned to 1.	It is assumed that spacing of cells 
 * in x direction is equal to 1. Side to side 
 * cell size is therefore equal to 1 too.
 *
 * sjf - March 25, 2008
 */
{
  const double rx[8] = {  1.0,  1.0,  0.0, -1.0, -1.0, -1.0,  0.0,  1.0 };
  const double ry[8] = {  0.0, -1.0, -1.0, -1.0,  0.0,  1.0,  1.0,  1.0 };
  const double ex[8] = {  0.0, -1.0, -1.0,  0.0,  0.0,  1.0,  1.0,  0.0 };
  const double ey[8] = { -1.0,  0.0,  0.0,  1.0,  1.0,  0.0,  0.0, -1.0 };

  if(ns<=1)
    {
      *x = 0; 
      *y = 0; 
    }
  else
    {
      const unsigned iring = ((julery_isqrt(ns-1)-1)>>1)+1;
      const unsigned iringcell = (ns-(2*iring-1)*(2*iring-1))-1;
      const unsigned iside = iringcell/iring;
      const unsigned icell = iringcell%iring;
      *x = double(iring)*rx[iside] + double(icell)*ex[iside];
      *y = double(iring)*ry[iside] + double(icell)*ey[iside];
    }
  return;
}

#ifdef TESTMAIN1
#include<stdio.h>
int main(int argc, char** argv)
{
  double x = -4;
  
  for(double y = -4.0; y<=4.0; y+=0.1)
    {
      for(double x = -4.0; x<=4.0; x+=0.1)
	{
	  double rx = x;
	  double ry = y;
	  int ns;
	  xy_to_ns(&rx,&ry,&ns);
	  printf("%d ",ns);
	}
      printf("\n");
    }
 
}
#endif

#ifdef TESTMAIN2
#include<stdio.h>
int main(int argc, char** argv)
{
  for(unsigned i=0;i<100;i++)
    {
      double x;
      double y;
      ns_to_xy(i,&x,&y);
      printf("%d %f %f\n",i,x,y);
    }
}
#endif
