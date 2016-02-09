//-*-mode:c++; mode:font-lock;-*-

/*! \file VSReconstruction_fast.cpp

  Implementation of VVV's "in the sky chi^2" algorithm to calculate
  the RMS distance between the shower core (a line defined by four
  parameters) and the back-propagated photons detected by the
  instrument. Note that the value of chi^2 produced here is not a
  "chi^2" in the statistical sense, in that it is not normalized by
  the variances in the measured quantities (i.e. it is defined as
  "chi^2=Sum x_i^2" rather than "chi^2=Sigma x_i^2/sigma_x_i^2") so
  you should not test it statistically with the chi^2 distribution or
  take chi^2/DOF. The value calculated here does have the physically
  useful interpretation as the squared width of the emission region in
  space.... i.e. taking W=sqrt(chi^2/N_rays) gives a real width
  which can be selected against!

  The algorithm is simple when phrased correctly, but it cannot be
  explained in text-only comments! See the accompanying documentation.

  This file contains the code to calculate chi^2 for a single
  (hypothesized) shower direction. It is seperated from the "driver"
  functions which manage the overall minimization so that it can be
  compiled using a different compiler if desired. The Intel compiler
  will do a good job vectorizing this code if you let it.  There will
  be a *significant* increase in speed on a system with SSE2/SSE3
  extensions. GCC-4 also has vectorization, however it does not seem
  to work with this code. That does not surprise me, I had to make
  some modifications to get the "icc" to vectorize it, I suspect that
  I would only need to make minor changes to some of the pointer
  definitions to get "g++" to vectorize.

  The "#pragma ivdep" instruct the Intel compiler that the loops are
  vector operations.

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       04/16/2005

*/

// Intel/SSE3: icc -vec_report2 -xP -O3 -c VSReconstruction_fast.cpp
// Intel/SSE2: icc -vec_report2 -xN -O3 -c VSReconstruction_fast.cpp

// GNU/SSE3: g++4 -march=nocona -ftree-vectorize -ftree-vectorizer-verbose=5 -msse3 -I../Physics -m32 -O3 -c VSReconstruction_fast.cpp

//#define __VS_NO_IOSTREAM

#include<cmath>
#include<cfloat>

#include <Vec3D.hpp>

#include "VSReconstruction.hpp"

using namespace VERITAS;
using namespace Physics;

void VSReconstruction::
calculateChi2ForHypothesis(bool minimize_shower_core_location)
{
  for(unsigned iscope=0; iscope<nscopes;iscope++)
    {
      const unsigned nph = scope_nray[iscope];
      
      const double *const ws_p = ray_weight[iscope];
      const double *const es_x_p = es_x[iscope];
      const double *const es_y_p = es_y[iscope];
      const double *const es_z_p = es_z[iscope];
      
#pragma ivdep
      for(unsigned iph=0; iph<nph; iph++)
	es_dot_ec[iph] = 
	  ec_x*es_x_p[iph] + ec_y*es_y_p[iph] + ec_z*es_z_p[iph];
      
#pragma ivdep
      for(unsigned iph=0; iph<nph; iph++)
	one_minus_es_dot_ec_squared[iph] = 
	  1.0 - es_dot_ec[iph]*es_dot_ec[iph];
      
#if 0
      for(unsigned iph=0; iph<nph; iph++)
	if(one_minus_es_dot_ec_squared[iph] < DBL_EPSILON)
	  std::cout << iscope << '/' << iph << ' '
		    << ec_x << " * " << es_x_p[iph] << " + "
		    << ec_y << " * " << es_y_p[iph] << " + "
		    << ec_z << " * " << es_z_p[iph] << " = "
		    << es_dot_ec[iph] << ' ' 
		    << one_minus_es_dot_ec_squared[iph] << ' '
		    << 1.0/one_minus_es_dot_ec_squared[iph] 
		    <<std::endl;
#endif
      
#pragma ivdep
      for(unsigned iph=0; iph<nph; iph++)
	q_tau_x[iph] = 
	  (es_x_p[iph] - es_dot_ec[iph] * ec_x)/
	  one_minus_es_dot_ec_squared[iph];
      
#pragma ivdep
      for(unsigned iph=0; iph<nph; iph++)
	q_tau_y[iph] = 
	  (es_y_p[iph] - es_dot_ec[iph] * ec_y)/
	  one_minus_es_dot_ec_squared[iph];
#pragma ivdep
      for(unsigned iph=0; iph<nph; iph++)
	q_tau_z[iph] = 
	  (es_z_p[iph] - es_dot_ec[iph] * ec_z)/
	  one_minus_es_dot_ec_squared[iph];
      
#pragma ivdep
      for(unsigned iph=0; iph<nph; iph++)
	q_p_x[iph] = 
	  (ec_x - es_dot_ec[iph] * es_x_p[iph])/
	  one_minus_es_dot_ec_squared[iph];
#pragma ivdep
      for(unsigned iph=0; iph<nph; iph++)
	q_p_y[iph] = 
	  (ec_y - es_dot_ec[iph] * es_y_p[iph])/
	  one_minus_es_dot_ec_squared[iph];
#pragma ivdep
      for(unsigned iph=0; iph<nph; iph++)
	q_p_z[iph] = 
	  (ec_z - es_dot_ec[iph] * es_z_p[iph])/
	  one_minus_es_dot_ec_squared[iph];
      
#if 0
      for(unsigned iph=0; iph<nph; iph++)
	std::cout << "q: " << iscope << '/' << iph << ' ' 
		  << q_tau_x[iph] << ' ' 
		  << q_tau_y[iph] << ' ' 
		  << q_tau_z[iph] << ' ' 
		  << q_p_x[iph] << ' ' 
		  << q_p_y[iph] << ' ' 
		  << q_p_z[iph] << std::endl;
#endif
      
#pragma ivdep
      for(unsigned iph=0; iph<nph; iph++)
	if(one_minus_es_dot_ec_squared[iph] < DBL_EPSILON)
	  Q_xx[iph] = 
	    1.0 - ec_x*ec_x;
	else
	  Q_xx[iph] = 
	    1.0 + q_tau_x[iph]*q_tau_x[iph] + q_p_x[iph]*q_p_x[iph] 
	    - 2.0 * es_x_p[iph]*q_tau_x[iph]
	    - 2.0 * ec_x*q_p_x[iph]
	    + es_dot_ec[iph] * (2.0 * q_tau_x[iph] * q_p_x[iph]);
#pragma ivdep
      for(unsigned iph=0; iph<nph; iph++)
	if(one_minus_es_dot_ec_squared[iph] < DBL_EPSILON)
	  Q_yy[iph] = 
	    1.0 - ec_y*ec_y;
	else
	  Q_yy[iph] = 
	    1.0 + q_tau_y[iph]*q_tau_y[iph] + q_p_y[iph]*q_p_y[iph] 
	    - 2.0 * es_y_p[iph]*q_tau_y[iph]
	    - 2.0 * ec_y*q_p_y[iph]
	    + es_dot_ec[iph] * (2.0 * q_tau_y[iph] * q_p_y[iph]);
#pragma ivdep
      for(unsigned iph=0; iph<nph; iph++)
	if(one_minus_es_dot_ec_squared[iph] < DBL_EPSILON)
	  Q_zz[iph] = 
	    1.0 - ec_z*ec_z;
	else
	  Q_zz[iph] = 
	    1.0 + q_tau_z[iph]*q_tau_z[iph] + q_p_z[iph]*q_p_z[iph] 
	    - 2.0 * es_z_p[iph]*q_tau_z[iph]
	    - 2.0 * ec_z*q_p_z[iph]
		+ es_dot_ec[iph] * (2.0 * q_tau_z[iph] * q_p_z[iph]);
#pragma ivdep
      for(unsigned iph=0; iph<nph; iph++)
	if(one_minus_es_dot_ec_squared[iph] < DBL_EPSILON)
	  Q_xy[iph] = 
		-ec_x*ec_y;
	else
	  Q_xy[iph] = 
	    q_tau_x[iph]*q_tau_y[iph] + q_p_x[iph]*q_p_y[iph] 
	    - es_x_p[iph]*q_tau_y[iph] - es_y_p[iph]*q_tau_x[iph]
	    - ec_x*q_p_y[iph] - ec_y*q_p_x[iph]
	    + es_dot_ec[iph] * (q_tau_x[iph] * q_p_y[iph] +
				q_tau_y[iph] * q_p_x[iph]);
#pragma ivdep
      for(unsigned iph=0; iph<nph; iph++)
	if(one_minus_es_dot_ec_squared[iph] < DBL_EPSILON)
	  Q_xz[iph] = 
	    -ec_x*ec_z;
	else
	  Q_xz[iph] = 
	    q_tau_x[iph]*q_tau_z[iph] + q_p_x[iph]*q_p_z[iph] 
	    - es_x_p[iph]*q_tau_z[iph] - es_z_p[iph]*q_tau_x[iph]
	    - ec_x*q_p_z[iph] - ec_z*q_p_x[iph]
	    + es_dot_ec[iph] * (q_tau_x[iph] * q_p_z[iph] +
				q_tau_z[iph] * q_p_x[iph]);
#pragma ivdep
      for(unsigned iph=0; iph<nph; iph++)
	if(one_minus_es_dot_ec_squared[iph] < DBL_EPSILON)
	  Q_yz[iph] = 
	    -ec_y*ec_z;
	else
	  Q_yz[iph] = 
	    q_tau_y[iph]*q_tau_z[iph] + q_p_y[iph]*q_p_z[iph] 
	    - es_y_p[iph]*q_tau_z[iph] - es_z_p[iph]*q_tau_y[iph]
	    - ec_y*q_p_z[iph] - ec_z*q_p_y[iph]
	    + es_dot_ec[iph] * (q_tau_y[iph] * q_p_z[iph] +
				q_tau_z[iph] * q_p_y[iph]);
      
#if 0
      for(unsigned iph=0; iph<nph; iph++)
	std::cout << "Q: " << iscope << '/' << iph << ' '
		  << Q_xx[iph] << ' ' 
		  << Q_yy[iph] << ' ' 
		  << Q_zz[iph] << ' ' 
		  << Q_xy[iph] << ' ' 
		  << Q_xz[iph] << ' ' 
		  << Q_yz[iph] << std::endl;
#endif
      
      double* Q_scope_xx_p = Q_scope_xx+iscope;
      double* Q_scope_yy_p = Q_scope_yy+iscope;
      double* Q_scope_zz_p = Q_scope_zz+iscope;
      double* Q_scope_xy_p = Q_scope_xy+iscope;
      double* Q_scope_xz_p = Q_scope_xz+iscope;
      double* Q_scope_yz_p = Q_scope_yz+iscope;
      
      *Q_scope_xx_p=0;
      *Q_scope_yy_p=0;
      *Q_scope_zz_p=0;
      *Q_scope_xy_p=0;
      *Q_scope_xz_p=0;
      *Q_scope_yz_p=0;
      
#pragma ivdep
      for(unsigned iph=0; iph<nph; iph++)*Q_scope_xx_p += ws_p[iph]*Q_xx[iph];
#pragma ivdep
      for(unsigned iph=0; iph<nph; iph++)*Q_scope_yy_p += ws_p[iph]*Q_yy[iph];
#pragma ivdep
      for(unsigned iph=0; iph<nph; iph++)*Q_scope_zz_p += ws_p[iph]*Q_zz[iph];
#pragma ivdep
      for(unsigned iph=0; iph<nph; iph++)*Q_scope_xy_p += ws_p[iph]*Q_xy[iph];
#pragma ivdep
      for(unsigned iph=0; iph<nph; iph++)*Q_scope_xz_p += ws_p[iph]*Q_xz[iph];
#pragma ivdep
      for(unsigned iph=0; iph<nph; iph++)*Q_scope_yz_p += ws_p[iph]*Q_yz[iph];
    }

  unsigned N = 0;
  double S = 0;
  double Vx = 0;
  double Vy = 0;
  double Vz = 0;
  double Txx = 0;
  double Tyy = 0;
  double Tzz = 0;
  double Txy = 0;
  double Txz = 0;
  double Tyz = 0;
  
#if 0
  for(unsigned iscope=0; iscope<nscopes; iscope++)
    std::cout << "Tel: " << iscope << ' '
	      << Q_scope_xx[iscope] << ' ' 
	      << Q_scope_yy[iscope] << ' ' 
	      << Q_scope_zz[iscope] << ' ' 
	      << Q_scope_xy[iscope] << ' ' 
	      << Q_scope_xz[iscope] << ' ' 
	      << Q_scope_yz[iscope] << std::endl;
#endif
  
#pragma ivdep
  for(unsigned iscope=0; iscope<nscopes;iscope++)
    N += scope_nray[iscope];
  
#pragma ivdep
  for(unsigned iscope=0; iscope<nscopes;iscope++)
    S += 
      r_x[iscope] * Q_scope_xx[iscope] * r_x[iscope] +
      r_y[iscope] * Q_scope_yy[iscope] * r_y[iscope] +
      r_z[iscope] * Q_scope_zz[iscope] * r_z[iscope] +
      2.0 * r_x[iscope] * Q_scope_xy[iscope] * r_y[iscope] +
      2.0 * r_x[iscope] * Q_scope_xz[iscope] * r_z[iscope] +
      2.0 * r_y[iscope] * Q_scope_yz[iscope] * r_z[iscope];
  
#pragma ivdep
  for(unsigned iscope=0; iscope<nscopes;iscope++)
    Vx += 
      Q_scope_xx[iscope] * r_x[iscope] +
      Q_scope_xy[iscope] * r_y[iscope] +
      Q_scope_xz[iscope] * r_z[iscope];
  
#pragma ivdep
  for(unsigned iscope=0; iscope<nscopes;iscope++)
    Vy += 
      Q_scope_xy[iscope] * r_x[iscope] +
      Q_scope_yy[iscope] * r_y[iscope] +
      Q_scope_yz[iscope] * r_z[iscope];
  
#pragma ivdep
  for(unsigned iscope=0; iscope<nscopes;iscope++)
    Vz += 
      Q_scope_xz[iscope] * r_x[iscope] +
      Q_scope_yz[iscope] * r_y[iscope] +
      Q_scope_zz[iscope] * r_z[iscope];
  
#pragma ivdep
  for(unsigned iscope=0; iscope<nscopes;iscope++)Txx += Q_scope_xx[iscope];
#pragma ivdep
  for(unsigned iscope=0; iscope<nscopes;iscope++)Tyy += Q_scope_yy[iscope];
#pragma ivdep
  for(unsigned iscope=0; iscope<nscopes;iscope++)Tzz += Q_scope_zz[iscope];
#pragma ivdep
  for(unsigned iscope=0; iscope<nscopes;iscope++)Txy += Q_scope_xy[iscope];
#pragma ivdep
  for(unsigned iscope=0; iscope<nscopes;iscope++)Txz += Q_scope_xz[iscope];
#pragma ivdep
  for(unsigned iscope=0; iscope<nscopes;iscope++)Tyz += Q_scope_yz[iscope];

  if(minimize_shower_core_location)
    {
      // It can be shown that (V,ec)=0 and hence that ec is an
      // eigenvector of T with eigenvector of zero, so T cannot be
      // inverted. So transform to new basis when system becomes two
      // dimensional: 
      
      // ev = V/|V| 
      // et = ec^ev
      // ec = ec
      
      double Vnorm = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
      
      // Unitary transformation matrix
      double Uvx = Vx/Vnorm;
      double Uvy = Vy/Vnorm;
      double Uvz = Vz/Vnorm;
      double Ucx = ec_x;
      double Ucy = ec_y;
      double Ucz = ec_z;
      double Utx = Ucy*Uvz-Ucz*Uvy;
      double Uty = Ucz*Uvx-Ucx*Uvz;
      double Utz = Ucx*Uvy-Ucy*Uvx;
      
      // Transform T to new coordinates - Trot = U T U'
      
      // First calculate: Temp = T U'
      // perl -e 'for($i=0;$i<3;$i++){ for($j=0;$j<3;$j++) { printf("double Temp_%c%c = ",ord("x")+$i,ord("x")+$j); for($n=0;$n<3;$n++) { printf("T%c%c*U%c%c",ord("x")+(($i<=$n)?$i:$n),ord("x")+(($i<=$n)?$n:$i),ord("x")+$j,ord("x")+$n); if($n==2){print ";"}else{print " + ";}} printf("\n"); } }'
  
      double Temp_xv = Txx*Uvx + Txy*Uvy + Txz*Uvz;
      double Temp_xt = Txx*Utx + Txy*Uty + Txz*Utz;
//    double Temp_xc = Txx*Ucx + Txy*Ucy + Txz*Ucz;
      double Temp_yv = Txy*Uvx + Tyy*Uvy + Tyz*Uvz;
      double Temp_yt = Txy*Utx + Tyy*Uty + Tyz*Utz;
//    double Temp_yc = Txy*Ucx + Tyy*Ucy + Tyz*Ucz;
      double Temp_zv = Txz*Uvx + Tyz*Uvy + Tzz*Uvz;
      double Temp_zt = Txz*Utx + Tyz*Uty + Tzz*Utz;
//    double Temp_zc = Txz*Ucx + Tyz*Ucy + Tzz*Ucz;
  
      // Second calculate Trot = U Temp
      // perl -e 'for($i=0;$i<3;$i++){ for($j=0;$j<3;$j++) { printf("double T%c%c = ",ord("x")+$i,ord("x")+$j); for($n=0;$n<3;$n++) { printf("U%c%c*Temp_%c%c",ord("x")+$i,ord("x")+$n,ord("x")+$n,ord("x")+$j); if($n==2){print ";"}else{print " + ";}} printf("\n"); } }'
      
      double Tvv = Uvx*Temp_xv + Uvy*Temp_yv + Uvz*Temp_zv;
      double Tvt = Uvx*Temp_xt + Uvy*Temp_yt + Uvz*Temp_zt;
//    double Tvc = Uvx*Temp_xc + Uvy*Temp_yc + Uvz*Temp_zc;
//    double Ttv = Utx*Temp_xv + Uty*Temp_yv + Utz*Temp_zv;
      double Ttt = Utx*Temp_xt + Uty*Temp_yt + Utz*Temp_zt;
//    double Ttc = Utx*Temp_xc + Uty*Temp_yc + Utz*Temp_zc;
//    double Tcv = Ucx*Temp_xv + Ucy*Temp_yv + Ucz*Temp_zv;
//    double Tct = Ucx*Temp_xt + Ucy*Temp_yt + Ucz*Temp_zt;
//    double Tcc = Ucx*Temp_xc + Ucy*Temp_yc + Ucz*Temp_zc;

      // Transform V to new coordinates, easy to do in this basis 
//    double Vv = Vnorm;
//    double Vt = 0;
//    double Vc = 0;

      // Calculate determinant of 2x2 upper-left portion of Trot
      double Tdet = Tvv*Ttt - Tvt*Tvt;
      
      // Calculate core location: rc_rot = Trot(inv) Vrot
      double rc_v = Vnorm*Ttt/Tdet;
      double rc_t = -Vnorm*Tvt/Tdet;
//    double rc_c = 0;

      // Recover core location in orignal coordinate system: rc = U' rc_rot
      rc_x = Uvx*rc_v + Utx*rc_t;
      rc_y = Uvy*rc_v + Uty*rc_t;
      rc_z = Uvz*rc_v + Utz*rc_t;

      chi2 = S - ( Vx*rc_x + Vy*rc_y + Vz*rc_z );
    }
  else
    {
      // Watch out -- if you do not minimize rc then you cannot use
      // simplified value for chi2 (above)

      double Tx = Txx*rc_x + Txy*rc_y + Txz*rc_z;
      double Ty = Txy*rc_x + Tyy*rc_y + Tyz*rc_z;
      double Tz = Txz*rc_x + Tyz*rc_y + Tzz*rc_z;

      chi2 = 
	S - 2.0*(Vx*rc_x + Vy*rc_y + Vz*rc_z) + (Tx*rc_x + Ty*rc_y + Tz*rc_z);
    }
  
#if 0
  std::cout << "Theta:" << distanceToAxis(optical_axis, 
					  ec_x, ec_y, ec_z)/M_PI*180
	    << std::endl;
  
  std::cout << std::endl;
  
  std::cout << "T:    " << Txx << ' ' << Txy << ' ' << Txz << std::endl
	    << "      " << Txy << ' ' << Tyy << ' ' << Tyz << std::endl
	    << "      " << Txz << ' ' << Tyz << ' ' << Tzz << std::endl;
  
  std::cout << std::endl;
  
  std::cout << "U:    " << Uvx << ' ' << Uvy << ' ' << Uvz << std::endl
	    << "      " << Utx << ' ' << Uty << ' ' << Utz << std::endl
	    << "      " << Ucx << ' ' << Ucy << ' ' << Ucz << std::endl;
  
  std::cout << std::endl;
  
  std::cout << "Trot: " << Tvv << ' ' << Tvt << ' ' << Tvc << std::endl
	    << "      " << Tvt << ' ' << Ttt << ' ' << Ttc << std::endl
	    << "      " << Tvc << ' ' << Ttc << ' ' << Tcc << std::endl;
  
  std::cout << std::endl;
  
  std::cout << "Tdet: " << Tdet << std::endl;
  
  std::cout << std::endl;
  
  std::cout << "V:    " << Vx << std::endl
	    << "      " << Vy << std::endl
	    << "      " << Vz << std::endl;
  
  std::cout << std::endl;
  
  std::cout << "Vv: " << Vnorm << std::endl;
  
  std::cout << std::endl;
  
  std::cout << "rc_v: " << rc_v << std::endl
	    << "  _t  " << rc_t << std::endl;
  
  std::cout << std::endl;
  
  std::cout << "rc:   " << rc_x << std::endl
	    << "      " << rc_y << std::endl
	    << "      " << rc_z << std::endl;
  
  std::cout << std::endl;
  
  std::cout << "S:    " << S << std::endl;
  
  std::cout << std::endl;
  
  std::cout << "chi2: " << chi2 << std::endl;
  
  std::cout << std::endl;
#endif
}
