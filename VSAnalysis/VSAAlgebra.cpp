//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAAlgebra.cpp
  Various linear algebra elements

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       11/08/2005
*/

#include<cassert>
#include<stdexcept>
#include<cmath>
#include<cfloat>

#include<VSAMath.hpp>
#include<VSAAlgebra.hpp>

using namespace VERITAS;
using namespace VERITAS::VSAAlgebra;

//#define DEBUG
#ifdef DEBUG
#include<iostream>
#endif

// ============================================================================
// ============================================================================
// 
// Vec2D
//
// ============================================================================
// ============================================================================

void Vec2D::rotate(double theta)
{
  const double c = cos(theta);
  const double s = sin(theta);
  set(c*fCx-s*fCy, c*fCy+s*fCx);
}

void Vec2D::rotate(const Vec2D& point, double theta)
{
  set(x()-point.x(), y()-point.y());
  rotate(theta);
  set(x()+point.x(), y()+point.y());
}

void Vec2D::perp(Vec2D& p1) const
{
  p1.setNormalized(-y(),x());
}

// ============================================================================
// ============================================================================
// 
// Vec3D
//
// ============================================================================
// ============================================================================

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/// Operator &= of rotation composition
/*! \note 
  Composition is opposite to the sense of matrix
  multiplication. The composition r=r1&r2 is equivalent to 
  a rotation first by r1 then by r2.

  For example, for two rotations r1 and r2, and for any vector p

  Vec3D r1(?,?,?);
  Vec3D r2(?,?,?);
  Vec3D p(?,?,?);

  Vec3D r=r1&r2;

  p.Rotate(r1);   \/  is equivalent to p.Rotate(r)
  p.Rotate(r2);   /\ 
*/

/*
   Composition algorithm comes from consideration of rotation
   as when expressed as spin matrices... i.e. write the both 
   rotations in terms of Pauli matrices, multiply and then compare 
   terms. See any QM book or Goldstein chap 4.

   $\mathbold{R_1} = \mathbold{1} \cos\theta_1 / 2 +
    i \hat{n}_1 . \vec{\mathbold{sigma}} \sin\theta_1 / 2 $

   $\mathbold{R_2} = \mathbold{1} \cos\theta_2 / 2 +
    i \hat{n}_2 . \vec{\mathbold{sigma}} \sin\theta_2 / 2 $

    Multiply the matrices using the following, collect terms
    and compare with $R_1$ or $R_2$ to deduce the composite
    angle and axis of rotation.

   $(\vec{\mathbold{sigma}}.\vec{A})
    (\vec{\mathbold{sigma}}.\vec{B}) = \mathbold{1}\vec{A}.\vec{B}+
    i\vec{\mathbold{sigma}}(\vec{A}\cross\vec{B})$
*/

Vec3D& Vec3D::operator &= (const Vec3D& r)
{
#define R1 (*this)
#define R2 r

  const double r1_theta = R1.norm();
  const double r2_theta = R2.norm();

  if(r1_theta == 0)
    {
      // R1 is zero so set to R2
      *this = R2;
    }
  else if(r2_theta == 0)
    {
      // R2 is zero so set to R1
      *this = R1;
    }
  else
    {
      const double sin_half_r1 = sin(r1_theta/2.0);
      const double cos_half_r1 = cos(r1_theta/2.0);
      const double sin_half_r2 = sin(r2_theta/2.0);
      const double cos_half_r2 = cos(r2_theta/2.0);

      const Vec3D r1_hat = R1/r1_theta;
      const Vec3D r2_hat = R2/r2_theta;

      const double coshalftheta =
        cos_half_r1*cos_half_r2-
        sin_half_r1*sin_half_r2*(r1_hat*r2_hat);

      const Vec3D axis(r1_hat*sin_half_r1*cos_half_r2+
		       r2_hat*sin_half_r2*cos_half_r1-
		       (r1_hat^r2_hat)*sin_half_r1*sin_half_r2);

      const double sinhalftheta = axis.norm();
      if(sinhalftheta==0)
	{
	  *this = Vec3D(0,0,0);
	}
      else
	{
	  const double halftheta=atan(sinhalftheta/coshalftheta);
	  *this = axis/sinhalftheta*halftheta*2.0;
	}
    }
  return *this;
}

/** 
    Rotate vector around an arbitrary axis.
    \param axis: axis of rotation vector [rad]
    \note The modulus of rotation angle given in radians is 
    equal to the norm of the axis vector. The rotation of e_x
    around e_z by PI/2 is equal to e_y.
*/

void Vec3D::rotate(const Vec3D& axis)
{
  double angle = axis.norm();
  if(angle == 0)return;

  Vec3D e_axis(axis/angle);

  double par   = (*this)*e_axis;
  Vec3D p_axis(*this - par*e_axis);

  // "Rotation Formula" -- see Goldstein chap 4
  *this = par*e_axis + cos(angle)*p_axis + sin(angle)*(e_axis^p_axis);
}

void Vec3D::rotate(const RotationVec3D& axis)
{
  double par = (*this)*axis.axis();
  Vec3D p_axis(*this - par*axis.axis());
  *this = ( par*axis.axis() 
	    + axis.cosTheta()*p_axis 
	    + axis.sinTheta()*(axis.axis()^p_axis) );
}

/**
   Return two orthogonal vectors which define the plane perpendicular to
   this vector
**/

void Vec3D::perp(Vec3D& p1, Vec3D& p2) const
{
  const double cx = fCx;
  const double cy = fCy;
  const double cz = fCz;

  double min;
  unsigned cpt;

  // Find smallest component;
  min = fabs(cx), cpt = 0;
  if(fabs(cy) < min)min = fabs(cy), cpt=1;
  if(fabs(cz) < min)cpt = 2;

  // The first vector is formed by taking a cross product of
  // (cx,cy,cz) and the unit vector along the direction of smallest
  // component of c

  // The second vector is the cross product of (cx,cy,cz) and the first

  if(cpt==0) 
    { 
      // c1y = cz;  c1z = -cy; 
      p1.setNormalized(0,cz,-cy);
      p2.setNormalized(-cy*cy-cz*cz,cx*cy,cx*cz);
    }
  else if(cpt==1)
    { 
      // c1x = -cz; c1z = cx;  
      p1.setNormalized(-cz,0,cx);
      p2.setNormalized(cy*cx,-cz*cz-cx*cx,cy*cz);
    }
  else
    { 
      // c1x = cy;  c1y = -cx; 
      p1.setNormalized(cy,-cx,0);
      p2.setNormalized(cz*cx,cz*cy,-cx*cx-cy*cy);
    }

  //p1 = Vec3D(c1x,c1y,c1z);
  //p2 = Vec3D(cy*c1z-cz*c1y, cz*c1x-cx*c1z, cx*c1y-cy*c1x);
}

/**
   Scatter the direction of the vector with 2-dimensional Gaussian
   probability distribution of RMS (in the x and y directions) of
   sigma. Two uniform random deviates between 0 and 1 must be
   supplied.
**/

void Vec3D::scatter(double sigma, double uniform_dev_1, double uniform_dev_2)
{
  if(sigma>0)
    {
      const double vtheta = theta();
      const double vphi   = phi();
      const double vnorm  = norm();

      const double stheta = sigma*sqrt(-2*log(uniform_dev_1));
      const double sphi   = 2.0*M_PI*uniform_dev_2;
      const double cstheta = cos(stheta);
      const double sstheta = sin(stheta);

      fCx = vnorm*cos(sphi)*sstheta;
      fCy = vnorm*sin(sphi)*sstheta;
      fCz = vnorm*cstheta;
      
      rotate(Vec3D(0,vtheta,0));
      rotate(Vec3D(0,0,vphi));
    }
}

Vec3D Vec3D::makeRotation(const Vec3D& v1, const Vec3D& v2, 
			  Subspace2D subspace)
{
  Vec3D ea(v1);
  ea.normalize();
  Vec3D eb(v2);
  eb -= (ea*eb)*ea;
  eb.normalize();
  Vec3D ec;

  switch(subspace)
    {
    case SS_12:
      ec = ea^eb;
      break;
    case SS_23:
      ec = eb;
      eb = ea;
      ea = eb^ec;
      break;
    case SS_31:
      ec = ea;
      ea = eb;
      eb = ec^ea;
      break;
    }

  Vec3D r2 = Vec3D(0,0,1)^ec;
  double sin_theta = r2.norm();
  double cos_theta = ec.z();
  if(sin_theta!=0)r2 *= atan2(sin_theta,cos_theta)/sin_theta;

  Vec3D r2t_ea(ea); r2t_ea.rotate(-r2);
  //Vec r2t_eb(eb); r2t_eb.rotate(-r2);

  double sin_phi = r2t_ea.y();
  double cos_phi = r2t_ea.x();
  Vec3D r1 = Vec3D(0,0,atan2(sin_phi,cos_phi));

  Vec3D rc(r1); rc &= r2;
  return rc;
}

// ============================================================================
// ============================================================================
// 
// Symmetric2D
// 
// ============================================================================
// ============================================================================

#define DBL_INV_EPSILON 1.0/DBL_EPSILON

#define ASSIGN2(A,B) \
  e.val[1-A] = a11; e.vec[1-A] = Vec2D(u11,u21); \
  e.val[1-B] = a22; e.vec[1-B] = Vec2D(u12,u22)

unsigned Symmetric2D::eigen(Eigen2D& e) const
{
  if(fabs(fA12) <= DBL_EPSILON)
    {
      unsigned i = (fA11>fA22) ? 1 : 0;
      e.val[i] = fA11; e.vec[i] = Vec2D(1,0);
      i = 1-i;
      e.val[i] = fA22; e.vec[i] = Vec2D(0,1);
    }
  else
    {
      static const double thetamax = sqrt(DBL_INV_EPSILON);
  
      const double h = (fA22-fA11);
      const double theta = h/(2.0*fA12);
      const double fabs_theta = fabs(theta);
      const double t =
	(fabs_theta >= thetamax)?
	0.5/theta:
	((theta>0)?
	 1.0/(fabs(theta)+sqrt(1.0+theta*theta)):
	 -1.0/(fabs(theta)+sqrt(1.0+theta*theta)));
      const double c = 1.0/sqrt(1.0+t*t);
      const double s = t*c;
      //const double tau = s/(1.0+c);
  
      const double ta12 = t*fA12;

#ifdef DEBUG
      std::cerr << "Theta: " << atan(1/theta)/M_PI*180/2 << std::endl;
#endif

      const double a11 = fA11 - ta12;
      const double a22 = fA22 + ta12;

      const double u11 = c; //1 - s*tau;
      const double u21 = -s;
      const double u12 = s;
      const double u22 = c; // u + s*tau;

      if(a11>a22){ ASSIGN2(0,1); }
      else { ASSIGN2(1,0); }
    }

  return 2;
}

unsigned Symmetric2D::eigenOld(Eigen2D& e) const
{
  unsigned nroots = VSAMath::realRoots2(-1,trace(),-det(),e.val);

  // Regularize nearly equal and nearly zero eigenvalues (see SVD
  // section in Numerical Recipies)

  double scale = maxAbsEl();
  double split = fabs(e.val[1]-e.val[0])/scale;
  if((split>0)&&(split<DBL_EPSILON))
    {
      e.val[1]=e.val[0];
      nroots--;
    }
  
  if(fabs(e.val[0])/scale < DBL_EPSILON)e.val[0]=0;
  if(fabs(e.val[1])/scale < DBL_EPSILON)e.val[1]=0;

  if(nroots==2)
    {
      for(unsigned i=0; i<3; i++)
	{
	  const double l = e.val[i];
	  
	  const double a11 = fA11-l;
	  const double a22 = fA22-l;
	  const double a12 = fA12;

	  bool zero1 = ((a11!=0)&&(a12==0))||((a12!=0)&&(a22==0));
	  bool zero2 = ((a12!=0)&&(a11==0))||((a22!=0)&&(a12==0));

#ifdef DEBUG
	  std::cerr << std::fixed;
	  std::cerr << a11 << '\t' << a12 << std::endl
		    << a12 << '\t' << a22 << std::endl;
	  
	  std::cerr << zero1 << zero2 << std::endl;
#endif		
	  
	  if(zero1)
	    e.vec[i].set(0.0,1.0);
	  else if(zero2)
	    e.vec[i].set(1.0,0.0);
	  else
	    e.vec[i].setNormalized(1,-a11/a12);
	}
    }
  else
    {
      // nroots == 1 since all symmetric 2D matrices have at least one e-vec
      e.vec[0].set(1,0);
      e.vec[1].set(0,1);
    }
  return nroots;
}

const Symmetric2D Symmetric2D::kronecker(1.0,1.0,0.0);

// ============================================================================
// ============================================================================
// 
// Symmetric3D: Definitions of longer functions
//
// ============================================================================
// ============================================================================

static inline
void rotateUST(const double s, const double tau,
	       double& u11, double& u21, double& u31, 
	       double& u12, double& u22, double& u32) 
{
  const double _u11 = u11;
  const double _u21 = u21;
  const double _u31 = u31;
  const double _u12 = u12;
  const double _u22 = u22;
  const double _u32 = u32;

  u11 -= s*(_u12 + tau*_u11);
  u21 -= s*(_u22 + tau*_u21);
  u31 -= s*(_u32 + tau*_u31);
  u12 += s*(_u11 - tau*_u12);
  u22 += s*(_u21 - tau*_u22);
  u32 += s*(_u31 - tau*_u32);
}

static inline
void rotateUSC(const double s, const double c,
	       double& u11, double& u21, double& u31, 
	       double& u12, double& u22, double& u32) 
{
  const double _u11 = u11;
  const double _u21 = u21;
  const double _u31 = u31;
  const double _u12 = u12;
  const double _u22 = u22;
  const double _u32 = u32;

  u11 = c*_u11 - s*_u12;
  u21 = c*_u21 - s*_u22;
  u31 = c*_u31 - s*_u32;
  u12 = c*_u12 + s*_u11;
  u22 = c*_u22 + s*_u21;
  u32 = c*_u32 + s*_u31;
}

static inline
void jacobi3(double& a11, double& a22, double& a12, double& a13, double& a23,
	     double& u11, double& u21, double& u31, 
	     double& u12, double& u22, double& u32)
{
  // Eigenvalue/eignevector calculation using Jacobi rotation method.
  // See numerical recipes, chapter 11, section 2

  // if (theta > thetamax) then (theta^2 > 1/eps) and 
  // (theta^2+1 = theta^2 ( 1 + 1/theta^2 ) = theta^2 )
  static const double thetamax = sqrt(DBL_INV_EPSILON);

  const double h = (a22-a11);
  const double theta = h/(2.0*a12);
  const double fabs_theta = fabs(theta);
  const double t =
    (fabs_theta >= thetamax)?
    0.5/theta:
    ((theta>0)?
     1.0/(fabs(theta)+sqrt(1.0+theta*theta)):
     -1.0/(fabs(theta)+sqrt(1.0+theta*theta)));
  const double c = 1.0/sqrt(1.0+t*t);
  const double s = t*c;
  const double tau = s/(1.0+c);

  const double ta12 = t*a12;

  const double _a13 = a13;
  const double _a23 = a23;

#ifdef DEBUG
  std::cerr << "Theta: " << atan(1/theta)/M_PI*180/2 << std::endl;
#endif

  a12 = 0;
  a11 -= ta12;
  a22 += ta12;
  a13 -= s*(_a23 + tau*_a13);
  a23 += s*(_a13 - tau*_a23);

  rotateUST(s, tau, u11, u21, u31, u12, u22, u32);
}

#define ASSIGN3E(A,B,C) \
  e.val[2-A] = a11; e.vec[2-A] = Vec3D(u11,u21,u31); \
  e.val[2-B] = a22; e.vec[2-B] = Vec3D(u12,u22,u32); \
  e.val[2-C] = a33; e.vec[2-C] = Vec3D(-u13,-u23,-u33)

#define ASSIGN3O(A,B,C) \
  e.val[2-A] = a11; e.vec[2-A] = Vec3D(u11,u21,u31); \
  e.val[2-B] = a22; e.vec[2-B] = Vec3D(u12,u22,u32); \
  e.val[2-C] = a33; e.vec[2-C] = Vec3D(u13,u23,u33)

unsigned Symmetric3D::eigenJacobi(Eigen3D& e) const
{
  // Eigenvalue/eignevector calculation using Jacobi rotation method.
  // See numerical recipes, chapter 11, section 2

  double u11 = 1;
  double u12 = 0;
  double u13 = 0;
  double u21 = 0;
  double u22 = 1;
  double u23 = 0;
  double u31 = 0;
  double u32 = 0;
  double u33 = 1;

  double a11 = fA11;
  double a22 = fA22;
  double a33 = fA33;
  double a12 = fA12;
  double a23 = fA23;
  double a13 = fA13;

  for(unsigned i=0; i<50; i++)
    {
      const double sum = fabs(a12) + fabs(a23) + fabs(a13);
      if(sum == 0)break;

      const double threshold = (i<3)?0.0222*sum:0;
      
      double cutoff;
      
      // ##### 1-2 #####

#ifdef DEBUG
      Symmetric3D m1(a11,a22,a33,a12,a13,a23);
      std::cerr << m1 << std::endl << std::endl;
#endif

      cutoff = 100.0*DBL_INV_EPSILON*fabs(a12);
      if((i>=3)&&(fabs(a11) >= cutoff)&&(fabs(a22) >= cutoff))
	a12=0;
      else if(fabs(a12) > threshold)
	jacobi3(a11, a22, a12, a13, a23, u11, u21, u31, u12, u22, u32);

      // ##### 2-3 #####

#ifdef DEBUG
      Symmetric3D m2(a11,a22,a33,a12,a13,a23);
      std::cerr << m2 << std::endl << std::endl;
#endif

      cutoff = 100.0*DBL_INV_EPSILON*fabs(a23);
      if((i>=3)&&(fabs(a22) >= cutoff)&&(fabs(a33) >= cutoff))
	a23=0;
      else if(fabs(a23) > threshold)
	jacobi3(a22, a33, a23, a12, a13, u22, u32, u12, u23, u33, u13);
      
      // ##### 3-1 #####

#ifdef DEBUG
      Symmetric3D m3(a11,a22,a33,a12,a13,a23);
      std::cerr << m3 << std::endl << std::endl;
#endif

      cutoff = 100.0*DBL_INV_EPSILON*fabs(a13);;
      if((i>=3)&&(fabs(a33) >= cutoff)&&(fabs(a11) >= cutoff))
	a13=0;
      else if(fabs(a13) > threshold)
	jacobi3(a33, a11, a13, a23, a12, u33, u13, u23, u31, u11, u21);
    }

#ifdef VSA_ASSERT  

#ifdef DEBUG
  std::cerr << fabs(a12) << " + " << fabs(a23) << " + " << fabs(a13) << " = "
	    << fabs(a12) + fabs(a23) + fabs(a13) << std::endl;
#endif

  assert(fabs(a12) + fabs(a23) + fabs(a13) == 0);
#endif // VSA_ASSERT

#ifdef DEBUG
  Symmetric3D m(a11,a22,a33,a12,a13,a23);
  std::cerr << m << std::endl << std::endl;
#endif

  if(a11>a22)
    if(a11>a33)
      if(a22>a33)
	{ ASSIGN3E(0,1,2); }
      else
	{ ASSIGN3O(0,2,1); }
    else
      { ASSIGN3E(1,2,0); }
  else
    if(a22>a33)
      if(a11>a33)
	{ ASSIGN3O(1,0,2); }
      else
	{ ASSIGN3E(2,0,1); }
    else
      { ASSIGN3O(2,1,0); }
  
  return 3;
}

static inline double SQR(const double x)
{
  return x*x;
}

static inline double pythag(const double a, const double b)
{
  double fabs_a = fabs(a);
  double fabs_b = fabs(b);
  if(fabs_a > fabs_b) return fabs_a*sqrt(1.0+SQR(fabs_b/fabs_a));
  else return (fabs_b == 0.0 ? 0.0 : fabs_b*sqrt(1.0+SQR(fabs_a/fabs_b)));
}

static inline double SIGN(const double a, const double b)
{
  return (b>=0.0)?fabs(a):-fabs(a);
}

unsigned Symmetric3D::eigenQL(Eigen3D& e) const
{
  // Eigenvalue/eignevector calculation using QL method.
  // See numerical recipes, chapter 11, section 3

  // Strategy:
  // 1) On a 3-by-3 a Jacobi rotation in the 3-1 plane brings the matrix 
  //    into tri-diagonal form, suitable for QL
  // 2) Implement QL until one of the remaining off-diagonals is below 
  //    threshold
  // 3) One final Jacobi rotation in 1-2 or 2-3 plane, depending on results
  //    of step 2, to suppress final off-diagonal
  
  double u11 = 1;
  double u12 = 0;
  double u13 = 0;
  double u21 = 0;
  double u22 = 1;
  double u23 = 0;
  double u31 = 0;
  double u32 = 0;
  double u33 = 1;

  double a11 = fA11;
  double a22 = fA22;
  double a33 = fA33;
  double a12 = fA12;
  double a23 = fA23;
  double a13 = fA13;

  // Jacobi rotation to emiminate a13 -- matrix is now tri-diagonal

  double cutoff = 100.0*DBL_INV_EPSILON*fabs(a13);;
  if((fabs(a33) >= cutoff)&&(fabs(a11) >= cutoff))
    a13=0;
  else 
    jacobi3(a33, a11, a13, a23, a12, u33, u13, u23, u31, u11, u21);

  unsigned m;
  while(1)
    {
      //Symmetric3D m1(a11,a22,a33,a12,a13,a23);
      //std::cerr << m1 << std::endl << std::endl;

      // Test whether a12 or a23 is small enough. If so then exit the
      // iteration with the variable "m" set appropriately

      if((fabs(a11)+fabs(a22))>=DBL_INV_EPSILON*fabs(a12))
	{
	  m=0;
	  a12 = 0;
	  break;
	}
      else if((fabs(a22)+fabs(a33))>=DBL_INV_EPSILON*fabs(a23))
	{
	  m=1;
	  a23 = 0;
	  break;
	}
      else m=2;

      // Since the "i" loop in the NR algorithm is only over two values
      // in the 3-by-3 case, I have just unrolled it and eliminated the
      // unnecessary calculations from the first go through

      double g = (a22-a11)/(2.0*a12);
      double r = pythag(g,1.0);
      g=a33-a11+a12/(g+SIGN(r,g));

      // ===== I = 1 =====

      r = pythag(a23,g);
      if(r == 0.0)
	{
	  // what happens here -- how can i ever break out of the
	  // loop in this case?
	  assert(0);
	  continue;
	}

      double s = a23/r;
      double c = g/r;
      r = (a22-a33)*s+2.0*c*a23;
      double p = s*r;
      a33 += p;
      g = c*r-a23;

      rotateUSC(s, c, u12, u22, u32, u13, u23, u33);

      // ===== I = 0 =====

      double f = s*a12;
      double b = c*a12;
      r = pythag(f,g);
      a23 = r;

      if(r == 0.0)
	{
	  // this is ok i think since the i=1 will have modified things
	  a22 -= p;
	  continue;
	}

      s = f/r;
      c = g/r;
      g = a22-p;
      r = (a11-g)*s+2.0*c*b;
      p = s*r;
      a22 = g+p;
      g = c*r-b;

      rotateUSC(s, c, u11, u21, u31, u12, u22, u32);

      a11 -= p;
      a12 = g;
    }

  if(m==0)
    {
      cutoff = 100.0*DBL_INV_EPSILON*fabs(a23);
      if((fabs(a22) >= cutoff)&&(fabs(a33) >= cutoff))
	a23=0; // not strictly necessary since we do not use a23 any more
      else
	jacobi3(a22, a33, a23, a12, a13, u22, u32, u12, u23, u33, u13);
    }
  else if(m==1)
    {
      cutoff = 100.0*DBL_INV_EPSILON*fabs(a12);
      if((fabs(a11) >= cutoff)&&(fabs(a22) >= cutoff))
	a12=0; // not strictly necessary since we do not use a12 any more
      else
	jacobi3(a11, a22, a12, a13, a23, u11, u21, u31, u12, u22, u32);
    }
  else
    {
      // this never happens (hopefully!)
      assert(0);
    }

  if(a11>a22)
    if(a11>a33)
      if(a22>a33)
	{ ASSIGN3E(0,1,2); }
      else
	{ ASSIGN3O(0,2,1); }
    else
      { ASSIGN3E(1,2,0); }
  else
    if(a22>a33)
      if(a11>a33)
	{ ASSIGN3O(1,0,2); }
      else
	{ ASSIGN3E(2,0,1); }
    else
      { ASSIGN3O(2,1,0); }
  
  return 3;
}

const Symmetric3D Symmetric3D::kronecker(1.0,1.0,1.0,0.0,0.0,0.0);

// ============================================================================
// 
// Orthogonal3D: Definitions of longer functions
//
// ============================================================================

Orthogonal3D::Orthogonal3D(const Vec3D& v1, const Vec3D& v2, 
			   Subspace2D subspace):
  fRotation(),
  fAX(0), fAY(0), fAZ(0),
  fBX(0), fBY(1), fBZ(0),
  fCX(0), fCY(0), fCZ(0)
{
  fRotation = Vec3D::makeRotation(v1,v2,subspace);
  calcElements();
}

void Orthogonal3D::calcElements()
{
  Vec3D ex(1,0,0); ex.rotate(fRotation);
  Vec3D ey(0,1,0); ey.rotate(fRotation);
  Vec3D ez(0,0,1); ez.rotate(fRotation);

#ifdef DEYUG
  std::cerr << "ex: " << ex << std::endl;
  std::cerr << "ey: " << ey << std::endl;
  std::cerr << "ez: " << ez << std::endl;
#endif

  fAX = ex.a();  fAY = ey.a();  fAZ = ez.a();
  fBX = ex.b();  fBY = ey.b();  fBZ = ez.b();
  fCX = ex.c();  fCY = ey.c();  fCZ = ez.c();
}

// ============================================================================
// 
// MatrixND: Definitions of longer functions
//
// ============================================================================

void MatrixND::gaussJordan(MatrixND& a, MatrixND& b)
{
  const unsigned n = a.nrow();
  const unsigned m = b.ncol();
  assert(a.isSquare());
  assert(b.nrow() == n);

  std::vector<unsigned> indxc(n,0);
  std::vector<unsigned> indxr(n,0);
  std::vector<unsigned> ipiv(n,0);
  for(unsigned i=0;i<n;i++)
    {
      double big = 0.0;
      unsigned irow = 0;
      unsigned icol = 0;
      for(unsigned j=0;j<n;j++)
	if(ipiv[j] != 1)
	  for(unsigned k=0;k<n;k++)
	    if(ipiv[k] == 0)
	      if(fabs(a(j,k))>=big)
		big=fabs(a(j,k)), irow=j, icol=k;
      ipiv[icol]++;
      if(irow != icol)
	{
	  for(unsigned l=0;l<n;l++)std::swap(a(irow,l), a(icol,l));
	  for(unsigned l=0;l<m;l++)std::swap(b(irow,l), b(icol,l));
	}
      indxr[i] = irow;
      indxc[i] = icol;
      if(a(icol,icol) == 0.0)
	throw std::invalid_argument(std::string(__PRETTY_FUNCTION__)
				    + ": singular Matrix");

      double pivinv = 1.0/a(icol,icol);
      a(icol,icol) = 1.0;
      for(unsigned l=0;l<n;l++)a(icol,l) *= pivinv;
      for(unsigned l=0;l<m;l++)b(icol,l) *= pivinv;
      for(unsigned ll=0;ll<n;ll++)
	if(ll != icol)
	  {
	    double temp = a(ll,icol);
	    a(ll,icol) = 0;
	    for(unsigned l=0;l<n;l++)a(ll,l) -= a(icol,l)*temp;
	    for(unsigned l=0;l<m;l++)b(ll,l) -= b(icol,l)*temp;
	  }
    }
  
  for(int l=n-1;l>=0;l--)
    if(indxr[l] != indxc[l])
      for(unsigned k=0;k<n;k++)
	std::swap(a(k,indxr[l]),a(k,indxc[l]));
}

#ifdef TEST_MAIN_1

#include<iostream>
#include<iomanip>
#include<VSAMath.hpp>

#include<RandomNumbers.hpp>

int main(int argc, char**argv)
{
  // Test roots solver
  
  unsigned nroots;

#if 0
  double roots2[2];
  double roots3[3];

  double a,b,c,d;
  double A,B,C;

  // (x-A)*(x-B) = x*x - (A+B)*x + A*B

  A = -17;
  B = 10;
  a = -3;
  b = -a*(A+B);
  c = a*A*B;

  nroots = VSAMath::realRoots2(a,b,c,roots2);
  std::cout << nroots << ' ' << roots2[0] << ' ' << roots2[1] << std::endl;

  A = -5;
  B = -5;
  a = 2;
  b = -a*(A+B);
  c = a*A*B;

  nroots = VSAMath::realRoots2(a,b,c,roots2);
  std::cout << nroots << ' ' << roots2[0] << ' ' << roots2[1] << std::endl;

  // (x-A)*(x-B)*(x-C) = x*x*x - (A+B+C)*x*x + (A*B+A*C+B*C)*x - A*B*C

  A = -17;
  B = 10;
  C = 34;
  a = -2;
  b = -a*(A+B+C);
  c = a*(A*B+A*C+B*C);
  d = -a*A*B*C;

  nroots = VSAMath::realRoots3(a,b,c,d,roots3);
  std::cout << nroots << ' ' 
	    << roots3[0] << ' ' << roots3[1] << ' ' << roots3[2] << std::endl;

  A = 12;
  B = 12;
  C = 17;
  a = 2;
  b = -a*(A+B+C);
  c = a*(A*B+A*C+B*C);
  d = -a*A*B*C;

  nroots = VSAMath::realRoots3(a,b,c,d,roots3);
  std::cout << nroots << ' ' 
	    << roots3[0] << ' ' << roots3[1] << ' ' << roots3[2] << std::endl;

  A = 10;
  B = 10;
  C = 10;
  a = 2;
  b = -a*(A+B+C);
  c = a*(A*B+A*C+B*C);
  d = -a*A*B*C;

  nroots = VSAMath::realRoots3(a,b,c,d,roots3);
  std::cout << nroots << ' ' 
	    << roots3[0] << ' ' << roots3[1] << ' ' << roots3[2] << std::endl;
#endif

  //Symmetric3D matrix(1,2,3,4,5,6);
  //Symmetric3D matrix(1,0,0,0,0,0);
  //Symmetric3D matrix(0,0,0,0,0,0);
  //Symmetric3D matrix(0,0,0,1,1,0);
  //Symmetric3D matrix(0,0,0,1,0,0);
  //Symmetric3D matrix(0,0,0,0,1,0);
  //Symmetric3D matrix(0,0,0,0,0,1);
  //Symmetric3D matrix(1.0625,1.1875,1.75,-0.108253175473055,0.216506350946110,-0.375);

  RandomNumbers rng(RandomNumbers::defaultFilename().c_str());
  Eigen3D e;

#if 0
  for(unsigned i=0;i<10000000;i++)
    {
      Symmetric3D matrix(rng.Uniform(),rng.Uniform(),rng.Uniform(),
			 rng.Uniform(),rng.Uniform(),rng.Uniform());
      //nroots = matrix.eigenJacobi(e);
      nroots = matrix.eigenQL(e);
      //std::cout << e.val[0] << std::endl;      
    }
#endif


#if 1
  Symmetric3D matrix(rng.Uniform(),rng.Uniform(),rng.Uniform(),
		     rng.Uniform(),rng.Uniform(),rng.Uniform());
#else
  Symmetric3D matrix(0.641933, 0.587790, 0.858495, 
		     0.683663, 0.302988, 0.970878);
#endif
  std::cout << std::fixed
	    << matrix.a11() << '\t' << matrix.a12() << '\t' << matrix.a13() 
	    << std::endl
	    << matrix.a12() << '\t' << matrix.a22() << '\t' << matrix.a23() 
	    << std::endl
	    << matrix.a13() << '\t' << matrix.a23() << '\t' << matrix.a33() 
	    << std::endl << std::endl;

  nroots = matrix.eigenJacobi(e);
  std::cout << "N eigenvalues: " << nroots << std::endl;
  for(unsigned i=0;i<3;i++)
    std::cout << "l" << i << " = " << e.val[i]
	      << "\te" << i << " = " 
	      << e.vec[i].x() << '\t' << e.vec[i].y() << '\t' << e.vec[i].z()
	      << std::endl;
  std::cout << std::endl;

  nroots = matrix.eigenQL(e);
  std::cout << "N eigenvalues: " << nroots << std::endl;
  for(unsigned i=0;i<3;i++)
    std::cout << "l" << i << " = " << e.val[i]
	      << "\te" << i << " = " 
	      << e.vec[i].x() << '\t' << e.vec[i].y() << '\t' << e.vec[i].z()
	      << std::endl;

#if 0
  Orthogonal3D u(e.vec[0],e.vec[1]);
  Symmetric3D d = u.transFwd(matrix);
  std::cout << std::fixed
	    << d.a11() << '\t' << d.a12() << '\t' << d.a13() << std::endl
	    << d.a12() << '\t' << d.a22() << '\t' << d.a23() << std::endl
	    << d.a13() << '\t' << d.a23() << '\t' << d.a33() << std::endl 
	    << std::endl;

  Symmetric3D d2 = u.inverse().transBwd(matrix);
  std::cout << std::fixed
	    << d2.a11() << '\t' << d2.a12() << '\t' << d2.a13() << std::endl
	    << d2.a12() << '\t' << d2.a22() << '\t' << d2.a23() << std::endl
	    << d2.a13() << '\t' << d2.a23() << '\t' << d2.a33() << std::endl 
	    << std::endl;

  Vec3D phi(1,2,3);
  Vec3D varphi;
  Vec3D psi(3,4,5);

  Orthogonal3D u2(psi,phi,SS_31);

  psi.normalize();
  phi -= psi*(psi*phi);
  phi.normalize();
  varphi = psi^phi;

  std::cerr << phi.x() << '\t' << phi.y() << '\t' << phi.z() << std::endl
	    << varphi.x() << '\t' << varphi.y() << '\t' << varphi.z() << std::endl
	    << psi.x() << '\t' << psi.y() << '\t' << psi.z() << std::endl
	    << std::endl
	    << u2.ax() << '\t' << u2.ay() << '\t' << u2.az() << std::endl
	    << u2.bx() << '\t' << u2.by() << '\t' << u2.bz() << std::endl
	    << u2.cx() << '\t' << u2.cy() << '\t' << u2.cz() << std::endl;
#endif

  return 0;
}


#endif

#ifdef TEST_MAIN_2

#include<iostream>
#include<iomanip>
#include<VSAMath.hpp>

#include<RandomNumbers.hpp>

int main(int argc, char**argv)
{
  RandomNumbers rng("test.dat");
  MatrixND a(10);
  for(unsigned irow=0;irow<10;irow++)
    for(unsigned icol=0;icol<10;icol++)
      a(irow,icol)=rng.Normal();

  MatrixND ainv;
  a.inverse(ainv);

  MatrixND i(a);
  i *= ainv;

  std::cout << std::scientific << std::setprecision(9);
  //  std::cout << std::fixed << std::setprecision(2);

  for(unsigned irow=0;irow<10;irow++)
    {
      std::cout << a(irow,0);
      for(unsigned icol=1;icol<10;icol++)
	std::cout << ' ' << a(irow,icol);
      std::cout << '\n';
    }
  std::cout << '\n';

  for(unsigned irow=0;irow<10;irow++)
    {
      std::cout << ainv(irow,0);
      for(unsigned icol=1;icol<10;icol++)
	std::cout << ' ' << ainv(irow,icol);
      std::cout << '\n';
    }
  std::cout << '\n';

  for(unsigned irow=0;irow<10;irow++)
    {
      std::cout << i(irow,0);
      for(unsigned icol=1;icol<10;icol++)
	std::cout << ' ' << i(irow,icol);
      std::cout << '\n';
    }
  std::cout << '\n';

}

#endif

