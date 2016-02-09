//-*-mode:c++; mode:font-lock;-*-

/*! \file VSASVD.hpp

  Singular Value Decomposition

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       12/02/2007
*/

#ifndef VSASVD_HPP
#define VSASVD_HPP

#include <limits>

#include <VSAAlgebra.hpp>

#define VSA_ASSERT

namespace VERITAS
{
  namespace VSAAlgebra
  {
    //! Performs QR decomposition of an m x n matrix supplied as the
    //! constructor argument.
    //!
    //! A = Q*R
    //!
    //! where Q is an m x m orthogonal matrix and R is an m x n upper
    //! triangular matrix.
    class QRDecomposition
    {
    public:
      
      QRDecomposition(const VSAAlgebra::MatrixND& a);

      void decompose(VSAAlgebra::MatrixND& a,
		     VSAAlgebra::VecND& tau);

      //! Return matrix R of QR decomposition.
      const VSAAlgebra::MatrixND& R() { return m_R; }

      //! Return n x n upper sub-matrix of R.
      const VSAAlgebra::MatrixND& R1() { return m_R1; }

      //! Return matrix Q of QR decomposition.
      const VSAAlgebra::MatrixND& Q() { return m_Q; }
      const VSAAlgebra::VecND& tau() { return m_tau; }

    private:

      double householder(VSAAlgebra::VecND& v);
      void householder_hm(double tau, const VSAAlgebra::VecND& v, 
			  VSAAlgebra::MatrixND& a);
      double dnrm2(const VSAAlgebra::VecND& v);
      double pythag(const double a, const double b);

      VSAAlgebra::MatrixND m_a;
      VSAAlgebra::MatrixND m_R;
      VSAAlgebra::MatrixND m_R1;
      VSAAlgebra::MatrixND m_Q;
      VSAAlgebra::VecND m_tau;
    };

    //! Performs SVD decomposition of an m x n matrix supplied as the
    //! constructor argument.
    //!
    //! A = U*W*V'
    //!
    //! where U is m x m orthogonal matrix, W is an m x n diagonal
    //! matrix of singular values, and V is an n x n orthogonal
    //! matrix.
    class SVD
    {
    public:
      SVD(const VSAAlgebra::MatrixND& a);
    
      const VSAAlgebra::MatrixND& u() const { return m_u; }
      const VSAAlgebra::MatrixND& v() const { return m_v; }
      const VSAAlgebra::VecND& w() const { return m_w; }
      double thresh() const { return m_thresh; }

      //! Find the minimum norm solution of the matrix equation:
      //!
      //! A*x = b
      void solve(const VSAAlgebra::VecND& b, VSAAlgebra::VecND& x, 
		 double thresh = -1.0);
      void solve(const VSAAlgebra::MatrixND& b, VSAAlgebra::MatrixND& x, 
		 double thresh = -1.0);

      //! Return the pseudoinverse of the input matrix A.
      VSAAlgebra::MatrixND inverse() const; 
      VSAAlgebra::MatrixND inverse(double thresh) const;       

      static VSAAlgebra::MatrixND inverse(const VSAAlgebra::MatrixND& a)
      {
	SVD svd(a);
	return svd.inverse();
      }

      //! Return the norm (largest singular value) of the matrix.
      static double norm(const VSAAlgebra::MatrixND& a) 
      {
	SVD svd(a);
	return svd.norm();
      }

      double norm() { return m_w[0]; }
      unsigned rank(double thresh = -1);
      unsigned nullity(double thresh = -1);

      void range(VSAAlgebra::MatrixND& m, double thresh = -1);
      void nullSpace(VSAAlgebra::MatrixND& m, double thresh = -1);

      double invCondition() 
      { return (m_w[0]<=0 || m_w[m_ncol-1]<=0)?0.0:m_w[m_ncol-1]/m_w[0]; }

    private:
      void decompose();
      void reorder();
      double pythag(const double a, const double b);
      void setThresh(double thresh);

      int                    m_nrow;
      int                    m_ncol;
      VSAAlgebra::MatrixND   m_u;
      VSAAlgebra::MatrixND   m_v;
      VSAAlgebra::VecND      m_w;

      double                 m_eps;
      double                 m_thresh;
    };

    class GSVD
    {
    public:
      GSVD(const VSAAlgebra::MatrixND& a,
	   const VSAAlgebra::MatrixND& b);

      const VSAAlgebra::MatrixND& u() const { return m_u; }
      const VSAAlgebra::MatrixND& v() const { return m_v; }
      const VSAAlgebra::MatrixND& x() const { return m_x; }
      const VSAAlgebra::MatrixND& c() const { return m_c; }
      const VSAAlgebra::MatrixND& s() const { return m_s; }

    private:
      
      void decompose();
      void csd(const MatrixND& q1, const MatrixND& q2,
	       MatrixND& u, MatrixND& v, MatrixND& z, MatrixND& c, MatrixND& s);
      void flip(MatrixND& x);

      VSAAlgebra::MatrixND   m_a;
      VSAAlgebra::MatrixND   m_b;

      VSAAlgebra::MatrixND   m_u;
      VSAAlgebra::MatrixND   m_v;
      VSAAlgebra::MatrixND   m_x;
      VSAAlgebra::MatrixND   m_c;
      VSAAlgebra::MatrixND   m_s;

    };

  }
}

#endif // ifndef VSASVD_HPP
