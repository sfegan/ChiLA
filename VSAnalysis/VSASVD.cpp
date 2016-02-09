//-*-mode:c++; mode:font-lock;-*-

/*! \file VSASVD.cpp

  Singular Value Decomposition

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    0.1
  \date       12/02/2007
*/

#include <stdexcept>
#include <VSASVD.hpp>

using namespace VERITAS;
using namespace VERITAS::VSAAlgebra;

QRDecomposition::QRDecomposition(const VSAAlgebra::MatrixND& a):
  m_a(a), m_R(), m_R1(), m_tau()
{
  decompose(m_a,m_tau);

  const unsigned nrow = a.nrow();
  const unsigned ncol = a.ncol();

  m_R = m_a;

  for(unsigned irow = 0; irow < nrow; irow++)
    for(unsigned icol = 0; icol < ncol; icol++)
      if(irow > icol) m_R(irow,icol) = 0;

  m_Q = VSAAlgebra::MatrixND(nrow,nrow);
  m_R1 = m_R.subMatrix(0,0,ncol,ncol);

  for(unsigned i = 0; i < std::min(nrow,ncol); i++)
    {
      VSAAlgebra::VecND v(nrow);
      //      v[0] = 1;
      for(unsigned j = i; j < nrow; j++)
	{
	  if(i == j) v[j] = 1;
	  else v[j] = m_a(j,i);
	}

      // std::cout << "v = " << std::endl;
      // std::cout << v << std::endl;

      VSAAlgebra::MatrixND Q(nrow,nrow);

      for(unsigned irow = 0; irow < nrow; irow++)
	for(unsigned icol = 0; icol < nrow; icol++)
	  {
	    if(irow == icol) Q(irow,icol) += 1;
	    Q(irow,icol) -= m_tau(i)*v(irow)*v(icol);
	  }

      // std::cout << "Q = " << std::endl;
      // std::cout << Q << std::endl;

      // if(i == 0) m_Q = Q;
      // else m_Q = Q*m_Q;
      if(i == 0) m_Q = Q;
      else m_Q = m_Q*Q;
    }

}

void QRDecomposition::decompose(VSAAlgebra::MatrixND& a,
				VSAAlgebra::VecND& tau)
{
  const unsigned nrow = a.nrow();
  const unsigned ncol = a.ncol();
  const unsigned k = std::min(nrow,ncol);

  tau = VSAAlgebra::VecND(k);

  for(unsigned i = 0; i < k; i++)
    {
      VSAAlgebra::VecND c_full = a.columnVector(i);
      VSAAlgebra::VecND c = c_full.subVector(i,nrow-i);

      double tau_i = householder(c);

      // std::cout << "c = " << std::endl;
      // std::cout << c << std::endl;

      for(unsigned j = i; j < nrow; j++) a(j,i) = c(j-i);

      // std::cout << "a = " << std::endl;
      // std::cout << a << std::endl;

      // std::cout << "tau " << tau_i << " " << tau.ndim() << std::endl;

      tau[i] = tau_i;

      if(i+1 < ncol)
      	{
      	  VSAAlgebra::MatrixND m = a.subMatrix(i,i+1,nrow-i,ncol-(i+1));
	  // std::cout << "m = " << std::endl;
	  // std::cout << m << std::endl;

      	  householder_hm(tau_i,c,m);
      	  a.setSubMatrix(i,i+1,m);
      	}
    }
}

double QRDecomposition::householder(VSAAlgebra::VecND& v)
{
  const unsigned n = v.ndim();

  double alpha, beta, tau ;
      
  VSAAlgebra::VecND x = v.subVector(1,n-1); 
      
  double xnorm = dnrm2(x);
  
  if (xnorm == 0) return 0.0;
  
  alpha = v[0];
  beta = - (alpha >= 0.0 ? +1.0 : -1.0) * pythag(alpha, xnorm);
  tau = (beta - alpha) / beta ;
  
  v *= 1.0 / (alpha - beta);
  v[0] = beta;
  //  for(unsigned idim = 1; idim < n; idim++) v[idim] = x[idim-1];
  
  return tau;
}

void QRDecomposition::householder_hm(double tau, const VSAAlgebra::VecND& v, 
				    VSAAlgebra::MatrixND& a)
{
  const unsigned m = a.nrow();
  const unsigned n = a.ncol();

  for (unsigned j = 0; j < n; j++)
    {
      /* Compute wj = Akj vk */
      
      double wj = a(0,j);
      for (unsigned i = 1; i < m; i++)  /* note, computed for v(0) = 1 above */
	wj += a(i,j)*v(i);
        
        /* Aij = Aij - tau vi wj */
        
        /* i = 0 */
        {
          double A0j = a(0, j);
	  a(0,j) = A0j - tau *  wj;
        }
        
        /* i = 1 .. M-1 */
        
        for (unsigned i = 1; i < m; i++)
          {
            double Aij = a(i,j);
            double vi = v(i);
	    a(i,j) = Aij - tau * vi * wj;
          }
      }
}

double QRDecomposition::dnrm2(const VSAAlgebra::VecND& v)
{
  double scale = 0.0;
  double ssq = 1.0;
  unsigned i;
  unsigned ix = 0;

  const unsigned n = v.ndim();
  if (n == 0) return 0;
  else if (n == 1) return fabs(v[0]);

  for (i = 0; i < n; i++) 
    {
      const double x = v[ix];

      if (x != 0.0) 
	{
	  const double ax = fabs(x);

	  if (scale < ax) {
	    ssq = 1.0 + ssq * (scale / ax) * (scale / ax);
	    scale = ax;
	  } else {
	    ssq += (ax / scale) * (ax / scale);
	  }
	}
      
      ix++;
    }

  return scale * sqrt(ssq);
}

double QRDecomposition::pythag(const double a, const double b) 
{
  double fabsa=fabs(a), fabsb=fabs(b);
  return (fabsa > fabsb ? fabsa*sqrt(1.0+std::pow(fabsb/fabsa,2)) :
	  (fabsb == 0.0 ? 0.0 : fabsb*sqrt(1.0+std::pow(fabsa/fabsb,2))));
}


// ============================================================================

SVD::SVD(const MatrixND& a):
  m_nrow(a.nrow()), m_ncol(a.ncol()), m_u(a), m_v(m_ncol,m_ncol), m_w(m_ncol),
  m_eps(std::numeric_limits<double>::epsilon()), m_thresh()
{
  decompose();
  reorder();
  setThresh(-1.0);
}

void SVD::solve(const VecND& b, VecND& x, double thresh)
{
  if(b.ndim() != unsigned(m_nrow) || x.ndim() != unsigned(m_ncol))throw
    std::invalid_argument(std::string(__PRETTY_FUNCTION__)
			  + ": bad input vector sizes");

  VecND temp(m_ncol); //m_nrow);
  setThresh(thresh);
  
  for(int jcol=0;jcol<m_ncol;jcol++)
    {
      double s = 0.0;
      if(m_w[jcol] > m_thresh)
	{
	  for(int irow=0;irow<m_nrow;irow++)s += m_u(irow,jcol)*b[irow];
	  s /= m_w[jcol];
	}
      temp[jcol] = s;
    }
  
  for(int jcol=0;jcol<m_ncol;jcol++)
    {
      double s = 0.0;
      for(int jjcol = 0;jjcol<m_ncol;jjcol++)
	s += m_v(jcol,jjcol)*temp[jjcol];
      x[jcol] = s;
    }
}

void SVD::solve(const MatrixND& b, MatrixND& x, double thresh)
{
  if(b.nrow() != unsigned(m_ncol) || x.nrow() != unsigned(m_ncol)
     || b.ncol() != x.ncol())
    throw std::invalid_argument(std::string(__PRETTY_FUNCTION__)
				+ ": bad input matrix sizes");
  
  VecND xx(m_ncol);
  for(int jrow=0; jrow<m_nrow; jrow++)
    {
      for(int icol=0; icol<m_ncol; icol++)xx[icol] = b(icol,jrow);
      solve(xx,xx,thresh);
      for(int icol=0; icol<m_ncol; icol++)x(icol,jrow) = xx[icol];
    }
}

VSAAlgebra::MatrixND SVD::inverse() const
{
  return inverse(m_thresh);
}

VSAAlgebra::MatrixND SVD::inverse(double thresh) const
{
  MatrixND v = m_v; 
  MatrixND u = m_u.getTranspose();
  MatrixND w(m_ncol);

  for(int icol = 0; icol < m_ncol;icol++)
    if(fabs(m_w(icol)) > thresh) w(icol,icol) = 1/m_w(icol);
    else w(icol,icol) = 0;

  v*=w;
  v*=u;
  return v;
}

unsigned SVD::rank(double thresh)
{
  setThresh(thresh);
  unsigned n = 0;
  for(int jcol=0;jcol<m_ncol;jcol++)if(m_w[jcol] > m_thresh)n++;
  return n;
}

unsigned SVD::nullity(double thresh)
{
  setThresh(thresh);
  unsigned n = 0;
  for(int jcol=0;jcol<m_ncol;jcol++)if(m_w[jcol] <= m_thresh)n++;
  return n;
}

void SVD::range(VSAAlgebra::MatrixND& m, double thresh)
{
  m.resize(m_nrow,rank(thresh));
  unsigned n=0;
  for(int jcol=0;jcol<m_ncol;jcol++)
    if(m_w[jcol] > m_thresh)
      {
	for(int irow=0;irow<m_nrow;irow++)
	  m(irow,n) = m_u(irow,jcol);
	n++;
      }
}

void SVD::nullSpace(VSAAlgebra::MatrixND& m, double thresh)
{
  m.resize(m_nrow,nullity(thresh));
  unsigned n=0;
  for(int jcol=0;jcol<m_ncol;jcol++)
    if(m_w[jcol] <= m_thresh)
      {
	for(int irow=0;irow<m_nrow;irow++)
	  m(irow,n) = m_u(irow,jcol);
	n++;
      }
}

void SVD::setThresh(double thresh)
{
  if(thresh>=0.0)m_thresh = thresh;
  else m_thresh = 0.5*sqrt(double(m_nrow+m_ncol)+1.0)*m_w[0]*m_eps;
}

// ============================================================================
//
// INTERNAL FUNCTIONS: Straight from NR3 Webnote 2 with some REGEXP changes
//
// ============================================================================

inline static double SIGN(const double& a, const double& b)
{
  return b >= 0.0 ? fabs(a) : -fabs(a);
}

void SVD::decompose() 
{
  bool flag;
  int i,its,j,jj,k,l=0,nm=0;
  double anorm,c,f,g,h,s,scale,x,y,z;
  VecND rv1(m_ncol);
  g = scale = anorm = 0.0;
  for (i=0;i<m_ncol;i++) 
    {
      l=i+2;
      rv1[i]=scale*g;
      g=s=scale=0.0;
      if (i < m_nrow) 
	{
	  for (k=i;k<m_nrow;k++) scale += fabs(m_u(k,i));
	  if (scale != 0.0) 
	    {
	      for (k=i;k<m_nrow;k++) 
		{
		  m_u(k,i) /= scale;
		  s += m_u(k,i)*m_u(k,i);
		}
	      f=m_u(i,i);
	      g = -SIGN(sqrt(s),f);
	      h=f*g-s;
	      m_u(i,i)=f-g;
	      for (j=l-1;j<m_ncol;j++) 
		{
		  for (s=0.0,k=i;k<m_nrow;k++) s += m_u(k,i)*m_u(k,j);
		  f=s/h;
		  for (k=i;k<m_nrow;k++) m_u(k,j) += f*m_u(k,i);
		}
	      for (k=i;k<m_nrow;k++) m_u(k,i) *= scale;
	    }
	}
      m_w[i]=scale *g;
      g=s=scale=0.0;
      if (i+1 <= m_nrow && i+1 != m_ncol) 
	{
	  for (k=l-1;k<m_ncol;k++) scale += fabs(m_u(i,k));
	  if (scale != 0.0) 
	    {
	      for (k=l-1;k<m_ncol;k++) 
		{
		  m_u(i,k) /= scale;
		  s += m_u(i,k)*m_u(i,k);
		}
	      f=m_u(i,l-1);
	      g = -SIGN(sqrt(s),f);
	      h=f*g-s;
	      m_u(i,l-1)=f-g;
	      for (k=l-1;k<m_ncol;k++) rv1[k]=m_u(i,k)/h;
	      for (j=l-1;j<m_nrow;j++) 
		{
		  for (s=0.0,k=l-1;k<m_ncol;k++) s += m_u(j,k)*m_u(i,k);
		  for (k=l-1;k<m_ncol;k++) m_u(j,k) += s*rv1[k];
		}
	      for (k=l-1;k<m_ncol;k++) m_u(i,k) *= scale;
	    }
	}
      anorm=std::max(anorm,(fabs(m_w[i])+fabs(rv1[i])));
    }

  for (i=m_ncol-1;i>=0;i--)
    { 
      if (i < m_ncol-1)
	{
	  if (g != 0.0) 
	    {
	      for (j=l;j<m_ncol;j++)
		m_v(j,i)=(m_u(i,j)/m_u(i,l))/g;
	      for (j=l;j<m_ncol;j++) 
		{
		  for (s=0.0,k=l;k<m_ncol;k++) s += m_u(i,k)*m_v(k,j);
		  for (k=l;k<m_ncol;k++) m_v(k,j) += s*m_v(k,i);
		}
	    }
	  for (j=l;j<m_ncol;j++) m_v(i,j)=m_v(j,i)=0.0;
	}
      m_v(i,i)=1.0;
      g=rv1[i];
      l=i;
    }

  for (i=std::min(m_nrow,m_ncol)-1;i>=0;i--) 
    {
      l=i+1;
      g=m_w[i];
      for (j=l;j<m_ncol;j++) m_u(i,j)=0.0;
      if (g != 0.0) 
	{
	  g=1.0/g;
	  for (j=l;j<m_ncol;j++) 
	    {
	      for (s=0.0,k=l;k<m_nrow;k++) s += m_u(k,i)*m_u(k,j);
	      f=(s/m_u(i,i))*g;
	      for (k=i;k<m_nrow;k++) m_u(k,j) += f*m_u(k,i);
	    }
	  for (j=i;j<m_nrow;j++) m_u(j,i) *= g;
	} 
      else for (j=i;j<m_nrow;j++) m_u(j,i)=0.0;
      ++m_u(i,i);
    }

  for (k=m_ncol-1;k>=0;k--) 
    {
      for (its=0;its<30;its++)
	{
	  flag=true;
	  for (l=k;l>=0;l--) 
	    {
	      nm=l-1;
	      if (l == 0 || fabs(rv1[l]) <= m_eps*anorm)
		{
		  flag=false;
		  break;
		}
	      if (fabs(m_w[nm]) <= m_eps*anorm) break;
	    }
	  if (flag) 
	    {
	      c=0.0;
	      s=1.0;
	      for (i=l;i<k+1;i++) 
		{
		  f=s*rv1[i];
		  rv1[i]=c*rv1[i];
		  if (fabs(f) <= m_eps*anorm) break;
		  g=m_w[i];
		  h=pythag(f,g);
		  m_w[i]=h;
		  h=1.0/h;
		  c=g*h;
		  s = -f*h;
		  for (j=0;j<m_nrow;j++) 
		    {
		      y=m_u(j,nm);
		      z=m_u(j,i);
		      m_u(j,nm)=y*c+z*s;
		      m_u(j,i)=z*c-y*s;
		    }
		}
	    }
	  z=m_w[k];
	  if (l == k) 
	    { 
	      if (z < 0.0) 
		{
		  m_w[k] = -z;
		  for (j=0;j<m_ncol;j++) m_v(j,k) = -m_v(j,k);
		}
	      break;
	    }
	  if (its == 29) 
	    throw std::overflow_error(std::string(__PRETTY_FUNCTION__)
				  + ": no convergence after 30 iterations");
	  x=m_w[l];
	  nm=k-1;
	  y=m_w[nm];
	  g=rv1[nm];
	  h=rv1[k];
	  f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
	  g=pythag(f,1.0);
	  f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
	  c=s=1.0; 
	  for (j=l;j<=nm;j++) 
	    {
	      i=j+1;
	      g=rv1[i];
	      y=m_w[i];
	      h=s*g;
	      g=c*g;
	      z=pythag(f,h);
	      rv1[j]=z;
	      c=f/z;
	      s=h/z;
	      f=x*c+g*s;
	      g=g*c-x*s;
	      h=y*s;
	      y *= c;
	      for (jj=0;jj<m_ncol;jj++) 
		{
		  x=m_v(jj,j);
		  z=m_v(jj,i);
		  m_v(jj,j)=x*c+z*s;
		  m_v(jj,i)=z*c-x*s;
		}
	      z=pythag(f,h);
	      m_w[j]=z;
	      if (z) 
		{
		  z=1.0/z;
		  c=f*z;
		  s=h*z;
		}
	      f=c*g+s*y;
	      x=c*y-s*g;
	      for (jj=0;jj<m_nrow;jj++) 
		{
		  y=m_u(jj,j);
		  z=m_u(jj,i);
		  m_u(jj,j)=y*c+z*s;
		  m_u(jj,i)=z*c-y*s;
		}
	    }
	  rv1[l]=0.0;
	  rv1[k]=f;
	  m_w[k]=x;
	}
    }
}

void SVD::reorder() 
{
  int i,j,k,s,inc=1;
  double sw;
  VecND su(m_nrow);
  VecND sv(m_ncol);
  do { inc *= 3; inc++; } while (inc <= m_ncol);
  do {
    inc /= 3;
    for (i=inc;i<m_ncol;i++) {
      sw = m_w[i];
      for (k=0;k<m_nrow;k++) su[k] = m_u(k,i);
      for (k=0;k<m_ncol;k++) sv[k] = m_v(k,i);
      j = i;
      while (m_w[j-inc] < sw) {
	m_w[j] = m_w[j-inc];
	for (k=0;k<m_nrow;k++) m_u(k,j) = m_u(k,j-inc);
	for (k=0;k<m_ncol;k++) m_v(k,j) = m_v(k,j-inc);
	j -= inc;
	if (j < inc) break;
      }
      m_w[j] = sw;
      for (k=0;k<m_nrow;k++) m_u(k,j) = su[k];
      for (k=0;k<m_ncol;k++) m_v(k,j) = sv[k];
    }
  } while (inc > 1);
  for (k=0;k<m_ncol;k++) {
    s=0;
    for (i=0;i<m_nrow;i++) if (m_u(i,k) < 0.) s++;
    for (j=0;j<m_ncol;j++) if (m_v(j,k) < 0.) s++;
    if (s > (m_nrow+m_ncol)/2) {
      for (i=0;i<m_nrow;i++) m_u(i,k) = -m_u(i,k);
      for (j=0;j<m_ncol;j++) m_v(j,k) = -m_v(j,k);
    }
  }
}

double SVD::pythag(const double a, const double b) 
{
  double fabsa=fabs(a), fabsb=fabs(b);
  return (fabsa > fabsb ? fabsa*sqrt(1.0+std::pow(fabsb/fabsa,2)) :
	  (fabsb == 0.0 ? 0.0 : fabsb*sqrt(1.0+std::pow(fabsa/fabsb,2))));
}

// ============================================================================

GSVD::GSVD(const MatrixND& a,const MatrixND& b):
  m_a(a), m_b(b)
{
  vsassert(a.ncol() == b.ncol());

  decompose();

}

// returns unitary matrices U and V,
//   a (usually) square matrix X, and nonnegative diagonal matrices
//   C and S so that
//
//       A = U*C*X'
//       B = V*S*X'
//       C'*C + S'*S = I 

void GSVD::decompose()
{
  const unsigned ncol = m_a.ncol();
  const unsigned nrow1 = m_a.nrow();
  const unsigned nrow2 = m_b.nrow();

  MatrixND ab = m_a;
  ab.concatenate(m_b);

  QRDecomposition qrd(ab);
  
  std::cout << ab << std::endl;

  std::cout << qrd.Q() << std::endl;


  MatrixND q1 = qrd.Q().subMatrix(0,0,nrow1,ncol);
  MatrixND q2 = qrd.Q().subMatrix(nrow1,0,nrow2,ncol);

  std::cout << "Q1: " << std::endl;
  std::cout << q1 << std::endl;
  std::cout << "Q2: " << std::endl;
  std::cout << q2 << std::endl;

  MatrixND z;
  csd(q1,q2,m_u,m_v,z,m_c,m_s);

  std::cout << "U = " << std::endl << m_u << std::endl;
  std::cout << "V = " << std::endl << m_v << std::endl;
  std::cout << "Z = " << std::endl << z << std::endl;
  std::cout << "C = " << std::endl << m_c << std::endl;
  std::cout << "S = " << std::endl << m_s << std::endl;

  


  m_x = qrd.R1().getTranspose()*z;
  
  std::cout << "X = " << std::endl << m_x << std::endl;

}

// ---------------------------------------------------------------------------
// Given Q1 and Q2 such that Q1'*Q1 + Q2'*Q2 = I, the
// C-S Decomposition is a joint factorization of the form
//    Q1 = U*C*Z' and Q2=V*S*Z'
// where U, V, and Z are orthogonal matrices and C and S
// are diagonal matrices (not necessarily square) satisfying
//    C'*C + S'*S = I
// ---------------------------------------------------------------------------
void GSVD::csd(const MatrixND& q1, const MatrixND& q2,
	       MatrixND& u, MatrixND& v, MatrixND& z, MatrixND& c, MatrixND& s)
{
  const unsigned m = q1.nrow();
  const unsigned n = q2.nrow();
  const unsigned p = q1.ncol();

  vsassert(q1.ncol() == q1.ncol());

  SVD svd(q1);

  const unsigned q = std::min(m,p);

  u = svd.u();
  c = MatrixND::makeDiagonal(svd.w());
  z = svd.v();
  z.transpose();

  std::cout << u*c*z << std::endl;
  std::cout << "u = " << std::endl << u << std::endl;
  std::cout << "c = " << std::endl << c << std::endl;
  std::cout << "z = " << std::endl << z << std::endl;

  VSAAlgebra::MatrixND c2 = c;
  for(unsigned i = 0; i < z.nrow(); i++)
    for(unsigned j = 0; j < z.ncol(); j++)
      c2(i,j) = c(c.nrow()-i-1,c.ncol()-j-1);

  c = c2;
  flip(u);
  flip(z);

  std::cout << "u = " << std::endl << u << std::endl;
  std::cout << "c = " << std::endl << c << std::endl;
  std::cout << "z = " << std::endl << z << std::endl;

  s = q2*z;

  std::cout << "s = " << std::endl << s << std::endl;

  unsigned k = 0;

  if(q == 1) k = 0;
  else if(m < p) k = n;
  else
    {
      
      while(k < c.nrow())
	{
	  std::cout << k << " " << c(k,k) << std::endl;
	  if(c(k,k) < 1/sqrt(2.)) break;
	  else k++;
	}
      k++;
    }

  //  VSAAlgebra::MatrixND s2 = s.subMatrix(0,0,s.nrow(),k);

  std::cout << "k = " << k << std::endl;

  QRDecomposition qrd(s.subMatrix(0,0,s.nrow(),k));
  
  v = qrd.Q();
  s = v.getTranspose()*s;

  std::cout << "s = " << std::endl << s << std::endl;
  
  unsigned r = std::min(k,m);

  std::cout << "r = " << r << std::endl;

  for(unsigned i = 0; i < z.nrow(); i++)
    for(unsigned j = 0; j < r; j++)
      if(i != j) s(i,j) = 0;

  std::cout << "s = " << std::endl << s << std::endl;

  if(m == 1 && p > 1) s(0,0) = 0;

  std::cout << "s = " << std::endl << s << std::endl;

  if(k < std::min(n,p))
    {
      r = std::min(n,p);
      
      std::cout << "here" << std::endl;
    }

  if(n < p)
    {
      for(unsigned i = 0; i < s.nrow(); i++)
	for(unsigned j = n+1; j < p; j++)
	  s(i,j) = 0;
    }

  std::cout << "u = " << std::endl << u << std::endl;
  std::cout << "c = " << std::endl << c << std::endl;
  std::cout << "v = " << std::endl << v << std::endl;
  std::cout << "s = " << std::endl << s << std::endl;

  std::cout << "ucz' = " << std::endl << u*c*z.getTranspose() << std::endl;
  std::cout << "vsz' = " << std::endl << v*s*z.getTranspose() << std::endl;
  std::cout << "c'c = " << std::endl << c.getTranspose()*c << std::endl;
  std::cout << "s's = " << std::endl << s.getTranspose()*s << std::endl;

}

void GSVD::flip(MatrixND& x)
{
  VSAAlgebra::MatrixND x2 = x;

  for(unsigned i = 0; i < x.nrow(); i++)
    for(unsigned j = 0; j < x.ncol(); j++)
      x2(i,j) = x(i,x.ncol()-j-1);

  x = x2;
}

// double diagk(const MatrixND& x, unsigned k)
// {

// }

// void diagp()
// {

// }

#ifdef TEST_MAIN

// g++ -O3 -g -Wall -D__STDC_LIMIT_MACROS -D__STDC_CONSTANT_MACROS -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE -DUSEALLOCA -DNOHDF -fPIC -I/usr/include/mysql -I../VSUtility -I../VSSimDB -I../Physics -I../VSShower -I../VSOptics -I../VSAnalysis -I../VSNSpace -I../VSDataReduction -I../SEphem -I. -I/opt/local/include -L/opt/local/lib -I/Users/mdwood/veritas/include -DTEST_MAIN -o test VSASVD.cpp VSAMath.o VSAAlgebra.o -lpthread 

#include<iostream>
#include<iomanip>
#include<VSAMath.hpp>
#include<RandomNumbers.hpp>

int main(int argc, char**argv)
{

}

#endif
