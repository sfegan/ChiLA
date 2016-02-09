//-*-mode:c++; mode:font-lock;-*-

/*! \file RandomNumbers_TNG.cpp

  The next generation random number generator. Features NR3 generators,
  RanLux and old NR2 generator.

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    $Revision: 1.3 $
  \date       10/31/2007

  $Id: RandomNumbers_TNG.cpp,v 1.3 2008/06/27 21:53:09 matthew Exp $

*/

#include <vector>
#include <climits>
#include <unistd.h>
#include <errno.h>

#include<RandomNumbers_TNG.hpp>

// ============================================================================
// ============================================================================
//
// RandomNumbersBase
//
// ============================================================================
// ============================================================================

#ifndef NOTHREADS
pthread_once_t 
RandomNumbersBase::s_td_once_control = PTHREAD_ONCE_INIT;
pthread_mutex_t 
RandomNumbersBase::s_td_mutex = PTHREAD_MUTEX_INITIALIZER;
std::set<std::string> RandomNumbersBase::s_td_locks;

void RandomNumbersBase::initializeThreadOnce(void)
{
  pthread_mutex_init(&s_td_mutex,NULL);
}
#endif

/** 
 *  Internal function. Lock the RNG state file using a pthreads
 *  construct and an external lock file
 */
bool RandomNumbersBase::lockState(std::string& filename)
{
  if(m_lock_file_fd!=-1)return false;

#ifndef NOTHREADS
  pthread_once(&s_td_once_control, initializeThreadOnce);
#endif

  std::string file = filename;
  std::string::size_type ifind;

  // Extract path component of filename
  std::string path;
  ifind = file.length();
  for(std::string::size_type ichar = 0; ichar<file.length(); ichar++)
    if(file[ichar]=='/')ifind=ichar;
  if(ifind!=file.length())
    {
      path=file.substr(0,ifind+1);
      file=file.substr(ifind+1);
    }
  
  // Extract extension component of filename
  std::string ext;
  ifind = file.length();
  for(std::string::size_type ichar = 0; ichar<file.length(); ichar++)
    if(file[ichar]=='.')ifind=ichar;
  if(ifind!=file.length())
    {
      ext=file.substr(ifind);
      file=file.substr(0,ifind);
    }

  // Loop to find an unlocked file of type /path/filename_N.extension
  // where N=0,1,... are tested until one is found
  unsigned instance = 0;
  bool continue_looping = true;

  while(continue_looping)
    {
      // Construct the lock file name

      char buffer[20];
      sprintf(buffer,"_%u",instance);
      m_lock_file_name = 
	path+std::string(".lock_")+file+std::string(buffer)+ext;

#ifndef NOTHREADS
      // If we are using pthreads then make sure the file is not
      // locked from a different thread, since file locking call will
      // allow multiple locks to be acquired from one process
      pthread_mutex_lock(&s_td_mutex);
      if(s_td_locks.find(m_lock_file_name) != s_td_locks.end())
	{
	  pthread_mutex_unlock(&s_td_mutex);
	  instance++;
	  continue;
	}
      s_td_locks.insert(m_lock_file_name);
      pthread_mutex_unlock(&s_td_mutex);
#endif

      // Create/Open the lock file
      m_lock_file_fd = open(m_lock_file_name.c_str(),
			    O_WRONLY|O_CREAT,S_IRUSR|S_IWUSR);

      if(m_lock_file_fd<0)
	{
	  std::ostringstream stream;
	  stream 
	    << "RandomNumbersBase::lockState: could not open lockfile: " 
	    << m_lock_file_name << '\n'
	    << "RandomNumbersBase::lockState: " << strerror(errno)
	    << '\n';
	  throw stream.str();
	}

      // Try to acquire the lock - use fcntl() rather than flock()
      // since it claims to work over NFS, which might be important on
      // a cluster.
      m_lock.l_type   = F_WRLCK;
      m_lock.l_whence = SEEK_SET;
      m_lock.l_start  = 0;
      m_lock.l_len    = 0;
      m_lock.l_pid    = 0;
      
      if(fcntl(m_lock_file_fd, F_SETLK, &m_lock)<0)
	{
	  if((errno==EACCES)||(errno==EAGAIN))
	    {
#ifndef NOTHREADS
	      // If the lock failed then we must delete the entry
	      // in the threads table
	      pthread_mutex_lock(&s_td_mutex);
	      s_td_locks.erase(m_lock_file_name);
	      pthread_mutex_unlock(&s_td_mutex);
#endif
	      // Close the lock file and try the next iteration
	      close(m_lock_file_fd);
	      instance++;
	    }
	  else
	    {
	      std::ostringstream stream;
	      stream 
		<< "RandomNumbersBase::lockState: could not lock lockfile: " 
		<< m_lock_file_name << '\n'
		<< "RandomNumbersBase::lockState: " << strerror(errno)
		<< '\n';
	      throw stream.str();
	    }
	}
      else
	{

	  filename = path+file+std::string(buffer)+ext;
	  continue_looping=false;
	}
    }
  return true;
}

/** 
 *  Internal function. Unlock the RNG state file using a pthreads
 *  construct and an external lock file
 */
void RandomNumbersBase::unlockState()
{
  if(m_lock_file_fd==-1)return;

  m_lock.l_type   = F_UNLCK;
  m_lock.l_whence = SEEK_SET;
  m_lock.l_start  = 0;
  m_lock.l_len    = 0;
  m_lock.l_pid    = 0;
  fcntl(m_lock_file_fd, F_SETLKW, &m_lock);

  close(m_lock_file_fd);
  m_lock_file_fd=-1;

#ifndef NOTHREADS
  pthread_mutex_lock(&s_td_mutex);
  s_td_locks.erase(m_lock_file_name);
  pthread_mutex_unlock(&s_td_mutex);
#endif  
}

/**
 *  Return the value of ln[Factorial(ix)] using cache and gammln(ix+1).
 */
double RandomNumbersBase::factln(unsigned ix)
{
  if(ix<m_factln.size())
    {
      if(m_factln[ix] < 0)m_factln[ix]=gammln(ix+1.0);
      return m_factln[ix];
    }
  return gammln(ix+1.0);
}

/**
 *  Returns the value ln[Gamma(xx)] for xx > 0. 
 */
double RandomNumbersBase::gammln(double xx)
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146, -86.50532032941677,
                        24.01409824083091, -1.231739572450155,
                    0.1208650973866179e-2, -0.5395239384953e-5};
  int j;

  if (xx <= 0) 
    {
      std::ostringstream stream;
      stream 
	<< "RandomNumbersBase::gammln: argument " << xx << " is not > 0\n"; 
      throw stream.str();
    }

  y=x=xx;
  tmp=x+5.5;
  tmp-=(x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser+=cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

// ============================================================================
// ============================================================================
//
// RNGCore::NR2Ran2
//
// ============================================================================
// ============================================================================

const double RNGCore::NR2Ran2::RNMX=(1.0-FLT_EPSILON);

RNGCore::NR2Ran2::
NR2Ran2(uint64_t seed, const Options& opt):
  ran2_idum1(), ran2_idum2(), ran2_iy(), ran2_iv()
{
  ran2_init(int32_t(seed&0x7FFFFFFF));
}
 
void RNGCore::NR2Ran2::saveCoreState(std::ostream& stream) const
{
  stream << ran2_idum1 << '\n' << ran2_idum2 << '\n' << ran2_iy << '\n';
  for(unsigned j=0;j<RANDOMNUMBERS_NTAB;j++)stream << ran2_iv[j] << '\n';
}

bool RNGCore::NR2Ran2::loadCoreState(std::istream& stream)
{
  int32_t my_ran2_idum1;
  int32_t my_ran2_idum2;
  int32_t my_ran2_iy;
  int32_t my_ran2_iv[RANDOMNUMBERS_NTAB];
  
  unsigned ngood = 0;
  if(stream >> my_ran2_idum1)ngood++;
  if(stream >> my_ran2_idum2)ngood++;
  if(stream >> my_ran2_iy)ngood++;
  for(unsigned j=0;j<RANDOMNUMBERS_NTAB;j++)if(stream >> my_ran2_iv[j])ngood++;
  
  if(ngood!=RANDOMNUMBERS_NTAB+3)return false;

  ran2_idum1 = my_ran2_idum1;
  ran2_idum2 = my_ran2_idum2;
  ran2_iy = my_ran2_iy;
  for(unsigned j=0;j<RANDOMNUMBERS_NTAB;j++)ran2_iv[j] = my_ran2_iv[j];
  return true;
}
     
double RNGCore::NR2Ran2::ran2()
{
  int j;
  int32_t   k;
  double    am=(1.0/RAN2_M1);
  int32_t   imm1=(RAN2_M1-1);
  int32_t   ndiv=(1+(RAN2_M1-1)/RANDOMNUMBERS_NTAB);
  double    tmp;

  k=ran2_idum1/RAN2_Q1;
  ran2_idum1=RAN2_A1*(ran2_idum1-k*RAN2_Q1)-RAN2_R1*k; 
  if (ran2_idum1<0) ran2_idum1+=RAN2_M1;
  k=ran2_idum2/RAN2_Q2;
  ran2_idum2=RAN2_A2*(ran2_idum2-k*RAN2_Q2)-RAN2_R2*k;
  if (ran2_idum2<0) ran2_idum2+=RAN2_M2;
  j=ran2_iy/ndiv;
  ran2_iy=ran2_iv[j]-ran2_idum2;
  ran2_iv[j]=ran2_idum1;
  if(ran2_iy<1) ran2_iy+=imm1;
  tmp=(double)(am*ran2_iy);
  if(tmp>RNMX) return RNMX;
  else return tmp;
}

void RNGCore::NR2Ran2::ran2_init(int32_t idum)
{
  int32_t k;

  if(idum<0)idum=-idum;
  if(idum==0)idum=1;

  ran2_idum2 = idum;
  for(unsigned j=RANDOMNUMBERS_NTAB+7;j;j--) {
    k=idum/RAN2_Q1;
    idum = RAN2_A1*(idum-k*RAN2_Q1)-RAN2_R1*k; 
    if(idum<0)idum+=RAN2_M1;
    if(j<RANDOMNUMBERS_NTAB)ran2_iv[j]= idum;
  }
  ran2_iy=ran2_iv[0];
  ran2_idum1=idum;
  return;
}
    
// ============================================================================
// ============================================================================
//
// RNGCore::RanLuxV32
//
// ============================================================================
// ============================================================================

// Adapted from "ranlxd" with header below

/******************************************************************************
*
* File ranlxd.c
*
* Copyright (C) 2005 Martin Luescher
*
* This software is distributed under the terms of the GNU General Public
* License (GPL)
*
* Random number generator "ranlxd". See the notes 
*
*   "User's guide for ranlxs and ranlxd v3.2" (December 2005)
*
*   "Algorithms used in ranlux v3.0" (May 2001)
*
* for a detailed description
*
* The externally accessible functions are 
*
*   void ranlxd(double r[],int n)
*     Computes the next n double-precision random numbers and 
*     assigns them to the elements r[0],...,r[n-1] of the array r[]
* 
*   void rlxd_init(int level,int seed)
*     Initialization of the generator
*
*   int rlxd_size(void)
*     Returns the number of integers required to save the state of
*     the generator
*
*   void rlxd_get(int state[])
*     Extracts the current state of the generator and stores the 
*     information in the array state[N] where N>=rlxd_size()
*
*   void rlxd_reset(int state[])
*     Resets the generator to the state defined by the array state[N]
*
******************************************************************************/

RNGCore::RanLuxV32::RanLuxV32(uint64_t seed, const Options& opt):
  pr(), prm(), ir(), jr(), is(), is_old(), next(), one_bit(), carry(), x()
{
  int iseed = seed&0x7FFFFFFF;
  if(iseed==0)iseed=1;
  rlxd_init(opt, iseed);
}

void RNGCore::RanLuxV32::saveCoreState(std::ostream& stream) const
{
  std::vector<int> state(rlxd_size());
  rlxd_get(&state.front());
  std::copy(state.begin(), state.end(), 
	    std::ostream_iterator<int>(stream,"\n"));
}

bool RNGCore::RanLuxV32::loadCoreState(std::istream& stream)
{
  std::istream_iterator<int> istream(stream);
  std::vector<int> state;
  state.reserve(rlxd_size());
  for(int istate=0;stream&&istate<rlxd_size();istate++)
    state.push_back(*istream++);
  if(state.size()!=unsigned(rlxd_size()))return false;
  rlxd_reset(&state.front());
  return true;
}

#define RANLUXV32_STEP(pi,pj) \
      d=(*pj).c1.c1-(*pi).c1.c1-carry.c1; \
      (*pi).c2.c1+=(d<0); \
      d+=BASE; \
      (*pi).c1.c1=d&MASK; \
      d=(*pj).c1.c2-(*pi).c1.c2-carry.c2; \
      (*pi).c2.c2+=(d<0); \
      d+=BASE; \
      (*pi).c1.c2=d&MASK; \
      d=(*pj).c1.c3-(*pi).c1.c3-carry.c3; \
      (*pi).c2.c3+=(d<0); \
      d+=BASE; \
      (*pi).c1.c3=d&MASK; \
      d=(*pj).c1.c4-(*pi).c1.c4-carry.c4; \
      (*pi).c2.c4+=(d<0); \
      d+=BASE; \
      (*pi).c1.c4=d&MASK; \
      d=(*pj).c2.c1-(*pi).c2.c1; \
      carry.c1=(d<0); \
      d+=BASE; \
      (*pi).c2.c1=d&MASK; \
      d=(*pj).c2.c2-(*pi).c2.c2; \
      carry.c2=(d<0); \
      d+=BASE; \
      (*pi).c2.c2=d&MASK; \
      d=(*pj).c2.c3-(*pi).c2.c3; \
      carry.c3=(d<0); \
      d+=BASE; \
      (*pi).c2.c3=d&MASK; \
      d=(*pj).c2.c4-(*pi).c2.c4; \
      carry.c4=(d<0); \
      d+=BASE; \
      (*pi).c2.c4=d&MASK


void RNGCore::RanLuxV32::error(int no) const
{
  std::ostringstream stream;
  switch(no)
    {
    case 0:
      stream 
	<< "Error in rlxd_init\n"
	<< "Arithmetic on this machine is not suitable for ranlxd\n";
      break;
    case 1:
      stream 
	<< "Error in subroutine rlxd_init\n"
	<< "Bad choice of luxury level (should be 1 or 2)\n";
      break;
    case 2:
      stream 
	<< "Error in subroutine rlxd_init\n"
	<< "Bad choice of seed (should be between 1 and 2^31-1)\n";
      break;
    case 3:
      stream 
	<< "Error in rlxd_get\n"
	<< "Undefined state (ranlxd is not initialized)\n";
      break;
    case 4:
      stream 
	<< "Error in rlxd_reset\n"
	<< "Arithmetic on this machine is not suitable for ranlxd\n";
      break;
    case 5:
      stream 
	<< "Error in rlxd_reset\n"
	<< "Unexpected input data\n";
      break;
    }
  throw(stream.str());
}
  
void RNGCore::RanLuxV32::update()
{
   int k,kmax,d;
   dble_vec_t *pmin,*pmax,*pi,*pj;

   kmax=pr;
   pmin=&x.vec[0];
   pmax=pmin+12;
   pi=&x.vec[ir];
   pj=&x.vec[jr];
      
   for (k=0;k<kmax;k++) 
   {
      RANLUXV32_STEP(pi,pj);
      pi+=1;
      pj+=1;
      if (pi==pmax)
         pi=pmin;      
      if (pj==pmax)
         pj=pmin; 
   }

   ir+=prm;
   jr+=prm;
   if (ir>=12)
      ir-=12;
   if (jr>=12)
      jr-=12;
   is=8*ir;
   is_old=is;
}


void RNGCore::RanLuxV32::define_constants()
{
   int k;

   one_bit=ldexp(1.0,-24);

   for (k=0;k<96;k++)
   {
      next[k]=(k+1)%96;
      if ((k%4)==3)
         next[k]=(k+5)%96;
   }   
}


void RNGCore::RanLuxV32::rlxd_init(int level,int seed)
{
   int i,k,l;
   int ibit,jbit,xbit[31];
   int ix,iy;

   if ((INT_MAX<2147483647)||(FLT_RADIX!=2)||(FLT_MANT_DIG<24)||
       (DBL_MANT_DIG<48))
      error(0);

   define_constants();
   
   if (level==1)
      pr=202;
   else if (level==2)
      pr=397;
   else
      error(1);
   
   i=seed;

   for (k=0;k<31;k++) 
   {
      xbit[k]=i%2;
      i/=2;
   }

   if ((seed<=0)||(i!=0))
      error(2);

   ibit=0;
   jbit=18;

   for (i=0;i<4;i++)
   {
      for (k=0;k<24;k++)
      {
         ix=0;

         for (l=0;l<24;l++) 
         {
            iy=xbit[ibit];
            ix=2*ix+iy;
         
            xbit[ibit]=(xbit[ibit]+xbit[jbit])%2;
            ibit=(ibit+1)%31;
            jbit=(jbit+1)%31;
         }

         if ((k%4)!=i)
            ix=16777215-ix;

         x.num[4*k+i]=ix;
      }
   }

   carry.c1=0;
   carry.c2=0;
   carry.c3=0;
   carry.c4=0;

   ir=0;
   jr=7;
   is=91;
   is_old=0;
   prm=pr%12;
}


void RNGCore::RanLuxV32::ranlxd(double r[],int n)
{
   int k;

   for (k=0;k<n;k++) 
   {
      is=next[is];
      if (is==is_old)
         update();
      r[k]=one_bit*((double)(x.num[is+4])+one_bit*(double)(x.num[is]));      
   }
}


int RNGCore::RanLuxV32::rlxd_size(void) const
{
   return(105);
}


void RNGCore::RanLuxV32::rlxd_get(int state[]) const
{
   int k;

   state[0]=rlxd_size();

   for (k=0;k<96;k++)
      state[k+1]=x.num[k];

   state[97]=carry.c1;
   state[98]=carry.c2;
   state[99]=carry.c3;
   state[100]=carry.c4;

   state[101]=pr;
   state[102]=ir;
   state[103]=jr;
   state[104]=is;
}


void RNGCore::RanLuxV32::rlxd_reset(int state[])
{
   int k;

   if ((INT_MAX<2147483647)||(FLT_RADIX!=2)||(FLT_MANT_DIG<24)||
       (DBL_MANT_DIG<48))
      error(4);

   define_constants();

   if (state[0]!=rlxd_size())
      error(5);

   for (k=0;k<96;k++)
   {
      if ((state[k+1]<0)||(state[k+1]>=167777216))
         error(5);

      x.num[k]=state[k+1];
   }

   if (((state[97]!=0)&&(state[97]!=1))||
       ((state[98]!=0)&&(state[98]!=1))||
       ((state[99]!=0)&&(state[99]!=1))||
       ((state[100]!=0)&&(state[100]!=1)))
      error(5);
   
   carry.c1=state[97];
   carry.c2=state[98];
   carry.c3=state[99];
   carry.c4=state[100];

   pr=state[101];
   ir=state[102];
   jr=state[103];
   is=state[104];
   is_old=8*ir;
   prm=pr%12;
   
   if (((pr!=202)&&(pr!=397))||
       (ir<0)||(ir>11)||(jr<0)||(jr>11)||(jr!=((ir+7)%12))||
       (is<0)||(is>91))
      error(5);
}
