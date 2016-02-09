//-*-mode:c++; mode:font-lock;-*-

/*! \file VSTestTimer.hpp

  Class to make testing of execution times simpler

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       09/01/2005
*/

#ifndef VSTESTTIMER_HPP
#define VSTESTTIMER_HPP

#include <cctype>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <sys/times.h>
#include <stdint.h>

namespace VERITAS
{

  class VSTimerCoreIA32
  {
  public:
    VSTimerCoreIA32(): 
      fCPUCycleTimeNano(), fCPUCycleTimeAtto(), fStartCycle() 
    { 
      std::ifstream stream("/proc/cpuinfo");
      std::string line;
      while(stream)
	{
	  std::getline(stream,line);

	  std::string::size_type icolon = line.find(':');
	  if(icolon==std::string::npos)continue;
	  
	  std::string::size_type ifind;
	  
	  ifind=icolon-1;
	  while((ifind!=0)&&(isspace(line[ifind])))ifind--;
	  std::string key;
	  if(ifind!=std::string::npos)key=line.substr(0,ifind+1);

	  if(key!="cpu MHz")continue;

	  ifind=icolon+1;
	  while((ifind!=line.size())&&(isspace(line[ifind])))ifind++;
	  std::string value;
	  if(ifind!=line.size())value=line.substr(ifind);
	  
	  std::istringstream mhz_stream(value);
	  double mhz;
	  mhz_stream >> mhz;

	  uint64_t atto = uint64_t(round((1/mhz)*1e12));
	  fCPUCycleTimeAtto = atto%1000000000ULL;
	  fCPUCycleTimeNano = atto/1000000000ULL;

	  break;
	}
      /* nothing to see here */ 
    }

    void start() 
    { 
      nanotime_ia32(fStartCycle); 
    }

    void stop(uint64_t& time_ns)
    {
      uint64_t end_cycle; 
      nanotime_ia32(end_cycle);

      uint64_t cycles_hi = (end_cycle-fStartCycle)/1000000000ULL;
      uint64_t cycles_lo = (end_cycle-fStartCycle)%1000000000ULL;

      uint64_t diff = 
	cycles_hi*fCPUCycleTimeNano*1000000000ULL
	+cycles_lo*fCPUCycleTimeNano
	+cycles_hi*fCPUCycleTimeAtto
	+(cycles_lo*fCPUCycleTimeAtto)/1000000000ULL;

      time_ns+=diff;
    }
    
    static void nanotime_ia32(uint64_t & val)
    {
      // This code snippet from /usr/include/asm-x86_64/msr.h
      //unsigned int a,d;
      //__asm__ __volatile__("rdtsc" : "=a" (a), "=d" (d));
      //val = uint64_t(a) | (uint64_t(d)<<32);

      __asm__ __volatile__("rdtsc" : "=A" (val) : );
    }

  private:
    uint64_t fCPUCycleTimeNano;
    uint64_t fCPUCycleTimeAtto;
    uint64_t fStartCycle;    
  };

  class VSTimerCoreUNIX
  {
  public:
    VSTimerCoreUNIX(): fStartTime() { /* nothing to see here */ }

    void start() 
    {
      gettimeofday(&fStartTime,0); 
    }

    void stop(uint64_t& time_ns)
    {
      timeval end_time; 
      gettimeofday(&end_time,0); 
      time_ns+=uint64_t(end_time.tv_sec-fStartTime.tv_sec)*1000000000;
      if(end_time.tv_usec>fStartTime.tv_usec)
	time_ns+=uint64_t(end_time.tv_usec-fStartTime.tv_usec)*1000; 
      else
	time_ns+=uint64_t(fStartTime.tv_usec-end_time.tv_usec)*1000; 
    }

  private:
    struct timeval fStartTime;    
  };    

  template<typename TimerCore>
  class VSTestTimer: protected TimerCore
  {
  public:
    VSTestTimer(const std::string& name = std::string("")):
      TimerCore(),
      fName(name), fEnabled(true), fRunning(false), fTimeNS() { }

    void reset() 
    {
      fEnabled=true;
      fRunning=false;
      fTimeNS=0; 
    }

    void start() 
    {
      if(!fEnabled)return;
      if(fRunning)
	{
	  std::cerr << "VSTestTimer::start(): timer ";
	  if(!fName.empty())std::cerr << '"' << fName << '"' << ' ';
	  std::cerr << "already running!" << std::endl;
	}
      else 
	{
	  TimerCore::start();
	  fRunning=true;
	}
    }

    void stop() 
    {
      if(!fEnabled)return;
      if(!fRunning)
	{
	  std::cerr << "VSTestTimer::stop(): timer ";
	  if(!fName.empty())std::cerr << '"' << fName << '"' << ' ';
	  std::cerr << "not running!" << std::endl;
	}
      else
	{
	  TimerCore::stop(fTimeNS);
	  fRunning=false;
	} 
    }

    void enable() { fEnabled=true; }
    void disable() { fEnabled=false; }
    const std::string& name() const { return fName; }
    uint64_t timeNS() const { return fTimeNS; }
  private:
    std::string    fName;
    bool           fEnabled;
    bool           fRunning;
    uint64_t       fTimeNS;
  };

}

template<typename TimerCore>
inline std::ostream& operator << (std::ostream& stream,
				  const VERITAS::VSTestTimer<TimerCore>& timer)
{
  std::ostringstream tmpstream;
  if(!timer.name().empty())tmpstream << timer.name() << ": ";
  tmpstream <<  timer.timeNS()/1000000000 << '.' 
	    << std::setfill('0') << std::setw(9) << timer.timeNS()%1000000000;
  stream << tmpstream.str();
  return stream;
}

/*

From: http://www.ncsa.uiuc.edu/UserInfo/Resources/Hardware/IA32LinuxCluster/Doc/timing.html#cputime

In order to use the Linux ASM timer on the IA64 platform (titan), you will need to compile the routine using the GNU C compiler. The routine is:

unsigned long long int nanotime_ia64(void)
{
    unsigned long long int val;
    __asm__ __volatile__("mov %0=ar.itc" : "=r"(val) :: "memory");
    return(val);
}
IA32
On the IA32 platform (platinum) you can use either the Intel or GNU compiler. The routine is:

unsigned long long int nanotime_ia32(void)
{
     unsigned long long int val;
    __asm__ __volatile__("rdtsc" : "=A" (val) : );
     return(val);
}

*/

#endif // VSTESTTIMER_HPP
