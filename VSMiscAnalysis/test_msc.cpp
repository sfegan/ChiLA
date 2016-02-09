/*! \file optimize_msc.cpp

  Test a set of cuts, calculated with optimize_msc or otherwise, on a
  bunch of files
  
  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       11/06/2007

  $Id: test_msc.cpp,v 1.10 2008/02/24 17:40:48 sfegan Exp $

*/

#include <iomanip>
#include <fstream>
#include <limits>

#include <VSOptions.hpp>
#include <VSOctaveH5Reader.hpp>
#include <VSFileUtility.hpp>
#include <VSEventData.hpp>
#include <VSSimulationData.hpp>
#include <VSNSpace.hpp>
#include <VSGenSpace.hpp>
#include <VSScaledParameterLibrary.hpp>
#include <VSNSpaceOctaveH5IO.hpp>
#include <Angle.h>

using namespace VERITAS;
using namespace SEphem;

#define A(x) x+=o.x;

void liandma(double Non, double Noff, double alpha,
             double* Nsig, double* Sig5, double* Sig9, double* Sig17)
{
  double alphasq;
  double oneplusalpha;
  double oneplusalphaoveralpha;

  double Ntot;

  alphasq=alpha*alpha;
  oneplusalpha=1.0+alpha;
  oneplusalphaoveralpha=oneplusalpha/alpha;

  *Nsig  = Non - alpha*Noff;
  Ntot   = Non + Noff;

  *Sig5  = *Nsig/sqrt(Non+alphasq*Noff);
  *Sig9  = *Nsig/sqrt(alpha*Ntot);
  *Sig17 = sqrt(2*( Non *log(oneplusalphaoveralpha*(Non/Ntot)) +
                    Noff*log(oneplusalpha*(Noff/Ntot)) ));
  if(*Nsig<0)*Sig17=-*Sig17;

  return;
}

class Counts
{
public:
  // Members
  Counts(): m_on(), m_off(), m_gps_lt(), m_l3_lt() { }
  Counts& operator+=(const Counts& o)
  { A(m_on); A(m_off); A(m_gps_lt); A(m_l3_lt); return *this; }
  void accumulateLT(double gps, double l3) { m_gps_lt+=gps; m_l3_lt+=l3; }
  void accumulateOn() { m_on++; }
  void accumulateOff() { m_off++; }

  void print(std::ostream& stream, double alpha)
  { 
    double background = m_off*alpha;
    double signal;
    double sig5;
    double sig9;
    double sig17;
    liandma(m_on,m_off,alpha,&signal,&sig5,&sig9,&sig17);
    double s2b = signal*signal/background;

    stream << m_on << ' '
	   << m_off << ' '
	   << signal/m_l3_lt*60 << ' '
	   << background/m_l3_lt*60 << ' '
	   << s2b << ' '
	   << 1.0/(s2b/m_l3_lt*60) << ' '
	   << sig17 << ' '
	   << sig17/sqrt(m_l3_lt/3600);
  }

  // Data
  unsigned m_on;
  unsigned m_off;
  double m_gps_lt;
  double m_l3_lt;
};

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname 
         << " [options] stage2_file [stage2_file..]" 
         << std::endl
         << std::endl;
  stream << "Options:" << std::endl;
  options.printUsage(stream);
}

int main(int argc, char** argv)
{
  std::string progname(*argv);
  std::vector<std::string> command_line(argc);
  for(unsigned iarg=0;iarg<unsigned(argc);iarg++)command_line[iarg]=argv[iarg];

  std::cerr << "Command line:";
  for(unsigned iarg=0;iarg<unsigned(argc);iarg++)
    std::cerr << ' ' << argv[iarg];
  std::cerr << std::endl;
 
  VSOptions options(argc, argv, true);

  bool print_usage = false;
  if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;

  bool dump_events = false;
  options.findBoolValue("dump_events", dump_events, true,
			"Dump paremeters of events that pass cuts.");

  unsigned itheta_on = 1;
  options.findWithValue("theta_on", itheta_on,
			"Index into theta array to use for ON source counts.");
  

  std::vector<unsigned> itheta_off;
  itheta_off.push_back(2);
  options.findWithValue("theta_off", itheta_off,
			"Set of indexes into theta array to use for OFF "
			"source counts.");

  unsigned size_cut_nscope=0;
  options.findWithValue("size_cut_nscope", size_cut_nscope, 
			"Number of telescope images that must pass size "
			"cut. Value of zero implies that all images present "
			"must pass.");

  double sizeN_hi = std::numeric_limits<double>::infinity();
  double sizeN_lo = 0;
  
  options.findWithValue("sizeN_hi", sizeN_hi,
			"Set the upper size cut that must be passed by the "
			"N'th brightest telescope. N is set with the "
			"\"size_cut_nscope\" option.");

  options.findWithValue("sizeN_lo", sizeN_lo,
			"Set the lower size cut that must be passed by the "
			"N'th brightest telescope. N is set with the "
			"\"size_cut_nscope\" option.");

  double theta_hi = std::numeric_limits<double>::infinity();
  
  options.findWithValue("theta_hi", theta_hi,
			"Set the theta cut.");

  double width_hi = std::numeric_limits<double>::infinity();
  double width_lo = -2;

  options.findWithValue("width_hi", width_hi,
			"Set the upper mean scaled width cut.");
  options.findWithValue("width_lo", width_lo,
			"Set the lower mean scaled width cut.");

  double length_hi = std::numeric_limits<double>::infinity();
  double length_lo = -2;

  options.findWithValue("length_hi", length_hi,
			"Set the upper mean scaled length cut.");
  options.findWithValue("length_lo", length_lo,
			"Set the lower mean scaled length cut.");

  if(!options.assertNoOptions())
    {
      std::cerr << progname << ": unknown options: ";
      for(int i=1;i<argc;i++)
        if(*(argv[i])=='-') std::cerr << ' ' << argv[i];
      std::cerr << std::endl;
      usage(progname, options, std::cerr);
      exit(EXIT_FAILURE);
    }

  argv++,argc--;

  if(print_usage)
    {
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  int arg_req = 1;
  if(argc < arg_req)
    {
      std::cerr << progname << ": need at least " << arg_req
		<< " arguments, got " << argc << std::endl;
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  double alpha = 1.0;
  if(itheta_off.size())alpha = 1.0/double(itheta_off.size());

  std::ostream* stream = &std::cout;
  if(dump_events)stream = &std::cerr;

  Counts cts_all;

  while(argc)
    {
      std::string stage2_file = *argv;
      argv++,argc--;
      VSFileUtility::expandFilename(stage2_file);

      try
	{
	  VSOctaveH5Reader reader(stage2_file);

	  Counts cts_one;
	  double gps_lt = 0;
	  double l3_lt = 0;
	  reader.readScalar("diagnostics.gps_livetime_sec", gps_lt);
	  reader.readScalar("diagnostics.l3_livetime_sec", l3_lt);
	  
	  cts_one.accumulateLT(gps_lt, l3_lt);

	  VSEventDataReader ed(reader.readStruct("events"));
	  unsigned nscope = ed.nscope();
	  unsigned nrow = ed.rows();
	  
	  std::cerr << stage2_file << ": " << nrow << std::endl;

	  for(unsigned irow=0;irow<nrow;irow++)
	    {
	      VSEventArrayDatum d;
	      ed.element(d,irow);

	      double sizeN;
	      std::vector<double> sizes;
	      for(unsigned iscope=0;iscope<nscope;iscope++)
		if((d.scope[iscope])
		   &&(d.scope[iscope]->used_in_reconstruction))
		  sizes.push_back(d.scope[iscope]->fp_N);
	      std::sort(sizes.begin(),sizes.end());
	      if(sizes.empty())continue;
	      else if(size_cut_nscope==0)sizeN=sizes.front();
	      else if(size_cut_nscope<=sizes.size())
		sizeN=sizes[sizes.size()-size_cut_nscope];
	      else continue;

#if 0
	      std::cout << sizeN << ' ' << d.msc_width << ' ' << d.msc_length
			<< '\n';
	      std::cout << (sizeN < sizeN_lo) << ' '
			<< (sizeN > sizeN_hi) << ' '
			<< (d.msc_width < width_lo) << ' '
			<< (d.msc_width > width_hi) << ' '
			<< (d.msc_length < length_lo) << ' '
			<< (d.msc_length > length_hi) << '\n';
#endif

	      if(sizeN < sizeN_lo || sizeN > sizeN_hi
		 || d.msc_width < width_lo || d.msc_width > width_hi
		 || d.msc_length < length_lo || d.msc_length > length_hi)
		continue;

#if 0
	      std::cout << "PASSED" << std::endl;
#endif

	      assert(itheta_on < d.theta.size());
	      double theta = d.theta[itheta_on];
	      if(theta <= theta_hi)cts_one.accumulateOn();

	      if(dump_events && theta <= theta_hi)
		{
		  std::cout << 0 << ' ' << irow << ' ' << d.event_num << ' ' 
			    << sizeN << ' ' << d.msc_width << ' ' 
			    << d.msc_length << ' ' << d.msc_disp << ' '
			    << theta;
		  for(unsigned iscope=0;iscope<nscope;iscope++)
		if((d.scope[iscope])
		   &&(d.scope[iscope]->used_in_reconstruction))
			  std::cout << ' ' << d.scope[iscope]->fp_N;
	else
			  std::cout << ' ' << 0.0;
		  std::cout << '\n';
		}

	      for(unsigned itheta=0;itheta<itheta_off.size();itheta++)
		{
		  assert(itheta_off[itheta] < d.theta.size());
		  theta = d.theta[itheta_off[itheta]];

		  if(dump_events && theta <= theta_hi)
		    {
		      std::cout << itheta+1 << ' ' << irow << ' ' << d.event_num << ' ' 
				<< sizeN << ' ' << d.msc_width << ' ' 
				<< d.msc_length << ' ' << d.msc_disp << ' '
				<< theta;
		      for(unsigned iscope=0;iscope<nscope;iscope++)
		if((d.scope[iscope])
		   &&(d.scope[iscope]->used_in_reconstruction))
			  std::cout << ' ' << d.scope[iscope]->fp_N;
	else
			  std::cout << ' ' << 0.0;
		      std::cout << '\n';
		    }

		  if(theta <= theta_hi)cts_one.accumulateOff();
		}
	    }

	  (*stream) << stage2_file << ": ";
	  cts_one.print(*stream,alpha);
	  (*stream) << std::endl;

	  cts_all += cts_one;
	}
      catch(const VSOctaveH5Exception& e)
	{
	  std::cerr << "Caught instance of VSOctaveH5Exception" << std::endl
		    << e.message() << std::endl
		    << "Skipping file: " << stage2_file << std::endl;
	}
      catch(const std::exception& x)
	{
	  std::cerr << "Caught instance of " << x.what() << std::endl
		    << "Skipping file: " << stage2_file << std::endl;
	}
      catch(...)
	{
	  std::cerr << "Caught some exception" << std::endl
		    << "Skipping file: " << stage2_file << std::endl;
	}
    }

  (*stream) << "total: ";
  cts_all.print(*stream,alpha);
  (*stream) << std::endl;
}
