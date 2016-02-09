//-*-mode:c++; mode:font-lock;-*-

/*! \file VSEventData.hpp

  Data structures for event analysis

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       04/13/2007

  $Id: VSEventData.hpp,v 3.20 2008/07/29 22:34:22 matthew Exp $

*/

#ifndef VSEVENTDATA_HPP
#define VSEVENTDATA_HPP

#include<vector>

#include<VSTime.hpp>
#include<VSOctaveIO.hpp>

namespace VERITAS
{
  
  struct VSEventScopeDatum
  {
    VSEventScopeDatum():
      has_image(), used_in_reconstruction(), 
      N(), nimage(), ntrig(), R(), d1(), d2(), 
      theta1(), theta2(), 
      delta1(), delta2l(), delta2m(), G(), t01(), t02(), lambdad(), lambdac(),
      fp_trigger_ichan(), fp_N(), fp_xc(), fp_yc(), fp_dist(), 
      fp_length(), fp_width(), fp_psi(), fp_ex(), fp_ey(), fp_disp(), 
      intrinsic_width(), intrinsic_length(),
      sc_width(), sc_length(), sc_disp(),
      lt_log10_energy(0), lt_log10_energy_err(0)
    { /* nothing to see here */ }

    bool         has_image;
    bool         used_in_reconstruction;
    double       N;
    uint16_t     nimage;
    uint16_t     ntrig;
    double       R;
    double       d1;
    double       d2;
    double       theta1;
    double       theta2;
    double       delta1;
    double       delta2l;
    double       delta2m;
    double       G;
    double       t01;
    double       t02;
    double       lambdad;
    double       lambdac;
    unsigned     fp_trigger_ichan;
    double       fp_N;
    double       fp_xc;
    double       fp_yc;
    double       fp_dist;
    double       fp_length;
    double       fp_width;
    double       fp_psi;
    double       fp_ex;
    double       fp_ey;
    double       fp_disp;
    double       intrinsic_width;
    double       intrinsic_length;
    double       sc_width;
    double       sc_length;
    double       sc_disp;
    double       lt_log10_energy;
    double       lt_log10_energy_err;
   
    static void _compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDMEMBER(c,VSEventScopeDatum,has_image);
      H5_ADDMEMBER(c,VSEventScopeDatum,used_in_reconstruction);
      H5_ADDMEMBER(c,VSEventScopeDatum,N);
      H5_ADDMEMBER(c,VSEventScopeDatum,nimage);
      H5_ADDMEMBER(c,VSEventScopeDatum,ntrig);
      H5_ADDMEMBER(c,VSEventScopeDatum,R);
      H5_ADDMEMBER(c,VSEventScopeDatum,d1);
      H5_ADDMEMBER(c,VSEventScopeDatum,d2);
      H5_ADDMEMBER(c,VSEventScopeDatum,theta1);
      H5_ADDMEMBER(c,VSEventScopeDatum,theta2);
      H5_ADDMEMBER(c,VSEventScopeDatum,delta1);
      H5_ADDMEMBER(c,VSEventScopeDatum,delta2l);
      H5_ADDMEMBER(c,VSEventScopeDatum,delta2m);
      H5_ADDMEMBER(c,VSEventScopeDatum,G);
      H5_ADDMEMBER(c,VSEventScopeDatum,t01);
      H5_ADDMEMBER(c,VSEventScopeDatum,t02);
      H5_ADDMEMBER(c,VSEventScopeDatum,lambdad);
      H5_ADDMEMBER(c,VSEventScopeDatum,lambdac);
      H5_ADDMEMBER(c,VSEventScopeDatum,fp_trigger_ichan);
      H5_ADDMEMBER(c,VSEventScopeDatum,fp_N);
      H5_ADDMEMBER(c,VSEventScopeDatum,fp_xc);
      H5_ADDMEMBER(c,VSEventScopeDatum,fp_yc);
      H5_ADDMEMBER(c,VSEventScopeDatum,fp_dist);
      H5_ADDMEMBER(c,VSEventScopeDatum,fp_length);
      H5_ADDMEMBER(c,VSEventScopeDatum,fp_width);
      H5_ADDMEMBER(c,VSEventScopeDatum,fp_psi);
      H5_ADDMEMBER(c,VSEventScopeDatum,fp_ex);
      H5_ADDMEMBER(c,VSEventScopeDatum,fp_ey);
      H5_ADDMEMBER(c,VSEventScopeDatum,fp_disp);
      H5_ADDMEMBER(c,VSEventScopeDatum,intrinsic_width);
      H5_ADDMEMBER(c,VSEventScopeDatum,intrinsic_length);
      H5_ADDMEMBER(c,VSEventScopeDatum,sc_width);
      H5_ADDMEMBER(c,VSEventScopeDatum,sc_length);
      H5_ADDMEMBER(c,VSEventScopeDatum,sc_disp);
      H5_ADDMEMBER(c,VSEventScopeDatum,lt_log10_energy);
      H5_ADDMEMBER(c,VSEventScopeDatum,lt_log10_energy_err);
    }
  };

  struct VSEventArrayDatum
  {
    VSEventArrayDatum(const std::vector<unsigned>& nchan = 
		      std::vector<unsigned>(),
		      unsigned ntarget = 0):
      event_num(), 
      abs_event_time(), mjd(), event_time(), elapsed_ticks(), live_ticks(), 
      trigger_mask(), l3_sent_mask(), has_datum_mask(),
      has_image_mask(), used_in_reconstruction_mask(),
      mean_array_zn(), mean_array_az(), nscope_image(), chi2e(), chi2R(),
      mean_fov_x(), mean_fov_y(), 
      mean_derotated_fov_x(), mean_derotated_fov_y(),
      zn(), az(), ra(), dec(), ra_J2000(), dec_J2000(), 
      R(), Rx(), Ry(), deltael(), deltaew(), 
      deltaRl(), deltaRw(), theta0(), theta1(),
      N2(), msc_width(), msc_length(), msc_disp(), 
      mlt_log10_energy(0), mlt_log10_energy_chi2(0),
      theta(ntarget), scope(nchan.size())
    { 
      for(unsigned iscope=0;iscope<nchan.size();iscope++)
	if(nchan[iscope])scope[iscope]=new VSEventScopeDatum;
    }
    
    ~VSEventArrayDatum()
    {
      for(unsigned iscope=0;iscope<scope.size();iscope++)delete scope[iscope];
    }

    unsigned     event_num;
    VSTime       abs_event_time;
    double       mjd;
    double       event_time;
    uint64_t     elapsed_ticks;
    uint64_t     live_ticks;
    uint32_t     trigger_mask;
    uint32_t     l3_sent_mask;
    uint32_t     has_datum_mask;
    uint32_t     has_image_mask;
    uint32_t     used_in_reconstruction_mask;
    double       mean_array_zn;
    double       mean_array_az;
    unsigned     nscope_image;
    double       chi2e;
    double       chi2R;
    double       mean_fov_x;
    double       mean_fov_y;
    double       mean_derotated_fov_x;
    double       mean_derotated_fov_y;
    double       zn;
    double       az;
    double       ra;
    double       dec;
    double       ra_J2000;
    double       dec_J2000;
    double       R;
    double       Rx;
    double       Ry;
    double       deltael;
    double       deltaew;
    double       deltaRl;
    double       deltaRw;
    double       theta0;
    double       theta1;
    double       N2;
    double       msc_width;
    double       msc_length;
    double       msc_disp;
    double       mlt_log10_energy;
    double       mlt_log10_energy_chi2;

    std::vector<double> theta;
    std::vector<VSEventScopeDatum*> scope;

    static void _compose(VSOctaveH5CompositeDefinition& c)
    {
      H5_ADDMEMBER(c,VSEventArrayDatum,event_num);
      H5_ADDSIMPLECOMPOSITE(c,VSEventArrayDatum,abs_event_time);
      H5_ADDMEMBER(c,VSEventArrayDatum,mjd);
      H5_ADDMEMBER(c,VSEventArrayDatum,event_time);
      H5_ADDMEMBER(c,VSEventArrayDatum,elapsed_ticks);
      H5_ADDMEMBER(c,VSEventArrayDatum,live_ticks);
      H5_ADDMEMBER(c,VSEventArrayDatum,trigger_mask);
      H5_ADDMEMBER(c,VSEventArrayDatum,l3_sent_mask);
      H5_ADDMEMBER(c,VSEventArrayDatum,has_datum_mask);
      H5_ADDMEMBER(c,VSEventArrayDatum,has_image_mask);
      H5_ADDMEMBER(c,VSEventArrayDatum,used_in_reconstruction_mask);
      H5_ADDMEMBER(c,VSEventArrayDatum,mean_array_zn);
      H5_ADDMEMBER(c,VSEventArrayDatum,mean_array_az);
      H5_ADDMEMBER(c,VSEventArrayDatum,nscope_image);
      H5_ADDMEMBER(c,VSEventArrayDatum,chi2e);
      H5_ADDMEMBER(c,VSEventArrayDatum,chi2R);
      H5_ADDMEMBER(c,VSEventArrayDatum,mean_fov_x);
      H5_ADDMEMBER(c,VSEventArrayDatum,mean_fov_y);
      H5_ADDMEMBER(c,VSEventArrayDatum,mean_derotated_fov_x);
      H5_ADDMEMBER(c,VSEventArrayDatum,mean_derotated_fov_y);
      H5_ADDMEMBER(c,VSEventArrayDatum,zn);
      H5_ADDMEMBER(c,VSEventArrayDatum,az);
      H5_ADDMEMBER(c,VSEventArrayDatum,ra);
      H5_ADDMEMBER(c,VSEventArrayDatum,dec);
      H5_ADDMEMBER(c,VSEventArrayDatum,ra_J2000);
      H5_ADDMEMBER(c,VSEventArrayDatum,dec_J2000);
      H5_ADDMEMBER(c,VSEventArrayDatum,R);
      H5_ADDMEMBER(c,VSEventArrayDatum,Rx);
      H5_ADDMEMBER(c,VSEventArrayDatum,Ry);
      H5_ADDMEMBER(c,VSEventArrayDatum,deltael);
      H5_ADDMEMBER(c,VSEventArrayDatum,deltaew);
      H5_ADDMEMBER(c,VSEventArrayDatum,deltaRl);
      H5_ADDMEMBER(c,VSEventArrayDatum,deltaRw);
      H5_ADDMEMBER(c,VSEventArrayDatum,theta0);
      H5_ADDMEMBER(c,VSEventArrayDatum,theta1);
      H5_ADDMEMBER(c,VSEventArrayDatum,N2);
      H5_ADDMEMBER(c,VSEventArrayDatum,msc_width);
      H5_ADDMEMBER(c,VSEventArrayDatum,msc_length);
      H5_ADDMEMBER(c,VSEventArrayDatum,msc_disp);
      H5_ADDMEMBER(c,VSEventArrayDatum,mlt_log10_energy);
      H5_ADDMEMBER(c,VSEventArrayDatum,mlt_log10_energy_chi2);
    }

    VSEventArrayDatum(const VSEventArrayDatum&);
    VSEventArrayDatum& operator=(const VSEventArrayDatum&);
  };

  class VSEventDataWriter
  {
  public:
    VSEventDataWriter(VSOctaveH5WriterStruct* s,
		      std::vector<unsigned> nchan, unsigned ntarget);
    ~VSEventDataWriter();

    bool append(const VSEventArrayDatum& x);

  private:
    VSEventDataWriter(const VSEventDataWriter&);
    VSEventDataWriter& operator=(const VSEventDataWriter&);

    typedef VSOctaveH5WriterCompositeVector<VSEventArrayDatum> ArrayDataWriter;
    typedef VSOctaveH5WriterCompositeVector<VSEventScopeDatum> ScopeDataWriter;

    ArrayDataWriter*                              m_array_writer;
    std::vector<VSOctaveH5WriterVector<double>*>  m_theta_writer;
    std::vector<ScopeDataWriter*>                 m_scope_writer;
  };

  class VSEventDataReader
  {
  public:
    typedef VSOctaveH5ReaderBase::MemberSubset MemberSubset;

    VSEventDataReader(VSOctaveH5ReaderStruct* s,
		      const MemberSubset& array_subset = MemberSubset(),
		      const MemberSubset& scope_subset = MemberSubset());
    ~VSEventDataReader();

    bool element(VSEventArrayDatum& x, unsigned index);

    inline unsigned rows() const { return m_array_reader->rows(); }

    inline unsigned nscope() const { return m_scope_reader.size(); }
    inline bool hasIScope(unsigned iscope) const 
    { return (iscope<nscope())?(m_scope_reader[iscope]!=0):false; }
    
    VSEventArrayDatum at(unsigned index) 
    { 
      VSEventArrayDatum x; 
      if(!element(x,index))throw std::out_of_range(__PRETTY_FUNCTION__); 
      return x; 
    }

    VSEventArrayDatum operator[] (unsigned index) 
    {
      VSEventArrayDatum x; 
      element(x,index); 
      return x;
    }

    static bool loadAllEvents(VSOctaveH5ReaderStruct* s, 
			      std::vector<VSEventArrayDatum>& x);

  private:
    VSEventDataReader(const VSEventDataReader&);
    VSEventDataReader& operator=(const VSEventDataReader&);

    typedef VSOctaveH5ReaderCompositeVector<VSEventArrayDatum> ArrayDataReader;
    typedef VSOctaveH5ReaderCompositeVector<VSEventScopeDatum> ScopeDataReader;

    ArrayDataReader*                              m_array_reader;
    std::vector<VSOctaveH5ReaderVector<double>*>  m_theta_reader;
    std::vector<ScopeDataReader*>                 m_scope_reader;
    bool                                          m_no_has_intrinsic;
  };

}

#endif // defined VSEVENTDATA_HPP
