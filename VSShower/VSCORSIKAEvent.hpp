//-*-mode:c++; mode:font-lock;-*-

/*! \file VSCORSIKAEvent.hpp
  Simple CORSIKA Cerenkov File Dispatcher and Vistor

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       03/17/2005
*/

#ifndef VSCORSIKAEVENT_HPP
#define VSCORSIKAEVENT_HPP

#include <iostream>

#define VSCORSIKAEVENT_MAX_BUNCHES  10000000
#define VSCORSIKAEVENT_MAX_PE       10000000
#define VSCORSIKAEVENT_MAX_TEL      1000
#define VSCORSIKAEVENT_MAX_ARRAY    100

namespace VERITAS
{

  struct VSCORSIKARunParameters
  {
    VSCORSIKARunParameters():
      run_num(), run_date(), run_version(), obs_height(), sample_e_min(), 
      sample_e_max(), sample_slope() { }
    unsigned run_num;          /**< Run number                               */
    unsigned run_date;         /**< Run date                                 */
    float    run_version;      /**< Version of CORSIKA                       */
    float    obs_height;       /**< Height of observation level          [m] */
    float    sample_e_min;     /**< Lower limit of simulated energies  [TeV] */
    float    sample_e_max;     /**< Upper limit of simulated energies  [TeV] */
    float    sample_slope;     /**< Spectral index of power-law spectrum     */
  };

  struct VSCORSIKAExtraRunParameters
  {
    VSCORSIKAExtraRunParameters():
      pri_start_alt(), magnetic_bx(), magnetic_bz(), sample_theta_lo(),
      sample_theta_hi(), sample_phi_lo(), sample_phi_hi(), sample_cone_lo(),
      sample_cone_hi(), cerenk_bunch(), cerenk_event_use(), cerenk_lambda_lo(),
      cerenk_lambda_hi(), cerenk_atmo() { }      
    float    pri_start_alt;    /**< Primary starting altitiude      [g/cm^2] */
    float    magnetic_bx;      /**< Magnetic field strength X-dir       [uT] */
    float    magnetic_bz;      /**< Magnetic field strength Z-dir       [uT] */
    float    sample_theta_lo;  /**< Lower bound on sampled zenith ang  [deg] */
    float    sample_theta_hi;  /**< Upper bound on sampled zenith ang  [deg] */
    float    sample_phi_lo;    /**< Lower bound on sampled azimuth ang [deg] */
    float    sample_phi_hi;    /**< Upper bound on sampled azimuth ang [deg] */
    float    sample_cone_lo;   /**< Lower bound on sampled view cone   [deg] */
    float    sample_cone_hi;   /**< Upper bound on sampled view cone   [deg] */
    unsigned cerenk_bunch;     /**< Cerenkov bunch size                      */
    unsigned cerenk_event_use; /**< Number of times each event is used       */
    float    cerenk_lambda_lo; /**< Cerenkov emission band lo           [nm] */
    float    cerenk_lambda_hi; /**< Cerenkov emission band hi           [nm] */
    bool     cerenk_atmo;      /**< Cerenkov atmospheric absorption on       */
  };

  struct VSCORSIKAArraySpec
  {
    VSCORSIKAArraySpec(): 
      array_max_tel(), array_num_tel() { }
    unsigned array_max_tel;    /**< Maximum number of telescopes             */
    unsigned array_num_tel;    /**< Number of telescopes in run              */
  };
  
  struct VSCORSIKATelescopeSpec
  {
    VSCORSIKATelescopeSpec():
      scope_num(), scope_x(), scope_y(), scope_z(), scope_r() { }
    unsigned scope_num;        /**< Telescope number                         */
    float    scope_x;          /**< X positions of telescope (north)     [m] */
    float    scope_y;          /**< Y positions of telescope (west)      [m] */
    float    scope_z;          /**< Z positions of telescope (up)        [m] */
    float    scope_r;          /**< Radius of sphere enclosing telescope [m] */
  };
  
  struct VSCORSIKAEvent
  {
    VSCORSIKAEvent():
      event_num(), pri_type(), pri_energy(), pri_interact_alt(), pri_px(),
      pri_py(), pri_pz(), pri_azimuth(), pri_elevation() { }
    unsigned event_num;        /**< Event number                             */
    int      pri_type;         /**< Primary particle num                     */
    float    pri_energy;       /**< Primary energy                     [TeV] */
    float    pri_interact_alt; /**< Primary interaction altitude         [m] */
    float    pri_px;           /**< Primary X-momentum               [TeV/c] */
    float    pri_py;           /**< Primary Y-momentum               [TeV/c] */
    float    pri_pz;           /**< Primary Z-momentum               [TeV/c] */
    float    pri_azimuth;      /**< Shower direction azimuth           [deg] */
    float    pri_elevation;    /**< Shower direction altitude          [deg] */
  };

  struct VSCORSIKAEventUse
  {
    VSCORSIKAEventUse():
      use_num(), pri_xcore(), pri_ycore(), pri_impact() { }
    unsigned use_num;          /**< Event (re)use number                     */
    float    pri_xcore;        /**< Core location of event               [m] */
    float    pri_ycore;        /**< Core location of event               [m] */
    float    pri_impact;       /**< Core impact distance                 [m] */
  };

  struct VSCORSIKATelescopeEvent
  {
    VSCORSIKATelescopeEvent():
      scope_num(), rel_xcore(), rel_ycore(), rel_zcore(), rel_impact(), 
      num_bunch(), num_ph() { }
    unsigned scope_num;        /**< Telescope Number                         */
    float    rel_xcore;        /**< Relative core location of event      [m] */
    float    rel_ycore;        /**< Relative core location of event      [m] */
    float    rel_zcore;        /**< Relative core location of event      [m] */
    float    rel_impact;       /**< Relative core impact distance        [m] */
    unsigned num_bunch;        /**< Number of photon bunches                 */
    double   num_ph;           /**< Number of photons                        */
  };

  struct VSCORSIKAPhotonBunch
  {
    VSCORSIKAPhotonBunch():
      ph_count(), ph_rel_x(), ph_rel_y(), ph_rel_impact(), 
      ph_cosine_x(), ph_cosine_y(), ph_cosine_z(), 
      ph_lambda(), ph_time(), ph_height() { }
    float    ph_count;         /**< Number of photons in the bunch           */
    float    ph_rel_x;         /**< X-position of photon at z[scope]     [m] */
    float    ph_rel_y;         /**< Y-position of photon at z[scope]     [m] */
    float    ph_rel_impact;    /**< Impact distance rel. to scope center [m] */
    float    ph_cosine_x;      /**< Direction cosine WRT X-axis              */
    float    ph_cosine_y;      /**< Direction cosine WRT Y-axis              */
    float    ph_cosine_z;      /**< Direction cosine WRT Z-axis              */
    float    ph_lambda;        /**< Wavelength                          [nm] */
    float    ph_time;          /**< Ground impact time from projected 
				    arrival of primary at ground level  [ns] */
    float    ph_height;        /**< Height of production of bunch        [m] */
  };

  class VSCORSIKAEventDispatcherStop
  {
  public:
    virtual ~VSCORSIKAEventDispatcherStop();
    virtual void stopProcessingFile() = 0;
  };

  class VSCORSIKAEventVisitor
  {
  public:
    VSCORSIKAEventVisitor(): fDispatchStop() { }
    virtual ~VSCORSIKAEventVisitor();

    // Register
    virtual void registerDispatcher(VSCORSIKAEventDispatcherStop* dispatch_stop);
    
    // File
    virtual void visitFile(const char* filename);
    virtual void leaveFile();

    // CORSIKA run info
    virtual void visitRun(const VSCORSIKARunParameters& param);
    virtual void visitRunExtra(const VSCORSIKAExtraRunParameters& param);
    virtual void leaveRun();

    // CORSIKA configuarion -- Steering Card Settings and Telescope Positions
    virtual void visitInputConfigEntry(const char* line);
    virtual void visitArraySpec(const VSCORSIKAArraySpec& arrayspec);
    virtual void visitTelescopeSpec(const VSCORSIKATelescopeSpec& scopespec);

    // Event setup
    virtual void visitEvent(const VSCORSIKAEvent& event, bool& veto);
    virtual void leaveEvent(bool veto);
    
    // Event as seen by sampled array
    virtual void visitEventUse(const VSCORSIKAEventUse& use, bool& veto);
    virtual void leaveEventUse(bool veto);

    // Event as seen from each telescope
    virtual void visitTelescopeEvent(const VSCORSIKATelescopeEvent& scope,
				     bool& veto);
    virtual void leaveTelescopeEvent(bool veto);

    // Photon generated in telescope
    virtual void visitPhotonBunch(const VSCORSIKAPhotonBunch& bunch);

  protected:
    VSCORSIKAEventDispatcherStop* fDispatchStop;
  };

  class VSCORSIKAEventDispatcher: private VSCORSIKAEventDispatcherStop
  {
  public:
    VSCORSIKAEventDispatcher(VSCORSIKAEventVisitor* visitor);
    virtual ~VSCORSIKAEventDispatcher();
    unsigned processFile(const char* filename);
    virtual void stopProcessingFile();
    void setDispatchEmptyBunches(bool deb=true) { fDispatchEmptyBunches=deb; }
      private:
    VSCORSIKAEventVisitor*        fVisitor;
    bool                          fStopProcessingFlag;
    bool                          fExtraParametersDispatched;
    bool                          fDispatchEmptyBunches;
  };

  class VSCORSIKAFileDumper: public VSCORSIKAEventVisitor
  {
  public:
    VSCORSIKAFileDumper(std::ostream& s): 
      VSCORSIKAEventVisitor(), fStream(s) { }
    
    virtual void visitFile(const char* filename);
    virtual void visitRun(const VSCORSIKARunParameters& param);
    virtual void visitRunExtra(const VSCORSIKAExtraRunParameters& param);
    virtual void visitInputConfigEntry(const char* line);
    virtual void visitArraySpec(const VSCORSIKAArraySpec& arrayspec);
    virtual void visitTelescopeSpec(const VSCORSIKATelescopeSpec& scopespec);
    virtual void visitEvent(const VSCORSIKAEvent& event, bool& veto);
    virtual void visitEventUse(const VSCORSIKAEventUse& use, bool& veto);
    virtual void visitTelescopeEvent(const VSCORSIKATelescopeEvent& scope, 
				     bool& veto);
    virtual void visitPhotonBunch(const VSCORSIKAPhotonBunch& bunch);
  private:
    std::ostream& fStream;
  };

}

std::ostream& operator <<(std::ostream& s, 
			  const VERITAS::VSCORSIKARunParameters& x);
std::ostream& operator <<(std::ostream& s, 
			  const VERITAS::VSCORSIKAExtraRunParameters& x);
std::ostream& operator <<(std::ostream& s, 
			  const VERITAS::VSCORSIKAArraySpec& x);
std::ostream& operator <<(std::ostream& s, 
			  const VERITAS::VSCORSIKATelescopeSpec& x);
std::ostream& operator <<(std::ostream& s, 
			  const VERITAS::VSCORSIKAEvent& x);
std::ostream& operator <<(std::ostream& s, 
			  const VERITAS::VSCORSIKAEventUse& x);
std::ostream& operator <<(std::ostream& s, 
			  const VERITAS::VSCORSIKATelescopeEvent& x);
std::ostream& operator <<(std::ostream& s, 
			  const VERITAS::VSCORSIKAPhotonBunch& x);

#endif // VSCORSIKAEVENT_HPP
