//-*-mode:c++; mode:font-lock;-*-

/*! \file VSReconstruction.hpp
  Implementation of VVV "in the sky" chi-2 shower core reconstruction
  algorithm

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n

  \version    1.0
  \date       04/16/2005
*/

/*
   Coordinate system:

   As defined in CORSIKA 
   +X - NORTH
   +Y - WEST
   +Z - UP
*/

#ifndef VSRECONSTRUCTION_HPP
#define VSRECONSTRUCTION_HPP

#include <Vec3D.hpp>

#include <vector>
#include <float.h>

#include "VSCORSIKAEvent.hpp"

namespace VERITAS
{

  class VSReconstruction
  {
  public:
    typedef VSCORSIKATelescopeSpec                         ScopeSpec;
    typedef std::vector<ScopeSpec>                         ScopeSpecs;

    struct Axis
    {
      Axis(double _cx, double _cy): cx(_cx), cy(_cy) { }
      double cx;
      double cy;
    };

    struct Ray: public Axis
    {
      Ray(double _cx, double _cy, double _w=1.0): 
	Axis(_cx,_cy), weight(_w) { }
      double weight;
    };

    struct Chi2MapElement
    {
      Chi2MapElement(double _chi2, 
		     double _ec_x, double _ec_y, double _ec_z,
		     double _rc_x, double _rc_y, double _rc_z):
	chi2(_chi2), 
	ec_x(_ec_x), ec_y(_ec_y), ec_z(_ec_z),
	rc_x(_rc_x), rc_y(_rc_y), rc_z(_rc_z) { }
      double chi2;
      double ec_x;
      double ec_y;
      double ec_z;
      double rc_x;
      double rc_y;
      double rc_z;      
    };

    struct MinimizationResults
    {
      MinimizationResults():
	nrays(), weight(), chi2(), 
	ec_x(), ec_y(), ec_z(), rc_x(), rc_y(), rc_z(),
	chi2_curv_l(), chi2_curv_w(), chi2_curv_ang(),
	ps(), ps_mean(), ps_var(), scope_theta_stat(), scope_ps_stat() { }
      MinimizationResults(unsigned _nrays, double _chi2,
			  double _ec_x, double _ec_y, double _ec_z,
			  double _rc_x, double _rc_y, double _rc_z,
			  double _chi2_curv_l, double _chi2_curv_w,
			  double _chi2_curv_ang):
	nrays(_nrays), chi2(_chi2),  
	ec_x(_ec_x), ec_y(_ec_y), ec_z(_ec_z),
	rc_x(_rc_x), rc_y(_rc_y), rc_z(_rc_z),
	chi2_curv_l(_chi2_curv_l), chi2_curv_w(_chi2_curv_w),
	chi2_curv_ang(_chi2_curv_ang), ps(), ps_mean(0.0), ps_var(0.0),
	scope_theta_stat(), scope_ps_stat() { }
      unsigned nrays;
      double weight;
      double chi2;
      double ec_x;
      double ec_y;
      double ec_z;
      double rc_x;
      double rc_y;
      double rc_z;
      double chi2_curv_l;
      double chi2_curv_w;
      double chi2_curv_ang;
      std::vector<std::pair<double, double> > ps;
      double ps_mean;
      double ps_var;

      typedef std::pair<double,double> Stat;

      std::vector<unsigned> scope_id;
      std::vector<double>   scope_ray_weight;
      std::vector<Stat>     scope_theta_stat;
      std::vector<Stat>     scope_ps_stat;
    };

    struct GridSpecs
    {
      GridSpecs(double _radius, double _step): 
	radius(_radius), step_size(_step) { }
      double radius;
      double step_size;
    };
    
    typedef std::vector<Ray>                               SingleScopeRays;
    typedef std::vector<SingleScopeRays>                   AllScopeRays;
    
    VSReconstruction(const ScopeSpecs& scopes, double curvature_distance=0.01);
    virtual ~VSReconstruction();

    // ------------------------------------------------------------------------
    // Calculate MinimizationResults for given axis (direction & core loc)
    // ------------------------------------------------------------------------

    void calculateChi2ForAxis(const AllScopeRays& rays,
			      const MinimizationResults& axis,
			      MinimizationResults& res);

    // ------------------------------------------------------------------------
    // Gridsearch method to find minimum of Chi^2
    // ------------------------------------------------------------------------

    template<class OutputIterator> OutputIterator 
    gridSearch(const AllScopeRays& rays, const Axis& optical_axis,
	       std::vector<GridSpecs> grid, MinimizationResults& res,
	       OutputIterator map_output);
    
    void gridSearch(const AllScopeRays& rays, const Axis& optical_axis,
		    std::vector<GridSpecs> grid, MinimizationResults& res);
    		     
  protected:
    Physics::Vec3D findRotationVector(const Axis& optical_axis) const;
    double distanceToAxis(const Axis& optical_axis, 
			  double ec_x, double ec_y, double ec_z) const;

    ScopeSpecs                  fScopes;
    double                      fCurvatureDistance;

    // ------------------------------------------------------------------------
    // Internal chi^2 calculation routine and data
    // ------------------------------------------------------------------------
    
    void setupArray(const ScopeSpecs& scopes);
    void deleteArray();

    void setupChi2Workspace(const AllScopeRays& rays);
    void deleteChi2Workspace();

    void calculateChi2ForHypothesis(bool minimize_shower_core_location);
    void calculateCurvature(MinimizationResults& res);
    void calculatePProfile(MinimizationResults& res);

    // These variables hold the INPUT and OUTPUT of the calculation
    double ec_x;
    double ec_y;
    double ec_z;
    double rc_x;
    double rc_y;
    double rc_z;
    double chi2;

    // These variables must be setup properly prior to calling calculate_chi2
    unsigned nscopes;
    unsigned nrays;
    double weight;

    // These arrays hold the array information (initialized by setArray)
    unsigned scope_array_size;
    unsigned* scope_nray;
    double** ray_weight;
    double** es_x;
    double** es_y;
    double** es_z;
    unsigned* scope_id;
    double* r_x;
    double* r_y;
    double* r_z;
    double* scope_ray_weight;
    double* scope_weighting;
    double* Q_scope_xx;
    double* Q_scope_yy;
    double* Q_scope_zz;
    double* Q_scope_xy;
    double* Q_scope_xz;
    double* Q_scope_yz;

    // These arrays hold the temporary ray workspace
    unsigned ray_array_size;
    double* es_dot_ec;
    double* one_minus_es_dot_ec_squared;
    double* q_tau_x;
    double* q_tau_y;
    double* q_tau_z;
    double* q_p_x;
    double* q_p_y;
    double* q_p_z;
    double* Q_xx;
    double* Q_yy;
    double* Q_zz;
    double* Q_xy;
    double* Q_xz;
    double* Q_yz;
  };

  template<class OutputIterator> OutputIterator VSReconstruction::
  gridSearch(const AllScopeRays& rays, const Axis& optical_axis,
	     std::vector<GridSpecs> grid, MinimizationResults& res,
	     OutputIterator map_output)
  {
    setupChi2Workspace(rays);
    Physics::Vec3D oa_rot(findRotationVector(optical_axis));

    bool found_one = false;
    MinimizationResults min_res;

    // Center initial search on optical axis
    min_res.ec_x = optical_axis.cx;
    min_res.ec_y = optical_axis.cy;    

    for(std::vector<GridSpecs>::iterator igrid = grid.begin(); 
	igrid != grid.end(); igrid++)
      {
	//found_one = false;
	
	double step = igrid->step_size;
	double r_max = ceil(igrid->radius/step)*step;

	Physics::Vec3D 
	  oa_rot(findRotationVector(Axis(min_res.ec_x,min_res.ec_y)));

	const int grid_start = int(r_max/step);
	const int grid_count = 2*grid_start+1;
	
	for(int x_bin = 0; x_bin<grid_count; x_bin++)
	  {
	    double x = (double(x_bin)-double(grid_start))*step;
	    for(int y_bin = 0; y_bin<grid_count; y_bin++)
	      {
		double y = (double(y_bin)-double(grid_start))*step;
		double r = sqrt(x*x+y*y);
		if(r>r_max)continue;
		
		Physics::Vec3D ec(0,0,-1);
		ec.Rotate(Physics::Vec3D(r*M_PI/180.0,0,0));
		ec.Rotate(Physics::Vec3D(0,0,atan2(y,x)));
		ec.Rotate(oa_rot);
		
		// Transfer to calculateChi2() variables
		ec_x = ec.x;
		ec_y = ec.y;
		ec_z = ec.z;
		
		calculateChi2ForHypothesis(true);
		
		if((!found_one)||(chi2<min_res.chi2))
		  {
		    found_one = true;
		    min_res.chi2 = chi2;
		    min_res.ec_x = ec_x;
		    min_res.ec_y = ec_y;
		    min_res.ec_z = ec_z;
		    min_res.rc_x = rc_x;
		    min_res.rc_y = rc_y;
		    min_res.rc_z = rc_z;
		  }
		
		// Send this point to the output iterator
		*(map_output++)=
		  Chi2MapElement(chi2,ec_x,ec_y,ec_z,rc_x,rc_y,rc_z);
	      }
	  }
      }

    min_res.nrays=nrays;
    min_res.weight=weight;
    res=min_res;
    calculateCurvature(res);
    calculatePProfile(res);
    
    return map_output;
  }

} // end namespace

#endif // VSRECONSTRUCTION_HPP
