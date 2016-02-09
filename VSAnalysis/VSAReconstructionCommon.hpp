//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAReconstructionCommon.hpp
  Common data structures for all reconstruction methods

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       11/12/2005
*/

#ifndef VSARECONSTRUCTIONCOMMON_HPP
#define VSARECONSTRUCTIONCOMMON_HPP

#define VSA_ASSERT

#include<vector>
#include<cmath>
#include<cfloat>

#include<VSAAlgebra.hpp>

namespace VERITAS 
{
  namespace VSAReconstructionInternals
  {

    struct ScopeInfo
    {
      ScopeInfo(unsigned npixel=0)
	: ri(), pij(npixel), Di() { /* nothing to see here */ }

      //! Position of telescope in global coordinate system
      //! (x=east, y=north, z=up) in meters
      VSAAlgebra::Vec3D ri; 

      //! Focal plane position of pixel in focal plane in radians. With
      //! the telescope at stow a photon with propagation vector
      //! in the global coordinates frame e=(a,-sqrt(1-a^2-b^2),b)
      //! would hit a pixel at pij=(a,b)
      std::vector<VSAAlgebra::Vec2D> pij;

      //! Diameter of the reflector
      double Di;
    };

    struct ArrayInfo
    {
      ArrayInfo(unsigned nscope=0)
	: elevation(), scopes(nscope) { /* nothing to see here */ }
	
      //! Array reference elevation
      double elevation;

      //! Telescope info
      std::vector<ScopeInfo> scopes;
    };

    struct PixelImage
    {
      PixelImage(unsigned _j=0, double _nij=0)
	: j(_j), nij(_nij), tij() { /* nothing to see here */ }

      //! Pixel number (starting with zero)
      unsigned j;
      //! Signal recorded by pixel
      double nij;
      //! Detection time (ns)
      double tij;
    };

    struct ScopeImage
    {
      ScopeImage(unsigned npixel=0)
	: zenith(0), azimuth(0), camera_rotation(0),
	  pixels(npixel), use_in_reconstruction(false)
      { /* nothing to see here */ }

      //! Zenith angle in radians
      double zenith;
      //! Azimuth angle (from north going east) in radians
      double azimuth;
      //! Camera rotation angle in radians. Right handed about the z-axis
      double camera_rotation;
      //! Pixels in the image in this telescope
      std::vector<PixelImage> pixels;
      //! Flag allowing the inclusion of this image in the reconstruction
      bool use_in_reconstruction;
    };

    struct ScopeMoments
    {
      ScopeMoments()
	: Ni(0), Wi(0), Vi(), Ti(), HillasParameters(), Vi_cam(), Ti_cam(), 
	  U_cam_to_event(), ei(), ri(), use_in_reconstruction(false)
      { /* nothing to see here */ }
      
      double                  Ni;
      double                  Wi;
      VSAAlgebra::Vec3D       Vi;
      VSAAlgebra::Symmetric3D Ti;
      VSAAlgebra::Eigen3D     HillasParameters;

      //! Moments in 2-D frame of camera (for convenience)
      VSAAlgebra::Vec2D       Vi_cam;
      VSAAlgebra::Symmetric2D Ti_cam;

      //! Orthogonal transformation from Camera frame to event Frame
      VSAAlgebra::Orthogonal3D U_cam_to_event;

      VSAAlgebra::Vec3D       ei; //! Optical axis in event frame
      VSAAlgebra::Vec3D       ri; //! Radius vector in event frame

      //! Flag allowing the inclusion of this image in the reconstruction
      bool use_in_reconstruction;
    };

    struct ArrayMoments
    {
      ArrayMoments(unsigned nscopes=0)
	: nuse_in_reconstruction(0), N(0), scopes(nscopes), 
	  U_earth_to_event() 
      { /* nothing to see here */ }
      
      unsigned                  nuse_in_reconstruction;
      double                    N;      
      std::vector<ScopeMoments> scopes;
      
      //! Orthogonal transformation from Earth frame to event Frame
      VSAAlgebra::Orthogonal3D U_earth_to_event; 
    };
    
    struct GridSearchSpecs
    {
      GridSearchSpecs(double _radius, double _step): 
	radius(_radius), step_size(_step) { /* nothing to see here */ }

      double                  radius;
      double                  step_size;
    };

    struct GridSearchElement
    {
      GridSearchElement(double _x, double _y, double _val)
	: x(_x), y(_y), val(_val) { /* nothing to see here */ }

      double                  x;
      double                  y;
      double                  val;
    };

    struct GridSearchResults
    {
      GridSearchResults()
	: min_x(0), min_y(0), min_at_edge(false), min_val(0), min_step()
      { /* nothing to see here */ }

      double                  min_x;
      double                  min_y;
      bool                    min_at_edge;
      double                  min_val;
      double                  min_step;
    };

    struct GridCurvatureResults
    {
      GridCurvatureResults():
	min_e(), min_R(),
	curv_grid_l(0), curv_grid_w(0), curv_grid_theta(0),
	curv_Y_l(0), curv_Y_w(0), curv_Y_axis() 
      { /* nothing to see here */ }

      VSAAlgebra::Vec3D       min_e;
      VSAAlgebra::Vec3D       min_R;
      double                  curv_grid_l;
      double                  curv_grid_w;
      double                  curv_grid_theta;
      double                  curv_Y_l;
      double                  curv_Y_w;
      VSAAlgebra::Vec3D       curv_Y_axis;
    };

    class NOOPGridSearchElementOutputIterator
    {
    public:
      NOOPGridSearchElementOutputIterator& operator++() { return *this; }
      NOOPGridSearchElementOutputIterator  operator++(int) { return *this; }
      NOOPGridSearchElementOutputIterator& operator*() { return *this; }
      NOOPGridSearchElementOutputIterator& operator=(const GridSearchElement&) 
      { return *this; }
    };

    template<typename OutputIterator, typename Calculator> 
    OutputIterator runGridSearch(GridSearchResults& results,
				 Calculator& calculator,
				 std::vector<GridSearchSpecs> grid, 
				 OutputIterator map_output)
    {
      bool found_one = false;

      if(grid.size() == 0)return map_output;
      double min_step = grid[0].step_size;

      double grid_center_x = 0;
      double grid_center_y = 0;

      for(std::vector<GridSearchSpecs>::iterator igrid = grid.begin(); 
	  igrid != grid.end(); igrid++)
	{
	  double step = igrid->step_size;
	  double r_max = ceil(igrid->radius/step)*step;
	  double r_max2 = r_max*r_max;
	  if(step<min_step)min_step=step;

	  const int grid_start = int(r_max/step);
	  const int grid_count = 2*grid_start+1;
	
	  for(int x_bin = 0; x_bin<grid_count; x_bin++)
	    {
	      double x = (double(x_bin)-double(grid_start))*step;
	      for(int y_bin = 0; y_bin<grid_count; y_bin++)
		{
		  double y = (double(y_bin)-double(grid_start))*step;
		  double r2 = x*x+y*y;
		  if(r2>r_max2)continue;
		
		  const double val = 
		    calculator(x+grid_center_x, y+grid_center_y);
		
		  if(((!found_one)&&(!std::isnan(val)))||(val<results.min_val))
		    {
		      found_one = true;
		      bool edge=false;
		      if(r2 + step*(step+2*x) > r_max2)edge=true;
		      if(r2 + step*(step-2*x) > r_max2)edge=true;
		      if(r2 + step*(step+2*y) > r_max2)edge=true;
		      if(r2 + step*(step-2*y) > r_max2)edge=true;
		      if(r2 + 2*step*(step+x+y) > r_max2)edge=true;
		      if(r2 + 2*step*(step-x+y) > r_max2)edge=true;
		      if(r2 + 2*step*(step+x-y) > r_max2)edge=true;
		      if(r2 + 2*step*(step-x-y) > r_max2)edge=true;
		      
		      results.min_val      = val;
		      results.min_x        = x+grid_center_x;
		      results.min_y        = y+grid_center_y;
		      results.min_at_edge  = edge;
		    }
		  
		  // Send this point to the output iterator
		  *(map_output++) = GridSearchElement(x,y,val);
		}
	    }
	  
	  grid_center_x = results.min_x;
	  grid_center_y = results.min_y;	  
	}
      results.min_step = min_step;

      return map_output;
    }


    template<typename Calculator> 
    void calculateGridCurvature(GridCurvatureResults& curv,
				const GridSearchResults& results,
				Calculator& calculator)
    {
      const double grid_center_x = results.min_x;
      const double grid_center_y = results.min_y;

      double grid[9];

      VSAAlgebra::Vec3D e;
      VSAAlgebra::Vec3D R;
      VSAAlgebra::Symmetric3D Y;

      grid[4] = 
	calculator.calculate(grid_center_x, grid_center_y, e, R, Y, true);

#if 0
      std::cerr << 4 << ' ' << grid[4] << std::endl;
#endif

      if(grid[4]!=results.min_val)
	{
	  std::cerr << "grid[4]!=results.min_val: " 
		    << grid[4] << "!=" << results.min_val << std::endl;
	}

#ifdef VSA_ASSERT
      assert(grid[4]==results.min_val);
#endif

      curv.min_e = e;
      curv.min_R = R;

      VSAAlgebra::Eigen3D Yeigen;
      Y.eigen(Yeigen);
      curv.curv_Y_l = sqrt(grid[4]/Yeigen.val[1]);
      curv.curv_Y_w = sqrt(grid[4]/Yeigen.val[2]);
      curv.curv_Y_axis = Yeigen.vec[1];
      
      unsigned i=0;
      for(int x_bin = -1; x_bin<=1; x_bin++)
	{
	  double x = double(x_bin)*results.min_step;
	  for(int y_bin = -1; y_bin<=1; y_bin++)
	    {	
	      if((y_bin==0)&&(x_bin==0)){ i++; continue; }
	      double y = double(y_bin)*results.min_step;

	      grid[i] = 
		calculator.calculate(x+grid_center_x, y+grid_center_y,
				     e, R, Y, false);
	      
#if 0
	      std::cerr << i << ' ' << grid[i] << std::endl;
#endif
	      i++;
	    }
	}

      double min_val = results.min_val;
      double norm = 1.0/(min_val*results.min_step*results.min_step);
      double min_val_pp_xy = (grid[0] - grid[2] - grid[6] + grid[8])*norm/4.0;
      double min_val_pp_xx = (grid[3] + grid[5] - 2.0*grid[4])*norm;
      double min_val_pp_yy = (grid[1] + grid[7] - 2.0*grid[4])*norm;

      VSAAlgebra::Symmetric2D 
	min_val_pp(min_val_pp_xx,min_val_pp_yy,min_val_pp_xy);
      VSAAlgebra::Eigen2D min_val_pp_eigen;
      min_val_pp.eigen(min_val_pp_eigen);
      
#if 0
      std::cerr << norm << std::endl;
      std::cerr << "xy " << min_val_pp_xy << std::endl;
      std::cerr << "xx " << min_val_pp_xx << std::endl;
      std::cerr << "yy " << min_val_pp_yy << std::endl;
      std::cerr << "e0 " << min_val_pp_eigen.val[0] << std::endl;
      std::cerr << "e1 " << min_val_pp_eigen.val[1] << std::endl;
#endif

      curv.curv_grid_l = 1/sqrt(min_val_pp_eigen.val[0]);
      curv.curv_grid_w = 1/sqrt(min_val_pp_eigen.val[1]);
      curv.curv_grid_theta = atan2(min_val_pp_eigen.vec[0].y(),
				   min_val_pp_eigen.vec[0].x());
    }

  } //  namespace VSAReconstructionInternals

} //namespace VERITAS 

#endif // defined VSARECONSTRUCTIONCOMMON_HPP
