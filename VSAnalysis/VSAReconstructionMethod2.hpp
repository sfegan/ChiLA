//-*-mode:c++; mode:font-lock;-*-

/*! \file VSAReconstructionMethod2.hpp

  Implementation of reconstruction 2, 2D simultaneously in focal
  planes of telescope

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       11/11/2005
*/

#ifndef VSARECONSTRUCTIONMETHOD2_HPP
#define VSARECONSTRUCTIONMETHOD2_HPP

#define VSA_ASSERT

#include<vector>

#include<VSAAlgebra.hpp>
#include<VSAReconstructionCommon.hpp>

namespace VERITAS 
{
  namespace VSAReconstructionInternals
  {

    class Method2
    {
    public:
      Method2() { /* nothing to see here */ }
      virtual ~Method2();
      
      double calculateChi2ForHypothesis(VSAAlgebra::Vec3D& e,
					const VSAAlgebra::Vec3D& R,
					VSAAlgebra::Symmetric3D& Y,
					const ArrayMoments& event,
					bool minimize_e = true) const;

      struct MinimizationResults
      {
	MinimizationResults()
	  : chi2(0), e(), R(), at_edge(false), 
	    chi2_e_curv_l(), chi2_e_curv_w(), chi2_e_curv_axis(),
	    chi2_R_curv_l(), chi2_R_curv_w(), chi2_R_curv_axis() { }
	double              chi2;
	VSAAlgebra::Vec3D   e;
	VSAAlgebra::Vec3D   R;
	bool                at_edge;
	double              chi2_e_curv_l;
	double              chi2_e_curv_w;
	VSAAlgebra::Vec3D   chi2_e_curv_axis;
	double              chi2_R_curv_l;
	double              chi2_R_curv_w;
	VSAAlgebra::Vec3D   chi2_R_curv_axis;
      };

      template<class OutputIterator> OutputIterator 
      gridSearch(MinimizationResults& results, const ArrayMoments& event,
		 const std::vector<GridSearchSpecs>& grid,
		 OutputIterator map_output);

      void gridSearch(MinimizationResults& results,
		      const ArrayMoments& event,
		      const std::vector<GridSearchSpecs>& grid);

    private:

      class GridCalculator
      {
      public:
	GridCalculator(const Method2& reconstruction, 
		       const ArrayMoments& event)
	  : fReconstruction(reconstruction), fEvent(event) { }
	  
	double operator () (double x, double y) const;
	double calculate(double x, double y,
			 VSAAlgebra::Vec3D& e, VSAAlgebra::Vec3D& R,
			 VSAAlgebra::Symmetric3D& Y,
			 bool minimize = true) const;

      private:
	const Method2& fReconstruction;
	const ArrayMoments& fEvent;
      };

    };

    template<class OutputIterator> OutputIterator Method2::
    gridSearch(MinimizationResults& results, const ArrayMoments& event,
	       const std::vector<GridSearchSpecs>& grid,
	       OutputIterator map_output)
    {
      GridCalculator calculator(*this,event);
      GridSearchResults search_results;
      runGridSearch(search_results,calculator,grid,map_output);
      GridCurvatureResults curv;
      calculateGridCurvature(curv,search_results,calculator);

      results.chi2             = search_results.min_val;
      results.e                = curv.min_e;
      results.R                = curv.min_R;
      results.at_edge          = search_results.min_at_edge;
      results.chi2_e_curv_l    = curv.curv_Y_l;
      results.chi2_e_curv_w    = curv.curv_Y_w;
      results.chi2_e_curv_axis = curv.curv_Y_axis;
      results.chi2_R_curv_l    = curv.curv_grid_l;
      results.chi2_R_curv_w    = curv.curv_grid_w;
      results.chi2_R_curv_axis.set(cos(curv.curv_grid_theta),
				   sin(curv.curv_grid_theta), 0);
      return map_output;
    }

  } // namespace VSAReconstructionInternals


} // namespace VERITAS

#endif // VSARECONSTRUCTIONMETHOD2_HPP
