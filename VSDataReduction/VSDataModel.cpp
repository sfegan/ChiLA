#include <VSDataModel.hpp>

using namespace VERITAS;

void VSOffRegion::getRegions(const VSAAlgebra::Vec2D& src_xy,
			     const VSAAlgebra::Vec2D& obs_xy,
			     double theta, unsigned max_nregion,
			     std::vector< VSAAlgebra::Vec2D >& off_coords)
{
  double rcamera = src_xy.d(obs_xy);
  double dphi = 2*asin(theta/rcamera);
  unsigned nregion = (unsigned)((2*M_PI-2*dphi)/dphi);

  if(max_nregion > 0) nregion = std::min(nregion,max_nregion);
  //  dphi = (2*M_PI-2*dphi)/(double)nregion;
  
  for(unsigned iregion = 0; iregion < nregion; iregion++)
    {
      VSAAlgebra::Vec2D coords = src_xy;      
     
      coords.rotate(obs_xy,M_PI+dphi*((double)(nregion-1)/2.-(double)iregion));
      vsassert(coords.d(src_xy) > 2*theta);
      off_coords.push_back(coords);
    }
}

void VSExclusionRegion::save(VSOctaveH5WriterStruct* writer) const
{
  writer->writeStructCellVector("exclusion_regions",m_regions);
}

bool VSExclusionRegion::load(VSOctaveH5ReaderStruct* reader) 
{
  return reader->readStructCellVector("exclusion_regions",m_regions);
}

VSDataModel::VSDataModel(VSSourceModel* source_fn,
			 VSAcceptanceModel* bkgnd_fn): 
  VSAFunction::CompositeSumFn<VSModelCoord>(),
  m_source_fn(source_fn->clone()), 
  m_bkgnd_fn(bkgnd_fn->clone())
{ 
  setFn1(m_bkgnd_fn);
  setFn2(m_source_fn);
}

VSDataModel::~VSDataModel()
{
  
}

VSDataModel::VSDataModel(const VSDataModel& o)
{
  m_source_fn = o.m_source_fn->clone();
  m_bkgnd_fn = o.m_bkgnd_fn->clone();
  setFn1(m_bkgnd_fn);
  setFn2(m_source_fn);
}

