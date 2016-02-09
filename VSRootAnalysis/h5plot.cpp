#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>

// ----------------------------------------------------------------------------
// ChiLA Includes
// ----------------------------------------------------------------------------
#include <VSOptions.hpp>
#include <VSOctaveIO.hpp>
#include <VSOctaveH5Reader.hpp>
#include <VSLineTokenizer.hpp>
#include <Vec3D.hpp>
#include <VSTimeDepPedData.hpp>
#include <VSRunInfoData.hpp>
#include <VSLaserData.hpp>
#include <VSEventData.hpp>
#include <VSNSpace.hpp>
#include <VSNSpaceOctaveH5IO.hpp>
#include <VSPointing.hpp>

// ----------------------------------------------------------------------------
// ROOT Includes
// ----------------------------------------------------------------------------
#include <TROOT.h>
#include <TH1F.h>
#include <TRint.h>
#include <TPDF.h>
#include <TStyle.h>
#include <TH2F.h>
#include <TF1.h>
#include <THStack.h>
#include <TPaveText.h>

// ----------------------------------------------------------------------------
// Local Includes
// ----------------------------------------------------------------------------
#include "VSRHistogram.hpp"
#include "VSRHistogram1D.hpp"
#include "VSRHistogram2D.hpp"
#include "VSRH5Loader.hpp"
#include "VSRCanvas.hpp"

using namespace VERITAS;
using namespace Physics;

template <class T>
inline std::string to_string (const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}

struct Options
{
  Options():
    normalize(),
    cumulative(),
    residual(false),
    ratio(false),
    log(false),
    rebin(1),
    projectx(false),
    projecty(false)
 { }

  std::vector< std::pair< double, double > > point;
  std::vector< std::pair< unsigned, double > > nspace_slice;
  std::vector< triple< unsigned, unsigned, unsigned > > nspace_marginalize;
  std::vector< unsigned > nspace_project;
  triple< unsigned, double,double > hist_range;
  std::string normalize;
  std::string cumulative;
  bool residual;
  bool ratio;
  bool log;
  std::vector<double> rescale;
  unsigned rebin;

  bool projectx;
  std::vector< std::pair< unsigned, unsigned > > projectx_range;
  bool projecty;
  std::vector< std::pair< unsigned, unsigned > > projecty_range;
};

void load_hdf5(const std::string& file,  std::string path,
	       Options& opt,
	       std::vector< VSRHistogram1D* >& hists_1d,
	       std::vector< VSRHistogram2D* >& hists_2d,
	       std::vector< VSRGraph1D* >& graphs,
	       std::vector< VSRCanvasElement* >& elements);

void load_txt(const std::string& file,  
	      Options& opt,
	      std::vector< VSRHistogram1D* >& hists_1d,
	      std::vector< VSRHistogram2D* >& hists_2d,
	      std::vector< VSRGraph1D* >& graphs,
	      std::vector< VSRCanvasElement* >& elements);

void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname << " [options] <h5 file> <path>" << std::endl
         << std::endl;
  stream << "Description: " 
	 << "Plot a histogram, graph, or histogram in an H5 file." 
	 << std::endl;
  stream << "Options:" << std::endl;
  options.printUsage(stream);
}



int main(int argc, char **argv) 
{
  // --------------------------------------------------------------------------
  // Parse options
  // --------------------------------------------------------------------------
  std::string progname(*argv);
  VSOptions options(argc, argv, true);

  bool print_usage = false;
  if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;

  VSRCanvas::configure(options);



  std::vector< std::string > labels;
  std::vector< std::string > bin_labels;
  std::string title;

  std::vector< std::string > h5_paths;

  std::pair< double, double > limits_y;
  std::pair< double, double > limits_x;
  std::pair< double, double > limits_z;


  Options opt;




  bool dump = false;

  options.findWithValue("normalize", opt.normalize, 
                        "");

  options.findWithValue("point", opt.point,"");

  options.findWithValue("nspace_slice", opt.nspace_slice, 
                        "Extract a slice of an nspace along one or more "
			"dimensions.  The option should "
			"specify a comma-separated list which gives the index "
			"of each dimension and a coordinate along that "
			"dimension (e.g. 0/0.5,1/2.5).");

  options.findWithValue("nspace_marginalize", opt.nspace_marginalize, 
			"");

  options.findWithValue("nspace_project", opt.nspace_project, 
			"When plotting an nspace with 2 or more dimensions "
			"dimensions, specify the dimensions along which "
			"you want to project the nspace.  The number of "
			"projecting dimensions cannot exceed two.  Note "
			"that projection is performed after selection "
			"of nspace slices with -nspace_slice.");

  options.findWithValue("rebin", opt.rebin, 
                        "Rebinning factor used for merging histogram bins.");

  options.findWithValue("cumulative", opt.cumulative, 
                        "When plotting histograms, show the cumulative "
			"distribution.");

  options.findBoolValue("residual", opt.residual, true,
                        "When plotting histograms, show the residual "
			"between the first histogram and each subsequent "
			"histogram.");

  options.findBoolValue("ratio", opt.ratio, true,
                        "When plotting 1D histograms or graphs show the ratio "
			"of each with respect to the first.");

  options.findBoolValue("log", opt.log, true,
                        "");

  options.findBoolValue("dump", dump, true,
                        "");
  
  options.findWithValue("rescale", opt.rescale,
                        "");

  options.findWithValue("ylim", limits_y,
                        "Specify range for Y axis.");

  options.findWithValue("xlim", limits_x,
                        "Specify range for X axis.");

  options.findWithValue("zlim", limits_z,
                        "Specify range for X axis.");
  
  options.findWithValue("path", h5_paths,
                        "Path within the H5 for the object to be plotted.");

  options.findWithValue("labels", labels,
                        "");

  options.findWithValue("hist_range", opt.hist_range,
                        "");


  options.findWithValue("bin_labels", bin_labels,
                        "");

  options.findWithValue("title", title,
                        "Set the canvas title.");

  if(options.findWithValue("projectx", opt.projectx_range, "")!=
     VSOptions::FS_NOT_FOUND) 
    opt.projectx = true;

  if(options.findWithValue("projecty", opt.projecty_range, "")!=
     VSOptions::FS_NOT_FOUND) 
    opt.projecty = true;

  if(!options.assertNoOptions())
    {
      std::cerr << progname << ": unknown options: ";
      for(int i=1;i<argc;i++)
	if(*(argv[i])=='-') std::cerr << ' ' << argv[i];
      std::cerr << std::endl;
      usage(progname, options, std::cerr);
      exit(EXIT_FAILURE);
    }  

  if(print_usage)
    {
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }

  argv++,argc--;


  std::vector< std::pair< std::string, std::string > > file_paths;

  if(h5_paths.size() == 0 && argc == 2)
    {
      std::vector<std::string> paths;
      VSDataConverter::fromString(paths,argv[1]);

      for(std::vector<std::string>::iterator itr = paths.begin();
	  itr != paths.end(); ++itr)
	file_paths.push_back(std::make_pair(argv[0],*itr));
    }

  for(int iarg = 0; iarg < argc; iarg++)
    {
      if(iarg < (int)h5_paths.size())
	file_paths.push_back(std::make_pair(argv[iarg],h5_paths[iarg]));
      else if(h5_paths.size() == 1)
	file_paths.push_back(std::make_pair(argv[iarg],h5_paths.front()));
    }
  
  // --------------------------------------------------------------------------
  // Main program
  // --------------------------------------------------------------------------
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);

  TApplication *app = new TApplication("app",NULL,NULL,false);




  std::vector< VSRHistogram1D* > hists_1d;
  std::vector< VSRHistogram2D* > hists_2d;
  std::vector< VSRGraph1D* > graphs;
  std::vector< VSRCanvasElement* > elements;

  if(opt.point.size() > 0)
    {
      VSRGraph1D* g = new VSRGraph1D;
      for(std::vector<std::pair<double,double> >::iterator itr =
	    opt.point.begin(); itr != opt.point.end(); ++itr)
	{
	  std::cout << itr->first << " " << itr->second << std::endl;
	  g->set(itr->first, itr->second);
	}

      graphs.push_back(g);
      elements.push_back(g);
    }



  for(std::vector< std::pair< std::string, std::string > >::iterator itr =
	file_paths.begin(); itr != file_paths.end(); )
    {
      std::string file = itr->first;
      std::string path = itr->second;
      std::vector<unsigned> index;
      std::string::size_type pos1 = path.find_first_of("{")+1;
      std::string::size_type pos2 = path.find_first_of("}");
      std::string index_string = path.substr(pos1,pos2-pos1);
      VSDatumConverter< std::vector<unsigned> >::
	fromString(index,index_string.c_str(),"/");

      while(index.size() == 1 && pos1 != std::string::npos)
	{
	  pos1 = path.find_first_of("{",pos1)+1;
	  pos2 = path.find_first_of("}",pos2+1);
	  index_string = path.substr(pos1,pos2-pos1);

	  VSDatumConverter< std::vector<unsigned> >::
	    fromString(index,index_string.c_str(),"/");
	}
      

      if(index.size() > 1)
	{

	  file_paths.erase(itr);
	  for(unsigned i = 0; i < index.size(); i++)
	    {
	      std::string new_path = 
		path.substr(0,pos1) + 
		VSDataConverter::toString(index[i]) +
		path.substr(pos2,path.size()-pos2);
	      file_paths.push_back(std::make_pair(file,new_path));
	    }
	  itr = file_paths.begin();
	}
      else ++itr;
    }

  const unsigned nfile = file_paths.size();
  for(unsigned ifile = 0; ifile < nfile; ifile++)
    {
      std::string file = file_paths[ifile].first;
      std::string path = file_paths[ifile].second;

      
      

      std::cout << "Opening " << file << "/" << path << std::endl;

      if(VSOctaveH5ReaderStruct::isHDF5(file))
	load_hdf5(file,path,opt,hists_1d,hists_2d,graphs,elements);
      else
	load_txt(file,opt,hists_1d,hists_2d,graphs,elements);
    }

  for(unsigned ihist = 0; ihist < hists_1d.size(); ihist++)
    {
      if(opt.rebin > 1) hists_1d[ihist]->rebin(opt.rebin);

      if(opt.normalize == "area") hists_1d[ihist]->normalize();
      else if(opt.normalize == "max") hists_1d[ihist]->normalizeMax();

      if(!opt.cumulative.empty()) hists_1d[ihist]->cumulative(opt.cumulative);
    }


  if(opt.normalize == "area")
    {
      for(unsigned ihist = 0; ihist < hists_2d.size(); ihist++)
	hists_2d[ihist]->normalize();
    }
  else if(opt.normalize == "max")
    {
      for(unsigned ihist = 0; ihist < hists_2d.size(); ihist++)
	hists_2d[ihist]->normalizeMax();
    }

  if(opt.residual && hists_1d.size() > 1)
    {
      VSRHistogram1D hist = *hists_1d[0];

      for(unsigned ihist = 0; ihist < hists_1d.size(); ihist++)
	*hists_1d[ihist] -= hist;
    }
  else if(opt.ratio && hists_1d.size() > 1)
    {
      VSRHistogram1D hist = *hists_1d[0];

      for(unsigned ihist = 0; ihist < hists_1d.size(); ihist++)
	*hists_1d[ihist] /= hist;
    }
  else if(opt.ratio && graphs.size() > 1)
    {
      VSRGraph1D graph = *graphs[0];

      const unsigned ngraph = graphs.size();
      for(unsigned igr = 0; igr < ngraph; igr++)
	*graphs[igr] /= graph;

      graphs[0]->setShowErrors(false);
    }
 

  if(opt.residual && hists_2d.size() > 1)
    {
      VSRHistogram2D hist = *hists_2d[0];

      for(unsigned ihist = 0; ihist < hists_2d.size(); ihist++)
	*hists_2d[ihist] -= hist;
    }
  else if(opt.ratio && hists_2d.size() > 1)
    {
      VSRHistogram2D hist = *hists_2d[0];

      for(unsigned ihist = 0; ihist < hists_2d.size(); ihist++)
	*hists_2d[ihist] /= hist;
    }

  VSRCanvas canvas("test");
  unsigned ipad = 0;

  const unsigned nelements = elements.size();
  for(unsigned iel = 0; iel < nelements; iel++)
    {
      if(dump) elements[iel]->dump();

      if(iel <  opt.rescale.size()) 
	elements[iel]->multiply(opt.rescale[iel]);

      std::string label;
      if(iel < labels.size()) label = labels[iel];

      canvas.add(ipad,elements[iel],label);

      if(limits_x.first != limits_x.second)
	canvas.setLimitsX(ipad,limits_x.first,limits_x.second);

      if(limits_y.first != limits_y.second)
	canvas.setLimitsY(ipad,limits_y.first,limits_y.second);

      if(limits_z.first != limits_z.second)
	canvas.setLimitsZ(ipad,limits_z.first,limits_z.second);

      canvas.setBinLabels(ipad, bin_labels);

      //if(elements[iel]->is2D()) ipad++;
    }

  canvas.setTitle(0,title);
//   for(unsigned ihist = 0; ihist < hists_2d.size(); ihist++)
//     {
//       canvas.setTitle(ihist,title);
//       canvas.setGrid(ihist,true);
//       canvas.add(ihist,hists_2d[ihist]);


//     }
//   if(limits_y.first != limits_y.second)
//     canvas.setLimitsY(0,limits_y.first,limits_y.second);

//   if(limits_x.first != limits_x.second)
//     canvas.setLimitsX(0,limits_x.first,limits_x.second);

//   if(limits_z.first != limits_z.second)
//     canvas.setLimitsZ(0,limits_z.first,limits_z.second);

  canvas.draw();

  app->Run(true);

}

void load_txt(const std::string& file,  
	      Options& opt,
	      std::vector< VSRHistogram1D* >& hists_1d,
	      std::vector< VSRHistogram2D* >& hists_2d,
	      std::vector< VSRGraph1D* >& graphs,
	      std::vector< VSRCanvasElement* >& elements)
{
  
  std::ifstream inFile(file.c_str());

  double x, y;

  VSRGraph1D* graph = new VSRGraph1D;
  graphs.push_back(graph);
  elements.push_back(graph);
  
  while(inFile >> x >> y)
    {
      x = log10(x)-3;

      std::cout << x << " " << y << std::endl;
      graph->set(x,y);
    }
}

void load_hdf5(const std::string& file,  std::string path,
	       Options& opt,
	       std::vector< VSRHistogram1D* >& hists_1d,
	       std::vector< VSRHistogram2D* >& hists_2d,
	       std::vector< VSRGraph1D* >& graphs,
	       std::vector< VSRCanvasElement* >& elements)
{
  VSNSpaceOctaveH5IO nspace_io;

  VSOctaveH5Reader *reader = new VSOctaveH5Reader(file);
  vsassert(reader);

  if(!reader->isValid(path))
    {
      std::cout << "Invalid path: " 
		<< path << std::endl << std::endl;

      std::vector<std::string> fields;
      while(1)
	{
	  std::string::size_type pos = path.find_last_of(".");
	  if(pos == std::string::npos) 
	    {
	      path = "";
	      fields = reader->variables();
	      break;
	    }

	  path = path.substr(0,pos);
	      
	  if(reader->isValid(path))
	    {
	      VSOctaveH5ReaderStruct* s = reader->readStruct(path);
	      vsassert(s);
	      fields = s->variables();
	      delete s;
	      break;
	    }
	}

      std::cout << "Contents of " << file << "/" 
		<< path << ":" << std::endl;
      for(std::vector<std::string>::iterator ifield = 
	    fields.begin(); ifield != fields.end(); ifield++)
	std::cout << *ifield << std::endl;

      exit(EXIT_FAILURE);
    }

  
  if(reader->isVector(path))
    {
      std::vector<double> v;
      reader->readVector(path,v);
      std::sort(v.begin(),v.end());

      VSRHistogram1D* hist;
      
      double lo = 0;
      double hi = 0;

      if(opt.log)
	{
	  std::vector<double>::iterator itr = v.begin();
	  while(*itr <= 0)
	    {
	      itr++;
	      if(itr+1 == v.end()) break;
	    }
	  lo = log10(*itr);
	  hi = log10(*v.rbegin());
	}
      else
	{
	  lo = *v.begin();
	  hi = *v.rbegin();
	}

      if(opt.hist_range.first > 0)
	hist = new VSRHistogram1D(opt.hist_range.first,
				  opt.hist_range.second,
				  opt.hist_range.third);
      else
	hist = new VSRHistogram1D(200,lo,hi);
      
      for(unsigned i = 0; i < v.size(); i++)
	{
	  if(opt.log) hist->fill(log10(v[i]));
	  else hist->fill(v[i]);
	}

      hists_1d.push_back(hist);
      elements.push_back(hist);
      
      delete reader;
      return;
    }
  
  VSOctaveH5ReaderStruct* s = reader->readStruct(path);
  vsassert(s);

  VSLimitedHist<double,double> limited_hist;
  VSLimitedErrorsHist<double,double> limited_errors_hist;
  VSSimple2DHist<double,double> simple2d_hist;
  VSSimpleGraph<double,double> simple_graph;

  if(nspace_io.isHistogram(s) || nspace_io.isVolume(s))
    {
      VSNSpace nspace;

      if(nspace_io.isVolume(s))
	{
	  VSNSpace::Volume volume;
	  nspace_io.readVolume(s,volume);
	  nspace = VSNSpace(volume);
	}
      else
	nspace_io.readHistogram(s,nspace);

      std::vector< std::pair< unsigned, double > > slice = opt.nspace_slice;

      const unsigned nslice = slice.size();
      for(unsigned islice = 0; islice < nslice; islice++)
	{
	  vsassert(nspace.slice(slice[islice].first,
				slice[islice].second));

	  for(unsigned jslice = islice+1; jslice < nslice; jslice++)
	    if(slice[jslice].first > slice[islice].first) 
	      slice[jslice].first--;
	}

      for(std::vector<triple<unsigned,unsigned,unsigned> >::iterator itr = 
	    opt.nspace_marginalize.begin(); 
	  itr != opt.nspace_marginalize.end(); ++itr)
	{
	  VSNSpace::Volume marginal_volume(nspace.space());
	  marginal_volume.setInvert();
	  marginal_volume.clearInsideRange(itr->first,itr->second,
					   itr->third);
	  nspace.marginalize(itr->first,marginal_volume);
	}

      std::set< unsigned > dims;

      for(std::vector<unsigned>::iterator itr = opt.nspace_project.begin();
	  itr != opt.nspace_project.end(); ++itr)
	dims.insert(*itr);

      nspace.project(dims);
	  
      if(nspace.space().ndim >= 2)
	{
	  VSRHistogram2D* hist = VSRHistogram2D::create(nspace);
	  elements.push_back(hist);
	  hists_2d.push_back(hist);
	}
      else if(nspace.space().ndim == 1)
	{
	  VSRHistogram1D* hist = VSRHistogram1D::create(nspace);
	  elements.push_back(hist);
	  hists_1d.push_back(hist);
	}
    }
  else if(limited_hist.load(s))
    {	  
      VSRHistogram1D* hist = VSRHistogram1D::create(&limited_hist);
      hists_1d.push_back(hist);
      elements.push_back(hist);
    }
  else if(limited_errors_hist.load(s))
    {
      VSRHistogram1D* hist = VSRHistogram1D::create(&limited_errors_hist);
      hists_1d.push_back(hist);
      elements.push_back(hist);
    }
  else if(simple2d_hist.load(s))
    {
      VSRHistogram2D* hist = VSRHistogram2D::create(&simple2d_hist);

      for(std::vector< std::pair<unsigned, unsigned> >::iterator itr =
	    opt.projectx_range.begin(); itr != opt.projectx_range.end();
	  ++itr)
	{

//       if(opt.projectx)
// 	{
	  VSRHistogram1D* histp = hist->projectX(itr->first,itr->second);

	  if(opt.rebin > 1) histp->rebin(opt.rebin);
	  elements.push_back(histp);	  
	}
	  
      for(std::vector< std::pair<unsigned, unsigned> >::iterator itr =
	    opt.projecty_range.begin(); itr != opt.projecty_range.end();
	  ++itr)
	{
	  VSRHistogram1D* histp = hist->projectY(itr->first,itr->second);
	  if(opt.rebin > 1) histp->rebin(opt.rebin);
	  elements.push_back(histp);	  
	}

//       if(opt.projecty)
// 	{
// 	  VSRHistogram1D* histp = 
// 	    hist->projectY(opt.projecty_range.first,opt.projecty_range.second);

// 	  if(opt.rebin > 1) histp->rebin(opt.rebin);
// 	  elements.push_back(histp);	  
// 	}

      if(!opt.projectx && !opt.projecty)
	{
	  elements.push_back(hist);	  
	  hists_2d.push_back(hist);
	}
    }
  else if(simple_graph.load(s))
    {
      VSRGraph1D* graph = VSRGraph1D::create(&simple_graph);
      graphs.push_back(graph);
      elements.push_back(graph);
    }
 
  delete reader;
}
