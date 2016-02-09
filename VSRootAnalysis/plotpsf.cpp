#include <iostream>

// ----------------------------------------------------------------------------
// ChiLA Includes
// ----------------------------------------------------------------------------
#include <VSOptions.hpp>
#include <VSSimpleStat.hpp>

// ----------------------------------------------------------------------------
// ROOT Includes
// ----------------------------------------------------------------------------
#include <TROOT.h>
#include <TMath.h>
#include <TH1F.h>
#include <TApplication.h>
#include <TPDF.h>
#include <TStyle.h>
#include <TH2F.h>
#include <TF1.h>
#include <TF2.h>
#include <THStack.h>
#include <TPaveText.h>
#include <Math/ParamFunctor.h>

// ----------------------------------------------------------------------------
// Local Includes
// ----------------------------------------------------------------------------
#include "VSRHistogram.hpp"
#include "VSRHistogram1D.hpp"
#include "VSRHistogram2D.hpp"
#include "VSRH5Loader.hpp"
#include "VSRCanvasBuilder.hpp"
#include "VSRSpectrumCalc.hpp"


#include "fitsio.h"

using namespace VERITAS;


class BkgndFunction : public ROOT::Math::ParamFunctor
{
public:
  BkgndFunction(double x, double y):m_x(x), m_y(y), m_exclude(true)
  { }

  virtual double operator() (const double* x, const double* p) const
  {
    double r = sqrt(std::pow(x[0]-m_x,2) + std::pow(x[1]-m_y,2));

    if (r < 0.2 && m_exclude) 
    {
      TF1::RejectPoint();
      return 0;
    }

    return p[0] + x[0]*p[1] + x[0]*x[0]*p[2] + x[1]*p[3] + x[1]*x[1]*p[4];
  }

  virtual BkgndFunction* Clone() const
  {
    return new BkgndFunction(m_x,m_y);
  }

  void setExclude(bool exclude) { m_exclude = exclude; }

private:
  double m_x;
  double m_y;
  bool   m_exclude;

};

bool open_fits_grab_vector(const std::string& filename, 
			   std::vector < std::vector < float > > &H, 
			   std::vector < std::string > &header, 
			   bool printheader = false);

void print_header(const std::vector< std::string >& header);
void find_peak(VSRHistogram2D* h, std::pair< unsigned, unsigned >& peak);
void usage(const std::string& progname, const VSOptions& options,
           std::ostream& stream)
{
  stream << "Usage: " << progname 
	 << " [options] stg3_file" 
	 << std::endl
         << std::endl;
  stream << "Options:" << std::endl;
  options.printUsage(stream);
}


int main(int argc, char **argv)
{
  std::string progname(*argv);
    

  VSOptions options(argc, argv, true);

  bool print_usage = false;
  if(options.find("h","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;
  if(options.find("help","Print this message.")!=VSOptions::FS_NOT_FOUND)
    print_usage=true;


  argv++,argc--;

  int arg_req = 1;
  if(argc != arg_req)
    {
      std::cerr << progname << ": need " << arg_req
		<< " arguments, got " << argc << std::endl;
      usage(progname, options, std::cout);
      exit(EXIT_SUCCESS);
    }


  std::string fits_file = argv[0];

  // --------------------------------------------------------------------------
  // Main program
  // --------------------------------------------------------------------------
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  TApplication *app = new TApplication("app",NULL,NULL);

  std::vector< std::vector< float > > fits_data;
  std::vector< std::string > fits_header;

  open_fits_grab_vector(fits_file,fits_data,fits_header);

  const unsigned ny = fits_data.size();
  const unsigned nx = fits_data[0].size();

  const double scale = 10.8/3600.0;

  double xmin = -scale*(nx/2);
  double xmax = scale*(nx/2);
  double ymin = -scale*(ny/2);
  double ymax = scale*(ny/2);



  VSRHistogram2D *h_raw = new VSRHistogram2D(nx,xmin,xmax,ny,ymin,ymax);
  VSRHistogram2D *h_signal = new VSRHistogram2D(nx,xmin,xmax,ny,ymin,ymax);
  VSRHistogram2D *h_smooth = new VSRHistogram2D(nx,xmin,xmax,ny,ymin,ymax);
  VSRHistogram2D *h_bkgnd = new VSRHistogram2D(nx,xmin,xmax,ny,ymin,ymax);
	



  for(unsigned ix = 0; ix < nx; ix++)
    {
      for(unsigned iy = 0; iy < ny; iy++)
	{
	  h_raw->set(ix,iy,fits_data[iy][ix]);
	}
    }
	

  *h_smooth = *h_raw;
  *h_signal = *h_raw;

  h_smooth->smooth(5);

  std::pair< unsigned, unsigned > peak_coord;

  find_peak(h_smooth,peak_coord);

  double peak_x = h_smooth->getValueX(peak_coord.first);
  double peak_y = h_smooth->getValueY(peak_coord.second);
  double xlo = peak_x-0.4;
  double xhi = peak_x+0.4;
  double ylo = peak_y-0.4;
  double yhi = peak_y+0.4;


//   std::cout << "PEAK X " << peak_coord.first << " " 
// 	    << h_smooth->getValueX(peak_coord.first) 
// 	    << " Y " << peak_coord.second << " "
// 	    << h_raw->getValueY(peak_coord.second) 
// 	    << std::endl;

  BkgndFunction bkgnd_functor(peak_x,peak_y);

  TF2 *bkgnd_func = new TF2("bkgnd_func",bkgnd_functor,xlo,xhi,ylo,yhi,5);


  h_raw->fit("bkgnd_func");

  
  double a = bkgnd_func->GetParameter(0);
  double b = bkgnd_func->GetParameter(1);
  double c = bkgnd_func->GetParameter(2);

  VSSimpleStat2<double,double> stat_x;
  VSSimpleStat2<double,double> stat_y;
  VSSimpleStat2<double,double> stat_xx;
  VSSimpleStat2<double,double> stat_xy;
  VSSimpleStat2<double,double> stat_yy;

  bkgnd_functor.setExclude(false);

  TF2 *bkgnd_func2 = new TF2("bkgnd_func2",bkgnd_functor,xlo,xhi,ylo,yhi,5);

  bkgnd_func2->SetParameter(0,bkgnd_func->GetParameter(0));
  bkgnd_func2->SetParameter(1,bkgnd_func->GetParameter(1));
  bkgnd_func2->SetParameter(2,bkgnd_func->GetParameter(2));
  bkgnd_func2->SetParameter(3,bkgnd_func->GetParameter(3));
  bkgnd_func2->SetParameter(4,bkgnd_func->GetParameter(4));

  for(unsigned ix = 0; ix < h_signal->getNBinsX(); ix++)
    {
      for(unsigned iy = 0; iy < h_signal->getNBinsY(); iy++)
	{
	  double x = h_signal->getValueX(ix);
	  double y = h_signal->getValueY(iy);
	  double z = h_signal->get(ix,iy);
	  double bkgnd = bkgnd_func2->Eval(x,y);

	  h_signal->set(ix,iy,z-bkgnd);

	  if(x < xlo || x > xhi)
	    continue;
	  else if(y < ylo || y > yhi)
	    continue;

	  stat_x.accumulate(x,z-bkgnd);
	  stat_y.accumulate(y,z-bkgnd);

	  stat_xx.accumulate(x*x,z-bkgnd);
	  stat_yy.accumulate(y*y,z-bkgnd);
	  stat_xy.accumulate(x*y,z-bkgnd);
	}
    }

  double mx = stat_x.mean();
  double my = stat_y.mean();
  double mxx = stat_xx.mean() - mx*mx;
  double myy = stat_yy.mean() - my*my;
  double sum =  stat_x.count();

  std::cout << "SUM:    " << sum << std::endl;
  std::cout << "MEAN X: " << mx << std::endl;
  std::cout << "MEAN Y: " << my << std::endl;
  //  std::cout << "RMS:    " << sqrt(mxx+myy) << std::endl;
  //  std::cout << "FWHM:   " << 1.67*sqrt(mxx+myy) << std::endl;

  VSRHistogram1D* h_th = new VSRHistogram1D(200,0,0.3);

  for(unsigned ix = 0; ix < h_signal->getNBinsX(); ix++)
    {
      for(unsigned iy = 0; iy < h_signal->getNBinsY(); iy++)
	{
	  double x = h_signal->getValueX(ix);
	  double y = h_signal->getValueY(iy);

	  double r1 = sqrt(std::pow(x-mx,2)+std::pow(y-my,2));

	  for(int ir = 0;ir < h_th->getNBins();ir++)
	    {
	      double r2 = h_th->getX(ir);
	      
	      if(r1 < r2) 
		h_th->set(ir,h_th->get(ir)+h_signal->get(ix,iy));

	    }

	}
    }
  

  double s1 = h_th->get(h_th->getNBins()-1);

  //  h_th->normalize();
  h_th->multiply(1./s1);

  TGraph *gr_th = new TGraph;

  for(unsigned ix = 0; ix < h_th->getNBins(); ix++)
    {
      gr_th->SetPoint(ix,h_th->get(ix),h_th->getX(ix));
    }
  
  std::cout << "68% Containment: " << 2*gr_th->Eval(0.68) << std::endl;
  std::cout << "80% Containment: " << 2*gr_th->Eval(0.80) << std::endl;


  VSRCanvasBuilder canvas_raw;
  VSRCanvasBuilder canvas_signal;
  VSRCanvasBuilder canvas_th;

  canvas_signal.addElement(0,h_signal);
  canvas_raw.addElement(0,h_raw);
  canvas_th.addElement(0,h_th);
  //  canvas.addElement(1,h_signal);
  // canvas.addElement(2,h_signal);

  h_raw->setRange(xlo,xhi,ylo,yhi);
  h_signal->setRange(xlo,xhi,ylo,yhi);



  canvas_signal.createCanvas("c1")->Draw();

  gPad->SetGrid();

  canvas_raw.createCanvas("c2")->Draw();
  canvas_th.createCanvas("c3")->Draw();

  

  //  print_header(fits_header);


  app->Run(true);
}

void find_peak(VSRHistogram2D* h, std::pair< unsigned, unsigned >& peak)
{
  double max = 0;

  for(unsigned ix = 0; ix < h->getNBinsX(); ix++)
    for(unsigned iy = 0; iy < h->getNBinsY(); iy++)
      {
	if(h->get(ix,iy) > max)
	  {
	    max = h->get(ix,iy);
	    peak.first=ix;
	    peak.second=iy;
	  }
      }

}

void print_header(const std::vector< std::string >& header)
{
  for(std::vector< std::string >::const_iterator itr = header.begin();
      itr != header.end(); ++itr)
    std::cout << *itr << std::endl;
}

bool open_fits_grab_vector(const std::string& filename, 
			   std::vector < std::vector < float > > &H, 
			   std::vector < std::string > &header, 
			   bool printheader)
{
  bool ioerror = false;

  //Initialisations
  H.clear();
  header.clear();

  int mydex = 0;
  std::vector < float > h;
  fitsfile *fptr;
  int status,  nfound, anynull;
  long naxes[2], fpixel, nbuffer, npixels, iii;
  int ii, nkeys;
  status = 0;
  int buffsize =  1000;
  float datamin, datamax, nullval, buffer[buffsize];
  char card[FLEN_CARD]; 
  std::string temp = ""; 
  status = 0;

  //Open fits file
  if ( fits_open_file(&fptr,filename.c_str(), READONLY, &status) )
    {
      fits_report_error(stderr, status);
      ioerror = true;
    }

  //Grab header
  fits_get_hdrspace(fptr, &nkeys, NULL, &status);
  
  for (ii = 1; ii <= nkeys; ii++)
    { 
      fits_read_record(fptr, ii, card, &status); /* read keyword */
      if(printheader){std::cout << card << std::endl;}
      temp = std::string(card);
      header.push_back(temp);
    }
  
  /* read the NAXIS1 and NAXIS2 keyword to get image size */
  if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
    fits_report_error(stderr, status);

  npixels  = naxes[0] * naxes[1];         /* number of pixels in the image */
  fpixel   = 1;
  nullval  = 0;                /* don't check for null values in the image */
  datamin  = 1.0E30;
  datamax  = -1.0E30;
  std::vector < std::vector < float > > G;

  while (npixels > 0)
    {
      nbuffer = npixels;
      if (npixels > buffsize)
	nbuffer = buffsize;     // read as many pixels as will fit in buffer

      // Note that even though the FITS images contains unsigned integer 
      // pixel values (or more accurately, signed integer pixels with    
      // a bias of 32768),  this routine is reading the values into a    
      // float array.   Cfitsio automatically performs the datatype      
      // conversion in cases like this.                                  

      if ( fits_read_img(fptr, TFLOAT, fpixel, nbuffer, &nullval,
			 buffer, &anynull, &status) )
	fits_report_error(stderr, status);

      for (iii = 0; iii < nbuffer; iii++)  
	{       
	  if(mydex < naxes[0] -1)
	    {
	      h.push_back(buffer[iii]);
	      mydex++;
	    }
	  else
	    {
	      G.push_back(h);
	      h.clear();
	      mydex = 0;
	    }	
      }
      npixels -= nbuffer;    // increment remaining number of pixels
      fpixel  += nbuffer;    // next pixel to be read in image
    }

  for(unsigned int j = 0; j < G.size(); j++)
    H.push_back(G[G.size() - 1 - j]);
	
  fits_close_file(fptr, &status);

  if(status)fits_report_error(stderr, status);

  return ioerror;
}
