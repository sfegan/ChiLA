#include <iostream>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cmath>

// ----------------------------------------------------------------------------
// ChiLA Includes
// ----------------------------------------------------------------------------
#include <WhippleCams.h>

// ----------------------------------------------------------------------------
// ROOT Includes
// ----------------------------------------------------------------------------
#include <TEllipse.h>
#include <TPaletteAxis.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TBox.h>
#include <TGaxis.h>
#include <TPaveText.h>
#include <TMath.h>

// ----------------------------------------------------------------------------
// Local Includes
// ----------------------------------------------------------------------------
#include "VSRHistogramCamera.hpp"

using namespace std;

VSRHistogramCamera::VSRHistogramCamera(double lo, double hi,
				       const std::string& title):
  VSRHistogram(), m_lo(lo), m_hi(hi),
  m_title(title), m_pixelValues(500,0), m_isSuppressed(500,true)
{

}

VSRHistogramCamera::~VSRHistogramCamera() {

}

void VSRHistogramCamera::draw()
{
  if(m_lo == m_hi) 
    {
      m_lo = *min_element(m_pixelValues.begin(),m_pixelValues.end());
      m_hi = *max_element(m_pixelValues.begin(),m_pixelValues.end());
//       double median = getMedian(m_pixelValues,m_isSuppressed);
//       m_lo = median*0.6;
//       m_hi = median*1.4;
    }

  gStyle->SetPalette(1);
  int ncolors = gStyle->GetNumberOfColors();
  const double alpha = 4.0;
  const unsigned nchan = 499;
  unsigned numSuppressed = 0;

  for(unsigned int ichan = 0; ichan < nchan; ichan++) 
    {
      double x = VC499Xcoord[ichan]/alpha + 0.45;
      double y = VC499Ycoord[ichan]/alpha + 0.45;
      double r = VC499Radius[ichan]/alpha*1.2;
      
      TEllipse *e1 = new TEllipse(x,y,r);
      double value;
      
      if(m_pixelValues[ichan] > m_hi)
	value = 1;
      else if(m_pixelValues[ichan] < m_lo)
	value = 0;
      else
	value = (m_pixelValues[ichan]-m_lo)/(m_hi-m_lo);
      
      int theColor = int(0.99*value*double(ncolors));
      
      if(m_isSuppressed[ichan]) 
	{
	  numSuppressed++;
	  e1->SetFillColor(0);
	} 
      else
	e1->SetFillColor(gStyle->GetColorPalette(theColor));
      
      e1->SetLineWidth(1);
      e1->Draw();
      
      ostringstream os;
      os << ichan+1 << endl;
     
      if(m_options.showChannelLabels)
	{
	  TText *label = new TText(x,y-r/10.,os.str().c_str());
	  //      label->SetTextFont(51);
	  label->SetTextSize(0.015);
	  label->SetTextAlign(22);
	  label->Draw();
	}      
    }


//   double ymin = 0.55;
//   double ymax = 0.95;
//   double xmin = 0.88;
//   double xmax = 0.92;
  
//   double y1,y2;
//   for(unsigned icolor = 0; icolor < ncolors; icolor++) 
//     {
//       y1 = ymin + (float)icolor*(ymax-ymin)/(float)ncolors;
//       y2 = ymin + (float)(icolor+1)*(ymax-ymin)/(float)ncolors;
//       TBox *box = new TBox(xmin,y1,xmax,y2);
//       box->SetFillColor(gStyle->GetColorPalette(icolor));
//       box->Draw();
//     }

//   TBox *box = new TBox(xmin,ymin,xmax,ymax);
//   box->SetFillStyle(0);
//   box->Draw();

//   TGaxis *axis = new TGaxis(xmax,ymin,xmax,ymax,m_lo,m_hi,510,"+L");
//   //axis->SetLabelOffset(0.05);
  
//   axis->SetLineWidth((Width_t)1);
//   axis->SetLabelSize(0.03);
//   axis->Draw();

//   TPaveText *text = new TPaveText(0.05,0.9,0.3,0.95);
//   text->AddText(m_title.c_str());
//   text->SetBorderSize(0);
//   text->SetFillColor(0);

//   text->Draw();

//   TPaveText *text2 = new TPaveText(0.03,0.1,0.22,0.22);
//   ostringstream os;
//   os << "Median: " << setprecision(4) << getMedian(m_pixelValues,
// 						   m_isSuppressed);
//   text2->AddText(os.str().c_str());
//   os.str("");
//   os<< "Mean: " << setprecision(4) << getMean(m_pixelValues,
// 					      m_isSuppressed);
//   text2->AddText(os.str().c_str());
//   os.str("");
//   os<< "Num Suppressed: " << numSuppressed;
//   text2->AddText(os.str().c_str());
//   os.str("");
  

//   text2->SetBorderSize(0);
//   text2->SetFillColor(0);
//   text2->SetTextAlign(13);
//   text2->SetTextSize(0.03);
  //  text2->Draw();
}

VSRHistogramCamera* VSRHistogramCamera::clone() const
{
  VSRHistogramCamera* hist = new VSRHistogramCamera(*this);
  return hist;
}

double VSRHistogramCamera::getMedian(vector< double > &pixelValues,
				     vector< bool > &isSuppressed) 
{
  const unsigned nchan = pixelValues.size();

  double values[nchan];

  unsigned n = 0;
  for(unsigned ichan = 0; ichan < nchan; ichan++) 
    {
      if(!isSuppressed[ichan]) 
	{      
	  values[n] = pixelValues[ichan];
	  n++;
	}
    }
  
  return TMath::Median(n,values);
}

double VSRHistogramCamera::getMean(vector< double > &pixelValues,
				   vector< bool > &isSuppressed) 
{
  const unsigned nchan = pixelValues.size();

  double values[nchan];

  unsigned n = 0;
  for(unsigned ichan = 0; ichan < nchan; ichan++) {

    if(!isSuppressed[ichan]) {      
      values[n] = pixelValues[ichan];
      n++;
    }

  }
  
  return TMath::Mean(n,values);
}

void VSRHistogramCamera::smooth(double radius)
{

  double median = getMedian(m_pixelValues,m_isSuppressed);

  const unsigned nchan = 499;

  std::vector< double > channel_values(nchan);

  for(unsigned ichan1 = 0; ichan1 < nchan; ichan1++)
    {
      double n = 0;
      for(unsigned ichan2 = 0; ichan2 < nchan; ichan2++)
	{
	  double dr = 
	    sqrt(std::pow(VC499Xcoord[ichan1]-VC499Xcoord[ichan2],2) +
		 std::pow(VC499Ycoord[ichan1]-VC499Ycoord[ichan2],2));

	  double diff = 
	    fabs(m_pixelValues[ichan2]-median)/median;

	  if(ichan1 == ichan2 ||
	     (dr < radius*VC499Radius[ichan1] && !m_isSuppressed[ichan2] &&
	      diff < 0.5))
	    {
	      channel_values[ichan1] += m_pixelValues[ichan2];
	      n++;
	    }
	}

      channel_values[ichan1]/=n;
    }

  m_pixelValues = channel_values;
}

// TPaletteAxis* VHistogramCamera::createPalette(float min, float max) {
//   TCanvas *ctemp = new TCanvas("ctemp","ctemp");
//   TH2F *htemp = new TH2F("htemp","htemp",100,0,1,100,0,1);
//   htemp->SetMinimum(min);
//   htemp->SetMaximum(max);
//   htemp->Fill(0.5,0.5);
//   htemp->SetContour(50);
//   htemp->Draw("COLZ");
//   ctemp->Paint();
//   TPaletteAxis *palette = 
//     (TPaletteAxis*)htemp->GetListOfFunctions()->FindObject("palette")->Clone();
//   htemp->SetDirectory(0);
//   delete htemp;
//   ctemp->Close();
//   delete ctemp;
//   return palette;
// }
