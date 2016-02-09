#include <iostream>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cassert>

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
#include <THStack.h>
#include <TMultiGraph.h>
#include <TLegend.h>

// ----------------------------------------------------------------------------
// Local Includes
// ----------------------------------------------------------------------------
#include "VSRText.hpp"

using namespace std;


VSRText::VSRText(): m_text()
{

}

VSRText::~VSRText()
{
  for(std::vector< TText* >::iterator itr = m_text.begin(); 
      itr != m_text.end(); ++itr)
    delete *itr;
}

void VSRText::addText(const std::string& text)
{
  double y = 0.9-0.06*m_text.size();
  m_text.push_back(new TText(0.1,y,text.c_str()));
}

void VSRText::draw()
{
  for(std::vector< TText* >::iterator itr = m_text.begin(); 
      itr != m_text.end(); ++itr)
    (*itr)->Draw();
}
