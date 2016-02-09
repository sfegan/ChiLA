//-*-mode:c++; mode:font-lock;-*-
#ifndef VSRTEXT_HPP
#define VSRTEXT_HPP

#include <vector>

// ----------------------------------------------------------------------------
// ChiLA Includes
// ----------------------------------------------------------------------------
#include <VSSimpleHist.hpp>
#include <VSSimpleErrorsHist.hpp>
#include <VSNSpace.hpp>

// ----------------------------------------------------------------------------
// ROOT Includes
// ----------------------------------------------------------------------------
#include <TCanvas.h>
#include <TPad.h>
#include <TH1F.h>
#include <TText.h>

// ----------------------------------------------------------------------------
// Local Includes
// ----------------------------------------------------------------------------
#include "VSRHistogram.hpp"

class VSRText
{
public:
  VSRText();
  ~VSRText();

  void addText(const std::string& text);

  void draw();

private:

  std::vector< TText* > m_text;

};


#endif // VSRTEXT_HPP
