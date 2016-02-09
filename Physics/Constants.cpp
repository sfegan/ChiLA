#include"Constants.hpp"

using namespace Physics;

const double Constants::num_pi;
const double Constants::num_e;
const double Constants::num_ln_10;
const double Constants::num_log10_e;
const double Constants::num_root_two;
const double Constants::num_root_half;

const double Constants::cgs_L;
const double Constants::cgs_M;
const double Constants::cgs_T;

const double Constants::cgs_L2;
const double Constants::cgs_L3;
const double Constants::cgs_M2;
const double Constants::cgs_M3;
const double Constants::cgs_T2;
const double Constants::cgs_T3;
const double Constants::cgs_V;
const double Constants::cgs_F;
const double Constants::cgs_E;
const double Constants::cgs_Q;
const double Constants::cgs_FLUX;
const double Constants::cgs_E_FLUX;

const double Constants::num_avogadro;
const double Constants::num_avogadro_err;
const double Constants::num_alpha;
const double Constants::num_alpha_err;

const double Constants::si_c;
const double Constants::si_mu0;
const double Constants::si_epsilon0;
const double Constants::si_c_err;
const double Constants::si_h;
const double Constants::si_h_err;
const double Constants::si_hbar;
const double Constants::si_hbar_err;
const double Constants::si_G;
const double Constants::si_G_err;
const double Constants::si_e;
const double Constants::si_e_err;
const double Constants::si_k;
const double Constants::si_k_err;
const double Constants::si_sigma;
const double Constants::si_sigma_err;
    
const double Constants::cgs_c;
const double Constants::cgs_c_err;
const double Constants::cgs_h;
const double Constants::cgs_h_err;
const double Constants::cgs_hbar;
const double Constants::cgs_hbar_err;
const double Constants::cgs_G;
const double Constants::cgs_G_err;
const double Constants::cgs_e;
const double Constants::cgs_e_err;
const double Constants::cgs_k;
const double Constants::cgs_k_err;
const double Constants::cgs_sigma;
const double Constants::cgs_sigma_err;

const double Constants::si_amu;
const double Constants::si_amu_err;
const double Constants::si_eV;
const double Constants::si_eV_err;

const double Constants::cgs_amu;
const double Constants::cgs_amu_err;
const double Constants::cgs_eV;
const double Constants::cgs_eV_err;
    
const double Constants::si_e_chg;
const double Constants::si_e_chg_err;
const double Constants::si_e_mass;
const double Constants::si_e_mass_err;
const double Constants::si_p_chg;
const double Constants::si_p_chg_err;
const double Constants::si_p_mass;
const double Constants::si_p_mass_err;
 
const double Constants::cgs_e_chg;
const double Constants::cgs_e_chg_err;
const double Constants::cgs_e_mass;
const double Constants::cgs_e_mass_err;
const double Constants::cgs_p_chg;
const double Constants::cgs_p_chg_err;
const double Constants::cgs_p_mass;
const double Constants::cgs_p_mass_err;

const double Constants::si_tropical_year;
const double Constants::si_julian_year;

const double Constants::cgs_tropical_year;
const double Constants::cgs_julian_year;

const double Constants::si_au;
const double Constants::si_light_year;
const double Constants::si_parsec;
const double Constants::si_sun_mass;
const double Constants::si_sun_radius;
const double Constants::si_sun_luminosity;

const double Constants::cgs_au;
const double Constants::cgs_light_year;
const double Constants::cgs_parsec;
const double Constants::cgs_sun_mass;
const double Constants::cgs_sun_radius;
const double Constants::cgs_sun_luminosity;

template<> std::string Physics::convertToString<>(const std::string& x)
{
  std::ostringstream stream;
  stream << '"';
  for(std::string::const_iterator i = x.begin(); i!=x.end(); i++)
    {
      if(*i == '"')stream << '\\';
      stream << *i;
    }
  stream << '"';
  return stream.str();
}

template<> std::string Physics::convertToString<>(const float& x)
{
  std::ostringstream stream;
  stream << std::setprecision(std::numeric_limits<float>::digits10) 
	 << std::scientific << x;
  return stream.str();
}

template<> std::string Physics::convertToString<>(const double& x)
{
  std::ostringstream stream;
  stream << std::setprecision(std::numeric_limits<double>::digits10) 
	   << std::scientific << x;
  return stream.str();
}

template<> std::string Physics::convertToString<>(const long double& x)
{
    std::ostringstream stream;
    stream << std::setprecision(std::numeric_limits<long double>::digits10) 
	   << std::scientific << x;
    return stream.str();
}

enum ParserState { PS_SKIP_SPACE, PS_NORMAL, PS_ESCAPED, PS_DOUBLE_QUOTED };

template<> bool Physics::convertFromString<>(const std::string& str, std::string& x)
{
  x.clear();

  std::vector<ParserState> state;
  state.push_back(PS_SKIP_SPACE);

  for(std::string::const_iterator i = str.begin(); !state.empty() && i!=str.end(); i++)
    {
      switch(state[state.size()-1])
	{
	case PS_SKIP_SPACE:
	  if(isspace(*i))continue;
	  state.pop_back();
	  state.push_back(PS_NORMAL);
	  // FALL THROUGH
	case PS_NORMAL:
	  if(*i=='\\')state.push_back(PS_ESCAPED);
	  else if(*i=='"')state.push_back(PS_DOUBLE_QUOTED);
	  else if(isspace(*i))state.pop_back();
	  else x.push_back(*i);
	  break;
	case PS_ESCAPED:
	  if(*i=='n')x.push_back('\n');
	  else if(*i=='r')x.push_back('\r');
	  else x.push_back(*i);
	  state.pop_back();
	  break;
	case PS_DOUBLE_QUOTED:
	  if(*i=='\\')state.push_back(PS_ESCAPED);
	  else if(*i=='"')state.pop_back();
	  else x.push_back(*i);
	  break;
	}
    }
  return true;
}

#ifdef TEST_MAIN

#include<sstream>
#include<iostream>
#include<iomanip>
#include<string>

std::string fmt(double val, const std::string& units="")
{
  std::ostringstream s;
  s << std::scientific << std::setprecision(10) << val << ' ' << units;
  std::ostringstream s2;
  s2 << std::setw(40) << std::left << s.str();
  return s2.str();
}

int main()
{
  int w=20;

  std::cout << std::setw(w) << std::left << "pi"              << fmt(Constants::num_pi) << std::endl
	    << std::setw(w) << std::left << "e"               << fmt(Constants::num_e) << std::endl
	    << std::setw(w) << std::left << "ln10"            << fmt(Constants::num_ln_10) << std::endl
	    << std::setw(w) << std::left << "log10(e)"        << fmt(Constants::num_log10_e) << std::endl
	    << std::setw(w) << std::left << "sqrt(2)"         << fmt(Constants::num_root_two) << std::endl
	    << std::setw(w) << std::left << "sqrt(0.5)"       << fmt(Constants::num_root_half) << std::endl
	    << std::endl
#if 1
	    << std::setw(w) << std::left << "cgs L"           << fmt(Constants::cgs_L) << std::endl
	    << std::setw(w) << std::left << "cgs M"           << fmt(Constants::cgs_M) << std::endl
	    << std::setw(w) << std::left << "cgs T"           << fmt(Constants::cgs_T) << std::endl
	    << std::setw(w) << std::left << "cgs L2"          << fmt(Constants::cgs_L2) << std::endl
	    << std::setw(w) << std::left << "cgs L3"          << fmt(Constants::cgs_L3) << std::endl
	    << std::setw(w) << std::left << "cgs M2"          << fmt(Constants::cgs_M2) << std::endl
 	    << std::setw(w) << std::left << "cgs M3"          << fmt(Constants::cgs_M3) << std::endl
	    << std::setw(w) << std::left << "cgs T2"          << fmt(Constants::cgs_T2) << std::endl
	    << std::setw(w) << std::left << "cgs T3"          << fmt(Constants::cgs_T3) << std::endl
	    << std::setw(w) << std::left << "cgs V"           << fmt(Constants::cgs_V) << std::endl
	    << std::setw(w) << std::left << "cgs F"           << fmt(Constants::cgs_F) << std::endl
	    << std::setw(w) << std::left << "cgs E"           << fmt(Constants::cgs_E) << std::endl
	    << std::setw(w) << std::left << "cgs Q"           << fmt(Constants::cgs_Q,"g cm^3 / A s^3") << std::endl
	    << std::endl
#endif
	    << std::setw(w) << std::left << "Avogadro"        << fmt(Constants::num_avogadro) << std::endl
	    << std::setw(w) << std::left << "Avogadro err"    << fmt(Constants::num_avogadro_err) << std::endl
	    << std::setw(w) << std::left << "Fine struct"     << fmt(Constants::num_alpha) << std::endl
	    << std::setw(w) << std::left << "Fine struct err" << fmt(Constants::num_alpha_err) << std::endl
	    << std::endl

	    << std::setw(w) << std::left << "mu0"             << fmt(Constants::si_mu0,"kg m / A^2 s^2") << std::endl
	    << std::setw(w) << std::left << "epsilon0"        << fmt(Constants::si_epsilon0,"A^2 s^4 / kg m^3") << std::endl
    
	    << std::setw(w) << std::left << "c"             
	    << fmt(Constants::si_c,"m / s") << ' ' << fmt(Constants::cgs_c,"cm / s") << std::endl

	    << std::setw(w) << std::left << "h"             
	    << fmt(Constants::si_h,"kg m^2 / s") << ' ' << fmt(Constants::cgs_h,"g cm^2 / s") << std::endl

	    << std::setw(w) << std::left << "h err"             
	    << fmt(Constants::si_h_err,"kg m^2 / s") << ' ' << fmt(Constants::cgs_h_err,"g cm^2 / s") << std::endl

	    << std::setw(w) << std::left << "h bar"             
	    << fmt(Constants::si_hbar,"kg m^2 / s") << ' ' << fmt(Constants::cgs_hbar,"g cm^2 / s") << std::endl

	    << std::setw(w) << std::left << "h bar err"             
	    << fmt(Constants::si_hbar_err,"kg m^2 / s") << ' ' << fmt(Constants::cgs_hbar_err,"g cm^2 / s") << std::endl

	    << std::setw(w) << std::left << "G"             
	    << fmt(Constants::si_G,"m^3 / kg / s^2") << ' ' << fmt(Constants::cgs_G,"cm^3 / g / s^2") << std::endl

	    << std::setw(w) << std::left << "G err"             
	    << fmt(Constants::si_G_err,"m^3 / kg / s^2") << ' ' << fmt(Constants::cgs_G_err,"cm^3 / g / s^2") << std::endl

	    << std::setw(w) << std::left << "e"             
	    << fmt(Constants::si_e,"A s") << ' ' << fmt(Constants::cgs_e,"esu (g^1/2 cm^3/2 / s)") << std::endl

	    << std::setw(w) << std::left << "e err"             
	    << fmt(Constants::si_e_err,"A s") << ' ' << fmt(Constants::cgs_e_err,"esu (g^1/2 cm^3/2 / s)") << std::endl

	    << std::setw(w) << std::left << "k"             
	    << fmt(Constants::si_k,"kg m^2 / K / s^2") << ' ' << fmt(Constants::cgs_k,"g cm^2 / K / s^2") << std::endl

	    << std::setw(w) << std::left << "k err"             
	    << fmt(Constants::si_k_err,"kg m^2 / K / s^2") << ' ' << fmt(Constants::cgs_k_err,"g cm^2 / K / s^2") << std::endl

	    << std::setw(w) << std::left << "stef-boltz"             
	    << fmt(Constants::si_sigma,"kg / K^4 / s^3") << ' ' << fmt(Constants::cgs_sigma,"g / K^4 / s^3") << std::endl

	    << std::setw(w) << std::left << "stef-boltz err"
	    << fmt(Constants::si_sigma_err,"kg / K^4 / s^3") << ' ' << fmt(Constants::cgs_sigma_err,"g / K^4 / s^3") << std::endl
	    << std::endl

	    << std::setw(w) << std::left << "A.M.U."             
	    << fmt(Constants::si_amu,"kg") << ' ' << fmt(Constants::cgs_amu,"g") << std::endl

	    << std::setw(w) << std::left << "A.M.U. err"             
	    << fmt(Constants::si_amu_err,"kg") << ' ' << fmt(Constants::cgs_amu_err,"g") << std::endl

	    << std::setw(w) << std::left << "eV"             
	    << fmt(Constants::si_eV,"kg m^2 / s^2") << ' ' << fmt(Constants::cgs_eV,"g cm^2 / s^2") << std::endl

	    << std::setw(w) << std::left << "eV err"             
	    << fmt(Constants::si_eV_err,"kg m^2 / s^2") << ' ' << fmt(Constants::cgs_eV_err,"g cm^2 / s^2") << std::endl
	    << std::endl

	    << std::setw(w) << std::left << "e-charge"             
	    << fmt(Constants::si_e_chg,"A s") << ' ' << fmt(Constants::cgs_e_chg,"esu (g^1/2 cm^3/2 / s)") << std::endl

	    << std::setw(w) << std::left << "e-charge err"             
	    << fmt(Constants::si_e_chg_err,"A s") << ' ' << fmt(Constants::cgs_e_chg_err,"esu (g^1/2 cm^3/2 / s)") << std::endl

	    << std::setw(w) << std::left << "e-mass"             
	    << fmt(Constants::si_e_mass,"kg") << ' ' << fmt(Constants::cgs_e_mass,"g") << std::endl

	    << std::setw(w) << std::left << "e-mass err"             
	    << fmt(Constants::si_e_mass_err,"kg") << ' ' << fmt(Constants::cgs_e_mass_err,"g") << std::endl

	    << std::setw(w) << std::left << "p-charge"             
	    << fmt(Constants::si_p_chg,"A s") << ' ' << fmt(Constants::cgs_p_chg,"esu (g^1/2 cm^3/2 / s)") << std::endl

	    << std::setw(w) << std::left << "p-charge err"             
	    << fmt(Constants::si_p_chg_err,"A s") << ' ' << fmt(Constants::cgs_p_chg_err,"esu (g^1/2 cm^3/2 / s)") << std::endl

	    << std::setw(w) << std::left << "p-mass"             
	    << fmt(Constants::si_p_mass,"kg") << ' ' << fmt(Constants::cgs_p_mass,"g") << std::endl

	    << std::setw(w) << std::left << "p-mass err"             
	    << fmt(Constants::si_p_mass_err,"kg") << ' ' << fmt(Constants::cgs_p_mass_err,"g") << std::endl
	    << std::endl

	    << std::setw(w) << std::left << "tropical year"   << fmt(Constants::si_tropical_year,"s") << std::endl
	    << std::setw(w) << std::left << "julian year"   << fmt(Constants::si_julian_year,"s") << std::endl
	    << std::endl

	    << std::setw(w) << std::left << "astronomical unit"             
	    << fmt(Constants::si_au,"m") << ' ' << fmt(Constants::cgs_au,"cm") << std::endl

	    << std::setw(w) << std::left << "light year"             
	    << fmt(Constants::si_light_year,"m") << ' ' << fmt(Constants::cgs_light_year,"cm") << std::endl

	    << std::setw(w) << std::left << "parsec"             
	    << fmt(Constants::si_parsec,"m") << ' ' << fmt(Constants::cgs_parsec,"cm") << std::endl

	    << std::setw(w) << std::left << "solar mass"             
	    << fmt(Constants::si_sun_mass,"kg") << ' ' << fmt(Constants::cgs_sun_mass,"g") << std::endl

	    << std::setw(w) << std::left << "solar radius"             
	    << fmt(Constants::si_sun_radius,"m") << ' ' << fmt(Constants::cgs_sun_radius,"cm") << std::endl

	    << std::setw(w) << std::left << "solar luminosity"             
	    << fmt(Constants::si_sun_luminosity,"kg m^2 / s^3") << ' ' << fmt(Constants::cgs_sun_luminosity,"g cm^2 / s^3") << std::endl
    ;
}

#endif 
