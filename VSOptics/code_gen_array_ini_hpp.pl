#!/usr/bin/perl -w

#
# Program: code_gen_array_ini_hpp.pl
# Author:  Stephen Fegan
#          Physics and Astronomy Department, UCLA
# E-mail:  sfegan@astro.ucla.edu
# Date:    2004/12/06
#
# The code generates reads in the file ARRAY_INI_TEMPLATE and generates
# the VSOArrayParameters.hpp header file.
#
# These two code generators read in the entries in an "array.ini" file
# (or the ARRAY_INI_TEMPLATE) and produce the C++ code required to
# read and write an array configuration to/from the "array.ini" file
# and transfer these values to/from the database. Initially the C++
# code was written by hand but it became very difficult to add values
# to or modify the structure of the "array.ini" file because the C++
# code would have to be changed by hand at the same time. This was
# very inconvenient and left the code prone to errors, so this code
# generator was developed. Hopefully this system will be relatively
# easy to maintain, although it does require the knowledge of PERL.
#
# The code proceeds as follows.
#
# (1) The file is read and the comments/members/units/types stored
#
# (2) The .hpp and .cpp files are generated using the stored data
#
# The ARRAY_INI_TEMPLATE should have entries as follows:
#
#   Comment comment comment comment comment comment comment comment
#   comment comment comment comment comment comment comment comment
# @ NameOfVariable [Units] <Variable Type>
#   Value
#
# For example:
#
#   Vector to intersection of rotation axes from the center of
#   reflector in reflector reference frame (+y points along optical axis,
#   +z is up, and +x points East when telescope is in Home position,
#   the origin is in the center of the reflector).
# @ ScopeTranslationX [cm] <double>
#   10.0
#

use strict;

my @Comments;
my @Members;
my @Types;
my @Units;
my @Values;

my @comment_bits;
my $ampersand=0;

foreach (<ARGV>)
  {
    chomp;
    s/\s+$//;

    if($ampersand==1)
      {
	$ampersand=0;
	s/^\s+//;
	s/\s+$//;
	push @Values,$_;
	next;
      }

    next unless $_;

    push(@comment_bits,$_);

    if(/^\@/)
      {
	my ( $amp, $var, $rest ) = split /\s+/,$_,3;
	my $unit="";
	my $type="double";
	
	$unit = $1 if($rest =~ /\[(.*)\]/);
	$type = $1 if($rest =~ /<(.*)>/);

	push @Comments,join("\\n",@comment_bits);
	push @Members,$var;
	push @Types,$type;
	push @Units,$unit;

	undef @comment_bits;

	$ampersand = 1;
      }
  }

print << 'END';
//-*-mode:c++; mode:font-lock;-*-

/*! \file VSOArrayParameters.hpp

  Array parameters class header file.

  AUTOMATICALLY GENERATED BY code_gen_array_ini_hpp.pl DO NOT EDIT!!!

  \author  Stephen Fegan               \n
           UCLA                        \n
	   sfegan@astro.ucla.edu       \n

  \date    11/30/2004
  \version 0.2

*/

#ifndef VSOARRAYPARAMETERS_HPP
#define VSOARRAYPARAMETERS_HPP

#include <iostream>
#include <string>

#include <Constants.hpp>
#include <Vec3D.hpp>

#include <VSDatabase.hpp>
#include <VSOptions.hpp>

namespace VERITAS
{

  //! Class encapsulates the parameters from which an array is generated.
  class VSOArrayParameters
  {
   public:
    // ************************************************************************
    // Array Parameters
    // ************************************************************************
END

my @members;
my @types;
my @units;
my $member;

@members=@Members;
@types=@Types;
@units=@Units;
while($member = shift @members)
  {
    my $type = shift @types;
    my $unit = shift @units;

    printf("    %-17s %-35s // %s\n",$type,$member.";",$unit);
  }

print << 'END';

    // ************************************************************************
    // Simple Member Functions
    // ************************************************************************
    VSOArrayParameters(bool use_canonical_values=false);
    virtual ~VSOArrayParameters();
    void dump(std::ostream& stream);
    void reset(bool use_canonical_values=false);

    // ************************************************************************
    // Command Line Member Functions
    // ************************************************************************
    static void zeroCanonicalValues();
    static void resetCanonicalValues();
    static void setCanonicalValuesFromArrayParameters(const VSOArrayParameters& o);
    static void setCanonicalValuesFromOptions(VSOptions& options);
    static void printOptions(std::ostream& stream);

    // ************************************************************************
    // ArrayINI Member Functions
    // ************************************************************************
    bool readFromArrayINIFile(std::istream& stream);
    bool readFromArrayINIFile(const std::string& filename);
    void writeToArrayINIFile(std::ostream& stream) const;
    void writeToArrayINIFile(const std::string& filename) const;

    // ************************************************************************
    // Database Member Functions
    // ************************************************************************
    bool readFromDatabase(VSDatabase* db, uint32_t optics_id);
    void writeToDatabase(VSDatabase* db, uint32_t optics_id) const;
    static void createSimulationParametersTable(VSDatabase* db);

    static const std::string scCollection;

   private:
END

@members=@Members;
@types=@Types;
@units=@Units;
while($member = shift @members)
  {
    my $type = shift @types;
    my $unit = shift @units;

    printf("    static %-10s %-35s // %s\n",$type,"sCanonical".$member.";",$unit);
  }

print << 'END';
  };
}

#endif // VSOARRAYPARAMETERS_HPP
END
