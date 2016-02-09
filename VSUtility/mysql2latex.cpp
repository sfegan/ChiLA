//-*-mode:c++; mode:font-lock;-*-

#include<iostream>
#include<iomanip>
#include<vector>
#include<string>
#include<algorithm>
#include<cstring>

#include<VSOptions.hpp>
#include<VSDBMySQL3X.hpp>
#include<VSDBMySQL3X.hpp>
#include<VSDBFactory.hpp>

using namespace VERITAS;

class DBWriter
{
public:
  virtual ~DBWriter();
  virtual void startDocument(const std::string& database) = 0;
  virtual void endDocument() = 0;
  virtual void startTableListingSection() = 0;
  virtual void endTableListingSection() = 0;
  virtual void startTableListing() = 0;
  virtual void listTable(const std::string& table) = 0;
  virtual void endTableListing() = 0;
  virtual void startTableDescriptionSection() = 0;
  virtual void endTableDescriptionSection() = 0;
  virtual void startTableDescription(const std::string& tablename) = 0;
  virtual void listTableEntry(const std::string& field, 
			      const std::string& type,
			      const std::string& null, 
			      const std::string& key,
			      const std::string& def, 
			      const std::string& extra) = 0;
  virtual void endTableDescription() = 0;
};

DBWriter::~DBWriter()
{
  // nothing to see here
}

class DBLaTeXWriter: public DBWriter
{
public:
  DBLaTeXWriter(std::ostream& stream): 
    fStream(stream), fTableMembers(), fTableName(), fTableRestart(),
    fAbbrev() { }
  virtual ~DBLaTeXWriter();
  virtual void startDocument(const std::string& database);
  virtual void endDocument();
  virtual void startTableListingSection();
  virtual void endTableListingSection();
  virtual void startTableListing();
  virtual void listTable(const std::string& table);
  virtual void endTableListing();
  virtual void startTableDescriptionSection();
  virtual void endTableDescriptionSection();
  virtual void startTableDescription(const std::string& tablename);
  virtual void listTableEntry(const std::string& field, 
			      const std::string& type,
			      const std::string& null, 
			      const std::string& key,
			      const std::string& def, 
			      const std::string& extra);
  virtual void endTableDescription();
protected:
  std::string Ec(const std::string& s) const;
  std::string E(const std::string& s) const;
  std::string abbrev(const std::string& s, unsigned count);
public:
  std::ostream& fStream;
  unsigned      fTableMembers;
  std::string   fTableName;
  bool          fTableRestart;

  std::map<std::string,unsigned> fAbbrev;
};

DBLaTeXWriter::~DBLaTeXWriter()
{
  // nothing to see here
}

void DBLaTeXWriter::startDocument(const std::string& database)
{
  fStream 
    << "\\documentclass[letterpaper,10pt]{article}\n"
    << "\\usepackage{newcent}\n"
    << "\n"
    << "\\newlength{\\sjfhmargin}\n"
    << "\\newlength{\\sjfvmargin}\n"
    << "\\setlength{\\sjfhmargin}{0.75in}\n"
    << "\\setlength{\\sjfvmargin}{0.75in}\n"
    << "\n"
    << "\\setlength{\\voffset}{0in}\n"
    << "\\setlength{\\hoffset}{0in}\n"
    << "\n"
    << "\\setlength{\\headheight}{0ex}\n"
    << "\\setlength{\\headsep}{0in}\n"
    << "\\setlength{\\topskip}{0in}\n"
    << "%\\setlength{\\footskip}{0in}\n"
    << "\n"
    << "\\setlength{\\textwidth}{\\paperwidth}\n"
    << "\\addtolength{\\textwidth}{-2\\sjfhmargin}\n"
    << "\\setlength{\\evensidemargin}{-1in}\n"
    << "\\addtolength{\\evensidemargin}{\\sjfhmargin}\n"
    << "\\setlength{\\oddsidemargin}{\\evensidemargin}\n"
    << "\n"
    << "\\setlength{\\textheight}{\\paperheight}\n"
    << "\\addtolength{\\textheight}{-2\\sjfvmargin}\n"
    << "\\addtolength{\\textheight}{-\\headheight}\n"
    << "\\addtolength{\\textheight}{-\\headsep}\n"
    << "\\addtolength{\\textheight}{-\\footskip}\n"
    << "\n"
    << "\\setlength{\\topmargin}{-1in}\n"
    << "\\addtolength{\\topmargin}{\\sjfvmargin}\n"
    << "\n"
    << "\\setlength{\\parindent}{0in}\n"
    << "\\setlength{\\parskip}{1.5ex plus 0.5ex minus 0.2ex}\n"
    << "\n"
    << "\\ifx\\pdfoutput\\undefined \\newcount\\pdfoutput \\fi\n"
    << "\\ifcase\\pdfoutput\n"
    << "  \\special{papersize=11in,8.5in}\n"
    << "\\else\n"
    << "  \\usepackage[pdftex]{hyperref}\n"
    << "  \\pdfpagewidth 8.5truein\n"
    << "  \\pdfpageheight 11truein\n"
    << "  \\pdfhorigin 1truein\n"
    << "  \\pdfvorigin 1truein\n"
    << "  \\hypersetup{\n"
    << "    pdfstartpage=1,\n"
    << "    pdftitle={Description of database: " << Ec(database) << "},\n"
    << "    pdfcreator={\\LaTeX\\ with package \\flqq hyperref\\frqq}\n"
    << "  }\n"
    << "\\fi\n"
    << "\n"
    << "\\begin{document}\n";
    //<< "\\title{Description of database " << Ec(database) << "}\n"
    //<< "\\maketitle\n";
}

void DBLaTeXWriter::endDocument()
{
  fStream 
    << "\\end{document}\n";
}

void DBLaTeXWriter::startTableListingSection()
{
  // nothing to see here
}

void DBLaTeXWriter::endTableListingSection()
{
  // nothing to see here
}

void DBLaTeXWriter::startTableListing()
{
  fStream 
    << "\\addtocounter{table}{1}\n"
    << "\\begin{minipage}{\\textwidth}\n"
    << "\\begin{center}\n"
    << "Table~\\thetable";
  if(fTableMembers)fStream << " -- {\\it continued}";
  fStream 
    << ". "
    << "Listing of tables\n"
    << "\\end{center}\n"
    << "\\begin{center}\n"
    << "\\begin{tabular}{ll}\\hline\n"
    << "NUM & TABLE\\\\\\hline\n";
}

void DBLaTeXWriter::listTable(const std::string& table)
{
  if(fTableMembers)
    {
      if(fTableMembers%50 == 0)
	{
	  unsigned remember = fTableMembers;
	  fTableRestart=true;
	  endTableListing();
	  fStream << "\\addtocounter{table}{-1}\n";
	  fTableMembers = remember;
	  startTableListing();
	  fTableRestart=false;
	}
      else fStream << "\\\\\n";
    }
  fStream 
    << fTableMembers+1 << " & {\\tt " << E(table) << " }";
  fTableMembers++;
}

void DBLaTeXWriter::endTableListing()
{
  fTableMembers=0;
  fStream << "\\\\";
  if(!fTableRestart)fStream << "\\hline\n";
  else fStream << '\n';
  fStream 
    << "\\end{tabular}\n"
    << "\\end{center}\n"
    << "\\end{minipage}\n\n";
}

void DBLaTeXWriter::startTableDescriptionSection()
{
  // nothing to see here
}

void DBLaTeXWriter::endTableDescriptionSection()
{
  if(!fAbbrev.empty())
    {
      fStream 
	<< "\\subsubsection*{Abbreviations:}\n";

      // Re-sort the abbreviation entries by number
      std::map<unsigned, std::string> abbrev;      
      for(std::map<std::string,unsigned>::iterator i = fAbbrev.begin();
	  i!=fAbbrev.end(); i++)abbrev[i->second]=i->first;
      
      for(std::map<unsigned,std::string>::iterator i = abbrev.begin();
	  i!=abbrev.end(); i++)
	{
	  char buffer[20];
	  sprintf(buffer,"%u",i->first);
	  std::string tex;
	  tex+=std::string("$^{");
	  tex+=std::string(buffer);
	  tex+=std::string("}$ {\\tt ");
	  while(!i->second.empty())
	    {
	      tex+=E(i->second.substr(0,80));
	      if(i->second.length()>80)
		{
		  i->second=i->second.substr(80);
		  tex+="\\\\";
		}
	      else i->second.erase();
	    }
	  tex+=std::string("}\n\n");
	  fStream << tex;
	}
    }
}

void DBLaTeXWriter::startTableDescription(const std::string& tablename)
{
  fTableName = tablename;
  fStream 
    << "\\addtocounter{table}{1}\n"
    << "\\begin{minipage}{\\textwidth}\n"
    << "\\begin{center}\n"
    << "Table~\\thetable";
  if(fTableMembers)fStream << " -- {\\it continued}";
  fStream 
    << ". "
    << "Description of table: { \\tt " << Ec(tablename) << "}\n"
    << "\\end{center}\n"
    << "\\begin{center}\n"
    << "\\begin{tabular}{llllll}\\hline\n"
    << "FIELD NAME & TYPE & NULL & KEY & DEF. & EXT.\\\\\\hline\n";
}

void DBLaTeXWriter::listTableEntry(const std::string& field, 
				   const std::string& type,
				   const std::string& null, 
				   const std::string& key,
				   const std::string& def, 
				   const std::string& extra)
{
  
  std::string type_upper;
  for(std::string::const_iterator c=type.begin();c!=type.end();c++)
    type_upper += toupper(*c);

  std::string extra_upper;
  for(std::string::const_iterator c=extra.begin();c!=extra.end();c++)
    extra_upper += toupper(*c);


  if(fTableMembers)
    {
      if(fTableMembers%50 == 0)
	{
	  unsigned remember = fTableMembers;
	  fTableRestart=true;
	  endTableDescription();
	  fStream << "\\addtocounter{table}{-1}\n";
	  fTableMembers = remember;
	  startTableDescription(fTableName);
	  fTableRestart=false;
	}
      else fStream << "\\\\\n";
    }
  fStream 
    << abbrev(field,25) << " & "
    << abbrev(type_upper,20) << " & "
    << abbrev(null,5) << " & "
    << abbrev(key,5) << " & "
    << abbrev(def,10) << " & "
    << abbrev(extra,5);
  fTableMembers++;
}

void DBLaTeXWriter::endTableDescription()
{
  fTableMembers=0;
  fStream << "\\\\";
  if(!fTableRestart)fStream << "\\hline\n";
  else fStream << '\n';
  fStream 
    << "\\end{tabular}\n"
    << "\\end{center}\n"
    << "\\end{minipage}\n\n";
}

std::string DBLaTeXWriter::E(const std::string& s) const
{
  std::string o;
  char verb='\0';
  for(std::string::const_iterator c=s.begin(); c!=s.end(); c++)
    {
      if(*c==verb)
	{ 
	  o+=verb; verb='\0';
	}
      if(verb=='\0')
	{ 
	  if(*c=='+')verb='-'; 
	  else verb='+'; 
	  o+="\\verb";
	  o+=verb; 
	}
      o+=*c;
    }
  if(verb!='\0')o+=verb;
  return o;
}

std::string DBLaTeXWriter::Ec(const std::string& s) const
{
  std::string o;
  for(std::string::const_iterator c=s.begin(); c!=s.end(); c++)
    {
      
      switch(*c)
	{
	case '$':
	case '&':
	case '%':
	case '#':
	case '_':
	case '{':
	case '}':
	  o+='\\';
	  o+=*c;
	  break;
	default:
	  o+=*c;
	  break;
	}
    }
  return o;
}

std::string DBLaTeXWriter::abbrev(const std::string& s, unsigned count)
{
  std::string o;
  o+=std::string("{\\tt ");
  o+=E(s.substr(0,count));
  o+=std::string("}");
  if(s.size()>count)
    {
      std::map<std::string,unsigned>::iterator f = fAbbrev.find(s);
      unsigned mark;
      char buffer[20];
      if(f!=fAbbrev.end())mark=f->second;
      else fAbbrev[s]=mark=fAbbrev.size();
      sprintf(buffer,"%u",mark);
      o+=std::string("...$^{");
      o+=std::string(buffer);
      o+=std::string("}$");
    }
  return o;
}

int main(int argc, char** argv)
{
  // --------------------------------------------------------------------------
  // Process command line
  // --------------------------------------------------------------------------

  VSOptions options(argc,argv);
  VSDBFactory::configure(&options);
  std::string progname(*argv);
  argc--,argv++;
  if(!argc)
    {
      std::cerr << "Usage: " << progname << " database" << std::endl;
      exit(EXIT_FAILURE);
    }
  std::string database(*argv);

  // --------------------------------------------------------------------------
  // Create the database
  // --------------------------------------------------------------------------

  VSDatabase* db = VSDBFactory::getInstance()->createVSDB();
  db->useDatabase(database);

  // --------------------------------------------------------------------------
  // Create the writer
  // --------------------------------------------------------------------------

  DBWriter* writer = new DBLaTeXWriter(std::cout);
  writer->startDocument(database);
  
  // --------------------------------------------------------------------------
  // Query tables
  // --------------------------------------------------------------------------

  std::vector<std::string> tables;
  
  VSDBStatement* stmt = db->createQuery("SHOW TABLES");
  std::string table;
  stmt->bindToResult(table);
  stmt->execute();
  while(stmt->retrieveNextRow())tables.push_back(table);
  delete stmt;
  std::sort(tables.begin(),tables.end());
  
  writer->startTableListingSection();
  writer->startTableListing();
  for(std::vector<std::string>::const_iterator t=tables.begin(); 
      t!=tables.end(); t++)writer->listTable(*t);
  writer->endTableListing();
  writer->endTableListingSection();

  writer->startTableDescriptionSection();
  for(std::vector<std::string>::const_iterator t=tables.begin(); 
      t!=tables.end(); t++)
    {
      stmt = db->createQuery("DESCRIBE "+*t);
      std::string field;
      std::string type;
      std::string null;
      std::string key;
      std::string def;
      std::string extra;
      stmt->bindToResult(field);
      stmt->bindToResult(type);
      stmt->bindToResult(null);
      stmt->bindToResult(key);
      stmt->bindToResult(def);
      stmt->bindToResult(extra);
      stmt->execute();

      writer->startTableDescription(*t);
      while(stmt->retrieveNextRow())
	writer->listTableEntry(field, type, null, key, def, extra);
      writer->endTableDescription();
      
      delete stmt;
    }
  writer->endTableDescriptionSection();

  // --------------------------------------------------------------------------
  // Clean up
  // --------------------------------------------------------------------------
  
  writer->endDocument();
  delete writer;
}
