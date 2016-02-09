//-*-mode:c++; mode:font-lock;-*-

/*! \file VSDatabase.cpp
  Base class for database input output

  \author     Stephen Fegan               \n
              UCLA                        \n
              sfegan@astro.ucla.edu       \n        

  \version    0.1
  \date       03/02/2005
  \note
*/

#include <stack>

#include "VSDatabase.hpp"

// ****************************************************************************
// ****************************************************************************
// **                                                                        **
// ** VSDatabase                                                             **
// **                                                                        **
// ****************************************************************************
// ****************************************************************************

/*! \class VERITAS::VSDatabase   

  VSDatabase is provided the API to a simple database system. The API
  is very limited, supporting only a small fraction of the
  functionality available with modern relational database management
  systems (RDBMS). Ideally it is completely independent of any
  specific RDBMS but it was designed with the MySQL 4.1 prepared
  statement (MySQL-PS) interface in mind. That is to say that the
  functionality of this class can be implemented with the MySQL-PS API
  with a minimum of effort. It can also be implemented relatively
  easily with a atandard text-based database API such as MySQL 3.XX
  and 4.0X, and probably with the text based interface to any SQL
  database. Other PS interfaces, such as PostgreSQL probably will not
  map as naturally to this API.

  This class is not meant as a fully featured, robust, mission
  critical interface to a database. It is missing the functionality of
  transactions and rollback. It supports only input and output of the
  usual numeric types, strings (C-type and STL) and a number of date and
  time types. Other types may work if the native RDMS-API can convert
  them to and from one of the supported types internally. For example,
  MySQL will convert almost anything to a string which can then be
  used with this class.

  In standard usage, the functions createDatabase() and useDatabase()
  are called to create and connect to a database. Tables can be
  created in the database using createTable(). Tables can be deleted
  using dropDatabase() and data can be removed from tables using
  deleteFromTable().
  
  Input and ouput on tables is done by calling createQuery() which
  parses the query and returns a "prepared statement" which can be
  executed, and re-executed as desired. The query should be standard
  SQL with optional parameter markers. These markers, which are '?' 
  characters embedded in the SQL, are bound to variables in the user's
  code and substituted with real data when the query is executed. The
  results of the query are returned in a row-by-row manner through
  variables which have been previously bound to the results. Binding
  of variables to the data is done using the
  VSDBStatement::bindToParam() and VSDBStatement::bindToResult()
  functions.

  \par Example: Inserting data to a table using bound parameters.
  \code
  VSDatabase *db = VSDBFactory::instantiateDB();
  db->createDatabase("example");
  db->useDatabase("example");
  db->dropTable("test");
  db->createTable("test","channel INTEGER UNSIGNED, mean DOUBLE, dev DOUBLE");
  VSDBStatement *stmt = db->createQuery("INSERT INTO test VALUES (?,?,?)");
  unsigned c;
  double mean;
  double dev;
  stmt->bindToParam(c);
  stmt->bindToParam(mean);
  stmt->bindToParam(dev);
  for(c=0;c<SOME_CAMERA_INSTANCE->numChannels();c++)
    {
      SOME_PEDESTAL_CALCULATOR_INSTANCE->getPed(c,mean,dev);
      stmt->execute();
    }
  delete stmt;
  delete db;
  \endcode

  \par Example: Extracting data from a table using bound results.
  \code
  VSDatabase *db = VSDBFactory::instantiateDB();
  db->useDatabase("example");
  VSDBStatement *stmt = db->createQuery("SELECT * FROM test");
  unsigned c;
  double mean;
  double dev;
  stmt->bindToResult(c);
  stmt->bindToResult(mean);
  stmt->bindToResult(dev);
  stmt->execute();
  while(stmt->retreiveNextRow())
    {
      std::cout << std::setw(3) << c <<
                << std::setw(10) << mean <<
                << std::setw(10) << dev << std::endl;
    }
  delete stmt;
  delete db;
  \endcode

*/

VERITAS::VSDatabase::~VSDatabase()
{
  // nothing to see here
}

/*! \fn VERITAS::VSDatabase::getErrorMessage() throw()

  \return Message which describes what went wrong with last query
*/

/*! \fn VERITAS::VSDatabase::createDatabase(const std::string& dbname) throw()
  
  \param dbname The name of the database to create, if it does not
  exist already.
  
  \return Negative value on failure. Zero or positive on success.
*/


/*! \fn VERITAS::VSDatabase::dropDatabase(const std::string& dbname) throw()

  \param dbname The name of the database to drop, if it exists.

  \return Negative value on failure. Zero or positive on success.
*/

/*! \fn VERITAS::VSDatabase::useDatabase(const std::string& dbname) throw()

  \param dbname The name of the database to select.

  \return Negative value on failure. Zero or positive on success.
*/

/*! \fn VERITAS::VSDatabase::createTable(const std::string& table, const std::string& spec) throw()

  \param table The name of the table to create, if it does not
  exist already.

  \param spec SQL specification of the columns in the table.

  \return Negative value on failure. Zero or positive on success.
*/

/*! \fn VERITAS::VSDatabase::dropTable(const std::string& table) throw()
  \param table The name of the table to drop, if it exists.
  \return Negative value on failure. Zero or positive on success.
*/

/*! \fn VERITAS::VSDatabase::deleteFromTable(const std::string& table, const std::string& cond) throw()

  \param table The name of the database to select.

  \param cond SQL condition defining what to delete. If it is empty then 
  all rows are deleted.

  \return Number of rows deleted, negative value on failure. 
*/

/*! \fn VERITAS::VSDatabase::createQuery(const std::string& query, uint32_t flags = 0) throw()
  
  \param query The SQL for the query with optional parameter markers.

  \param flags Flags which modify the behavior of the query. Multiple
  flags can be or'd together. See, for example, FLAG_NO_BUFFER.

  \return A prepared query which can be exeuted multiple times with
  different paramaters. The statement is created on the heap and
  should should be deleted when it is no longer being used.
*/

/*! \fn VERITAS::VSDatabase::createSelectQuery(const std::string& table, const std::string& condition, const std::string& columns, uint32_t flags) throw()

  Utility function to create a prepared statement from a SELECT query

  \param table Name of the table to within which to execute the SELECT query.

  \param condition Condition to select the data. The condition is used to
  generate the WHERE clause to the SELECT. If the condition is is empty
  then no WHERE clause is included.

  \param columns List of column names, or "*".

  \param Flags which modify the behavior of the query. Multiple
  flags can be or'd together. See, for example, FLAG_NO_BUFFER.

  \return A statement which can be used to execute the SELECT query
*/

/*! \fn VERITAS::VSDatabase::createInsertQuery(const std::string& table, unsigned num_columns, const std::string& columns, uint32_t flags) throw()

  Utility function to create a prepared statement from an INSERT query

  \param table Name of the table to within which to execute the INSERT query.

  \param num_columns Number of columns to INSERT data into, which 
  corresponds to the number of parameters which must be bound to the 
  query.

  \param columns List of column names to insert data into or an empty 
  string, meaning insert into all columns.

  \param Flags which modify the behavior of the query. Multiple
  flags can be or'd together. See, for example, FLAG_NO_BUFFER.

  \return A statement which can be used to execute the INSERT query
*/
  
VERITAS::VSDBStatement* 
VERITAS::VSDatabase::createSelectQuery(const std::string& table, 
				       const std::string& condition,
				       const std::string& columns, 
				       uint32_t flags) throw()
{
  std::ostringstream qstream;
  qstream << "SELECT ";
  if(!columns.empty())qstream << columns;
  else qstream << "*";
  qstream << " FROM " << table;
  if(!condition.empty())qstream << " WHERE " << condition;
  return createQuery(qstream.str(), flags);
}

#warning Finish comments

/*!
  This convenience makes it easier to assemble the SQL which specifies
  the columns of a table. The output is intended to be used in createTable()
  to create a table.
  
  The function simply assembles the components to one of the following
  - <code>"NAME TYPE"</code> if <i>first</i> is true and <i>post</i> is empty
  - <code>"NAME TYPE POST"</code> if <i>first</i> is true
  - <code>",NAME TYPE"</code> if <i>first</i> is false and <i>post</i> is empty
  - <code>",NAME TYPE POST"</code> if <i>first</i> is false
 
  \param name Name of the column
  \param type SQL for type of column
  \param first Set to true if this is the first item in the specification 
  list. Otherwise set it to false.
  \param post Any SQL to append after the specification. For example
  "NOT NULL"
  \return SQL of this column
*/
std::string VERITAS::VSDatabase::
sqlSpec(const std::string& name, const std::string& type,
	bool first, const std::string& post) const throw()
{
  std::string o;
  if(!first)o=",";
  o+=name;
  o+=' ';
  o+=type;
  if(!post.empty())
    {
      o+=' ';
      o+=post;
    }
  return o;
}

/*! \fn VERITAS::VSDatabase::sqlSpecOf(const std::string& name, const T& x, bool first, const std::string& post="") const throw()

  This function is essentially just a call to sqlSpec() with the
  <i>type</i> parameter set by the appropriate typeOf() function.  The
  SQL specification required to construct a column with the type of
  <i>x</i> is returned in a database independent manner. This is
  achieved through the typeOf() functions, which are virtual in
  VSDatabase and implemented in each derived database class. Each such
  class knows the SQL for each of the C++ types for which a typeOf()
  function exists. Since only MySQL is implemeted at this time, there
  is nothing to be gained from this approach at present.

  The usefulness of this function may not be obvious. From the example
  below one could argue that it would be much easier to just write the
  SQL specification for the table by hand. Indeed, it ise easier to do
  so. However, this function allows the SQL which constructs a table
  to be tied more closely to the specification of the C++ type that
  will interface with that table. If it is decided to change the type
  of a variable in the struct (for example from <code>int</code> to
  <code>float</code>) then the code constructing the table will change
  at the same time. In the example you can see that a temporary
  variable may be needed to use this function. Ideally this ideal
  would be taken further, the specification of the table and all I/O
  code would be generated directly from the C++ code using a code
  generator.

  \param name Name of the column
  \param x Value is unused. The parameter is only for the purpose of 
  overloading.
  \param first Set to true if this is the first item in the specification 
  list. Otherwise set it to false.
  \param post Any SQL to append after the specification. For example
  "NOT NULL"
  \return SQL of this column
  
  \section Example
  \code
  struct Pedestal
  {
    std::string type;
    VSDBDataTime time;
    unsigned channel;
    double mean;
    double rms;
  };
  Pedestal temp;
  VSDatabase* db = VSDBFactory::instantiateDB();
  db->useDatabase("Analysis");
  db->dropTable("pedestals");
  db->createTable("pedestals",
                  db->sqlSpecOf("type",temp.type,true,"NOT NULL")+
                  db->sqlSpecOf("time",temp.time,false,"NOT NULL")+
                  db->sqlSpecOf("channel",temp.channel,false,"NOT NULL")+
                  db->sqlSpecOf("mean",temp.mean,false,"NOT NULL")+
                  db->sqlSpecOf("rms",temp.rms,false,"NOT NULL")+
		  std::string(" PRIMARY KEY (type,time,channel)");
  delete db;
  \endcode 
*/
  

// ****************************************************************************
// ****************************************************************************
// **                                                                        **
// ** VSDBStatement                                                          **
// **                                                                        **
// ****************************************************************************
// ****************************************************************************

/*! \class VERITAS::VSDBStatement

  VSDBStatement encapsulates the concept of a prepared SQL
  statement. A statement is prepared through the VSDatabase class
  using the VSDatabase::createQuery function. This class then allows
  user-space variables to be bound to the parameters and to results of
  the query. When this is done, the statement can be executed, either
  one time or repeatedly, by setting the bound input variables as
  desired and calling execute(). It the query is succesful and returns
  data, the results can be transferred to the bound output variables
  by calling retrieveNextRow().

  If the API to the database manager supports the concept of prepared
  statements natively then concrete implementations of this class can 
  naturally take advantage of this ability. If they do not then this
  functionality must be simulated.

*/

VERITAS::VSDBStatement::~VSDBStatement()
{
  // nothing to see here
}

// ****************************************************************************
// ****************************************************************************
// **                                                                        **
// ** parseQueryForParameters                                                **
// **                                                                        **
// ****************************************************************************
// ****************************************************************************

enum SQLPQFP_State { PS_ACTIVE, PS_STRING, PS_ESCAPE };

/*! 
  
  A very simple parser of SQL which searches for parameter markers
  ('?')  and splits the query into pieces around the markers. The
  parser correctly handles '?' characters which are in strings.

  For example, the query 

  <code>"INSERT INTO test VALUES (?,'How are you ?',?)"</code>

  returns the following elements in <i>query_bits</i>

  -# <code>"INSERT INTO test VALUES ("</code>
  -# <code>",'How are you ?',"</code>
  -# <code>",)"</code>
  
  The number of pieces of in the vector on return always equals the
  number of parameter markers plus one. Therefore, if the first (last)
  character in the query is a '?' then the first (last) element in the
  vector is a null string.

  \param query SQL query to parse.
  \param query_bits Vector of STL strings in which to return the
  components of the query which lie between the parameter markers
  
  \return <i>true</i> on success, <i>false</i> on failure.
  
*/

bool VERITAS::parseSQLForParameters(const std::string& query, 
				    std::vector<std::string>& query_bits)
{
  query_bits.clear();
  std::string bit;
  std::stack<SQLPQFP_State> states;
  states.push(PS_ACTIVE);
  std::string::const_iterator c=query.begin();
  while(c!=query.end())
    {
      switch(states.top())
	{
	case PS_ACTIVE:
	  switch(*c)
	    {
	    case '\'':
	      bit += *c;
	      c++;
	      states.push(PS_STRING);
	      break;
	    case '\\':
	      bit += *c;
	      c++;
	      states.push(PS_ESCAPE);
	      break;
	    case '?':
	      c++;
	      query_bits.push_back(bit);
	      bit.erase();
	      break;
	    default:
	      bit += *c;
	      c++;
	      break;	      
	    }
	  break;
	case PS_STRING:
	  switch(*c)
	    {
	    case '\'':
	      bit += *c;
	      c++;
	      states.pop();
	      break;
	    case '\\':
	      bit += *c;
	      c++;
	      states.push(PS_ESCAPE);
	      break;
	    default:
	      bit += *c;
	      c++;
	      break;	      
	    }
	  break;
	case PS_ESCAPE:
	  states.pop();
	  bit += *c;
	  c++;
	  break;
	}
    }
  query_bits.push_back(bit);
  return states.size()==1;
}

std::ostream& operator<<(std::ostream& stream, const VERITAS::VSDBDate& x)
{
  std::ostringstream ts;
  ts << std::setfill('0') 
     << std::setw(4) << x.year << '-'
     << std::setw(2) << x.month << '-'
     << std::setw(2) << x.day;
  stream << ts.str();
  return stream;
}

std::ostream& operator<<(std::ostream& stream, const VERITAS::VSDBTime& x)
{
  std::ostringstream ts;
  ts << std::setfill('0') 
     << std::setw(2) << x.hour << ':'
     << std::setw(2) << x.minute << ':'
     << std::setw(2) << x.second;
  stream << ts.str();
  return stream;
}

std::ostream& operator<<(std::ostream& stream, const VERITAS::VSDBDateTime& x)
{
  std::ostringstream ts;
  ts << std::setfill('0') 
     << std::setw(4) << x.year << '-'
     << std::setw(2) << x.month << '-'
     << std::setw(2) << x.day << ' '
     << std::setw(2) << x.hour << ':'
     << std::setw(2) << x.minute << ':'
     << std::setw(2) << x.second;
  stream << ts.str();  
  return stream;
}

std::ostream& operator<<(std::ostream& stream, const VERITAS::VSDBTimestamp& x)
{
  std::ostringstream ts;
  ts << std::setfill('0') 
     << std::setw(4) << x.year << '-'
     << std::setw(2) << x.month << '-'
     << std::setw(2) << x.day << ' '
     << std::setw(2) << x.hour << ':'
     << std::setw(2) << x.minute << ':'
     << std::setw(2) << x.second;
  stream << ts.str();  
  return stream;
}
