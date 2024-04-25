// tsv2insert.cpp

/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE                          
*               National Center for Biotechnology Information
*                                                                          
*  This software/database is a "United States Government Work" under the   
*  terms of the United States Copyright Act.  It was written as part of    
*  the author's official duties as a United States Government employee and 
*  thus cannot be copyrighted.  This software/database is freely available 
*  to the public for use. The National Library of Medicine and the U.S.    
*  Government have not placed any restriction on its use or reproduction.  
*                                                                          
*  Although all reasonable efforts have been taken to ensure the accuracy  
*  and reliability of the software and data, the NLM and the U.S.          
*  Government do not and cannot warrant the performance or results that    
*  may be obtained by using this software or data. The NLM and the U.S.    
*  Government disclaim all warranties, express or implied, including       
*  warranties of performance, merchantability or fitness for any particular
*  purpose.                                                                
*                                                                          
*  Please cite the author in any work or product based on this material.   
*
* ===========================================================================
*
* Author: Vyacheslav Brover
*
* File Description:
*   Group columns of a table
*
*/

#undef NDEBUG

#include "../common.hpp"
#include "tsv.hpp"
using namespace Common_sp;
#include "../version.inc"

#include "../common.inc"



namespace
{
  
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Convert a tsv-table into SQL insert statements")
  	{
      version = VERSION;
  	  addPositional ("tsv_table", ".tsv-table");
  	  addKey ("sql_table", "SQL table name (default = <tsv_table>");
  	  addKey ("file_name_col", "Column to store the file name");  	  
  	  addFlag ("add_go", "Add 'go' at the end of the SQL statements");
  	}
  	
  	
 
	void body () const final
	{
		const string pathName    = getArg ("tsv_table");
		      string sqlTable    = getArg ("sql_table");
		const string fileNameCol = getArg ("file_name_col");
		const bool add_go        = getFlag ("add_go");


    TextTable tt (pathName);
    tt. qc ();
    
    const string fName (trimExtension (getFileName (tt. name)));
    
    if (! fileNameCol. empty ())
    {
      TextTable::Header h (fileNameCol);
      h. numeric = false;
      tt. header << std::move (h);
      for (StringVector& row : tt. rows)            
        row << fName;
      tt. qc ();
    }
    
    if (sqlTable. empty ())
      sqlTable = fName;

    string s ("insert into " + sqlTable + " (");
    bool first = true;    
    for (const TextTable::Header& h : tt. header)
    {
      if (! first)
        s += ", ";
      s += h. name;
      first = false;
    }
    s += ") values (";
    for (const StringVector& row : tt. rows)            
    {
      cout << s;
      first = true;    
      FFOR (size_t, i, row. size ())
      {
        if (! first)
          cout << ", ";
        const string& f = row [i];
        if (f. empty ())
          cout << "null";
        else
        {
          const bool numeric = tt. header [i]. numeric;
          if (! numeric)
            cout << "'";
          cout << f;  // ' -> '' ??
          if (! numeric)
            cout << "'";
        }
        first = false;
      }
      cout << ")" << endl;      
    	if (add_go)
    	  cout << "go" << endl;
    }
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}

