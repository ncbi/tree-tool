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
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../version.inc"



namespace
{
  
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Convert a tsv-table into SQL insert statements")
  	{
      version = VERSION;
  	  addPositional ("table", "tsv-table file with a header");
  	}
  	
  	
 
	void body () const final
	{
		const string fName  = getArg ("table");


    TextTable tt (fName);
    tt. qc ();

    string s ("insert into " + trimExtension (getFileName (tt. name)) + " (");
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
      cout << ");" << endl;
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



