// tsv_group.cpp

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
    : Application ("Group columns of a tsv-table")
  	{
      version = VERSION;
  	  addPositional ("table", "tsv-table file with a header");
  	  addKey ("by", "Comma-separated list of columns to group by"); 
  	  addKey ("count", "name of the added \"count\" column");
  	  addKey ("sum", "Comma-separated list of columns to sum"); 
  	  addKey ("aggr", "Comma-separated list of columns to aggregate: make unique and sort"); 
  	}
  	
  	
 
	void body () const final
	{
		const string fName  = getArg ("table");
		const string byS    = getArg ("by");
		const string countS = getArg ("count");
		const string sumS   = getArg ("sum");
		const string aggrS  = getArg ("aggr");


    TextTable tt (fName);
    tt. qc ();
    if (verbose ())
      tt. printHeader (cout);      
    
    const StringVector by   (byS,   ',', true);
          StringVector sum  (sumS,  ',', true);    
    const StringVector aggr (aggrS, ',', true);    


    // QC
    for (const string& s : by)
      if (! tt. hasColumn (s))
        throw runtime_error ("Table has no by-column " + strQuote (s));
        
    if (! countS. empty ())
    { 
      if ( tt. hasColumn (countS))
        throw runtime_error ("Table already has the column " + strQuote (countS));
      tt. header << TextTable::Header (countS);
      for (StringVector& row : tt. rows)
        row << "1";
      tt. qc ();
      sum << countS;
    }
    
    for (const string& s : sum)
    {
      if (! tt. hasColumn (s))
        throw runtime_error ("Table has no sum-column " + strQuote (s));
      if (by. contains (s))
        throw runtime_error ("Same by-column and sum-column: " + strQuote (s));
      if (! tt. header [tt. col2num (s)]. numeric)
        throw runtime_error ("Column " + strQuote (s) + " is not numeric");
    }
        
    for (const string& s : aggr)
    {
      if (! tt. hasColumn (s))
        throw runtime_error ("Table has no aggregation column " + strQuote (s));
      if (by. contains (s))
        throw runtime_error ("Same by-column and aggregation column: " + strQuote (s));
    }
        

    tt. group (by, sum, aggr);
    tt. qc ();
    
    tt. saveText (cout);
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



