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
    : Application ("Group columns of a tsv-table")
  	{
      version = VERSION;
  	  addPositional ("table", "tsv-table");
  	  addKey ("by", "Comma-separated list of columns to group by"); 
  	  addKey ("count", "name of the added \"count\" column");
  	  addKey ("sum", "Comma-separated list of columns to sum"); 
  	  addKey ("min", "Comma-separated list of columns to find minimum"); 
  	  addKey ("max", "Comma-separated list of columns to find maximum"); 
  	  addKey ("aggr", "Comma-separated list of columns to aggregate: make unique and sort"); 
  	}
  	
  	
 
	void body () const final
	{
		const string fName  = getArg ("table");
		const string byS    = getArg ("by");
		const string countS = getArg ("count");
		const string sumS   = getArg ("sum");
		const string minS   = getArg ("min");
		const string maxS   = getArg ("max");
		const string aggrS  = getArg ("aggr");


    TextTable tt (fName);
    tt. qc ();
    if (verbose ())
      tt. printHeader (cout);      
    
    const StringVector by   (byS,   ',', true);
          StringVector sum  (sumS,  ',', true);    
          StringVector minV (minS, ',', true);    
          StringVector maxV (maxS, ',', true);    
    const StringVector aggr (aggrS, ',', true);    

    for (string& s : sum)
      tt. substitueColumn (s, s + "_sum");
    for (string& s : minV)
      tt. substitueColumn (s, s + "_min");
    for (string& s : maxV)
      tt. substitueColumn (s, s + "_max");

    if (! countS. empty ())
    { 
      if ( tt. hasColumn (countS))
        throw runtime_error ("Table already has column " + strQuote (countS));
      tt. header << TextTable::Header (countS);  // numeric
      for (StringVector& row : tt. rows)
        row << "1";
      tt. qc ();
      sum << countS;
    }
    
    tt. group (by, sum, minV, maxV, aggr);
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



