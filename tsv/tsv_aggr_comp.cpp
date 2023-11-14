// tsv_aggr_comp.cpp

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
*   Compare two aggregated columns of a tsv-table
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
    : Application ("Compare two aggregated columns (values are comma-separated) of a tsv-table")
  	{
      version = VERSION;
      addPositional ("tsv", "tsv-file name");
  	  addPositional ("col1", "First aggregated column");
  	  addPositional ("col2", "Second aggregated column");
  	  addPositional ("diff1", "New column containing col1 \\ col2");
  	  addPositional ("diff2", "New column containing col2 \\ col1");
  	}
  	
  	
 
	void body () const final
	{
		const string fName = getArg ("tsv");
		const string col1  = getArg ("col1");
		const string col2  = getArg ("col2");
		const string diff1 = getArg ("diff1");
		const string diff2 = getArg ("diff2");

    {		
  		const StringVector vec {{col1, col2, diff1, diff2}};
  		if (! vec. checkUniq ())
  		  throw runtime_error ("Parameter columns must be different");
    }
    				
		
		TextTable tt (fName);
		tt. qc ();
		
		const size_t c1 = tt. col2num (col1);
		const size_t c2 = tt. col2num (col2);
		if (tt. header [c1]. numeric)
		  throw runtime_error ("Column " + strQuote (col1) + " is numeric");
		if (tt. header [c2]. numeric)
		  throw runtime_error ("Column " + strQuote (col2) + " is numeric");
		if (tt. hasColumn (diff1))
		  throw runtime_error ("Column " + strQuote (diff1) + " already exists");
		if (tt. hasColumn (diff2))
		  throw runtime_error ("Column " + strQuote (diff2) + " already exists");
		  
		TextTable::Header h1 (diff1);
		h1. numeric = false;
		TextTable::Header h2 (diff2);
		h2. numeric = false;
    tt. header << std::move (h1) << std::move (h2);
		
		const string sep (1, TextTable::aggr_sep);
    for (StringVector& row : tt. rows)
    {
      const StringVector val1 (TextTable::aggr2values (row [c1]));
      const StringVector val2 (TextTable::aggr2values (row [c2]));
      StringVector res1 (val1);
      res1. setMinus (val2);
      StringVector res2 (val2);
      res2. setMinus (val1);
      row << res1. toString (sep);
      row << res2. toString (sep);
    }
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



