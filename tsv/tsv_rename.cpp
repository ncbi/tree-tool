// tsv_rename.cpp

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
*   Rename a column in a tsv-table
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
    : Application ("Print a tsv-table with a renamed column")
  	{
      version = VERSION;
      addPositional ("tsv", "tsv-file name");
  	  addPositional ("col_num", "1-based column number");
  	  addPositional ("col_name", "New column name");
  	}
  	
  	
 
	void body () const final
	{
		const string fName   = getArg ("tsv");
		      size_t colNum  = (size_t) arg2uint ("col_num");
		      string colName = getArg ("col_name");
				
		
		TextTable tt (fName);
		tt. qc ();
		
		QC_ASSERT (colNum > 0);
		colNum--;
		if (colNum >= tt. header. size ())
		  throw runtime_error ("Column number is too large");

    trim (colName);		
		QC_ASSERT (! colName. empty ());
		const size_t oldNum = tt. col2num_ (colName);
		if (oldNum != no_index && oldNum != colNum)
		  throw runtime_error ("Column " + strQuote (colName) + " already exists");
    tt. header [colNum]. name = colName;
    
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



