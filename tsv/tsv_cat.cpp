// tsv_cat.cpp

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
*   Concatenate tsv-tables
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
    : Application ("Concatenate tsv-tables. File names are added as the new last column")
  	{
      version = VERSION;
  	  addPositional ("list", "List of tsv-tables with the same headers");
  	  addPositional ("col_name", "Column name of the added file names; may be empty");
  	}
  	
  	
 
	void body () const final
	{
		const string listFName = getArg ("list");
		const string colName   = getArg ("col_name");
		
		
		TextTable total;
		total. pound = true;
		{
  		LineInput li (listFName, 1);  // PAR
  		while (li. nextLine ())
  		{
        TextTable tab (li. line);
        tab. qc ();
        
        if (! colName. empty ())
    		  if (tab. hasColumn (colName))
    		    throw runtime_error ("Files already have column " + strQuote (colName));

        // header
        QC_ASSERT (! tab. header. empty ());
  		  if (total. header. empty ())
  		    total. header = tab. header;
  		  else
  		  {
  		    if (total. header. size () != tab. header. size ())
  		      throw runtime_error ("Different number of columns in files");
  		    FFOR (size_t, i, total. header. size ())
  		      if (total. header [i]. name != tab. header [i]. name)
  		        throw runtime_error ("Different column #" + to_string (i + 1) + " in files");
  		  }
  		  
        if (! colName. empty ())
        {
    		  const string fName = getFileName (li. line);
    		  for (StringVector& row : tab. rows)
    		    row << fName;
    		}
  		  
  		  total. rows << tab. rows;		  
  		}
    }
    if (! colName. empty ())
  		total. header << TextTable::Header (colName);
		total. qc ();


    total. saveText (cout);
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



