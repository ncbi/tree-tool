// tsv_shift.cpp

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
*   Shift a tsv-table 
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
    : Application ("<row_i> <nl> <row_i+1> <row_i+2> <nl> ... -> <row_i> <tab> <row_i+1> <nl> <row_i+1> <tab> <row_i+2> <nl> ...")
  	{
      version = VERSION;
  	  addPositional ("table", "tsv-table");
  	}
  	
  	
 
	void body () const final
	{
		const string fName  = getArg ("table");


    const TextTable tIn (fName);
    tIn. qc ();
    if (tIn. rows. empty ())
      throw runtime_error ("Table " + strQuote (fName) + " is empty");

    TextTable tOut;
    FOR_START (int, i, 1, 3)
      for (TextTable::Header h : tIn. header)
      {
        h. name += "_" + to_string (i);
        tOut. header << std::move (h);
      }
    FFOR (TextTable::RowNum, i, tIn. rows. size () - 1)
    {
      StringVector row (tIn. rows [i]);
      row << tIn. rows [i + 1];
      tOut. rows << std::move (row);
    }  
    tOut. qc ();
    tOut. saveText (cout);
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



