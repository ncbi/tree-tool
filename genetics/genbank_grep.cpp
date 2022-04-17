// asn_grep.cpp

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
*   Print accession numbers of the sequence GenBank records containing a keyword
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "genbank_text.hpp"



namespace 
{


struct ThisApplication : Application
{
	ThisApplication ()
	  : Application ("Print accession numbers of the sequence GenBank records containing a keyword")
{
  addPositional ("in", "Genbank file");
  addPositional ("keyword", "Keyword to find (case-insensitive)");
  addKey ("field", "Field name");
}

  
  
	void body () const final
  {
	  const string inFName = getArg ("in");
	  const string keyword = getArg ("keyword");
	  const string field   = getArg ("field");
	        
	  
	  LineInput f (inFName, 100000);  // PAR
	  Progress prog;
	  while (! f. eof)
	  {
  	  const Asn_sp::GenbankText gt (f);
  	  if (verbose ())
  	    gt. saveText (cout);  
  	    
  	  const string acc_ver (gt. name2value ("VERSION"));
  	  
	    prog (acc_ver);
	    
	    string content;
  		const string recordName (gt. keyword2name (keyword, true, content));  // PAR
  		if (   ! recordName. empty () 
  		    && (field. empty () || field == recordName)
  		   )
  		  cout << acc_ver << '\t' << content << endl;
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



