// extractPairs.cpp

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
*   Extract pairs out of a specifically formatted file
*
*/


#undef NDEBUG

#include "common.hpp"
using namespace Common_sp;
#include "version.inc"

#include "common.inc"



namespace 
{
  
  
List<string> getList (string &s,
                      const string &delim)
{
  trim (s);
  replaceStr (s, delim, " ");
  replaceStr (s, "  ", " ");
  return str2list (s);  
}

  
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("A1<delim>A2<delim>A3...\\tB1<delim>B2<delim>B3... -> A1\\tB1\\nA1\\tB2\\nA1\\tB3\\nA2\\tB1\\nA2\\tB2\\nA2\\tB3\\nA3\\tB1\\nA3\\tB2\\nA3\\tB3...")
  	{
      version = VERSION;
  	  addPositional ("in", "Input file");
  	  addKey ("delim", "Delimiter", " ");
  	}



	void body () const final
	{
		const string fName = getArg ("in");
		const string delim = getArg ("delim");
		ASSERT (! delim. empty ());

    LineInput f (fName, 1000);  // PAR
  	while (f. nextLine ())
  	{ 
  	  trim (f. line);
  	  if (f. line. empty ())
  	    continue;
  	  string lhs (findSplit (f. line, '\t'));
  	  const List<string> lhsVec (getList (lhs, delim));
  	  const List<string> rhsVec (getList (f. line, delim));
  	  for (const string& a : lhsVec)
    	  for (const string& b : rhsVec)
    	    cout << a << '\t' << b << endl;
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
