// mergePairs.cpp

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
*   Merge pairs with identical first element
*
*/


#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;



namespace 
{
  
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Merge pairs with identical first element:\nA1\tB1\nA1\tB2\nA2\tB3\nA4\tB4\n->\nA1\tB1\tB2\nA2\tB3\tB4")
  	{
  	  addPositional ("in", "Input file");
  	}



	void body () const final
	{
		const string fName = getArg ("in");

    LineInput f (fName, 100 * 1024, 1000);
    string lhs, lhs_prev;
  	while (f. nextLine ())
  	{ 
  	  trim (f. line);
  	  if (f. line. empty ())
  	    continue;
  	  lhs = findSplit (f. line, '\t');
  	  QC_ASSERT (! lhs. empty ());
  	  const string& rhs = f. line;
  	  QC_ASSERT (! rhs. empty ());
  	  if (lhs != lhs_prev)
  	  {
  	    if (! lhs_prev. empty ())
  	      cout << endl;
  	    cout << lhs;  	    
  	  }
  	  cout << '\t' << rhs;
  	  lhs_prev = move (lhs);
  	}
    if (! lhs_prev. empty ())
    	cout << endl;
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}
