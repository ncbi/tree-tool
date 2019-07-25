// setMinus.cpp

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
*   Print all combinations of pairs out of a list of objects
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
    : Application ("Print the pairs of words from a list. In each pair: word1 < word2")
  	{
  	  addPositional ("list", "File with different words");
  	}



	void body () const final
	{
		const string listFName = getArg ("list");


    StringVector words;  
    {
      LineInput f (listFName);
      words = f. getVector ();
    }
    words. sort ();
    ASSERT (words. isUniq ());
    
    {
    	Progress prog (words. size ());
	    for (const string& w1 : words)
	    {
	    	prog (w1);
	    	for (const string& w2 : words)
	    		if (& w1 == & w2)
	    			break;
	    	  else
	    	  {
	    	  	ASSERT (w1 != w2);
	    	  	const string* p1 = & w1;
	    	  	const string* p2 = & w2;
	    	  	if (*p1 > *p2)
	    	  		swap (p1, p2);
	    	  	cout << *p1 << '\t' << *p2 << endl;
	    	  }
	    }
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
