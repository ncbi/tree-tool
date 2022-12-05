// replace_dict.cpp

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
*   Replace words using a dictionary
*
*/

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;
#include "version.inc"



namespace
{
  
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Replace words in a text using a dictionary")
  	{
      version = VERSION;
  	  addPositional ("file", "Text file");
  	  addPositional ("dictionary", "Pairs of words: <word1> <word2>");
  	  addFlag ("reverse", "Replace <word2> by <word1>, otherwise <word1> by <word2>");
  	}
  	
  	
 
	void body () const final
	{
		const string fName = getArg ("file");
		const string dict  = getArg ("dictionary");
		const bool   rev   = getFlag ("reverse");

    
    Vector<Pair<string>> pairs;
    {
      LineInput f (dict);
      Istringstream iss;
      while (f. nextLine ())
      {
        iss. reset (f. line);
        string a, b;
        iss >> a >> b;
        QC_ASSERT (! b. empty ());
        Pair<string> p (move (a), move (b));
        if (rev)
          p. swap ();
        pairs << move (p);;
      }
    }

    LineInput f (fName);
    while (f. nextLine ())
    {
      for (const Pair<string>& p : pairs)
        replaceStr (f. line, p. first, p. second);
      cout << f. line << endl;
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



