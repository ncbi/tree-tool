// ramdom_words.cpp

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
*   Generate distinct random words
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
    : Application ("Generate distinct ASCII-printable random words")
  	{
      version = VERSION;
  	  addPositional ("num", "Number of random words");
  	  addPositional ("len_max", "Max. word length");
  	  addKey ("exclude", "Characters to exclude");
  	  addKey ("exclude_first", "Characters in the beginning of a word to exclude");
  	}
  	
  	
 
	void body () const final
	{
		const size_t num      = str2<size_t> (getArg ("num"));
		const size_t len_max  = str2<size_t> (getArg ("len_max"));
		const string exclude  =               getArg ("exclude");
		const string exclude1 =               getArg ("exclude_first");
		

    Rand rand (seed_global);
    
    string alphabet;
    for (char c = ' ' + 1; c < 127; c++)
      if (! contains (exclude, c))
        alphabet += c;    
      
    string alphabet1;
    for (const char c : alphabet)
      if (! contains (exclude1, c))
        alphabet1 += c;    
        
    Set<string> words;
    while (words. size () < num)
    {
      string s;
      FFOR (size_t, i, rand. get (len_max) + 1)
        s += i ? alphabet  [rand. get (alphabet.  size ())]
               : alphabet1 [rand. get (alphabet1. size ())];
      words << s;
    }
    
    for (const string& s : words)
      cout << s << endl;
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



