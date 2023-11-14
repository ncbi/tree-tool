// file2hash.cpp

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
*   Print file hash class
*
*/


#undef NDEBUG

#include "common.hpp"
using namespace Common_sp;
#include "version.inc"

#include "common.inc"



namespace 
{
  
  
void process (const string &name,
              bool append)
{
  cout << str2hash_class (name);
  if (append)
    cout << '\t' << name;
  cout << endl;
}
  
  
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print file hash class: 0.." + to_string (hash_class_max - 1))
    {
      version = VERSION;
      addPositional ("name", "File name");
      addFlag ("file", "<file> is a file with names");
      addFlag ("append", "Append the file name");
    }



	void body () const final
	{
	  const string name = getArg ("name");
	  const bool isFile = getFlag ("file");
	  const bool append = getFlag ("append");
	  
	  
	  if (isFile)
	  {
	    LineInput f (name);
	    while (f. nextLine ())
	      if (! strBlank (f. line))
    	    process (f. line, append);
	  }
	  else	  
	    process (name, append);
  }  
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



