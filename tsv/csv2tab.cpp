// csv2tab.cpp

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
*   Convert a .csv-file to a tab-delimited file
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
    : Application ("Convert a .csv-file to a tab-delimited file")
    {
      version = VERSION;
  	  addPositional ("in", "Text file");
  	  addFlag ("tsv", "First line is a header");
  	}


	void body () const final
	{
		const string in  = getArg ("in");
		const bool   tsv = getFlag ("tsv");


    LineInput f (in);  
    while (f. nextLine ())
    {
      if (contains (f. line, '\t'))
        throw runtime_error ("File contains tabs");
      if (tsv && f. tp. lineNum == 1)
        cout << '#';
      bool quote = false;
      for (const char c : f. line)
        if (c == '\"')
          toggle (quote);
        else
          if (! quote && c == ',')
            cout << '\t';
          else
            cout << c;
      cout << endl;
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
