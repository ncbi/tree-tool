// tsv_test_key.cpp

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
*   Test a key of a tsv-file
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
  
  
struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Test a key of a tsv-file, print column numbers (1-based)")
  	{
      version = VERSION;
  	  addPositional ("tsv", "tsv-table");
  	  addPositional ("keys", "Object columns in <tsv>, comma-delimited");
  	}
  	
  	
 
	void body () const final
	{
		const string fName = getArg ("tsv");
		const string keysS = getArg ("keys");

    
    const TextTable tab (fName, noString, true, 10000);  // PAR
    tab. qc ();

    const StringVector keysVec (keysS, ',', true);
    const Vector<TextTable::ColNum> keyCols (tab. columns2nums (keysVec));
    
    {
      bool first = true;
      for (const TextTable::ColNum c : keyCols)
      {
        if (! first)
          cout << ',';
        cout << c + 1;
        first = false;
      }
      cout << endl;
    }
      
    Set<string> keys;
    for (const StringVector& row : tab. rows)
    {
      string key;
      bool first = true;
      for (const size_t i : keyCols)
      {
        if (first)
          first = false;
        else
          key += '\t';
        key += row [i];
      }
      if (! keys. insert (key). second)
        throw runtime_error ("Duplicate key: " + key);
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




