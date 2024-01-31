// tsv2triple.cpp

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
*   Convert a tsv-file into triples 
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
    : Application ("Convert a tsv-file into an object-attribute-value table")
  	{
      version = VERSION;
  	  addPositional ("tsv", "tsv-table");
  	  addKey ("keys", "Object columns in <tsv>, comma-delimited");
  	}
  	
  	
 
	void body () const final
	{
		const string fName = getArg ("tsv");
		const string keysS = getArg ("keys");

    
    const TextTable tab (fName);
    tab. qc ();

    Vector<TextTable::ColNum> keys;
    {
      const StringVector keysVec (keysS, ',', true);
      for (const string& s : keysVec)
        keys << tab. col2num (s);
    }
      
    for (const StringVector& row : tab. rows)
    {
      string keyS;
      bool first = true;
      for (const size_t i : keys)
      {
        if (first)
          first = false;
        else
          keyS += '\t';
        keyS += row [i];
      }
      if (! keyS. empty ())
        keyS += '\t';
      FFOR (size_t, i, row. size ())
        if (   ! keys. contains (i)
            && ! row [i]. empty ()
           )
          cout << keyS << tab. header [i]. name << '\t' << row [i] << endl;
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




