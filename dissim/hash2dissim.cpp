// hash2dissim.cpp  

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
*   Convert hashes to a dissimilarity 
*
*/


#undef NDEBUG

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
#include "../dm/dataset.hpp"
using namespace DM_sp;
#include "evolution.hpp"
using namespace DM_sp;
#include "../version.inc"

#include "../common.inc"



namespace 
{


const string attrName = "cons";



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Convert hashes to a dissimilarity named " + strQuote (attrName) + " and print a " + dmSuff + "-file")
    {
      version = VERSION;
    	// Input
  	  addPositional ("objects", "File with a list of objects");
  	  addPositional ("hash_dir", "Directory with hashes for each object");
  	  addKey ("intersection_min", "Min. number of common hashes to compute distance", "50");
  	  addKey ("ratio_min", "Min. ratio of hash sizes (0..1)", "0.5");
  	  // Output
  	  addPositional ("out", "Output " + dmSuff + "-file without " + dmSuff);
  	}



	void body () const final
	{
		const string objectsFName     = getArg  ("objects");
		const string hash_dir         = getArg  ("hash_dir");
		const size_t intersection_min = str2<size_t> (getArg ("intersection_min"));
		const Prob   hashes_ratio_min = str2real (getArg ("ratio_min"));
		const string out              = getArg ("out");
		ASSERT (isProb (hashes_ratio_min));
		ASSERT (! out. empty ());
		
		
    Dataset ds;
    {
  		Set<string> objNames;
      {
        LineInput f (objectsFName);
        while (f. nextLine ())
        {
          trim (f. line);
          objNames << std::move (f. line);
        }
      }
      cerr << "# Objects: " << objNames. size () << endl;  
      for (const string& name : objNames)
        ds. appendObj (name);
    }        
    ds. qc ();
    
    
    Vector<Hashes> obj2hashes;  obj2hashes. reserve (ds. objs. size ());
    {
      Progress prog (ds. objs. size ());
      FFOR (size_t, objNum, ds. objs. size ())
      {
        prog (ds. objs [objNum] -> name);
        obj2hashes << std::move (Hashes (hash_dir + "/" + ds. objs [objNum] -> name));
      }
    }
    
    auto attr = new PositiveAttr2 (attrName, ds, 6);  // PAR
    {
      Progress prog (ds. objs. size ());
      FFOR (size_t, row, ds. objs. size ())
      {
        prog (ds. objs [row] -> name);
        const Hashes& hash1 = obj2hashes [row];
        FOR (size_t, col, row)
        {
          const Hashes& hash2 = obj2hashes [col];
          const Real dissim = hash1. getDissim (hash2, intersection_min, hashes_ratio_min);
          attr->putSymm (row, col, dissim);
        }
        attr->put (row, row, 0);
      }
    }
    
    ds. qc ();
    {
      OFStream f (out + dmSuff);
      ds. saveText (f);    
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



