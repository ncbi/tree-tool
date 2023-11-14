// hash_request2dissim.cpp

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
*   Compute hash dissimilarities
*
*/


#undef NDEBUG

#include "../common.hpp"
using namespace Common_sp;
#include "evolution.hpp"
using namespace DM_sp;
#include "../version.inc"

#include "../common.inc"



namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Compute hash dissimilarities for pairs of hash files")
    {
      version = VERSION;
    	// Input
  	  addPositional ("pairs", "File with pairs of files");
  	  addKey ("intersection_min", "Min. number of common hashes to compute distance", "50");
  	  addKey ("ratio_min", "Min. ratio of hash sizes (0..1)", "0.5");
  	  // Output
  	  addPositional ("out", "Output file with lines: <obj1> <obj2> <dissimlarity>; <obj1> < <obj2>");
  	}



	void body () const final
	{
		const string pairsFName       = getArg  ("pairs");
		const size_t intersection_min = str2<size_t> (getArg ("intersection_min"));
		const Prob   hashes_ratio_min = str2real (getArg ("ratio_min"));
		const string out              = getArg  ("out");
		ASSERT (isProb (hashes_ratio_min));
		ASSERT (! out. empty ());
		
		
    // Cf. hash2dissim.cpp
    OFStream output (out);
    ONumber on (output, 6, true);  // PAR
    map<string/*fName*/,Hashes> name2hashes;
    PairFile input (pairsFName, false, false);
    while (input. next ())
    {
      const string fName1 (input. name1);
      const string fName2 (input. name2);
      if (! contains (name2hashes, fName1))  name2hashes [fName1] = std::move (Hashes (fName1));
      if (! contains (name2hashes, fName2))  name2hashes [fName2] = std::move (Hashes (fName2));
      const Hashes& h1 = name2hashes [fName1];
      const Hashes& h2 = name2hashes [fName2];
      const double dissim = h1. getDissim (h2, intersection_min, hashes_ratio_min);
      output << input. name1 << '\t' << input. name2 << '\t' << dissim << endl;
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



