// distTree_refresh_dissim.cpp

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
*   Compute files with new dissimilarities
*
*/


#undef NDEBUG

#include "../common.hpp"
using namespace Common_sp;
#include "distTree.hpp"
using namespace DistTree_sp;
#include "../version.inc"

#include "../common.inc"



namespace 
{
  
 
struct ThisApplication : Application
{
	ThisApplication ()
	: Application ("Create files with representative dissimilarities")
	{
	  version = VERSION;
		// Input
	  addPositional ("dir", "Directory with data for an incremental tree ending with '/'");
    // Output
	  addPositional ("dissim_request", "Output file with requests to compute needed dissimilarities, tab-delimited line format: <obj1> <obj2>");
	  addPositional ("output_dissim", "Output file with dissimilarities used in the tree, tab-delimited line format: <obj1> <obj2> <dissimilarity>");
	}
	
	
	
	void body () const final
  {
	  const string dirName        = getArg ("dir");
		const string dissim_request = getArg ("dissim_request");
		const string output_dissim  = getArg ("output_dissim");

    QC_ASSERT (! dirName. empty ());
    QC_ASSERT (! dissim_request. empty ());
    QC_ASSERT (! output_dissim. empty ());

		
    const DistTree tree (dirName + "/tree", string (), string (), string ());
    tree. qc ();
    
    const Vector<LeafPair> leafPairs (tree. getMissingLeafPairs_ancestors (sparsingDepth, true));
    const Vector<DissimLine> dissimLines (tree. getDissimLines (dirName + "/dissim", 100000000));  // PAR
    
    OFStream fRequest (dissim_request);
    OFStream fDissim (output_dissim);
    const ONumber on (fDissim, dissimDecimals, true);
    {
      section ("Saving indiscernibles", false);
      Progress prog (dissimLines. size (), dissim_progress);
      for (const DissimLine &dl : dissimLines)
      {
        prog ();
        if (dl. dissim == 0.0)
          fDissim << dl. name1 << '\t' << dl. name2 << '\t' << dl. dissim << endl;
      }
    }
    {
      section ("Saving requests", false);
      Progress prog (leafPairs. size (), dissim_progress);
      for (const LeafPair& lp : leafPairs)
      {
        prog ();
        const DissimLine bait ( lp. first ->name
                              , lp. second->name
                              );
        const size_t i = dissimLines. binSearch (bait);  // use unordered_map ??
        if (i == no_index)
          fRequest << bait. name1 << '\t' << bait. name2 << endl;
        else
        {
          const DissimLine& fish = dissimLines [i];
          if (fish. dissim != 0.0)
            fDissim << fish. name1 << '\t' << fish. name2 << '\t' << fish. dissim << endl;
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


