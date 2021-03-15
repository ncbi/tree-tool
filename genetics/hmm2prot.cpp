// hmm2prot.cpp

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
*   Find new seeds for a set of HMMs
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
using namespace DM_sp;
#include "hmm.hpp"
using namespace Hmm_sp;
#include "../version.inc"




namespace 
{
  

struct Hit : Named
// name: hmmName
{
  double score1 {0};
    // > 0

  Hit ()
    {}
  Hit (const string &hmmName,
       double score1_arg
      )
    : Named (hmmName)
    , score1 (score1_arg)
    { ASSERT (score1 > 0); }
};



struct Hmm : Root
{
  StringVector sseqids;  
};

  

  
struct ThisApplication : Application
{
	ThisApplication ()
	  : Application ("Find new seeds for a set of HMMs")
	  {
      version = VERSION;
		  addPositional ("hmmsearch", "Output of hmmsearch");
		  addPositional ("dir", "Directory for output files with ordered seed ids");
	  }


	
	void body () const final
  {
	  const string hmmsearch = getArg ("hmmsearch");  
	  const string dir       = getArg ("dir");  
	  
	  
	  map<string/*prot_name*/,Hit> prot2bestHmm;
	  {
	    Hmmsearch hs (hmmsearch);
	    while (hs. next ())
	    {
	      const Hit hit (prot2bestHmm [hs. prot_name]);
	      if (    hs. score1 >  hit. score1
	          || (hs. score1 == hit. score1 && hs. hmm_name < hit. name)
	         )
	        prot2bestHmm [hs. prot_name] = Hit (hs. hmm_name, hs. score1);
	    }
	  }

    MeanVar mv (1); 	// PAR  
	  map<string/*hmm_name*/,Hmm> hmms;  
	  for (const auto &p : prot2bestHmm)
	  {
	    const Hit& hit = p. second;
	    hmms [hit. name] = Hmm ();
	    mv << hit. score1;
	  }
	  cout << mv. getMean () << endl;
	    
	  for (const auto &p : prot2bestHmm)
	  {
	    const Hit& hit = p. second;
	    Hmm& hmm = hmms [hit. name];
	    hmm. sseqids << p. first;
	  }

	  for (auto &it : hmms)
	  {
	    Hmm& hmm = it. second;
  	  sort (hmm. sseqids);
  	  OFStream ofs (dir, it. first, string ());
  	  save (ofs, hmm. sseqids, '\n');
  	  ofs << endl;
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



