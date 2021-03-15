// hmmsearch2besthits.cpp

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
*   Select a best protein hit for each HMM. Print: <protein> <HMM> <score> <ali_from> <ali_to>
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "hmm.hpp"
using namespace Hmm_sp;
#include "../version.inc"



namespace
{
  
  
struct Match
{
	string prot_name;
	double score {0.0};
	size_t ali_from {(size_t) -1};
	size_t ali_to {0};
	
	Match (const string &prot_name_arg,
	       double score_arg)
	  : prot_name (prot_name_arg)
	  , score (score_arg)
	  { ASSERT (score > 0.0); }
	Match ()
	  {}
	  
	bool valid () const
	  { return ali_to > 0; }
};
  
  
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Select a best protein hit for each HMM. Print: <protein> <HMM> <score> <ali_from> <ali_to>")
  	{
      version = VERSION;
  	  addPositional ("tblout", "hmmsearch --tblout result");
  	  addKey ("domtblout", "hmmsearch --domtblout result");
  	}
  


	void body () const final
	{
		const string tbloutFName    = getArg ("tblout");
		const string domtbloutFName = getArg ("domtblout");
		
    
    map <string/*hmm*/, Match> hmm2match;  
    {
	    Hmmsearch hs (tbloutFName);
      while (hs. next ())
        if (hmm2match [hs. hmm_name]. score < hs. score1)
        	hmm2match [hs. hmm_name] = Match (hs. prot_name, hs. score1);
    }
    
    if (! domtbloutFName. empty ())
    {
	    HmmDom hd (domtbloutFName);
      while (hd. next ())
        if (const Match* m = findPtr (hmm2match, hd. hmm_name))
          if (m->prot_name == hd. prot_name)
        	{
        		minimize (var_cast (m) -> ali_from, hd. ali_from);
        		maximize (var_cast (m) -> ali_to,   hd. ali_to);
        	}
    }

    for (const auto& it : hmm2match)
    {
      const Match& m = it. second;
    	if (m. valid ())
	    	cout         << m. prot_name 
	    	     << '\t' << it. first 
    	  	   << '\t' << m. score 
    	  	   << '\t' << m. ali_from + 1
    	  	   << '\t' << m. ali_to
	    	     << endl;
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


