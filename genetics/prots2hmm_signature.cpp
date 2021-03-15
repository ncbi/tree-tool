// prot2hmm_signature.cpp

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
*   Print HMM signature for each protein
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
	
	
struct ThisApplication : Application
{
	ThisApplication ()
	  : Application ("Print HMM signature for each protein")
	  {
      version = VERSION;
		  addPositional ("hmmsearch", "Output of hmmsearch");
	  }


	
	void body () const final
  {
	  const string hmmsearch = getArg ("hmmsearch");  

	  	  
  	map <string/*prot_name*/, StringVector/*HMMs*/> prot2hmms;
  	{
	    Hmmsearch hs (hmmsearch);
	    while (hs. next ())
	      prot2hmms [hs. prot_name] << hs. hmm_accession;
	  }
	  
	  for (auto& it : prot2hmms)
	  {
			StringVector& hmms = it. second;
			ASSERT (! hmms. empty ());
			hmms. sort ();
			hmms. uniq ();
	  	cout << it. first << '\t' << hmms. toString (",") << endl;
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



