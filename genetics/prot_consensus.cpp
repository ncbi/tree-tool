// prot_consensus.cpp

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
*   Print a protein alignment consensus
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "seq.hpp"
using namespace Seq_sp;
#include "../version.inc"



namespace 
{
	
	
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print a protein alignment consensus")
    {
      version = VERSION;
  	  addPositional ("Multi_FASTA", "Multi-FASTA protein file with aligned sequences");
    }


	
	void body () const final
  {
	  const string inFName = getArg ("Multi_FASTA");


	  size_t len = 0;
	  Vector<Peptide> prots;
	  {
  	  Multifasta faIn (inFName, true);
  	  while (faIn. next ())
  	  {
  	    Peptide prot (faIn, 1000, true);  // PAR
  	    prot. qc ();
  	    QC_ASSERT (! prot. seq. empty ());
  	    if (len)
  	    {
  	      if (prot. seq. size () != len)
  	        throw runtime_error ("Unequal sequence lengths");
  	    }
  	    else 
  	      len = prot. seq. size ();
  	    prots << std::move (prot);
  	  }
  	}
  	

    Peptide consensus (getFileName (inFName), len, true);
    FOR (size_t, i, len)
    {
      char c = '\0';
      for (const Peptide& prot : prots)
      {      
  	    if (! c)
  	      c = prot. seq [i];
  	    if (c != prot. seq [i])
  	      c = 'X';
  	  }
  	  consensus. seq [i] = c;
  	}
  	consensus. qc ();

    consensus. saveText (cout);    
  }
};


}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



