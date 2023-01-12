// alignment2seq.cpp

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
*   Print a DNA alignment consensus
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
    : Application ("Print a DNA alignment consensus")
    {
      version = VERSION;
  	  addPositional ("Multi_FASTA", "Multi-FASTA file with aligned DNA sequences");
    }


	
	void body () const final
  {
	  const string inFName = getArg ("Multi_FASTA");


	  size_t len = 0;
	  Vector<Dna> dnas;
	  {
  	  Multifasta faIn (inFName, false);
  	  while (faIn. next ())
  	  {
  	    Dna dna (faIn, 1000, true);  // PAR
  	    dna. qc ();
  	    QC_ASSERT (! dna. seq. empty ());
  	    if (len)
  	    {
  	      if (dna. seq. size () != len)
  	        throw runtime_error ("Unequal sequence lengths");
  	    }
  	    else 
  	      len = dna. seq. size ();
  	    dnas << move (dna);
  	  }
  	}
  	

    Dna consensus (getFileName (inFName), len, true);
    FOR (size_t, i, len)
    {
      constexpr size_t n = 5;
      bool acgtb_sum [n] {false, false, false, false, false};
      for (const Dna& dna : dnas)
      {      
        bool acgtb [n] {false, false, false, false, false};
  	    wild2nucleotides (dna. seq [i], acgtb);
  	    FOR (size_t, j, n)
  	      if (acgtb [j])
  	        acgtb_sum [j] = true;  	    
  	  }
  	  consensus. seq [i] = nucleotides2wild (acgtb_sum);
  	}

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



