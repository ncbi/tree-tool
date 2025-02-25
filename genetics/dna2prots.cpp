// dna2prots.cpp

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
*   Print all longest proteins between start and stop codons
*
*/


#undef NDEBUG

#include "../common.hpp"
using namespace Common_sp;
#include "seq.hpp"
using namespace Seq_sp;
#include "../version.inc"

#include "../common.inc"




namespace 
{


struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Print all longest proteins between start and stop codons")
	  {
      version = VERSION;
	  	addPositional ("in", "Input DNA multi-FASTA file");
	  	addPositional ("gencode", "NCBI genetic code");
	  	addPositional ("len_min", "Min. protein length in aa");
		  addKey ("complexity_min", "Min. protein sequence complexity", "0");
		  addFlag ("no_x", "Suppress sequences with ambiguous aa");
	  }



  void body () const final
  {
    const string  inFName       = getArg ("in");
    const Gencode gencode       = (Gencode) str2<uint> (getArg ("gencode"));
    const size_t len_min        = str2<size_t> (getArg ("len_min"));
	  const double complexity_min = arg2double ("complexity_min");
	  const bool   allowX         = ! getFlag ("no_x");


    Vector<Peptide> peps;  peps. reserve (10000); // PAR
	  {
		  Multifasta faIn (inFName, false, 1);
		  while (faIn. next ())
		  {
		    Dna dna (faIn, 1e6/*PAR*/, false);
	    	strLower (dna. seq);
			  dna. qc ();
			  for (const int frame : {-3, -2, -1, 1, 2, 3})
			    peps << dna. getPeptides ((Frame) frame, gencode, len_min);
			}
		}		 
		    
    for (const Peptide& pep : peps)
    {
      pep. qc ();
      ASSERT (pep. seq. size () >= len_min);
      if (   pep. getComplexity () >= complexity_min
          && (allowX || ! pep. getXs ())
         )
        pep. saveText (cout);
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



