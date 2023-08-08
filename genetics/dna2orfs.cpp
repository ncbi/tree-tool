// dna2stat.cpp

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
*   Print all ORFs between stop codons
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
    : Application ("Print all ORFs between stop codons")
	  {
      version = VERSION;
	  	addPositional ("in", "Input DNA multi-FASTA file");
	  	addPositional ("gencode", "NCBI genetic code");
	  	addPositional ("len_min", "Min. ORF length in aa without 'X'");
	  }



  void body () const final
  {
    const string  inFName = getArg ("in");
    const Gencode gencode = (Gencode) str2<uint> (getArg ("gencode"));
    const size_t len_min  = str2<size_t> (getArg ("len_min"));


    Vector<Peptide> peps;  peps. reserve (10000); // PAR
	  {
		  Multifasta faIn (inFName, false);
		  while (faIn. next ())
		  {
		    Dna dna (faIn, 1e6/*PAR*/, false);
	    	strLower (dna. seq);
			  dna. qc ();
			  for (const int frame : {-3, -2, -1, 1, 2, 3})
			    peps << dna. getOrfs ((Frame) frame, gencode, len_min);
			}
		}		     
    for (const Peptide& pep : peps)
    {
      pep. saveText (cout);
      pep. qc ();
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



