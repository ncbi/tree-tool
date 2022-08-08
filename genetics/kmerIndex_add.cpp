// kmerIndex_add.cpp

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
*   Add DNA sequences to a DNA k-mer index
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
    : Application ("Add DNA sequences to a DNA k-mer index.\n\
Time: O(n), where n is the number of sequences, assuming constant average length of sequences and sequence identifiers\
")
    {
      version = VERSION;
  	  addPositional ("kmer_index", "K-mer index file name");
  	  addPositional ("Multi_FASTA", "Multi-FASTA file with DNA sequences");  
    }


	
	void body () const final
  {
	  const string kmerFName = getArg ("kmer_index");
	  const string inFName   = getArg ("Multi_FASTA");


    KmerIndex kmi (kmerFName);
    kmi. qc ();
    
	  size_t n = 0;
	  size_t kmers_total = 0;
	  size_t kmersRejected_total = 0;
    {
		  Multifasta faIn (inFName, false, 100);  // PAR
		  while (faIn. next ())
		  {
		    const Dna dna (faIn, 10000, false);  // PAR
		    dna. qc ();
		    size_t kmers = 0;
        size_t kmersRejected = 0;
		    kmi. add (dna, kmers, kmersRejected);
		    ASSERT (kmers >= kmersRejected);
  	    if (kmers == kmersRejected)
  	      throw runtime_error ("DNA is not indexed: " + dna. getId ());
  	    kmers_total         += kmers;
  	    kmersRejected_total += kmersRejected;
  	    n++;
		  }
		}
    kmi. qc ();
		cout << "# Sequences: " << n << endl;
		cout << "Average # k-mers per sequence: "          << (double) kmers_total         / (double) n << endl;
		cout << "Average # rejected k-mers per sequence: " << (double) kmersRejected_total / (double) n << endl;
  }
};


}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



