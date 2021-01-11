// kmerIndex_find.cpp

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
*   Find closest sequences using a DNA k-mer index
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "seq.hpp"
using namespace Seq_sp;



namespace 
{
  
  

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Find close sequences using a DNA k-mer index.\n\
Time: O(log n), where n is the number of DNA sequences in the index\
")
    {
  	  addPositional ("kmer_index", "DNA k-mer index file name");
  	  addPositional ("FASTA", "Input DNA sequence");  
  	  addPositional ("top", "Number of closest sequence identifiers to print");
  	  addFlag ("common_kmers", "Print number of common k-mers");
  	  addFlag ("self", "Do not remove the input DNA identifier");
    }


	
	void body () const final
  {
	  const string kmerFName  = getArg ("kmer_index");
	  const string inFName    = getArg ("FASTA");
	  const size_t top        = (size_t) arg2uint ("top");
	  const bool common_kmers = getFlag ("common_kmers");
	  const bool self         = getFlag ("self");


    KmerIndex kmi (kmerFName);
    kmi. qc ();
    
    LineInput li (inFName);
    EXEC_ASSERT (li. nextLine ());
    const Dna dna (li, 10000, false);  // PAR
    dna. qc ();

    const Vector<KmerIndex::NumId> numIds (kmi. find (dna));

  #ifndef NDEBUG
    size_t num_prev = numeric_limits<size_t>::max ();
  #endif
    const string id (dna. getId ());
    size_t n = 0;
    for (const KmerIndex::NumId& numId : numIds)
    {
      if (n >= top)
        break;
      if (self || numId. id != id)
      {
        cout << numId. id;
        if (common_kmers)
          cout << '\t' << numId. n;
        cout << endl;
        n++;
      }
      ASSERT (num_prev >= numId. n);
    #ifndef NDEBUG
      num_prev = numId. n;
    #endif
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



