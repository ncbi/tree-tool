// symbet.cpp  

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
*   Dissimilarity by BLASTP symmetric best hits
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
using namespace DM_sp;
#include "evolution.hpp"
using namespace DM_sp;
#include "../genetics/seq.hpp"
using namespace Seq_sp;
#include "nw/sm_blosum62.c"
#include "../version.inc"



namespace 
{
  
  
unordered_set<string> readSeqs (const string &fName,
                                size_t len_min)
{
  unordered_set<string> seqs;  seqs. rehash (100000);  // PAR
  Multifasta fa (fName, true, 16 * 1024 * 1024);  // PAR
  while (fa. next ())
  {
    Peptide pep (fa, Peptide::stdAveLen, false);  
    pep. qc ();
      // Convert pep.seq to "positives" ??
    if (pep. seq. size () >= len_min)
      seqs. insert (move (pep. seq));
  }
  return seqs;
}



double kmer2weight (const string &kmer)
{
  ASSERT (! kmer. empty ());
  double sum = 0.0;
  for (const char c : kmer)
  {
    const char* c_pos = strchr (NCBISM_Blosum62. symbols, c);
    ASSERT (c_pos);
    const size_t i = (size_t) (c_pos - NCBISM_Blosum62. symbols);
    const double w = s_Blosum62PSM [i * 25 + i];  // Different matrix ??
    ASSERT (w > 0.0);
    sum += w;
  }
  return sum;
}



Vector<size_t> getBests (const unordered_set<string> &seqs1,
                         const unordered_set<string> &seqs2,
                         size_t k)
// Return: size() = seqs1.size()
//         values are indexes of seqs2 or no_index
{
  ASSERT (k);
  
  unordered_map<string/*kmer*/,Vector<size_t/*id*/>> kmer2ids;  kmer2ids. rehash (seqs2. size () * 100);  // PAR
  {
    size_t id = 0;
    for (const string& seq : seqs2)
    {
      const size_t seqSize = seq. size ();
      FOR (size_t, i, seqSize)
        if (i + k <= seqSize)
        {
          const Peptide kmer ("x", seq. substr (i, k), false);
          kmer. qc ();
          ASSERT (kmer. seq. size () == k);
          if (kmer. getXs ())
            continue;
          kmer2ids [kmer. seq] << id;
        }
      id++;
    }
  }

  Vector<size_t> bests;  bests. reserve (seqs1. size ());
  {
    unordered_map<size_t/*id*/,double> id2weight;  id2weight. rehash (100000);  // PAR
    for (const string& seq : seqs1)
	  {
	    id2weight. clear ();
      const size_t seqSize = seq. size ();
      FOR (size_t, i, seqSize)
        if (i + k <= seqSize)
        {
          const Peptide kmer ("x", seq. substr (i, k), false);
          kmer. qc ();
          ASSERT (kmer. seq. size () == k);
          if (kmer. getXs ())
            continue;
          if (const Vector<size_t>* ids = findPtr (kmer2ids, kmer. seq))
            for (const size_t id : *ids)
            {
              const double weight = kmer2weight (kmer. seq);
              ASSERT (weight > 0.0);
              id2weight [id] += weight;
            }
        }
      size_t id_best = no_index;
      double weight_max = 0;
      for (const auto& it : id2weight)
        // ignore small it.second ??
        if (maximize (weight_max, it. second))
          id_best = it. first;
      bests << id_best;
	  }
	}
	ASSERT (bests. size () == seqs1. size ());
	
	return bests;
}




struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print dissimilarity by k-mer symmetric best hits: >= 0 or nan")
    {
      version = VERSION;
  	  addPositional ("fasta1", "Protein FASTA file 1");
  	  addPositional ("fasta2", "Protein FASTA file 2");
  	  addKey ("k", "k-mer size, >= 3", "5");
  	  addKey ("min_prot_len", "Min. protein length", "0");
  	}



	void body () const final
	{
		const string fName1  = getArg ("fasta1");
		const string fName2  = getArg ("fasta2");
		const size_t k       = (size_t) arg2uint ("k");
		const size_t len_min = (size_t) arg2uint ("min_prot_len");
		QC_ASSERT (k >= 3);  // PAR
		
		
		const unordered_set<string> seqs1 (readSeqs (fName1, len_min));
		const unordered_set<string> seqs2 (readSeqs (fName2, len_min));
		
		const Vector<size_t> bests1 (getBests (seqs1, seqs2, k));
		const Vector<size_t> bests2 (getBests (seqs2, seqs1, k));
		
		// For approximate k-mer symbets run Needleman-Wunsch ??!
		
		size_t symbet = 0;
		FFOR (size_t, i, bests1. size ())
		  if (bests1 [i] != no_index)
		    if (bests2 [bests1 [i]] == i)
		      symbet++;


	  if (verbose ())
	  {
  	  PRINT (seqs1. size ());
  	  PRINT (seqs2. size ());
  	  PRINT (symbet);
  	}
	  cout << intersection2dissim ( (Real) seqs1. size ()
	                              , (Real) seqs2. size ()
	                              , (Real) symbet
	                              , 50.0   // PAR
	                              , 0.5    // PAR
	                              , true)  
	       << endl; 
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



