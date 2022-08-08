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
  Multifasta fa (fName, true);  
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



struct TopMatches
{
  struct Item
  {
    size_t id {no_index};
    double weight {NaN};
    bool operator< (const Item &other) const
      { return weight > other. weight; }
  };
  size_t size_max {0};
  Vector<Item> items;
    // size() <= size_max
    // Orderted by Item::weight descending
  

  explicit TopMatches (size_t size_max_arg)
    : size_max (size_max_arg)
    {
      ASSERT (size_max);
      items. reserve (size_max);
    }


  void add (size_t id,
            double weight)
    { 
      ASSERT (id != no_index);
      ASSERT (weight > 0.0);
      if (! items. empty () && greaterReal (items. back (). weight, weight))
        return;
      ASSERT (items. size () <= size_max);
      if (items. size () == size_max)
        items. pop_back ();
      items << Item {id, weight};
      items. sortBubble ();
    }
  Vector<size_t> getIds () const
    // Sorted
    {
      Vector<size_t> ids;  ids. reserve (items. size ());
      for (const Item& item : items)
        ids << item. id;
      ids. sort ();
      return ids;
    }
};



typedef  Vector<Vector<size_t>>  Bests;



Bests getBests (const unordered_set<string> &seqs1,
                const unordered_set<string> &seqs2,
                size_t k,
                size_t ploidy)
// Return: size() = seqs1.size()
//         values are indexes of seqs2 or no_index
{
  ASSERT (k);
  ASSERT (ploidy);
  
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

  Bests bests;  bests. reserve (seqs1. size ());
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
      TopMatches tm (ploidy); 
      for (const auto& it : id2weight)
        tm. add (it. first, it. second);
          // Skip if it.second is too small ??
      bests. push_back (tm. getIds ());
	  }
	}
	ASSERT (bests. size () == seqs1. size ());
	
	return bests;
}



size_t getMaps (const Bests &bests1,
                const Bests &bests2)
{
	size_t maps = 0;
	FFOR (size_t, i, bests1. size ())
	  for (const size_t j : bests1 [i])
	    if (bests2 [j]. containsFast (i))
	    {
	      maps++;
	      break;
	    }
	return maps;
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
  	  addKey ("ploidy", "Number of chromosome copies", "1");
  	}



	void body () const final
	{
		const string fName1   = getArg ("fasta1");
		const string fName2   = getArg ("fasta2");
		const size_t k        = (size_t) arg2uint ("k");
		const size_t len_min  = (size_t) arg2uint ("min_prot_len");
		const size_t ploidy = (size_t) arg2uint ("ploidy");
		QC_ASSERT (k >= 3);  // PAR
		QC_ASSERT (ploidy >= 1);
		
		
		const unordered_set<string> seqs1 (readSeqs (fName1, len_min));
		const unordered_set<string> seqs2 (readSeqs (fName2, len_min));		
		const Real size1 = (Real) seqs1. size ();
		const Real size2 = (Real) seqs2. size ();
		
		Real dissim = NaN;
		constexpr Real sizes_ratio_min = 0.5;  // PAR
    if (  min (size1, size2) 
    	  / max (size1, size2)
  	    >= sizes_ratio_min   // Cf. maps2dissim()
  	   )
  	{
  		const Bests bests1 (getBests (seqs1, seqs2, k, ploidy));
  		const Bests bests2 (getBests (seqs2, seqs1, k, ploidy));
  		// For approximate k-mer symbets run Needleman-Wunsch ??!
      const size_t maps1 = getMaps (bests1, bests2);
      const size_t maps2 = getMaps (bests2, bests1);
  	  if (verbose ())
  	  {
    	  PRINT (seqs1. size ());
    	  PRINT (seqs2. size ());
    	  PRINT (maps1);
    	  PRINT (maps2);
    	}
    	dissim = maps2dissim ( size1
                           , size2
                           , (Real) maps1
                           , (Real) maps2
                           , 50.0   // PAR
                           , sizes_ratio_min
                           , true);
    }
  	
	  cout << dissim << endl; 
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



