// prot_clust.cpp

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
*   Cluster proteins by hash-based Jaccard index
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
  

hash<string> str_hash;

  
  
struct HashPep : Named
{
  // Sorted, unique
  Vector<size_t> jaccardHashes;
  Vector<size_t> indexHashes;
  
  
  HashPep (const Peptide &p,
           size_t jaccard_k,
           size_t index_k)
    : Named (p. getId ())
    { 
      ASSERT (jaccard_k);
      ASSERT (jaccard_k <= index_k);
      setHashes (p, jaccardHashes, jaccard_k);
      setHashes (p, indexHashes,   index_k);  
    }
private:
  void setHashes (const Peptide &p,
                  Vector<size_t> &hashes,
                  size_t k)
    {
      ASSERT (k);
      ASSERT (hashes. empty ());
      ASSERT (p. seq. size () >= k);
      hashes. reserve (p. seq. size () - k + 1);
      string s;
      FFOR (size_t, i, p. seq. size () - k + 1)
      {
        s = p. seq. substr (i, k);
        bool ambig = false;
        for (const char c : s)
          if (isAmbigAa (c))
          {
            ambig = true;
            break;
          }
        if (ambig)
          continue;
        hashes << str_hash (s);
      }
      hashes. sort ();
      hashes. uniq ();
    }
};



double getJaccard (const HashPep &p1,
                   const HashPep &p2)
// Time: O (|p1| + |p2|)
// Symmetric
{
  const size_t size1 = p2. jaccardHashes. size ();
  const size_t size2 = p2. jaccardHashes. size ();
  ASSERT (size1);
  ASSERT (size2);
  
  size_t common = 0;
  {
    size_t i = 0;
    for (const size_t h1 : p1. jaccardHashes)
    {
      while (   i < size2 
             && p2. jaccardHashes [i] < h1
            )
        i++;
      if (i == size2)
        break;
      if (p2. jaccardHashes [i] == h1)
      {
        common++;
        i++;
      }
    }
  }
  ASSERT (size1 >= common);
  ASSERT (size2 >= common);
  return (double) common / (double) (size1 + (size2 - common));
}

	
	
struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Print: <protein 1> <protein 2> <Jaccard index of protein sequence hashes>. Time: 3.5 hours/5.5 M sequences")
    {
      version = VERSION;
  	  addPositional ("in", "Protein FASTA file");
  	  addPositional ("jaccard_min", "Min. <Jaccard index>. 0.5 matches 85% identity, 0.4 - 50% identity");
  	  addPositional ("out", "Matches: <id1> <id2> <Jaccard index>");
  	  addKey ("jaccard_k", "K-mer length for Jaccard index", "4");
  	  addKey ("index_k", "K-mer length for indexing, larger than <jaccard_k>", "20");
  	  addKey ("freq_max", "Max. relative frequency of an index k-mer. Set to 1 for remote matches", "0.001");
    }


	
	void body () const final
  {
	  const string inFName     = getArg ("in");
	  const double jaccard_min = arg2double ("jaccard_min");
	  const string outFName    = getArg ("out");
	  const size_t jaccard_k   = (size_t) arg2uint ("jaccard_k");
	  const size_t index_k     = (size_t) arg2uint ("index_k");
	  const double relFreq_max = arg2double ("freq_max");
	  
	  QC_ASSERT (jaccard_min >= 0.0);
	  QC_ASSERT (jaccard_k >= 3);  // PAR
	  QC_ASSERT (jaccard_k <= index_k);
	  QC_ASSERT (relFreq_max > 0.0);
	  QC_ASSERT (relFreq_max <= 1.0);
	  
    
    Vector<HashPep> hashPeps;  hashPeps. reserve (10000);  // PAR
    map<size_t/*index k-mer hash*/,VectorPtr<HashPep>> kmer2hps;
    {
      OFStream fOut (outFName);

      // hashPeps
      {
    	  Multifasta faIn (inFName, true);
    	  while (faIn. next ())
    	  {
      	  const Peptide peptide (faIn, 1024 * 1024, false);    	  
    	    peptide. qc ();
    	    fOut         << peptide. getId () 
    	         << '\t' << peptide. getId () 
    	         << '\t' << 1
    	         << '\n';
    	    if (peptide. seq. size () < index_k)
    	      continue;  
    	    HashPep hp (peptide, jaccard_k, index_k);
    	    if (hp. indexHashes. empty ())
    	      continue;  
    	    ASSERT (! hp. jaccardHashes. empty ());
    	    hashPeps << std::move (hp);
      	}
      }
      cout << "# Proteins: " << hashPeps. size () << endl;
        
      section ("Indexing", false);
      // kmer2hps
      {
        Progress prog (hashPeps. size (), 1000);  // PAR
        for (const HashPep& hp : hashPeps)
        {
          prog ();
          for (const size_t h : hp. indexHashes)
            kmer2hps [h] << & hp;
        }
      }
      
    #if 0
      {
        OFStream f ("hashes");
        for (const auto& it : kmer2hps)
          f << it. second. size () << '\n';
      }
    #endif
        
      section ("Searching", false);
      {
        Progress prog (hashPeps. size (), 1000);  // PAR
        VectorPtr<HashPep> neighbors;  neighbors. reserve (1000);  // PAR
        const size_t freq_max = size_t (relFreq_max * (double) hashPeps. size ());
        for (const HashPep& hp1 : hashPeps)
        {
          prog ();
          neighbors. clear ();
          for (const size_t h : hp1. indexHashes)
            if (kmer2hps [h]. size () < freq_max)  // Speed improvement!
              for (const HashPep* hp2 : kmer2hps [h])
                if (hp1. name < hp2->name)
                  neighbors << hp2;
          neighbors. sort ();
          neighbors. uniq ();
          if (verbose ())
            cerr << hp1. name << "\tneighbors=" << neighbors. size () << '\n';
          for (const HashPep* hp2 : neighbors)
          {
            const double jaccard = getJaccard (hp1, *hp2);
            if (jaccard >= jaccard_min)
              fOut         << hp1. name 
                   << '\t' << hp2->name 
                   << '\t' << jaccard
                   << '\n';
          }
        }
      }
    }
      
    quick_exit (0);  // Freeing memory takes too much time
  }
};


}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



