// splitFasta.cpp

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
*   Split a multi-fasta file into a set of fasta files 
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
  
  
const string splitSuffix ("-split");  // PAR
string seqSuffix;

  
  
struct Group
{
  const size_t seqCount_max;
  const size_t seqSize_max;
  const string dir;
  VectorOwn<const Seq> seqs;  
  size_t seqSize {0};
  size_t group {0};
  
  
  Group (size_t seqCount_max_arg,
         size_t seqSize_max_arg,
         const string &dir_arg)
    : seqCount_max (seqCount_max_arg ? seqCount_max_arg : numeric_limits<size_t>::max ())
    , seqSize_max  (seqSize_max_arg  ? seqSize_max_arg  : numeric_limits<size_t>::max ())
    , dir (dir_arg)
    { ASSERT (! dir. empty ()); }
 ~Group ()
    { save (); }
  
  
  void add (Seq* seq)
    { ASSERT (seq);
      ASSERT (! seq->seq. empty ());
      ASSERT (seqs. size () <= seqCount_max);
      ASSERT (seqSize <= seqSize_max);
      if (   seqs. size () == seqCount_max
          || (   seqSize + seq->seq. size () > seqSize_max
              && seqSize >= seqSize_max / 2
             )
         )
        save ();
      size_t n = 0;
      while (seqSize + seq->seq. size () > seqSize_max)
      {
        Seq* prefix = seq->copy ();
        n++;
        prefix->name = seq->getId () + ":" + to_string (n) + splitSuffix;  
        ASSERT (seqSize_max >= seqSize);
        const size_t len = seqSize_max - seqSize;
        ASSERT (len < seq->seq. size ());
        ASSERT (len >= seqSize_max / 2);
        prefix->seq. erase (len);
        seqs << prefix;
        save ();
        seq->seq. erase (0, len);
      }
      seqs << seq;
      seqSize += seq->seq. size ();
    }
  void save ()
    { if (seqs. empty ())
        return;
      group++;
      OFStream f (dir + "/" + to_string (group));
      for (const Seq* seq : seqs)
      {
        var_cast (seq) -> appendId (seqSuffix);
        seq->saveText (f);
      }
      seqs. deleteData ();
      seqSize = 0;
    }
};
  
  
  
struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Split a multi-fasta file into a set of fasta files each containing one sequence")
  	{
      version = VERSION;
  	  addPositional ("in", "Input multi-FASTA file");
  	  addPositional ("out_dir", "Output directory");
  	  addFlag ("aa", "Multi-FASTA file contains protein sequences, otherwise DNA sequences");
  	  addFlag ("sparse", "Sparse sequence");
  	  addKey ("len_min", "Minimum sequence length", "0");
		  addFlag ("whole", "Sequence identifiers are whole strings which should not be split by '|'");
  	  addFlag ("pseudo", "Pseudo-proteins are allowed");
  	  addKey ("mono_nuc_max", "Mask mononucleotide repeats longer than <mono_nuc_max>. 0 = infinity", "0");
    #ifndef _MSC_VER
  	  addFlag ("large", "Create files in subdirectories \"0\" .. \"" + to_string (small_hash_class_max - 1) + "\" which are the hashes of file names");
  	#endif
  	  addKey ("group_count", "Group by <group_count> number of sequences. Group names are sequential numbers. 0 - no grouping.", "0");
  	  addKey ("group_size", "Make groups of sequences with total sequence size <= <group_size> with possible splitting sequences. Group names are sequential numbers. 0 - no grouping.", "0");
  	  addKey ("extension", "Add file extension (not compatible with -group_count or -group_size)", "");
  	  addKey ("suffix", "Suffix to add to each sequence separated by dash", "");
  	}



	void body () const final
  {
		const string in           = getArg ("in");
		const string out_dir      = getArg ("out_dir");
		const bool   aa           = getFlag ("aa");
		const bool   sparse       = getFlag ("sparse");
		const size_t len_min      = (size_t) arg2uint ("len_min");
	  const bool   whole        = getFlag ("whole");
		const bool   pseudo       = getFlag ("pseudo");
		const size_t mono_nuc_max = (size_t) arg2uint ("mono_nuc_max");
  #ifndef _MSC_VER
		const bool   large        = getFlag ("large");
  #endif
    const size_t group_count  = (size_t) arg2uint ("group_count");
    const size_t group_size   = (size_t) arg2uint ("group_size");
    const string ext          = getArg ("extension");
                 seqSuffix    = getArg ("suffix");

  
    QC_IMPLY (pseudo, aa);
    QC_IMPLY (mono_nuc_max, ! aa);


    const bool groupP = group_count || group_size;

    QC_ASSERT (! out_dir. empty ());    
  #ifndef _MSC_VER
    QC_IMPLY (large, groupP);
  #endif
    QC_IMPLY (groupP, ext. empty ());


    { // For ~Progress()      
  	  Multifasta fa (in, aa); 
      Group group (group_count, group_size, out_dir);
  	  while (fa. next ())
  	  {
  	    unique_ptr<Seq> seq;
  	    if (aa)
  	    {
  	      auto pep = new Peptide (fa, Peptide::stdAveLen, sparse);
  	      seq. reset (pep);  
  	      if (pseudo && pep->hasInsideStop ())
  	        pep->pseudo = true;
  	    }
  	    else
  	    {
  	      auto dna = new Dna (fa, 128 * 1024, sparse);
  	      dna->qc ();
  	      if (mono_nuc_max)
  	        dna->monoNuc2n (mono_nuc_max + 1);
  	      seq. reset (dna);  
  	    }
  	    seq->qc ();
  	    ASSERT (! seq->name. empty ());
  	    if (isRight (seq->name, splitSuffix))
  	      throw runtime_error ("Suffix " + strQuote (splitSuffix) + " is used in the file " + strQuote (in));
		    if (seq->seq. size () < len_min)
		    	continue;
        string s (seq->name);
        s = findSplit (s);
        if (! whole)
          s = findSplit (s, '|');
        string dir (out_dir);
      #ifndef _MSC_VER
        if (large)
        {
          dir += "/" + to_string (str2hash_class (s, false));
          Dir (dir). create ();
        }
      #endif
        if (groupP)
          group. add (seq. release ());
        else
        {
          seq->appendId (seqSuffix);
          seq->saveFile (dir + "/" + s + (ext. empty () ? "" : ("." + ext)));
        }
  	  }
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



