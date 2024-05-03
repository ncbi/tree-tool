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
*   Split a multi-fasta file into a set of fasta files each containing one sequence
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
  
  
struct Group
{
  const size_t size_max;
  const string dir;
  VectorOwn<const Seq> seqs;
  size_t group {0};
  
  
  Group (size_t size_max_arg,
         const string &dir_arg)
    : size_max (size_max_arg)
    , dir (dir_arg)
    { ASSERT (size_max >= 1); 
      ASSERT (! dir. empty ());
    }
 ~Group ()
    { save (); }
  
  
  void add (const Seq* seq)
    { ASSERT (size_max > 1);
      ASSERT (seq);
      seqs << seq;
      ASSERT (seqs. size () <= size_max);
      if (seqs. size () == size_max)
        save ();
    }
  void save ()
    { if (seqs. empty ())
        return;
      group++;
      OFStream f (dir + "/" + to_string (group));
      for (const Seq* seq : seqs)
        seq->saveText (f);
      seqs. deleteData ();
    }
};
  
  
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Split a multi-fasta file into a set of fasta files each containing one sequence")
  	{
      version = VERSION;
  	  addPositional ("in", "Input multi-FASTA file");
  	  addPositional ("out_dir", "Output directory");
  	  addFlag ("aa", "Multi-FASTA file contains protein sequences, otherwise DNA sequences");
  	  addFlag ("sparse", "Sparse sequence");
  	  addFlag ("pseudo", "Pseudo-proteins are allowed");
  	  addKey ("len_min", "Minimum sequence length", "0");
		  addFlag ("whole", "Sequence identifiers are whole strings which should not be split by '|'");
    #ifndef _MSC_VER
  	  addFlag ("large", "Create files in subdirectories \"0\" .. \"" + to_string (hash_class_max - 1) + "\" which are the hashes of file names");
  	#endif
  	  addKey ("group", "Group by <group> number of sequences. Group names are sequential numbers. 1 - no grouping.", "1");
  	  addKey ("extension", "Add file extension (not compatible with -group)", "");
  	}



	void body () const final
  {
		const string in          = getArg ("in");
		const string out_dir     = getArg ("out_dir");
		const bool   aa          = getFlag ("aa");
		const bool   sparse      = getFlag ("sparse");
		const bool   pseudo      = getFlag ("pseudo");
		const size_t len_min     = (size_t) arg2uint ("len_min");
	  const bool   whole       = getFlag ("whole");
  #ifndef _MSC_VER
		const bool   large       = getFlag ("large");
  #endif
    const size_t group_size = (size_t) arg2uint ("group");
    const string ext        = getArg ("extension");

    QC_ASSERT (! out_dir. empty ());    
    QC_ASSERT (group_size >= 1);
  #ifndef _MSC_VER
    QC_IMPLY (large, group_size == 1);
  #endif
    QC_IMPLY (group_size > 1, ext. empty ());


    { // For ~Progress()      
  	  Multifasta fa (in, aa); 
      Group group (group_size, out_dir);
  	  while (fa. next ())
  	  {
  	    unique_ptr<const Seq> seq;
  	    if (aa)
  	    {
  	      auto pep = new Peptide (fa, Peptide::stdAveLen, sparse);
  	      seq. reset (pep);  
  	      if (pseudo && pep->hasInsideStop ())
  	        pep->pseudo = true;
  	    }
  	    else
  	      seq. reset (new Dna (fa, 128 * 1024, sparse));  
  	    seq->qc ();
  	    ASSERT (! seq->name. empty ());
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
          dir += "/" + to_string (str2hash_class (s));
          Dir (dir). create ();
        }
      #endif
        if (group_size == 1)
          seq->saveFile (dir + "/" + s + (ext. empty () ? "" : ("." + ext)));
        else
          group. add (seq. release ());
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



