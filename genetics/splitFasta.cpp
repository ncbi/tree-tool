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
    : Application ("Split a multi-fasta file into a set of fasta files each containing one sequence")
  	{
  	  addPositional ("in", "Input multi-FASTA file");
  	  addPositional ("out_dir", "Output directory");
  	  addFlag ("aa", "Multi-FASTA file contains protein sequences, otherwise DNA sequences");
  	  addFlag ("sparse", "Sparse sequence");
  	  addKey ("len_min", "Minimum sequence length", "0");
    #ifndef _MSC_VER
  	  addFlag ("large", "Create files in subdirectories \"0\" .. \"" + to_string (hash_class_max - 1) + "\" which are the hashes of file names");
  	#endif
  	}



	void body () const final
  {
		const string in      = getArg ("in");
		const string out_dir = getArg ("out_dir");
		const bool   aa      = getFlag ("aa");
		const bool   sparse  = getFlag ("sparse");
		const size_t len_min = (size_t) arg2uint ("len_min");
  #ifndef _MSC_VER
		const bool   large   = getFlag ("large");
  #endif


    {
      // For ~Progress()
  	  Multifasta fa (in, aa, 16 * 1024 * 1024);  // PAR
  	  while (fa. next ())
  	  {
  	    unique_ptr<const Seq> seq;
  	    if (aa)
  	      seq. reset (new Peptide (fa, Peptide::stdAveLen, sparse));  
  	    else
  	      seq. reset (new Dna (fa, 128 * 1024, sparse));  
  	    seq->qc ();
  	    ASSERT (! seq->name. empty ());
		    if (seq->seq. size () < len_min)
		    	continue;
  	    string s (seq->name);
		    s = findSplit (s);
		    s = findSplit (s, '|');
		    string dir (out_dir);
      #ifndef _MSC_VER
		    if (large)
		    {
		      dir += "/" + to_string (str2hash_class (s));
		      createDirectory (dir, true);
		    }
		  #endif
  	    seq->saveFile (dir + "/" + s);
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



