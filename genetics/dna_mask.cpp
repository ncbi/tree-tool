// dna_mask.cpp

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
*   replace segments of DNA by Ns
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
  
  
struct Mask 
{
  size_t from;
  size_t to;
};

  
  
struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Replace segments of DNA by Ns")
  	{
      version = VERSION;
  	  addPositional ("in", "Input multi-FASTA file");
  	  addPositional ("mask_file", "File with lines: contig id\\tfrom\\to (1-based)");
  	}



	void body () const final
  {
		const string in         = getArg ("in");
		const string maskFName  = getArg ("mask_file");
	  

    map<string/*contig*/,Vector<Mask>> contig2masks;
    {
      LineInput f (maskFName);
      while (f. nextLine ())
      {
        istringstream iss (f. line);
        string contig;
        size_t from;
        size_t to {no_index};
        iss >> contig >> from >> to;
        QC_ASSERT (to != no_index);
        if (from > to)
          swap (from, to);
        QC_ASSERT (from);
        from--;
        contig2masks [contig] << Mask {from, to};
      }
    }


    { // For ~Progress()      
  	  Multifasta fa (in, false); 
  	  while (fa. next ())
  	  {
	      Dna dna (fa, 16 * 1024 * 1024, false);
	      dna. qc ();
	      for (const Mask& m : contig2masks [dna. getId ()])
	      {
	        ASSERT (m. from < m. to);
	        ASSERT (m. to <= dna. seq. size ());
	        FFOR_START (size_t, i, m. from, m. to)
	          dna. seq [i] = 'n';
	      }
	      dna. qc ();
	      dna. saveText (cout);
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



