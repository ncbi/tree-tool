// fragmentize.cpp

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
*   Split a DNA multi-fasta file into a set of contigs
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
  
  
constexpr double x_frac_max = 0.25;  // PAR

  
  
struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Split a DNA multi-fasta file a set of contigs.\n\
Contigs with ambiguities > " + to_string (x_frac_max * 100) + " % are skipped")
  	{
      version = VERSION;
  	  addPositional ("in", "Input multi-FASTA file");
  	  addPositional ("contig_len", "Contig length");
  	}



	void body () const final
  {
		const string in         = getArg ("in");
		const size_t contig_len = (size_t) arg2uint ("contig_len");
	  
	  QC_ASSERT (contig_len);


    { // For ~Progress()      
  	  Multifasta fa (in, false); 
  	  while (fa. next ())
  	  {
	      const Dna dna (fa, 16 * 1024 * 1024, false);
	      dna. qc ();          
        size_t start = 0; 
        bool stopP = false;
        while (! stopP)
        {
          ASSERT (start <= dna. seq. size ());
          size_t stop = start + contig_len; 
          if (minimize (stop, dna. seq. size ()))
            stopP = true;
          if (start < stop)
          {
            const size_t start_human = start + 1;
            const size_t stop_human = stop;
            const string name (dna. getId () + ":" + to_string (start_human) + "-" + to_string (stop_human) + " " + dna. getDescription (false));
            Dna frag (name, dna. seq. substr (start, stop - start), false);    
            frag. trimAmbiguous ();
            frag. qc ();
            if (   frag. seq. size () >= contig_len / 2  // PAR
                && (double) frag. getXs () <= (double) frag. seq. size () * x_frac_max
               ) 
              frag. saveText (cout);
            start = stop;
          }
          else
          {             
            ASSERT (stopP);
          }
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



