// dna_cut.cpp

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
*   Update a DNA sequence by excising a segment out of it
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
    : Application ("Cut a segment out of a DNA sequence, make it to be in a positive strand")
    {
      version = VERSION;
  	  addPositional ("in", "DNA FASTA file with one sequence");
  	  addPositional ("start", "Start position of a segment, 1-based");
  	  addPositional ("stop", "Stop position of a segment, >= start position");
  	  addKey ("flank", "length of flanking sequence", "0");
  	  addFlag ("excise", "Excise the segment, otherwise leave the segment");
  	  addKey ("strand", "Strand (0 or - / 1 or +)", "1");
    }


	
	void body () const final
  {
	  const string inFName = getArg ("in");
	        size_t start   = str2<size_t> (getArg ("start"));
	        size_t stop    = str2<size_t> (getArg ("stop"));
	  const size_t flank   = str2<size_t> (getArg ("flank"));
	  const bool excise    = getFlag ("excise");
	  const string strandS = getArg ("strand");
	  
	  QC_ASSERT (start >= 1);
	  QC_ASSERT (stop >= start);
	  start--;  // To make 0-based
	  
	  bool strand = false;
	  if (strandS == "0" || strandS == "-")
	    ;
	  else if (strandS == "1" || strandS == "+")
	    strand = true;
	  else
	    throw runtime_error ("Wrong \"strand\" parameter");

	  
	  Dna* dna = nullptr;  // Not delete'd
	  {
  	  constexpr size_t size = 1024 * 1024;  // PAR
  	  LineInput li (inFName, size);
  	  QC_ASSERT (li. nextLine ());
      dna = new Dna (li, size, false);  
    }
    ASSERT (dna);
    dna->qc ();    
    const size_t len = dna->seq. size ();
    QC_ASSERT (stop <= len);
    
    if (start > flank)
      start -= flank;
    else
      start = 0;
      
    if (stop + flank < len)
      stop += flank;
    else
      stop = len;
      
    dna->seq = (excise 
                  ? dna->seq. substr (0, start) + dna->seq. substr (stop)
                  : dna->seq. substr (start, stop - start)
               );
    if (! strand)
      dna->reverse ();
    dna->saveText (cout);
  }
};


}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



