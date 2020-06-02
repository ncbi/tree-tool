// dna_trim.cpp

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
*   Update a DNA sequence by trimming head and tail
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
    : Application ("Update a DNA sequence by trimming head and tail")
    {
  	  addPositional ("in", "DNA FASTA file with a DNA sequence");
  	  addPositional ("head", "Head length");
  	  addPositional ("tail", "Tail length");
    }


	
	void body () const final
  {
	  const string inFName = getArg ("in");
	  const size_t head    = str2<size_t> (getArg ("head"));
	  const size_t tail    = str2<size_t> (getArg ("tail"));

    if (! head && ! tail)
      return;	  

	  
	  Dna* dna = nullptr;  // Not delete'd
	  {
  	  constexpr size_t size = 1024 * 1024;  // PAR
  	  LineInput li (inFName, size);
  	  QC_ASSERT (li. nextLine ());
      dna = new Dna (li, size, false);  
    }
    ASSERT (dna);
    dna->qc ();    
    
    if (head + tail >= dna->seq. size ())
      throw runtime_error ("Too short DNA");
    
    dna->seq. erase (dna->seq. size () - tail);
    dna->seq. erase (0, head);
    dna->qc ();
    
    dna->saveFile (inFName);
  }
};


}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



