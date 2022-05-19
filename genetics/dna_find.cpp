// dna_find.cpp

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
*   Find DNA subsequences
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
	
	

void search (const Dna &hay,
             const Dna &needle,
             bool strand)
{
  size_t pos = 0;
  for (;;)
  {
    pos = hay. seq. find (needle. seq, pos);
    if (pos == string::npos)
      break;
    cout << hay. getId () << '\t' << pos + 1 << '\t' << pos + needle. seq. size () << '\t' << strand << endl;
    pos++;
  }  	    	    
}



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print 1-based coordinates and strand of a DNA needle in a DNA hay")
	  {
      version = VERSION;
		  addPositional ("hay",  "Multi-FASTA DNA file");
		  addPositional ("needle", "FASTA file with a target DNA sequence to be found");
	  }

  
  
	void body () const final
  {
	  const string inFName     = getArg ("hay");
	  const string targetFName = getArg ("needle");


    const FastaDna target (targetFName, 1000, false);  // PAR
    target. qc ();
    
    Dna targetRev (target);
    targetRev. reverse ();
    targetRev. qc();
    

	  Multifasta fa (inFName, false, 128 * 1024 * 1024, 0);  // PAR 
	  while (fa. next ())
	  {
	    const Dna dna (fa, 1000/*PAR*/, false);
	    dna. qc ();
	    search (dna, target, true);
	    search (dna, targetRev, false);
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



