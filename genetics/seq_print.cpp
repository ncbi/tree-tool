// seq_print.cpp

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
*   Print a sequence with spacing
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


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print a sequence with spacing")
	  {
      version = VERSION;
	  	addPositional ("in", "Input DNA multi-FASTA file");
  	  addFlag ("aa", "Sequence is a protein, otherwise DNA");
	  }



  void body () const final
  {
    const string inFName = getArg ("in");
    const bool   aa      = getFlag ("aa");


	  Multifasta fa (inFName, aa, 0); 
	  bool first = true;
  	while (fa. next ())
  	{
      unique_ptr<const Seq> seq;
      if (aa)
        seq. reset (new Peptide (fa, Peptide::stdAveLen, false));  
      else
        seq. reset (new Dna (fa, 128 * 1024, false));  
      seq->qc ();
      
      if (! first)
        cout << endl;
      cout << seq->name << endl;
      for (size_t i = 0; i < seq->seq. size (); i += 10)
      {
        if (i && ! (i % 60))
          cout << ' ' << i << endl;
        cout << seq->seq. substr (i, 10) << ' ';
      }
      cout << ' ' << seq->seq. size () << endl;
      
      first = false;
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



