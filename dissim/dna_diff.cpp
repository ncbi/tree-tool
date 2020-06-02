// dna_diff.cpp

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
*   Print the number of nucleotide differences between two DNA sequences, case-insensitive
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../genetics/seq.hpp"
using namespace Seq_sp;




namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print the number of nucleotide differences between two DNA sequences, case-insensitive")
	  {
	  	addPositional ("in", "Input file with two DNA sequences");
	  }



  void body () const final
  {
    const string inFName = getArg ("in");
    
    
    string seq1, seq2;
    {
      LineInput f (inFName);
      EXEC_ASSERT (f. nextLine ());
      seq1 = f. line;
      EXEC_ASSERT (f. nextLine ());
      seq2 = f. line;
    }
    trim (seq1);
    trim (seq2);
    strLower (seq1);
    strLower (seq2);
    QC_ASSERT (! seq1. empty ());
    QC_ASSERT (! seq2. empty ());
    QC_ASSERT (seq1. size () == seq2. size ());
    
    const Dna dna1 ("1", seq1, true);    
    const Dna dna2 ("2", seq2, true);    
    dna1. qc ();
    dna2. qc ();
    ASSERT (dna1. seq. size () == dna2. seq. size ());
    
    size_t n = 0;
    FFOR (size_t, i, dna1. seq. size ())
      if (! nucleotideMatch ( dna1. seq [i]
                            , dna2. seq [i]
                            )
         )
        n++;
        
    cout << n << endl;
  }  
};


}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



