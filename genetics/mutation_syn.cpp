// mutation_syn.cpp

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
*   Print synonymous mutations
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


void process (const string &gene,
              bool insertion,
              string gapS,
              const string &fullS,
              size_t start,
              size_t gapLen)
{
  ASSERT (gapS. size () == fullS. size ());
  ASSERT (start + gapLen <= gapS. size ());
  ASSERT (! contains (fullS, '-'));
  ASSERT (gapLen);
  ASSERT (gapS [start] == '-');
  IMPLY (start, gapS [start - 1] != '-');

  if (verbose ())
  {
    cout << gapS << endl;  
    cout << fullS << endl;  
  }

  bool fixed = false;
  for (;;)
  {
    ASSERT (gapS [start] == '-');
    start--;
    if (start == (size_t) -1)
      break;
    ASSERT (gapS [start] != '-');
    ASSERT (gapS [start] == fullS [start]);
    if (fullS [start] != fullS [start + gapLen])
      break;
    ASSERT (gapS [start + gapLen] == '-');
    gapS [start + gapLen] = gapS [start];
    gapS [start] = '-';
    fixed = true;
  }
  if (fixed)
  {
    cout << "Fix:\t" << gene + Mutation::delimiter;
    if (insertion)
    {
      const size_t start_ = start == (size_t) -1 ? 0 : start;
      cout << fullS [start_] + to_string (start_ + 1) + fullS. substr (start_, gapLen + 1);
    }
    else
      cout << fullS. substr (start, gapLen) << to_string (start + 1) << "del";
    cout << endl;
  }

  for (;;)
  {
    start++;
    ASSERT (gapS [start] == '-');
    IMPLY (start, gapS [start - 1] != '-');
    ASSERT (gapS [start + gapLen] == fullS [start + gapLen]);
    if (fullS [start] != fullS [start + gapLen])
      break;
    gapS [start] = gapS [start + gapLen];
    gapS [start + gapLen] = '-';
    cout << "Synonym:\t" << gene + Mutation::delimiter;
    if (insertion)
      cout << fullS [start] + to_string (start + 1) + fullS. substr (start, gapLen + 1);
    else
      cout << fullS. substr (start, gapLen) << to_string (start + 1) << "del";
    cout << endl;
  }
}

	

	
// ThisApplication

struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Print synonymous mutations")
    {
      version = VERSION;
      addPositional ("seq", "Sequence file");
      addPositional ("pos", "Real mutation position (1-based)");
      addPositional ("mut", "Mutation symbol");
  	  addFlag ("aa", "Sequence is a protein, otherwise DNA");
    }



  void body () const final
  {
    const string seqFName = getArg ("seq");
          size_t pos_real = (size_t) arg2uint ("pos");  
    const string mutS     = getArg ("mut");  
		const bool   aa       = getFlag ("aa");


    const Mutation mut (aa, pos_real, mutS, mutS);
    mut. qc ();

    Multifasta mf (seqFName, aa, 0);
    unique_ptr<Seq> seq (makeSeq (mf, false));  
    seq->qc ();         
    
    const string orig (seq->seq);
    mut. apply (seq->seq);
    
    if (mut. isInsertion ())
    {
      ASSERT (! mut. isDeletion ());
      ASSERT (mut. insertionHasRef);
      ASSERT (mut. reference. size () == 1);
      ASSERT (mut. allele. size () > 1);
      const size_t diff = mut. allele. size () - 1;
      process ( mut. gene
              , true
              , orig. substr (0, mut. pos_real + 1) + string (diff, '-') + orig. substr (mut. pos_real + 1)
              , seq->seq
              , mut. pos_real + 1
              , diff
              );
    }
    else if (mut. isDeletion ())
    {
      ASSERT (mut. allele. empty ());
      const size_t diff = mut. reference. size ();
      process ( mut. gene
              , false
              , orig. substr (0, mut. pos_real) + string (diff, '-') + orig. substr (mut. pos_real + diff)
              , orig
              , mut. pos_real
              , diff
              );
    }
    else
      throw runtime_error ("Insertion or deletion is expected");
  }
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



