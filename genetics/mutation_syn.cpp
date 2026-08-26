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
*   Print synonymous mutations for indels
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


void process (const Mutation &mut,
              string gapS,
              const string &fullS,
              size_t start)
{
  ASSERT (gapS. size () == fullS. size ());
  ASSERT (! contains (fullS, '-'));
  ASSERT (gapS [start] == '-');
  IMPLY (start, gapS [start - 1] != '-');
  
  const size_t gapLen = (size_t) abs (mut. getInsertionSize ());
  ASSERT (gapLen);
  ASSERT (start + gapLen <= gapS. size ());

  if (verbose ())
  {
    cout << gapS << endl;  
    cout << fullS << endl;  
  }

  bool toFix = false;
  for (;;)
  {
    ASSERT (gapS [start] == '-');
    IMPLY (start, gapS [start - 1] != '-');
    ASSERT (gapS [start + gapLen] == fullS [start + gapLen]);
    if (fullS [start] != fullS [start + gapLen])
      break;
    gapS [start] = gapS [start + gapLen];
    gapS [start + gapLen] = '-';
    toFix = true;
    start++;
  }
  if (toFix)
  {
    cout << "Fix:\t" << mut. gene + Mutation::delimiter;
    const int pos = mut. getOffset () + (int) start + 1;  // pos = 0 is impossible ??
    if (mut. isInsertion ())
      cout << str2upper (string (1, fullS [start]), ! mut. prot) 
           << pos 
           << str2upper (fullS. substr (start, gapLen + 1), ! mut. prot);
    else
      cout << str2upper (fullS. substr (start, gapLen), ! mut. prot)
           << pos 
           << "del";
    cout << endl;
  }
  
  for (;;)
  {
    ASSERT (gapS [start] == '-');
    IMPLY (start, gapS [start - 1] != '-');
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
    cout << "Synonym:\t" << mut. gene + Mutation::delimiter;
    const int pos = mut. getOffset () + (int) start + 1;  // pos = 0 is impossible ??
    if (mut. isInsertion ())
    {
      const size_t start_ = start ? start - 1 : start;  
      const int pos_      = start ? pos - 1 : pos;
      cout << str2upper (string (1, fullS [start_]), ! mut. prot) 
           << pos_
           << str2upper (fullS. substr (start_, gapLen + 1), ! mut. prot);
    }
    else
      cout << str2upper (fullS. substr (start, gapLen), ! mut. prot)
           << pos 
           << "del";
    cout << endl;
  }
}

	

	
// ThisApplication

struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Print synonymous mutations for indels")
    {
      version = VERSION;
      addPositional ("seq", "Sequence file");
      addPositional ("pos", "Real mutation position (1-based)");
      addPositional ("mut", "Indel mutation symbol");
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
    QC_ASSERT (! mut. ambig);

    Multifasta mf (seqFName, aa, 0);
    unique_ptr<Seq> seq (makeSeq (mf, false));  
    seq->qc ();    
    
    const string orig (seq->seq);
    mut. apply (seq->seq);
    
    const size_t gapLen = (size_t) abs (mut. getInsertionSize ());
    if (mut. isInsertion ())
    {
      ASSERT (! mut. isDeletion ());
      ASSERT (mut. insertionHasRef);
      ASSERT (mut. allele. size () > 1);
      ASSERT (mut. reference == mut. allele. substr (0, 1));
      process ( mut
              , orig. substr (0, mut. pos_real + 1) + string (gapLen, '-') + orig. substr (mut. pos_real + 1)
              , seq->seq
              , mut. pos_real + 1
              );
    }
    else if (mut. isDeletion ())
    {
      ASSERT (mut. allele. empty ());
      process ( mut
              , orig. substr (0, mut. pos_real) + string (gapLen, '-') + orig. substr (mut. pos_real + gapLen)
              , orig
              , mut. pos_real
              );
    }
    else
      throw runtime_error ("Insertion or deletion mutation is expected");
  }
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



