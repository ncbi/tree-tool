// align2html.cpp

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
    : Application ("Print aligned sequences as HTML")
	  {
      version = VERSION;
	  	addPositional ("in", "Input DNA multi-FASTA file, alignment");
  	  addFlag ("aa", "Sequence is a protein, otherwise DNA");
	  }



  void body () const final
  {
    const string inFName = getArg ("in");
    const bool   aa      = getFlag ("aa");


    VectorOwn<Seq> seqs;
    size_t len = 0;
    size_t name_max = 0;
    {
  	  Multifasta fa (inFName, aa, 0); 
    	while (fa. next ())
    	{
    	  Seq* seq = nullptr;
        if (aa)
          seq = new Peptide (fa, Peptide::stdAveLen, true);  
        else
          seq = new Dna (fa, 128 * 1024, true);  
        ASSERT (seq);
        seq->qc ();
        seqs << seq;
        if (len)
        {
          QC_ASSERT (seq->seq. size () == len);
        }
        else
          len = seq->seq. size ();
        seq->seq. reserve (len * 2);  // PAR
        maximize (name_max, seq->name. size ());
      }
    }
    
    if (seqs. empty ())
      return;

    string consensus (seqs. front () -> seq);
    FFOR (size_t, i, consensus. size ())
    {
      map<char,size_t> char2n;
      for (const Seq* seq : seqs)
        char2n [seq->seq [i]] ++;
      char c = '\0';
      size_t n_max = 0;
      for (const auto& it : char2n)
        if (maximize (n_max, it. second))
          c = it. first;
      ASSERT (c);
      consensus [i] = c;
    }
    
    cout << "<pre>";
    for (const Seq* seq : seqs)
    {
      cout << pad (seq->name, name_max, efalse) << "  ";
      FOR_REV (size_t, i, consensus. size ())
        if (seq->seq [i] != consensus [i])
        {
          const string by (string ("<font color=red><b>") + seq->seq [i] + "</b></font>");
          var_cast (seq) -> seq = seq->seq. substr (0, i) + by + seq->seq. substr (i + 1);
        }
      cout << seq->seq << endl;
    }
    cout << "</pre>";
  }  
};


}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



