// uniqSeq.cpp

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
*   Print unique sequences
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
    : Application ("Print unique sequences")
	  {
		  addPositional ("in", "Input FASTA file");
		  addFlag ("aa", "Sequences are proteins, otherwise DNA");
		  addKey ("pair", "Replacement id pairs: redundant id <tab> main id");
	  }

  
  
	void body () const final
  {
	  const string inFName   = getArg ("in");
	  const bool   aa        = getFlag ("aa");
	  const string pairFName = getArg ("pair");


    OFStream* pairF = nullptr;
    if (! pairFName. empty ())  
      pairF = new OFStream (pairFName);
    unordered_map<size_t/*hash*/,Seq*> hash2seq;  hash2seq. rehash (1000000);  // PAR
      // Not delete'd
		{
		  Multifasta fa (inFName, aa);
		  while (fa. next ())
		  {
		    unique_ptr<Seq> seq;
		    if (aa)
		      seq. reset (new Peptide (fa, 1000/*PAR*/, false));
		    else
		      seq. reset (new Dna (fa, 1000/*PAR*/, false));  
		    ASSERT (seq. get ());
		    seq->qc ();
		    
		    const size_t h = str_hash (seq->seq);
		    if (const Seq* const* other = findPtr (hash2seq, h))
		    {
		      if ((*other)->getId () > seq->getId ())
		      {
		        if (pairF)
		          *pairF << (*other)->getId () << '\t' << seq->getId () << endl;
		        var_cast (*other) -> name = seq->name;
		      }
		      else
		      {
		        if (pairF)
		          *pairF << seq->getId () << '\t' << (*other)->getId () << endl;
		      }
		    }
		    else
		      hash2seq [h] = seq. release ();
		  }
		}


    for (const auto& it : hash2seq)
    {
      const Seq* seq = it. second;
      ASSERT (seq);
		  seq->saveText (cout);
      if (pairF)
        *pairF << seq->getId () << '\t' << seq->getId () << endl;
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



