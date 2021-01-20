// interSeq.cpp

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
*   Find identical (intersection) sequences in 2 FASTA files
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
	
	
	
struct UniqSeq
{
  const Seq* seq {nullptr};
    // !nullptr
    // Not delete'd
  StringVector ids;
};
	
	

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Find identical (intersection) sequences in 2 FASTA files")
	  {
		  addPositional ("in1", "Input FASTA file 1");
		  addPositional ("in2", "Input FASTA file 2");
		  addFlag ("aa", "Sequences are proteins, otherwise DNA");
	  }

  
  
	void body () const final
  {
	  const string inFName1  = getArg ("in1");
	  const string inFName2  = getArg ("in2");
	  const bool   aa        = getFlag ("aa");


    unordered_map<size_t/*hash*/,string> hash2id;  hash2id. rehash (1000000);  // PAR
		{
		  Multifasta fa (inFName1, aa);
		  while (fa. next ())
		  {
		    unique_ptr<Seq> seq;
		    if (aa)
		      seq. reset (new Peptide (fa, 1000/*PAR*/, false));
		    else
		      seq. reset (new Dna (fa, 1000/*PAR*/, false));  
		    ASSERT (seq. get ());
		    seq->qc ();
		    QC_ASSERT (! seq->getId (). empty ());
		    
		    const size_t h = str_hash (seq->seq);
		    hash2id [h] = seq->name;
		  }
		}

		{
		  Multifasta fa (inFName2, aa);
		  while (fa. next ())
		  {
		    unique_ptr<Seq> seq;
		    if (aa)
		      seq. reset (new Peptide (fa, 1000/*PAR*/, false));
		    else
		      seq. reset (new Dna (fa, 1000/*PAR*/, false));  
		    ASSERT (seq. get ());
		    seq->qc ();
		    QC_ASSERT (! seq->getId (). empty ());
		    
		    const size_t h = str_hash (seq->seq);
		    if (const string* s = findPtr (hash2id, h))
		      cout << seq->name << '\t' << *s << endl;
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



