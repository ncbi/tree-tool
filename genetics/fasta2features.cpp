// fasta2features.cpp

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
*   Save position-based features of an alignment
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
    : Application ("Save position-based features of an alignment")
    {
      version = VERSION;
  	  addPositional ("Multi_FASTA", "Multi-FASTA file");
  	  addPositional ("out_dir", "Output directory");
  	  addFlag ("aa", "Multi-FASTA file contains protein sequences, otherwise DNA sequences");
    }


	
	void body () const final
  {
	  const string inFName = getArg ("Multi_FASTA");
	  const string outDir  = getArg ("out_dir");
	  const bool   aa      = getFlag ("aa");


    size_t len = 0;
    { // For ~Progress()      
  	  Multifasta fa (inFName, aa); 
  	  while (fa. next ())
  	  {
  	    const Seq* seq = nullptr;
  	    if (aa)
  	      seq = new Peptide (fa, Peptide::stdAveLen, true);  
  	    else
  	      seq = new Dna (fa, 128 * 1024, true);  
  	    ASSERT (seq);
  	    QC_ASSERT (! seq->name. empty ());
  	    seq->qc ();

  	    if (len)
  	    {
  		    QC_ASSERT (len == seq->seq. size ());
  	    }
  	    else
  	      len = seq->seq. size ();  	      
  	      
  	    OFStream f (outDir + "/" + seq->name);
      	FOR (size_t, i, len)
      	{
    	    const char c = seq->seq [i];
    	    if (! seq->isAmbiguous (c))
    	      f << (i + 1) << ":" << c << endl;
    	  }
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



