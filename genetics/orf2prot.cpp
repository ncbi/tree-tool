// orf2prot.cpp

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
*   Convert ORFs into proteins
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
    : Application ("Convert ORFs into proteins")
    {
      version = VERSION;
  	  addPositional ("Multi_FASTA", "Multi-FASTA file with ORFs");
  	  addFlag ("no_stop", "It is allowed to miss the stop codon");
  	  addKey ("error", "File with errors", "");
    }


	
	void body () const final
  {
	  const string inFName  = getArg ("Multi_FASTA");
	  const bool   no_stop  = getFlag ("no_stop");
	  const string errFName = getArg ("error");


    const Gencode gencode = 11;  // PAR

    unique_ptr<OFStream> fOut;
    if (! errFName. empty ())
      fOut. reset (new OFStream (errFName));
    
	  Multifasta faIn (inFName, false);
	  while (faIn. next ())
	  {
	    const Dna dna (faIn, 1024 * 1024, false);
	    dna. qc ();
	    size_t translationStart = 0;
	    Peptide pep (dna. makePeptide (1, gencode, true, true, translationStart));
	    ASSERT (translationStart == 0);
	    pep. name = dna. name;
	  //trimSuffix (pep. name, ".fr1");
	    pep. saveText (cout);
	    pep. pseudo = true;  // For qc()
	    pep. qc ();
	    if (fOut. get ())
	    {
	      string err;
	      if (pep. seq. empty ())
	        err += "\tno sequence";
	      if (! trimSuffix (pep. seq, "*"))
	        if (! no_stop)
	          err += "\tno stop codon";
	      if (pep. hasInsideStop ())
	        err += "\tin-frame stop codon";
	      if (pep. seq [0] != 'm')
	        err += "\tno start codon";
	      if (! err. empty ())
          *fOut << pep. name << err << endl;
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



