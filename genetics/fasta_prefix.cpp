// fasta_prefix.cpp

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
*   Replace sequence identifiers with <prefix><number>
*
*/


#undef NDEBUG

#include "../common.hpp"
using namespace Common_sp;
#include "seq.hpp"
using namespace Seq_sp;

#include "../common.inc"



namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Replace sequence identifiers with <prefix><number>")
    {
  	  addPositional ("in", "Input FASTA file");
  	  addPositional ("match", "Output pairs: <old id>\t<new id>");
  	  addPositional ("prefix", "Prefix of new identifiers");
  	  addFlag ("aa", "FASTA is a protein file");
    }


	
	void body () const final
  {
	  const string inFName    = getArg ("in");
	  const string matchFName = getArg ("match");
	  const string prefix     = getArg ("prefix");
	  const bool   aa         = getFlag ("aa");
  

    Multifasta fa (inFName, aa);
    size_t n = 1;
    OFStream matchF (matchFName);
    while (fa. next ())
    {
      unique_ptr<const Seq> seq;
      if (aa)
        seq. reset (new Peptide (fa, 1000, false));
      else
        seq. reset (new Dna (fa, 100000, false));
      seq->qc ();
        
      const string newId = prefix + to_string (n);
      matchF << seq->getId () << '\t' << newId << endl;
      var_cast (seq. get ()) -> name = newId + seq->name. substr (seq->getIdSize ());
      seq->qc ();
      n++;
      
      if (aa)
        seq->asPeptide () -> saveText (cout);
      else
        seq->asDna () -> saveText (cout);
    }
  }
};



}  // namespace



int main(int argc, 
         const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



