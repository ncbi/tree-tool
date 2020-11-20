// extractFastaProt.cpp

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
*   Print unique protein sequences
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
    : Application ("Print unique protein sequences")
	  {
		  addPositional ("in", "Input FASTA file with proteins");
	  }

  
  
	void body () const final
  {
	  const string inFName = getArg ("in");

  
    unordered_map<size_t/*hash*/,Peptide> hash2pep;  hash2pep. rehash (100000);  // PAR
		{
		  Multifasta fa (inFName, true);
		  while (fa. next ())
		  {
		    Peptide pep (fa, 1000/*PAR*/, false);
		    pep. qc ();
		    const size_t h = str_hash (pep. seq);
		    if (! contains (hash2pep, h))
		      hash2pep [h] = move (pep);
		  }
		}


    for (const auto& it : hash2pep)
		  it. second. saveText (cout);
  }
};


}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



