// fasta2GC.cpp

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
*   Print a Data Master file with GC% of DNA sequences
*
*/


#undef NDEBUG

#include "../common.hpp"
using namespace Common_sp;
#include "seq.hpp"
using namespace Seq_sp;
#include "../dm/dataset.hpp"
using namespace DM_sp;
#include "../version.inc"

#include "../common.inc"



namespace 
{
	
	
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print " + dmSuff + "-file with GC% of DNA sequences")
    {
      version = VERSION;
  	  addPositional ("Multi_FASTA", "Multi-FASTA file with DNA sequences");
    }


	
	void body () const final
  {
	  const string inFName = getArg ("Multi_FASTA");


    Dataset ds;  
    auto gcAttr = new RealAttr1 ("GC", ds, 6);
    const size_t len_min = 10000;  // PAR
    {
		  Multifasta faIn (inFName, false);
		  while (faIn. next ())
		  {
		    const Dna dna (faIn, len_min, false);
		    dna. qc ();
		    if (dna. seq. size () < len_min)  
		      continue;
		    
		    const map<char,size_t> chars (dna. getCharCount ());
		    const size_t a = chars. at ('a');
		    const size_t c = chars. at ('c');
		    const size_t g = chars. at ('g');
		    const size_t t = chars. at ('t');
		    const size_t all = a + c + g + t;
		    if (all < len_min)
		  //if ((Real) all / (Real) dna. seq. size () < 0.99)  // PAR
		      continue;
		      
		    const size_t objNum = ds. appendObj (dna. name);
		    const_cast <Obj*> (ds. objs [objNum]) -> mult = (Real) all;
		    (*gcAttr) [objNum] = (Real) (c + g) / (Real) all;
		  }
		}
		
		ds. saveText (cout);
  }
};


}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



