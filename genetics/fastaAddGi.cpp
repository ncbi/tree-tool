// fastaAddGi.cpp

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
*   "Create a FASTA file <fasta>.out where each identifier is prepended by a new gi starting from <min_gi>
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


//const size_t len = 20;  // PAR



struct ThisApplication : Application
{
	ThisApplication ()
	  : Application ("Create a FASTA file <fasta>.out where each identifier is prepended by a new gi starting from <min_gi>")
           //"Create the map file <fasta>.map"; 
	  {
      version = VERSION;
		  addPositional ("fasta", "Input FASTA file");
		  addPositional ("min_gi", "Min. gi");
		  addFlag ("aa", "Amino acid sequence, otherwise nucleotide");
	  }

	

	void body () const final
  {
	  const string inFName =             getArg ("fasta");
	  const long minGi     = str2<long> (getArg ("min_gi"));
	  const bool aa = getFlag ("aa");
	  ASSERT (minGi > 0);
    
    
    OFStream fOut ("", inFName, "out");
  //OFStream fMap ("", inFName, "map");
    long /*CSeq_id::TGi*/ gi = minGi;
	  Multifasta faIn (inFName, aa);
	  while (faIn. next ())
	  {
	  	unique_ptr<Seq> seq (makeSeq (faIn, false));

	    ASSERT (! strBlank (seq->name));
	    ASSERT (! contains (seq->getId (), "gi|"));
	    		    
	  //fMap << gi << " " << seq->getId () << endl;
	    
	    seq->name = "gi|" + toString (gi) + "|" + seq->name;
	    seq->printCase (fOut, aa);

	    ASSERT (gi < numeric_limits<int>::max());
	    gi++;
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



