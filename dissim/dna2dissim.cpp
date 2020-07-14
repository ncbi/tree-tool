// dna2dissim.cpp

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
*   Align and compute dissimilarity for a pair of DNAs.
*   Print: ref_match: <target start>-<target stop> (human coordinates)
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../genetics/seq.hpp"
using namespace Seq_sp;
#include "align.hpp"



namespace 
{
	

Dna readDna (const string &fName)
{
  LineInput in (fName);
  EXEC_ASSERT (in. nextLine ());
  const Dna dna (in, 1024 * 1024, false);  // PAR
  dna. qc ();
  
  return dna;
}



struct ThisApplication : Application
{
  
  ThisApplication ()
    : Application ("Align and compute dissimilarity for a pair of DNAs.\nPrint: ref_match: <target start>-<target stop> (human coordinates)")
    {
  	  addPositional ("FASTA1", "FASTA file 1 with a DNA sequence (target)");
  	  addPositional ("FASTA2", "FASTA file 1 with a DNA sequence (reference)");
  	  addFlag ("global", "Global alignment, otherwise semiglobal");
  	  addKey ("match_len_min", "Min. match length. Valid for semiglobal alignment", "60");
  	  addKey ("mutation", "file for mutations");
    }


	
	void body () const final
  {
	  const string inFName1   = getArg ("FASTA1");
	  const string inFName2   = getArg ("FASTA2");
	  const bool global       = getFlag ("global");
	        size_t match_len_min = str2<size_t> (getArg ("match_len_min"));
	  const string mutFName   = getArg ("mutation");
    
    if (global)
    	match_len_min = 0;
    else
    	if (match_len_min == 0)
    		throw runtime_error ("match_len_min cannot be 0 for a semiglobal alignment");

    const Dna dna1 = readDna (inFName1);
    const Dna dna2 = readDna (inFName2);

		Align_sp::Align align (dna1, dna2, ! global, match_len_min, /*false,*/ 0);
		
		if (verbose ())
		  align. print (cout);
		
		cout << "Identity = " << align. matches << '/' << align. tr. size () 
		     << " (" << ((double) align. matches / (double) align. tr. size () * 100) << "%)" << endl;
		cout << "Min. edit distance = " << align. getMinEditDistance () << endl;
		cout << endl;

		
		align. setAlignment (dna1. seq, dna2. seq);


 		size_t targetStart = 0;
		size_t targetStop = dna1. seq. size ();
    {		
  		while (align. sparse2 [targetStart] == '-')
  	  {
  	    ASSERT (align. sparse1 [targetStart] != '-');
  		  targetStart++;
  		}
  		
  		size_t i = align. sparse2. size ();
  		while (align. sparse2 [i - 1] == '-')
  	  {
  	    ASSERT (align. sparse1 [i] != '-');
  	    ASSERT (i);
  	    i--;
  	    ASSERT (targetStop);
  		  targetStop--;
  		}
    }
		if (targetStart >= targetStop)
		  throw runtime_error ("No alignment");
		cout << "ref_match: " << targetStart + 1 << ' ' << targetStop << endl;
  		

		if (! mutFName. empty ())
		{
		  OFStream f (mutFName);
		  size_t refPos = 0;
  		size_t mismatchStart = no_index;  
  		size_t refStart = no_index;
  		FFOR (size_t, i, align. sparse1. size ())
  	  {
  		  if (align. sparse1 [i] == align. sparse2 [i])
  		  {
  		    if (mismatchStart != no_index)
  		    {
  		      ASSERT (mismatchStart < i);
  		      const size_t len = i - mismatchStart;
  		      string from (align. sparse2. substr (mismatchStart, len));
  		      string to   (align. sparse1. substr (mismatchStart, len));
  		      replaceStr (from, "-", "");
  		      replaceStr (to,   "-", "");
  		      ASSERT (! (from. empty () && to. empty ()));
  		      if (from. empty ())
  		        from = "INS";
  		      if (to. empty ())
  		        to = "DEL";
  		      f << from << refStart + 1 << to << endl;
  		      mismatchStart = no_index;
  		      refStart = no_index;
  		    }
  		  }
  		  else
  		    if (mismatchStart == no_index)
  		    {
  		      mismatchStart = i;
  		      refStart = refPos;
  		    }
  		  if (align. sparse2 [i] != '-')
  		    refPos++;
  		}
    }


		align. printAlignment (60);  // PAR
  }
};


}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



