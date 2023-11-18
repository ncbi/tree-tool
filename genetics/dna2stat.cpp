// dna2stat.cpp

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
*   Compute statistics for a DNA multi-FASTA file: length, ambiguities fraction, GC fraction, octamers fraction
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
    : Application ("Compute statistics for a DNA multi-FASTA file: length, ambiguities fraction, GC fraction, octamers fraction")
	  {
      version = VERSION;
	  	addPositional ("in", "Input DNA multi-FASTA file");
	  	addFlag ("acgt", "Add A, C, G and T fractions");
	  }



  void body () const final
  {
    const string inFName = getArg ("in");
    const bool acgtP     = getFlag ("acgt");


    size_t len = 0;
	  size_t ambiguities = 0;
	  size_t gc = 0;
    constexpr size_t repeatLen_max = 15;  // PAR
    Vector<size_t> acgt (128, 0);
    Vector<size_t> repeat (repeatLen_max, 0);
	  {
		  Multifasta faIn (inFName, false);
		  while (faIn. next ())
		  {
		    Dna dna (faIn, 1e6/*PAR*/, false);
	    	strLower (dna. seq);
			  dna. qc ();
	
			  len += dna. seq. size ();
			  
			  char prev = ' ';
			  size_t repeatLen = 1;
			  for (const char c : dna. seq)
			    if (charInSet (c, "acgt"))
			    {
			      if (   c == 'g'
			          || c == 'c'
			         )
			        gc++;
			      if (c == prev)
			        repeatLen++;
			      else
			      {
			        if (repeatLen < repeatLen_max)
			          repeat [repeatLen] ++;
			        repeatLen = 1;
			      }
			      acgt [(size_t) c] ++;
			      prev = c;
			    }
			    else
			    {
			      ambiguities++;
			      prev = ' ';
			    }
	      if (repeatLen < repeatLen_max)
	        repeat [repeatLen] ++;
	    }
      ASSERT (! repeat [0]);
		}
		
		const size_t len_pure = len - ambiguities;

	  ONumber on (cout, 6, false);  // PAR
    cout << inFName
	       << '\t' << len 
	       << '\t' << (double) ambiguities / (double) len 
	       << '\t' << (double) gc / (double) len_pure;
	#if 0
	  FOR_START (size_t, i, 2, repeatLen_max)
	    cout << '\t' << (double) (repeat [i] * i) / (double) len_pure;
	#else
	    cout << '\t' << (double) (repeat [8] * 8) / (double) len_pure;
	#endif
	  if (acgtP)
	    cout << '\t' << (double) acgt ['a'] / (double) len_pure
	         << '\t' << (double) acgt ['c'] / (double) len_pure
	         << '\t' << (double) acgt ['g'] / (double) len_pure
	         << '\t' << (double) acgt ['t'] / (double) len_pure;
	  cout << endl;
  }  
};


}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



