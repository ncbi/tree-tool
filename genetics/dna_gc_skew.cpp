// dna_gc_skew.cpp

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
*   DNA GC skew (G - C) / (G + C) distribution
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
	  : Application ("DNA GC skew (G - C) / (G + C) distribution")
		{
      version = VERSION;
		  addPositional ("in", "Input FASTA file with a DNA sequence");
		}



	void body () const final
  {
		const string inFName = getArg ("in");

    LineInput in (inFName);
    EXEC_ASSERT (in. nextLine ());
    Dna dna (in, 1024 * 1024, false);  // PAR
    ASSERT (! dna. name. empty ());
    dna. qc ();
    
    constexpr size_t window = 100000;  // PAR
    
    dna. seq += dna. seq. substr (0, window);
    
    // CpG ??
    
    size_t cg = 0;
    int g = 0; 
    FFOR (size_t, i, dna. seq. size ())
    {
      if (i % window == 0)
      {
        if (cg)
          cout << i - window << '\t' << (double) g / (double) cg << endl;
        cg = 0;
        g = 0;
      }
      const char nuc = dna. seq [i];
      if (   nuc == 'c' 
          || nuc == 'g'
         )
      {
        cg++;
        if (nuc == 'g')
          g++;
        else
          g--;
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



