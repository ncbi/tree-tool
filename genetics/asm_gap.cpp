// asm_gap.cpp

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
*   Compute the distribution of gaps in a DNA multi-FASTA file; print: # contigs, N50; save contigs split by gaps
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "seq.hpp"
using namespace Seq_sp;
#include "../dm/dataset.hpp"
using namespace DM_sp;
#include "../version.inc"




namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Compute the distribution of gaps in a DNA multi-FASTA file; print: # contigs, N50; save contigs split by gaps")
	  {
      version = VERSION;
	  	addPositional ("in", "Input DNA multi-FASTA file");
	  	addPositional ("gap_min", "Minimum gap length");
	  	addKey ("gap_distr", "Output file with gap distribution");
	  	addKey ("out", "Output multi-FASTA file split into separate contigs using gaps as separators");
	  }



  void body () const final
  {
    const string inFName    = getArg ("in");
    const size_t gap_min    = (size_t) arg2uint ("gap_min");
    const string distrFName = getArg ("gap_distr");
    const string outFName   = getArg ("out");


    unique_ptr<OFStream> outF;
    if (! outFName. empty ())
      outF. reset (new OFStream (outFName));

    Vector<size_t> contigs;  contigs. reserve (1000);  // PAR
    Dataset ds;  
    auto gapAttr = new IntAttr1 ("gap_len", ds);
	  {
		  Multifasta faIn (inFName, false);
		  while (faIn. next ())
		  {
		    Dna dna (faIn, 1e6/*PAR*/, false);
	    	strLower (dna. seq);
			  dna. qc ();
	
			  size_t gapLen = 0;
			  size_t start = 0;
			  FFOR (size_t, i, dna. seq. size ())
			  {
			    if (charInSet (dna. seq [i], "acgt"))
			    {
			      if (gapLen)
			      {			        
      		    const size_t objNum = ds. appendObj ();
      		    (*gapAttr) [objNum] = (int) gapLen;
      		    if (gapLen >= gap_min)
      		    {
      		      ASSERT (i >= start + gapLen);
      		      const size_t contigLen = i - start - gapLen;
    		        contigs << contigLen;
    		        if (outF. get ())
  		          {
  		            const Dna contig ( dna. getId () + ":" + to_string (start + 1) + "-" + to_string (i)
  		                             , dna. seq. substr (start, contigLen)
  		                             , false
  		                             );
  		            contig. saveText (*outF);
  		          }
    		        start = i;
      		    }
			      }
			      gapLen = 0;
			    }
			    else
			      gapLen++;
			  }
	      if (gapLen)
	      {			        
  		    const size_t objNum = ds. appendObj ();
  		    (*gapAttr) [objNum] = (int) gapLen;
  		    if (gapLen >= gap_min)
  		    {
  		      const size_t dnaSize = dna. seq. size ();
  		      ASSERT (dnaSize >= start + gapLen);		        
  		      const size_t contigLen = dnaSize - start - gapLen;
		        contigs << contigLen;
		        if (outF. get ())
	          {
	            const Dna contig ( dna. getId () + ":" + to_string (start + 1) + "-" + to_string (dnaSize)
	                             , dna. seq. substr (start, contigLen)
	                             , false
	                             );
	            contig. saveText (*outF);
	          }
  		    }
	      }
	    }
		}
		
		
		cout << contigs. size () << endl;
		cout << contigs [contigs. size () / 2] << endl;
		
		
		if (! distrFName. empty ())
		{
		  OFStream f (distrFName + dmSuff);
		  ds. saveText (f);
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



