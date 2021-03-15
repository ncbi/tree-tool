// hmm_tc1.cpp

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
*   Find TC1 for an HMM
*
*/


#undef NDEBUG
#include "../common.inc"

//#include <objects/seqloc/Seq_id.hpp>

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
using namespace DM_sp;
#include "seq.hpp"
using namespace Seq_sp;
#include "hmm.hpp"
using namespace Hmm_sp;
#include "../version.inc"



namespace 
{
  
  
Prob error_max = NaN;

  

struct Hit : Named
// name: sseqid
{
  Real score1;
  bool seed;
  
  Hit (const string &sseqid,
       Real score1_arg,
       bool seed_arg)
    : Named (sseqid)
    , score1 (score1_arg)
    , seed (seed_arg)
    { ASSERT (score1 > 0); }
  void saveText (ostream &os) const
    { os << name << ' ' << score1 << ' ' << seed; }
    
  static bool compare (const Hit &h1,
                       const Hit &h2)
    { return h1. score1 < h2. score1; }
};



struct Hmm : Named
// name: hmmName
{
  // PAR
  static constexpr uint decimals = 0;  
  static constexpr Real lenFrac = 0.9;
  
  // Input
  Set<string> seedNames;
  Vector<Hit> hits;
  MeanVar len {0};
  // Output
  bool separating {false};
  Real tc1 {NaN};
  MeanVar score {decimals + 1};
  size_t wrong_min {0};

  
  Hmm () = default;
  void qc () const override
    { if (! qc_on)
        return;
      Named::qc ();
      len. qc ();
      score. qc ();
      QC_ASSERT (! seedNames. empty ());
      QC_ASSERT_EQ ((Real) seedNames. size (), len. n, 1e-6);  // PAR
      QC_ASSERT (wrong_min <= hits. size ());
		  QC_ASSERT (len. getMean () > 0);
      if (! isNan (tc1))
      { QC_ASSERT (tc1 > 0);
        QC_ASSERT (geReal ((Real) hits. size (), score. n));
        QC_ASSERT (geReal ((Real) seedNames. size (), score. n));
        QC_ASSERT (! isNan (score. getMean ()));
		    if (score. getMean () <= 0)
		    { saveText (cout);
		      ERROR;
		    }
		    if (! separating)
		      { QC_ASSERT_EQ ((Real) hits. size (), score. n, 1e-6); }  // PAR
      }
    }
  void saveText (ostream &os) const override
    { os. width (4);
      { ONumber on (os, (int) decimals, false);  
  	    os << "TC1: " << floor (tc1) << " " << name;
  	  }
 	    os << "  good:" << good ();
  	  if (separating)
  	  { ONumber on (os, 1, false); 
  	    os << "  error = " << wrong_min
   	       << " (" << getError () * 100 << " %)";
  	  }
 	    os << "  ";
 	    { ONumber on (os, (int) decimals + 1, false);  
  	    score. saveText (os);
  	  }
  	  os << "  len: ";
  	  len. saveText (os);
    }
  

  bool getSeparating () const
    { for (const auto& hit : hits)
        if (! hit. seed)
          return true;
      return false;
    }  
  void find_tc1 ()  
    // Input: separating
    // Update: hits (sort)
    // Output: wrong_min, tc1, score
    { score. clear ();
      if (separating)
      { sort (hits, Hit::compare);
        size_t wrong [2 /*bool*/] = {0, 0};
        Real prev = hits [0]. score1 * 0.5;  // PAR
        for (const auto& hit : hits)
          if (! hit. seed)
            wrong [false] ++;
        ASSERT (wrong [false] <= hits. size ());
        wrong_min = hits. size () + 1;
        tc1 = NaN;
        for (const auto& hit : hits)
        { if (minimize (wrong_min, wrong [false] + wrong [true]))
            tc1 = (prev + hit. score1) / 2;
          if (hit. seed)
          { wrong [true] ++;
            score << hit. score1;
          }
          else
            wrong [false] --;
          prev = hit. score1;
        }
        if (minimize (wrong_min, wrong [false] + wrong [true]))
          tc1 = prev * 1.5;  // PAR
        size_t missedPositives = seedNames. size ();  
        missedPositives -= (size_t) score. n;
        wrong_min += missedPositives;
      }
      else
      { for (const auto& hit : hits)
          score << hit. score1;
        tc1 = len. getMean () * lenFrac;
      }
    }
  Prob getError () const
    { return (Real) wrong_min / (Real) seedNames. size (); }
  bool good () const
    { if (separating)
        return    getError () <= error_max
             /*&& geReal (tc1 / len. getMean (), lenFrac)*/;  // PAR
  	  else
        return    fabs ((Real) seedNames. size () - score. n) < 1e-6  // PAR
               && tc1 < score. v_min
               && (score. getMean () - tc1) / score. getSD () >= 3.0;  // PAR  
    }  
};

  

  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Find TC1 for an HMM")
	  {
      version = VERSION;
		  addKey ("hmmsearch", "Output of hmmsearch of <seqFName> against one HMM");
		  addKey ("error_max", "Max. error fraction w.r.t. seeds to be \"good\"", "0");
		  addKey ("seqFName", "File with seed protein sequences");
	  }


	
	void body () const final
  {
	  const string hmmsearch =           getArg ("hmmsearch");  
	  error_max              = str2real (getArg ("error_max"));  
	  const string seqFName  =           getArg ("seqFName");  
	  ASSERT (isProb (error_max));
	  
	  
	  Hmm hmm;  
	  
	  // hmm.{seedNames,len}
    { 
		  Multifasta fa (seqFName, true);
		  while (fa. next ())
		  {
	      const Peptide seq (fa, Peptide::stdAveLen, false);
	      const string name (seq. getId ());
		    ASSERT (! name. empty ());
		    ASSERT (! contains (seq. seq, '-'));
	      hmm. seedNames << name;
	      hmm. len << (Real) seq. seq. size ();
		  }
		// ~Progress()
		}
	  	  
	  // hmm.hits
	  // Suppress repeated domains ??
	  {
	    Hmmsearch hs (hmmsearch);
	    while (hs. next ())
	    {
	      if (hmm. name. empty ())
	        hmm. name = hs. hmm_name;
	      else
	        ASSERT (hmm. name == hs. hmm_name);
	      const Hit hit (hs. prot_name, hs. score1, hmm. seedNames. contains (hs. prot_name));
	      hmm. hits << hit;
	      if (verbose ())
	      {
	        hit. saveText (cout);
	        cout << endl;
	      }
	    }
	  }
	  
	  hmm. separating = hmm. getSeparating ();
    hmm. find_tc1 ();
    hmm. qc ();

    hmm. saveText (cout);
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



