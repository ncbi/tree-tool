// tblastn2disruption.cpp

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
*   Find disruptions in genes
*
*/
   
   
#undef NDEBUG 

#include "../common.hpp"
#include "../graph.hpp"
using namespace Common_sp;
#include "seq.hpp"  
using namespace Seq_sp;
#include "../version.inc"

#include "../common.inc"



namespace 
{
  
void process (const VectorPtr<Hsp> &hsps,
            //const SubstMat sm,
              bool bacteria,
              double ident_min,
              double coverage_min)
{
  Hsp::Merge merge (hsps, nullptr/*sm*/, 20, bacteria);  // PAR
  for (;;)
  { 
    const Hsp* origHsp = nullptr;
    AlignScore score = - score_inf;    
    const Hsp hsp (merge. get (origHsp, score)); 
    if (hsp. empty ()) 
      break; 
    if (   hsp. relIdentity ()  >= ident_min
        && hsp. qRelCoverage () >= coverage_min
       )
      for (const Disruption& disr : hsp. disrs)
      {
        cout         << hsp. sseqid
             << '\t' << hsp. qseqid
             << '\t' << disr. genesymbol_raw ()
             << '\t' << hsp. qInt. start + 1
             << '\t' << hsp. qInt. stop
             << '\t' << hsp. sInt. start + 1
             << '\t' << hsp. sInt. stop
             << '\t' << strand2char (hsp. sInt. strand)
             << '\t' << hsp. length
             << '\t' << hsp. nident;
        if (verbose ())
           cout << '\t' << hsp. str ();
        cout << '\n';
      }
  }
}
  
}



// ThisApplication

struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Find gene disruptions using tblastn output.\nOutput: qseqid sseqid <Disruption::genesymbol_raw()>")
    {
      version = VERSION;
      // Input
      addPositional ("tblastn", "tblastn output in the format: qseqid sseqid qstart qend sstart send qseq sseq. Ordered by qseqid, sseqid"); 
  	//addKey ("matrix", "Protein matrix", "BLOSUM62");
  	  addFlag ("stop_codon", "Query proteins end with stop codons");
  	  addFlag ("bacteria", "Bacterial genome");
  	  addKey ("ident_min", "Min. identity of the reference protein", "0.8");
  	  addKey ("coverage_min", "Min. coverage of the reference protein", "0.5");
    }



  void body () const final
  {
    const string blastFName   = getArg ("tblastn");
	//const string matrix       = getArg ("matrix");
	  const bool qStopCodon     = getFlag ("stop_codon");
	  const bool bacteria       = getFlag ("bacteria");
	  const double ident_min    = arg2double ("ident_min");
	  const double coverage_min = arg2double ("coverage_min");

    QC_ASSERT (between (ident_min, 0.0, 1.0));
    QC_ASSERT (between (coverage_min, 0.0, 1.0));
    

  #if 0
    const SubstMat sm (execDir + "/matrix/" + matrix);  
    sm. qc ();
  #endif
  
    VectorOwn<Hsp> hsps;
	  {
  	  LineInput in (blastFName, 10000);  // PAR
  	  string qseqid_prev;
  	  string sseqid_prev;
  	  while (in. nextLine ())
  	    try 
  	    { 
  	      const auto hsp = new Hsp (in. line, true/*qProt*/, false/*sProt*/, true/*aProt*/, qStopCodon, bacteria);
  	      hsp->qc ();
  	      if (! (   qseqid_prev == hsp->qseqid
  	             && sseqid_prev == hsp->sseqid
  	            )
  	         )
  	      {
  	        process (hsps, /*sm,*/ bacteria, ident_min, coverage_min);
  	        hsps. deleteData ();
            if (qseqid_prev > hsp->qseqid)
              throw runtime_error ("Non-sorted qseqid");
            if (   qseqid_prev == hsp->qseqid 
                && sseqid_prev >  hsp->sseqid
               )
              throw runtime_error ("Non-sorted sseqid");
  	      }
  	      qseqid_prev = hsp->qseqid;
  	      sseqid_prev = hsp->sseqid;
  	      hsps << hsp;
  	    }
	      catch (const exception &e)
	        { throw logic_error (string (e. what ()) + "\n" + blastFName + ": " + in. lineStr ()); }
  	}
  	process (hsps, /*sm,*/ bacteria, ident_min, coverage_min);
  }
};



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



