// tblastn2marker_euk.cpp

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
*   Find best gene matches and print proteins
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
  
  
//unique_ptr<OFStream> gapF;
char delim {'\0'};
double threshold = NaN;

// PAR
constexpr double complexity_min = 3.0;  
static_assert (complexity_min >= Peptide::stdMinComplexity);



string qseqid2gene (const string &qseqid)
{
  if (   ! qseqid. empty () 
      && delim
     )
  {
    const size_t pos = qseqid. find (delim);
    if (pos == string::npos)
      throw runtime_error ("Protein identifier " + strQuote (qseqid) + " has no delimiter");
    QC_ASSERT (pos);
    return qseqid. substr (0, pos);
  }
  
  return qseqid;
}



void report (const Hsp &hsp,
             AlignScore score)
{
  if (hsp. empty ())
    return;
    
  const string name (qseqid2gene (hsp. qseqid) + " " + hsp. qseqid + " " + hsp. sseqid + " strand=" + to_string (hsp. sInt. strand == 1 ? 1 : 0) 
                     + " score=" + to_string (int (score)) + " len=" + to_string (hsp. cleanSseq (). size ()));
  Peptide pep (name, hsp. cleanSseq (), false);
  pep. pseudo = true;
  pep. qc ();
  if (pep. getComplexity () >= complexity_min)  
    pep. saveText (cout);
}



//OFStream f1 ("aa.log");  



void process (const VectorPtr<Hsp> &origHsps,
              const SubstMat &sm,
              Hsp &hsp_best,
              AlignScore &score_best)
{
  Hsp::Merge merge (origHsps, & sm, 0, false/*bacteria*/);  
  const Hsp* origHsp = nullptr;
  AlignScore score = - score_inf;
  Hsp hsp (merge. get (origHsp, score)); 
  if (hsp. empty ())
    return;
  if (score / (double) hsp. cleanSseq (). size () < threshold)
    return;
//f1 << hsp. qseqid << ' ' << hsp. sseqid << ' ' << score << '\n';  
    
  if (qseqid2gene (hsp_best. qseqid) == qseqid2gene (hsp. qseqid))
  {
    if (score_best < score)
    {
      hsp_best = std::move (hsp);
      score_best = score;
    }
    return;
  }

  report (hsp_best, score_best);
    
  hsp_best = std::move (hsp);
  score_best = score;
}



struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Find best protein matches and print proteins")
    {
      version = VERSION;
  	  addPositional ("tblastn", "tblastn output in format: " + string (Hsp::format [true]) + ". Ordered by qseqid, sseqid. If -delimeter then qseqid = <gene><delimiter><variant>, where <gene> has no <delimiter>.");
  	  addKey ("delimiter", "Delimiter character separating gene names from variant names in <tblastn>");
  	  addKey ("threshold", "Min. protein score to length ratio to be reported", "4");
  	  addKey ("matrix", "Protein matrix", "BLOSUM62");
  	//addKey ("gap_stats", "Output file with gap lengths");
    }



	void body () const final
  {
	  const string tblastnFName = getArg ("tblastn");
	  const string delimiterS   = getArg ("delimiter"); 
	               threshold    = arg2double ("threshold");
	  const string matrix       = getArg ("matrix");
	//const string gap_stats    = getArg ("gap_stats");
  
  
    if (! delimiterS. empty ())
    {
      QC_ASSERT (delimiterS. size () == 1);
      delim = delimiterS [0];
    }
  
    const SubstMat sm (execDir + "/matrix/" + matrix);  
    sm. qc ();
  
  #if 0
	  if (! gap_stats. empty ())
      gapF. reset (new OFStream (gap_stats));
  #endif
  
    VectorOwn<Hsp> hsps;
	  LineInput in (tblastnFName, 10000);  // PAR
	  string qseqid_prev;
	  string sseqid_prev;
	  Hsp hsp_best;
    AlignScore score_best = - score_inf;
	  while (in. nextLine ())
    { 
      // qseqid, sseqid
      const Hsp* hsp = nullptr;
	    try 
	    { 
 	      unique_ptr<Hsp> hsp_ (new Hsp (in. line, true/*qProt*/, false/*sProt*/, true/*aProt*/, false/*qStopCodon*/, false/*bacteria*/));
	      hsp_->qc ();
        hsp = hsp_. release ();
	    }
	    catch (const exception &e)
	      { throw logic_error (string (e. what ()) + "\n" + tblastnFName + ": " + in. lineStr ()); }
	    ASSERT (hsp);

      if (! (   qseqid_prev == hsp->qseqid
             && sseqid_prev == hsp->sseqid
            )
         )
      {
    	  process (hsps, sm, hsp_best, score_best);
        hsps. deleteData ();
        if (qseqid_prev > hsp->qseqid)
          throw runtime_error ("Non-ordered qseqid");
        if (   qseqid_prev == hsp->qseqid 
            && sseqid_prev >  hsp->sseqid
           )
          throw runtime_error ("Non-ordered sseqid");
      }
      hsps << hsp;
      
      qseqid_prev = hsp->qseqid;
      sseqid_prev = hsp->sseqid;
    }
 	  process (hsps, sm, hsp_best, score_best);
 	  report (hsp_best, score_best);
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



