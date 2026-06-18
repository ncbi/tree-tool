// assignment.cpp

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
*   Max. assignment problem
*
*/


#undef NDEBUG

#include "common.hpp"
#include "numeric.hpp"
using namespace Common_sp;
#include "bipartite.hpp"
using namespace Bipartite_sp;
#include "version.inc"

#include "common.inc"



namespace 
{
  
  
struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Maximum assignment problem. Print maximum total score")
    {
      version = VERSION;
      addPositional ("in", "File with bipartite arcs: <job>\\t<worker>\\t<score>. Missing arcs assume to have score 0");
    }

  
	void body () const final
	{
	  const string fName = getArg ("in");
	  
	  
	  Bipartite gr; 
	  Real total_hi = 0.0;
	  {
  	  LineInput f (fName);
  	  map<string,BpNode*> name2job;  
  	  map<string,BpNode*> name2worker;  
  	  while (f. nextLine ())
  	  {
  	    string jobName    (findSplit (f. line, '\t'));
  	    string workerName (findSplit (f. line, '\t'));
  	    const Real score = str2<Real> (f. line);
  	    QC_ASSERT (! isNan (score));
  	    if (nullReal (score))
  	      continue;
  	    total_hi += score;
  	    trim (jobName);
  	    QC_ASSERT (! jobName. empty ());
  	    trim (workerName);
  	    QC_ASSERT (! workerName. empty ());
  	    
  	    BpNode* &job = name2job [jobName];
  	    if (! job)
  	      job = new BpNode (gr, false, jobName);
  	      
  	    BpNode* &worker = name2worker [workerName];
  	    if (! worker)
  	      worker = new BpNode (gr, true, workerName);
  	      
  	    if (const DiGraph::Arc* arc = job->incident (worker, true))
  	      var_cast (static_cast <const BpArc*> (arc)) -> score += score;
  	    else
  	      new BpArc (job, worker, score);
  	  }
  	}
    if (verbose ())
    	cout << "total_hi = " << total_hi << endl;
  	  	
    
    gr. assignmentSparse ();
        
      
    size_t matched = 0;
    const Real score = gr. getMatchScore (matched);
    cout << score << endl;
    if (qc_on)        
    {
      size_t n = 0;  // ??
      Real score_min = inf;  // ??
      for (const BpNode* j : gr. bpNodes [false])
        if (   ! j->match
            && positive (j->potential)
           )
        {
          n++;
          minimize (score_min, j->potential);
        //cout << *j << endl;
        //ERROR;
        }
      PRINT (n);
      PRINT (score_min);
      if (n)
        ERROR;
      for (const BpNode* w : gr. bpNodes [true])
      {
        ASSERT (w);
        IMPLY (! w->match, ! w->potential);
      }
    }
  //ASSERT (leReal (total_match_lo, score));
    if (verbose ())
      for (const BpNode* j : gr. bpNodes [false])
        if (j->match)
          cout << "Match: " << * j->match << endl;  
  }  
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



