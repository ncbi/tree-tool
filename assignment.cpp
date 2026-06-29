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
    : Application ("Maximum assignment problem in a sparse graph. Print maximum total score")
    {
      version = VERSION;
      addPositional ("in", "File with bipartite arcs: <job>\\t<worker>\\t<score>. Missing arcs are assumed to have score 0");
    }

  
	void body () const final
	{
	  const string fName = getArg ("in");
	  
	  
    map<const DiGraph::Node*, Bipartite> components;
	  {
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
      	      
      // --> graph.hpp ??
      gr. connectedComponents ();
      VectorPtr<BpNode> nodes;
      for (const DiGraph::Node* n : gr. nodes)
        nodes << static_cast <const BpNode*> (n);
      for (const BpNode* n : nodes)
      {
        ASSERT (n);
        Bipartite& bp = components [var_cast (n) -> getDisjointCluster ()];
        var_cast (n) -> moveTo (bp);
        bp. bpNodes [n->out] << n;
      }
      ASSERT (gr. nodes. empty ());
    }
  	  	

    // ??
    VectorPtr<Bipartite> componentsVec;  componentsVec. reserve (components. size ());
    for (const auto& it : components)
      componentsVec << & it. second;
    componentsVec. sortPtr ();


    Real total = 0.0;
    {
      Progress prog (components. size ());
    //for (auto& it : components)
      for (const Bipartite* bp_ : componentsVec)
      {
        prog ();
      //Bipartite& bp = it. second;
        Bipartite& bp = * var_cast (bp_);  // ??
        bp. qc ();
        bp. assignmentSparse ();
        size_t matched = 0;
        total += bp. getMatchScore (matched);
        if (verbose (-1))
          for (const BpNode* j : bp. bpNodes [false])
            if (j->match)
              cout << "Match: " << * j->match << endl;  
      }
    }
    cout << total << endl;
  }  
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



