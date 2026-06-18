// bipartite.hpp

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
*   Bipartite graph tools
*
*/


#undef NDEBUG

#include "bipartite.hpp"

#include "common.inc"



namespace Bipartite_sp
{
  

// BpNode  

BpNode::BpNode (Bipartite &bp,
                bool out_arg,
                const string &name_arg)
: DiGraph::Node (bp)
, out (out_arg)
, name (name_arg)
{
  bp. bpNodes [out] << this;
}
  


void BpNode::qc () const 
{ 
  if (! qc_on)
    return;
  DiGraph::Node::qc ();
    
  QC_ASSERT (! name. empty ());
  QC_ASSERT (Common_sp::finite (potential));
//QC_ASSERT (! negative (potential));
  if (match)
  {
    QC_ASSERT (match->node [out] == this);
    QC_ASSERT (static_cast <const BpNode*> (match->node [! out]) -> match == match);
  }
  QC_IMPLY (prev, prev->node [out] == this);
}



void BpNode::saveText (ostream &os) const 
{ 
  os << name << endl;
  for (const DiGraph::Arc* arc : arcs [! out])
    os << "  " << *arc << (static_cast <const BpArc*> (arc) == match ? " M" : "") << endl;
}


  
// BpArc  

void BpArc::qc () const 
{ 
  if (! qc_on)
    return;
  DiGraph::Arc::qc ();  
    
  QC_ASSERT (! getNode (false) -> out);
  QC_ASSERT (getNode (true) -> out);
  QC_ASSERT (node [false] != node [true]);
  
  QC_ASSERT (! nullReal (score));
  QC_ASSERT (Common_sp::finite (score));
    
  QC_ASSERT (! negative (reducedScore ()));
  QC_ASSERT (Common_sp::finite (reducedScore ()));
    
  QC_IMPLY (matched (), zero ());
}




// Bipartite

void Bipartite::qc () const
{
  if (! qc_on)
    return;
  DiGraph::qc ();
    
  QC_ASSERT (bpNodes [false]. size () + bpNodes [true]. size () == nodes. size ());
  for (const bool out : {false, true})
    for (const BpNode* n : bpNodes [out])
      QC_ASSERT (n->out == out);
}


  
void Bipartite::expand (const VectorPtr<BpNode> &start) const
{
  for (const DiGraph::Node* n : nodes)
     var_cast (static_cast <const BpNode*> (n)) -> prev = nullptr;
     
  if (start. empty ())
    return;
    
  bool out = start [0] -> out;
  if (qc_on)
    for (const BpNode* n : start)
      QC_ASSERT (n->out == out);
  
  array<VectorPtr<BpNode>,2/*out*/> flagged;
  flagged [out] = start;
  while (! flagged [out]. empty ())
  {
    flagged [! out]. clear ();
    for (const BpNode* n : flagged [out])
      for (const DiGraph::Arc* arc : n->arcs [! out])
      {
        const BpArc* sa = static_cast <const BpArc*> (arc);
        ASSERT (sa);
        if (   sa->zero ()
            && out == sa->matched ()
           )
        {
          const BpNode* other = sa->getNode (! out);
          ASSERT (other);
          var_cast (other) -> prev = sa;
          flagged [! out] << other;
        }
      }
    toggle (out);
  }
}        



const BpArc* Bipartite::findAugmentingPath () const
{
  VectorPtr<BpNode> start;  
  for (const BpNode* j : bpNodes [false])
    if (! j->match)
      start << j;
  expand (start);

  const BpArc* goal_arc = nullptr;
  for (const BpNode* w : bpNodes [true])
    if (   w->prev
        && ! w->match
       )
    {
      goal_arc = w->prev;
      break;
    }      
  
  return goal_arc;
}



void Bipartite::rematch (const BpArc* goal_arc)
{
  const BpArc* sa = goal_arc;
  for (;;)
  {
    ASSERT (sa);
    ASSERT (! sa->matched ());
    const BpNode* j = sa->getNode (false);
    ASSERT (j);
    const BpArc* sa1 = j->prev;
    IMPLY (sa1, sa1->matched ());
    var_cast (sa) -> use ();
    if (verbose ())
      cout << "New match: " << *sa << endl;  
    if (! sa1)
      break;
    const BpNode* w = sa1->getNode (true);
    ASSERT (w);    
    sa = w->prev;
  }
  ASSERT (sa);
}



void Bipartite::assignmentComplete ()
{
	VectorPtr<BpArc> arcs;  
	for (const DiGraph::Arc* arc : arcs)
	  arcs << static_cast <const BpArc*> (arc);  	  
  arcs. sortPtr ();  	
	    	
  Real total_match_lo = 0.0;  
  for (const BpArc* ca : arcs)
  {
    ASSERT (ca);
    if (   ! ca->getNode (false) -> match
        && ! ca->getNode (true)  -> match
        && ca->zero ()
       )
      var_cast (ca) -> use ();
    if (   ! ca->getNode (false) -> naive_match
        && ! ca->getNode (true)  -> naive_match
       )
    {
      total_match_lo += ca->score;
      var_cast (ca->getNode (false)) -> naive_match = true;
      var_cast (ca->getNode (true))  -> naive_match = true;
    }
  }
  if (verbose ())
    cout << "total_lo = " << total_match_lo << endl;   


  for (;;)
  {
    qc (); 
    maximumMatching ();
    qc ();  
    
  #if 0
    size_t matched = 0;
    const Real score = getMatchScore (matched);
    prog (to_string (matched) + " / " + to_string (gr. jobs. size ()) + " (" + to_string (score) + ")");
  #endif
  
    const BpArc* arc_best = nullptr;
    Real score_min = inf;
    for (const BpNode* j : bpNodes [false])
      if (   ! j->match
          || j->prev
         )
        for (const DiGraph::Arc* arc : j->arcs [true])
        {
          const BpArc* sa = static_cast <const BpArc*> (arc);
          ASSERT (sa);
          BpNode* w = var_cast (sa->getNode (true));
          ASSERT (w);
          if (   ! w->prev
              && minimize (score_min, sa->reducedScore ())
             )
            arc_best = sa;
        }
    ASSERT (! negative (score_min));
  //if (nullReal (score_min))
    //break;          
    if (! arc_best)
    {
      ASSERT (score_min == inf);
      break;
    }
      
    if (verbose ())
      cout << "To make zero(): " << *arc_best << endl;  
    ASSERT (Common_sp::finite (score_min));
    for (const BpNode* j : bpNodes [false])
      if (   ! j->match
          || j->prev
         )
        var_cast (j) -> potential -= score_min;
    for (const BpNode* w : bpNodes [true])
      if (w->prev)
        var_cast (w) -> potential += score_min;
    ASSERT (arc_best->zero ());
  }
} 



void Bipartite::assignmentSparse ()
{
  assignmentComplete ();
  
  
  for (;;)
  {
    VectorPtr<BpNode> start;  
    for (const BpNode* j : bpNodes [false])
      if (   ! j->match
          && ! nullReal (j->potential)
         )
        start << j;
    expand (start);
  
    const BpArc* goal_arc = nullptr;
    for (const BpNode* w : bpNodes [false])
      if (   w->prev
          && nullReal (w->potential)
         )
      {
        ASSERT (w->match);
        goal_arc = w->prev;
        break;
      }
      
    if (! goal_arc)
      break;
    
    ASSERT (goal_arc->matched ());
    const BpNode* j = goal_arc->getNode (false);
    ASSERT (j);
    ASSERT (nullReal (j->potential));
    var_cast (j) -> match = nullptr;
    const BpNode* w = goal_arc->getNode (true);
    ASSERT (w);
    rematch (w->prev);
  }
}



}  // namespace

