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
    
//QC_ASSERT (! name. empty ());
  QC_ASSERT (Common_sp::finite (potential));
//QC_ASSERT (! negative (potential));
  if (match)
  {
    QC_ASSERT (match->node [out] == this);
    QC_ASSERT (static_cast <const BpNode*> (match->node [! out]) -> match == match);
  }
  QC_IMPLY (prev, prev->node [out] == this);
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
  
//QC_ASSERT (! nullReal (score));
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
          if (! other->prev)
          {
            var_cast (other) -> prev = sa;
            flagged [! out] << other;
          }
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
    if (verbose (-1))
      cout << "New match: " << *sa << endl;  
    IMPLY (sa1, sa1->matched ());
    var_cast (sa) -> use ();  // Destroys sa1->matched()
    if (! sa1)
      break;
    const BpNode* w = sa1->getNode (true);
    ASSERT (w);    
    sa = w->prev;
  }
  ASSERT (sa);
}



void Bipartite::assignment ()
{
	VectorPtr<BpArc> arcs;  
	for (const DiGraph::Node* n : nodes)
  	for (const DiGraph::Arc* arc : n->arcs [true])
  	  arcs << static_cast <const BpArc*> (arc);  	  
  if (arcs. empty ())
    return;
  if (arcs. size () == 1)
  {
    var_cast (arcs [0]) -> use ();
    if (verbose ())
      cout         << "total_lo = "   << arcs [0] -> score 
           << '\t' << "total_init = " << arcs [0] -> score
           << endl;   
    return;
  }

  arcs. sortPtr ();  		    	
  Real total_match_lo = 0.0;  
  Real total_match_init = 0.0;
  for (const BpArc* ca : arcs)
  {
    ASSERT (ca);
    if (   ! ca->getNode (false) -> match
        && ! ca->getNode (true)  -> match
        && ca->zero ()
       )
    {
      var_cast (ca) -> use ();
      total_match_init += ca->score;
    }
    if (   ! ca->getNode (false) -> naive_match
        && ! ca->getNode (true)  -> naive_match
       )
    {
      total_match_lo += ca->score;
      var_cast (ca->getNode (false)) -> naive_match = true;
      var_cast (ca->getNode (true))  -> naive_match = true;
    }
  }
  ASSERT (leReal (total_match_init, total_match_lo));
  if (verbose ())
    cout         << "total_lo = "   << total_match_lo 
         << '\t' << "total_init = " << total_match_init
         << endl;   
         
         
  // Auxiliary BpNode's
  {
    const size_t diff = difference ( bpNodes [false]. size ()
                                   , bpNodes [true].  size ()
                                   );
    const bool out = bpNodes [false]. size () > bpNodes [true]. size ();
    FFOR (size_t, i, diff)
      new BpNode (*this, out, noString);
    ASSERT (bpNodes [false]. size () == bpNodes [true]. size ());
  }
      
  // Auxiliary BpArc's
  for (const BpNode* a : bpNodes [false])
  {
    ASSERT (a);
    VectorPtr<DiGraph::Node> neighborhood (a->getNeighborhood (true));
    neighborhood. sort ();
    for (const BpNode* b : bpNodes [true])
      if (! neighborhood. containsFast (b))
        new BpArc (var_cast (a), var_cast (b), 0.0);
  }


  for (;;)
  {
    qc (); 
    maximumMatching ();
    qc ();  
    
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
    if (! arc_best)
    {
      ASSERT (score_min == inf);
      break;
    }
      
    if (verbose (-1))
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
  
  
  if (qc_on)      
  {
    size_t matched = 0;
    const Real total = getMatchScore (matched);
    ASSERT (leReal (total_match_lo, total));    
    for (const bool out : {false, true})
      for (const BpNode* w : bpNodes [out])
      {
        QC_ASSERT (w);
        QC_ASSERT (w->match);
      }
  }
} 



}  // namespace

