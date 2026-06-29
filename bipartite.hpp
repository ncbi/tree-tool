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


#ifndef BIPARTITE_HPP_75943  // random number
#define BIPARTITE_HPP_75943


#include "common.hpp"
#include "graph.hpp"
#include "numeric.hpp"
using namespace Common_sp;



namespace Bipartite_sp
{
  
  
struct Bipartite;
struct BpArc;



struct BpNode : DiGraph::Node
{
  const bool out;
	  // false: "job"
	  // true:  "worker"
  const string name;
  bool naive_match {false};
  const BpArc* match {nullptr};
    // Incident to *this
  Real potential {0.0};
  const BpArc* prev {nullptr};    
    // For augmenting paths

          
  BpNode (Bipartite &bp,
          bool out_arg,
          const string &name_arg);
  void qc () const override;
  void saveText (ostream &os) const override
    { os << name << endl;
      for (const DiGraph::Arc* arc : arcs [! out])
        os << "  " << *arc << endl;
    }
};



struct BpArc : DiGraph::Arc
// From BpNode where !out to BpNode where out
{
  Real score {NaN};
    // Greater is better
  
  
  BpArc (BpNode* start,
         BpNode* end,
         Real score_arg)
    : DiGraph::Arc (start, end)
    , score (score_arg)
    { maximize (start->potential, score); }
  void qc () const override;
  void saveText (ostream &os) const override
    { os << getNode (false) -> name;
      if (getNode (false) -> match)
        os << "(M)";
      os << '-' << getNode (true) -> name;
      if (getNode (true) -> match)
        os << "(M)";
      os << ": " << reducedScore () 
         << " (" << getNode (false) -> potential << " + " << getNode (true) -> potential << " - " << score << ')';
      if (matched ())
        os << " M";
    }

    
  const BpNode* getNode (bool out) const
    { return static_cast <const BpNode*> (node [out]); }
  Real reducedScore () const
    { return   getNode (false) -> potential
             + getNode (true)  -> potential
             - score;
    }
  bool matched () const
    { return getNode (false) -> match == this; }
  bool zero () const
    { return nullReal (reducedScore ()); }
  bool operator< (const BpArc &other) const
    { return score > other. score; }
  void use ()
    { for (const bool out : {false, true})
        var_cast (static_cast <const BpNode*> (node [out])) -> match = this;
    }
};



struct Bipartite : DiGraph
// Of BpNode* and BpArc*
{
	array<VectorPtr<BpNode>,2/*bool out*/> bpNodes;


  Bipartite () = default;
  void qc () const override;

 
  // Time ??
private:
  void expand (const VectorPtr<BpNode> &start,
               bool from) const;
    // Output: BpNode::prev
  const BpArc* findAugmentingPath (bool from) const;
    // Find an augmenting path on zero()-BpArc's
    // Return: end BpArc* of the path
    //         may be nullptr
    // Invokes: expand()
  void rematch (const BpArc* goal_arc,
                bool from);
    // Input: BpNode::prev
  void maximumMatching (bool from)
    { while (const BpArc* goal_arc = findAugmentingPath (from))
        rematch (goal_arc, from);
    }
    // Find a maximum matching on zero()-BpArc's
    // Fulkerson-Ford
    // Output: BpNode::prev define a minimum vertex cut    
public:
  bool assignmentComplete (Real &total_match_lo);
    // Hungarian algorithm 
    // Optimal for a complete graph
    // Return: maximumMatching( )is invoked
    // Invokes: maximumMatching()
  void assignmentSparse ();
    // Invokes: assignmentComplete(), expand()
  Real getMatchScore (size_t &matched) const
    { Real s = 0.0;    
      matched = 0;
      for (const BpNode* j : bpNodes [false])
        if (j->match)
        { s += j->match->score;
          matched++;
        }
      return s;
    }

  // ??
  bool operator< (const Bipartite &other) const
    { return nodes. size () < other. nodes. size (); }
};



}  // namespace



#endif
