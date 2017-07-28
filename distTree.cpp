// distTree.cpp

#undef NDEBUG
#include "../cpp/common.inc"

#include "distTree.hpp"

#include "prediction.hpp"



namespace DistTree_sp
{


// VarianceType

const StringVector varianceTypeNames {"lin", "exp", "linExp"};
VarianceType varianceType = varianceType_linExp;




// DTNode::Closest

void DTNode::Closest::qc () const
{
  ASSERT (node);
  ASSERT (node->len > 0);
  ASSERT (! node->inDiscernable ());
  ASSERT (node != node->getTree (). root);
  ASSERT (absCriterion_delta >= 0);
  ASSERT (leafLen >= 0);
  ASSERT (arcLen >= 0);
  ASSERT (arcLen <= node->len);
}



Steiner* DTNode::Closest::insert () 
{
  auto st = new Steiner ( const_cast <DistTree&> (node->getDistTree ())
                        , const_static_cast <Steiner*> (node->getParent ())
                        , node->len - arcLen
                        );
  const_cast <DTNode*> (node) -> setParent (st);
  const_cast <DTNode*> (node) -> len = arcLen;
  st->subtreeLeaves = node->subtreeLeaves;
  return st;
}




// DTNode

DTNode::DTNode (DistTree &tree,
                Steiner* parent_arg,
    	          Real len_arg)
: TreeNode (tree, parent_arg)  // DistTree must be declared
, len (len_arg)
{}



void DTNode::qc () const
{ 
	TreeNode::qc ();

  if ((bool) getParent () != ! isNan (len))
  {
    cout << getName () << endl;
    getTree (). print (cout);
    ERROR;
  }
  IMPLY (! isNan (len), len >= 0);	
  IMPLY (! childrenDiscernable (), ! inDiscernable ());
  IMPLY (getDistTree (). optimizable (), attr);
  
  IMPLY (paths, errorDensity >= 0);
}



void DTNode::saveContent (ostream& os) const
{ 
  ONumber oNum (os, 6, true);  // PAR
  {
		os << "len=" << len;
		if (subtreeLen. weights)
		  os << "  len_mean=" << getHeight () << "  len_weight=" << subtreeLen. weights;
	}
  if (paths)
	  os << "  err_density=" << errorDensity << "  paths=" << paths;
}



const DistTree& DTNode::getDistTree () const
{
  return static_cast <const DistTree&> (getTree ());
}



void DTNode::addAttr () 
{
  ASSERT (! attr);
    
  DistTree& tree = const_cast <DistTree&> (getDistTree ());
  tree. attrNum_max++;
  ASSERT (tree. attrNum_max);
  const string name ("v" + toString (tree. attrNum_max));
  attr = new CompactBoolAttr1 (name, tree. ds);
}



const Leaf* DTNode::inDiscernable () const
{ 
  const Leaf* g = asLeaf ();
  return g && ! g->discernable ? g : nullptr; 
}



void DTNode::saveFeatureTree (ostream &os,
                              size_t offset) const
{
  FOR (size_t, i, offset)
    os << ' ';
  ONumber on (os, 6, true);
  string s = (asLeaf () ? "s" : "") + getName ();
  os << s << ": t=" << (isNan (len) ? 0 : len) << "  C=0  dC=+0-0" << endl;
  if (asLeaf ())
  {
    FOR (size_t, i, offset + 2)
      os << ' ';
    os << 'g' << getName () << ": " << endl;
  }
  else
  	for (const DiGraph::Arc* arc : arcs [false])
  	  static_cast <const DTNode*> (arc->node [false]) -> saveFeatureTree (os, offset + 2);
}



void DTNode::setSubtreeLenUp ()
{
  subtreeLen. clear ();
  if (! isLeaf ())
  	for (const DiGraph::Arc* arc : arcs [false])
  	{
  	  DTNode* child = static_cast <DTNode*> (arc->node [false]);
  	  child->setSubtreeLenUp ();
  	  subtreeLen. add (child->getHeight () + child->len, child->subtreeLen. weights + child->len);
  	}
}



void DTNode::setGlobalLenDown (DTNode* &bestDTNode,
                               Real &bestDTNodeLen_new,
                               WeightedMeanVar &bestGlobalLen)
{
  if (const DTNode* parent_ = static_cast <const DTNode*> (getParent ())) 
  {
    WeightedMeanVar parentSubtreeLen (parent_->subtreeLen /*global len*/);  
    {
      WeightedMeanVar subtreeLen1;
      subtreeLen1. add (getHeight () + len, subtreeLen. weights + len);
      parentSubtreeLen. subtract (subtreeLen1);  
    }
      // Subtree goes upward from *parent_
    Real lenDelta =   len / 2 
                    + (  (parentSubtreeLen. getMean () + parentSubtreeLen. weights)
                       - (      subtreeLen. getMean () +       subtreeLen. weights)
                      ) / 4;
    maximize (lenDelta, 0.0);
    minimize (lenDelta, len);
    WeightedMeanVar globalLen;
    globalLen. add (      subtreeLen. getMean () + lenDelta,               subtreeLen. weights + lenDelta);
    globalLen. add (parentSubtreeLen. getMean () + (len - lenDelta), parentSubtreeLen. weights + (len - lenDelta));
    if (   ! bestDTNode 
        || bestGlobalLen. getMean () > globalLen. getMean ()
       )
    {
      bestDTNode = this;
      bestDTNodeLen_new = lenDelta;
      bestGlobalLen = globalLen;
    }
    subtreeLen. add (parentSubtreeLen. getMean () + len, parentSubtreeLen. weights + len);
      // global len
  }
	for (const DiGraph::Arc* arc : arcs [false])
	  static_cast <DTNode*> (arc->node [false]) -> setGlobalLenDown (bestDTNode, bestDTNodeLen_new, bestGlobalLen);
}



void DTNode::setSubtreeLeaves (/*const Set<const Leaf*>* superSet*/)
{
  for (const DiGraph::Arc* arc : arcs [false])
  {
    const DTNode* child = static_cast <const DTNode*> (arc->node [false]);
    const_cast <DTNode*> (child) -> setSubtreeLeaves ();
    ASSERT (subtreeLeaves. size () == child->subtreeLeaves. size ());
    FOR (size_t, i, subtreeLeaves. size ())
      if (child->subtreeLeaves [i])
        subtreeLeaves [i] = true;
  }
}



void DTNode::findClosestNode (const Leaf2dist &leaf2dist,
                              Leaf2dist &leaf2hat_dist,
                              Closest &closest) const
{
  ASSERT (! leaf2dist. empty ());
  ASSERT (leaf2dist. size () == subtreeLeaves. size ());
  ASSERT (leaf2hat_dist. empty () || leaf2hat_dist. size () == leaf2dist.size ());


  bool foundSubLeaf   = false;
  bool foundSuperLeaf = false;
	FOR (size_t, i, leaf2dist. size ())
	  if (! isNan (leaf2dist [i]))
	  {
	    if (subtreeLeaves [i])
	      foundSubLeaf = true;
	    if (! subtreeLeaves [i])
	      foundSuperLeaf = true;
	  }
	if (! foundSubLeaf)
	  return;


  // leaf2hat_dist
  // Depth-first order
  if (leaf2hat_dist. empty ())
  {
    if (arcs [false]. empty ())
    {
      leaf2hat_dist. resize (leaf2dist. size (), NAN);
      const Leaf* leaf = asLeaf();
      ASSERT (leaf);
      ASSERT (leaf->index != NO_INDEX);
      ASSERT (! isNan (leaf2dist [leaf->index]));
      for (const DiGraph::Node* node : getTree (). nodes)
        if (const Leaf* g = static_cast <const DTNode*> (node) -> asLeaf ())
          if (   g->index != NO_INDEX
              && ! isNan (leaf2dist [g->index])
             )
            leaf2hat_dist [g->index] = DistTree::path2prediction (DistTree::getPath (leaf, g));
    }
  }
  else
    // *getParent() --> *this
  	FOR (size_t, i, leaf2dist. size ())
  	{
  	  const Real u = subtreeLeaves [i] ? 1 : -1;
  	  leaf2hat_dist [i] -= len * u;
  	}

#ifndef NDEBUG
  Leaf2dist leaf2hat_dist_ (leaf2hat_dist);
#endif
  for (const DiGraph::Arc* arc : arcs [false])
  {
    const DTNode* child = static_cast <const DTNode*> (arc->node [false]);
    if (const Leaf* g = child->asLeaf ())
      if (g->index == NO_INDEX)
        continue;
    const_cast <DTNode*> (child) -> findClosestNode (leaf2dist, leaf2hat_dist, closest);
  #ifndef NDEBUG
    if (leaf2hat_dist_. empty ())
      leaf2hat_dist_ = leaf2hat_dist;
    else
    {
      ASSERT (leaf2hat_dist. size () == leaf2hat_dist_. size ());
    	FOR (size_t, i, leaf2hat_dist. size ())
      {
        ASSERT (isNan (leaf2hat_dist [i]) == isNan (leaf2hat_dist_ [i]));
        if (! isNan (leaf2hat_dist [i]))
          { ASSERT_EQ (leaf2hat_dist [i], leaf2hat_dist_ [i], 1e-3); } // PAR
      }
    }
  #endif
  }
  ASSERT (! leaf2hat_dist. empty ());
	  

	if (this == getTree (). root)
	  return;
	  

	if (len > 0 && foundSuperLeaf)
	{
  	Real mult_sum = 0;
  	Real u_avg = 0;	
  	Real delta_avg = 0;
  	FOR (size_t, i, leaf2dist. size ())
  	{
  	  const Real dist = leaf2dist [i];
	    if (isNan (dist))
	      continue;
  	  ASSERT (positive (dist));
  	  const Real mult = dissim2mult (dist);
  	  mult_sum += mult;
  	  const Real u = subtreeLeaves [i] ? 1 : -1;
  	  const Real delta = dist - leaf2hat_dist [i];
  	  u_avg     += mult * u;
  	  delta_avg += mult * delta;
  	}
  	ASSERT (positive (mult_sum));
  	u_avg     /= mult_sum;
  	delta_avg /= mult_sum;
  	
  	Real u_var = 0;
  	Real delta_u_cov = 0;
  	FOR (size_t, i, leaf2dist. size ())
  	{
  	  const Real dist = leaf2dist [i];
	    if (isNan (dist))
	      continue;
  	  const Real mult = dissim2mult (dist);
  	  const Real u = subtreeLeaves [i] ? 1 : -1;
  	  const Real delta = dist - leaf2hat_dist [i];
  	  delta_u_cov += mult * (delta - delta_avg) * (u - u_avg);
  	  u_var       += mult * sqr (u - u_avg);
  	}
  	u_var       /= mult_sum;
  	delta_u_cov /= mult_sum;  	
  	ASSERT (u_var > 0);

  	const Real arcLen = min (max (0.0, delta_u_cov / u_var), len);
  	const Real leafLen = max (0.0, delta_avg - u_avg * arcLen);
  	
  	Real absCriterion_delta = 0;
  	FOR (size_t, i, leaf2dist. size ())
  	{
  	  const Real dist = leaf2dist [i];
	    if (isNan (dist))
	      continue;
  	  const Real mult = dissim2mult (dist);
  	  const Real u = subtreeLeaves [i] ? 1 : -1;
  	  const Real delta = dist - leaf2hat_dist [i];
  	  absCriterion_delta += mult * sqr (delta - arcLen * u - leafLen);
  	}
  	ASSERT (absCriterion_delta >= 0);
  	
  	if (absCriterion_delta < closest. absCriterion_delta)
  	{
  	  closest = Closest (this, absCriterion_delta, leafLen, arcLen);
  	  closest. qc ();
  	}
  }


  // *this --> *getParent()
	FOR (size_t, i, leaf2dist. size ())
	{
	  const Real u = subtreeLeaves [i] ? 1 : -1;
	  leaf2hat_dist [i] += len * u;
	}
}




// Steiner


Steiner::Steiner (DistTree &tree,
      	          Steiner* parent_arg,
      	          Real len_arg)
: DTNode (tree, parent_arg, len_arg)  // DistTree must be declared
{
  for (const bool b : Bool)
    bootstrap [b] = 0;
}



void Steiner::qc () const
{
	DTNode::qc ();
	  
	ASSERT (! isLeaf ());	
//IMPLY (reprLeaf, getChildren (). contains (reprLeaf));
}



void Steiner::reverseParent (const Steiner* target, 
                             Steiner* child)
{
  ASSERT (target);
  ASSERT (descendentOf (target));
  IMPLY (child, child->getParent () == this);  
  
  const DTNode* oldParent = static_cast <const DTNode*> (getParent ());
  if (this != target)
  {
    ASSERT (oldParent);
    const_cast <Steiner*> (oldParent -> asSteiner ()) -> reverseParent (target, this);
  }
  
  if (child)
  {
    setParent (child);
    // Arc-specific data
    len  = child->len;
    attr = child->attr;
  }
}



void Steiner::makeRoot (Steiner* ancestor2descendant)
{
  ASSERT (ancestor2descendant);
  ASSERT (descendentOf (ancestor2descendant));

  const DTNode* parent_new = static_cast <const DTNode*> (ancestor2descendant->getParent ());
  // Arc-specific data
  const Real len_new               = ancestor2descendant->len;
  const CompactBoolAttr1* attr_new = ancestor2descendant->attr;
  
  reverseParent (ancestor2descendant, nullptr);

  setParent (const_cast <DTNode*> (parent_new));
  // Arc-specific data
  len = len_new;
  attr = attr_new;
}




const Steiner* Steiner::makeDTRoot ()
{ 
  const Steiner* root = static_cast <const DTNode*> (getTree (). root) -> asSteiner ();
  ASSERT (root);
  makeRoot (const_cast <Steiner*> (root)); 
  ASSERT (this == getTree (). root);
  
  return root;
}




// Leaf

Leaf::Leaf (DistTree &tree,
  	        Steiner* parent_arg,
  	        Real len_arg,
  	        const string &name_arg)
: DTNode (tree, parent_arg, len_arg)  // DistTree must be declared
, name (name_arg)
{}



void Leaf::qc () const
{ 
	DTNode::qc ();
	  
  ASSERT (getParent ());
	ASSERT (isLeaf ());	
  ASSERT (! name. empty());
  ASSERT (! isLeft (name, "0x"));
  if (! isNan (len) && ! discernable && len != 0)
  {
    cout << getName () << " " << len << endl;
    ERROR;
  }
  if (getParent () && discernable != static_cast <const DTNode*> (getParent ()) -> childrenDiscernable ())
  {
    cout << getName () << " " << discernable << " " << getParent () -> getName () 
             << " " << static_cast <const DTNode*> (getParent ()) -> childrenDiscernable () << endl;
    ERROR;
  }
  
  IMPLY (reprLeaf, reprLeaf == this);
}



const DTNode* Leaf::getDiscernable () const
{ 
  if (const Leaf* leaf = inDiscernable ())
    return static_cast <const DTNode*> (leaf->getParent ());
  return this;
}



Real Leaf::getRelLenError () const
{ 
  return getLenError () / getDistTree (). getLenError ();
}



Steiner* Leaf::collapse (Leaf* other)
{
  ASSERT (other);
  ASSERT (this != other);
  ASSERT (other->graph == graph);
  ASSERT (! positive (len));
  
  discernable = false;
  len = 0;

  Steiner* otherParent = const_static_cast <Steiner*> (other->getParent ());
  ASSERT (otherParent);
  Steiner* st = other->discernable
                  ? new Steiner (const_cast <DistTree&> (getDistTree ()), otherParent, other->len)
                  : nullptr;
  if (st)
  {
    this ->setParent (st);
    other->setParent (st);
    st->subtreeLeaves = other->subtreeLeaves;
  }
  else
    setParent (otherParent);
  
  other->discernable = false;
  other->len = 0;
  
  return st;
}




// Change

Change::~Change ()
{
  ASSERT (status != eApplied);
}



void Change::qc () const
{
  Root::qc ();
    
	ASSERT (valid ());
	IMPLY (! isNan (improvement), positive (improvement));
	ASSERT (! targets. empty ());
}



bool Change::apply ()
{
  ASSERT (status == eInit);
  
  * const_cast <RealAttr1*> (tree. prediction_old) = * tree. prediction;
  fromLen = from->len;
  const bool ok = apply_ ();  
//IMPLY (tree. ds. objs. size () == tree. dissimSize_max (), ok);  
  status = eApplied;
  
  return ok;
}



void Change::restore ()
{
  ASSERT (status == eApplied);
  status = eInit;

  restore_ ();
  * const_cast <RealAttr1*> (tree. prediction) = * tree. prediction_old;
}



void Change::commit ()
{
  ASSERT (status == eApplied);
	status = eDone;

  commit_ ();
}



bool Change::compare (const Change* a, 
	                    const Change* b)
{ 
  ASSERT (a);
  ASSERT (! isNan (a->improvement));
  
  if (! positive (a->improvement))
    return false;
  
  if (a == b)  
    return false;
  if (! b)  
    return positive (a->improvement);
  ASSERT (positive (b->improvement));
    
  if (a->improvement > b->improvement)  return true;  
  if (a->improvement < b->improvement)  return false;
  const int cmp = strcmp (a->type (), b->type ());
  if (cmp < 0)  return true; 
  if (cmp > 0)  return false; 
  if (a < b)  return true;  // non-stable result ??
    
  return false;
}




// ChangeToSibling

bool ChangeToSibling::apply_ ()
{
  ASSERT (targets. size () == 2);
  ASSERT (! oldParent);
  ASSERT (! arcEnd);
  ASSERT (! inter);
  
  oldParent = const_static_cast <Steiner*> (from->getParent ());
  ASSERT (oldParent);
  // May be: oldParent == to
  arcEnd = const_static_cast <Steiner*> (to->getParent ());
    // May be nullptr
  ASSERT (from != arcEnd);

  // targets  
  targets << oldParent;
  if (arcEnd)
    targets << arcEnd;

  toLen = to->len;

  // Topology
  inter = new Steiner (const_cast <DistTree&> (tree), arcEnd, 0);
  const_cast <DTNode*> (from) -> setParent (inter); 
  const_cast <DTNode*> (to)   -> setParent (inter);

  
  // DTNode::len
  const_cast <DTNode*> (from) -> len = 0;
  const_cast <DTNode*> (to)   -> len = 0;
  ASSERT (inter->len == 0);

  // tree.{fromAttr_new,toAttr_new,interAttr,target_new}
  CompactBoolAttr1* fromAttr_new = const_cast <CompactBoolAttr1*> (tree. fromAttr_new);
  CompactBoolAttr1* toAttr_new   = const_cast <CompactBoolAttr1*> (tree. toAttr_new);
  CompactBoolAttr1* interAttr    = const_cast <CompactBoolAttr1*> (tree. interAttr);
  RealAttr1* target_new = const_cast <RealAttr1*> (tree. target_new);
  fromAttr_new->setAll (false);
  toAttr_new  ->setAll (false);
  interAttr   ->setAll (false);
  FOR (size_t, objNum, tree. ds. objs. size ())
  {
    const VectorPtr<Tree::TreeNode> path (tree. getPath ( tree. obj2leaf1 [objNum]
                                                        , tree. obj2leaf2 [objNum]
                                                        )
                                         );
    (*target_new) [objNum] = (* tree. target) [objNum] - tree. path2prediction (path);  
      // redundant <= map paths in the old tree to path in the new tree ??
    const bool fromUsed = path. contains (from);
    fromAttr_new->setCompactBool (objNum, fromUsed);
    if (fromUsed)
    {
      const bool toVia    = path. contains (to);
      const bool interVia = path. contains (inter);
      ASSERT (toVia != interVia);
      toAttr_new->setCompactBool (objNum, toVia);
      interAttr ->setCompactBool (objNum, interVia);
    }
    else
    {
      const bool toUsed = path. contains (to);
      ASSERT (toUsed == path. contains (inter));
      toAttr_new->setCompactBool (objNum, toUsed);
      interAttr ->setCompactBool (objNum, toUsed);
    }
  }

  Space1<NumAttr1> sp (tree. ds, false);  sp. reserve (3);
  sp << fromAttr_new;  // 0
  if (arcEnd)
    sp << toAttr_new   // 1
       << interAttr;   // 2
  L2LinearNumPrediction lr (tree. dsSample, sp, *target_new);
  lr. solveUnconstrained ();
  if (verbose ())
    lr. qc ();  
  if (isNan (lr. absCriterion))
    return false;

  FOR (size_t, i, lr. beta. size ())
    maximize (lr. beta [i], 0.0);

  if (arcEnd)
  {
    const_cast <DTNode*> (from) -> len = lr. beta [0];
    const_cast <DTNode*> (to)   -> len = lr. beta [1];
    inter                       -> len = lr. beta [2];
  }
  else
  {
    const_cast <DTNode*> (from) -> len = lr. beta [0] / 2;
    const_cast <DTNode*> (to)   -> len = from->len;
    inter                       -> len = NAN;
  }

  // tree.prediction
  FOR (size_t, objNum, tree. ds. objs. size ())
    (* const_cast <RealAttr1*> (tree. prediction)) [objNum] = (* tree. target) [objNum] - lr. getResidual (objNum);
    
  return true;
}



void ChangeToSibling::restore_ ()
{
  ASSERT (oldParent);
//ASSERT (arcEnd);
  ASSERT (inter);

  // targets  
  targets. pop ();
  if (arcEnd)
    targets. pop ();

  // DTNode::len
  const_cast <DTNode*> (from) -> len = fromLen;
  const_cast <DTNode*> (to)   -> len = toLen;

  // Topology
  const_cast <DTNode*> (from) -> setParent (oldParent);  
  const_cast <DTNode*> (to)   -> setParent (arcEnd);
  delete inter;
  
  oldParent = nullptr;
  arcEnd = nullptr;
  inter = nullptr;
}



void ChangeToSibling::commit_ ()
{
  ASSERT (oldParent);
//ASSERT (arcEnd);
  ASSERT (inter);
  
  inter->addAttr ();  
  if (oldParent->isTransient ())
    const_cast <DistTree&> (tree). delayDeleteRetainArcs (oldParent);
}




#if 0
// ChangeToChild

bool ChangeToChild::apply_ ()
{
  ASSERT (targets. size () == 2);
  ASSERT (! oldParent);
  ASSERT (! arcEnd);
  ASSERT (! inter);
  
  Steiner* from_ = const_cast <Steiner*> (from->asSteiner ());
  ASSERT (from_);
  
  oldParent = const_static_cast <Steiner*> (from->getParent ());
    // May be nullptr
  arcEnd = const_static_cast <Steiner*> (to->getParent ());
  ASSERT (arcEnd);
  // May be: arcEnd == from
  ASSERT (oldParent != arcEnd);

  // targets  
  if (oldParent)
    targets << oldParent;
  targets << arcEnd;

  toLen = to->len;

  // Topology
  inter = new Steiner (const_cast <DistTree&> (tree), arcEnd, 0);
  const_cast <DTNode*> (to) -> setParent (inter);
  ASSERT (! inter->attr);
  inter->makeRoot (from_);
  ASSERT (to->attr);
  ASSERT (inter->attr);
  ASSERT ((bool) from->attr == (from != arcEnd));
  
    
  // DTNode::len
  const_cast <DTNode*> (to) -> len = 0;
  ASSERT (arcEnd->len == 0);
  inter->len = 0;

  // tree.{fromAttr_new,toAttr_new,interAttr,target_new}
  CompactBoolAttr1* arcEnd_new = const_cast <CompactBoolAttr1*> (tree. fromAttr_new);
  CompactBoolAttr1* toAttr_new = const_cast <CompactBoolAttr1*> (tree. toAttr_new);
  CompactBoolAttr1* interAttr  = const_cast <CompactBoolAttr1*> (tree. interAttr);
  RealAttr1* target_new        = const_cast <RealAttr1*>        (tree. target_new);
  arcEnd_new->setAll (false);
  toAttr_new->setAll (false);
  interAttr ->setAll (false);
  FOR (size_t, objNum, tree. ds. objs. size ())
  {
    const VectorPtr<Tree::TreeNode> path (tree. getPath ( tree. obj2leaf1 [objNum]
                                                        , tree. obj2leaf2 [objNum]
                                                        )
                                         );
    (*target_new) [objNum] = (* tree. target) [objNum] - tree. path2prediction (path);
    const bool interUsed = path. contains (inter);
    interAttr->setCompactBool (objNum, interUsed);
    if (interUsed)
    {
      const bool toVia     = path. contains (to);
      const bool arcEndVia = path. contains (arcEnd);
      ASSERT (toVia != arcEndVia);
      toAttr_new->setCompactBool (objNum, toVia);
      arcEnd_new->setCompactBool (objNum, arcEndVia);
    }
    else
    {
      const bool toUsed = path. contains (to);
      ASSERT (toUsed == path. contains (arcEnd));
      toAttr_new->setCompactBool (objNum, toUsed);
      arcEnd_new->setCompactBool (objNum, toUsed);
    }
  }

  Space1<NumAttr1> sp (tree. ds, false);  sp. reserve (3);
  sp << toAttr_new;   // 0
  if (oldParent)
    sp << arcEnd_new  // 1
       << interAttr;  // 2
  L2LinearNumPrediction lr (tree. dsSample, sp, *target_new);
  lr. solveUnconstrained ();
  if (verbose ())
    lr. qc ();  
  if (isNan (lr. absCriterion))
    return false;

  FOR (size_t, i, lr. beta. size ())
    maximize (lr. beta [i], 0.0);

  if (oldParent)
  {
    const_cast <DTNode*> (to) -> len = lr. beta [0];
    arcEnd                    -> len = lr. beta [1];
    inter                     -> len = lr. beta [2];
  }
  else
  {
    const_cast <DTNode*> (to) -> len = lr. beta [0] / 2;
    arcEnd                    -> len = to->len;
    inter                     -> len = NAN;
  }

  // tree.prediction
  FOR (size_t, objNum, tree. ds. objs. size ())
    (* const_cast <RealAttr1*> (tree. prediction)) [objNum] = (* tree. target) [objNum] - lr. getResidual (objNum);

  return true;
}



void ChangeToChild::restore_ ()
{
  ASSERT (arcEnd);
  ASSERT (inter);

  if (oldParent)
    targets. pop ();
  targets. pop ();

  // DTNode::len
  const_cast <DTNode*> (to) -> len = toLen;  
  inter->len = fromLen;

  // Topology
  const_cast <Steiner*> (from->asSteiner ()) -> makeRoot (inter);  // from->len = inter->len
  const_cast <DTNode*> (to) -> setParent (arcEnd);
  ASSERT (! inter->attr);
  delete inter;
    
  ASSERT (sameReal (from->len, fromLen));
  
  oldParent = nullptr;
  arcEnd = nullptr;
  inter = nullptr;
}



void ChangeToChild::commit_ ()
{
  ASSERT (arcEnd);
  ASSERT (inter);

  arcEnd->addAttr ();
  if (from->isTransient ())
    const_cast <DistTree&> (tree). delayDeleteRetainArcs (const_cast <DTNode*> (from));
}




// Swap

bool Swap::apply_ ()
{
  Steiner* fromParent = const_static_cast <Steiner*> (from->getParent ());
  ASSERT (fromParent);
  Steiner* toParent = const_static_cast <Steiner*> (to->getParent ());
  ASSERT (toParent);

  toLen = to->len;

  // Topology
  const_cast <DTNode*> (from) -> setParent (toParent); 
  const_cast <DTNode*> (to)   -> setParent (fromParent); 


  // DTNode::len
  const_cast <DTNode*> (from) -> len = 0;
  const_cast <DTNode*> (to)   -> len = 0;

  // tree.{fromAttr_new,toAttr_new,target_new}
  CompactBoolAttr1* fromAttr_new = const_cast <CompactBoolAttr1*> (tree. fromAttr_new);
  CompactBoolAttr1* toAttr_new   = const_cast <CompactBoolAttr1*> (tree. toAttr_new);
  RealAttr1* target_new          = const_cast <RealAttr1*>        (tree. target_new);
  fromAttr_new->setAll (false);
  toAttr_new  ->setAll (false);
  FOR (size_t, objNum, tree. ds. objs. size ())
  {
    const VectorPtr<Tree::TreeNode> path (tree. getPath ( tree. obj2leaf1 [objNum]
                                                        , tree. obj2leaf2 [objNum]
                                                        )
                                         );
    fromAttr_new->setCompactBool (objNum, path. contains (from));
    toAttr_new  ->setCompactBool (objNum, path. contains (to));
    (*target_new) [objNum] = (* tree. target) [objNum] - tree. path2prediction (path);
  }

  Space1<NumAttr1> sp (tree. ds, false);  sp. reserve (2);
  sp << fromAttr_new   // 0
     << toAttr_new     // 1
     ;
  L2LinearNumPrediction lr (tree. dsSample, sp, *target_new);
  lr. solveUnconstrained ();
  if (verbose ())
    lr. qc ();  
  if (isNan (lr. absCriterion))
    return false;

  FOR (size_t, i, lr. beta. size ())
    maximize (lr. beta [i], 0.0);

  const_cast <DTNode*> (from) -> len = lr. beta [0];
  const_cast <DTNode*> (to)   -> len = lr. beta [1];

  // tree.prediction
  FOR (size_t, objNum, tree. ds. objs. size ())
    (* const_cast <RealAttr1*> (tree. prediction)) [objNum] = (* tree. target) [objNum] - lr. getResidual (objNum);

  return true;
}



void Swap::restore_ ()
{
  // DTNode::len
  const_cast <DTNode*> (from) -> len = fromLen;
  const_cast <DTNode*> (to)   -> len = toLen;

  // Topology
  Steiner* fromParent = const_static_cast <Steiner*> (to->getParent ());
  ASSERT (fromParent);
  Steiner* toParent = const_static_cast <Steiner*> (from->getParent ());
  ASSERT (toParent);
  const_cast <DTNode*> (from) -> setParent (fromParent); 
  const_cast <DTNode*> (to)   -> setParent (toParent); 
}
#endif




// DistTree

namespace
{
  // PAR
  // The greater then better DistTree::absCriterion
  constexpr uint areaRadius_std = 4;  // >= 4 <= ChangeToCousin can be applied  
}



DistTree::DistTree (const string &treeFName,
	                  const string &dissimFName,
	                  const string &attrName,
	                  bool sparse)
: dsSample (ds)
{
	ASSERT (dissimFName. empty () == attrName. empty ());

  // Initial tree topology
  loadTreeFile (treeFName);
  ASSERT (root);
  ASSERT (nodes. front () == root);
  ASSERT (static_cast <const DTNode*> (root) -> asSteiner ());

  setName2leaf ();
        
  if (dissimFName. empty ())
    setReprLeaves ();
  else
  {
  //Verbose verb;
    loadDissimDs (dissimFName, attrName);
    dissimDs2ds (sparse);  
    if (! getConnected ())
      throw runtime_error ("Disconnected objects");
  	setDiscernable ();  
  	
    topology2attrs (nodes);
    setPrediction ();
    setAbsCriterion (); 
  }
}



DistTree::DistTree (const string &dirName,
	                  const string &dissimFName,
	                  const string &attrName)
: dsSample (ds)
{
  // Initial tree topology
  loadTreeDir (dirName);
  ASSERT (root);
  ASSERT (nodes. front () == root);
  ASSERT (static_cast <const DTNode*> (root) -> asSteiner ());

  setName2leaf ();
        
  loadDissimDs (dissimFName, attrName);
  dissimDs2ds (false);  
  if (! getConnected ())
    throw runtime_error ("Disconnected objects");
	setDiscernable ();  
	
  setGlobalLen ();    

  topology2attrs (nodes);
  setPrediction ();
  setAbsCriterion (); 
}



DistTree::DistTree (const string &dissimFName,
	                  const string &attrName)
: dsSample (ds)
{
  loadDissimDs (dissimFName, attrName);
  ASSERT (dissimDs. get ());

  // Initial tree topology: star topology of indiscernable objects
  ASSERT (! root);
  auto root_ = new Steiner (*this, nullptr, NAN);
  for (const Obj* obj : dissimDs->objs)
	  new Leaf (*this, root_, 0, obj->name);
  ASSERT (root);
  ASSERT (nodes. front () == root);
  ASSERT (static_cast <const DTNode*> (root) -> asSteiner ());

  setName2leaf ();
        
  dissimDs2ds (false);  
  if (! getConnected ())
    throw runtime_error ("Disconnected objects");
	setDiscernable (); 
	 
  neighborJoin ();
  
  topology2attrs (nodes);
  setPrediction ();
  setAbsCriterion (); 
}



DistTree::DistTree (Prob branchProb,
                    size_t leafNum_max)
: dsSample (ds)
{
	ASSERT (isProb (branchProb));
	ASSERT (branchProb < 1);
	ASSERT (leafNum_max);


  const Real len = 1;  // PAR

  Set<Steiner*> open;
  open << new Steiner (*this, nullptr, NAN);
  Bernoulli bernoulli;
  bernoulli. setParam (branchProb);
  size_t leafIndex = 0;
  size_t leafNum = 1;  // isLeaf()
  while (! open. empty ())
  {
    auto it = open. begin ();
    ASSERT (it != open. end ());
    Steiner* st = *it;
    ASSERT (st);
    bernoulli. randVariable ();
    if (leafNum < leafNum_max && bernoulli. variable)
    {
      if (! st->isLeaf ())
        leafNum++;
      open << new Steiner (*this, st, len);
    }
    else
    {
      open. erase (it);
      if (st->isLeaf ())
      {
        leafIndex++;
	      new Leaf (*this, const_static_cast <Steiner*> (st->getParent ()), st->len, toString (leafIndex));
	      delete st;
      }
    }
  }	

  if (root->isTransient ())
    delayDeleteRetainArcs (const_static_cast <DTNode*> (root));
  deleteTransients ();  

	ASSERT (root);
  ASSERT (nodes. front () == root);
  if (! static_cast <const DTNode*> (root) -> asSteiner ())
    throw runtime_error ("One-node tree");

  setName2leaf ();
        
  setReprLeaves ();
}



DistTree::DistTree (const DTNode* center,
                    uint areaRadius,
                    VectorPtr<TreeNode> &area,
                    const DTNode* &area_root,
                    Node2Node &newLeaves2boundary)
: dsSample (ds)
{
  ASSERT (center);
  ASSERT (center->graph != this);
  ASSERT (! center->inDiscernable ());
  ASSERT (areaRadius >= 2);
  ASSERT (area. empty ());
  ASSERT (newLeaves2boundary. empty ());
  
  const DistTree& wholeTree = center->getDistTree ();
  
  // area
  VectorPtr<TreeNode> boundary;
  center->getArea (areaRadius, area, boundary);  
  ASSERT (! area. empty ());
  if (area. size () == 1)
    throw runtime_error ("Singleton tree");
  // Remove !discernable
  for (Iter<VectorPtr<TreeNode>> iter (area); iter. next (); )
    if (static_cast <const DTNode*> (*iter) -> inDiscernable ())
      iter. erase ();

  Set<const TreeNode*> areaSet;
  areaSet. insertAll (area);
  ASSERT (areaSet. contains (center));

  // boundary: remove !discernable
  {
    Set<const TreeNode*> newBoundary; 
    for (Iter<VectorPtr<TreeNode>> iter (boundary); iter. next ();)
    {
      const DTNode* node = static_cast <const DTNode*> (*iter);
      if (node->inDiscernable ())
      {
        iter. erase ();
        newBoundary << node->getParent ();
      }
    }
    ASSERT (areaSet. contains (newBoundary));
  #ifndef NDEBUG
    for (const TreeNode* node : boundary)
      ASSERT (! newBoundary. contains (node));
  #endif
    for (const TreeNode* node : newBoundary)
      boundary << node;    
  }
#ifndef NDEBUG
  for (const TreeNode* node : boundary)
    ASSERT (! static_cast <const DTNode*> (node) -> inDiscernable ());
#endif


  // nodes, area_root
  Node2Node old2new;  // 1-1
  area_root = nullptr;  
  size_t leafNum = 0;
  for (const TreeNode* node_old : area)
  {
    ASSERT (node_old);
    if (! areaSet. contains (node_old->getParent ()))
    {
      ASSERT (! area_root);
      area_root = static_cast <const DTNode*> (node_old);
    }
    Set<const TreeNode*> children_old;
    for (const DiGraph::Node* child : node_old->getChildren ())
      children_old << static_cast <const TreeNode*> (child);
    DTNode* node_new = nullptr;
    const Real len = static_cast <const DTNode*> (node_old) -> len;
    if (children_old. intersects (areaSet))
      node_new = new Steiner (*this, nullptr, len);
    else
    {
      leafNum++;
      node_new = new Leaf (*this, nullptr, len, "L" + toString (leafNum));
    }
    ASSERT (node_new);
    old2new [node_old] = node_new;
  }
  ASSERT (area_root);
  ASSERT (areaSet. contains (area_root));
  ASSERT (! areaSet. contains (area_root->getParent ()));
  ASSERT (nodes. size () == area. size ());
  
  borrowArcs (old2new, true);  // Arc's are not parallel

  setRoot ();
  ASSERT (static_cast <const DTNode*> (root) -> asSteiner ());
  ASSERT (findPtr (old2new, static_cast <const DiGraph::Node*> (area_root)) == root);
  
  {
    const VectorPtr<DiGraph::Node> children (root->getChildren ());
    ASSERT (! children. empty ());
    if (children. size () == 1)
    {
      // old2new[area_root] --> Leaf
      DTNode* child = const_static_cast <DTNode*> (children. front ());
      ASSERT (child);
      const TreeNode* root_ = root;                            // Old root in *this
      auto inter = new Steiner (*this, nullptr, NAN);  // New root in *this
      child->setParent (inter);
      child->len /= 2;
      ASSERT (root_->getChildren (). empty ());
      delete root_;
      leafNum++;  
      auto leaf = new Leaf (*this, inter, child->len, "L" + toString (leafNum));
      old2new [area_root] = leaf;  
      ASSERT (inter == root);
      ASSERT (nodes. size () == area. size () + 1/*inter*/);
      if (verbose ())
        cout << "area_root -> " << leaf->getName () << endl;
      ASSERT (static_cast <const DTNode*> (findPtr (old2new, static_cast <const DiGraph::Node*> (area_root))) -> asLeaf ());
    }
    else
    { 
      ASSERT (static_cast <const DTNode*> (findPtr (old2new, static_cast <const DiGraph::Node*> (area_root))) -> asSteiner ());
      ASSERT (area_root == wholeTree. root);
      ASSERT (children. size () == area_root->getChildren (). size ()); 
    }
  }
  ASSERT (isNan (static_cast <const DTNode*> (root) -> len));  

  
  setName2leaf ();

  Set<const TreeNode*> boundarySet;
  boundarySet. insertAll (boundary);
  ASSERT (boundarySet. size () == boundary. size ());
  ASSERT (boundarySet. size () == name2leaf. size ());

  // newLeaves2boundary
  {
    const Node2Node new2old (DiGraph::reverse (old2new));
    ASSERT (new2old. size () == old2new. size ());
    for (const DiGraph::Node* node : nodes)
    {
      const DiGraph::Node* node_old = findPtr (new2old, node);  // May be nullptr
      ASSERT ((bool) static_cast <const DTNode*> (node) -> asLeaf () == boundarySet. contains (static_cast <const TreeNode*> (node_old)));
      if (static_cast <const DTNode*> (node) -> asLeaf ())
      {
        ASSERT (node_old);
        newLeaves2boundary [node] = node_old;
      }
    }
  }
  ASSERT (newLeaves2boundary. size () == boundary. size ());
  

  // ds.objs[]->mult, *target, obj2leaf1[], obj2leaf2[]
  // For some leaf pairs the dissimilarity may be missing
  ds. objs. reserve (name2leaf. size () * (name2leaf. size () - 1) / 2);
  target = new RealAttr1 ("Target", ds, wholeTree. target->decimals); 
  obj2leaf1. reserve (ds. objs. capacity ());
  obj2leaf2. reserve (ds. objs. capacity ());  
  for (const auto it2 : name2leaf)
    for (const auto it1 : name2leaf)
    {
      if (it1. first == it2. first)
        break;
      ASSERT (it1. first < it2. first);
      const size_t objNum = ds. appendObj (getObjName (it1. first, it2. first));  
      const_cast <Obj*> (ds. objs [objNum]) -> mult = 0;
      (* const_cast <RealAttr1*> (target)) [objNum] = 0;
      obj2leaf1. resize (objNum + 1);
      obj2leaf2. resize (objNum + 1);
      obj2leaf1 [objNum] = it1. second;
      obj2leaf2 [objNum] = it2. second;
    }
  ds. setName2objNum ();


  const VectorOwn<Obj>& wholeObjs = wholeTree. ds. objs;
  FOR (size_t, objNum_whole, wholeObjs. size ())
  {
    if (wholeObjs [objNum_whole] -> mult == 0)
      continue;
      
    const VectorPtr<TreeNode> path (wholeTree. getPath ( wholeTree. obj2leaf1 [objNum_whole]
                                                       , wholeTree. obj2leaf2 [objNum_whole]
                                                       )
                                   );
    const Real dist_hat_whole = path2prediction (path);
    Set<const TreeNode*> pathSet (path);
    pathSet. intersect (areaSet);
    pathSet. erase (area_root);    
    if (pathSet. empty ())
      continue;

    VectorPtr<TreeNode> path_new;  path_new. reserve (pathSet. size ());
    insertAll (path_new, pathSet);
    const Real dist_hat_sub = path2prediction (path_new);

    Set<const TreeNode*> extremes (pathSet);
    extremes. intersect (boundarySet);
    ASSERT (! extremes. empty ());
    ASSERT (extremes. size () <= 2);
    if (extremes. size () == 1)
      extremes << area_root;
    ASSERT (extremes. size () == 2);
    
    const Real dist_whole = (* wholeTree. target) [objNum_whole] - (dist_hat_whole - dist_hat_sub);
      // May be < 0
    const Real mult_whole = wholeObjs [objNum_whole] -> mult;
    ASSERT (mult_whole > 0);
    
    // objNum  
    const DTNode* node0 = static_cast <const DTNode*> (findPtr (old2new, static_cast <const DiGraph::Node*> (extremes. front ())));
    const DTNode* node1 = static_cast <const DTNode*> (findPtr (old2new, static_cast <const DiGraph::Node*> (extremes. back ())));
    ASSERT (node0);
    ASSERT (node1);
    ASSERT (node0 != node1);
    const Leaf* leaf0 = node0->asLeaf ();
    const Leaf* leaf1 = node1->asLeaf ();
    ASSERT (leaf0);
    ASSERT (leaf1);
    const size_t objNum = ds. getName2objNum (getObjName (leaf0->name, leaf1->name));
    ASSERT (objNum != NO_INDEX);
    
    const_cast <Obj*> (ds. objs [objNum]) -> mult += mult_whole;
    (* const_cast <RealAttr1*> (target)) [objNum] += mult_whole * dist_whole;
  }


  // ds.objs[]->mult, *target: finish
  // dissim2_sum
  FOR_REV (size_t, objNum, ds. objs. size ())
  {
    const Real mult = ds. objs [objNum] -> mult;
    if (mult == 0)
      const_cast <RealAttr1*> (target) -> setMissing (objNum);
    else
      (* const_cast <RealAttr1*> (target)) [objNum] /= mult;
    dissim2_sum += mult * sqr ((*target) [objNum]);  // max (0.0, dist);
  }
    

  loadDissimFinish ();
//ASSERT (setDiscernable () == 0);
  topology2attrs (nodes);
  setPrediction ();
  setAbsCriterion (); 
}



void DistTree::loadTreeDir (const string &dir)
{
	ASSERT (! dir. empty ());
  ASSERT (isRight (dir, "/"));

  Vector<string> fileNames;  
  {
    const string outFName (dir + ".list");
    EXEC_ASSERT (system (("ls " + dir + " > " + outFName). c_str ()) == 0);
    LineInput f (outFName);
    fileNames = f. getVector ();
    EXEC_ASSERT (system (("rm " + outFName). c_str ()) == 0);
  }
  ASSERT (! fileNames. empty ());

  Name2steiner name2steiner;
  for (string name : fileNames)
  {
    EXEC_ASSERT (trimSuffix (name, dmSuff));    
    const Dataset leafDs (dir + name);
    if (leafDs. objs. empty ())
      continue;
    Steiner* steiner = getName2steiner (name, name2steiner);
    for (const Obj* obj : leafDs. objs)
      new Leaf (*this, steiner, NAN, obj->name);
  }

  deleteTransients ();  
}



Steiner* DistTree::getName2steiner (const string &name,
                                    Name2steiner &name2steiner) 
{
  if (name. empty ())
    return nullptr;
  if (Steiner* s = name2steiner [name])
    return s;

  string prefix (name);
  rfindSplit (prefix, '.');
  auto s = new Steiner (*this, getName2steiner (prefix, name2steiner), NAN);
  name2steiner [name] = s;
  
  return s;
}



void DistTree::loadTreeFile (const string &fName)
{
	ASSERT (! fName. empty ());

  Vector<string> lines;
  try
  {
    LineInput in (fName, 10000);  // PAR
    lines = in. getVector ();
  }
  catch (const exception &e)
  {
    throw runtime_error ("Loading tree file \"" + fName + "\"\n" + e. what ());
  }
  ASSERT (! lines. empty ());

	size_t lineNum = 0; 
  Steiner* steiner = nullptr;  	
	EXEC_ASSERT (loadLines (lines, lineNum, steiner, 0));
}



namespace 
{
	Real token2real (const string &s,
	                 const string &token)
	{
	  const string s1 (" " + s);
		const string token1 (" " + token + "=");
		const size_t pos = s1. find (token1);
		if (pos == string::npos)
			return NAN;  
		string valueS (s1. substr (pos + token1. size ()));
		return str2real (findSplit (valueS));
	}
}



bool DistTree::loadLines (const Vector<string>& lines,
						              size_t &lineNum,
						              Steiner* parent,
						              size_t expectedOffset)
{ 
	ASSERT (lineNum <= lines. size ());

	if (lineNum == lines. size ())
		return false;

	const string& line = lines [lineNum];
	string s (line);
	trimLeading (s);
	const size_t offset = line. size () - s. size ();
	if (offset < expectedOffset)
		return false;
	if (offset != expectedOffset)
	{
		cout << "Line " << lineNum + 1 << ": " << line << endl;
		ERROR;
	}
	
	lineNum++;

  string idS (findSplit (s));
	ASSERT (isRight (idS, ":"));
	idS. erase (idS. size () - 1);
  const Real len          = token2real (s, "len");
  const Real paths        = token2real (s, "paths");
  const Real errorDensity = token2real (s, "err_density");
  const Real leafError    = token2real (s, "leaf_error");
  DTNode* dtNode = nullptr;
	if (isLeft (idS, "0x"))
	{
	  ASSERT (isNan (leafError));
    Steiner* steinerParent = nullptr;
    if (parent)
    {
      steinerParent = const_cast <Steiner*> (parent->asSteiner ());
      ASSERT (steinerParent);
    }
		auto steiner = new Steiner (*this, steinerParent, len);
		while (loadLines (lines, lineNum, steiner, expectedOffset + Offset::delta))
			;
	  dtNode = steiner;
	}
	else
	{
	  ASSERT (parent);
		auto leaf = new Leaf (*this, parent, len, idS);
		leaf->relLenError = leafError;
		dtNode = leaf;
	}
	ASSERT (dtNode);
	
	if (! isNan (paths))
	{
  	dtNode->paths = (size_t) DM_sp::round (paths);
  	dtNode->errorDensity = errorDensity;
  }
	
	return true;
}



void DistTree::setName2leaf ()
{
  ASSERT (name2leaf. empty ());
  for (const DiGraph::Node* node : nodes)
    if (const Leaf* g = static_cast <const DTNode*> (node) -> asLeaf ())
      name2leaf [g->name] = g;
}



void DistTree::loadDissimDs (const string &dissimFName,
                             const string &attrName)
{
	ASSERT (! dissimFName. empty ());
	ASSERT (! dissimDs. get ());
	ASSERT (! dissimAttr);
	ASSERT (! optimizable ());


  // dissimAttr, dissimDs
  {
    dissimDs. reset (new Dataset (dissimFName));
    const Attr* attr = dissimDs->name2attr (attrName);
    ASSERT (attr);
    dissimAttr = attr->asPositiveAttr2 ();
  }
  ASSERT (dissimAttr);
  ASSERT (dissimAttr->matr. isSquare ());
  ASSERT (dissimAttr->matr. rowsSize (false) == dissimDs->objs. size ());
  {
    Real maxCorrection;
    size_t row_bad, col_bad;
    const_cast <PositiveAttr2*> (dissimAttr) -> matr. symmetrize (maxCorrection, row_bad, col_bad);
    if (maxCorrection > 2 * pow (10, - (Real) dissimAttr->decimals))
      cout << "maxCorrection = " << maxCorrection 
           << " at " << dissimDs->objs [row_bad] -> name 
           << ", "   << dissimDs->objs [col_bad] -> name 
           << endl;
  }
  dissimDs->setName2objNum ();
}



void DistTree::dissimDs2ds (bool sparse)
{
	ASSERT (dissimDs. get ());
	ASSERT (dissimAttr);
	ASSERT (! optimizable ());


  const Set<string> objNames (dissimDs->getObjNames ());
  restrictLeaves (objNames, true);
#ifndef NDEBUG
  {
    const Set<string> leafNames (name2leaf);
    ASSERT (objNames. contains (leafNames));
  }
#endif
  ASSERT (name2leaf. size () <= dissimDs->objs. size ());

  setReprLeaves ();
  
  // Leaf::comment
  for (const Obj* obj : dissimDs->objs)
  {
    const Leaf* g = findPtr (name2leaf, obj->name);
    IMPLY (! extraObjs (), g);
    if (g)
      const_cast <Leaf*> (g) -> comment = obj->comment;
  }

  Set<string> selected;
  if (sparse)
    for (const DiGraph::Node* node : nodes)
      if (const Leaf* leaf = static_cast <const DTNode*> (node) -> asLeaf ())
      {
        ASSERT (leaf == leaf->reprLeaf);
        const TreeNode* ancestor = leaf;
        while (ancestor)
        {
        	for (const DiGraph::Arc* arc : ancestor->arcs [false])
        	{
        	  const DTNode* child = static_cast <const DTNode*> (arc->node [false]);
        	  ASSERT (child->reprLeaf);
        	  if (child->reprLeaf != leaf)
          	  selected << getObjName (leaf->name, child->reprLeaf->name);
        	}
          ancestor = ancestor->getParent ();
        }
      }      

  // ds.objs, *target, obj2leaf1, obj2leaf2, dissim2_sum
  if (verbose ())
    cout << "Leaf pairs -> data objects ..." << endl;
  ds. objs. reserve (name2leaf. size () * (name2leaf. size () - 1) / 2);
  target = new RealAttr1 ("Target", ds, dissimAttr->decimals); 
  obj2leaf1. reserve (ds. objs. capacity ());
  obj2leaf2. reserve (ds. objs. capacity ());
  FOR (size_t, row, dissimDs->objs. size ())
  {
    const string name1 = dissimDs->objs [row] -> name;
    if (! findPtr (name2leaf, name1))
      continue;
    FOR (size_t, col, row)  // dissimAttr is symmetric
    {
      string name2 = dissimDs->objs [col] -> name;
      if (! findPtr (name2leaf, name2))
        continue;
      const Real dissim = dissimAttr->get (row, col);
      if (positive (dissim) && sparse && ! selected. contains (getObjName (name1, name2)))
        continue;
      addDissim (name1, name2, dissim);
    }
  }
  
  if (! extraObjs ())
  {
    dissimDs. reset (nullptr);
    dissimAttr = nullptr;
  }


  loadDissimFinish ();      
}



bool DistTree::addDissim (const string &name1,
                          const string &name2,
                          Real dissim)
{
  if (   isNan (dissim)  // prediction must be large ??
      || ! DM_sp::finite (dissim)  
     )
  {
    cout << name1 << '-' << name2 << ": " << dissim << endl;
    return false;  
  }
    
  const Real mult = dissim2mult (dissim);
  if (nullReal (mult))
  {
    if (nullReal (dissim))
      ;  // Needed for setDiscernable();
    else
    {
      cout << name1 << '-' << name2 << ": " << dissim << endl;
      return false;
    }
  }
  
  dissim2_sum += mult * sqr (dissim);  // max (0.0, dissim);
  
  const size_t objNum = ds. appendObj (getObjName (name1, name2));
  
  const_cast <Obj*> (ds. objs [objNum]) -> mult = mult;
  (* const_cast <RealAttr1*> (target)) [objNum] = dissim;

  obj2leaf1. resize (objNum + 1);
  obj2leaf2. resize (objNum + 1);
  EXEC_ASSERT (obj2leaf1 [objNum] = findPtr (name2leaf, name1));
  EXEC_ASSERT (obj2leaf2 [objNum] = findPtr (name2leaf, name2));
  if (name1 > name2)
    swap ( obj2leaf1 [objNum]
         , obj2leaf2 [objNum]
         );
  
  return true;
}



namespace
{
  struct Neighbors
  {
    // !nullptr, discernable
    array <const DTNode*, 2> nodes;
      // nodes[0] < nodes[1]
    Real dissim {NAN};
      // !isNan()
    Real mult {NAN};  

    Neighbors ()
      { nodes [0] = nullptr;
        nodes [1] = nullptr;
      }
    Neighbors (const Leaf* leaf1,
               const Leaf* leaf2,
               Real dissim_arg,
               Real mult_arg)
      : dissim (dissim_arg)
      , mult (mult_arg)
      { ASSERT (leaf1);
        ASSERT (leaf2);
        nodes [0] = leaf1->getDiscernable ();
        nodes [1] = leaf2->getDiscernable ();
        orderNodes ();
      }
      
    void print (ostream &os) const
      { os << nodes [0] << ' ' << nodes [1] << ' ' << dissim << ' ' << mult << endl; }
    bool operator== (const Neighbors& other) const
      { return    nodes [0] == other. nodes [0]
               && nodes [1] == other. nodes [1];
      }
    void orderNodes ()
      { swapGreater (nodes [0], nodes [1]); }
    bool same () const
      { return nodes [0] == nodes [1]; }
    bool merge (const Neighbors &from)
      { if (! (*this == from))
          return false;
        dissim = (dissim * mult + from. dissim * from. mult) / (mult + from. mult);
        mult += from. mult;
        return true;
      }
    void setNodeDissimAve (bool first,
                           bool add)
      { const_cast <DTNode*> (nodes [first]) -> subtreeLen. add (dissim, mult * getSign (add)); }
      // Requires: !same()
    Real getParentDissim (bool first) const
      { return 0.5 * max (0.0, dissim + (  nodes [0] -> subtreeLen. getMean () 
                                         - nodes [1] -> subtreeLen. getMean ()
                                        ) 
                                        * getSign (first)
                         ); 
      }
    static bool compare (const Neighbors& n1,
                         const Neighbors& n2)
      { LESS_PART (n1, n2, nodes [0]);
        LESS_PART (n1, n2, nodes [1]);
        return false;  
      }
  };
}



void DistTree::neighborJoin ()
{
  ASSERT (ds. objs. size () >= 2);
  
  cout << "Neighbor joining..." << endl;
  
  // DTNode:subtreeLen: dissimilarity from other leaves

  clearSubtreeLen ();  
    
  Vector<Neighbors> neighborsVec;  neighborsVec. reserve (ds. objs. size ());
  FOR (size_t, objNum, ds. objs. size ())
  {
    Neighbors neighbors ( obj2leaf1 [objNum]
                        , obj2leaf2 [objNum]
                        , (*target) [objNum]
                        , 1  // better than Obj::mult ??
                        );
    if (neighbors. same ())
      continue;
    ASSERT (positive (neighbors. dissim));
    ASSERT (positive (neighbors. mult));
    for (const bool first : Bool)
      neighbors. setNodeDissimAve (first, true);
    neighborsVec << neighbors;
  }
  
  if (neighborsVec. empty ())
    return;
  

  Progress prog ((uint) nodes. size ());
  Neighbors neighbors_best;
  for (;;)
  {
    prog ();
    
    // Remove duplicate Neighbors 
    {
      Common_sp::sort (neighborsVec, Neighbors::compare);
      size_t toRemove = 0;
      FOR (size_t, i, neighborsVec. size ())
      {
        const size_t j = i - toRemove;
        if (j != i)
          neighborsVec [j] = neighborsVec [i];
        if (   neighborsVec [j] == neighbors_best 
            || (j > 0 && neighborsVec [j - 1]. merge (neighborsVec [j]))
           )
          toRemove++;
      }
      while (toRemove)
      {
        neighborsVec. pop_back ();
        toRemove--;
      }
    }
    
    if (neighborsVec. size () == 1)
      break;
      
    bool first_best = false;
    Real dissim_min = INF;
    {
      size_t i_best = NO_INDEX;
      FOR (size_t, i, neighborsVec. size ())
        for (const bool first : Bool)
          if (minimize (dissim_min, neighborsVec [i]. getParentDissim (first)))
          {
            i_best = i;
            first_best = first;
          }
      ASSERT (i_best != NO_INDEX);
      neighbors_best = neighborsVec [i_best];
    }
    ASSERT (dissim_min >= 0);
    ASSERT (dissim_min < neighbors_best. dissim);
    
  #if 0
    cout << neighborsVec. size () << " " << nodes. size () << " " << dissim_min << endl;  
    if (nodes. size () > 600)
    {
      neighbors_best. print (cout);
      cout << endl;
      {
        OFStream of ("", "neighbors", "");
        for (const auto& n : neighborsVec)
          n. print (of);
      }
      ERROR;
    }
  #endif

    auto newNode = new Steiner (*this, const_static_cast <Steiner*> (neighbors_best. nodes [first_best] -> getParent ()), 0);
    {
      DTNode* a = const_cast <DTNode*> (neighbors_best. nodes [0]);
      DTNode* b = const_cast <DTNode*> (neighbors_best. nodes [1]);
      if (! first_best)
        swap (a, b);
      a->setParent (newNode);
      b->setParent (newNode);
      a->len = dissim_min;
      b->len = neighbors_best. dissim - dissim_min;
      a->subtreeLen. addValue (- a->len);
      b->subtreeLen. addValue (- b->len);
      newNode->subtreeLen = a->subtreeLen;
      newNode->subtreeLen. add (b->subtreeLen);
    }
    
    // neighborsVec
    for (Neighbors& neighbors : neighborsVec)
      if (! (neighbors == neighbors_best))
        for (const bool first : Bool)
        {
          bool found = false;
          for (const bool best_first : Bool)
            if (neighbors. nodes [first] == neighbors_best. nodes [best_first])
            {
              neighbors. nodes [first] = newNode;
              neighbors. setNodeDissimAve (! first, false);
              neighbors. dissim -= neighbors_best. nodes [best_first] -> len;
              maximize (neighbors. dissim, 0.0);
              neighbors. setNodeDissimAve (! first, true);
              found = true;
            }
          if (found)
            neighbors. orderNodes ();
          ASSERT (! neighbors. same ());
        }
  }
  ASSERT (neighborsVec. size () == 1);
  
  
  Neighbors& neighbors = neighborsVec [0];
  ASSERT (! neighbors. same ());
  for (const bool first : Bool)  
    const_cast <DTNode*> (neighbors. nodes [first]) -> len = neighbors. dissim / 2;  

  clearSubtreeLen ();  
    
  for (DiGraph::Node* node : nodes)
  {
    DTNode* dtNode = static_cast <DTNode*> (node);
    if (! dtNode->attr)
      dtNode->addAttr ();
  }
  
  finishChanges ();
}



void DistTree::loadDissimFinish ()
{
  ASSERT (! prediction);
  ASSERT (! fromAttr_new);
  ASSERT (! toAttr_new);
  ASSERT (! interAttr);
  ASSERT (! target_new);
  ASSERT (! prediction_old);
  

  dsSample = Sample (ds);
  absCriterion_delta = dsSample. multSum * 1e-5;  // PAR

  // ds.attrs, DTNode::attr
  if (verbose ())
    cout << "Arcs -> data attributes ..." << endl;
  {
    Progress prog ((uint) nodes. size (), (uint) verbose () * 100);  // PAR
    for (const DiGraph::Node* node : nodes)
    {
      prog ();
      const_static_cast <DTNode*> (node) -> addAttr ();
    }
  }
  
  prediction = new RealAttr1 ("Prediction", ds, target->decimals); 
  
  fromAttr_new = new CompactBoolAttr1 ("from_new", ds);
  toAttr_new   = new CompactBoolAttr1 ("to_new", ds);
  interAttr    = new CompactBoolAttr1 ("inter", ds);
  const_cast <CompactBoolAttr1*> (fromAttr_new) -> setAll (false);
  const_cast <CompactBoolAttr1*> (toAttr_new)   -> setAll (false);
  const_cast <CompactBoolAttr1*> (interAttr)    -> setAll (false);

  target_new     = new RealAttr1 ("target_new",    ds, target->decimals);
  prediction_old = new RealAttr1 ("PredictionOld", ds, target->decimals); 


	ASSERT (optimizable ());
}



size_t DistTree::setDiscernable ()
{
  ASSERT (optimizable ());


  map <const DisjointCluster*, VectorPtr<Leaf>> cluster2leaves;
  {
    // Leaf::DisjointCluster
   	for (DiGraph::Node* node : nodes)
   	{
   	  const DTNode* dtNode = static_cast <DTNode*> (node);
   	  if (Leaf* leaf = const_cast <Leaf*> (dtNode->asLeaf ()))
   	  {
   	    leaf->discernable = true;
   	    leaf->DisjointCluster::init ();
   	  }
   	}
    FOR (size_t, objNum, ds. objs. size ())
      if (! positive ((*target) [objNum]))
      {
        ASSERT (! target->isMissing (objNum));
        const_cast <Leaf*> (obj2leaf1 [objNum]) -> merge (* const_cast <Leaf*> (obj2leaf2 [objNum]));
      }
    // cluster2leaves
   	for (DiGraph::Node* node : nodes)
   	{
   	  const DTNode* dtNode = static_cast <DTNode*> (node);
   	  if (const Leaf* leaf = dtNode->asLeaf ())
        cluster2leaves [const_cast <Leaf*> (leaf) -> getDisjointCluster ()] << leaf;
   	}
  }
 
  size_t n = 0;   
  for (const auto it : cluster2leaves)
  {
    const VectorPtr<Leaf>& clusterNodes = it. second;
    ASSERT (! clusterNodes. empty ());
    if (clusterNodes. size () == 1)
      continue;
    const Leaf* first = clusterNodes [0];
    ASSERT (first);
    Steiner* parent = const_cast <Steiner*> (static_cast <const DTNode*> (first->getParent ()) -> asSteiner ());
    ASSERT (parent);
    auto steiner = new Steiner (*this, parent, 0);
    steiner->addAttr ();  
    for (const Leaf* leaf_ : clusterNodes)
    {
      Leaf* leaf = const_cast <Leaf*> (leaf_);
      leaf->setParent (steiner);
      leaf->discernable = false;  
      leaf->len = 0;  
      n++;
    }
    ASSERT (! steiner->isTransient ());
    ASSERT (n);
  }
  
  
  if (n)
  {
   	Vector<DiGraph::Node*> nodeVec;  nodeVec. reserve (nodes. size ());
  
   	// delete isLeaf() && Steiner
   	insertAll (nodeVec, nodes);
   	for (DiGraph::Node* node : nodeVec)  
  	  if (const Steiner* s = static_cast <DTNode*> (node) ->asSteiner ())
  	    while (s && s->isLeaf ())
  	    {
  	      ASSERT (s->graph); 
     		  ASSERT (! s->inDiscernable ());
     		  ASSERT (s->childrenDiscernable ());  // Actually no children
     		  const Steiner* parent = static_cast <const DTNode*> (s->getParent ()) -> asSteiner ();
     			delayDeleteRetainArcs (const_cast <Steiner*> (s));
     			s = parent;
     	  }
  
  	// delete isTransient()
   	nodeVec. clear ();
   	insertAll (nodeVec, nodes);
   	for (DiGraph::Node* node : nodeVec)  
   	{
   	  const DTNode* dtNode = static_cast <DTNode*> (node);
   		if (dtNode->isTransient ())
   		{
   		  const Steiner* s = dtNode->asSteiner ();
   		  ASSERT (s);
   		  ASSERT (s->childrenDiscernable ());
   		  ASSERT (s->arcs [false]. size () == 1);
   			delayDeleteRetainArcs (const_cast <Steiner*> (s));
   	  }
   	}
  
    toDelete. deleteData ();
  }

  
  return n;
}



bool DistTree::getConnected () 
{
  ASSERT (optimizable ());


  map <const DisjointCluster*, VectorPtr<Leaf>> cluster2leaves;
  {
    // Leaf::DisjointCluster
   	for (DiGraph::Node* node : nodes)
   	{
   	  const DTNode* dtNode = static_cast <DTNode*> (node);
   	  if (Leaf* leaf = const_cast <Leaf*> (dtNode->asLeaf ()))
   	    leaf->DisjointCluster::init ();
   	}
    FOR (size_t, objNum, ds. objs. size ())
      if (! nullReal (ds. objs [objNum] -> mult))
        const_cast <Leaf*> (obj2leaf1 [objNum]) -> merge (* const_cast <Leaf*> (obj2leaf2 [objNum]));
    // cluster2leaves
   	for (DiGraph::Node* node : nodes)
   	{
   	  const DTNode* dtNode = static_cast <DTNode*> (node);
   	  if (const Leaf* leaf = dtNode->asLeaf ())
        cluster2leaves [const_cast <Leaf*> (leaf) -> getDisjointCluster ()] << leaf;
   	}
  }
  ASSERT (! cluster2leaves. empty ());
  
  if (cluster2leaves. size () == 1)
    return true;
    
  const DisjointCluster* cluster_main = nullptr;
  size_t size_max = 0;
  for (const auto it : cluster2leaves)
    if (maximize (size_max, it. second. size ()))
      cluster_main = it. first;
  ASSERT (cluster_main);
 
  for (const auto it : cluster2leaves)
    if (it. first != cluster_main)
    {
      cout << endl;
      cout << "Cluster:" << endl;
      for (const Leaf* leaf : it. second)
        cout << leaf->name << endl;
    }

  return false;
}



void DistTree::qc () const
{ 
	Tree::qc ();

  ASSERT (nodes. size () >= 2);
		
	ASSERT (root);
	ASSERT (root->graph == this);
	const Steiner* root_ = static_cast <const DTNode*> (root) -> asSteiner ();
  ASSERT (root_);

 	for (DiGraph::Node* node : nodes)
 	{
 	  const DTNode* dtNode = static_cast <DTNode*> (node);
 		ASSERT (! dtNode->isTransient ());
 	}

  const size_t leaves = root->getLeavesSize ();
  ASSERT (name2leaf. size () == leaves);
  
  if (! optimizable ())
    return;
    
  ASSERT ((bool) dissimDs. get () == (bool) dissimAttr);
  if (dissimAttr)
  {
    dissimDs->qc ();
    ASSERT (& dissimAttr->ds == dissimDs. get ());
  }

  ds. qc ();
  ASSERT (ds. objs. size () <= dissimSize_max ());
	ASSERT (dsSample. ds == & ds);
	ASSERT (dsSample. size () == ds. objs. size ());
 	ASSERT (positive (dsSample. multSum));
 	
  
  ASSERT (target);
  ASSERT (prediction); 	
 	ASSERT (& target    ->ds == & ds);
 	ASSERT (& prediction->ds == & ds);
 	ASSERT (! target    ->existsMissing ());
 	
  ASSERT (obj2leaf1. size () == obj2leaf2. size ());
  ASSERT (obj2leaf1. size () == ds. objs. size ());
  Set<pair<const Leaf*,const Leaf*> > leafSet;
  FOR (size_t, objNum, ds. objs. size ())
  {
    const Leaf* g1 = obj2leaf1 [objNum];
    const Leaf* g2 = obj2leaf2 [objNum];
    ASSERT (g1);
    ASSERT (g2);
    ASSERT (g1->name < g2->name);
    leafSet. checkUnique (pair<const Leaf*,const Leaf*> (g1, g2));
  //ASSERT ((*target) [objNum] >= 0);
    IMPLY (ds. objs [objNum] -> mult && (*target) [objNum] == 0, 
           ! g1->discernable && ! g2->discernable && g1->getParent () == g2->getParent ()
          );
    if (! isNan ((*prediction) [objNum]))
    {
      ASSERT ((*prediction) [objNum] >= 0);
    /*
      if (topologyOptimized && nullReal ((*prediction) [objNum]) && ! nullReal ((*target) [objNum]))
      {
        cout << g1->name << " " << g2->name << endl;
        ERROR;
      }
    */
    }
  }

  ASSERT (ds. objs. size () <= (leaves * (leaves - 1)) / 2);

 	ASSERT (fromAttr_new);
 	ASSERT (toAttr_new);
 	ASSERT (interAttr);
 	ASSERT (target_new);
 	ASSERT (prediction_old);
 	
 	if (! isNan (absCriterion))
 	{
   	ASSERT (absCriterion >= 0);
   	ASSERT (eqReal (absCriterion, getAbsCriterion (), absCriterion_delta));
 	  ASSERT (! prediction->existsMissing ());
 	}

  ASSERT (absCriterion_delta > 0);
}



void DistTree::qcAttrs () const
{
  Set<const Attr*> nodeAttrs;
 	for (DiGraph::Node* node : nodes)
 	{
 	  const DTNode* dtNode = static_cast <DTNode*> (node);
 		if (dtNode->attr)
 		{
 		  nodeAttrs. checkUnique (dtNode->attr);
 		//const Sample sample (ds);
 		  if (dtNode != root && ! positive (dtNode->attr->getProb (dsSample)))
 		  {
 		    cout << dtNode->getName () << endl;
 		    cout << endl;
 		    FOR (size_t, objNum, ds. objs. size ())
 		      cout << ds. objs [objNum] -> name << '\t' << ds. objs [objNum] -> mult << '\t' << (*target) [objNum] << endl;
 		    ERROR;
 		  }
 		}
 	}

  Set<const Attr*> pathAttrs;
  for (const Attr* attr : ds. attrs)
    if (   attr != target
        && attr != prediction
        && attr != fromAttr_new
        && attr != toAttr_new
        && attr != interAttr
        && attr != target_new
        && attr != prediction_old
       )
      pathAttrs. checkUnique (attr);
  ASSERT (nodeAttrs == pathAttrs);

  FOR (size_t, objNum, ds. objs. size ())
  {
    const Leaf* g1 = obj2leaf1 [objNum];
    const Leaf* g2 = obj2leaf2 [objNum];
    ASSERT (g1);
    ASSERT (g2);
    Unverbose unv;
    if (verbose ())
    {
      const VectorPtr<TreeNode> path (getPath (g1, g2));
      for (const DiGraph::Node* node : nodes)
      {
        const DTNode* dtNode = static_cast <const DTNode*> (node);
        ASSERT ((*dtNode->attr) [objNum] == path. contains (dtNode));
      }
    }
  }

 	if (! isNan (absCriterion))
    for (const Attr* attr : ds. attrs)
      if (   attr != fromAttr_new
          && attr != toAttr_new
          && attr != interAttr
          && attr != target_new
          && attr != prediction_old
         )
      {
        ASSERT (! attr->existsMissing ());
      }
}



void DistTree::deleteLeaf (TreeNode* node,
                           bool deleteTransientAncestor)
{
  ASSERT (! optimizable ());
  
  ASSERT (node);
  ASSERT (& node->getTree () == this);
  ASSERT (node != root);
  
  const Leaf* leaf = static_cast <DTNode*> (node) -> asLeaf ();
  ASSERT (leaf);

  if (verbose ())
    cout << "Deleting: " << leaf->name << endl;
  const Steiner* parent = static_cast <const DTNode*> (leaf->getParent ()) -> asSteiner ();
  ASSERT (parent);
  delayDeleteRetainArcs (const_cast <Leaf*> (leaf));
  if (deleteTransientAncestor && parent->isTransient ())
    delayDeleteRetainArcs (const_cast <Steiner*> (parent));
	
  toDelete. deleteData ();
}



string DistTree::getObjName (const string &name1,
                             const string &name2)
{
  ASSERT (! name1. empty ());
  ASSERT (! name2. empty ());
  ASSERT (name1 != name2);  
  const string* p1 = & name1;
  const string* p2 = & name2;
  if (name1 > name2)
    swap (p1, p2);        
  return *p1 + "-" + *p2;
}



const DTNode* DistTree::lcaName2node (const string &lcaName) const
{
  ASSERT (! lcaName. empty ());

  string s (lcaName);
  const string name1 = findSplit (s, '-');
  const Leaf* leaf1 = findPtr (name2leaf, name1);
  ASSERT (leaf1);
  const Leaf* leaf2 = findPtr (name2leaf, s);
  ASSERT (leaf2);
  const DTNode* node = static_cast<const DTNode*> (getLowestCommonAncestor (leaf1, leaf2));
  ASSERT (node);
  
  return node;
}



Set<const DTNode*> DistTree::getDiscernables () const
{
  Set<const DTNode*> s;
  for (const DiGraph::Node* node : nodes)
    if (const Leaf* leaf = static_cast <const DTNode*> (node) -> asLeaf ())
      s << leaf->getDiscernable ();
  return s;
}

  

void DistTree::printInput (ostream &os) const
{
  os << "# Objects: " << root->getLeavesSize () << endl;
  os << "# Discernable objects: " << getDiscernables (). size () << endl;
  os << "# Nodes: " << nodes. size () << endl;
  if (! optimizable ())
    return;
  os << "# Dissimilarities: " << ds. objs. size () << " (" << (Real) ds. objs. size () / (Real) dissimSize_max () * 100 << " %)" << endl; 
  os << "# Binary attributes: " << ds. attrs. size () - 2 - 5 << endl;  // less target, prediction, for Change
  os << "Dissimilarity variance: " << varianceTypeNames [varianceType] << endl;
  os << "# Objects weighted: " << dsSample. multSum << endl;
  os << "# Dissimilarity average: " << getDissim_ave () << endl;
}



Real DistTree::path2prediction (const VectorPtr<TreeNode> &path) 
{
  Real dHat = 0;
  for (const TreeNode* node : path)
    dHat += static_cast <const DTNode*> (node) -> len;
  ASSERT (! isNan (dHat));
  ASSERT (dHat >= 0);
  return dHat;
}



void DistTree::saveFeatureTree (const string &fName) const
{
  if (fName. empty ())
    return;
  
  OFStream f ("", fName, "");
  static_cast <const DTNode*> (root) -> saveFeatureTree (f, 0);
}



void DistTree::topology2attrs (const List<DiGraph::Node*>& nodes_arg)
{
  ASSERT (optimizable ());

  if (verbose ())
    cout << "Data values ..." << endl;

  Progress prog ((uint) ds. objs. size (), (uint) verbose () * 10000);  // PAR
  FOR (size_t, objNum, ds. objs. size ())
  {
    prog ();
    const VectorPtr<TreeNode> path (getPath ( obj2leaf1 [objNum]
                                            , obj2leaf2 [objNum]
                                            )
                                   );
    for (const DiGraph::Node* node : nodes_arg)
    {
      const DTNode* dtNode = static_cast <const DTNode*> (node);
      CompactBoolAttr1* attr = const_cast <CompactBoolAttr1*> (dtNode->attr);
      ASSERT (attr);
      attr->setCompactBool (objNum, path. contains (dtNode));
    }
  }
}



void DistTree::clearSubtreeLen ()
{
  for (DiGraph::Node* node : nodes)
    static_cast <DTNode*> (node) -> subtreeLen. clear ();
}



void DistTree::setGlobalLen ()
{
  // DTNode::subtreeLen
  for (DiGraph::Node* node : nodes)
  {
    DTNode* dtNode = static_cast <DTNode*> (node);
    dtNode->subtreeLen. clear ();
    if (Leaf* leaf = const_cast <Leaf*> (dtNode->asLeaf ()))
      leaf->subtreeLen. add (0);  
  }
  FOR (size_t, objNum, obj2leaf1. size ())
  {
    const Leaf* g1 = obj2leaf1 [objNum];
    const Leaf* g2 = obj2leaf2 [objNum];
    const Real d = max (0.0, (*target) [objNum]);
    const TreeNode* ancestor = getLowestCommonAncestor (g1, g2);
    ASSERT (ancestor);
    Steiner* s = const_cast <Steiner*> (static_cast <const DTNode*> (ancestor) -> asSteiner ());
    ASSERT (s);
    s->subtreeLen. add (d / 2);
  }

  // DTNode::len
  for (DiGraph::Node* node : nodes)
  {
    DTNode* dtNode = static_cast <DTNode*> (node);
    if (dtNode->inDiscernable ())
      { ASSERT (dtNode->len == 0); }
    else
    {
      if (const DTNode* parent = static_cast <const DTNode*> (dtNode->getParent ()))
        dtNode->len = max (0.0, parent->subtreeLen. getMean () - dtNode->subtreeLen. getMean ());
    }
  }
}



void DistTree::setPrediction () 
{
  FOR (size_t, objNum, obj2leaf1. size ())
  {
    const VectorPtr<TreeNode> path (getPath ( obj2leaf1 [objNum]
                                            , obj2leaf2 [objNum]
                                            )
                                   );
    (* const_cast <RealAttr1*> (prediction)) [objNum] = path2prediction (path);  
  }
}



void DistTree::checkPrediction () const
{
  ASSERT (optimizable ());

  VectorPtr<DTNode> dtNodes;  dtNodes. reserve (nodes. size ());
  for (DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <DTNode*> (node);
    if (dtNode != root)
      dtNodes << dtNode;
  }

  Space1<NumAttr1> sp (ds, false);  sp. reserve (dtNodes. size ());
  for (const DTNode* dtNode : dtNodes)
    sp << dtNode->attr;
  
  // lr
  L2LinearNumPrediction lr (dsSample, sp, *target);
  ASSERT (lr. beta. size () == dtNodes. size ());
  FOR (size_t, attrNum, dtNodes. size ())
    lr. beta [attrNum] = dtNodes [attrNum] -> len;

  FOR (size_t, objNum, obj2leaf1. size ())
    if (! eqReal ((* const_cast <RealAttr1*> (prediction)) [objNum], lr. predict (objNum), getDissim_ave () * 1e-2))  // PAR 
    {
      ONumber on (cout, 6, true);
      cout << (* const_cast <RealAttr1*> (prediction)) [objNum] << " " << lr. predict (objNum) << "  objNum = " << objNum << endl;
      cout << obj2leaf1 [objNum] -> getName () << " " << obj2leaf2 [objNum] -> getName () << endl;
      cout << "Path: ";
      FOR (size_t, attrNum, dtNodes. size ())
        if ((* dtNodes [attrNum] -> attr) [objNum])
          cout << " " << dtNodes [attrNum] -> getName ();
      cout << endl;
      print (cout);
      ERROR;
    }
}



Real DistTree::getAbsCriterion () const
{
  ASSERT (optimizable ());

  Real absCriterion_ = 0;
  FOR (size_t, objNum, ds. objs. size ())
    if (const Real mult = ds. objs [objNum] -> mult)
    {
      const Real d = (*target) [objNum];
    //ASSERT (d >= 0);
      const Real dHat = (*prediction) [objNum];
      ASSERT (! isNan (dHat));
      absCriterion_ += mult * sqr (dHat - d);
    }
  ASSERT (DM_sp::finite (absCriterion_));
  
  return absCriterion_;
}



void DistTree::checkAbsCriterion (const string &title) const
{
  cout << "checkAbsCriterion: " << title << endl;
  qc ();
	const Real absCriterion_new = getAbsCriterion ();
	if (! eqReal (absCriterion, absCriterion_new, absCriterion_delta)) 
	{
		cout << absCriterion << " " << absCriterion_new << endl;
		ERROR;
	}
	ASSERT (leReal (absCriterion, dissim2_sum));
}



void DistTree::printAbsCriterion_halves () const
{
  ASSERT (optimizable ());
  
  const Real dissim_ave = getDissim_ave ();
  size_t pairs = 0;
  size_t size_half       [2/*bool*/] = {0, 0};
  Real absCriterion_half [2/*bool*/] = {0, 0};
  Real dissim2_half      [2/*bool*/] = {0, 0};
  FOR (size_t, objNum, ds. objs. size ())
    if (const Real mult = ds. objs [objNum] -> mult)
    {
      const Real d = (*target) [objNum];
      const Real dHat = (*prediction) [objNum];
      ASSERT (! isNan (dHat));
      const bool half2 = d > dissim_ave;
      pairs++;
      size_half [half2] ++;
      absCriterion_half [half2] += mult * sqr (dHat - d);
      dissim2_half        [half2] += mult * sqr (d);
    }
    
  for (const bool half2 : Bool)
  {
    ONumber on (cout, 4, false);  // PAR
    cout << "absCriterion [" << half2 << "] = " << absCriterion_half [half2];
    const Prob unexplainedFrac = absCriterion_half [half2] / dissim2_half [half2];
    const Prob r = 1 - unexplainedFrac;
    const Real arcLenError = sqrt ((1 - r) / r); 
    cout << "  Error density [" << half2 << "] = " << arcLenError * 100 << " %";
    cout << "  Fraction = " << (Real) size_half [half2] / (Real) pairs * 100 << " %" << endl;
  }
}



void DistTree::setLeafAbsCriterion () 
{
  ASSERT (optimizable ());
  
  for (DiGraph::Node* node : nodes)
    if (const Leaf* leaf = static_cast <DTNode*> (node) -> asLeaf ())
      const_cast <Leaf*> (leaf) -> absCriterion. clear ();

  FOR (size_t, objNum, ds. objs. size ())
    if (const Real mult = ds. objs [objNum] -> mult)
    {
      const Real d = (*target) [objNum];
    //ASSERT (d >= 0);
      const Real dHat = (*prediction) [objNum];
      ASSERT (! isNan (dHat));
      const Real residual = sqr (dHat - d);      
      ASSERT (DM_sp::finite (residual));
      const_cast <Leaf*> (obj2leaf1 [objNum]) -> absCriterion. add (residual, mult);
      const_cast <Leaf*> (obj2leaf2 [objNum]) -> absCriterion. add (residual, mult);
    }

  for (DiGraph::Node* node : nodes)
    if (const Leaf* leaf = static_cast <DTNode*> (node) -> asLeaf ())
      const_cast <Leaf*> (leaf) -> relLenError = leaf->getRelLenError ();
}



bool DistTree::optimizeLen ()
{
  ASSERT (optimizable ());

  DTNode* toSkip = nullptr;  // toSkip->attr is redundant
  DTNode* toRetain = nullptr; 
  {
    const VectorPtr<DiGraph::Node> rootChildren (root->getChildren ());
    ASSERT (rootChildren. size () >= 2);
    if (rootChildren. size () == 2)  // root is transient in an undirected tree
    {
      toSkip   = const_static_cast <DTNode*> (rootChildren [0]);
      toRetain = const_static_cast <DTNode*> (rootChildren [1]);
      toRetain->len += toSkip->len;
      toSkip->len = 0;
    }
  }
  ASSERT ((bool) toSkip == (bool) toRetain);

  VectorPtr<DTNode> dtNodes;  dtNodes. reserve (nodes. size ());
  for (DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <DTNode*> (node);
    if (   dtNode != root
        && dtNode != toSkip
        && ! dtNode->inDiscernable ()
       )
      dtNodes << dtNode;
    IMPLY (dtNode->inDiscernable (), dtNode->len == 0);
  }

  Space1<NumAttr1> sp (ds, false);  sp. reserve (dtNodes. size ());
  for (const DTNode* dtNode : dtNodes)
    sp << dtNode->attr;
  
  // lr
  if (verbose ())
    cout << "Linear regression ..." << endl;
  L2LinearNumPrediction lr (dsSample, sp, *target);
  ASSERT (lr. beta. size () == dtNodes. size ());
    // lr.beta must be equal to DTNode::len
  FOR (size_t, attrNum, dtNodes. size ())
    lr. beta [attrNum] = dtNodes [attrNum] -> len;
  /*const bool solved =*/ lr. solveUnconstrainedFast (prediction, true, 10, 0.01);  // PAR  // Time = O(leaves^3) 
  if (verbose ())
    lr. qc ();
  if (/*! solved ||*/ isNan (lr. absCriterion))
    return false;

  // DTNode::len
  FOR (size_t, attrNum, dtNodes. size ())
    const_cast <DTNode*> (dtNodes [attrNum]) -> len = lr. beta [attrNum];
  if (toSkip)
  {
    toSkip->len = toRetain->len / 2;
    toRetain->len = toSkip->len;
  }
  
  setPrediction ();
  setAbsCriterion ();  
  
  return true;
}



void DistTree::optimizeLenLocal ()  
{
  if (verbose (1))
    cout << "optimizeLenLocal ..." << endl;

  ASSERT (optimizable ());
  
//cout << "absCriterion (before optimizeLenLocal) = " << absCriterion << endl;  
  Progress prog ((uint) nodes. size ());
  for (const DiGraph::Node* node : nodes)
  {
    prog ();
    const DTNode* dtNode = static_cast <const DTNode*> (node);
    if (   dtNode == root
        || dtNode->asLeaf () 
      //|| dtNode->inDiscernable ()
        || ! dtNode->childrenDiscernable ()
       )
      continue;

    ASSERT (! dtNode->isTransient ());
    
    const Real absCriterion_old = absCriterion;  
    
    VectorPtr<DiGraph::Node> dtNodes (dtNode->getChildren ());
    dtNodes << dtNode;  

    // lenOld, *target_new
    Vector<Real> lenOld;  lenOld. reserve (dtNodes. size ());    
    FOR (size_t, attrNum, dtNodes. size ())
    {
      const DTNode* n = static_cast <const DTNode*> (dtNodes [attrNum]);
      ASSERT (! n->inDiscernable ());
      lenOld << n->len;
      const_cast <DTNode*> (n) -> len = 0;
    }
    ASSERT (lenOld. size () == dtNodes. size ());
    FOR (size_t, objNum, ds. objs. size ())
    {
      const VectorPtr<TreeNode> path (getPath ( obj2leaf1 [objNum]
                                              , obj2leaf2 [objNum]
                                              )
                                     );
      (* const_cast <RealAttr1*> (target_new)) [objNum] = (*target) [objNum] - path2prediction (path);
    }
    
    Space1<NumAttr1> sp (ds, false);  sp. reserve (dtNodes. size ());
    for (const DiGraph::Node* dtNode1 : dtNodes)
      sp << static_cast <const DTNode*> (dtNode1) -> attr;

    // lr
    if (verbose ())
      cout << "Linear regression ..." << endl;
    L2LinearNumPrediction lr (dsSample, sp, *target_new);
    ASSERT (lr. beta. size () == dtNodes. size ());
    bool solved = true;
    lr. solveUnconstrained ();
    FOR (size_t, i, lr. beta. size ())
      if (maximize (lr. beta [i], 0.0))
        solved = false;
    if (verbose ())  
      lr. qc ();
    if (isNan (lr. absCriterion))
    {
    //ASSERT (ds. objs. size () < dissimSize_max ());  // Obj::mult = 0 must be skipped
      solved = false;
    }
  
    // DTNode::len
    FOR (size_t, attrNum, dtNodes. size ())
      const_static_cast <DTNode*> (dtNodes [attrNum]) -> len = solved ? lr. beta [attrNum] : lenOld [attrNum];
    if (verbose ())
    {
      setPrediction ();
      setAbsCriterion ();  
      ASSERT (leReal (absCriterion, absCriterion_old));
    }
  #if 0
    cout << "absCriterion (optimizeLenLocal) = " << absCriterion << endl;  
    if (absCriterion > absCriterion_old)  
    {
      for (const DiGraph::Node* node : dtNodes)
        cout << " " << node->getName ();
      cout << endl;
    }
  #endif
  }

  
  setPrediction ();
  setAbsCriterion ();  
}



void DistTree::optimize2 () 
{
  ASSERT (optimizable ());
  ASSERT (obj2leaf1. size () == 1);

  VectorPtr<Leaf> leaves;
  for (DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <DTNode*> (node);
    if (! dtNode->isLeaf ())
      continue;
    const Leaf* leaf = dtNode->asLeaf ();
    ASSERT (leaf);
    leaves << leaf;
  }
  ASSERT (leaves. size () == 2);  

  const Real t = (*target) [0];
  for (const Leaf* leaf : leaves)
    const_cast <Leaf*> (leaf) -> len = t / 2;
  
  absCriterion = 0;
  * const_cast <RealAttr1*> (prediction) = *target;
}



void DistTree::optimize3 () 
{
  ASSERT (optimizable ());
  ASSERT (obj2leaf1. size () == 3);

  VectorPtr<Leaf> leaves;
  for (DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <DTNode*> (node);
    if (! dtNode->isLeaf ())
      continue;
    const Leaf* leaf = dtNode->asLeaf ();
    ASSERT (leaf);
    leaves << leaf;
  }
  ASSERT (leaves. size () == 3);  

  for (const Leaf* leaf_ : leaves)
  {
    Leaf* leaf = const_cast <Leaf*> (leaf_);
    leaf->len = 0;
    leaf->setParent (const_static_cast <DTNode*> (root));
  }
  FOR (size_t, i, 3)
  {
    const Real t = (*target) [i];
    for (const Leaf* leaf_ : leaves)
    {
      Leaf* leaf = const_cast <Leaf*> (leaf_);
      if (   leaf == obj2leaf1 [i]
          || leaf == obj2leaf2 [i]
         )
        leaf->len += t;
      else
        leaf->len -= t;
    }
  }
  for (const Leaf* leaf_ : leaves)
  {
    Leaf* leaf = const_cast <Leaf*> (leaf_);
    leaf->len /= 2;  // >= 0 <= *target is a distance and triangle inequality 
    maximize (leaf->len, 0.0);
  }

  absCriterion = 0;
  * const_cast <RealAttr1*> (prediction) = *target;
}



bool DistTree::optimize () 
{ 
  ASSERT (optimizable ());
  ASSERT (ds. objs. size () >= 2 * root->getLeavesSize () - 2);
//{(obj2leaf1[obj],obj2leaf2[obj]):obj} must make up a connected graph ??
  
	VectorOwn<Change> changes;  changes. reserve (256);  // PAR
	{
	 	Vector<DTNode*> nodeVec;  nodeVec. reserve (nodes. size ());
    for (DiGraph::Node* node : nodes)
    {
      DTNode* dtNode = static_cast <DTNode*> (node);
      if (! dtNode->stable)
        nodeVec << dtNode;
    }
		Progress prog ((uint) nodeVec. size ());
	 	for (const DTNode* node : nodeVec)  
	 	{
	 	  static Chronometer chr;
	 	  chr. start ();
	   	prog ();
		 	if (const Change* bestChange = getBestChange (node)) 
		  { 
		  	ASSERT (positive (bestChange->improvement));
		  	changes << bestChange;
		  }
		  chr. stop ();
		  if (Chronometer::enabled)
		    cerr << ' ' << chr << endl;
 	  }
  }
    
  return applyChanges (changes);
}



void DistTree::optimizeIter (const string &output_tree)
{
  if (verbose ())
  {
    checkPrediction ();
    checkAbsCriterion ("optimizeIter");
  }
    
  uint iter = 0;
  for (;;)
  {
    iter++;
    if (verbose (1))
      cout << endl << "Topology optimization iter = " << iter << endl;
    if (! optimize ())
      break;	    
    saveFile (output_tree);
    if (verbose ())   
    {
      checkPrediction ();
      checkAbsCriterion ("optimize");
    }
    if (verbose (1))
      cout. flush ();
  }
}



void DistTree::optimizeAdd (bool sparse,
                            const string &output_tree)
{
  ASSERT (! nodes. empty ());
  ASSERT (dissimDs. get ());
  ASSERT (dissimAttr);
  ASSERT (optimizable ());


  {
    if (verbose (1))
      cout << endl << "New data objects ..." << endl;  
    setSubtreeLeaves ();
    Progress prog ((uint) extraObjs ());
    FOR (size_t, row, dissimDs->objs. size ())
    {
      Steiner* root_ = const_cast <Steiner*> (static_cast <const DTNode*> (root) -> asSteiner ());
      ASSERT (root_);
  
      const string name1 = dissimDs->objs [row] -> name;
      if (findPtr (name2leaf, name1))
        continue;
      prog (name1);
      
  		auto g = new Leaf (*this, root_, NAN, name1);  // Temporary parent
  		g->comment = dissimDs->objs [row] -> comment;
      name2leaf [g->name] = g;
  
      const Leaf* closestLocal = nullptr;  // other node
      Real dissim_min = INF;
      Leaf2dist leaf2dist (root_->subtreeLeaves. size (), NAN);
      FOR (size_t, col, dissimDs->objs. size ())
      {
        string name2 = dissimDs->objs [col] -> name;
        const Leaf* other = findPtr (name2leaf, name2);
        if (! other)
          continue;
        if (row == col)
          continue;
        const Real dissim = dissimAttr->get (row, col);
        if (   isNan (dissim)    // prediction must be large ??
            || ! DM_sp::finite (dissim)  
           )
          continue;
        if (minimize (dissim_min, dissim))
          closestLocal = other;
        leaf2dist [other->index] = dissim;
      }
      
      if (! closestLocal)
      {
        if (verbose (1))
          cout << "Deleting singleton: " << g->name << endl;
        EXEC_ASSERT (name2leaf. erase (g->name) == 1);
        delete g;  
        continue;
      }
      ASSERT (DM_sp::finite (dissim_min));
      if (verbose (1))
      {
        ONumber on (cerr, 3, false);  // PAR
        cerr << " dissim_min = " << dissim_min << "  ";
      }
  
    #if 0
      // closestLocal
      Real dist_left = max (0.0, dissim_min / 2);
      g->len = dist_left;
      while (closestLocal != root && greaterReal (dist_left, closestLocal->len))  
      {
        dist_left -= closestLocal->len;
        EXEC_ASSERT (closestLocal = static_cast <const DTNode*> (closestLocal->getParent ()));
      }
      maximize (dist_left, 0.0);
      ASSERT (closestLocal);
      IMPLY (g->len == 0, closestLocal->asLeaf ());
      IMPLY (g->len != 0, ! closestLocal->inDiscernable ());
    #endif
      
      Steiner* st = nullptr;
      if (! positive (dissim_min))
        st = g->collapse (const_cast <Leaf*> (closestLocal));
      else    
      {
        // Best Steiner point on the tree arcs
        ASSERT (g->index == NO_INDEX);
        ASSERT (root_->subtreeLeaves. size () == leaf2dist. size ()); 
        DTNode::Closest closestGlobal;
        Leaf2dist leaf2hat_dist;
        root_->findClosestNode (leaf2dist, leaf2hat_dist, closestGlobal);
        closestGlobal. qc ();
        if (verbose ()) 
        {
          closestGlobal. print (cout);
          cout << "  dist(local,global) = " << path2prediction (getPath (closestLocal, closestGlobal. node)) << "  ";
        }
        // May be: DTNode::len = 0
        EXEC_ASSERT (st = closestGlobal. insert ());
        g->setParent (st);
        g->len = closestGlobal. leafLen;
      }
      addSubtreeLeaf (g);
      
      Set<const Leaf*> representatives;  // size = O(log leaves)
      if (sparse)
      {
        if (! positive (dissim_min))
          representatives << closestLocal;
        sort ();    
        const DTNode* node = g;
        while (const DTNode* ancestor = static_cast <const DTNode*> (node->getParent ()))
        {
          for (const DiGraph::Arc* arc : ancestor->arcs [false])
          {
            const DTNode* child = static_cast <const DTNode*> (arc->node [false]);
            if (child != node)
              if (const Leaf* leaf = child->selectRepresentative (leaf2dist))
                representatives << leaf;
          }
          node = ancestor;
        }
      }
      
      // ds
      const size_t objNum_start = ds. objs. size (); 
      if (sparse)
        for (const Leaf* leaf : representatives)
        {
          ASSERT (leaf);
          ASSERT (leaf != g);
          ASSERT (leaf->index != NO_INDEX);
          EXEC_ASSERT (addDissim (name1, leaf->name, leaf2dist [leaf->index]));  
        }
      else
        FOR (size_t, col, dissimDs->objs. size ())
        {
          string name2 = dissimDs->objs [col] -> name;
          const Leaf* other = findPtr (name2leaf, name2);
          if (! other)
            continue;
          if (row == col)
            continue;
          const Real dissim = dissimAttr->get (row, col);
          if (   isNan (dissim)    // prediction must be large ??
              || ! DM_sp::finite (dissim)  
             )
            continue;
          EXEC_ASSERT (addDissim (name1, name2, dissim));  
        }
      ASSERT (objNum_start < ds. objs. size ());

      dsSample = Sample (ds);
      absCriterion_delta = dsSample. multSum * 1e-5;  // PAR

      g->addAttr ();
      if (st)
        st->addAttr ();
      
      {
        FOR (size_t, objNum, objNum_start)
        {
          const VectorPtr<TreeNode> path (getPath ( obj2leaf1 [objNum]
                                                  , obj2leaf2 [objNum]
                                                  )
                                         );
          ASSERT (! path. contains (g));
          const_cast <CompactBoolAttr1*> (g->attr) -> setCompactBool (objNum, false);
          if (st)
            const_cast <CompactBoolAttr1*> (st->attr) -> setCompactBool (objNum, path. contains (st));
        //(* const_cast <RealAttr1*> (prediction)) [objNum] = path2prediction (path);  
        }
    
        FOR_START (size_t, objNum, objNum_start, ds. objs. size ())
        {
          const VectorPtr<TreeNode> path (getPath ( obj2leaf1 [objNum]
                                                  , obj2leaf2 [objNum]
                                                  )
                                         );
          ASSERT (path. contains (g));
          for (const DiGraph::Node* node : nodes)
          {
            const DTNode* dtNode = static_cast <const DTNode*> (node);
            CompactBoolAttr1* attr = const_cast <CompactBoolAttr1*> (dtNode->attr);
            attr->setCompactBool (objNum, path. contains (dtNode));
          }
          (* const_cast <RealAttr1*> (prediction)) [objNum] = path2prediction (path);  
        }
      }
      
      FOR_START (size_t, objNum, objNum_start, ds. objs. size ())
      {
        const_cast <CompactBoolAttr1*> (fromAttr_new) -> setCompactBool (objNum, false);
        const_cast <CompactBoolAttr1*> (toAttr_new)   -> setCompactBool (objNum, false);
        const_cast <CompactBoolAttr1*> (interAttr)    -> setCompactBool (objNum, false);
      }
  
      setAbsCriterion (); 
      finishChanges (); 
      if (verbose ())
      {
        qc ();  
        qcAttrs ();
      }
  
      setLeafAbsCriterion (); 
      const Real leafRelCriterion = g->getRelLenError ();
      {
        ONumber on (cerr, 1, false);
        cerr << "rel.absCriterion = " << leafRelCriterion << "  ";
      }
      
      if (leafRelCriterion > 2)  // PAR
      {
        Unverbose unv;
        Set<const DTNode*> boundary;  
        optimizeSubtree (const_static_cast <Steiner*> (g->getParent ()), boundary);  // slow ??
        setSubtreeLeaves ();
      }
      if (verbose (1))
        reportErrors (cerr);
  
    #if 0
      const DTNode* boundaryNode = g->discernable ? g : static_cast <const DTNode*> (g->getParent ());
      ASSERT (boundary. contains (boundaryNode));
      const Steiner* optimNode = static_cast <const Steiner*> (boundaryNode->getParent ());
      ASSERT (optimNode);
      ASSERT (! boundary. contains (optimNode));
      Real len_min = INF;
      if (boundary. contains (static_cast <const DTNode*> (optimNode->getParent ())))
        minimize (len_min, optimNode->len);
      const VectorPtr<DiGraph::Node> children (optimNode->getChildren ());
      for (const DiGraph::Node* child : children)
        if (child != boundaryNode && boundary. contains (static_cast <const DTNode*> (child)))
          minimize (len_min, static_cast <const DTNode*> (child)->len);      
      if (verbose (1))
        cerr << "  " << leafLen_min << endl;
    #endif
      if (verbose (1))
        cerr << endl;
        
      if (prog. n % 10 == 0)  // PAR 
        saveFile (output_tree);
    }
    ASSERT (! extraObjs ());
    dissimDs. reset (nullptr);
    dissimAttr = nullptr;
    
  }


  saveFile (output_tree);
    
  if (verbose (1))
    cout << "Optimizing arc lengths..." << endl;
  EXEC_ASSERT (optimizeLen ());   
  finishChanges ();
}



void DistTree::setSubtreeLeaves ()
{
	const size_t leavesSize = root->getLeavesSize ();
	size_t index = 0;
  for (DiGraph::Node* node : nodes)
  {
    DTNode* dtNode = static_cast <DTNode*> (node);
    dtNode->subtreeLeaves. clear ();
    dtNode->subtreeLeaves. resize (leavesSize, false);
    if (Leaf* leaf = const_cast <Leaf*> (dtNode->asLeaf ()))
    {
      leaf->index = index;
      leaf->subtreeLeaves [index] = true;
      index++;
    }
  }
  ASSERT (index == leavesSize);
  
  const_static_cast <DTNode*> (root) -> setSubtreeLeaves ();
}




void DistTree::addSubtreeLeaf (Leaf* leaf)
{
  ASSERT (leaf);
  ASSERT (leaf->graph == this);
  ASSERT (leaf->subtreeLeaves. empty ());
  ASSERT (leaf->index == NO_INDEX);
  
	const size_t leavesSize = static_cast <const DTNode*> (root) -> subtreeLeaves. size ();  
	
	leaf->index = leavesSize;

  for (DiGraph::Node* node : nodes)
  {
    DTNode* dtNode = static_cast <DTNode*> (node);
    IMPLY (dtNode != leaf, dtNode->subtreeLeaves. size () == leavesSize);
    dtNode->subtreeLeaves. resize (leavesSize + 1, false);
  }
    
  // Use in setSubtreeLeaves() ??
  DTNode* node = leaf;
  while (node)
  {
    node->subtreeLeaves [leaf->index] = true;
    node = const_static_cast <DTNode*> (node->getParent ());
  }
}




Real DistTree::optimizeSubtree (const Steiner* center,
                                Set<const DTNode*> &boundary)  
{
  ASSERT (center);
  ASSERT (& center->getTree () == this);
  ASSERT (boundary. empty ());

  VectorPtr<TreeNode> area;
  const DTNode* area_root = nullptr;
  Node2Node new2old;  // Initially: newLeaves2boundary
  DistTree tree (center, areaRadius_std + 1/*PAR*/, area, area_root, new2old);
  if (verbose ())
  {
    tree. qc ();
    tree. qcAttrs ();
  }
  ASSERT (area_root);
  ASSERT (area_root->graph == this);
  ASSERT (area. contains (center));
  ASSERT (area. contains (area_root));
#ifndef NDEBUG
  for (const auto it : new2old)
  {
    // New
    ASSERT (it. first);
    ASSERT (it. first->graph == & tree);
    ASSERT (static_cast <const DTNode*> (it. first) -> asLeaf ());
    ASSERT (! static_cast <const DTNode*> (it. first) -> inDiscernable ());
    // Old = boundary
    ASSERT (it. second);
    ASSERT (it. second->graph == this);
    ASSERT (! static_cast <const DTNode*> (it. second) -> inDiscernable ());
  }
#endif

  // boundary
  for (const auto it : new2old)
    boundary << static_cast <const DTNode*> (it. second);

  const size_t leaves = tree. root->getLeavesSize ();
  {
    Unverbose unv;
    if (leaves > 3)
    {
      EXEC_ASSERT (tree. optimizeLen ());
      tree. finishChanges (); 
      tree. optimizeLenLocal ();  
      tree. finishChanges (); 
      tree. optimizeIter (string ());
    }
    else if (leaves == 3)
      tree. optimize3 ();
    else if (leaves == 2)
      tree. optimize2 ();
    else
    {
      if (verbose (1))
        cout << "Singleton" << endl;
      return INF;
    }
  }
  if (verbose ())
  {
    tree. qc ();
    tree. qcAttrs ();
  }
  if (verbose ())
  {
    cout << "Subtree: ";
    tree. reportErrors (cout);
    cout << endl;
  }
  
  const Real leafLen_min = tree. getMinLeafLen ();  


  Node2Node boundary2new (DiGraph::reverse (new2old));
  ASSERT (boundary2new. size () == new2old. size ());

  const Real absCriterion_old = absCriterion;

#if 0
  cout << "center = " << center << endl;
  cout << "Area:";
  for (const TreeNode* node : area)
    cout << ' ' << node->getName ();
  cout << endl;
  cout << "area_root = " << area_root << endl;
  cout << endl;
#endif

  if (const DTNode* dtNode = static_cast <const DTNode*> (findPtr (boundary2new, static_cast <const DiGraph::Node*> (area_root))))
    if (const Leaf* leaf = dtNode->asLeaf ())
    {
      // tree
      if (verbose ())
        cout << "Re-rooting..." << endl;
      auto st = new Steiner (tree, const_static_cast <Steiner*> (leaf->getParent ()), leaf->len);
      EXEC_ASSERT (new2old. erase (leaf) == 1);
      tree. delayDeleteRetainArcs (const_cast <Leaf*> (leaf));
      ASSERT (isNan (static_cast <const DTNode*> (tree. root) -> len));
      const Steiner* root_ = st->makeDTRoot ();
      boundary2new [area_root] = st;
      new2old [st] = area_root;
      ASSERT (root_ != st);
      ASSERT (! isNan (root_->len));
      if (root_->isTransient ())
        tree. delayDeleteRetainArcs (const_cast <Steiner*> (root_));
      tree. toDelete. deleteData ();
      ASSERT (st == tree. root);
      ASSERT (isNan (st->len));
      // For DTNode::qc()
      tree. ds. deleteAttrs ();
      ASSERT (! tree. optimizable ());
    }
  ASSERT (new2old. size () == boundary2new. size ());
  if (verbose ())
    tree. Tree::qc ();

  // nodes, new2old[]
  // center may be delete'd
  List<DiGraph::Node*> addedNodes;
  for (const TreeNode* node : area)
    if (! contains (boundary2new, static_cast <const DiGraph::Node*> (node)))
    {
      DTNode* dtNode = const_static_cast <DTNode*> (node);
      delete dtNode->attr;
      delete dtNode;
    }
    // Between boundary TreeNode's there are no Arc's
  for (const DiGraph::Node* node_new : tree. nodes)
  {
    DTNode* node = const_static_cast <DTNode*> (findPtr (new2old, node_new));
    if (! node)
    {
      node = new Steiner (*this, nullptr, NAN);
      node->addAttr ();
      new2old [node_new] = node;
      addedNodes << node;
    }
    ASSERT (node);
    if (node_new != tree. root)
      node->len = static_cast <const DTNode*> (node_new) -> len;
  }
  ASSERT (new2old. size () == tree. nodes. size ());

  borrowArcs (new2old, true);  // Arc's are not parallel

  setRoot ();

#if 0  
  cout << "Boundary:";
  for (const auto it : boundary2new)
    cout << ' ' << it. first->getName ();
  cout << endl << endl;
  Tree::qc ();  
#endif
    
  topology2attrs (addedNodes); 
  setPrediction ();
  finishChanges ();
  setAbsCriterion ();
  if (verbose ())
  {
    cout << "absCriterion = " << absCriterion << endl;
    checkPrediction ();
    checkAbsCriterion ("optimizeSubtree");
  }

  if (verbose (1))
    reportErrors (cout);
  	
  ASSERT (absCriterion < absCriterion_old);
  
  return leafLen_min;
}



const Change* DistTree::getBestChange (const DTNode* from) 
{
	ASSERT (from);
	
 	const Change* bestChange = nullptr;
 	
 	#define TRY_CHANGE(T,P)  if (T::valid_ P) tryChange (new T P, bestChange)

  VectorPtr<TreeNode> area;      area.     reserve (nodes. size ()); 
  VectorPtr<TreeNode> boundary;  boundary. reserve (nodes. size ()); 
  from->getArea (areaRadius_std, area, boundary);  
  if (verbose (1))
    cerr << " area=" << area. size () << " ";
    
 	for (const TreeNode* node : area)  
 	{
 		DTNode* to = const_static_cast <DTNode*> (node);
 	  TRY_CHANGE (ChangeToSibling, (from, to));  
	  if (verbose ())
	  {
  	  ASSERT (to->graph);
	    to->qc ();
	  }
  }
  
#if 0
  if (   ! bestChange 
      && area. size () < 100  // PAR
     )
  {
    if (verbose (1))
      cerr << "more ";
   	for (const Tree::TreeNode* node : area)  
   	{
   		DTNode* to = const_static_cast <DTNode*> (node);
 	    // Applied rarely
 	    TRY_CHANGE (ChangeToChild,  (from, to));  
   	  TRY_CHANGE (ChangeToUncle,  (from, to));  
      TRY_CHANGE (ChangeToCousin, (from, to));  
  	  if (verbose ())
  	  {
    	  ASSERT (to->graph);
  	    to->qc ();
  	  }
    }
  }
#endif
  
  #undef TRY_CHANGE
      
  if (bestChange)
  {
    if (verbose (1))
      cerr << "found ";
  	ASSERT (positive (bestChange->improvement));
    return bestChange;	
	}

  return nullptr;
}



bool DistTree::applyChanges (VectorOwn<Change> &changes)
{ 
	ASSERT (toDelete. empty ());
	
	
  const Real absCriterion_init = absCriterion;


  // DTNode::stable: init
  for (DiGraph::Node* node : nodes)
    static_cast <DTNode*> (node) -> stable = true;


  if (verbose (1))
    cout << "# Changes: " << changes. size () << endl;
  size_t commits = 0;
  {
    Unverbose un;
    Common_sp::sort (changes, Change::compare);	
    size_t nChange = 0;
    Progress prog ((uint) changes. size ());
  	for (const Change* ch_ : changes)
  	{
  	  prog ();
  	  nChange++;
  	  
  		Change* ch = const_cast <Change*> (ch_);
  		ASSERT (ch);
  		ASSERT (ch->improvement > 0);

  	  if (! ch->valid ())
    		continue;
  
      if (verbose ())
  		  ch->qc ();  
  
  	#ifndef NDEBUG
  	  const bool first = ch_ == changes. front ();  
  	#endif
  
  	  if (verbose ())
  	    checkAbsCriterion ("Before");
  	    
   	  const Real absCriterion_old = absCriterion;
  
      if (verbose (1))
      {
  	    cout << "Apply " << nChange << "/" << changes. size () << ": ";
     	  ch->print (cout);  
     	}
  	  const bool success = ch->apply ();    
  	  const Real absCriterion_new = getAbsCriterion ();
      if (verbose ())
    	  cout << "absCriterion_new = " << absCriterion_new << endl;
  	  IMPLY (first, success && eqReal (absCriterion_new, absCriterion - ch->improvement));
  	  if (! success || geReal (absCriterion_new, absCriterion))
  	  {
  	  //ASSERT (! first);
        if (verbose (1))
  	      cout << "Restore" << endl;
  	  	ch->restore ();
  	  	ASSERT (eqReal (absCriterion, absCriterion_old));
        if (verbose ())
    	  	checkAbsCriterion ("restore");
  	  }
  	  else
  	  {
       	ch->commit ();
       	commits++;
        setAbsCriterion ();
    	  if (verbose ())
    	  {
     	    cout << "absCriterion = " << absCriterion << endl;
    	    checkAbsCriterion ("commit");
    	    Unverbose unv;
    	    if (verbose ())
    	      qc ();
    	  }
        ASSERT (absCriterion < absCriterion_old);
        // DTNode::stable
        for (const TreeNode* node : ch->targets)
    	    if (node->graph)  
    	    {
            VectorPtr<TreeNode> area, boundary;
            node->getArea (areaRadius_std, area, boundary);  
            for (const TreeNode* areaNode : area)
              const_static_cast <DTNode*> (areaNode) -> stable = false;
          }
      }
  	}
  }
	

  const Real improvement = max (0.0, absCriterion_init - absCriterion);
  if (verbose ())
    cout << "Improvement = " << improvement << endl;
  ASSERT (geReal (improvement, - absCriterion_delta));
  ASSERT ((bool) commits == (bool) improvement);

  if (commits)
  {
    topology2attrs (nodes); 
    EXEC_ASSERT (optimizeLen ()); 
    finishChanges ();
    optimizeLenLocal ();  
    finishChanges ();
    if (verbose (1))
    {
      reportErrors (cout);
      cout << endl;
    }
  }
  
  if (verbose ())
  {
    qc ();
    qcAttrs ();
  }
  
  return positive (improvement);
}



void DistTree::tryChange (Change* ch,
	                        const Change* &bestChange)
{ 
  ASSERT (! isNan (absCriterion));
	ASSERT (ch);
	ASSERT (ch->from->graph == this);
	ASSERT (isNan (ch->improvement));
  
  if (verbose ())
  {
    ch->print (cout); 
    ch->qc ();
  }

  ch->improvement = ch->apply () ? absCriterion - getAbsCriterion () : 0;
  ASSERT (! isNan (ch->improvement));
  ch->restore ();

  Unverbose unv;
  if (verbose ())
  {
    ch->print (cout); 
  	checkAbsCriterion ("tryChange");
  }
  
	if (Change::compare (ch, bestChange))
	{
		delete bestChange;
  	bestChange = ch;
  }
  else
  	delete ch;
}



void DistTree::delayDeleteRetainArcs (DTNode* s)
{
	ASSERT (s);
	ASSERT (! s->inDiscernable ());

  if (verbose ())
    cout << "To delete: " << s->getName () << endl;
    
  if (const Leaf* leaf = s->asLeaf ())
    name2leaf. erase (leaf->name);

  const VectorPtr<DiGraph::Node> children (s->getChildren ());
  for (const DiGraph::Node* child : children)
    const_static_cast <DTNode*> (child) -> len += (s == root ? NAN : s->len);

	delete s->attr;
	s->attr = nullptr;

	s->isolateChildrenUp ();
	ASSERT (! s->graph);

	toDelete << s;
}



void DistTree::finishChanges ()
{
  if (const size_t n = deleteLenZero ())
    if (verbose ())
      cout << "# Nodes with zero arcs deleted = " << n << endl;
  toDelete. deleteData ();
}



size_t DistTree::deleteLenZero ()
{
  size_t n = 0;
 	Vector<DiGraph::Node*> nodeVec;  nodeVec. reserve (nodes. size ());
 	insertAll (nodeVec, nodes);
 	for (DiGraph::Node* node :  nodeVec)  
 		if (const Steiner* s = static_cast <DTNode*> (node) -> asSteiner ())
 			if (   s->getParent () 
 				  && s->len == 0
 				  && s->childrenDiscernable ()
 				 )
 			{
 				delayDeleteRetainArcs (const_cast <Steiner*> (s));
 				n++;
 			}
 			
  return n;  
}



void DistTree::setReprLeaves ()
{
  for (const DiGraph::Node* node : nodes)
    const_static_cast <DTNode*> (node) -> reprLeaf = nullptr;
  sort ();
  for (const DiGraph::Node* node : nodes)
  {
    DTNode* dtNode = const_static_cast <DTNode*> (node);
    EXEC_ASSERT (dtNode->reprLeaf = static_cast <const DTNode*> (dtNode->getLeftmostDescendent ()) -> asLeaf ());
  }
}



void DistTree::reroot (DTNode* underRoot,
                       Real arcLen) 
{
  ASSERT (underRoot);
  ASSERT (& underRoot->getTree () == this);
  ASSERT (underRoot != root);  
  ASSERT (! underRoot->inDiscernable ());

  ASSERT (! isNan (arcLen));
  ASSERT (arcLen >= 0);
  ASSERT (arcLen <= underRoot->len);
  
  DTNode* root_ = const_static_cast<DTNode*> (root);
    
  const Steiner* newRootParent = static_cast <const DTNode*> (underRoot->getParent ()) -> asSteiner ();
  auto newRoot = new Steiner (*this, const_cast <Steiner*> (newRootParent), underRoot->len - arcLen);
  if (optimizable ())
  {
    newRoot->addAttr ();
    * const_cast <CompactBoolAttr1*> (newRoot->attr) = * underRoot->attr;
  }
  underRoot->setParent (newRoot); 
  underRoot->len = arcLen;
  
  newRoot->makeDTRoot ();
  ASSERT (newRoot == root);
  ASSERT (root_ != root);
  
  if (root_->isTransient ())
    delayDeleteRetainArcs (root_);

  finishChanges ();
}



void DistTree::reroot ()
{
  DTNode* root_ = const_static_cast<DTNode*> (root);
  
  if (! root_->childrenDiscernable ())
    return;
  
  root_->setSubtreeLenUp ();
  
  DTNode* bestDTNode = nullptr;
  Real bestDTNodeLen_new = NAN;
  WeightedMeanVar bestGlobalLen;
  root_->setGlobalLenDown (bestDTNode, bestDTNodeLen_new, bestGlobalLen);
  ASSERT (bestDTNode);
  
  if (verbose ())
  {
    cout << bestDTNode << " " << bestDTNodeLen_new << endl;
    cout << bestGlobalLen. getMean () << endl;
    print (cout);
  }
  
  reroot (bestDTNode, bestDTNodeLen_new);
}



Real DistTree::getMeanResidual () const
{
  ASSERT (optimizable ());

  Real s = 0;
  FOR (size_t, objNum, ds. objs. size ())
    if (const Real mult = ds. objs [objNum] -> mult)
    {
      const Real d = (*target) [objNum];
    //ASSERT (d >= 0);
      const Real dHat = (*prediction) [objNum];
      ASSERT (! isNan (dHat));
      s += mult * (dHat - d);
    }
  
  return s / dsSample. multSum;
}



Real DistTree::getMinLeafLen () const
{
  Real len_min = INF;
  for (DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <DTNode*> (node);
    if (   (  dtNode->isLeaf () && ! dtNode->inDiscernable ())
        || (! dtNode->isLeaf () && ! dtNode->childrenDiscernable ())
       )
      minimize (len_min, dtNode->len);
  }
  return len_min;
}



Real DistTree::getSqrResidualCorr () const
{
  ASSERT (optimizable ());

  Correlation corr;
  FOR (size_t, objNum, ds. objs. size ())
  {
    const Real d = (*target) [objNum];
  //ASSERT (d >= 0);
    const Real dHat = (*prediction) [objNum];
    ASSERT (! isNan (dHat));
    const Real residual2 = sqr (dHat - d);
    corr. add (d, residual2);
  }
  
  return corr. getCorrelation ();
}
  
  

Real DistTree::setErrorDensities () 
{
  ASSERT (optimizable ());

  for (const DiGraph::Node* node : nodes)
  {
    DTNode* dtNode = const_static_cast <DTNode*> (node);
    dtNode->paths = 0;
    dtNode->errorDensity = 0;
  }
  
  Real epsilon2_0 = 0;
  FOR (size_t, objNum, obj2leaf1. size ())
  {
    const Real mult = ds. objs [objNum] -> mult;
    const Real d = (*target) [objNum];
  //ASSERT (d >= 0);
    const Real dHat = (*prediction) [objNum];
    ASSERT (! isNan (dHat));
    if (nullReal (dHat))
    {
      epsilon2_0 += mult * sqr (d);
      continue;
    }

    ASSERT (positive (d));
    const Real a = mult * sqr (dHat - d) / dHat;
    if (! DM_sp::finite (a))
    {
      cout << obj2leaf1 [objNum] -> getName () << " " << obj2leaf2 [objNum] -> getName () << ": " << d << " " << dHat << endl;
      ERROR;
    }
    
    const VectorPtr<TreeNode> path (getPath ( obj2leaf1 [objNum]
                                            , obj2leaf2 [objNum]
                                            )
                                   );
    for (const TreeNode* node : path)
    {
      DTNode* dtNode = const_static_cast <DTNode*> (node);
      dtNode->paths++;
      dtNode->errorDensity += a;
    }
  }

  for (const DiGraph::Node* node : nodes)
  {
    DTNode* dtNode = const_static_cast <DTNode*> (node);
    if (! dtNode->getParent ())
      continue;
    ASSERT (dtNode->paths);
    dtNode->errorDensity = sqrt (dtNode->errorDensity / (Real) dtNode->paths);
  }
  
  return epsilon2_0;
}



size_t DistTree::printLeafRelLenErros (ostream &os,
                                       Real relErr_min) const
{
  ASSERT (os. good ());
  ASSERT (relErr_min >= 0);
  
  ASSERT (optimizable ());

  size_t n = 0;
  for (const DiGraph::Node* node : nodes)
    if (const Leaf* leaf = static_cast <const DTNode*> (node) -> asLeaf ())
    {
      const Real relErr = leaf->getRelLenError ();
      if (lessReal (relErr, relErr_min)) 
        continue;
      os << leaf->getName () << '\t' << relErr << endl;
      n++;
    }
    
  return n;
}



RealAttr1* DistTree::getResiduals2 () 
{
  ASSERT (optimizable ());

  Common_sp::AutoPtr<RealAttr1> resid2Attr (new RealAttr1 ("resid2", ds, target->decimals + 1));
  FOR (size_t, objNum, ds. objs. size ())
    if (ds. objs [objNum] -> mult)
    {
      const Real d = (*target) [objNum];
    //ASSERT (d >= 0);
      const Real dHat = (*prediction) [objNum];
      ASSERT (! isNan (dHat));
      const Real residual = abs (dHat - d);      
      ASSERT (DM_sp::finite (residual));
      (*resid2Attr) [objNum] = sqr (residual);
    }
    
  return resid2Attr. release ();
}



RealAttr1* DistTree::getLogPredictionDiff () 
{
  ASSERT (optimizable ());

  Common_sp::AutoPtr<RealAttr1> logDiffAttr (new RealAttr1 ("logDiff", ds, target->decimals + 1));
  FOR (size_t, objNum, ds. objs. size ())
    if (ds. objs [objNum] -> mult)
    {
      const Real d = (*target) [objNum];
    //ASSERT (d >= 0);
      const Real dHat = (*prediction) [objNum];
      ASSERT (! isNan (dHat));
      if (   positive (d)
          && positive (dHat)
         )
        (*logDiffAttr) [objNum] = log (d) - log (dHat);
    }
    
  return logDiffAttr. release ();
}



void DistTree::pairResiduals2dm (const RealAttr1* resid2Attr,
                                 const RealAttr1* logDiffAttr,
                                 ostream &os) const
{
  ASSERT (optimizable ());
  ASSERT (resid2Attr);
  ASSERT (& resid2Attr->ds == & ds);
  ASSERT (logDiffAttr);
  ASSERT (& logDiffAttr->ds == & ds);
    
  VectorPtr<Attr> attrs;
  attrs << target << prediction << resid2Attr << logDiffAttr;
  dsSample. save (attrs, os);
}




}





/* TO DO ??

not solved LinearRegression
  
Subtree for optimizeSubtree():
  many short arcs
  big error in the subgraph
  
optimize():
  cover by optimizeSubtree() if nodes.size() is large: should improve speed without decreasing the criterion

optimizeAdd():
  sparse() each new N genomes  
    Dataset::sparse()

non-stability of results

*/
