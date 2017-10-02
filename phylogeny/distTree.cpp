// distTree.cpp

#undef NDEBUG
#include "../common.inc"

#include "distTree.hpp"

#include "../dm/prediction.hpp"



namespace DistTree_sp
{


// VarianceType

const StringVector varianceTypeNames {"lin", "exp", "linExp"};
VarianceType varianceType = varianceType_linExp;




// DTNode

DTNode::DTNode (DistTree &tree,
                Steiner* parent_arg,
    	          Real len_arg)
: TreeNode (tree, parent_arg)  // DistTree must be declared
, len (len_arg)
{}



void DTNode::qc () const
{ 
  if (! qc_on)
    return;
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
  {
    const ONumber oLen (os, dissimDecimals, true);
    os << "len=" << len;
		if (subtreeLen. weights)
		  os << "  len_mean=" << getHeight ();
  }
  const ONumber oNum (os, 6, true);  // PAR
	if (subtreeLen. weights)
	  os << "  len_weight=" << subtreeLen. weights;
  if (paths)
	  os << "  err_density=" << errorDensity << "  paths=" << paths;
//os << "  " << dissimSum << ' ' << dissimWeightedSum;  
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
  const ONumber on (os, dissimDecimals, true);
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



#if 0
void DTNode::setSubtreeLeaves ()
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



void DTNode::getCenters (size_t depth,
                         VectorPtr<Steiner> &centers) const
{
  static_assert (areaRadius_std > 1, "areaRadius_std > 1");
  ASSERT (depth < areaRadius_std);
  
  if (depth == 0)
  {
    const Steiner* st = asSteiner ();
    if (! st)
      st = static_cast <const Steiner*> (getParent ());
    centers << st;
  }

  depth++;
  if (depth == areaRadius_std)
    depth = 0;
    
  for (const DiGraph::Arc* arc : arcs [false])
  {
    const DTNode* child = static_cast <const DTNode*> (arc->node [false]);
    child->getCenters (depth, centers);
  }
}
#endif




// Steiner


Steiner::Steiner (DistTree &tree,
      	          Steiner* parent_arg,
      	          Real len_arg)
: DTNode (tree, parent_arg, len_arg)  // DistTree must be declared
{
#if 0
  for (const bool b : Bool)
    bootstrap [b] = 0;
#endif
}



void Steiner::qc () const
{
  if (! qc_on)
    return;
	DTNode::qc ();
	  
	ASSERT (! isLeaf ());	
//IMPLY (reprLeaf, getChildren (). contains (reprLeaf));
}



void Steiner::getDescendants (VectorPtr<DTNode> &descendants,
                              size_t depth) const
{
  if (depth)
    if (childrenDiscernable ())
      for (const DiGraph::Arc* arc : arcs [false])
      {
        const DTNode* child = static_cast <const DTNode*> (arc->node [false]);
        child->getDescendants (descendants, depth - 1);
      }
    else
      descendants << this;  
  else
    descendants << this;  
}



void Steiner::reverseParent (const Steiner* target, 
                             Steiner* child)
{
  ASSERT (target);
  ASSERT (descendantOf (target));
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
  ASSERT (descendantOf (ancestor2descendant));

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

const string Leaf::non_discernable ("non-discernable");



Leaf::Leaf (DistTree &tree,
  	        Steiner* parent_arg,
  	        Real len_arg,
  	        const string &name_arg)
: DTNode (tree, parent_arg, len_arg)  // DistTree must be declared
, name (name_arg)
{}



void Leaf::qc () const
{ 
  if (! qc_on)
    return;
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



#if 0
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
#endif




// Change

Change::~Change ()
{
  ASSERT (status != eApplied);
}



void Change::qc () const
{
  if (! qc_on)
    return;
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



bool Change::strictlyLess (const Change* a, 
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
        
  if (! dissimFName. empty ())
  {
    loadDissimDs (dissimFName, attrName);
    
    if (! getConnected ())
      throw runtime_error ("Disconnected objects");
  	setDiscernable ();  

    dissimDs2ds (sparse);    	
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
  // Initial tree topology, no DTNode::len
  loadTreeDir (dirName);
  ASSERT (root);
  ASSERT (nodes. front () == root);
  ASSERT (static_cast <const DTNode*> (root) -> asSteiner ());

  setName2leaf ();
        
  loadDissimDs (dissimFName, attrName);
  
  if (! getConnected ())
    throw runtime_error ("Disconnected objects");
	setDiscernable ();  

  setGlobalLen ();  

  dissimDs2ds (false);  
  topology2attrs (nodes);
  setPrediction ();
  setAbsCriterion (); 
}



DistTree::DistTree (const string &dissimFName,
	                  const string &attrName,
	                  bool sparse)
: dsSample (ds)
{
  loadDissimDs (dissimFName, attrName);

  // Initial tree topology: star topology of indiscernable objects
  ASSERT (! root);
  auto root_ = new Steiner (*this, nullptr, NAN);
  for (const Obj* obj : dissimDs->objs)
	  new Leaf (*this, root_, 0, obj->name);
  ASSERT (root);
  ASSERT (nodes. front () == root);
  ASSERT (static_cast <const DTNode*> (root) -> asSteiner ());

  setName2leaf ();
        
  if (! getConnected ())
    throw runtime_error ("Disconnected objects");
	setDiscernable ();

  neighborJoin ();
  
  dissimDs2ds (sparse);    	 
  topology2attrs (nodes);
  setPrediction ();
  setAbsCriterion (); 
}



DistTree::DistTree (const string &dataDirName,
 	                  bool loadDissim)
: completeDs (false)
, dsSample (ds)
{
	ASSERT (! dataDirName. empty ());
	ASSERT (dataDirName. back () == '/');

  // Initial tree topology
  loadTreeFile (dataDirName + "tree");  
  ASSERT (root);
  ASSERT (nodes. front () == root);
  ASSERT (static_cast <const DTNode*> (root) -> asSteiner ());
  
  setName2leaf ();

  VectorPtr<Leaf> newLeaves;  newLeaves. reserve (nodes. size ());
  {
    LineInput f (dataDirName + "leaf", 10 * 1024, 1);  // PAR
    string leafName, anchorName;
    Real leafLen, arcLen;
    while (f. nextLine ())
    {
      istringstream iss (f. line);
      iss >> leafName >> anchorName >> leafLen >> arcLen;
      ASSERT (iss. eof ());
      DTNode* anchor = const_cast <DTNode*> (lcaName2node (anchorName));
      ASSERT (anchor);
      Leaf* leaf = nullptr;
      if (! anchor->childrenDiscernable () && leafLen == 0 && arcLen == 0)
      {
        const Steiner* anchorSt = anchor->asSteiner ();
        ASSERT (anchorSt);
        leaf = new Leaf (*this, const_cast <Steiner*> (anchorSt), 0, leafName);
        leaf->discernable = false;
      }
      else
      {
        auto st = new Steiner ( *this
                              , const_static_cast <Steiner*> (anchor->getParent ())
                              , anchor->len - arcLen
                              );
        anchor->setParent (st);
        anchor->len = arcLen;
        leaf = new Leaf (*this, st, leafLen, leafName);
        if (const Leaf* anchorLeaf = anchor->asLeaf ())
          if (leafLen == 0 && arcLen == 0)
          {
            const_cast <Leaf*> (anchorLeaf) -> discernable = false;
            leaf->discernable = false;
          }
      }
      ASSERT (leaf);
      name2leaf [leaf->name] = leaf;
      newLeaves << leaf;
    }
  }


  if (loadDissim)
  {
    loadDissimPrepare (name2leaf. size () * (size_t) log ((Real) name2leaf. size () * 10), dissimDecimals);  // PAR
    {
      const string fName (dataDirName + "dissim");
      cout << "Loading " << fName << " ..." << endl;
      LineInput f (fName, 10 * 1024 * 1024, 100000);  // PAR
      while (f. nextLine ())
      {
        const string name1 = findSplit (f. line, '\t');
        const string name2 = findSplit (f. line, '\t');
        if (! name2leaf [name1])
          throw runtime_error ("Tree has no object " + name1);
        if (! name2leaf [name2])
          throw runtime_error ("Tree has no object " + name2);
        const Real dissim = str2<Real> (f. line);
      #if 0
        ASSERT (! f. line. empty ());
        if (   contains (f. line, 'e')
            || contains (f. line, 'E')
           )
          throw runtime_error ("Dissimilarities must be in a fixed point format");
        const size_t dot = f. line. find ('.');
        if (dot != string::npos)
          maximize (dissimDecimals, (streamsize) (f. line. size () - 1 - dot));
      #endif
        EXEC_ASSERT (addDissim (name1, name2, dissim));
      }
    }
  
    loadDissimFinish ();      
  
    topology2attrs (nodes);
    setPrediction ();
    setAbsCriterion (); 

    if (verbose ())
      qc ();     

    {  
      cout << "Optimizing new leaves ..." << endl;
      Progress prog ((uint) newLeaves. size ());
      for (const Leaf* leaf : newLeaves)
      {
        prog (real2str (absCriterion, 6));
        Unverbose unv;
        optimizeSubgraph (static_cast <const Steiner*> (leaf->getParent ()));
      }
    }
  }  


  ASSERT (! dissimDs. get ());
  ASSERT (! dissimAttr);
}



DistTree::DistTree (const string &newickFName)
: completeDs (false)
, dsSample (ds)
{ 
  {
    ifstream f (newickFName);
    newick2node (f, nullptr);
  }

	ASSERT (root);
  ASSERT (nodes. front () == root);
  if (! static_cast <const DTNode*> (root) -> asSteiner ())
    throw runtime_error ("One-node tree");
    
  finishChanges ();

  setName2leaf ();        
}



DistTree::DistTree (Prob branchProb,
                    size_t leafNum_max)
: completeDs (false)
, dsSample (ds)
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
  ASSERT (areaRadius >= 1);  
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
    DTNode* node_new = nullptr;
    {
      const Real len = static_cast <const DTNode*> (node_old) -> len;
      if (intersects (node_old->getChildren (), areaSet))
        node_new = new Steiner (*this, nullptr, len);
      else  
      {
        leafNum++;
        node_new = new Leaf (*this, nullptr, len, "L" + toString (leafNum));
      }
    }
    ASSERT (node_new);
    old2new [node_old] = node_new;
  }
  ASSERT (area_root);
  ASSERT (areaSet. contains (area_root));
  ASSERT (! areaSet. contains (area_root->getParent ()));
  ASSERT (nodes. size () == area. size ());
  ASSERT (nodes. size () == old2new. size ());
  
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
      const TreeNode* root_ = root;                    // Old root in *this
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
  for (const auto it : old2new)
    if (static_cast <const DTNode*> (it. second) -> asLeaf ())
      newLeaves2boundary [it. second] = it. first;
  ASSERT (newLeaves2boundary. size () == boundary. size ());
  

  // ds.objs[]->mult, *target, obj2leaf1[], obj2leaf2[]: init
  // For some leaf pairs the dissimilarity may be missing
  ds. objs. reserve (dissimSize_max ());
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
      obj2leaf1 << it1. second;
      obj2leaf2 << it2. second;
    }
  ds. setName2objNum ();


  // ds.objs[]->mult, *target: sum
  const VectorOwn<Obj>& wholeObjs = wholeTree. ds. objs;
  VectorPtr<TreeNode> extremes (2);  // temporary
  FOR (size_t, objNum_whole, wholeObjs. size ())
  {
    if (wholeObjs [objNum_whole] -> mult == 0)
      continue;
      
    VectorPtr<TreeNode> path (wholeTree. getPath ( wholeTree. obj2leaf1 [objNum_whole]
                                                 , wholeTree. obj2leaf2 [objNum_whole]
                                                 )
                             );
    const Real dist_hat_whole = path2prediction (path);

    path. filter ([&] (size_t i) { return ! areaSet. contains (path [i]) || path [i] == area_root; });    
    if (path. empty ())
      continue;

    const Real dist_hat_sub = path2prediction (path);

    // extremes
    extremes. clear ();
    for (const TreeNode* node : path)
      if (boundarySet. contains (node))
        extremes << node;
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


  // *target: finish
  // dissim2_sum
  FOR_REV (size_t, objNum, ds. objs. size ())
  {
    const Real mult = ds. objs [objNum] -> mult;
    if (mult == 0)
      const_cast <RealAttr1*> (target) -> setMissing (objNum);
    else
    {
      (* const_cast <RealAttr1*> (target)) [objNum] /= mult;
      dissim2_sum += mult * sqr ((*target) [objNum]);  // max (0.0, dist);
    }
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

  StringVector fileNames;  
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
    Unverbose unv;
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

  StringVector lines;
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



bool DistTree::loadLines (const StringVector &lines,
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
  const bool indiscernable = contains (s, Leaf::non_discernable);
  IMPLY (parent, ! isNan (len));
  DTNode* dtNode = nullptr;
	if (isLeft (idS, "0x"))
	{
	  ASSERT (isNan (leafError));
	  ASSERT (! indiscernable);
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
		leaf->discernable = ! indiscernable;
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
    if (const Leaf* leaf = static_cast <const DTNode*> (node) -> asLeaf ())
      name2leaf [leaf->name] = leaf;
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
    Unverbose unv;
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
  
//dissimDecimals = dissimAttr->decimals;
}



void DistTree::newick2node (ifstream &f,
                            Steiner* parent)
{
  char c;
  f >> c;
  DTNode* node = nullptr;
  if (c == '(')
  {
    auto st = new Steiner (*this, parent, NAN);
    for (;;)
    {
      newick2node (f, st);
      f >> c;
      if (c == ')')
        break;
      if (c != ',')
        throw runtime_error ("Comma expected");
    }
    f >> c;
    node = st;
  }
  else
  {
    string name;
    while (c != ':' && c != ';')
    {
      name += c;
      f >> c;
    }
    if (name. empty ())
      throw runtime_error ("Empty name");
    node = new Leaf (*this, parent, NAN, name);
  }
  ASSERT (node);
  
  if (c == ';')
    return;
  
  if (c != ':')
    throw runtime_error ("':' expected");

  Real len;
  f >> len;
  node->len = max (0.0, len);
}



bool DistTree::getConnected () 
{
	ASSERT (dissimDs. get ());
	ASSERT (dissimAttr);
	ASSERT (! optimizable ());


  map <const DisjointCluster*, VectorPtr<Leaf>> cluster2leaves;
  {
    // Leaf::DisjointCluster
   	for (DiGraph::Node* node : nodes)
   	{
   	  const DTNode* dtNode = static_cast <DTNode*> (node);
   	  if (Leaf* leaf = const_cast <Leaf*> (dtNode->asLeaf ()))
   	    leaf->DisjointCluster::init ();
   	}
    FOR (size_t, row, dissimDs->objs. size ())
    {
      const Leaf* leaf1 = name2leaf [dissimDs->objs [row] -> name];
      ASSERT (leaf1);
      FOR (size_t, col, row)  // dissimAttr is symmetric
        if (dissimAttr->get (row, col) < INF)
        {
          const Leaf* leaf2 = name2leaf [dissimDs->objs [col] -> name];
          ASSERT (leaf2);
          const_cast <Leaf*> (leaf1) -> merge (* const_cast <Leaf*> (leaf2));
        }
    }
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



size_t DistTree::setDiscernable ()
{
	ASSERT (! optimizable ());
  ASSERT (dissimDs. get ());
  ASSERT (dissimAttr);


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
    FOR (size_t, row, dissimDs->objs. size ())
    {
      const Leaf* leaf1 = name2leaf [dissimDs->objs [row] -> name];
      ASSERT (leaf1);
      FOR (size_t, col, row)  // dissimAttr is symmetric
        if (! positive (dissimAttr->get (row, col)))
        {
          const Leaf* leaf2 = name2leaf [dissimDs->objs [col] -> name];
          ASSERT (leaf2);
          const_cast <Leaf*> (leaf1) -> merge (* const_cast <Leaf*> (leaf2));
        }
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



void DistTree::setGlobalLen ()
{
	ASSERT (dissimDs. get ());
	ASSERT (dissimAttr);
	ASSERT (! optimizable ());


  // DTNode::subtreeLen
  for (DiGraph::Node* node : nodes)
  {
    DTNode* dtNode = static_cast <DTNode*> (node);
    dtNode->subtreeLen. clear ();
    if (Leaf* leaf = const_cast <Leaf*> (dtNode->asLeaf ()))
      leaf->subtreeLen. add (0);  
  }
  FOR (size_t, row, dissimDs->objs. size ())
  {
    const Leaf* leaf1 = name2leaf [dissimDs->objs [row] -> name];
    ASSERT (leaf1);
    FOR (size_t, col, row)  // dissimAttr is symmetric
    {
      const Leaf* leaf2 = name2leaf [dissimDs->objs [col] -> name];
      ASSERT (leaf2);
      const Real d = dissimAttr->get (row, col);
      if (isNan (d))
        continue;
      const TreeNode* ancestor = getLowestCommonAncestor (leaf1, leaf2);
      ASSERT (ancestor);
      Steiner* s = const_cast <Steiner*> (static_cast <const DTNode*> (ancestor) -> asSteiner ());
      ASSERT (s);
      s->subtreeLen. add (max (0.0, d) / 2);  
    }
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

  finishChanges ();
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

    Neighbors ()
      { nodes [0] = nullptr;
        nodes [1] = nullptr;
      }
    Neighbors (const Leaf* leaf1,
               const Leaf* leaf2,
               Real dissim_arg)
      : dissim (dissim_arg)
      { ASSERT (leaf1);
        ASSERT (leaf2);
        nodes [0] = leaf1->getDiscernable ();
        nodes [1] = leaf2->getDiscernable ();
        orderNodes ();
      }
      
    void print (ostream &os) const
      { os        << nodes [0] -> getLcaName () 
           << ' ' << nodes [1] -> getLcaName () 
           << ' ' << dissim 
           << endl; 
      }
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
        dissim = (dissim + from. dissim) / 2;
        return true;
      }
    Real getCriterion (size_t n) const
      { return dissim - (nodes [0] -> len + nodes [1] -> len) / (Real) (n - 2); }
      // James A. Studier, Karl J. Keppler, A Note on the Neighbor-Joining Algorithm of Saitou and Nei
    Real getParentDissim (size_t n) const
      { return min (dissim, max (0.0, 0.5 * (dissim + (nodes [0] -> len - nodes [1] -> len) / (Real) (n - 2)))); }
    static bool strictlyLess (const Neighbors& n1,
                              const Neighbors& n2)
      { LESS_PART (n1, n2, nodes [0]);
        LESS_PART (n1, n2, nodes [1]);
        return false;  
      }
  };
}



void DistTree::neighborJoin ()
{
	ASSERT (dissimDs. get ());
	ASSERT (dissimAttr);
	ASSERT (! optimizable ());
  ASSERT (dissimDs->objs. size () >= 2);
  ASSERT (completeDs);
    
  cout << "Neighbor joining ..." << endl;
  
  // DTNode::len: sum of dissimilarities from other objects (dissim_sum)

  size_t n = 0;
	for (const DiGraph::Arc* arc : root->arcs [false])
	{
	  const_static_cast <DTNode*> (arc->node [false]) -> len = 0;
	  n++;
	}
    
  Vector<Neighbors> neighborsVec;  neighborsVec. reserve (dissimSize_max ());
  FOR (size_t, row, dissimDs->objs. size ())
  {
    const Leaf* leaf1 = name2leaf [dissimDs->objs [row] -> name];
    FOR (size_t, col, row)  // dissimAttr is symmetric
    {
      const Leaf* leaf2 = name2leaf [dissimDs->objs [col] -> name];
      const Neighbors neighbors (leaf1, leaf2, dissimAttr->get (row, col));
      if (neighbors. same ())
        continue;
      if (! positive (neighbors. dissim))
        throw runtime_error ("No distance for " + leaf1->name + " - " + leaf2->name);
      if (neighbors. dissim == INF)
        throw runtime_error ("Infinite distance for " + leaf1->name + " - " + leaf2->name);
      neighborsVec << neighbors;
    }
  }  
  if (neighborsVec. empty ())
    return;
  

  Progress prog ((uint) n - 1);
  Neighbors neighbors_best;
  for (;;)
  {
    prog ();
    
    // Remove duplicate Neighbors 
    Common_sp::sort (neighborsVec, Neighbors::strictlyLess);
    neighborsVec. filter ([&] (size_t i) 
                            { return    neighborsVec [i] == neighbors_best 
                                     || (i && neighborsVec [i - 1]. merge (neighborsVec [i]));
                            }
                         );    

    if (neighborsVec. size () == 1)
      break;
    ASSERT (n > 2);
    
    if (! neighbors_best. nodes [0])  // First iteration
      for (Neighbors& neighbors : neighborsVec)
        for (const bool first : Bool)
          const_cast <DTNode*> (neighbors. nodes [first]) -> len += neighbors. dissim;
          
    if (verbose (-2))  
    {
      cout << endl << "Nodes and sum:" << endl;
    	for (const DiGraph::Arc* arc : root->arcs [false])
    	{
    	  const DTNode* node = static_cast <const DTNode*> (arc->node [false]);
    	  cout << node->getLcaName () << ": " << node->len << endl;
    	}    	
    	cout << endl << "Pairs:" << endl;
      for (const Neighbors& neighbors : neighborsVec)
        neighbors. print (cout);
    }
      
      
    Steiner* newNode = nullptr;
    {
      // neighbors_best
      Real dissim_min = NAN;
      {
        Real criterion = INF;
        size_t i_best = NO_INDEX;
        FOR (size_t, i, neighborsVec. size ())
          // P (criterion1 < criterion2) ??
          if (minimize (criterion, neighborsVec [i]. getCriterion (n)))            
            i_best = i;
        ASSERT (i_best != NO_INDEX);
        neighbors_best = neighborsVec [i_best];
        if (verbose (-1))
        {
          cout << endl << "Best: ";
          neighbors_best. print (cout);
        }
        dissim_min = neighbors_best. getParentDissim (n);
      }
      ASSERT (dissim_min >= 0);
      ASSERT (dissim_min <= neighbors_best. dissim);
      {
        DTNode* a = const_cast <DTNode*> (neighbors_best. nodes [0]);
        DTNode* b = const_cast <DTNode*> (neighbors_best. nodes [1]);
        const Real dissim_sum_a = a->len;
        const Real dissim_sum_b = b->len;
        a->len = dissim_min;
        b->len = neighbors_best. dissim - dissim_min;
        // dissim_sum
        const Real dissim_sum = (  dissim_sum_a - neighbors_best. dissim - (Real) (n - 2) * a->len
                                 + dissim_sum_b - neighbors_best. dissim - (Real) (n - 2) * b->len
                                ) / 2;
        newNode = new Steiner (*this, const_static_cast <Steiner*> (neighbors_best. nodes [0] -> getParent ()), dissim_sum);
        a->setParent (newNode);
        b->setParent (newNode);
        if (verbose (-1))
          cout << "New: " << newNode->getLcaName () << " " << a->len << " " << b->len << " " << newNode->len << endl;
      }
    }
    ASSERT (newNode);
    
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
              const_cast <DTNode*> (neighbors. nodes [! first]) -> len -= neighbors. dissim;
              neighbors. dissim -= neighbors_best. nodes [best_first] -> len;
              maximize (neighbors. dissim, 0.0);
              const_cast <DTNode*> (neighbors. nodes [! first]) -> len += neighbors. dissim / 2;  
                // Done twice for neighbors. nodes [! first]
              found = true;
            }
          if (found)
            neighbors. orderNodes ();
          ASSERT (! neighbors. same ());
        }

    n--;
  }
  ASSERT (neighborsVec. size () == 1);
  ASSERT (n == 2);
  
  
  Neighbors& neighbors = neighborsVec [0];
  ASSERT (! neighbors. same ());
  for (const bool first : Bool)  
  {
    DTNode* node = const_cast <DTNode*> (neighbors. nodes [first]);
    ASSERT (node->getParent () == root);
    node->len = neighbors. dissim / 2;  
  }

  finishChanges ();
  
  reroot ();  // Reduces ds.objs.size() if ds is sparse
}



void DistTree::dissimDs2ds (bool sparse)
{
	ASSERT (dissimDs. get ());
	ASSERT (dissimAttr);
	ASSERT (! optimizable ());
	ASSERT (! ds. name2objNumSet ());


  const Set<string> objNames (dissimDs->getObjNames ());
  restrictLeaves (objNames, true);
#ifndef NDEBUG
  {
    const Set<string> leafNames (name2leaf);
    ASSERT (objNames. contains (leafNames));
  }
#endif
  ASSERT (name2leaf. size () <= dissimDs->objs. size ());

  // Leaf::comment
  for (const Obj* obj : dissimDs->objs)
  {
    const Leaf* leaf = findPtr (name2leaf, obj->name);
    IMPLY (! extraObjs (), leaf);
    if (leaf)
      const_cast <Leaf*> (leaf) -> comment = obj->comment;
  }

  // Sparsing
  Set<string> selectedPairs;
  if (sparse)
  {
    ASSERT (ds. objs. empty ());  // => ds does not affect selectedPairs
    selectedPairs = selectPairs ();
  }

  // ds.objs, *target, obj2leaf1, obj2leaf2, dissim2_sum
  loadDissimPrepare (dissimSize_max (), dissimAttr->decimals);
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
      if (sparse && ! selectedPairs. contains (getObjName (name1, name2)))
        continue;
          // The graph of comparable objects remains connected (getConnected())
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



void DistTree::loadDissimPrepare (size_t pairs_max,
                                  streamsize target_decimals)
{
  ASSERT (! target);
  ASSERT (ds. objs. empty ());
  
  if (verbose ())
    cout << "Leaf pairs -> data objects ..." << endl;

  ds. objs. reserve (pairs_max);
  target = new RealAttr1 ("Target", ds, target_decimals); 
  obj2leaf1. reserve (ds. objs. capacity ());
  obj2leaf2. reserve (ds. objs. capacity ());
}
  


bool DistTree::addDissim (const string &name1,
                          const string &name2,
                          Real dissim)
{
  if (   isNan (dissim)  // prediction must be large ??
      || ! DM_sp::finite (dissim)  
     )
  {
    cout << name1 << " - " << name2 << ": " << dissim << endl;
    completeDs = false;
    return false;  
  }
    
  const Real mult = dissim2mult (dissim);
  if (nullReal (mult))
  {
    if (! nullReal (dissim))
      cout << name1 << " - " << name2 << ": " << dissim << endl;
    return false;
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
    Unverbose unv;
    Progress prog ((uint) nodes. size (), 100);  // PAR
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



void DistTree::qc () const
{ 
  if (! qc_on)
    return;
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

//ASSERT (ds. objs. size () <= (leaves * (leaves - 1)) / 2);
  const size_t discernables = getDiscernables (). size ();
  ASSERT (discernables <= leaves);
  const size_t pairs_max = (discernables * (discernables - 1)) / 2;
  ASSERT (ds. objs. size () <= pairs_max);
  IMPLY (completeDs, ds. objs. size () == pairs_max);

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
  if (! qc_on)
    return;
  
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
  return *p1 + objNameSeparator + *p2;
}



const DTNode* DistTree::lcaName2node (const string &lcaName) const
{
  ASSERT (! lcaName. empty ());

  string s (lcaName);
  const string name1 = findSplit (s, objNameSeparator);
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
  os << "INPUT:" << endl;
  os << "# Leaves: " << root->getLeavesSize () << endl;
  os << "# Discernable leaves: " << getDiscernables (). size () << endl;
  os << "# Nodes: " << nodes. size () << endl;
  if (! optimizable ())
    return;
  os << "# Dissimilarities: " << ds. objs. size () << " (" << (Real) ds. objs. size () / (Real) dissimSize_max () * 100 << " %)" << endl; 
  os << "Ave. dissimilarity = " << getDissim_ave () << endl;
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
  
  OFStream f (fName);
  static_cast <const DTNode*> (root) -> saveFeatureTree (f, 0);
}



void DistTree::topology2attrs (const List<DiGraph::Node*>& nodes_arg)
{
  ASSERT (optimizable ());

  if (verbose ())
    cout << "Data values ..." << endl;

  for (DiGraph::Node* node : nodes_arg)
    node->inStack = false;

  Progress prog ((uint) ds. objs. size (), 10000);  // PAR
  FOR (size_t, objNum, ds. objs. size ())
  {
    prog ();
    const VectorPtr<TreeNode> path (getPath ( obj2leaf1 [objNum]
                                            , obj2leaf2 [objNum]
                                            )
                                   );
    for (const TreeNode* node : path)
      const_cast<TreeNode*> (node) -> inStack = true;
    for (const DiGraph::Node* node : nodes_arg)
    {
      const DTNode* dtNode = static_cast <const DTNode*> (node);
      CompactBoolAttr1* attr = const_cast <CompactBoolAttr1*> (dtNode->attr);
      ASSERT (attr);
      attr->setCompactBool (objNum, dtNode->inStack);
    }
    for (const TreeNode* node : path)
      const_cast<TreeNode*> (node) -> inStack = false;
  }
}



void DistTree::clearSubtreeLen ()
{
  for (DiGraph::Node* node : nodes)
    static_cast <DTNode*> (node) -> subtreeLen. clear ();
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

  const Real dissim_ave = getDissim_ave ();
  FOR (size_t, objNum, obj2leaf1. size ())
    if (! eqReal ((* const_cast <RealAttr1*> (prediction)) [objNum], lr. predict (objNum), dissim_ave * 1e-2))  // PAR 
    {
      const ONumber on (cout, dissimDecimals, true);
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
  if (! optimizable ())
    return;
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
      dissim2_half      [half2] += mult * sqr (d);
    }
    
  for (const bool half2 : Bool)
  {
    const ONumber on (cout, 6, false);  // PAR
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



void DistTree::getSkipRetain (DTNode* &toSkip,
                              DTNode* &toRetain)
{
  toSkip   = nullptr;  
  toRetain = nullptr; 
  const VectorPtr<DiGraph::Node> rootChildren (root->getChildren ());
  ASSERT (rootChildren. size () >= 2);
  if (rootChildren. size () == 2)  // root is transient in an undirected tree
  {
    toSkip   = const_static_cast <DTNode*> (rootChildren [0]);
    toRetain = const_static_cast <DTNode*> (rootChildren [1]);
    toRetain->len += toSkip->len;
    toSkip->len = 0;
  }
  ASSERT ((bool) toSkip == (bool) toRetain);
}



#if 0
bool DistTree::optimizeLenAll ()  
// Global optimization of DTNode::len
// Update: DTNode::len
// Output: prediction, absCriterion
// Invokes: setPrediction(), setAbsCriterion()
// To be followed by: finishChanges()
// Time: O(p log^2(n) + n^3)
{
	ASSERT (optimizable ());

  if (verbose (1))
    cout << "Optimizing arc lengths ..." << endl;

  DTNode* toSkip = nullptr;  
  DTNode* toRetain = nullptr; 
  getSkipRetain (toSkip, toRetain);

  // DTNode::index, dtNodes, arcs
  for (DiGraph::Node* node : nodes)
    static_cast <DTNode*> (node) -> index = NO_INDEX;
  VectorPtr<DTNode> dtNodes;  dtNodes. reserve (nodes. size ());
  size_t arcs = 0;
  for (DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <DTNode*> (node);
    if (   dtNode != root
        && dtNode != toSkip
        && ! dtNode->inDiscernable ()
       )
    {
      dtNodes << dtNode;
      const_cast <DTNode*> (dtNode) -> index = arcs;
      arcs++;
    }
    IMPLY (dtNode->inDiscernable (), dtNode->len == 0);
  }
  ASSERT (arcs == dtNodes. size ());

  // DTNode::dissimSum, matr
  Matrix matr (false, arcs, arcs /*+ 1*/);
  matr. putAll (0);
  for (DiGraph::Node* node : nodes)
    static_cast <DTNode*> (node) -> dissimSum = 0;
  FOR (size_t, objNum, ds. objs. size ())
  {
    const Real dissim = (*target) [objNum];    
    const Real mult = ds. objs [objNum] -> mult;
    const VectorPtr<TreeNode> path (getPath ( obj2leaf1 [objNum]
                                            , obj2leaf2 [objNum]
                                            )
                                   );
    for (const TreeNode* node1 : path)
    {
      DTNode* dtNode1 = const_static_cast <DTNode*> (node1);
      if (dtNode1->index == NO_INDEX)
        continue;
      dtNode1->dissimSum += dissim * mult;
      for (const TreeNode* node2 : path)
      {
        const DTNode* dtNode2 = static_cast <const DTNode*> (node2);
        if (dtNode2->index != NO_INDEX)
          matr. putInc ( false
                       , dtNode1->index
                       , dtNode2->index
                       , mult
                       );
      }
    }
  }    
    
  // Cf. L2LinearNumPrediction::solveUnconstrained()
  Matrix betaCovariance (matr);
  if (   ! betaCovariance. inverse (). get ()
      || ! matr. checkInverse (betaCovariance, false)
     )
    return false;
        
  MVector xy (arcs);
  xy. putAll (NAN);
  for (const DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <const DTNode*> (node);
    if (dtNode->index != NO_INDEX)
    {
      xy [dtNode->index] = dtNode->dissimSum;
    //matr. put (false, dtNode->index, arcs, dtNode->dissimSum);
    }
  }
  ASSERT (xy. defined ());

  MVector beta (arcs);
#if 1
  beta. multiply ( false
                 , betaCovariance, false
                 , xy,             false
                 );
#else
  beta. solveSystem (false, 0, matr, false);
#endif
        
  FOR (size_t, i, dtNodes. size ())
    const_cast <DTNode*> (dtNodes [i]) -> len = max (0.0, beta [i]);
      
  if (toSkip)
  {
    toSkip->len = toRetain->len / 2;
    toRetain->len = toSkip->len;
  }

  setPrediction ();
  setAbsCriterion ();  
  
  return true;
}
#endif



#if 0
void DistTree::quartet2arcLen ()  
// Assumes: Obj::mult = 1
// Output: DTNode::dissimSum, DTNode::len
// Requires: completeDs, # leaves > 2
// Invokes: setLeaves()
// Time: O(p log(n) + n)
{
	ASSERT (optimizable ());
	ASSERT (completeDs);


  if (verbose (1))
    cout << "Optimizing arc lengths by quartets ..." << endl;


  setLeaves ();
  const size_t allLeaves = root->leaves;
  ASSERT (allLeaves == name2leaf. size ());
  ASSERT (allLeaves > 2);
  

  for (DiGraph::Node* node : nodes)
  {
    DTNode* dtNode = static_cast <DTNode*> (node);
    dtNode->dissimSum = 0;
    dtNode->dissimWeightedSum = 0;
  }
  FOR (size_t, objNum, ds. objs. size ())
  {
    const Real dissim = (*target) [objNum];    
    const VectorPtr<TreeNode> path (getPath ( obj2leaf1 [objNum]
                                            , obj2leaf2 [objNum]
                                            )
                                   );
    FOR (size_t, i, path. size ())
    {
      DTNode* dtNode = const_static_cast <DTNode*> (path [i]);
      ASSERT (dtNode != root);
      dtNode->dissimSum += dissim;
      
      const TreeNode* parent = dtNode->getParent ();
      ASSERT (parent);
      const TreeNode* prev = nullptr;
      const TreeNode* next = nullptr;
      {
        if (i >= 1)
          prev = path [i - 1];
        if (i + 1 < path. size ())
          next = path [i + 1];
        const TreeNode* prevPrev = i >= 2 ? path [i - 2] : nullptr;
        const TreeNode* nextNext = (i + 2 < path. size ()) ? path [i + 2] : nullptr;
        if (prev && prev->getParent () != dtNode)
        {
          swap (prev,     next);
          swap (prevPrev, nextNext);
        }
        if (next && parent != next)
        {
          ASSERT (parent == next->getParent ());
          if (parent == root && root->arcs [false]. size () == 2)
            next = nextNext;
        }
      }
      
      ASSERT (allLeaves > dtNode->leaves);
      
      size_t prevPairs = 2;
      if (prev)
      {
        ASSERT (prev->getParent () == dtNode || prev->getParent () == parent);
        ASSERT (prev->leaves);
        ASSERT (dtNode->leaves > prev->leaves);
        prevPairs = dtNode->leaves - prev->leaves;
      }
      
      size_t nextPairs = 2;
      if (next)
      {
        ASSERT (next->getParent ());
        ASSERT (allLeaves > next->leaves);
        if (parent == next)
        {
          ASSERT (next->leaves > dtNode->leaves);
          nextPairs = next->leaves - dtNode->leaves;
        }
        else
        {
          ASSERT (parent == next->getParent () || parent == next->getParent () -> getParent ());
          ASSERT (allLeaves - dtNode->leaves > next->leaves);
          nextPairs = (allLeaves - dtNode->leaves) - next->leaves;
        }
      }

      ASSERT (prevPairs);
      ASSERT (nextPairs);
      dtNode->dissimWeightedSum += dissim * (Real) prevPairs * (Real) nextPairs;
    }
  }


  for (DiGraph::Node* node : nodes)   
  {
    DTNode* a = static_cast <DTNode*> (node);
    if (a->inDiscernable ())
      continue;
    if (a == root)
      continue;
      
  	ASSERT (a->leaves);

    const DTNode* b = static_cast <const DTNode*> (a->getParent ());
    ASSERT (b);
    if (b == root && root->arcs [false]. size () == 2)
      b = static_cast <const DTNode*> (root->getOtherChild (a));
    ASSERT (b);
    ASSERT (a != b);

    // Estimatating the length of arc (a,b)
      
  	ASSERT (allLeaves > a->leaves);
  	const size_t bLeaves = allLeaves - a->leaves;
  	
    Real belowSum = 0;  // Times 2
    size_t belowPairs = 0;
  	for (const Arc* arc : a->arcs [false])
  	{
  	  const DTNode* child = static_cast <const DTNode*> (arc->node [false]);
  	  belowSum += child->dissimSum;
  	  ASSERT (child->leaves);
  	  ASSERT (a->leaves > child->leaves);
  	  belowPairs += child->leaves * (a->leaves - child->leaves);
  	}
  	if (a->leaves > 1)
  	{
  	  belowSum -= a->dissimSum;
    	ASSERT (even (belowPairs));
    	belowPairs /= 2;
  	}
  	else
  	  belowPairs = 1;
  	ASSERT (! negative (belowSum));
  	IMPLY (a->leaves == 1, belowSum == 0);
  	ASSERT (belowPairs);

    Real aboveSum = 0;  // Times 2
    size_t abovePairs = 0;
  	for (const Arc* arc : b->arcs [false])
  	{
  	  const DTNode* child = static_cast <const DTNode*> (arc->node [false]);
  	  if (child != a)
  	  {
  	    aboveSum += child->dissimSum;
    	  ASSERT (child->leaves);
    	  ASSERT (bLeaves > child->leaves);
    	  abovePairs += child->leaves * (bLeaves - child->leaves);
  	  }
  	}
	  if (b == a->getParent () && b != root)
	  {
  	  aboveSum += b->dissimSum;
  	  ASSERT (allLeaves > b->leaves);
  	  ASSERT (b->leaves > a->leaves);
  	  abovePairs += (allLeaves - b->leaves) * (b->leaves - a->leaves); 
  	}
  	aboveSum -= a->dissimSum;
  	ASSERT (! negative (aboveSum));
  	ASSERT (even (abovePairs));
  	abovePairs /= 2;
  	
    // quartet
    const Real len_old = a->len;
  	a->len = 0.25 * (  a->dissimWeightedSum / ((Real) belowPairs * (Real) abovePairs) 
  	                 - belowSum / (Real) belowPairs
  	                 - aboveSum / (Real) abovePairs
  	                );
  	maximize (a->len, 0.0);
  	
  	if (b != a->getParent ())
  	{
  	  ASSERT (a->getParent () == root); 
  	  ASSERT (b->getParent () == root);
  	  ASSERT (root->arcs [false]. size () == 2);
  	  a->len /= 2;
  	  const_cast <DTNode*> (b) -> len = a->len;
  	}

  #if 0
  	const ONumber on (cout, dissimDecimals, true); 
  	cout << a->getName () << ": " << len_old << " -> " << a->len << " " << (a->len - len_old) / len_old * 100 << " %" << endl;  
  #endif
  }
  

  setPrediction ();
  setAbsCriterion ();  
}
#endif



namespace
{

bool DTNode_len_strictlyLess (const DTNode* a,
                              const DTNode* b)
{
  ASSERT (a);
  ASSERT (b);
  ASSERT (! isNan (a->len));
  ASSERT (! isNan (b->len));
  return a->len > b->len;  
}

}



bool DistTree::optimizeLenArc ()
{
  ASSERT (optimizable ());
  
  if (verbose (1))
    cout << "Optimizing arc lengths at each arc ..." << endl;
  
  DTNode* toSkip = nullptr;  
  DTNode* toRetain = nullptr; 
  getSkipRetain (toSkip, toRetain);

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
  Common_sp::sort (dtNodes, DTNode_len_strictlyLess);

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
  /*const bool solved =*/ lr. solveUnconstrainedFast (prediction, true, 10, 0.01);  // PAR  // Time = O(p n) 
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



namespace
{
  
struct NodeStar
{
  const Steiner* st;
  VectorPtr<DiGraph::Node> star;
  Real lenSum;
  
  explicit NodeStar (const Steiner* st_arg)
    : st     (st_arg)
    , star   (st->getChildren ())
    , lenSum (0)
    {
      ASSERT (st);
      ASSERT (! st->isTransient ());
      star << st;
      ASSERT (star. size () >= 3);
      for (const DiGraph::Node* node : star)
        lenSum += static_cast <const DTNode*> (node) -> len;
      ASSERT (lenSum >= 0);
    }
   
  static bool strictlyLess (const NodeStar &a,
                            const NodeStar &b)
    {
      ASSERT (! isNan (a. lenSum));
      ASSERT (! isNan (b. lenSum));
      return a. lenSum > b. lenSum;  
    }
};
  
}



void DistTree::optimizeLenNode ()  
{
  ASSERT (optimizable ());
  
  if (verbose (1))
    cout << "Optimizing arc lengths at each node ..." << endl;

//cout << "absCriterion (before optimizeLenLocal) = " << absCriterion << endl;  


  Vector<NodeStar> dtNodes;  dtNodes. reserve (nodes. size ());
  for (const DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <const DTNode*> (node);
    if (   dtNode == root
        || dtNode->asLeaf () 
      //|| dtNode->inDiscernable ()
        || ! dtNode->childrenDiscernable ()
       )
      continue;
    dtNodes << NodeStar (dtNode->asSteiner ());
  }
  Common_sp::sort (dtNodes, NodeStar::strictlyLess);  


  Progress prog ((uint) countInteriorNodes ());  // Upper bound
  Space1<NumAttr1> sp (ds, false); 
  Vector<Real> lenOld;  
  for (const NodeStar& ns : dtNodes)
  {
  #ifndef NDEBUG    
    const Real absCriterion_old = absCriterion;  
  #endif
    
    const VectorPtr<DiGraph::Node>& star = ns. star;

    // lenOld, *target_new
    lenOld. clear (); lenOld. reserve (star. size ());    
    for (const DiGraph::Node* node : star)
    {
      const DTNode* n = static_cast <const DTNode*> (node);
      ASSERT (! n->inDiscernable ());
      lenOld << n->len;
      const_cast <DTNode*> (n) -> len = 0;
    }
    ASSERT (lenOld. size () == star. size ());
    FOR (size_t, objNum, ds. objs. size ())
    {
      const VectorPtr<TreeNode> path (getPath ( obj2leaf1 [objNum]
                                              , obj2leaf2 [objNum]
                                              )
                                     );
      (* const_cast <RealAttr1*> (target_new)) [objNum] = (*target) [objNum] - path2prediction (path);
    }
    
    sp. clear ();  sp. reserve (star. size ());
    for (const DiGraph::Node* dtNode1 : star)
      sp << static_cast <const DTNode*> (dtNode1) -> attr;

    // lr
    Unverbose unv;
    if (verbose ())
      cout << "Linear regression ..." << endl;
    L2LinearNumPrediction lr (dsSample, sp, *target_new);
    ASSERT (lr. beta. size () == star. size ());
    FOR (size_t, attrNum, lr. beta. size ())
      lr. beta [attrNum] = lenOld [attrNum];
    const bool solved = lr. solveUnconstrainedFast (nullptr, true, 10, 0.01);  // PAR
    if (verbose ())  
      lr. qc ();
  
    // DTNode::len
    FOR (size_t, attrNum, star. size ())
      const_static_cast <DTNode*> (star [attrNum]) -> len = solved ? lr. beta [attrNum] : lenOld [attrNum];
    if (qc_on && verbose ())
    {
      setPrediction ();
      setAbsCriterion ();  
      ASSERT (leReal (absCriterion, absCriterion_old));
    }

    if (solved)
      prog (real2str (lr. absCriterion, 6));  // PAR
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
  if (qc_on && verbose ())
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
    if (qc_on && verbose ())   
    {
      checkPrediction ();
      checkAbsCriterion ("optimize");
    }
    if (verbose (1))
      cout. flush ();
  }
}



#if 0
void DistTree::optimizeAdd (bool sparse,
                            const string &output_tree)
// Input: dissimDs, dissimAttr
// Requires: (bool)dissimAttr
// Invokes: addDissim(), optimizeSubgraph() if leafRelCriterion is large, root->findClosestNode(), DTNode::selectRepresentative()
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
        const ONumber on (cerr, dissimDecimals, true); 
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
        // Try NJ ??
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
      
      Set<const Leaf*> representatives;  // size = O(log n)
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
        }
    
        FOR_START (size_t, objNum, objNum_start, ds. objs. size ())
        {
          const VectorPtr<TreeNode> path (getPath ( obj2leaf1 [objNum]
                                                  , obj2leaf2 [objNum]
                                                  )
                                         );
          ASSERT (path. contains (g));
          // Use Node::inStack ??
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
        const ONumber on (cerr, 1, false);
        cerr << "rel.absCriterion = " << leafRelCriterion << "  ";
      }
      
      if (leafRelCriterion > 2)  // PAR
      {
        Unverbose unv;
      //VectorPtr<DTNode> boundary;  boundary. reserve (nodes. size ());
        optimizeSubgraph (const_static_cast <Steiner*> (g->getParent ()) /*, boundary*/);
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
    cout << "Optimizing arc lengths ..." << endl;
  EXEC_ASSERT (optimizeLenArc ());   
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
#endif



Real DistTree::optimizeSubgraph (const Steiner* center)  
{
  ASSERT (center);
  ASSERT (center->graph);
  ASSERT (& center->getTree () == this);

  VectorPtr<TreeNode> area;
  const DTNode* area_root = nullptr;
  Node2Node new2old;  // Initially: newLeaves2boundary
  DistTree tree (center, areaRadius_std, area, area_root, new2old);
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

  const size_t leaves = tree. name2leaf. size ();
  {
    Unverbose unv;
    if (leaves > 3)
    {
      EXEC_ASSERT (tree. optimizeLenArc ());
      tree. finishChanges (); 
      tree. optimizeLenNode ();  
      tree. finishChanges (); 
      //
      tree. optimizeIter (string ());
        // tree.neighborJoin(), tree.optimizeSubgraphs() if tree is large ??
          // allows using mdsTree.sh
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
    cout << "Subtree: ";
    tree. qc ();
    tree. qcAttrs ();
    tree. reportErrors (cout);
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

  // tree
  // new2old: newBoundary2oldBoundary
  if (const DTNode* dtNode = static_cast <const DTNode*> (findPtr (boundary2new, static_cast <const DiGraph::Node*> (area_root))))
    if (const Leaf* leaf = dtNode->asLeaf ())
    {
      if (verbose ())
        cout << "Re-rooting ..." << endl;
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

  // nodes: delete
  // center may be delete'd
  for (const TreeNode* node : area)
    if (! contains (boundary2new, static_cast <const DiGraph::Node*> (node)))
    {
      DTNode* dtNode = const_static_cast <DTNode*> (node);
      delete dtNode->attr;
      delete dtNode;
    }
    // Between boundary TreeNode's there are no Arc's

  // nodes, new2old[]: new
  List<DiGraph::Node*> addedNodes;
  for (const DiGraph::Node* node_new : tree. nodes)
  {
    DTNode* node = const_static_cast <DTNode*> (findPtr (new2old, node_new));
    if (! node)
    {
      node = new Steiner (*this, nullptr, NAN);
      node->addAttr ();
      new2old [node_new] = node;
      addedNodes << node;
      node->stable = true;
    }
    ASSERT (node);
    if (node_new != tree. root)
      node->len = static_cast <const DTNode*> (node_new) -> len;
  }
  ASSERT (new2old. size () == tree. nodes. size ());

  if (contains (boundary2new, static_cast <const DiGraph::Node*> (center)))
  {
    ASSERT (center->graph);
    ASSERT (& center->getTree () == this);
    const_cast <Steiner*> (center) -> stable = true;
  }
  
  borrowArcs (new2old, true);  // Arc's are not parallel

  setRoot ();

  topology2attrs (addedNodes); 
  setPrediction ();
  finishChanges ();
  setAbsCriterion ();
  if (qc_on && verbose ())
  {
    cout << "absCriterion = " << absCriterion << endl;
    checkPrediction ();
    checkAbsCriterion ("optimizeSubgraph");
  }

  if (verbose (1))
    reportErrors (cout);
  	
  if (greaterReal (absCriterion, absCriterion_old, 1e-5)) // PAR
  {
    const ONumber on (cout, 9, true);
    cout << absCriterion << " " << absCriterion_old << endl;
    ERROR;
  }
  
  return leafLen_min;
}



void DistTree::optimizeSubgraphs ()
{
  for (DiGraph::Node* node : nodes)
    static_cast <DTNode*> (node) -> stable = false;

  Progress prog;
  for (;;)
  {
    size_t steiners = 0;
    size_t stables  = 0;
    const Steiner* center = nullptr;
    for (const DiGraph::Node* node : nodes)
      if (const Steiner* newSt = static_cast <const DTNode*> (node) -> asSteiner ())
      {
        steiners++;
        if (newSt->stable)
          stables++;
        else
          if (! newSt->getParent () || static_cast <const DTNode*> (newSt->getParent ()) -> stable)
          {
            if (! center)
              center = newSt;
          }
      }
    if (! center)
      break;
    {
      Unverbose unv;
      optimizeSubgraph (center);
    }
    {
      ostringstream oss;
      const ONumber on (oss, 6, false);  // PAR
      oss << stables << '/' << steiners << ' ' << absCriterion;
      prog (oss. str ());  
    }
  }

#if 0
  // Almost no improvement
  EXEC_ASSERT (optimizeLenArc ()); 
  finishChanges ();
  optimizeLenNode ();  
  finishChanges ();
  if (verbose (1))
    reportErrors (cout);
#endif
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
    Common_sp::sort (changes, Change::strictlyLess);	
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
  
  	  if (qc_on && verbose ())
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
        if (qc_on && verbose ())
    	  	checkAbsCriterion ("restore");
  	  }
  	  else
  	  {
       	ch->commit ();
       	commits++;
        setAbsCriterion ();
    	  if (qc_on && verbose ())
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
  if (verbose (1))
    cout << "Improvement = " << improvement << endl;
  ASSERT (geReal (improvement, - absCriterion_delta));
  ASSERT ((bool) commits == (bool) improvement);

  if (commits)
  {
    topology2attrs (nodes); 
    EXEC_ASSERT (optimizeLenArc ()); 
    finishChanges ();
    optimizeLenNode ();  
    finishChanges ();
    if (verbose (1))
      reportErrors (cout);
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
  if (qc_on && verbose ())
  {
    ch->print (cout); 
  	checkAbsCriterion ("tryChange");
  }
  
	if (Change::strictlyLess (ch, bestChange))
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



size_t DistTree::finishChanges ()
{
  const size_t n = deleteLenZero ();
  if (n && verbose ())
    cout << "# Nodes with zero arcs deleted = " << n << endl;

  toDelete. deleteData ();
  
  return n;
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



Set<string> DistTree::selectPairs ()
{
  setReprLeaves ();  
  ds. setName2objNum ();
  
  Set<string> selectedPairs;
 	VectorPtr<DTNode> descendants;  descendants. reserve (2 * powInt (2, (uint) sparsingDepth));  // PAR
  for (const DiGraph::Node* node : nodes)
    if (const Leaf* leaf = static_cast <const DTNode*> (node) -> asLeaf ())
    {
      ASSERT (leaf == leaf->reprLeaf);
      const DTNode* ancestor = leaf;
      while (ancestor)
      {
      	descendants. clear ();
      	ancestor->getDescendants (descendants, sparsingDepth); 
      	for (const DTNode* descendant : descendants)
      	{
      	  ASSERT (descendant->reprLeaf);
      	  if (leaf->getDiscernable () != descendant->reprLeaf->getDiscernable ())
      	  {
      	    const string objName (getObjName (leaf->name, descendant->reprLeaf->name));
      	    if (ds. getName2objNum (objName) == NO_INDEX)
        	    selectedPairs << objName;
        	}
      	}
        ancestor = static_cast <const DTNode*> (ancestor->getParent ());
      }
    }          
  return selectedPairs;
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
    if (dtNode->inDiscernable ())
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
      if (geReal (relErr, relErr_min)) 
      {
        os << leaf->getName () << '\t' << relErr << endl;
        n++;
      }
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



void DistTree::printDissim (ostream &os) const
{
  ASSERT (optimizable ());

  const ONumber on (os, dissimDecimals, true);
  FOR (size_t, objNum, ds. objs. size ())
  {
    ASSERT (ds. objs [objNum] -> mult);
    os         << obj2leaf1 [objNum] -> name
       << '\t' << obj2leaf2 [objNum] -> name
       << '\t' << (*target) [objNum]
       << endl;
  }
}




//////////////////////////////////////////////////////////////////////

// NewLeaf::Location

void NewLeaf::Location::qc () const
{
  if (! qc_on)
    return;

  ASSERT (anchor);
  ASSERT (! anchor->inDiscernable ());
  if (anchor == anchor->getTree (). root)
  {
    ASSERT (isNan (anchor->len));
    ASSERT (arcLen == 0);
  }
  else
  {
    ASSERT (anchor->len > 0);
    ASSERT (anchor->len < INF);
    ASSERT (arcLen >= 0);
    ASSERT (arcLen <= anchor->len);
  }

  ASSERT (leafLen >= 0);
  ASSERT (leafLen < INF);
  
  ASSERT (absCriterion_delta >= 0);  
}


    

// NewLeaf::Leaf2dissim

NewLeaf::Leaf2dissim::Leaf2dissim (const Leaf* leaf_arg,
                                   Real dissim_arg,
                                   const DTNode* anchor)
: leaf (leaf_arg)
, dissim (dissim_arg)
, mult (dissim2mult (dissim_arg))
, dist_hat (0)
, leafIsBelow (false)
{ 
  ASSERT (anchor);
  
  ASSERT (leaf);
  
  ASSERT (dissim >= 0);
  ASSERT (dissim < INF);
  
  ASSERT (mult >= 0);
  
  // dist_hat, leafIsBelow
  { 
    const DTNode* node = leaf;
    for (;;)
    { 
      ASSERT (node);
      if (node == anchor)
      {
        ASSERT (! leafIsBelow);
        leafIsBelow = true;
        break;
      }
      if (node == node->getTree (). root)
        break;
      dist_hat += node->len;
      node = static_cast <const DTNode*> (node->getParent ());
    }
    if (! leafIsBelow)
    {
      node = anchor;
      for (;;)
      { 
        ASSERT (node);
        if (node == node->getTree (). root)
          break;
        dist_hat += node->len;
        node = static_cast <const DTNode*> (node->getParent ());
      }
    }
  }  
  ASSERT (dist_hat >= 0);
  ASSERT (dist_hat < INF);
}




// NewLeaf

NewLeaf::NewLeaf (const DistTree &tree_arg,
                  const string &dataDir_arg,
                  const string &name_arg,
                  bool init)
: Named (name_arg)
, tree (tree_arg)
, dataDir (dataDir_arg)
{
  ASSERT (isRight (dataDir, "/"));
  ASSERT (! name. empty ());
  
  if (fileExists (getRequestFName ()))
    throw runtime_error ("File \"" + getRequestFName () + "\" exists");
    
  location. anchor = static_cast <const DTNode*> (tree. root);
  location. leafLen = 0;
  location. arcLen = 0;

  if (init)
  { 
    if (fileExists (getLeafFName ()))
      throw runtime_error ("File \"" + getLeafFName () + "\" exists");
    if (fileExists (getDissimFName ()))
      throw runtime_error ("File \"" + getDissimFName () + "\" exists");
    { OFStream of (getDissimFName ()); }  // destructor creates an empty file
  }
  else
  {
  #if 0
    {
      ifstream f (getLeafFName ());
      string name1, lcaName, tmp;
      f >> name1 >> lcaName >> location. leafLen >> location. arcLen >> tmp;
      ASSERT (name1 == name);
      ASSERT (tmp. empty ());
      location. anchor = tree. lcaName2node (lcaName);
    }
  #endif
    if (! fileExists (getLeafFName ()))
      throw runtime_error ("File \"" + getLeafFName () + "\" does not exist");
    if (! fileExists (getDissimFName ()))
      throw runtime_error ("File \"" + getDissimFName () + "\" does not exist");
    {
      LineInput f (getDissimFName ());
      string name1, name2;
      Real dissim;
      while (f. nextLine ())
        try
        {
          istringstream iss (f. line);
          iss >> name1 >> name2 >> dissim;
          ASSERT (iss. eof ());
          ASSERT (name1 != name2);
          if (name1 != name)
            throw runtime_error ("First object in " + getDissimFName () + " must be " + name);
          if (dissim < 0)
            throw runtime_error ("Dissimilarity must be non-negative");
          if (! DM_sp::finite (dissim))
            throw runtime_error ("Dissimilarity must be finite");
          leaf2dissims << Leaf2dissim (findPtr (tree. name2leaf, name2), dissim, location. anchor);
        }          
        catch (...)
        {
          cout << f. line << endl;
          throw;
        }
      Common_sp::sort (leaf2dissims);
    #ifndef NDEBUG
      // Check uniqueness
      const Leaf2dissim* prev = nullptr;
      for (const Leaf2dissim& ld : leaf2dissims)
      {
        IMPLY (prev, prev->leaf < ld. leaf);
        prev = & ld;
      }
    #endif
      optimize ();
    }
  }

  saveLeaf ();
  saveRequest ();
}



void NewLeaf::qc () const
{
  if (! qc_on)
    return;
  Named::qc ();
    
  location. qc ();    
  ASSERT (& location. anchor->getDistTree () == & tree);
}



void NewLeaf::saveRequest () const
{	
  Set<const Leaf*> requested;
  {
    // Cf. DistTree::selectPairs()
    VectorPtr<DTNode> descendants;  descendants. reserve (2 * powInt (2, (uint) sparsingDepth));  // PAR
    const DTNode* ancestor = location. anchor;
    while (ancestor)
    {
    	descendants. clear ();
    	ancestor->getDescendants (descendants, sparsingDepth); 
    	for (const DTNode* descendant : descendants)
    	{
    	  ASSERT (descendant->reprLeaf);
    	  const Leaf2dissim ld (descendant->reprLeaf);
    	  if (leaf2dissims. binSearch (ld) == NO_INDEX)
    	    requested << descendant->reprLeaf;
    	}
      ancestor = static_cast <const DTNode*> (ancestor->getParent ());
    }
  }
  
  {
    OFStream of (getRequestFName ());
    for (const Leaf* leaf : requested)
      of << name << '\t' << leaf->name << endl;
  }
}



void NewLeaf::optimize ()
{
	for (const Leaf2dissim& ld : leaf2dissims)
	  if (ld. dissim == 0)
	  {
	    location. anchor = ld. leaf->getDiscernable ();
	    location. leafLen = 0;
	    location. arcLen = 0;
	    location. absCriterion_delta = 0;
	    return;
	  }


  setLocation ();
  Location location_best (location);
  Vector<Leaf2dissim> leaf2dissims_best (leaf2dissims);
  while (location. arcLen == 0 && location. anchor->childrenDiscernable ())
  {
    VectorPtr<DiGraph::Arc> arcs;  arcs. reserve (location. anchor->arcs [false]. size ());
    insertAll (arcs, location. anchor->arcs [false]);
    bool allChildrenTried = true;
  	for (const DiGraph::Arc* arc : arcs)
  	{
  	  const DTNode* child = static_cast <const DTNode*> (arc->node [false]);
      ASSERT (! child->inDiscernable ());
  	  const Location location_old (location);
  	  const Vector<Leaf2dissim> leaf2dissims_old (leaf2dissims);
  	  if (descend (child))
  	  {
    	  setLocation ();
    	  if (minimize (location_best. absCriterion_delta, location. absCriterion_delta))
    	  {
    	    location_best = location;
          leaf2dissims_best = leaf2dissims;
    	  }
    	}
    	else
    	  allChildrenTried = false;
  	  location = location_old;
  	  leaf2dissims = leaf2dissims_old;
  	}
    if (! allChildrenTried || location_best. anchor == location. anchor)
      break;
    location = location_best;
    leaf2dissims = leaf2dissims_best;
    if (verbose ())  
    {
  	  cout << "Location: "; 
  	  location. saveText (cout); 
  	  cout << '\t' << "anchor = " << location. anchor->len 
  	       << '\t' << "# dissims = " << leaf2dissims. size ()
  	       << '\t' << "criterion = " << location. absCriterion_delta 
  	       << endl;  
  	}
  }
}



void NewLeaf::setLocation ()
{
	Real mult_sum = 0;
	Real u_avg = 0;	
	Real delta_avg = 0;
	for (const Leaf2dissim& ld : leaf2dissims)
	{
	  ASSERT (ld. dissim > 0);
	  mult_sum  += ld. mult;
	  u_avg     += ld. mult * ld. getU ();
	  delta_avg += ld. mult * ld. getDelta ();
	  
	}	
	ASSERT (positive (mult_sum));
	u_avg     /= mult_sum;
	delta_avg /= mult_sum;
	
	Real u_var = 0;
	Real delta_u_cov = 0;
	for (const Leaf2dissim& ld : leaf2dissims)
	{
	  delta_u_cov += ld. mult * (ld. getDelta () - delta_avg) * (ld. getU () - u_avg);
	  u_var       += ld. mult * sqr (ld. getU () - u_avg);
	}
	u_var       /= mult_sum;
	delta_u_cov /= mult_sum;  	

  if (location. anchor == tree. root)
  {
  	ASSERT (u_var == 0);
    location. arcLen = 0;
  }
  else
  {
  	ASSERT (u_var > 0);
  	location. arcLen = min (max (0.0, delta_u_cov / u_var), location. anchor->len);
  }

	location. leafLen = max (0.0, delta_avg - u_avg * location. arcLen);
	
	location. absCriterion_delta = 0;
	for (const Leaf2dissim& ld : leaf2dissims)
	  location. absCriterion_delta += ld. mult * sqr (ld. getEpsilon (*this));
	ASSERT (location. absCriterion_delta >= 0);
}



bool NewLeaf::descend (const DTNode* anchor_new)
{
  ASSERT (anchor_new);
  ASSERT (anchor_new->getParent () == location. anchor);
  
  location. anchor = anchor_new;
  
  bool hasLeaves = false;
  for (Leaf2dissim& ld : leaf2dissims)
  {
    if (ld. leafIsBelow)
      ld. leafIsBelow = ld. leaf->descendantOf (location. anchor);
    ld. dist_hat -= location. anchor->len * ld. getU ();
    if (ld. leafIsBelow)
      hasLeaves = true;
  }
  
  return hasLeaves;
}



}





/* TO DO ??

not solved LinearRegression
  
non-stability of results

*/
