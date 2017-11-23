// distTree.cpp

#undef NDEBUG
#include "../common.inc"

#include "distTree.hpp"

#include "../dm/prediction.hpp"



namespace DistTree_sp
{



Chronometer chron_getBestChange ("getBestChange");
Chronometer chron_tree2subgraph ("tree2subgraph");
Chronometer chron_subgraphOptimize ("subgraphOptimize");
Chronometer chron_subgraph2tree ("subgraph2tree");



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

  if (graph)
  {
    if ((bool) getParent () != ! isNan (len))
    {
      cout << getName () << endl;
      getTree (). print (cout);
      ERROR;
    }
    if (! childrenDiscernable ())
    {
      ASSERT (! inDiscernable ());
    	for (const DiGraph::Arc* arc : arcs [false])
    	  ASSERT (static_cast <const DTNode*> (arc->node [false]) -> inDiscernable ());
    }
  //IMPLY (getDistTree (). optimizable (), attr);
    ASSERT (getDistTree (). nodeAttrExist == (bool) attr);
    if (attr)
      attr->qc ();
  }
  else
    ASSERT (! attr);
  
  IMPLY (! isNan (len), len >= 0);	  
  IMPLY (paths, errorDensity >= 0);
}



void DTNode::saveContent (ostream& os) const
{ 
  {
    const ONumber oLen (os, dissimDecimals, true);
    os << "len=" << len;
  }
  const ONumber oNum (os, 6, true);  // PAR
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
  if (! tree. nodeAttrExist)
    return;
  tree. attrNum_max++;
  ASSERT (tree. attrNum_max);
  const string name ("v" + toString (tree. attrNum_max));
  attr = new CompactBoolAttr1 (name, tree. ds);
  attr->qc ();
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
  os << s << ": " << /*" t=" << (isNan (len) ? 0 : len) <<*/ "  C=0  dC=+0-0" << endl;
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



void DTNode::setGlobalLenDown (bool topological,
                               DTNode* &bestDTNode,
                               Real &bestDTNodeLen_new,
                               WeightedMeanVar &bestGlobalLen)
{
  if (const DTNode* parent_ = static_cast <const DTNode*> (getParent ())) 
  {
    WeightedMeanVar parentSubtreeLen (parent_->subtreeLen /*global len*/);  
    {
      WeightedMeanVar subtreeLen1;
      if (topological)
      {
        subtreeLen1 = subtreeLen;
        subtreeLen1. addValue (1);
      }
      else
        subtreeLen1. add (subtreeLen. getMean () + len, subtreeLen. weights + len);
      parentSubtreeLen. subtract (subtreeLen1);  
    }
    // Subtree goes upward from *parent_
      
    WeightedMeanVar globalLen;
    Real lenDelta = NAN;
    if (topological)
    {
      lenDelta = len / 2;
      globalLen. add (      subtreeLen);
      globalLen. add (parentSubtreeLen);
      globalLen. addValue (1);
    }
    else
    {
      // Optimal
      lenDelta =   len / 2 
                 + (  (parentSubtreeLen. getMean () + parentSubtreeLen. weights)
                    - (      subtreeLen. getMean () +       subtreeLen. weights)
                   ) / 4;
      maximize (lenDelta, 0.0);
      minimize (lenDelta, len);
      globalLen. add (      subtreeLen. getMean () + lenDelta,               subtreeLen. weights + lenDelta);
      globalLen. add (parentSubtreeLen. getMean () + (len - lenDelta), parentSubtreeLen. weights + (len - lenDelta));
    }
    ASSERT (! isNan (lenDelta));
    
    if (   ! bestDTNode 
        || bestGlobalLen. getMean () > globalLen. getMean ()
       )
    {
      bestDTNode = this;
      bestDTNodeLen_new = lenDelta;
      bestGlobalLen = globalLen;
    }
    
    if (topological)
    {
      parentSubtreeLen. addValue (1);
      subtreeLen. add (parentSubtreeLen);
    }
    else
      subtreeLen. add (parentSubtreeLen. getMean () + len, parentSubtreeLen. weights + len);
    // global len
  }


	for (const DiGraph::Arc* arc : arcs [false])
	  static_cast <DTNode*> (arc->node [false]) -> setGlobalLenDown (topological, bestDTNode, bestDTNodeLen_new, bestGlobalLen);
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
#endif



void DTNode::setLca ()
{
  ASSERT (getDisjointCluster () == this);
  ASSERT (! tarjanLca);    

  tarjanLca = this;
}




// Steiner


Steiner::Steiner (DistTree &tree,
      	          Steiner* parent_arg,
      	          Real len_arg)
: DTNode (tree, parent_arg, len_arg)  // DistTree must be declared
{}



void Steiner::qc () const
{
  if (! qc_on)
    return;
	DTNode::qc ();
	  
	ASSERT (! isLeaf ());	
//IMPLY (reprLeaf, getChildren (). contains (reprLeaf));
}



void Steiner::setSubtreeLenUp (bool topological)
{
  subtreeLen. clear ();
	for (const DiGraph::Arc* arc : arcs [false])
	{
	  DTNode* child = static_cast <DTNode*> (arc->node [false]);
	  child->setSubtreeLenUp (topological);
	  if (topological)
	  {
	    WeightedMeanVar childSubtreeLen (child->subtreeLen);
	    childSubtreeLen. addValue (1);
  	  subtreeLen. add (childSubtreeLen);
  	}
	  else
  	  subtreeLen. add (child->subtreeLen. getMean () + child->len, child->subtreeLen. weights + child->len);
	}
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



void Steiner::setLca ()
{
  DTNode::setLca ();

  for (const DiGraph::Arc* arc : arcs [false])
  {
    DTNode* child = const_static_cast <DTNode*> (arc->node [false]);
    child->setLca ();  // DFS
    merge (*child);
    static_cast <DTNode*> (getDisjointCluster ()) -> tarjanLca = this;
  }
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
    pathObjNums = child->pathObjNums;
  }
}



void Steiner::makeRoot (Steiner* ancestor2descendant)
{
  ASSERT (ancestor2descendant);
  ASSERT (descendantOf (ancestor2descendant));

  const DTNode* parent_new = static_cast <const DTNode*> (ancestor2descendant->getParent ());
  // Arc-specific data
  const auto len_new         = ancestor2descendant->len;
  const auto attr_new        = ancestor2descendant->attr;
  const auto pathObjNums_new = ancestor2descendant->pathObjNums;
  
  reverseParent (ancestor2descendant, nullptr);

  setParent (const_cast <DTNode*> (parent_new));
  // Arc-specific data
  len         = len_new;
  attr        = attr_new;
  pathObjNums = pathObjNums_new;
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



Leaf::Leaf (DistTree &tree,
  	        Leaf* other,
  	        const string &name_arg)
: DTNode ( tree   // DistTree must be declared
         , const_static_cast <Steiner*> (checkPtr<Leaf> (other) -> getParent ())  // temporary
         , 0
         )  
, name (name_arg)
{
  collapse (other);
}



void Leaf::qc () const
{ 
  if (! qc_on)
    return;
	DTNode::qc ();
	  
  if (graph)
  {
    ASSERT (getParent ());
  	ASSERT (isLeaf ());	
    if (getParent () && discernable != static_cast <const DTNode*> (getParent ()) -> childrenDiscernable ())
    {
      cout << getName () << " " << discernable << " " << getParent () -> getName () 
               << " " << static_cast <const DTNode*> (getParent ()) -> childrenDiscernable () << endl;
      ERROR;
    }
    for (const size_t objNum : pathObjNums)
      ASSERT (getDistTree (). dissims [objNum]. hasLeaf (this));
  }

  ASSERT (! name. empty());
  ASSERT (! isLeft (name, "0x"));
  if (! isNan (len) && ! discernable && len != 0)
  {
    cout << getName () << " " << len << endl;
    ERROR;
  }  
  IMPLY (reprLeaf, reprLeaf == this);
}



void Leaf::setLca ()
{
  DTNode::setLca ();

  for (const size_t objNum : pathObjNums)
  {
    const Leaf* other = getDissimOther (objNum);
    ASSERT (other != this);
    if (! other->graph)
      continue;
    if (! other->tarjanLca)
      continue;
    const DTNode* lca = static_cast <DTNode*> (const_cast <Leaf*> (other) -> getDisjointCluster ()) -> tarjanLca;
    ASSERT (lca);
    const Steiner* st = lca->asSteiner ();
    ASSERT (st);
    const Dissim& dissim = getDistTree (). dissims [objNum];
    ASSERT (! dissim. lca);
    const_cast <Dissim&> (dissim). lca = st;
  }
}



const Leaf* Leaf::getDissimOther (size_t objNum) const
{ 
  const Dissim& dissim = getDistTree (). dissims [objNum];
  ASSERT (dissim. hasLeaf (this));
  return dissim. leaf1 == this ? dissim. leaf2 : dissim. leaf1; 
}



const DTNode* Leaf::getDiscernable () const
{ 
  if (const Leaf* leaf = inDiscernable ())
    return static_cast <const DTNode*> (leaf->getParent ());
  return this;
}



Real Leaf::getRelCriterion (bool strong) const
{ 
  const Real x = absCriterion / getDistTree (). getLeafAbsCriterion (); 
  if (isNan (x))
    return x;
  ASSERT (x >= 0);
  return strong 
           ? x == 0 
             ? NAN 
             : log (x)
           : x;
  
}



void Leaf::collapse (Leaf* other)
{
  ASSERT (this != other);
  ASSERT (len == 0);
//ASSERT (discernable);

  if (! other)
    return;

  ASSERT (graph);
  ASSERT (graph == other->graph);
  
  Vector<Leaf*> indiscernables (1, this);
  {
    const Steiner* oldParent = static_cast <const Steiner*> (getParent ()); 
    ASSERT (oldParent);
    if (! oldParent->childrenDiscernable ())
    {
      indiscernables. clear ();
      indiscernables. reserve (oldParent->arcs [false]. size ());
      for (const DiGraph::Arc* arc : oldParent->arcs [false])
        indiscernables << static_cast <Leaf*> (arc->node [false]);
      ASSERT (indiscernables. size () >= 2);
    }
  }
  ASSERT (! indiscernables. empty ());
    
  if (other->discernable)
  {
    auto st = new Steiner ( const_cast <DistTree&> (getDistTree ())
                          , const_static_cast <Steiner*> (other->getParent ())
                          , other->len
                          );
    other->setParent (st);
    for (Leaf* leaf : indiscernables)
      leaf->setParent (st);
    other->len = 0;
    other->discernable = false;
  }
  else
  { 
    ASSERT (other->len == 0);
    if (getParent () != other->getParent ())
      for (Leaf* leaf : indiscernables)
        leaf->setParent (const_static_cast <Steiner*> (other->getParent ()));
  }

  for (Leaf* leaf : indiscernables)
    leaf->discernable = false;
}




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
//IMPLY (tree. ds. objs. size () == tree. getDissimSize_max (), ok);  
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




// Change

bool Change::apply_ ()
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
  

#if 0
  Subgraph subgraph (const_cast <DistTree&> (tree));
  {
    const Tree::TreeNode* lca = nullptr;
    subgraph. changedLcas = DistTree::getPath (from, to, from->getAncestor (areaRadius_std), lca);
    ASSERT (subgraph. changedLcas. size () >= 2);
    ASSERT (lca);
    if (to == lca)    
    {
      subgraph. changedLcas << lca;
      subgraph. changedLcas_root = static_cast <const DTNode*> (lca->getParent ());
    }
    else
      subgraph. changedLcas_root = static_cast <const DTNode*> (lca);

    // ??
    cout << "root: " << tree. root << endl;
    cout << "from: ";
    from->printAncestors (lca);
    cout << "to: ";
    to->printAncestors (lca);
  }
#endif
  

  // Topology
  inter = new Steiner (const_cast <DistTree&> (tree), arcEnd, 0);
  const_cast <DTNode*> (from) -> setParent (inter); 
  const_cast <DTNode*> (to)   -> setParent (inter);
  
  // DTNode::len
  const_cast <DTNode*> (from) -> len = 0;
  const_cast <DTNode*> (to)   -> len = 0;
  ASSERT (inter->len == 0);


#if 0
  // subgraph
  {
  #if 0
    VectorPtr<Tree::TreeNode>& area     = subgraph. area;
    VectorPtr<Tree::TreeNode>& boundary = subgraph. boundary;
    const Tree::TreeNode* lca_ = nullptr;
    area = DistTree::getPath (from, to, from->getAncestor (areaRadius_std), lca_);
    ASSERT (area. size () >= 2);
    ASSERT (lca_);
    const Steiner* lca = static_cast <const DTNode*> (lca_) -> asSteiner ();
    ASSERT (lca);
    if (   area. front () == from
        || area. front () == to
       )
      area. eraseAt (0);
    if (   area. back () == from
        || area. back () == to
       )
      area. eraseAt (area. size () - 1);
    area << lca;
    area. sort ();
    ASSERT (area. isUniq ());  // Interior area
    boundary. reserve (2 * area. size ());  // PAR
    for (const Tree::TreeNode* node : area)
    {
      const VectorPtr<DiGraph::Node> children (node->getChildren ());
      for (const DiGraph::Node* child : children)
        if (! area. containsFast (static_cast <const Tree::TreeNode*> (child)))
          boundary << static_cast <const Tree::TreeNode*> (child);
    }
    if (const Tree::TreeNode* parent = lca->getParent ())
      boundary << parent;
    boundary. sort ();
    ASSERT (boundary. isUniq ());
    ASSERT (! area. intersectsFast (boundary));
    area << boundary;  // Real area
    // ??
    cout << "boundary:" << endl;
    for (const Tree::TreeNode* node : boundary)
      cout << ' ' << node;
    cout << endl;
    cout << "from = " << from << "  to = " << to << endl;
    
    ASSERT (subgraph. area_root == lca);
  #else
    subgraph. boundary << from << to;
    if (arcEnd)
      subgraph. boundary << arcEnd;
    subgraph. area << subgraph. boundary;
    subgraph. area << inter;
  #endif
    subgraph. finish ();
    subgraph. addSubPaths (from->pathObjNums);
    subgraph. addSubPaths (to  ->pathObjNums);
    subgraph. subPaths. sort ();
    subgraph. subPaths. uniq ();
    subgraph. finishSubPaths ();
    subgraph. qc ();
    if (qc_on)
      for (const SubPath& subPath : subgraph. subPaths)
        ASSERT (   subPath. contains (from)
                || subPath. contains (to)
               );
  }
#endif
  

  // ??
  {
    // tree.{fromAttr_new,toAttr_new,interAttr,target_new}
    CompactBoolAttr1* fromAttr_new = const_cast <CompactBoolAttr1*> (tree. fromAttr_new);
    CompactBoolAttr1* toAttr_new   = const_cast <CompactBoolAttr1*> (tree. toAttr_new);
    CompactBoolAttr1* interAttr    = const_cast <CompactBoolAttr1*> (tree. interAttr);
    RealAttr1* target_new = const_cast <RealAttr1*> (tree. target_new);
    fromAttr_new->setAll (false);
    toAttr_new  ->setAll (false);
    interAttr   ->setAll (false);
    for (Iterator it (tree. dsSample); it ();)  
    {
      const VectorPtr<Tree::TreeNode> path (tree. dissim2path (*it));
      (*target_new) [*it] = (* tree. target) [*it] - tree. path2prediction (path);  
        // redundant <= map paths in the old tree to path in the new tree ??
      if (path. contains (from))
      {
        fromAttr_new->setCompactBool (*it, true);
        const bool toVia    = path. contains (to);
        const bool interVia = path. contains (inter);
        ASSERT (toVia != interVia);
        if (toVia)
          toAttr_new->setCompactBool (*it, true);
        if (interVia)
          interAttr ->setCompactBool (*it, true);
      }
      else
      {
        const bool toUsed = path. contains (to);
        ASSERT (toUsed == path. contains (inter));
        if (toUsed)
        {
          toAttr_new->setCompactBool (*it, true);
          interAttr ->setCompactBool (*it, true);
        }
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
      ASSERT (inter == tree. root);
    }
  
    // tree.prediction
    for (Iterator it (tree. dsSample); it ();)  
      (* const_cast <RealAttr1*> (tree. prediction)) [*it] = (* tree. target) [*it] - lr. getResidual (*it);
  }
    

#if 0
  Dataset ds;
  ds. objs. reserve (subgraph. subPaths. size ());
  auto target = new RealAttr1 ("target", ds);
//for (const SubPath& subPath : subgraph. subPaths)
  FOR (size_t, objNum, subgraph. subPaths. size ())
    ds. appendObj ();
  auto fromAttr  = new ExtBoolAttr1 ("from",  ds);
  auto toAttr    = new ExtBoolAttr1 ("to",    ds);
  auto interAttr = new ExtBoolAttr1 ("inter", ds);
  fromAttr ->setAll (EFALSE);
  toAttr   ->setAll (EFALSE);
  interAttr->setAll (EFALSE);
  FOR (size_t, objNum, subgraph. subPaths. size ())
  {
    const SubPath& subPath = subgraph. subPaths [objNum];
    const Tree::TreeNode* lca_ = nullptr;
  //cout << "SubPath: " << subPath. node1 << ' ' << subPath. node2 << ' ' << subgraph. area_root << endl;  // ??
    VectorPtr<Tree::TreeNode> path (DistTree::getPath (subPath. node1, subPath. node2, subgraph. area_root, lca_));
    if (path. contains (from))
    {
      (*fromAttr) [objNum] = ETRUE;
      const bool toVia    = path. contains (to);
      const bool interVia = path. contains (inter);
      ASSERT (toVia != interVia);
      if (toVia)
        (*toAttr) [objNum] = ETRUE;
      if (interVia)
        (*interAttr) [objNum] = ETRUE;
    }
    else
    {
      const bool toUsed = path. contains (to);
      ASSERT (toUsed);
      ASSERT (toUsed == path. contains (inter));
      if (toUsed)
      {
        (*toAttr)    [objNum] = ETRUE;
        (*interAttr) [objNum] = ETRUE;
      }
    }
    const size_t wholeObjNum = subPath. objNum;
    const_cast <Obj*> (ds. objs [objNum]) -> mult = tree. ds. objs [wholeObjNum] -> mult; 
    (*target) [objNum] = (* tree. target) [wholeObjNum] - subPath. dist_hat_tails;
  }
  
  Space1<NumAttr1> sp (ds, false);  sp. reserve (3);
  sp << fromAttr;  // 0
  if (arcEnd)
    sp << toAttr      // 1
       << interAttr;  // 2
  const Sample sample (ds);
  L2LinearNumPrediction lr (sample, sp, *target);
  lr. solveUnconstrained ();
  if (verbose ())
    lr. qc ();  
  ASSERT (! isNan (lr. absCriterion));

  FOR (size_t, i, lr. beta. size ())
    maximize (lr. beta [i], 0.0);

  ONumber on (cout, 6, false);  // ??
  if (arcEnd)
  {
  #if 0
    ??
    const_cast <DTNode*> (from) -> len = lr. beta [0];
    const_cast <DTNode*> (to)   -> len = lr. beta [1];
    inter                       -> len = lr. beta [2];
  #else
    cout << from ->len << ' ' << lr. beta [0] << endl;
    cout << to   ->len << ' ' << lr. beta [1] << endl;
    cout << inter->len << ' ' << lr. beta [2] << endl;
  #endif
  }
  else
  {
  #if 0
    ??
    const_cast <DTNode*> (from) -> len = lr. beta [0] / 2;
    const_cast <DTNode*> (to)   -> len = from->len;
    inter                       -> len = NAN;
  #else
    cout << from ->len << ' ' << lr. beta [0] / 2 << endl;
    cout << to   ->len << ' ' << lr. beta [0] / 2 << endl;
    cout << inter->len << endl;
  #endif
    ASSERT (inter == tree. root);
  }

#if 0
  ??
  // tree.prediction
  for (Iterator it (tree. dsSample); it ();)  
    (* const_cast <RealAttr1*> (tree. prediction)) [*it] = (* tree. target) [*it] - lr. getResidual (*it);
#endif
#endif


  return true;
}



void Change::restore_ ()
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



void Change::commit_ ()
{
  ASSERT (oldParent);
//ASSERT (arcEnd);
  ASSERT (inter);
  
  inter->addAttr ();  
  if (oldParent->isTransient ())
    const_cast <DistTree&> (tree). delayDeleteRetainArcs (oldParent);
}




// Dissim

Dissim::Dissim (const Leaf* leaf1_arg,
                const Leaf* leaf2_arg)
: leaf1 (leaf1_arg)
, leaf2 (leaf2_arg)
{
  ASSERT (leaf1);
  ASSERT (leaf2);
  ASSERT (leaf1->graph);
  ASSERT (leaf1->graph == leaf2->graph);
  if (leaf1->name > leaf2->name)
    swap (leaf1, leaf2);
  ASSERT (leaf1->name < leaf2->name);
}



void Dissim::qc () const
{
  ASSERT (leaf1);
  ASSERT (leaf2);
  ASSERT (leaf1 != leaf2);
  ASSERT (mult >= 0);
}



string Dissim::getObjName () const
{ 
  return DistTree::getObjName (leaf1->name, leaf2->name); 
}




// SubPath

void SubPath::qc () const
{
  ASSERT (objNum != NO_INDEX);
  
  ASSERT (node1);
  ASSERT (node2);
  ASSERT (node1 != node2);
  ASSERT (node1->graph);
  ASSERT (node1->graph == node2->graph);
  
  ASSERT (dist_hat_tails >= 0);
}




// Subgraph

Subgraph::Subgraph (const DistTree &tree_arg)
: tree (tree_arg)
{ 
  subPaths. reserve (tree. dissims. size ()); 
}



void Subgraph::qc () const 
{
  if (! qc_on)
    return;
    
//IMPLY (changedLcas. empty (), ! changedLcas_root);
    
  // area
  ASSERT (area. searchSorted);
  ASSERT (area. size () >= 2);
  for (const Tree::TreeNode* node : area)
  {
    ASSERT (node);
    ASSERT (node->graph == & tree);
    ASSERT (! static_cast <const DTNode*> (node) -> inDiscernable ());
  }
    
  // boundary
  ASSERT (boundary. searchSorted);
  ASSERT (area. containsFastAll (boundary));

  // area_root
  ASSERT (area_root);
  ASSERT (area_root->graph == & tree);
  ASSERT (area. containsFast (area_root));
  ASSERT (! area. containsFast (area_root->getParent ()));
  
  // area_underRoot
  if (area_underRoot)
  {
    ASSERT (area_underRoot->getParent () == area_root);
    ASSERT (area. containsFast (area_underRoot));
  }
  ASSERT (boundary. containsFast (area_root) == (bool) area_underRoot);

  // subPaths
  ASSERT (! subPaths. empty ());
  {
    Vector<size_t> objNums;  objNums. reserve (subPaths. size ());
    for (const SubPath& subPath : subPaths)
    {
      subPath. qc ();
      ASSERT (boundary. containsFast (subPath. node1));
      ASSERT (boundary. containsFast (subPath. node2));
      ASSERT (tree. dissims [subPath. objNum]. valid ());
      objNums << subPath. objNum;
    }
    objNums. sort ();
    ASSERT (objNums. isUniq ());
  }
}  



void Subgraph::finish ()
{
  area. sort ();
  boundary. sort ();
//changedLcas. sort ();
  
  // area_root
  ASSERT (! area_root);
  for (const Tree::TreeNode* node : area)
  {
    ASSERT (node);
    if (! area. containsFast (node->getParent ()))
    {
      ASSERT (! area_root);
      area_root = static_cast <const DTNode*> (node);
    }
  }
  ASSERT (area_root);
  
  // area_underRoot
  ASSERT (! area_underRoot);
  if (boundary. containsFast (area_root))
  {
    const VectorPtr<DiGraph::Node> children (area_root->getChildren ());
    for (const DiGraph::Node* child_ : children)
    {
      const DTNode* child = static_cast <const DTNode*> (child_);
      if (area. containsFast (child))
      {
        ASSERT (! area_underRoot);
        area_underRoot = child;
      }
    }
    ASSERT (area_underRoot);
  }
}



void Subgraph::addSubPaths (const Vector<size_t> &objNums)
{ 
  for (const size_t objNum : objNums)
    if (tree. dissims [objNum]. valid ())
      subPaths << SubPath (objNum);
}



void Subgraph::area2subPaths ()
// Parameter sparse: Use O(boundary.size()) dissimilarities, use reprLeaf's like in sparsing ??
{
  ASSERT (subPaths. empty ());
  for (const Tree::TreeNode* node : boundary)
  {
    const DTNode* dtNode = static_cast <const DTNode*> (node);
    const Vector<size_t>& pathObjNums = dtNode == area_root /* && area_underRoot */
                                          ? area_underRoot->pathObjNums
                                          : dtNode        ->pathObjNums;
    addSubPaths (pathObjNums);
  }
}



void Subgraph::finishSubPaths ()
{
  ASSERT (isNan (subPathsAbsCriterion));
  
  subPathsAbsCriterion = 0;
  subPaths. sort ();
  subPaths. uniq ();

  // SubPath::{node1,node2}
  for (const Tree::TreeNode* node : boundary)
  {
    const DTNode* dtNode = static_cast <const DTNode*> (node);
    const Vector<size_t>& pathObjNums = dtNode == area_root /* && area_underRoot */
                                          ? area_underRoot->pathObjNums
                                          : dtNode        ->pathObjNums;
    for (const size_t objNum : pathObjNums)
    {
      if (! tree. dissims [objNum]. valid ())
        continue;
      const SubPath subPath_ (objNum);
      const size_t index = subPaths. binSearch (subPath_);
      if (index == NO_INDEX)  // impossible ??
        continue;
      SubPath& subPath = subPaths [index];
      ASSERT (subPath. objNum == objNum);
      if (! subPath. node1)
        subPath. node1 = dtNode;
      else if (! subPath. node2)
        subPath. node2 = dtNode;
      else
        { ERROR; }
    }
  }


  // subPathsAbsCriterion, SubPath::dist_hat_tails
  for (SubPath& subPath : subPaths)
  {
    const size_t wholeObjNum = subPath. objNum;

    ASSERT (! isNan ((* tree. prediction) [wholeObjNum]));
    subPathsAbsCriterion += tree. getAbsCriterion (wholeObjNum);

    ASSERT (subPath. node1);
    ASSERT (subPath. node2);
    const Dissim& dissim = tree. dissims [wholeObjNum];
    ASSERT (dissim. lca);
    const Tree::TreeNode* lca = nullptr;
    VectorPtr<Tree::TreeNode> path (DistTree::getPath (subPath. node1, subPath. node2, dissim. lca, lca));
    ASSERT (path. size () >= 2);

  #if 0
    const Real dist_hat_whole = DistTree::path2prediction (path1);
    ASSERT (eqReal (dist_hat_whole, (* tree. prediction) [wholeObjNum], 1e-4));  // PAR
    for (const Tree::TreeNode* node : path)
      ASSERT (path1. contains (node));
  #endif

    const Real dist_hat_sub = DistTree::path2prediction (path);
    ASSERT (isNan (subPath. dist_hat_tails));
    subPath. dist_hat_tails = max (0.0, (* tree. prediction) [wholeObjNum] - dist_hat_sub);
    ASSERT (! isNan (subPath. dist_hat_tails));

    subPath. qc ();
  }
}



void Subgraph::subPaths2tree ()
{
  chron_subgraph2tree. start ();  

  DistTree& tree_ = const_cast <DistTree&> (tree);
  
  const DTNode* area_root_ = area_root->graph ? area_root: nullptr;

  tree_. absCriterion -= subPathsAbsCriterion;
  for (const SubPath& subPath : subPaths)
  {
    subPath. qc ();
    ASSERT (subPath. node1->graph == & tree);
    ASSERT (subPath. node2->graph == & tree);
    const size_t objNum = subPath. objNum;
    const Tree::TreeNode* lca_ = nullptr;
    const VectorPtr<Tree::TreeNode> path (DistTree::getPath (subPath. node1, subPath. node2, area_root_, lca_));
    ASSERT (lca_);
    if (! viaRoot (subPath))
    {
      const Steiner* lca = static_cast <const DTNode*> (lca_) -> asSteiner ();
      ASSERT (lca);
      tree_. dissims [objNum]. lca = lca;
    }
    for (const Tree::TreeNode* node : path)
      if (! boundary. containsFast (node))
      {
        const Steiner* st = static_cast <const DTNode*> (node) -> asSteiner ();
        ASSERT (st);
        const_cast <Steiner*> (st) -> pathObjNums << objNum;
      }
    (* const_cast <RealAttr1*> (tree. prediction)) [objNum] = subPath. dist_hat_tails + DistTree::path2prediction (path);  
    ASSERT (! isNan ((* tree. prediction) [objNum]));
    tree_. absCriterion += tree. getAbsCriterion (objNum);
  }
  ASSERT (DM_sp::finite (tree. absCriterion));

  chron_subgraph2tree. stop ();
  
  if (verbose ())
    tree. qcPaths (); 
}




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

    nodeAttrExist = true;
    dissimDs2ds (sparse);    	
    topology2attrs ();
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

  setGlobalLen ();  // --> after dissim2Ds(), use dissims ??

  nodeAttrExist = true;
  dissimDs2ds (false);  
  topology2attrs ();
}



DistTree::DistTree (const string &dissimFName,
	                  const string &attrName,
	                  bool sparse)
: dsSample (ds)
{
  loadDissimDs (dissimFName, attrName);

  // Initial tree topology: star topology 
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
  
  nodeAttrExist = true;
  dissimDs2ds (sparse);    	 
  topology2attrs ();
}



DistTree::DistTree (const string &dataDirName,
                    bool loadNewLeaves,
 	                  bool loadDissim)
: dsSample (ds)
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
  if (loadNewLeaves)
  {
    LineInput f (dataDirName + "leaf", 10 * 1024, 1);  // PAR
    string leafName, anchorName;
    Real leafLen, arcLen;
    while (f. nextLine ())
    {
      istringstream iss (f. line);
      iss >> leafName >> anchorName >> leafLen >> arcLen;
      ASSERT (leafLen >= 0);
      ASSERT (arcLen >= 0);
    //ASSERT (iss. eof ());  // Extra fields
      const DTNode* anchor = lcaName2node (anchorName);
      ASSERT (anchor);
      // Use absCriterion_leaf in f to test if leafName is an outlier
      //   if leafName is an outlier then do not add leafName to *this, but add leafName to outliers ??
      // Attach a leaf
      Leaf* leaf = nullptr;
      if (const Leaf* anchorLeaf = anchor->asLeaf ())
        if (leafLen == 0 && arcLen == 0)
          leaf = new Leaf (*this, const_cast <Leaf*> (anchorLeaf), leafName);
      if (! leaf)
      {
        while (anchor != root && greaterReal (arcLen, anchor->len))
        {
          arcLen -= anchor->len;
          anchor = static_cast <const DTNode*> (anchor->getParent ());
        }
        if (anchor == root)  
          leaf = new Leaf (*this, const_static_cast <Steiner*> (root), arcLen + leafLen, leafName);
        else
        {
          auto st = new Steiner ( *this
                                , const_static_cast <Steiner*> (anchor->getParent ())
                                , max (0.0, anchor->len - arcLen)
                                );
          ASSERT (st);
          const_cast <DTNode*> (anchor) -> setParent (st);
          const_cast <DTNode*> (anchor) -> len = arcLen;
          leaf = new Leaf (*this, st, leafLen, leafName);
        }
      }
      ASSERT (leaf);
      name2leaf [leaf->name] = leaf;
      newLeaves << leaf;
    }
  }


  if (loadDissim)
  {
    StringVector outliers;
    {
      FileItemGenerator fig (0, true, dataDirName + "outlier");  
  	  string item;
  	  while (fig. next (item))
  	    outliers << item;
    }
    outliers. sort ();
    loadDissimPrepare (name2leaf. size () * (size_t) log ((Real) name2leaf. size () * 10), dissimDecimals);  // PAR
    {
      const string fName (dataDirName + "dissim");
      cout << "Loading " << fName << " ..." << endl;
      LineInput f (fName, 10 * 1024 * 1024, 100000);  // PAR
      while (f. nextLine ())
      {
        const string name1 = findSplit (f. line, '\t');
        const string name2 = findSplit (f. line, '\t');
        if (name2. empty ())
          throw runtime_error ("Empty name2");
        if (outliers. containsFast (name1))
          continue;
        if (outliers. containsFast (name2))
          continue;
        const Leaf* leaf1 = findPtr (name2leaf, name1);
        if (! leaf1)
          continue;
        //throw runtime_error ("Tree has no object " + name1);
        const_cast <Leaf*> (leaf1) -> paths ++;
        const Leaf* leaf2 = findPtr (name2leaf, name2);
        if (! leaf2)
          continue;
        //throw runtime_error ("Tree has no object " + name2);
        const_cast <Leaf*> (leaf2) -> paths ++;
        const Real dissim = str2real (f. line);
        if (dissim == 0 && ! leaf1->getCollapsed (leaf2))  // Only for new Leaf's
          const_cast <Leaf*> (leaf1) -> collapse (const_cast <Leaf*> (leaf2));
        EXEC_ASSERT (addDissim (name1, name2, dissim)); 
      }
    }
    // qc
    for (const DiGraph::Node* node : nodes)
      if (const Leaf* leaf = static_cast <const DTNode*> (node) -> asLeaf ())
      {
        if (! leaf->paths)
          throw runtime_error ("No dissimilarities for object " + leaf->name);
        const_cast <Leaf*> (leaf) -> paths = 0;
      }
    // qc: pairs of <leaf1,leaf2> must be unique in ds ??
  } 
  
  
  cleanTopology ();
   
  
  if (loadDissim)
  {
    loadDissimFinish ();        
    qc ();     

    {  
      const Chronometer_OnePass cop ("Optimizing new leaves");  
      cout << "Optimizing new leaves ..." << endl;
      Progress prog ((uint) newLeaves. size ());
      for (const Leaf* leaf : newLeaves)
      {
        prog (real2str (absCriterion, 6));
        Unverbose unv;
        optimizeSubgraph (static_cast <const Steiner*> (leaf->getAncestor (areaRadius_std)));  // PAR ??
      }
    }
  }  


  ASSERT (! dissimDs. get ());
  ASSERT (! dissimAttr);
}



DistTree::DistTree (const string &newickFName)
: dsSample (ds)
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
}



DistTree::DistTree (const DTNode* center,
                    uint areaRadius,
                    Subgraph& subgraph,
                    Node2Node &newLeaves2boundary)
: dsSample (ds)
{
  ASSERT (center);
  ASSERT (center->graph != this);
  ASSERT (! center->inDiscernable ());
  ASSERT (areaRadius >= 1);  
  ASSERT (newLeaves2boundary. empty ());
  
  const DistTree& wholeTree = center->getDistTree ();


  // subgraph
  ASSERT (subgraph. empty ());
  ASSERT (& wholeTree == & subgraph. tree);
  VectorPtr<TreeNode>& area     = subgraph. area;
  VectorPtr<TreeNode>& boundary = subgraph. boundary;
  const DTNode* &area_root      = subgraph. area_root;

  // area, boundary
  center->getArea (areaRadius, area, boundary);  
  ASSERT (! area. empty ());
  if (area. size () == 1)
    throw runtime_error ("Singleton tree");
  // area: remove !discernable
  for (Iter<VectorPtr<TreeNode>> iter (area); iter. next (); )
    if (static_cast <const DTNode*> (*iter) -> inDiscernable ())
      iter. erase ();
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
    if (qc_on)
      for (const TreeNode* node : boundary)
        ASSERT (! newBoundary. contains (node));
    for (const TreeNode* node : newBoundary)
      boundary << node;    
  }
  
  subgraph. finish ();
  ASSERT (area. containsFast (center));
  subgraph. area2subPaths ();    
  subgraph. finishSubPaths ();
  subgraph. qc ();


  // nodes
  Node2Node old2new;  // 1-1
  size_t leafNum = 0;
  for (const TreeNode* node_old : area)
  {
    ASSERT (node_old);
    DTNode* node_new = nullptr;
    {
      const Real len = static_cast <const DTNode*> (node_old) -> len;
      VectorPtr<DiGraph::Node> children (node_old->getChildren ());
      children. sort ();
      if (children. intersectsFast (area))
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
  ASSERT (nodes. size () == area. size ());
  ASSERT (nodes. size () == old2new. size ());
  
  borrowArcs (old2new, true);  // Arc's are not parallel
  setRoot ();
  ASSERT (static_cast <const DTNode*> (root) -> asSteiner ());
  ASSERT (findPtr (old2new, area_root) == root);
  
  // subgraph.area_root --> Leaf
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
      ASSERT (static_cast <const DTNode*> (findPtr (old2new, area_root)) -> asLeaf ());
    }
    else
    { 
      ASSERT (static_cast <const DTNode*> (findPtr (old2new, area_root)) -> asSteiner ());
      ASSERT (area_root == wholeTree. root);
      ASSERT (children. size () == area_root->getChildren (). size ()); 
    }
  }
  ASSERT (isNan (static_cast <const DTNode*> (root) -> len));  

  
  setName2leaf ();
  ASSERT (boundary. size () == name2leaf. size ());

  // newLeaves2boundary
  for (const auto it : old2new)
    if (static_cast <const DTNode*> (it. second) -> asLeaf ())
      newLeaves2boundary [it. second] = it. first;
  ASSERT (newLeaves2boundary. size () == boundary. size ());
    

  // ds.objs[]->mult, *target, dissims: init
  // For some leaf pairs the dissimilarity may be missing
  ds. objs. reserve (getDissimSize_max ());
  target = new RealAttr1 ("Target", ds, wholeTree. target->decimals); 
  dissims. reserve (ds. objs. capacity ());
  for (const auto it2 : name2leaf)
    for (const auto it1 : name2leaf)
    {
      if (it1. first == it2. first)
        break;
      ASSERT (it1. first < it2. first);
      const size_t objNum = ds. appendObj (getObjName (it1. first, it2. first));  
      const_cast <Obj*> (ds. objs [objNum]) -> mult = 0;
      (* const_cast <RealAttr1*> (target)) [objNum] = 0;
      EXEC_ASSERT (leaves2dissims (const_cast <Leaf*> (it1. second), const_cast <Leaf*> (it2. second)) == objNum);
    }
  ds. setName2objNum ();


  // ds.objs[]->mult, *target: sum
  for (const SubPath& subPath : subgraph. subPaths)
  {
    const size_t wholeObjNum = subPath. objNum;
    
    const Real dist_whole = (* wholeTree. target) [wholeObjNum] - subPath. dist_hat_tails;
    ASSERT (! isNan (dist_whole));
      // May be < 0
    const Real mult_whole = wholeTree. ds. objs [wholeObjNum] -> mult; 
    ASSERT (mult_whole > 0);
    
    // objNum  
    const DTNode* node1 = static_cast <const DTNode*> (findPtr (old2new, subPath. node1));
    const DTNode* node2 = static_cast <const DTNode*> (findPtr (old2new, subPath. node2));
    ASSERT (node1);
    ASSERT (node2);
    ASSERT (node1 != node2);
    const Leaf* leaf1 = node1->asLeaf ();
    const Leaf* leaf2 = node2->asLeaf ();
    ASSERT (leaf1);
    ASSERT (leaf2);
    const size_t objNum = ds. getName2objNum (getObjName (leaf1->name, leaf2->name));
    ASSERT (objNum != NO_INDEX);
    
    const_cast <Obj*> (ds. objs [objNum]) -> mult += mult_whole;
    (* const_cast <RealAttr1*> (target)) [objNum] += mult_whole * dist_whole;
  }


  // *target, dissim2_sum
  FOR (size_t, objNum, ds. objs. size ())
  {
    const Real mult = ds. objs [objNum] -> mult;
    if (mult == 0)
      const_cast <RealAttr1*> (target) -> setMissing (objNum);
    else
    {
      (* const_cast <RealAttr1*> (target)) [objNum] /= mult;
      dissim2_sum += mult * sqr ((*target) [objNum]);  
    }
    dissims [objNum]. mult = mult;
  }
    

  nodeAttrExist = true;
  loadDissimFinish ();
//ASSERT (setDiscernable () == 0);
  topology2attrs ();
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
    leaf->relCriterion = leafError;
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
      if (const Leaf* leaf1 = findPtr (name2leaf, dissimDs->objs [row] -> name))
        FOR (size_t, col, row)  // dissimAttr is symmetric
          if (dissimAttr->get (row, col) < INF)
            if (const Leaf* leaf2 = findPtr (name2leaf, dissimDs->objs [col] -> name))
              const_cast <Leaf*> (leaf1) -> merge (* const_cast <Leaf*> (leaf2));
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
   	  {
   	    leaf->discernable = true;
   	    leaf->DisjointCluster::init ();
   	  }
   	}
    FOR (size_t, row, dissimDs->objs. size ())
      if (const Leaf* leaf1 = findPtr (name2leaf, dissimDs->objs [row] -> name))
        FOR (size_t, col, row)  // dissimAttr is symmetric
          if (! positive (dissimAttr->get (row, col)))
            if (const Leaf* leaf2 = findPtr (name2leaf, dissimDs->objs [col] -> name))
              const_cast <Leaf*> (leaf1) -> merge (* const_cast <Leaf*> (leaf2));
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
    cleanTopology ();

  
  return n;
}



void DistTree::cleanTopology ()
{
 	Vector<DiGraph::Node*> nodeVec;  nodeVec. reserve (nodes. size ());

 	// delete isLeaf() && Steiner
 	insertAll (nodeVec, nodes);
 	for (DiGraph::Node* node : nodeVec)  
	  if (const Steiner* s = static_cast <DTNode*> (node) -> asSteiner ())
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
 		  ASSERT (s->arcs [false]. size () == 1);
 		  ASSERT (s->childrenDiscernable ());
 			delayDeleteRetainArcs (const_cast <Steiner*> (s));
 	  }
 	}

  toDelete. deleteData ();
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
    if (const Leaf* leaf1 = findPtr (name2leaf, dissimDs->objs [row] -> name))
      FOR (size_t, col, row)  // dissimAttr is symmetric
      {
        if (const Leaf* leaf2 = findPtr (name2leaf, dissimDs->objs [col] -> name))
        {
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

  clearSubtreeLen ();

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
  ASSERT (isStar ());
	ASSERT (dissimDs. get ());
	ASSERT (dissimAttr);
  ASSERT (dissimDs->objs. size () >= 2);
	ASSERT (! optimizable ());
    
  cout << "Neighbor joining ..." << endl;
  
  // DTNode::len: sum of dissimilarities from other objects (dissim_sum)

  size_t n = 0;
	for (const DiGraph::Arc* arc : root->arcs [false])
	{
	  const_static_cast <DTNode*> (arc->node [false]) -> len = 0;
	  n++;
	}
    
  Vector<Neighbors> neighborsVec;  neighborsVec. reserve (getDissimSize_max ());
  FOR (size_t, row, dissimDs->objs. size ())
  {
    const Leaf* leaf1 = findPtr (name2leaf, dissimDs->objs [row] -> name);
    ASSERT (leaf1);
    FOR (size_t, col, row)  // dissimAttr is symmetric
    {
      const Leaf* leaf2 = findPtr (name2leaf, dissimDs->objs [col] -> name);
      ASSERT (leaf2);
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
    neighborsVec. sort (Neighbors::strictlyLess);
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
        for (const bool first : {false, true})
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
        for (const bool first : {false, true})
        {
          bool found = false;
          for (const bool best_first : {false, true})
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
  for (const bool first : {false, true})  
  {
    DTNode* node = const_cast <DTNode*> (neighbors. nodes [first]);
    ASSERT (node->getParent () == root);
    node->len = neighbors. dissim / 2;  
  }

  finishChanges ();
  
  reroot (true);  // Reduces ds.objs.size() if ds is sparse
}



void DistTree::dissimDs2ds (bool sparse)
{
	ASSERT (dissimDs. get ());
	ASSERT (dissimAttr);
	ASSERT (! optimizable ());
	ASSERT (! ds. name2objNumSet ());
	ASSERT (detachedLeaves. empty ());


  const StringVector objNames (dissimDs->getObjNames ());
  restrictLeaves (objNames, true);
  if (qc_on)
  {
    const Set<string> leafNames (name2leaf);
    ASSERT (objNames. containsFastAll (leafNames));
  }
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

  // ds.objs, *target, dissims, dissim2_sum
  loadDissimPrepare (getDissimSize_max (), dissimAttr->decimals);
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
  dissims. reserve (ds. objs. capacity ());
}
  


size_t DistTree::leaves2dissims (Leaf* leaf1,
                                 Leaf* leaf2)
{ 
  const size_t objNum = dissims. size ();
  dissims << Dissim (leaf1, leaf2);
  leaf1->pathObjNums << objNum;
  leaf2->pathObjNums << objNum;      
  
  return objNum;
}



bool DistTree::addDissim (const string &name1,
                          const string &name2,
                          Real dissim)
{
  ASSERT (detachedLeaves. empty ());

  if (   isNan (dissim)  // prediction must be large ??
    //|| ! DM_sp::finite (dissim)  // For selectPairs()
     )
  {
    cout << name1 << " - " << name2 << ": " << dissim << endl;
    return false;  
  }
  
  ASSERT (dissim >= 0);
    
  const Real mult = dissim2mult (dissim);  
  if (mult)
    dissim2_sum += mult * sqr (dissim);  // max (0.0, dissim);
  
  const size_t objNum = ds. appendObj (getObjName (name1, name2));
  
  const_cast <Obj*> (ds. objs [objNum]) -> mult = mult;
  (* const_cast <RealAttr1*> (target)) [objNum] = dissim;

  const Leaf* leaf1 = findPtr (name2leaf, name1);
  const Leaf* leaf2 = findPtr (name2leaf, name2);
  ASSERT (leaf1);
  ASSERT (leaf2);
  IMPLY (dissim == 0, leaf1->getCollapsed (leaf2));
  EXEC_ASSERT (leaves2dissims (const_cast <Leaf*> (leaf1), const_cast <Leaf*> (leaf2)) == objNum);  
  
  dissims [objNum]. mult = mult;
  
  return true;
}



void DistTree::loadDissimFinish ()
{
	ASSERT (optimizable ());
  ASSERT (! prediction);
  ASSERT (! fromAttr_new);
  ASSERT (! toAttr_new);
  ASSERT (! interAttr);
  ASSERT (! target_new);
  ASSERT (! prediction_old);
  

  dsSample = Sample (ds);
  absCriterion_delta = dsSample. mult_sum * 1e-5;  // PAR

  // ds.attrs, DTNode::attr
  if (nodeAttrExist)
  {
    if (verbose ())
      cout << "Arcs -> data attributes ..." << endl;
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

  target_new     = new RealAttr1 ("target_new",    ds, target->decimals);
  prediction_old = new RealAttr1 ("PredictionOld", ds, target->decimals); 

  setLca ();  
  for (Iterator it (dsSample); it ();)  
  {
    ASSERT (dissims [*it]. valid ());
    const VectorPtr<TreeNode> path (dissimLca2path (*it));
    (* const_cast <RealAttr1*> (prediction)) [*it] = path2prediction (path);  
    for (const TreeNode* node : path)
      if (const Steiner* st = static_cast <const DTNode*> (node) -> asSteiner ())
        const_cast <Steiner*> (st) -> pathObjNums << *it;
  }
  setAbsCriterion (); 
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
  for (const auto it : name2leaf)
  {
    ASSERT (it. second);
    ASSERT (it. second->graph);
  }
  
  ASSERT ((bool) dissimDs. get () == (bool) dissimAttr);
  if (dissimAttr)
  {
    dissimDs->qc ();
    ASSERT (& dissimAttr->ds == dissimDs. get ());
  }

  ds. qc ();
	ASSERT (dsSample. ds == & ds);
	ASSERT (dsSample. size () == ds. objs. size ());
	
  ASSERT (dissims. size () == ds. objs. size ());
  
  IMPLY (nodeAttrExist, optimizable ());


  size_t leafDissims = 0;
	if (optimizable ())
	{
   	ASSERT (positive (dsSample. mult_sum));
    ASSERT (dsSample. nEffective <= getDissimSize_max ());
   	  
    ASSERT (target);
    ASSERT (prediction); 	
   	ASSERT (& target    ->ds == & ds);
   	ASSERT (& prediction->ds == & ds);
   	ASSERT (! target->existsMissing ());
   	
    Set<pair<const Leaf*,const Leaf*> > leafSet;
    for (Iterator it (dsSample); it ();)  
    {
      const Leaf* leaf1 = dissims [*it]. leaf1;
      const Leaf* leaf2 = dissims [*it]. leaf2;
      ASSERT (leaf1);
      ASSERT (leaf2);
      ASSERT (leaf1->graph);
      ASSERT (leaf2->graph);
      ASSERT (leaf1->name < leaf2->name);
      leafSet. checkUnique (pair<const Leaf*,const Leaf*> (leaf1, leaf2));
    //ASSERT ((*target) [*it] >= 0);
      IMPLY (ds. objs [*it] -> mult && (*target) [*it] == 0, 
             ! leaf1->discernable && ! leaf2->discernable && leaf1->getParent () == leaf2->getParent ()
            );
      if (! isNan ((*prediction) [*it]))
      {
        ASSERT ((*prediction) [*it] >= 0);
      /*
        if (topologyOptimized && nullReal ((*prediction) [*it]) && ! nullReal ((*target) [*it]))
        {
          cout << leaf1->name << " " << leaf2->name << endl;
          ERROR;
        }
      */
      }
      Unverbose unv;
      if (verbose ())
      {
        const VectorPtr<TreeNode> path (dissim2path (*it));
        for (const DiGraph::Node* node : nodes)
        {
          const DTNode* dtNode = static_cast <const DTNode*> (node);
          ASSERT ((*dtNode->attr) [*it] == path. contains (dtNode));
        }
      }
    }
    
    for (const Dissim &dissim : dissims)
      dissim. qc ();
  
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
   	
   	   	
   	// ds.attrs
    Set<const Attr*> nodeAttrs;
   	for (const DiGraph::Node* node : nodes)
   	{
   	  const DTNode* dtNode = static_cast <const DTNode*> (node);
   	  ASSERT (nodeAttrExist == (bool) dtNode->attr);
   		if (dtNode->attr)
   		{
   		  nodeAttrs. checkUnique (dtNode->attr);
   		  if (dtNode != root && ! positive (dtNode->attr->getProb (dsSample)))
   		  {
   		    cout << dtNode->getName () << endl;
   		    cout << endl;
   		    FOR (size_t, objNum, ds. objs. size ())
   		      cout << ds. objs [objNum] -> name << '\t' << ds. objs [objNum] -> mult << '\t' << (*target) [objNum] << endl;
   		    ERROR;
   		  }
   		  ASSERT (& dtNode->attr->ds == & ds);
   		}
    //ASSERT ((dtNode == root) == dtNode->pathObjNums. empty ());  ??
   		if (const Leaf* leaf = dtNode->asLeaf ())
   		  leafDissims += leaf->pathObjNums. size ();
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

    ASSERT (absCriterion_delta > 0);
  }
  

//ASSERT (ds. objs. size () <= (leaves * (leaves - 1)) / 2);
  const size_t discernables = getDiscernables (). size ();
  ASSERT (discernables <= leaves);

  for (const Leaf* leaf : detachedLeaves)
  {
    ASSERT (leaf);
    ASSERT (! leaf->graph);
    leaf->qc ();
	  leafDissims += leaf->pathObjNums. size ();
  }

 	ASSERT (leafDissims == 2 * ds. objs. size ());
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
  os << "# Leaves: " << name2leaf. size () << endl;
  os << "# Discernable leaves: " << getDiscernables (). size () << endl;
  os << "# Nodes: " << nodes. size () << endl;
  if (! optimizable ())
    return;
  os << "# Dissimilarities: " << ds. objs. size () << " (" << (Real) ds. objs. size () / (Real) getDissimSize_max () * 100 << " %)" << endl; 
  os << "Ave. dissimilarity = " << getDissim_ave () << endl;
  if (nodeAttrExist)
    os << "Arc length is optimizable globally" << endl;
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



void DistTree::qcPaths () const
{
  if (! qc_on)
    return;

  // slow ??

  size_t n = 0;
  FOR (size_t, objNum, dissims. size ())
    if (dissims [objNum]. valid ())
    {
      const VectorPtr<TreeNode> path (dissimLca2path (objNum));
      for (const Tree::TreeNode* node : path)
      {
        ASSERT (static_cast <const DTNode*> (node) -> pathObjNums. contains (objNum));
        n++;
      }
    }

  size_t m = 0;
  for (const DiGraph::Node* node : nodes)
    for (const size_t objNum : static_cast <const DTNode*> (node) -> pathObjNums)
      if (dissims [objNum]. valid ())
        m++;
        
  ASSERT (n == m);
}



void DistTree::resetAttrs ()
{ 
  ASSERT (nodeAttrExist);
  
  for (DiGraph::Node* node : nodes)
 	  const_cast <CompactBoolAttr1*> (static_cast <DTNode*> (node) -> attr) -> setAll (false);
}



void DistTree::setLca ()
{
  for (DiGraph::Node* node : nodes)
  {
    DTNode* dtNode = static_cast <DTNode*> (node);
    dtNode->tarjanLca = nullptr;
    dtNode->DisjointCluster::init ();
  }
  for (Dissim& dissim : dissims)
    dissim. lca = nullptr;
  const_static_cast <DTNode*> (root) -> setLca ();
}



VectorPtr<Tree::TreeNode> DistTree::dissimLca2path (size_t objNum) const  
{ 
  const Dissim& dissim = dissims [objNum];
  ASSERT (dissim. lca);
  const TreeNode* lca = nullptr;
  const VectorPtr<TreeNode> path (getPath (dissim. leaf1, dissim. leaf2, dissim. lca, lca));
  ASSERT (lca == dissim. lca);
  
  return path;
}



void DistTree::topology2attrs ()
{
  ASSERT (nodeAttrExist);

  if (verbose ())
    cout << "Data values ..." << endl;

  Progress prog ((uint) ds. objs. size (), 10000);  // PAR
  setLca ();
  for (Iterator it (dsSample); it ();)  
  {
    prog ();
    const VectorPtr<TreeNode> path (dissimLca2path (*it));
    for (const TreeNode* node : path)
    {
      const DTNode* dtNode = static_cast <const DTNode*> (node);
      CompactBoolAttr1* attr = const_cast <CompactBoolAttr1*> (dtNode->attr);
      ASSERT (attr);
      attr->setCompactBool (*it, true);
    }
  }
}



void DistTree::removeTopologyAttrs ()
{
 	if (! nodeAttrExist)
 	  return;

  for (DiGraph::Node* node : nodes)
  {
 	  const CompactBoolAttr1* attr = static_cast <DTNode*> (node) -> attr;
 	  delete attr;
 	  static_cast <DTNode*> (node) -> attr = nullptr;
 	}
 	
 	nodeAttrExist = false;
}



void DistTree::clearSubtreeLen ()
{
  for (DiGraph::Node* node : nodes)
    static_cast <DTNode*> (node) -> subtreeLen. clear ();
}



void DistTree::setPrediction () 
{
  for (Iterator it (dsSample); it ();)  
  {
    const VectorPtr<TreeNode> path (dissim2path (*it));
    (* const_cast <RealAttr1*> (prediction)) [*it] = path2prediction (path);  
  }
}



Real DistTree::getAbsCriterion (size_t objNum) const
{
  ASSERT (optimizable ());

  const Real d = (*target) [objNum];
//ASSERT (d >= 0);
  const Real dHat = (*prediction) [objNum];
  ASSERT (! isNan (dHat));
  return ds. objs [objNum] -> mult * sqr (dHat - d);
}



Real DistTree::getAbsCriterion () const
{
  Real absCriterion_ = 0;
  for (Iterator it (dsSample); it ();)  
    absCriterion_ += getAbsCriterion (*it);  // it. mult * sqr (dHat - d);
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
  for (Iterator it (dsSample); it ();)  
  {
    const Real d = (*target) [*it];
    const Real dHat = (*prediction) [*it];
    ASSERT (! isNan (dHat));
    const bool half2 = d > dissim_ave;
    pairs++;
    size_half [half2] ++;
    absCriterion_half [half2] += it. mult * sqr (dHat - d);
    dissim2_half      [half2] += it. mult * sqr (d);
  }
    
  for (const bool half2 : {false, true})
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
      const_cast <Leaf*> (leaf) -> absCriterion = 0;

  for (Iterator it (dsSample); it ();)  
  {
    const Real d = (*target) [*it];
  //ASSERT (d >= 0);
    const Real dHat = (*prediction) [*it];
    ASSERT (! isNan (dHat));
    const Real residual = sqr (dHat - d);      
    ASSERT (DM_sp::finite (residual));
    const Dissim& dissim = dissims [*it];
    const_cast <Leaf*> (dissim. leaf1) -> absCriterion += residual * it. mult;
    const_cast <Leaf*> (dissim. leaf2) -> absCriterion += residual * it. mult;
  }
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
  for (Iterator it (dsSample); it ();)  
  {
    const Real dissim = (*target) [*it];    
    const VectorPtr<TreeNode> path (dissim2path (*it));
    for (const TreeNode* node1 : path)
    {
      DTNode* dtNode1 = const_static_cast <DTNode*> (node1);
      if (dtNode1->index == NO_INDEX)
        continue;
      dtNode1->dissimSum += dissim * it. mult;
      for (const TreeNode* node2 : path)
      {
        const DTNode* dtNode2 = static_cast <const DTNode*> (node2);
        if (dtNode2->index != NO_INDEX)
          matr. putInc ( false
                       , dtNode1->index
                       , dtNode2->index
                       , it. mult
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
// Requires: complete *dissimAttr, # leaves > 2
// Invokes: setLeaves()
// Time: O(p log(n) + n)
{
	ASSERT (optimizable ());


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
  for (Iterator it (dsSample); it ();)  
  {
    const Real dissim = (*target) [*it];    
    const VectorPtr<TreeNode> path (dissim2path (*it));
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
  ASSERT (nodeAttrExist);
  
  if (verbose (1))
    cout << "Optimizing arc lengths at each arc ..." << endl;
    
#ifndef NDENUG
  const Real absCriterion_old = absCriterion;
#endif
  
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
  dtNodes. sort (DTNode_len_strictlyLess);

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
  ASSERT (leReal (absCriterion, absCriterion_old));
  
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
  ASSERT (nodeAttrExist);
  
  if (verbose (1))
    cout << "Optimizing arc lengths at each node ..." << endl;

#ifndef NDENUG
  const Real absCriterion_old1 = absCriterion;
#endif

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
  dtNodes. sort (NodeStar::strictlyLess);  


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
    for (Iterator it (dsSample); it ();)  // optimize ??
    {
      const VectorPtr<TreeNode> path (dissim2path (*it));
      (* const_cast <RealAttr1*> (target_new)) [*it] = (*target) [*it] - path2prediction (path);
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
  #ifndef NDEBUG
    if (verbose ())
    {
      setPrediction ();
      setAbsCriterion ();  
      ONumber on (cout, 6, false);
    //cout << "Optimizing " << ns. st->getName () << ": " << absCriterion_old << " -> " << absCriterion << endl;
      ASSERT (leReal (absCriterion, absCriterion_old));
    }
  #endif

    if (solved)
      prog (real2str (lr. absCriterion, 6));  // PAR
  }

  
  setPrediction ();
  setAbsCriterion ();  
#ifndef NDEBUG
  if (verbose ())
  {
    ONumber on (cout, 6, false);
    cout << absCriterion_old1 << " -> " << absCriterion << endl;
  }
  ASSERT (leReal (absCriterion, absCriterion_old1));
#endif
}



void DistTree::checkAttrPrediction () const
{
  if (! nodeAttrExist)
    return;

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
  for (Iterator it (dsSample); it ();)  
    if (! eqReal ((* const_cast <RealAttr1*> (prediction)) [*it], lr. predict (*it), dissim_ave * 1e-2))  // PAR 
    {
      const ONumber on (cout, dissimDecimals, true);
      cout << (* const_cast <RealAttr1*> (prediction)) [*it] << " " << lr. predict (*it) << "  objNum = " << *it << endl;
      cout << dissims [*it]. getObjName () << endl;
      cout << "Path: ";
      FOR (size_t, attrNum, dtNodes. size ())
        if ((* dtNodes [attrNum] -> attr) [*it])
          cout << " " << dtNodes [attrNum] -> getName ();
      cout << endl;
      print (cout);
      ERROR;
    }
}



void DistTree::optimize2 () 
{
  ASSERT (optimizable ());
  ASSERT (dissims. size () == 1);

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
  ASSERT (dissims. size () == 3);

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
      if (dissims [i]. hasLeaf (leaf))
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
  ASSERT (nodeAttrExist);
  ASSERT (ds. objs. size () >= 2 * name2leaf. size () - 2);

  
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
	   	prog ();
	 	  chron_getBestChange. start ();
		 	if (const Change* bestChange = getBestChange (node)) 
		  { 
		  	ASSERT (positive (bestChange->improvement));
		  	changes << bestChange;
		  }
		  chron_getBestChange. stop ();
 	  }
  }
    
  return applyChanges (changes);
}



void DistTree::optimizeIter (const string &output_tree)
{
  if (qc_on && verbose ())
  {
    checkAttrPrediction ();
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
      checkAttrPrediction ();
      checkAbsCriterion ("optimize");
    }
    if (verbose (1))
      cout. flush ();
  }
}



Real DistTree::optimizeSubgraph (const Steiner* center)  
{
  ASSERT (center);
  ASSERT (center->graph);
  ASSERT (& center->getTree () == this);
  ASSERT (optimizable ());
  ASSERT (! nodeAttrExist);


  chron_tree2subgraph. start ();
  Subgraph subgraph (*this);
  Node2Node new2old;  // Initially: newLeaves2boundary
  DistTree tree (center, areaRadius_std, subgraph, new2old);
  tree. qc ();
  ASSERT (subgraph. area_root->graph == this);
  ASSERT (subgraph. area. contains (center));
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
  ASSERT (tree. nodeAttrExist);
  chron_tree2subgraph. stop ();
  

  const TreeNode* root_old = root;
  const bool rootInArea = (subgraph. area_root == root);


  // Optimization
  chron_subgraphOptimize. start ();
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
        // Invoke: tree.neighborJoin(), tree.optimizeSubgraphs() if tree is large ??
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
      chron_subgraphOptimize. stop ();  
      return INF;
    }
  }
  tree. qc ();
  if (verbose ())
  {
    cout << "Subtree: ";
    tree. reportErrors (cout);
  }
  chron_subgraphOptimize. stop (); 
   

  const Real leafLen_min = tree. getMinLeafLen ();  

  Node2Node boundary2new (DiGraph::reverse (new2old));
  ASSERT (boundary2new. size () == new2old. size ());

  // tree.root
  // new2old: newBoundary2oldBoundary
  if (const DTNode* dtNode = static_cast <const DTNode*> (findPtr (boundary2new, subgraph. area_root)))
  {
  //ASSERT (! rootInArea);
    const Leaf* leaf = dtNode->asLeaf ();
    ASSERT (leaf);
    if (verbose ())
      cout << "Re-rooting subgraph ..." << endl;
    auto st = new Steiner (tree, const_static_cast <Steiner*> (leaf->getParent ()), leaf->len);
    EXEC_ASSERT (new2old. erase (leaf) == 1);
    tree. delayDeleteRetainArcs (const_cast <Leaf*> (leaf));
    ASSERT (isNan (static_cast <const DTNode*> (tree. root) -> len));
    const Steiner* root_ = st->makeDTRoot ();
    boundary2new [subgraph. area_root] = st;
    new2old [st] = subgraph. area_root;
    ASSERT (root_ != st);
    ASSERT (! isNan (root_->len));
    if (root_->isTransient ())
      tree. delayDeleteRetainArcs (const_cast <Steiner*> (root_));
    tree. toDelete. deleteData ();
    ASSERT (st == tree. root);
    ASSERT (isNan (st->len));
  }
  ASSERT (new2old. size () == boundary2new. size ());
  ASSERT (subgraph. boundary. size () == boundary2new. size ());
//tree. qc ();


  const Real absCriterion_old = absCriterion;
  

  // tree, subgraph --> topology

  // nodes: delete
  // center may be delete'd
  for (const TreeNode* node : subgraph. area)
    if (! contains (boundary2new, node))
    {
      DTNode* dtNode = const_static_cast <DTNode*> (node);
      ASSERT (! dtNode->attr);
      dtNode->isolate ();
      dtNode->detach ();
      toDelete << dtNode;
    }
    // Between boundary TreeNode's there are no Arc's

  // nodes, new2old[]: new
  for (const DiGraph::Node* node_new : tree. nodes)
  {
    DTNode* node = const_static_cast <DTNode*> (findPtr (new2old, node_new));
    if (! node)
    {
      node = new Steiner (*this, nullptr, NAN);
      new2old [node_new] = node;
      node->stable = true;
    }
    ASSERT (node);
    if (node_new != tree. root)
      node->len = static_cast <const DTNode*> (node_new) -> len;
  }
  ASSERT (new2old. size () == tree. nodes. size ());

  if (contains (boundary2new, center))
  {
    ASSERT (center->graph);
    ASSERT (& center->getTree () == this);
    const_cast <Steiner*> (center) -> stable = true;
  }
  
  borrowArcs (new2old, true);  // Arc's are not parallel
  
  // root
  if (rootInArea)
  {
    root = nullptr;
    for (const auto it : new2old)
    {
      const TreeNode* node = static_cast <const TreeNode*> (it. second);
      if (! node->getParent ())
      {
        ASSERT (! root);
        root = node;
      }
    }
  }
  else
  {
    ASSERT (root_old);
    ASSERT (! root_old->getParent ());
    root = root_old;
  }
  ASSERT (root);
   
  // Topology, absCriterion
  subgraph. subPaths2tree ();

  
  // deleteLenZero
  for (const auto it : new2old)
  {
    DTNode* node = const_static_cast <DTNode*> (it. second);
    deleteLenZero (node);
  }
  
  toDelete. deleteData ();


  if (qc_on && verbose ())
  {
    cout << "absCriterion = " << absCriterion << endl;
    checkAttrPrediction ();
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
  
  qc ();
  

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
  if (nodeAttrExist)
  {
    EXEC_ASSERT (optimizeLenArc ()); 
    finishChanges ();
    optimizeLenNode ();  
    finishChanges ();
    if (verbose (1))
      reportErrors (cout);
  }
#endif
}



const Change* DistTree::getBestChange (const DTNode* from) 
{
	ASSERT (from);
	
 	const Change* bestChange = nullptr;
 	
  VectorPtr<TreeNode> area;      area.     reserve (nodes. size ()); 
  VectorPtr<TreeNode> boundary;  boundary. reserve (nodes. size ()); 
  from->getArea (areaRadius_std, area, boundary);  
  if (verbose (1))
    cerr << " area=" << area. size () << " ";
    
 	for (const TreeNode* node : area)  
 	{
 		DTNode* to = const_static_cast <DTNode*> (node);
 	  if (Change::valid (from, to))
 	    tryChange (new Change (from, to), bestChange);
	  if (verbose ())
	  {
  	  ASSERT (to->graph);
	    to->qc ();
	  }
  }
  
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
    changes. sort (Change::strictlyLess);	
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
  {
    ONumber on (cout, 6, false);
    cout << "Improvement = " << improvement /*<< "  from: " << absCriterion_init << " to: " << absCriterion*/ << endl;
  }
  ASSERT (geReal (improvement, - absCriterion_delta));
  ASSERT ((bool) commits == (bool) improvement);

  if (commits && nodeAttrExist)
  {
    resetAttrs ();
    topology2attrs (); 
    setPrediction ();
    EXEC_ASSERT (optimizeLenArc ()); 
    finishChanges ();
    optimizeLenNode ();  
    finishChanges ();  
    if (verbose (1))
      reportErrors (cout);
  }
  
  qc ();
  
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



void DistTree::delayDeleteRetainArcs (DTNode* node)
{
	ASSERT (node);

  if (verbose ())
    cout << "To delete: " << node->getName () << endl;
    
  if (const Leaf* leaf = node->asLeaf ())
    name2leaf. erase (leaf->name);

  const VectorPtr<DiGraph::Node> children (node->getChildren ());
  for (const DiGraph::Node* child : children)
    const_static_cast <DTNode*> (child) -> len += (node == root ? NAN : node->len);

	delete node->attr;
	node->attr = nullptr;
	
	if (const Steiner* st = node->asSteiner ())
    if (const TreeNode* parent_ = st->getParent ())
    {
      const Steiner* parent = static_cast <const DTNode*> (parent_) -> asSteiner ();
      ASSERT (parent);
      for (Dissim &dissim : dissims)  // slow ??
        if (dissim. lca == st)
          dissim. lca = parent;
    }

	node->detachChildrenUp ();
	ASSERT (! node->graph);
	
	toDelete << node;
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
 	  if (deleteLenZero (static_cast <DTNode*> (node)))
 	    n++;
 			
  return n;  
}



bool DistTree::deleteLenZero (DTNode* node)
{
  if (const Steiner* s = node->asSteiner ())
  	if (   s->getParent () 
  		  && s->len == 0
  		  && s->childrenDiscernable ()
  		 )
  	{
  		delayDeleteRetainArcs (const_cast <Steiner*> (s));
  		return true;
  	}
 			
  return false;
}



void DistTree::removeLeaf (Leaf* leaf)
{
  ASSERT (leaf);
  
  if (! leaf->graph)
    return;
    
  ASSERT (! dissimDs. get ());
  ASSERT (! dissimAttr);
  
  const TreeNode* parent = leaf->getParent ();
  ASSERT (parent);
    
  // tree.detachedLeaves
  leaf->detachChildrenUp ();  
  ASSERT (! leaf->graph);
	delete leaf->attr;
	leaf->attr = nullptr;
	detachedLeaves << leaf;

  name2leaf. erase (leaf->name);

  // Clean topology, parent, Leaf::discernable consistency
  ASSERT (! parent->isLeaf ());
  if (const TreeNode* child = parent->isTransient ())
	{
    ASSERT (parent->graph); 
    if (const Leaf* childLeaf = static_cast <const DTNode*> (child) -> asLeaf ())
      const_cast <Leaf*> (childLeaf) -> discernable = true;
	  Steiner* st = const_static_cast <Steiner*> (parent);
	  parent = parent->getParent ();
		delayDeleteRetainArcs (st);
    if (! parent)
      parent = root;
  }  
    
  FOR (size_t, objNum, ds. objs. size ())
    if (dissims [objNum]. hasLeaf (leaf))
      const_cast <Obj*> (ds. objs [objNum]) -> mult = 0;
  dsSample = Sample (ds);

  if (verbose ())
    qcPaths (); 

  setAbsCriterion ();
  {
    Unverbose unv;
    ASSERT (parent);
    optimizeSubgraph (static_cast <const Steiner*> (parent->getAncestor (areaRadius_std)));  // PAR
  }

  toDelete. deleteData ();
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
  newRoot->pathObjNums = underRoot->pathObjNums;
  if (nodeAttrExist)
  {
    newRoot->addAttr ();
    ASSERT (newRoot->attr);
    ASSERT (underRoot->attr);
    * const_cast <CompactBoolAttr1*> (newRoot->attr) = * underRoot->attr;
  }
  underRoot->setParent (newRoot); 
  underRoot->len = arcLen;
  
  newRoot->makeDTRoot ();
  ASSERT (newRoot == root);
  ASSERT (root_ != root);
  
  setLca ();  

  if (root_->isTransient ())
    delayDeleteRetainArcs (root_);

  finishChanges ();

  if (verbose ())
    qcPaths ();  
}



Real DistTree::reroot (bool topological)
{
  DTNode* root_ = const_static_cast<DTNode*> (root);
  
  if (! root_->childrenDiscernable ())
    return 0;
  
  root_->setSubtreeLenUp (topological);
//cout << "Old radius: " << root_->getHeight_ave () << endl;  
  
  DTNode* bestDTNode = nullptr;
  Real bestDTNodeLen_new = NAN;
  WeightedMeanVar bestGlobalLen;
  root_->setGlobalLenDown (topological, bestDTNode, bestDTNodeLen_new, bestGlobalLen);
  ASSERT (bestDTNode);
  
  const Real height = bestGlobalLen. getMean ();
  
  if (verbose ())
  {
    cout << bestDTNode << " " << bestDTNodeLen_new << endl;
    cout << bestGlobalLen. getMean () << endl;
    print (cout);
  }
  
  reroot (bestDTNode, bestDTNodeLen_new);
  
  clearSubtreeLen ();
  
  return height;
}



Real DistTree::getMeanResidual () const
{
  ASSERT (optimizable ());

  Real s = 0;
  for (Iterator it (dsSample); it ();)  
  {
    const Real d = (*target) [*it];
  //ASSERT (d >= 0);
    const Real dHat = (*prediction) [*it];
    ASSERT (! isNan (dHat));
    s += it. mult * (dHat - d);
  }
  
  return s / dsSample. mult_sum;
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
  for (Iterator it (dsSample); it ();)  
  {
    const Real d = (*target) [*it];
  //ASSERT (d >= 0);
    const Real dHat = (*prediction) [*it];
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
  for (Iterator it (dsSample); it ();)  
  {
    const Real d = (*target) [*it];
  //ASSERT (d >= 0);
    const Real dHat = (*prediction) [*it];
    ASSERT (! isNan (dHat));
    if (nullReal (dHat))
    {
      epsilon2_0 += it. mult * sqr (d);
      continue;
    }

    ASSERT (positive (d));
    const Real a = it. mult * sqr (dHat - d) / dHat;
    if (! DM_sp::finite (a))
    {
      cout << dissims [*it]. getObjName () << ": " << d << " " << dHat << endl;
      ERROR;
    }
    
    const VectorPtr<TreeNode> path (dissim2path (*it));
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
    dtNode->errorDensity = dtNode->paths ? sqrt (dtNode->errorDensity / (Real) dtNode->paths) : INF;
  }
  
  return epsilon2_0;
}



VectorPtr<Leaf> DistTree::findOutliers (bool strong,
                                        Real &outlier_min) const
{
	Dataset leafDs;
  auto criterionAttr = new RealAttr1 ("LeafCriterion", leafDs);
  
  leafDs. objs. reserve (name2leaf. size ());
  for (const auto& it : name2leaf)
    if (it. second->graph)
    {
      const Real relErr = it. second->getRelCriterion (strong);
      if (! isNan (relErr))
      {
    	  const size_t index = leafDs. appendObj (it. first);
    	  (*criterionAttr) [index] = relErr;
      }
    }
  
  const Sample sample (leafDs);
  outlier_min = criterionAttr->normal2outlier (sample, 0.1);  // PAR  

  VectorPtr<Leaf> res;
  if (! isNan (outlier_min))
    for (const auto& it : name2leaf)
      if (it. second->graph)
        if (geReal (it. second->getRelCriterion (strong), outlier_min))
          res << it. second;
      
	return res;
}
  


RealAttr1* DistTree::getResiduals2 () 
{
  ASSERT (optimizable ());

  Common_sp::AutoPtr<RealAttr1> resid2Attr (new RealAttr1 ("resid2", ds, target->decimals + 1));
  for (Iterator it (dsSample); it ();)  
  {
    const Real d = (*target) [*it];
  //ASSERT (d >= 0);
    const Real dHat = (*prediction) [*it];
    ASSERT (! isNan (dHat));
    const Real residual = abs (dHat - d);      
    ASSERT (DM_sp::finite (residual));
    (*resid2Attr) [*it] = sqr (residual);
  }
    
  return resid2Attr. release ();
}



RealAttr1* DistTree::getLogPredictionDiff () 
{
  ASSERT (optimizable ());

  Common_sp::AutoPtr<RealAttr1> logDiffAttr (new RealAttr1 ("logDiff", ds, target->decimals + 1));
  for (Iterator it (dsSample); it ();)  
  {
    const Real d = (*target) [*it];
  //ASSERT (d >= 0);
    const Real dHat = (*prediction) [*it];
    ASSERT (! isNan (dHat));
    if (   positive (d)
        && positive (dHat)
       )
      (*logDiffAttr) [*it] = log (d) - log (dHat);
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



void DistTree::saveDissim (ostream &os) const
{
  ASSERT (optimizable ());

  const ONumber on (os, dissimDecimals, true);
  FOR (size_t, objNum, ds. objs. size ())
  {
    const Dissim& dissim = dissims [objNum];
    if (! dissim. leaf1->graph)
      continue;
    if (! dissim. leaf2->graph)
      continue;
    if (ds. objs [objNum] -> mult || (*target) [objNum] == 0)
      os         << dissim. leaf1->name
         << '\t' << dissim. leaf2->name
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
  
//IMPLY (leafLen == 0 && arcLen == 0, anchor->asLeaf ());
  
  ASSERT (absCriterion_leaf >= 0);  
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
  ASSERT (leaf);
  ASSERT (dissim >= 0);
  ASSERT (anchor);
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
  {
    cout << "File \"" << getRequestFName () << "\" exists" << endl;
    return;
  }
    
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
      string name1, name2, dissimS;
      while (f. nextLine ())
        try
        {
          istringstream iss (f. line);
          iss >> name1 >> name2 >> dissimS;
          ASSERT (iss. eof ());
          ASSERT (name1 != name2);
          if (name1 != name)
            throw runtime_error ("First object in " + getDissimFName () + " must be " + name);
          const Real dissim = str2real (dissimS);
          if (dissim < 0)
            throw runtime_error ("Dissimilarity must be non-negative");
          leaf2dissims << Leaf2dissim (findPtr (tree. name2leaf, name2), dissim, location. anchor);
        }          
        catch (...)
        {
          cout << f. line << endl;
          throw;
        }
      leaf2dissims. sort ();
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
    	  if (! leaf2dissims. containsFast (ld))
    	    requested << descendant->reprLeaf;
    	}
      ancestor = static_cast <const DTNode*> (ancestor->getParent ());
    }
  }
  
  {
    OFStream of (getRequestFName ());
    for (const Leaf* leaf : requested)
    {
      ASSERT (name != leaf->name);
      of << name << '\t' << leaf->name << endl;
    }
  }
}



void NewLeaf::optimize ()
{
  bool dissimExists = false;
	for (const Leaf2dissim& ld : leaf2dissims)
	  if (ld. dissim != INF)
  	{
  	  dissimExists = true;
  	  if (ld. dissim == 0)  // if all dissim = INF ??
  	  {
  	    location. anchor = ld. leaf;
  	    location. leafLen = 0;
  	    location. arcLen = 0;
  	    location. absCriterion_leaf = 0;
  	    return;
  	  }
  	}
  	
  if (! dissimExists)
  {
    location. anchor = static_cast <const DTNode*> (tree. root);
    location. leafLen = INF;
    location. arcLen = INF;
    location. absCriterion_leaf = 0;
    return;
  }

  Location location_best (location);
  Vector<Leaf2dissim> leaf2dissims_best (leaf2dissims);
  optimizeAnchor (location_best, leaf2dissims_best);
  location = location_best;
  leaf2dissims = leaf2dissims_best;
}



void NewLeaf::optimizeAnchor (Location &location_best,
                              Vector<Leaf2dissim> &leaf2dissims_best)
{
  setLocation ();
  if (minimize (location_best. absCriterion_leaf, location. absCriterion_leaf))
  {
    location_best = location;
    leaf2dissims_best = leaf2dissims;
  }

  if (! location. anchor->childrenDiscernable ())
    return;
    
  VectorPtr<DiGraph::Arc> arcs;  arcs. reserve (location. anchor->arcs [false]. size ());
  insertAll (arcs, location. anchor->arcs [false]);
	for (const DiGraph::Arc* arc : arcs)
	{
	  const DTNode* child = static_cast <const DTNode*> (arc->node [false]);
    ASSERT (! child->inDiscernable ());
    const Keep<Location> location_old (location);
	  const Keep<Vector<Leaf2dissim>> leaf2dissims_old (leaf2dissims);
	  if (descend (child))
	    optimizeAnchor (location_best, leaf2dissims_best);
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
	  if (ld. mult)
	  {
  	  mult_sum  += ld. mult;
  	  u_avg     += ld. mult * ld. getU ();
  	  delta_avg += ld. mult * ld. getDelta ();
	  }
	}	
	ASSERT (positive (mult_sum));
	u_avg     /= mult_sum;
	delta_avg /= mult_sum;
	
	Real u_var = 0;
	Real delta_u_cov = 0;
	for (const Leaf2dissim& ld : leaf2dissims)
	  if (ld. mult)
  	{
  	  delta_u_cov += ld. mult * (ld. getDelta () - delta_avg) * (ld. getU () - u_avg);
  	  u_var       += ld. mult * sqr (ld. getU () - u_avg);
  	}
	u_var       /= mult_sum;
	delta_u_cov /= mult_sum;  	

  location. arcLen = 0;
  if (location. anchor == tree. root)
    {	ASSERT (u_var == 0); }
  else
    if (u_var > 0)  // dissimilarity to all leaves may be INF
    	location. arcLen = min (max (0.0, delta_u_cov / u_var), location. anchor->len);

	location. leafLen = max (0.0, delta_avg - u_avg * location. arcLen);
	
	location. absCriterion_leaf = 0;
	for (const Leaf2dissim& ld : leaf2dissims)
	  if (ld. mult)
  	  location. absCriterion_leaf += ld. mult * sqr (ld. getEpsilon (*this));
	ASSERT (location. absCriterion_leaf >= 0);
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
	  if (ld. mult)
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
