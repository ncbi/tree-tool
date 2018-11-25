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
const StringVector varianceTypeNames {"lin", "sqr", "exp", "linExp"};
VarianceType varianceType = varianceType_linExp;

Real dissim_coeff = 1;
Real hybridness_min = 1.1;
Real dissim_boundary = NaN;




// Neighbor

Neighbor::Neighbor (const Leaf* leaf_arg,
                    Real target_arg)
: leaf (leaf_arg)
, target (target_arg)
{ 
	ASSERT (leaf);
	ASSERT (leaf->graph);
	ASSERT (target > 0); 
}



Neighbor::Neighbor (const Leaf* leaf_arg)
: leaf (leaf_arg)
{ 
	ASSERT (leaf);
	ASSERT (leaf->graph);
}




// Triangle

Triangle::Triangle (const Leaf* child_arg,
		  	  	        Real hybridness_arg,
								  	const Leaf* parent1,
								  	const Leaf* parent2,
								  	Real parent1_dissim,
								  	Real parent2_dissim)
: child (child_arg)
, hybridness (hybridness_arg)
{
	ASSERT (parent1);
	ASSERT (parent2);
	parents [0]. leaf = parent1;
	parents [1]. leaf = parent2;
	parents [0]. dissim = parent1_dissim;
	parents [1]. dissim = parent2_dissim;
	if (parents [0]. leaf->name > parents [1]. leaf->name)
		swap (parents [0], parents [1]);
}



void Triangle::qc () const
{
	if (! qc_on)
		return;
		
	ASSERT (child);
	ASSERT (child->graph);
	for (const bool i : {false, true})
	{
		const Parent& p = parents [i];
		ASSERT (p. leaf);	
		ASSERT (p. leaf->graph);
 		ASSERT (child->getDiscernible () != p. leaf->getDiscernible ());
		ASSERT (p. dissim > 0);
 	}
 	ASSERT (! (   child_hybrid 
 	           && parents [false]. hybrid 
 	           && parents [true].  hybrid 
 	          )
 	       );
  ASSERT (parents [0]. leaf->getDiscernible () != parents [1]. leaf->getDiscernible ());
  ASSERT (parents [0]. leaf->name < parents [1]. leaf->name);
	ASSERT (hybridness >= hybridness_min);
}



void Triangle::print (ostream &os) const
{ 
	TabDel td (3, false);  // PAR  
	const Parent& p1 = parents [0];
	const Parent& p2 = parents [1];
	td << child->name << hybridness 
	   << p1. leaf->name << p2. leaf->name << p1. dissim << p2. dissim 
	   << child_hybrid << p1. hybrid << p2. hybrid;
	os << td. str () << endl;	
}



void Triangle::qcMatchHybrids (const VectorPtr<Leaf> &hybrids) const
{
	#define QC_MATC_HHYBRIDS(leaf,isHybrid)  if ((isHybrid) != hybrids. containsFast (leaf)) { cout << (leaf)->name << endl; print (cout); ERROR; }	
	QC_MATC_HHYBRIDS (child, child_hybrid);
	QC_MATC_HHYBRIDS (parents [0]. leaf, parents [0]. hybrid);
	QC_MATC_HHYBRIDS (parents [1]. leaf, parents [1]. hybrid);
	#undef QC_MATC_HHYBRIDS
}




// TriangleParentPair

void TriangleParentPair::setTriangles (const DistTree &tree) 
{
	for (const bool i : {false, true})
    ASSERT (parents [i]. leaf);
  ASSERT (parents [0]. leaf->name < parents [1]. leaf->name);
  ASSERT (parentsDissim > 0);
  ASSERT (hybridness_min > 1);
  ASSERT (triangles. empty ());

  Vector<Neighbor> hybridParents;  hybridParents. reserve (parents [0]. leaf->pathObjNums. size ());
	for (const uint objNum : parents [0]. leaf->pathObjNums)
	{
		const Dissim& dissim = tree. dissims [objNum];
		if (! dissim. mult)
			continue;
		const Leaf* other = dissim. getOtherLeaf (parents [0]. leaf);
	  ASSERT (other->graph);
    hybridParents << Neighbor (other, dissim. target);
	}
	hybridParents. sort ();
	ASSERT (hybridParents. isUniq ());

	for (const uint objNum : parents [1]. leaf->pathObjNums)
	{
		const Dissim& dissim = tree. dissims [objNum];
		if (! dissim. mult)
			continue;
		const Leaf* other = dissim. getOtherLeaf (parents [1]. leaf);
    ASSERT (other->graph);
  	const size_t i = hybridParents. binSearch (Neighbor (other));
  	if (i == NO_INDEX)
  		continue;
  	const Real hybridness = parentsDissim / (dissim. target + hybridParents [i]. target);
  	if (hybridness >= hybridness_min)
  		triangles << Triangle (other, hybridness, parents [0]. leaf, parents [1]. leaf, hybridParents [i]. target, dissim. target);
	}
	ASSERT (! triangles. empty ());

  // Statistics
	MeanVar hybridness_mv;
	for (const Triangle& tr : triangles)
		hybridness_mv << tr. hybridness;
	hybridness_ave = hybridness_mv. getMean ();
	
	ASSERT (hybridness_ave >= hybridness_min);
}



void TriangleParentPair::finish (const DistTree &tree,
                                 const Set<const Leaf*> &hybrids)
{ 
	for (const bool i : {false, true})
	{
    ASSERT (parents [i]. leaf);
		ASSERT (! parents [i]. classSize);
  }
	ASSERT (! triangles. empty ());
	ASSERT (triangle_best_index == NO_INDEX);
	

 	triangles. filterValue ([&hybrids] (const Triangle &tr) { return hybrids. contains (tr. child); });	
	if (triangles. empty ())
		return;

	Triangle* best = nullptr;
	Real hybridness = 0;
	for (Triangle& tr : triangles)
		if (! best || maximize (hybridness, tr. hybridness))
			best = & tr;
	ASSERT (best);
	triangle_best_index = (size_t) (best - & triangles [0]);
  ASSERT (triangle_best_index < triangles. size ());
  ASSERT (& getBest () == best);

	for (const bool i : {false, true})
		parents [i]. classSize = child_parent2parents (tree, best->child, best->parents [! i]. leaf, best->parents [! i]. dissim);


  // Decision
	if (dissimError ())
		return;


  // ??
  // Remove a hybrid candidate class, compute criterion
  // Select the hybrid candidate class which minimizes average absCriterion
  

  // PAR
  constexpr Prob classSizeFrac_max = 0.2;  

	const Real child_classSize = (Real) triangles. size ();	
  const Real all = child_classSize + (Real) parents [0]. classSize + (Real) parents [1]. classSize;

	for (const bool i : {false, true})  
	  if (   best->parents [! i]. dissim > dissim_boundary
	  	  && best->parents [  i]. dissim < dissim_boundary
	  	 ) 
	  {
		  const Parent& p = parents [i];		
			if (child_classSize < all * classSizeFrac_max)
				setChildrenHybrid ();
			else if (p. classSize < all * classSizeFrac_max)
				best->parents [i]. hybrid = true;		
			else if (child_classSize + (Real) p. classSize < all * classSizeFrac_max)
			{
				setChildrenHybrid ();
				best->parents [i]. hybrid = true;		
			}
			return;
	  }
  
  {
		const Prob p = child_classSize / all;
		if (p <= classSizeFrac_max)
		{
			// aa bb ab
			setChildrenHybrid ();
			return;
		}
		if (p >= 1 - classSizeFrac_max)
		{
			// ab ac aa
			// common(ab,ac) = 1/4 aa <= random halves of aa
			best->parents [0]. hybrid = true;
			best->parents [1]. hybrid = true;
			return;  			
		}
	}
	
	// Singletons
  bool allSingletons = false;
	if (triangles. size () == 1)
	{
	  allSingletons = true;
  	for (const bool i : {false, true})  
  		if (parents [i]. classSize > 1)
  			allSingletons = false;
  }
  if (allSingletons)
  {
    size_t i_bad = NO_INDEX;
    Real badCriterion = triangles [0]. child->badCriterion;
  	for (const bool i : {false, true})  
  		if (maximize (badCriterion, parents [i]. leaf->badCriterion))
  		  i_bad = (size_t) i;
  	if (i_bad == NO_INDEX)
  	  setChildrenHybrid ();
  	else  
  	{
  	  ASSERT (i_bad < 2);
  	  best->parents [i_bad]. hybrid = true;	
  	}
  }
  else
  {
  	if (triangles. size () == 1)
  		setChildrenHybrid ();
  	for (const bool i : {false, true})  
  		if (parents [i]. classSize == 1)
  			best->parents [i]. hybrid = true;	
  }
			
	// Undecided cases ??
}



void TriangleParentPair::qc () const
{
	if (! qc_on)
		return;

	for (const bool i : {false, true})  
		ASSERT (parents [i]. leaf);
	ASSERT (parents [0]. leaf->name < parents [1]. leaf->name);
	ASSERT (parentsDissim > 0);
	
	if (triangles. empty ())
		return;
	
  bool child_hybrid;
  bool first = true;
	for (const Triangle& tr : triangles)
	{
		tr. qc ();
		for (const bool i : {false, true})  
			ASSERT (tr. parents [i]. leaf == parents [i]. leaf);
		ASSERT (eqReal (parentsDissim, tr. parentsDissim ()));
		if (first)
		{
		  child_hybrid = tr. child_hybrid;
		  first = false;
		}
		else
		  ASSERT (child_hybrid == tr. child_hybrid);
	}
	
	ASSERT (& getBest ());

	for (const bool i : {false, true})  
		ASSERT (parents [i]. classSize);

	ASSERT (hybridness_ave >= hybridness_min);
}



void TriangleParentPair::print (ostream &os) const
{ 
	TabDel td (3, false);  // PAR
	const Parent& p1 = parents [0];
	const Parent& p2 = parents [1];
	td << getBest (). child->name << p1. leaf->name << p2. leaf->name 
	   << triangles. size () << p1. classSize << p2. classSize
	   << hybridness_ave << getBest (). parents [0]. dissim << getBest (). parents [1]. dissim
	   << getBest (). child_hybrid << getBest (). parents [0]. hybrid << getBest (). parents [1]. hybrid;
	os << td. str () << endl;
}



size_t TriangleParentPair::child_parent2parents (const DistTree &tree,
					                                       const Leaf* child,
								                                 const Leaf* parent,
									                               Real parentDissim) const
{
  ASSERT (child);
  ASSERT (parent);
  ASSERT (child != parent);
  ASSERT (parentDissim > 0);
  ASSERT (hybridness_min > 1);
  	
  Vector<Neighbor> hybridParents;  hybridParents. reserve (child->pathObjNums. size ());
	for (const uint objNum : child->pathObjNums)
	{
		const Dissim& dissim = tree. dissims [objNum];
		if (! dissim. mult)
			continue;
		const Leaf* other = dissim. getOtherLeaf (child);
	  ASSERT (other->graph);
    hybridParents << Neighbor (other, dissim. target);
	}
	hybridParents. sort ();
	ASSERT (hybridParents. isUniq ());

  size_t n = 0;
	for (const uint objNum : parent->pathObjNums)
	{
		const Dissim& dissim = tree. dissims [objNum];
		if (! dissim. mult)
			continue;
		const Leaf* otherParent = dissim. getOtherLeaf (parent);
    ASSERT (otherParent->graph);
  	const size_t i = hybridParents. binSearch (Neighbor (otherParent));
  	if (i == NO_INDEX)
  		continue;
  	const Real hybridness = dissim. target / (parentDissim + hybridParents [i]. target);
  	if (hybridness >= hybridness_min)  // otherParent's may be too diverse to be all hybrid
  		n++;
	}

  return n;
}




//

static const string rel_errorS ("rel_error");



// DTNode

DTNode::DTNode (DistTree &tree,
                Steiner* parent_arg,
    	          Real len_arg)
: TreeNode (tree, parent_arg)  // DistTree must be declared
, len (isNan (len_arg) ? len_arg : max (0.0, len_arg))
{}



void DTNode::qc () const
{ 
  if (! qc_on)
    return;
	TreeNode::qc ();

  if (graph)
  {
    ASSERT ((bool) getParent () == ! isNan (len));
    if (! childrenDiscernible ())
    {
      ASSERT (! inDiscernible ());
    	for (const DiGraph::Arc* arc : arcs [false])
    	  ASSERT (static_cast <const DTNode*> (arc->node [false]) -> inDiscernible ());
    }
  }
  
  IMPLY (! isNan (len), len >= 0);	  
  IMPLY (paths, errorDensity >= 0);
  ASSERT (! contains (name, '\''));
}



void DTNode::saveContent (ostream& os) const
{ 
  {
    const ONumber oLen (os, dissimDecimals, true);
    os << "len=" << len;
  }
  const ONumber oNum (os, criterionDecimals, true);  // PAR
  if (paths)
	  os << "  err_density=" << errorDensity << "  paths=" << paths;
  if (! isNan (absCriterion_ave))
    os << "  " << rel_errorS << "=" << getRelCriterion ();
}



const DistTree& DTNode::getDistTree () const
{
  return static_cast <const DistTree&> (getTree ());
}



const Leaf* DTNode::inDiscernible () const
{ 
  const Leaf* g = asLeaf ();
  return g && ! g->discernible ? g : nullptr; 
}



Real DTNode::getRelCriterion () const
{ 
  return absCriterion_ave / getDistTree (). getAbsCriterion_ave (); 
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
    Real lenDelta = NaN;
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
    
    const Real zScore = 3;  // PAR
    if (! inDiscernible ())
      if (   ! bestDTNode 
          || bestGlobalLen. getOutlier_min (zScore) /*getMean ()*/ > globalLen. getOutlier_min (zScore) /*getMean ()*/
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
    FFOR (size_t, i, subtreeLeaves. size ())
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



Vector<uint> DTNode::getLcaObjNums () 
{
  Vector<uint> lcaObjNums;  
  if (getDistTree (). dissims. empty ())
    return lcaObjNums;
    
  lcaObjNums. reserve ((size_t) log (getDistTree (). dissims. size ()));  // PAR

#if 0
  Vector<bool> childDissims (tree. dissims. size (), false);  // Faster for small n
#endif

  const VectorPtr<DiGraph::Node> children (getChildren ());
  FOR_REV (size_t, i, children. size ())
  {
    Vector<uint>& childPathObjNums = const_static_cast <DTNode*> (children [i]) -> pathObjNums;
  #if 0
    // Faster for small n
    for (const size_t objNum : childPathObjNums)
      if (childDissims [objNum])
        lcaObjNums << objNum;
      else
        childDissims [objNum] = true;
  #else
    childPathObjNums. sort ();
    FOR (size_t, j, i)
      for (const uint objNum : static_cast <const DTNode*> (children [j]) -> pathObjNums)
        if (childPathObjNums. containsFast (objNum))
          lcaObjNums << objNum;
  #endif
  }
        
  return lcaObjNums;
}



VectorPtr<Leaf> DTNode::getSparseLeafMatches (size_t depth_max,
                                              bool subtractDissims,
                                              bool refreshDissims) const
{
  IMPLY (depth_max, depth_max >= sparsingDepth);
  IMPLY (subtractDissims, asLeaf ());

  map <const Steiner* /*lca*/, VectorPtr<Leaf>>  lca2leaves;  
  if (subtractDissims)
  {
  	// Time: ~ O(p/n) ~ O(log(n))
	  for (const uint objNum : pathObjNums)
	  {
	  	const Dissim& dissim = getDistTree (). dissims [objNum];
	  	if (! dissim. mult)
	  		continue;
	  	ASSERT (dissim. lca);
	  	lca2leaves [dissim. lca] << dissim. getOtherLeaf (asLeaf ());
	  }
  	// Time: ~ O(log(n) p/(n log(n))) = O(p/n) ~ O(log(n))
	  for (auto& it : lca2leaves)
	  {
	  	VectorPtr<Leaf>& otherLeaves = it. second;
	  	ASSERT (! otherLeaves. empty ());
	  	otherLeaves. sort ();
	  	ASSERT (otherLeaves. isUniq ());
	  }
	}

  VectorPtr<Leaf> matches;  matches. reserve (getDistTree (). getSparseDissims_size ());  // PAR
 	VectorPtr<DTNode> descendants;  descendants. reserve (2 * powInt (2, (uint) sparsingDepth));  // PAR
  size_t depth = 0;  
  const DTNode* ancestor = this;
  const DTNode* ancestor_prev = nullptr;
  // Time: ~ O(log^2(n))
  while (ancestor) 
  {
  	descendants. clear ();
  	{
	  	size_t searchDepth = sparsingDepth;
	  	if (depth_max && depth > depth_max)
	  	{
	  	  const size_t dec = depth - depth_max;
	  	  if (searchDepth > dec)
	  	    searchDepth -= dec;
	  	  else
	  	    searchDepth = 1;
	  	}
	    maximize<size_t> (searchDepth, 3);  // PAR
	  	ASSERT (searchDepth);
	  	ancestor->getDescendants (descendants, searchDepth, ancestor_prev);  
	  }
	  ASSERT (! descendants. empty ());
	  {
		  size_t added = 0;
		  const VectorPtr<Leaf>* otherLeaves = nullptr;  // having a dissimilarity with *this where LCA = ancestor
		  if (subtractDissims)
		  	if (const Steiner* st = ancestor->asSteiner ())
		  	  otherLeaves = findPtr (lca2leaves, st);
	  	for (const DTNode* descendant : descendants)
	  	{
	  	  if (! refreshDissims)
  	  		if (otherLeaves && otherLeaves->size () + added >= descendants. size ())  // maybe not through the right children
  	  			break;  
	  		const Leaf* repr = descendant->getReprLeaf ();
	  		ASSERT (repr);
	  		if (   ! otherLeaves 
	  			  || ! otherLeaves->containsFast (repr)
	  			 )
				{
	        matches << repr;
	        added++;
	      }
	 	  }
	 	}
	 	ancestor_prev = ancestor;
    ancestor = static_cast <const DTNode*> (ancestor->getParent ());
    depth++;
  }
  matches. sort ();
  matches. uniq ();

  return matches;
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
}



void Steiner::saveContent (ostream& os) const 
{ 
  DTNode::saveContent (os);
  if (! name. empty ())
    os << "  name='" << name << '\'';
#if 0
	if (frequentChild)
	  os << " " << "frequent";
	os << " fd=" << frequentDegree;
#endif
}



const Leaf* Steiner::getReprLeaf () const
{ 
	const size_t n = arcs [false]. size ();
	ASSERT (n);
	const size_t index = getDistTree (). rand. get (n);
	return static_cast <const DTNode*> (arcs [false]. at (index) -> node [false]) -> getReprLeaf ();
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
                              size_t depth,
                              const DTNode* exclude) const
{
  if (depth && childrenDiscernible ())
    for (const DiGraph::Arc* arc : arcs [false])
    {
      const DTNode* child = static_cast <const DTNode*> (arc->node [false]);
      if (child != exclude)
        child->getDescendants (descendants, depth - 1, exclude);
    }
  else  
    if (this != exclude)
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
  const auto pathObjNums_new = ancestor2descendant->pathObjNums;
  
  reverseParent (ancestor2descendant, nullptr);

  setParent (const_cast <DTNode*> (parent_new));
  // Arc-specific data
  len         = len_new;
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



LeafCluster Steiner::getIndiscernibles ()
{
  ASSERT (! getDistTree (). subDepth);
  ASSERT (getDistTree (). optimizable ());
  
  LeafCluster leafCluster;

  if (childrenDiscernible ())
    return leafCluster;

  VectorPtr<Leaf> children;  children. reserve (arcs [false]. size ());
  for (const DiGraph::Arc* arc : arcs [false])
  {
    const Leaf* leaf = static_cast <const DTNode*> (arc->node [false]) -> asLeaf ();
    ASSERT (leaf);
    var_cast (leaf) -> DisjointCluster::init ();
    children << leaf;
 	}
 	children. sort ();
 	ASSERT (children. isUniq ());

  for (const DiGraph::Arc* arc : arcs [false])
  {
    const Leaf* leaf = static_cast <const DTNode*> (arc->node [false]) -> asLeaf ();
    ASSERT (leaf);
    for (const uint objNum : leaf->pathObjNums)
    {
      const Dissim& dissim = getDistTree (). dissims [objNum];
   	  if (   dissim. valid ()
   	      && dissim. target <= 0
   	      && children. containsFast (dissim. leaf1)
   	      && children. containsFast (dissim. leaf2)
   	     )
        var_cast (dissim. leaf1) -> DisjointCluster::merge (* var_cast (dissim. leaf2));
    }
  }

  for (const DiGraph::Arc* arc : arcs [false])
  {
    const Leaf* leaf = static_cast <const DTNode*> (arc->node [false]) -> asLeaf ();
    ASSERT (leaf);
    leafCluster [var_cast (leaf) -> DisjointCluster::getDisjointCluster ()] << leaf;
 	}
  
  return leafCluster;
}



#if 0
void Steiner::setSubTreeWeight ()
{
	subTreeWeight = lcaNum;
	for (const DiGraph::Arc* arc : arcs [false])
	  if (const Steiner* child = static_cast <DTNode*> (arc->node [false]) -> asSteiner ())
	  {
		  const_cast <Steiner*> (child) -> setSubTreeWeight ();
		  subTreeWeight += child->subTreeWeight;
		}
}



void Steiner::threadNum2subTree (size_t threadNum_arg)
{
	ASSERT (threadNum_arg != NO_INDEX);
	
  threadNum = threadNum_arg;
	for (const DiGraph::Arc* arc : arcs [false])
	  if (Steiner* child = const_cast <Steiner*> (static_cast <DTNode*> (arc->node [false]) -> asSteiner ()))
		  child->threadNum2subTree (threadNum_arg);
}



void Steiner::threadNum2ancestors (size_t threadNum_arg)
{
	ASSERT (threadNum_arg != NO_INDEX);

	Steiner* n = this;
  while (Steiner* parent = const_static_cast<Steiner*> (n->getParent ()))
  {
  	parent->threadNum = threadNum_arg;
  	n = parent;
  }
}
#endif




// Leaf

const string Leaf::non_discernible ("non-discernible");



Leaf::Leaf (DistTree &tree,
  	        Steiner* parent_arg,
  	        Real len_arg,
  	        const string &name_arg)
: DTNode (tree, parent_arg, len_arg)  // DistTree must be declared
, index (getDistTree (). leafNum)
{
  ASSERT (name. empty ());
  name = name_arg;
  const_cast <DistTree&> (getDistTree ()). leafNum ++;  
}



#if 0
Leaf::Leaf (DistTree &tree,
  	        Leaf* other,
  	        const string &name_arg)
: DTNode ( tree   // DistTree must be declared
         , const_static_cast <Steiner*> (checkPtr<Leaf> (other) -> getParent ())  // temporary
         , 0
         )  
, name (name_arg)
, discernible (other->discernible)
, index (getDistTree (). leafNum)
{
  const_cast <DistTree&> (getDistTree ()). leafNum ++;  
  collapse (other);
}
#endif



void Leaf::qc () const
{ 
  if (! qc_on)
    return;
	DTNode::qc ();
	  
  if (graph)
  {
    ASSERT (getParent ());
  	ASSERT (isLeaf ());	
    if (getParent () && discernible != static_cast <const DTNode*> (getParent ()) -> childrenDiscernible ())
    {
      cout << getName () << " " << discernible << " " << getParent () -> getName () 
               << " " << static_cast <const DTNode*> (getParent ()) -> childrenDiscernible () << endl;
      ERROR;
    }
    for (const uint objNum : pathObjNums)
      ASSERT (getDistTree (). dissims [objNum]. hasLeaf (this));
    IMPLY (getDistTree (). subDepth, discernible);
  }

  ASSERT (! name. empty());
  ASSERT (! isLeft (name, "0x"));
  if (! isNan (len) && ! discernible && len != 0)
  {
    cout << getName () << " " << len << endl;
    ERROR;
  }  
  ASSERT (getReprLeaf () == this);
  
  ASSERT (index < NO_INDEX);
}



void Leaf::setLca ()
{
  DTNode::setLca ();

  for (const uint objNum : pathObjNums)
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



const DTNode* Leaf::getDiscernible () const
{ 
  if (const Leaf* leaf = inDiscernible ())
    return static_cast <const DTNode*> (leaf->getParent ());
  return this;
}



void Leaf::collapse (Leaf* other)
{
  ASSERT (! getDistTree (). subDepth);
  ASSERT (this != other);
//ASSERT (len == 0);
//ASSERT (discernible);

  if (! other)
    return;

  ASSERT (graph);
  ASSERT (graph == other->graph);
  
#if 0
  printAncestors (nullptr);
  other->printAncestors (nullptr);
  cout << discernible << ' ' << other->discernible << endl;
  cout << name << ' ' << other->name << endl;
  getParent () -> saveText (cout);
#endif
  
  Vector<Leaf*> indiscernibles (1, this);
  {
    const Steiner* oldParent = static_cast <const Steiner*> (getParent ()); 
    ASSERT (oldParent);
    if (! discernible)   //  discernible == oldParent->childrenDiscernible ()
    {
      indiscernibles. clear ();
      indiscernibles. reserve (oldParent->arcs [false]. size ());
      for (const DiGraph::Arc* arc : oldParent->arcs [false])
      {
        const Leaf* child = static_cast <const DTNode*> (arc->node [false]) -> asLeaf ();
        ASSERT (child);
        ASSERT (! child->discernible);
        indiscernibles << const_cast <Leaf*> (child);
      }
      ASSERT (indiscernibles. size () >= 2);
    }
  }
  ASSERT (indiscernibles. contains (this));
#if 0
  for (const Leaf* leaf : indiscernibles)
    cout << ' ' << leaf;
  cout << endl;
#endif
    
  if (other->discernible)
  {
    auto st = new Steiner ( const_cast <DistTree&> (getDistTree ())
                          , const_static_cast <Steiner*> (other->getParent ())
                          , other->len
                          );
    other->setParent (st);
    for (Leaf* leaf : indiscernibles)
      leaf->setParent (st);
    other->len = 0;
    other->discernible = false;
  }
  else
  { 
    ASSERT (other->len == 0);
    if (getParent () != other->getParent ())
      for (Leaf* leaf : indiscernibles)
        leaf->setParent (const_static_cast <Steiner*> (other->getParent ()));
  }
  
  len = 0;
  discernible = false;

  ASSERT (getParent () == other->getParent ());
#ifndef NDEBUG
  for (const DiGraph::Arc* arc : getParent () -> arcs [false])
  {
    const Leaf* child = static_cast <const DTNode*> (arc->node [false]) -> asLeaf ();
    ASSERT (child);
    ASSERT (! child->discernible);
  }
#endif
}



Real Leaf::getHybridness (const Leaf* &parent1,
	                        const Leaf* &parent2,
	                        Real &parent1_dissim,
	                        Real &parent2_dissim) const
{
  ASSERT (hybridness_min > 1);
	ASSERT (graph);
	
	parent1 = nullptr;
	parent2 = nullptr;
	parent1_dissim = NaN;
	parent2_dissim = NaN;
	
	
  const auto& dissims = getDistTree (). dissims;

  Vector<Neighbor> hybridParents;  hybridParents. reserve (pathObjNums. size ());
    // Neighbor::leaf is unique
	  // Neighbor::target: dissimilarity between *this and leaf
	for (const uint objNum : pathObjNums)
	{
		const Dissim& dissim = dissims [objNum];
	  if (dissim. mult)  
  	{
  		const Leaf* other = dissim. getOtherLeaf (this);
  	  ASSERT (other->graph);  
      if (dissim. prediction < dissim. target) 
		    hybridParents << Neighbor (other, dissim. target);
		}
	}
	hybridParents. sort ();
	ASSERT (hybridParents. isUniq ());
	
	Real hybridness_ = hybridness_min;
	for (const Neighbor& hybridParent1 : hybridParents)
	{
		const Leaf* parent1_ = hybridParent1. leaf;
		ASSERT (parent1_ != this);
		for (const uint objNum : parent1_->pathObjNums)
		{
			const Dissim& dissim = dissims [objNum];
			if (! dissim. mult)
				continue;
			if (dissim. target <= 0)
				continue;
      if (dissim. prediction >= dissim. target) 
      	continue;
  		const Leaf* parent2_ = dissim. getOtherLeaf (parent1_);
			if (parent2_ == this)
			  continue;
    	const size_t i = hybridParents. binSearch (Neighbor (parent2_));
    	if (i == NO_INDEX)
    		continue;
    	if (maximize (hybridness_, dissim. target / (hybridParent1. target + hybridParents [i]. target)))
      {
    		parent1 = parent1_;
    		parent2 = parent2_;
    		parent1_dissim = hybridParent1. target;
    		parent2_dissim = hybridParents [i]. target;
    	}
		}
	}
	ASSERT ((bool) parent1 == (bool) parent2);
	ASSERT ((bool) parent1 == ! isNan (parent1_dissim));
	ASSERT ((bool) parent2 == ! isNan (parent2_dissim));
	
	if (parent1)
	{
		ASSERT (hybridness_ > hybridness_min);
		return hybridness_;
	}
	

#if 0
  // Compute if dissim_requests are needed ??
  Tree::LcaBuffer buf;
	for (const Neighbor& hybridParent1 : hybridParents)
  	for (const Neighbor& hybridParent2 : hybridParents)
  	{
  		if (& hybridParent1 == & hybridParent2)
  			break;
	    const TreeNode* lca_ = nullptr;
	    const VectorPtr<TreeNode>& path = Tree::getPath ( hybridParent1. leaf
																								      , hybridParent2. leaf
																		                  , nullptr
																		                  , lca_
																		                  , buf
																		                  );
	    ASSERT (! path. empty ());
	    const Real hybridness_tree = DistTree::path2prediction (path) / (hybridParent1. target + hybridParent2. target);
	    if (maximize (hybridness_, hybridness_tree))
	    {
	    	parent1 = hybridParent1. leaf;
	    	parent2 = hybridParent2. leaf;
	    }
    }    
	ASSERT ((bool) parent1 == (bool) parent2);
	ASSERT (isNan (parentDissimilarity));  
#endif

	
  return NaN;
}



void Leaf::getHybridTriangles (Vector<Triangle> &hybrids,
				                       Set<const Leaf*> &tried) const
{
	ASSERT (hybridness_min > 1);
	
	if (tried. contains (this))
		return;
		
	tried << this;
		
	const Leaf* parent1 = nullptr;
  const Leaf* parent2 = nullptr;
  Real parent1_dissim = NaN;
  Real parent2_dissim = NaN;
	const Real hybridness_ = getHybridness (parent1, parent2, parent1_dissim, parent2_dissim);
	if (isNan (hybridness_))
		return;
	if (hybridness_ < hybridness_min)
		return;
				
	ASSERT (parent1);
	ASSERT (parent2);
	
	hybrids << Triangle (this, hybridness_, parent1, parent2, parent1_dissim, parent2_dissim);
	
	parent1->getHybridTriangles (hybrids, tried);
	parent2->getHybridTriangles (hybrids, tried);	
}




// SubPath

void SubPath::qc () const
{
  if (! qc_on)
  	return;
//ASSERT (objNum != NO_INDEX);
  
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
{}



void Subgraph::qc () const 
{
  if (! qc_on)
    return;
    
  if (empty ())
    return;
    
  // area
  ASSERT (area. searchSorted);
  for (const Tree::TreeNode* node : area)
  {
    ASSERT (node);
    ASSERT (node->graph == & tree);
    ASSERT (! static_cast <const DTNode*> (node) -> inDiscernible ());
  }
    
  // boundary
  ASSERT (boundary. size () >= 2);
  ASSERT (boundary. searchSorted);
  ASSERT (area. containsFastAll (boundary));
  ASSERT ((area. size () == 2) == (area. size () == boundary. size ()));
  for (const Tree::TreeNode* node : boundary)
  	ASSERT (! boundary. containsFast (node->getParent ()));

  ASSERT ((bool) area_root == (bool) area_underRoot);

  if (area_root)
  {
    // area_root
  //ASSERT (area_root->asSteiner ());
    ASSERT (area_root->graph == & tree);
    ASSERT (area. containsFast (area_root));
    ASSERT (! area. containsFast (area_root->getParent ()));
    ASSERT (boundary. containsFast (area_root));  
    // area_underRoot
    ASSERT (area_underRoot->getParent () == area_root);
    ASSERT (area. containsFast (area_underRoot));
  }

  // subPaths
  ASSERT (! subPaths. empty ());
  {
    Vector<uint> objNums;  objNums. reserve (subPaths. size ());
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



void Subgraph::reserve (size_t radius)
{
  ASSERT (radius);
  
 	const size_t boundary_size = 2 * (size_t) powInt (2, (uint) radius);  // PAR
  area.     reserve (boundary_size * 2); 
  boundary. reserve (boundary_size);
  
  const Real path_size_max = 2 * log2 ((Real) tree. name2leaf. size ());
  const Prob area_prob = min (1.0, (Real) boundary_size / (Real) tree. name2leaf. size ());
  const Prob p = 1 - pow (1 - area_prob, path_size_max);
  ASSERT (isProb (p));
  subPaths. reserve ((size_t) ((Real) tree. dissims. size () * p)); 
}



void Subgraph::removeIndiscernibles ()
{
  ASSERT (! area. empty ());
  
  for (Iter<VectorPtr<Tree::TreeNode>> iter (area); iter. next (); )
    if (static_cast <const DTNode*> (*iter) -> inDiscernible ())
      iter. erase ();

  Set<const Tree::TreeNode*> newBoundary; 
  for (Iter<VectorPtr<Tree::TreeNode>> iter (boundary); iter. next ();)
  {
    const DTNode* node = static_cast <const DTNode*> (*iter);
    if (node->inDiscernible ())
    {
      iter. erase ();
      newBoundary << node->getParent ();
    }
  }
  if (qc_on)
    for (const Tree::TreeNode* node : boundary)
      ASSERT (! newBoundary. contains (node));
  for (const Tree::TreeNode* node : newBoundary)
    boundary << node;    
}



void Subgraph::finish ()
{
	ASSERT (area. size () >= 2);
	
  area. sort ();
  boundary. sort ();
  
  // area_root
  ASSERT (! area_root);
  for (const Tree::TreeNode* node : area)
  {
    ASSERT (node);
    if (! area. containsFast (node->getParent ()))
    {
      ASSERT (! area_root);
      EXEC_ASSERT (area_root = static_cast <const DTNode*> (node) -> asSteiner ());
    }
  }
  ASSERT (area_root);
  if (! boundary. containsFast (area_root))
    area_root = nullptr;
  
  // area_underRoot
  ASSERT (! area_underRoot);
  if (area_root)
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



void Subgraph::finishSubPaths ()
{
  ASSERT (subPathsAbsCriterion == 0);
  

  subPaths. sort ();
  subPaths. uniq ();


  // SubPath::{node1,node2}
  for (const Tree::TreeNode* node : boundary)
  {
    const DTNode* dtNode = static_cast <const DTNode*> (node);
    const Vector<uint>& pathObjNums = boundary2pathObjNums (dtNode);
    for (const uint objNum : pathObjNums)
    {
      if (! tree. dissims [objNum]. valid ())
        continue;
      const SubPath subPath_ (objNum);
      const size_t index = subPaths. binSearch (subPath_);
      if (index == NO_INDEX)  // impossible for optimizeSmallSubgraph() ??
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
  Tree::LcaBuffer buf;
  for (SubPath& subPath : subPaths)
  {
    const Dissim& dissim = tree. dissims [subPath. objNum];

    subPathsAbsCriterion += dissim. getAbsCriterion ();

    const VectorPtr<Tree::TreeNode>& path = getPath (subPath, buf);
    ASSERT (! path. empty ());
    const Real dist_hat_sub = DistTree::path2prediction (path);

    ASSERT (isNan (subPath. dist_hat_tails));
    subPath. dist_hat_tails = max (0.0, dissim. prediction - dist_hat_sub);
    ASSERT (! isNan (subPath. dist_hat_tails));

    subPath. qc ();
  }
}



Real Subgraph::getImprovement (const DiGraph::Node2Node &boundary2new) const
{
  Real s = 0;
  Tree::LcaBuffer buf;
  for (const SubPath& subPath : subPaths)
  {
    const Dissim& dissim = tree. dissims [subPath. objNum];
    if (! dissim. mult)
      continue;      
    const Tree::TreeNode* lca_ = nullptr;
    const VectorPtr<Tree::TreeNode>& path = Tree::getPath ( static_cast <const Tree::TreeNode*> (findPtr (boundary2new, subPath. node1))
																				                  , static_cast <const Tree::TreeNode*> (findPtr (boundary2new, subPath. node2))
																				                  , nullptr
																				                  , lca_
																				                  , buf
																				                  );
    ASSERT (! path. empty ());
    const Real dist_hat_sub = DistTree::path2prediction (path);
    s += dissim. mult * sqr ((dissim. target - subPath. dist_hat_tails) - dist_hat_sub);
  }
  
  return subPathsAbsCriterion - s;
}



void Subgraph::subPaths2tree ()
{
  chron_subgraph2tree. start ();  

  DistTree& tree_ = const_cast <DistTree&> (tree);

  // Time: O(p + |subPaths|)
  Vector<bool> subPathDissims (tree. dissims. size (), false);  
  for (const SubPath& subPath : subPaths)
    subPathDissims [subPath. objNum] = true;
  
  // Time: O(|area| (log(|boundary|) + p/n log(n)))  
  for (const Tree::TreeNode* node : area)
    if (! boundary. containsFast (node))
      const_static_cast <DTNode*> (node) -> pathObjNums. filterValue 
      //([this] (uint objNum) { const SubPath subPath (objNum); return subPaths. containsFast (subPath); });
          // Faster if n >> 0 <= O(p) above is not needed
        ([&subPathDissims] (uint objNum) { return subPathDissims [objNum]; });
  
  tree_. absCriterion -= subPathsAbsCriterion;
  Tree::LcaBuffer buf;
  // Time: O(|subPaths| log(|area|) log(|boundary|))
  for (const SubPath& subPath : subPaths)
  {
    subPath. qc ();
    ASSERT (subPath. node1->graph == & tree);
    ASSERT (subPath. node2->graph == & tree);
    const uint objNum = subPath. objNum;
    const Tree::TreeNode* lca_ = nullptr;
    const VectorPtr<Tree::TreeNode>& path = Tree::getPath (subPath. node1, subPath. node2, area_root, lca_, buf);
    ASSERT (lca_);
    Dissim& dissim = tree_. dissims [objNum];
    if (! viaRoot (subPath))
    {
      const Steiner* lca = static_cast <const DTNode*> (lca_) -> asSteiner ();
      ASSERT (lca);
      dissim. lca = lca;
    }
    for (const Tree::TreeNode* node : path)
      if (! boundary. containsFast (node))
      {
        const Steiner* st = static_cast <const DTNode*> (node) -> asSteiner ();
        ASSERT (st);
        if (qc_on && verbose ())
        {
        //ASSERT (area. containsFast (st));  
          ASSERT (! st->pathObjNums. contains (objNum));  
        }
        const_cast <Steiner*> (st) -> pathObjNums << objNum;
      }
    dissim. prediction = subPath. dist_hat_tails + DistTree::path2prediction (path);
    tree_. absCriterion += dissim. getAbsCriterion ();
  }
  ASSERT (DM_sp::finite (tree. absCriterion));
  maximize (tree_. absCriterion, 0.0);

  chron_subgraph2tree. stop ();
  
  tree_. qcPaths (); 
}



void Subgraph::node2subPaths (const DTNode* node)
{ 
  ASSERT (node);
  ASSERT (node->graph == & tree);
  for (const uint objNum : node->pathObjNums)
    if (tree. dissims [objNum]. valid ())
      subPaths << SubPath (objNum);
}



void Subgraph::area2subPaths ()
// Approximation: parameter sparse: Use O(boundary.size()) dissimilarities, use getReprLeaf()'s like in sparsing ??
{
  ASSERT (subPaths. empty ());
  bool first = true;
  for (const Tree::TreeNode* node : boundary)
  {
    const DTNode* dtNode = static_cast <const DTNode*> (node);
    if (dtNode == area_root)
      dtNode = area_underRoot;
    if (first)
      first = false;  // One of dtNode's is redundant given the others
    else
      node2subPaths (dtNode);
  }
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
	subgraph. qc ();
}



bool Change::apply ()
{
  ASSERT (status == eInit);
  ASSERT (targets. size () == 2);
  ASSERT (! oldParent);
  ASSERT (! arcEnd);
  ASSERT (! inter);
  

  status = eApplied;
  const bool tried = ! isNan (improvement);
  improvement = 0;

  
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

  fromLen = from->len;
  toLen = to->len;
  

  // subgraph
  ASSERT (subgraph. empty ());
  {
    subgraph. reserve (areaRadius_std);
    VectorPtr<Tree::TreeNode>& area     = subgraph. area;
    VectorPtr<Tree::TreeNode>& boundary = subgraph. boundary;
    const Tree::TreeNode* lca_ = nullptr;
	  Tree::LcaBuffer buf;
    area = Tree::getPath (from, to, tried ? nullptr : from->getAncestor (areaRadius_std), lca_, buf);  
    ASSERT (area. size () >= 1);
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
      area. pop_back ();
    area << lca;
    area. sort ();
    ASSERT (area. isUniq ());  // Interior area
    for (const Tree::TreeNode* node : area)
    {
      const VectorPtr<DiGraph::Node> children (node->getChildren ());
      for (const DiGraph::Node* child : children)
        if (! area. containsFast (static_cast <const Tree::TreeNode*> (child)))
          boundary << static_cast <const Tree::TreeNode*> (child);
    }
    if (const Tree::TreeNode* lcaParent = lca->getParent ())
      boundary << lcaParent;
    boundary. sort ();
    ASSERT (boundary. isUniq ());
    ASSERT (! area. intersectsFast_merge (boundary));
    area << boundary;  // Real area
    subgraph. finish ();
    ASSERT (subgraph. boundary. containsFast (from));
    ASSERT (subgraph. area. containsFast (to));
    IMPLY (arcEnd, subgraph. area. containsFast (arcEnd));

    subgraph. node2subPaths (from);
    subgraph. node2subPaths (to);
    subgraph. finishSubPaths ();

    subgraph. qc ();
    if (qc_on)
      for (const SubPath& subPath : subgraph. subPaths)
        ASSERT (   subPath. contains (from)
                || subPath. contains (to)
                || subPath. contains (subgraph. area_root)
               );
  }
  

  // Topology
  inter = new Steiner (const_cast <DistTree&> (tree), arcEnd, 0);
  const_cast <DTNode*> (from) -> setParent (inter); 
  const_cast <DTNode*> (to)   -> setParent (inter);
  
  // DTNode::len
  const_cast <DTNode*> (from) -> len = 0;
  const_cast <DTNode*> (to)   -> len = 0;
  ASSERT (inter->len == 0);
  

  Dataset ds;
  ds. objs. reserve (subgraph. subPaths. size ());
  auto target = new RealAttr1 ("target", ds);
  FFOR (size_t, objNum, subgraph. subPaths. size ())
    ds. appendObj ();
  auto fromAttr  = new ExtBoolAttr1 ("from",  ds);
  auto toAttr    = new ExtBoolAttr1 ("to",    ds);
  auto interAttr = new ExtBoolAttr1 ("inter", ds);
  fromAttr ->setAll (EFALSE);
  toAttr   ->setAll (EFALSE);
  interAttr->setAll (EFALSE);
  Tree::LcaBuffer buf;
  FFOR (size_t, objNum, subgraph. subPaths. size ())
  {
    const SubPath& subPath = subgraph. subPaths [objNum];
    const VectorPtr<Tree::TreeNode>& path = subgraph. getPath (subPath, buf);
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
    const Dissim& dissim = tree. dissims [subPath. objNum];
    const_cast <Obj*> (ds. objs [objNum]) -> mult = dissim. mult; 
    (*target) [objNum] = dissim. target - subPath. dist_hat_tails - DistTree::path2prediction (path);
  }
  ds. qc ();
  
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
  if (isNan (lr. absCriterion))
    return false;

  FFOR (size_t, i, lr. beta. size ())
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
    inter                       -> len = NaN;
    ASSERT (inter == tree. root);
  }

  Real subPathsAbsCriterion = 0;
  FFOR (size_t, objNum, subgraph. subPaths. size ())
  {
    const SubPath& subPath = subgraph. subPaths [objNum];
    const Dissim& dissim = tree. dissims [subPath. objNum];
    const Real prediction = nonNegative (dissim. target - lr. getResidual (objNum));
    subPathsAbsCriterion += dissim. getAbsCriterion (prediction);
  }    
  improvement = max (0.0, subgraph. subPathsAbsCriterion - subPathsAbsCriterion);


  return true;
}



void Change::restore ()
{
  ASSERT (status == eApplied);
  ASSERT (oldParent);
//ASSERT (arcEnd);
  ASSERT (inter);

  status = eInit;

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
  
  subgraph. clear ();
}



void Change::commit ()
{
  ASSERT (status == eApplied);
  ASSERT (oldParent);
  ASSERT (inter);

	status = eDone;

  if (arcEnd)
    inter->pathObjNums = to->pathObjNums;
  else
  {
    ASSERT (to->pathObjNums. empty ());
    const_cast <DTNode*> (to) -> pathObjNums = from->pathObjNums;
  }

  subgraph. area << inter;
  subgraph. area. sort ();  
  subgraph. subPaths2tree ();
  
  if (oldParent->isTransient ())
    const_cast <DistTree&> (tree). delayDeleteRetainArcs (oldParent);

  subgraph. clear ();
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
  if (a < b)  return true;  // non-stable result ??
    
  return false;
}




// Dissim

Dissim::Dissim (const Leaf* leaf1_arg,
                const Leaf* leaf2_arg,
                Real target_arg,
                Real mult_arg)
: leaf1 (leaf1_arg)
, leaf2 (leaf2_arg)
, target (target_arg)
, mult (mult_arg)
{
  ASSERT (leaf1);
  ASSERT (leaf2);
  ASSERT (leaf1->graph);
  ASSERT (leaf1->graph == leaf2->graph);
  if (leaf1->name > leaf2->name)
    swap (leaf1, leaf2);
  ASSERT (leaf1->name < leaf2->name);
  ASSERT (mult >= 0);
}



void Dissim::qc () const
{
	if (! qc_on)
		return;
	
  ASSERT (leaf1);
  ASSERT (leaf2);
  ASSERT (leaf1 != leaf2);
  ASSERT (leaf1->name < leaf2->name);
  
  ASSERT (mult >= 0);
  ASSERT (prediction >= 0);
        
  if (! mult)
  	return;

  ASSERT (valid ());
  ASSERT (! isNan (target));
  ASSERT (& leaf1->getDistTree () == & leaf2->getDistTree ());
  IMPLY (! leaf1->getDistTree (). subDepth && target == 0, ! leaf1->discernible && ! leaf2->discernible && leaf1->getParent () == leaf2->getParent ());
  ASSERT (lca);
}



string Dissim::getObjName () const
{ 
  return DistTree::getObjName (leaf1->name, leaf2->name); 
}



VectorPtr<Tree::TreeNode>& Dissim::getPath (Tree::LcaBuffer &buf) const  
{ 
  ASSERT (lca);
  const Tree::TreeNode* lca_ = nullptr;
  VectorPtr<Tree::TreeNode>& path = Tree::getPath (leaf1, leaf2, lca, lca_, buf);
  ASSERT (lca_ == lca);
  return path;
}



Real Dissim::getAbsCriterion (Real prediction_arg) const
{
  ASSERT (mult >= 0);
  if (mult == 0)
    return 0;
  if (! (prediction_arg >= 0))
  {
    cout << prediction_arg << ' ' << prediction << ' ' << mult << ' ' << target << endl;
    ERROR;
  }
  return mult * sqr (prediction_arg - target);
}



Real Dissim::setPathObjNums (size_t objNum,
	                           Tree::LcaBuffer &buf)
{
  ASSERT (valid ());

  const VectorPtr<Tree::TreeNode>& path = getPath (buf);
  prediction = DistTree::path2prediction (path);  
  for (const Tree::TreeNode* node : path)
    if (const Steiner* st = static_cast <const DTNode*> (node) -> asSteiner ())
      const_cast <Steiner*> (st) -> pathObjNums << (uint) objNum;  
  qc ();

  return getAbsCriterion ();
}
	
	


// Image

Image::Image (const DistTree &mainTree)
: subgraph (mainTree)
{ 
	ASSERT (mainTree. optimizable ());
	subgraph. reserve (areaRadius_std); 
}
	        


Image::~Image ()
{ 
	delete tree; 
}



void Image::processSmall (const DTNode* center_arg,
			                    uint areaRadius)
{ 
  ASSERT (subgraph. empty ());
	ASSERT (! tree);    	   	  	
  ASSERT (! center);

  
	center = center_arg;	
  ASSERT (center);
  ASSERT (center->graph);
  ASSERT (! center->inDiscernible ());

  ASSERT (areaRadius >= 1);
//ASSERT (areaRadius <= areaRadius_std);	  

  const DistTree& wholeTree = center->getDistTree ();
  ASSERT (& wholeTree == & subgraph. tree);


  chron_tree2subgraph. start ();

  // subgraph: area, boundary
  // areaRadius
  VectorPtr<Tree::TreeNode>& area     = subgraph. area;
  VectorPtr<Tree::TreeNode>& boundary = subgraph. boundary;
  for (;;)
  {
    center->getArea (areaRadius, area, boundary);
    subgraph. removeIndiscernibles ();
    if (areaRadius == 1 || boundary. size () <= boundary_size_max_std)
      break;
    area.     clear ();
    boundary. clear ();    
    ASSERT (areaRadius > 1);
    areaRadius--;
  }

  tree = new DistTree (subgraph, new2old);
  tree->qc ();
  rootInArea = (! subgraph. area_root || subgraph. area_root == wholeTree. root);

  ASSERT (area. containsFast (center));
  ASSERT (boundary. size () <= boundary_size_max_std);
#ifndef NDEBUG
  for (const auto& it : new2old)
  {
    // New
    ASSERT (it. first);
    ASSERT (it. first->graph == tree);
    ASSERT (static_cast <const DTNode*> (it. first) -> asLeaf ());
    ASSERT (! static_cast <const DTNode*> (it. first) -> inDiscernible ());
    // Old = boundary  
    ASSERT (it. second);
    ASSERT (it. second->graph == & wholeTree);
    ASSERT (! static_cast <const DTNode*> (it. second) -> inDiscernible ());
  }
#endif

  chron_tree2subgraph. stop ();
  

  if (verbose ())
    cerr << "  " << areaRadius 
         << ' ' << boundary. size () 
         << ' ' << area. size () - boundary. size ();
  

  // Optimization
  chron_subgraphOptimize. start ();
  const size_t leaves = tree->name2leaf. size ();
  ASSERT (leaves >= 2);
  {
    Unverbose unv;
    if (leaves > 3)
    {
      if (subgraph. unresolved ())
      {
      //cerr << " NJ";  
        // A bifurcating tree is needed for speed
        tree->neighborJoin ();
        neighborJoinP = true;
      }
      tree->optimizeLenArc ();
      tree->optimizeLenNode ();  
      if (subgraph. large () && areaRadius > 1)
        tree->optimizeSmallSubgraphs (areaRadius / 2, false); 
      else
        tree->optimizeWholeIter (20, string ());  // PAR
    }
    else if (leaves == 3)
    {
      tree->optimize3 ();
      neighborJoinP = true;
    }
    else if (leaves == 2)
      tree->optimize2 ();
  }
  tree->qc ();
  if (verbose ())
  {
    cout << "Subtree: ";
    tree->reportErrors (cout);
  }
  chron_subgraphOptimize. stop (); 
}



void Image::processLarge (const Steiner* subTreeRoot,
	                        const VectorPtr<Tree::TreeNode> &possibleBoundary)
{
  ASSERT (subgraph. empty ());
	ASSERT (! tree);    	   	  	
  ASSERT (! center);
#ifndef NDEBUG
  for (const Tree::TreeNode* node : possibleBoundary)
  {
  	ASSERT (node->graph == & subgraph. tree);
  	ASSERT (node->isTransient ());
  }
#endif

  const DistTree& wholeTree = subgraph. tree;

  // subgraph: area, boundary
  VectorPtr<Tree::TreeNode>& area     = subgraph. area;
  VectorPtr<Tree::TreeNode>& boundary = subgraph. boundary;
	area.     reserve (2 * wholeTree. name2leaf. size ());
  boundary. reserve (    wholeTree. name2leaf. size ());
  if (subTreeRoot)
  {
	  ASSERT (& subTreeRoot->getDistTree () == & wholeTree);
	  ASSERT (subTreeRoot->isTransient ());
    subTreeRoot->getSubtreeArea (possibleBoundary, area, boundary);
    boundary << subTreeRoot;
  }
  else
    static_cast <const Steiner*> (wholeTree. root) -> getSubtreeArea (possibleBoundary, area, boundary);
  subgraph. removeIndiscernibles ();
//cout << "Subgraph size: " << subgraph. area. size () << endl;
  
  try
  {
    tree = new DistTree (subgraph, new2old);
    tree->qc ();
    rootInArea = (! subgraph. area_root || subgraph. area_root == wholeTree. root);

    // Optimization
  	tree->optimizeSmallSubgraphs (areaRadius_std, false);
  	
  	tree->qc ();
  }
  catch (const bad_alloc &)  
  { 
    delete tree;
    tree = nullptr;
    subgraph. clear ();
  }
}



bool Image::apply ()
{
	if (! tree)
	  return false;
	  
	subgraph. qc ();
	

  DistTree& wholeTree = const_cast <DistTree&> (subgraph. tree);
	

  DiGraph::Node2Node boundary2new (DiGraph::reverse (new2old));
  ASSERT (boundary2new. size () == new2old. size ());

  if (   neighborJoinP   // Subgraph::getImprovement() is slow
      && subgraph. getImprovement (boundary2new) <= 1e-5  // PAR
     )
  {
    for (const Tree::TreeNode* node : subgraph. area)
      if (! subgraph. boundary. containsFast (node))
        const_static_cast <DTNode*> (node) -> stable = true;
  //cerr << " no improvement" << endl;  
  }
  else
  {
	  // Optimization of tree is not needed any more
	  for (const DiGraph::Node* node_new : tree->nodes)
	    const_static_cast <DTNode*> (node_new) -> pathObjNums. clear ();  
	
	  const Tree::TreeNode* root_old = wholeTree. root;
	
	  // tree->root
	  // new2old: newBoundary2oldBoundary
	  if (const DTNode* dtNode = static_cast <const DTNode*> (findPtr (boundary2new, subgraph. area_root)))
	  {
	  //ASSERT (! rootInArea);
	    ASSERT (subgraph. area_root);
	    const Leaf* leaf = dtNode->asLeaf ();
	    ASSERT (leaf);
	    if (verbose ())
	      cout << "Re-rooting subgraph ..." << endl;
	    // Replace leaf by st
	    auto st = new Steiner (*tree, const_static_cast <Steiner*> (leaf->getParent ()), leaf->len);
	    EXEC_ASSERT (new2old. erase (leaf) == 1);
	    tree->delayDeleteRetainArcs (const_cast <Leaf*> (leaf));
	    ASSERT (isNan (static_cast <const DTNode*> (tree->root) -> len));
	    const Steiner* root_ = st->makeDTRoot ();  // Old root
	    boundary2new [subgraph. area_root] = st;
	    new2old [st] = subgraph. area_root;
	    ASSERT (root_ != st);
	    ASSERT (root_->len >= 0);
	    if (root_->isTransient ())
	      tree->delayDeleteRetainArcs (const_cast <Steiner*> (root_));
	    tree->toDelete. deleteData ();
	    ASSERT (st == tree->root);
	    ASSERT (isNan (st->len));
	  }
	  ASSERT (new2old. size () == boundary2new. size ());
	  ASSERT (subgraph. boundary. size () == boundary2new. size ());
	//tree->qc ();
	
	  // tree, subgraph --> topology of wholeTree
	
	  // nodes: delete
	  // center may be delete'd
	#ifndef NDEBUG
	  size_t deleted = 0;
	#endif
	  for (const Tree::TreeNode* node : subgraph. area)
	    if (! contains (boundary2new, node))
	    {
	      DTNode* dtNode = const_static_cast <DTNode*> (node);
	      dtNode->isolate ();
	      dtNode->detach ();
	      wholeTree. toDelete << dtNode;
		  #ifndef NDEBUG
	      deleted++;
	    #endif
	    }
	    // Between boundary Tree::TreeNode's there are no Arc's => borrowArcs() will not break
	
	  // nodes, new2old[]: new
	#ifndef NDEBUG
	  size_t created = 0;
	#endif
	  for (const DiGraph::Node* node_new : tree->nodes)
	  {
	    DTNode* node = const_static_cast <DTNode*> (findPtr (new2old, node_new));
	    if (! node)
	    {
	      node = new Steiner (wholeTree, nullptr, NaN);
	      new2old [node_new] = node;
	      node->stable = true;
			#ifndef NDEBUG
	      created++;
	    #endif
	    }
	    ASSERT (node);
	    if (node_new != tree->root)
	      node->len = static_cast <const DTNode*> (node_new) -> len;
	  }
	  ASSERT (new2old. size () == tree->nodes. size ());
	  
	  ASSERT ((bool) deleted == (bool) created);
	  ASSERT ((bool) deleted == (subgraph. area. size () > 2));

	  if (subgraph. area. size () > 2)
	    wholeTree. borrowArcs (new2old, true);  // Arc's are not parallel
	  
	  // wholeTree.root
	  if (rootInArea)
	  {
	    wholeTree. root = nullptr;
	    for (const auto& it : new2old)
	    {
	      const Tree::TreeNode* node = static_cast <const Tree::TreeNode*> (it. second);
	      if (! node->getParent ())
	      {
	        ASSERT (! wholeTree. root);
	        wholeTree. root = node;
	      }
	    }
	  }
	  else
	  {
	    ASSERT (root_old);
	    ASSERT (! root_old->getParent ());
	    wholeTree. root = root_old;
	  }
	  ASSERT (wholeTree. root);
	   
	  // Topology, absCriterion
	  subgraph. subPaths2tree ();
	    
	  // deleteLenZero
	  for (const auto& it : new2old)
	  {
	    DTNode* node = const_static_cast <DTNode*> (it. second);
	    if (! subgraph. boundary. containsFast (node))
	      wholeTree. deleteLenZero (node);  
	  }

    wholeTree. toDelete. deleteData ();
	}


  if (center && contains (boundary2new, center))
  {
    ASSERT (center->graph == & wholeTree);
    const_cast <DTNode*> (center) -> stable = true;
  }


  if (verbose (1))
    wholeTree. reportErrors (cout);

    
  return true;
}




// DistTree

DistTree::DistTree (const string &treeFName,
	                  const string &dissimFName,
	                  const string &dissimAttrName,
	                  bool sparse)
: rand (seed_global)
{
	ASSERT (dissimFName. empty () == dissimAttrName. empty ());

  // Initial tree topology
  loadTreeFile (treeFName);  
  ASSERT (root);
  ASSERT (nodes. front () == root);
  ASSERT (static_cast <const DTNode*> (root) -> asSteiner ());

  setName2leaf ();
        
  if (! dissimFName. empty ())
  {
    loadDissimDs (dissimFName, dissimAttrName);
    
    if (! getConnected ())
      throw runtime_error ("Disconnected objects");
  	setDiscernibles_ds ();  

    dissimDs2dissims (sparse);    	
  }
}



DistTree::DistTree (const string &dirName,
	                  const string &dissimFName,
	                  const string &dissimAttrName)
: rand (seed_global)
{
  // Initial tree topology, no DTNode::len
  loadTreeDir (dirName);
  ASSERT (root);
  ASSERT (nodes. front () == root);
  ASSERT (static_cast <const DTNode*> (root) -> asSteiner ());

  setName2leaf ();
        
  loadDissimDs (dissimFName, dissimAttrName);
  
  if (! getConnected ())
    throw runtime_error ("Disconnected objects");
	setDiscernibles_ds ();  

  setGlobalLen ();  // --> after dissim2Ds(), use dissims ??

  dissimDs2dissims (false);  
}



DistTree::DistTree (const string &dissimFName,
	                  const string &dissimAttrName,
	                  bool sparse)
: rand (seed_global)
{
  loadDissimDs (dissimFName, dissimAttrName);

  // Initial tree topology: star topology 
  ASSERT (! root);
  auto root_ = new Steiner (*this, nullptr, NaN);
  for (const Obj* obj : dissimDs->objs)
	  new Leaf (*this, root_, 0, obj->name);
  ASSERT (root);
  ASSERT (nodes. front () == root);
  ASSERT (static_cast <const DTNode*> (root) -> asSteiner ());

  setName2leaf ();
        
  if (! getConnected ())
    throw runtime_error ("Disconnected objects");
	setDiscernibles_ds ();

  neighborJoin ();
  
  dissimDs2dissims (sparse);    	 
}



struct DissimLine
{
	// Input
	string name1;
	string name2;
	Real dissim {NaN};
	// Output
  Leaf* leaf1 {nullptr};
	Leaf* leaf2 {nullptr};
	

  DissimLine ()
    {}
	DissimLine (string &line,
	            uint lineNum)
		{ replace (line, '\t', ' ');
      name1 = findSplit (line);
      name2 = findSplit (line);
      const string errorS ("Line " + toString (lineNum) + ": ");
      if (name2. empty ())
        throw runtime_error (errorS + "empty name2");
      if (name1 >= name2)
        throw runtime_error (errorS + "name1 >= name2");
    #if 0
      if (outliers. containsFast (name1))
        continue;
      if (outliers. containsFast (name2))
        continue;
    #endif
      dissim = str2real (line);
      if (isNan (dissim))
        throw runtime_error (errorS + "dissimilarity is NaN");
      if (dissim < 0)
        throw runtime_error (errorS + "dissimilarity is negative");
      if (! DM_sp::finite (dissim))
        throw runtime_error (errorS + "dissimilarity is infinite");
		}
		

	void process (const map<string/*Leaf::name*/,const Leaf*> &name2leaf)
	  {
			ASSERT (! leaf1);
			ASSERT (! leaf2);
			
		  leaf1 = const_cast <Leaf*> (findPtr (name2leaf, name1));
		  if (! leaf1)
		    return;
		  //throw runtime_error ("Tree has no object " + name1);
		  leaf1->paths = 1;  // temporary
		
		  leaf2 = const_cast <Leaf*> (findPtr (name2leaf, name2));
		  if (! leaf2)
		    return;
		  //throw runtime_error ("Tree has no object " + name2);
		  leaf2->paths = 1;  // temporary
	  }
	  
	  
	void apply (DistTree &tree) const
    { ASSERT (! isNan (dissim));
      if (! leaf1)
        return;
      if (! leaf2)
        return;
      if (dissim == 0 && ! leaf1->getCollapsed (leaf2))  // Only for new Leaf's
        leaf1->collapse (leaf2);
      if (! tree. addDissim (leaf1, leaf2, dissim))
      	throw runtime_error ("Cannot add dissimilarity: " + name1 + " " + name2 + " " + toString (dissim));
    }


  bool operator< (const DissimLine &other) const
    { LESS_PART (*this, other, name1);
    	LESS_PART (*this, other, name2);
    	return false;
    }
  bool operator== (const DissimLine &other) const
    { return    name1 == other. name1
    	       && name2 == other. name2;
    }
};



namespace
{

void processDissimLine (size_t from,
                        size_t to,
                        Notype /*&res*/,
                        Vector<DissimLine> &dls,
                        const map<string/*Leaf::name*/,const Leaf*> &name2leaf)
{
  Progress prog (to - from, 1000);  // PAR
	FOR_START (size_t, i, from, to)
	{
	  prog ();
		dls [i]. process (name2leaf);
	}
}



struct OptimizeSmallSubgraph
{
  Image image;
	const DTNode* center;
	const uint radius;
	
	OptimizeSmallSubgraph (DistTree &tree,
	                       const DTNode* center_arg,
	                       uint radius_arg)
	  : image (tree)
	  , center (center_arg)
	  , radius (radius_arg)
	  { ASSERT (center);
	  	ASSERT (radius);
	  	ASSERT (center->graph == & image. subgraph. tree);
	  }
	  
	bool close (const OptimizeSmallSubgraph &oss) const
	  { Tree::LcaBuffer buf; 
	  	const Tree::TreeNode* lca = nullptr;
   	  const size_t arcDist = Tree::getPath (center, oss. center, nullptr, lca, buf). size ();
   	  return arcDist <= radius + oss. radius + 1;  // ??
   	}
  void process ()
    { image. processSmall (center, radius); }
  bool apply ()
    { Unverbose unv;
    	return image. apply (); 
    }
};



void processOptimizeSmallSubgraph (OptimizeSmallSubgraph &oss)
{
	oss. process ();
}

}



DistTree::DistTree (const string &dataDirName,
                    bool loadNewLeaves,
 	                  bool loadDissim,
 	                  bool optimizeP)
: rand (seed_global)
{
	ASSERT (! dataDirName. empty ());
	ASSERT (dataDirName. back () == '/');
	IMPLY (optimizeP, loadDissim);
	

  // Initial tree topology
  loadTreeFile (dataDirName + "tree");  
  ASSERT (root);
  ASSERT (nodes. front () == root);
  ASSERT (static_cast <const DTNode*> (root) -> asSteiner ());
  
  setName2leaf ();  
  
  
  VectorPtr<Leaf> newLeaves;  newLeaves. reserve (name2leaf. size () / 10 + 1);  // PAR
  if (loadNewLeaves)
  {
    LineInput f (dataDirName + "leaf", 10 * 1024, 1);  // PAR
    string leafName, anchorName;
    Real leafLen, arcLen;
	  Tree::LcaBuffer buf;
    while (f. nextLine ())
    {
      istringstream iss (f. line);
      iss >> leafName >> anchorName >> leafLen >> arcLen;
      ASSERT (leafLen >= 0);
      ASSERT (arcLen >= 0);
    //ASSERT (iss. eof ());  // Extra fields
      const DTNode* anchor = lcaName2node (anchorName, buf);
      ASSERT (anchor);
      // Attach a leaf
      Leaf* leaf = nullptr;
      if (const Leaf* anchorLeaf = anchor->asLeaf ())
      {
      #if 0
        if (leafLen == 0 && arcLen == 0)
          leaf = new Leaf (*this, const_cast <Leaf*> (anchorLeaf), leafName);  // DissimLine::apply() will collapse() if needed
        else
      #endif
          if (! anchorLeaf->discernible)
            anchor = static_cast <const DTNode*> (anchor->getParent ());
      }
      ASSERT (anchor);
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
                                , anchor->len - arcLen
                                );
          ASSERT (st);
          const_cast <DTNode*> (anchor) -> setParent (st);
          const_cast <DTNode*> (anchor) -> len = max (0.0, arcLen);
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
  #if 0
    StringVector outliers;
    {
      FileItemGenerator fig (0, true, dataDirName + "outlier");  
  	  string item;
  	  while (fig. next (item))
  	    outliers << item;
    }
    outliers. sort ();
  #endif
    loadDissimPrepare (name2leaf. size () * getSparseDissims_size ()); 
    const string fName (dataDirName + "dissim");
    {
	    Vector<DissimLine> dissimLines;  dissimLines. reserve (dissims. capacity ());
	    constexpr uint displayPeriod = 100000;  // PAR
	    {
	      cout << "Loading " << fName << " ..." << endl;
	      {
		      LineInput f (fName, 10 * 1024 * 1024, displayPeriod);  
		      while (f. nextLine ())
		      	dissimLines << DissimLine (f. line, f. lineNum);
		      if (! f. lineNum)
		        throw runtime_error ("Empty " + fName);
		    }
	      dissimLines. sort ();
	      dissimLines. uniq ();
	    }
	    vector<Notype> notypes;
	    arrayThreads (processDissimLine, dissimLines. size (), notypes, ref (dissimLines), cref (name2leaf));
	    {
	    	Progress prog (dissimLines. size (), displayPeriod);
		    for (const DissimLine &dl : dissimLines)
		    {
		    	prog ();
		    	dl. apply (*this);
		    }
		  }
    }
    // qc
	  for (const auto& it : name2leaf)
	  {
	  	const Leaf* leaf = it. second;
      if (! leaf->paths)
        throw runtime_error ("No dissimilarities for object " + leaf->name);
      const_cast <Leaf*> (leaf) -> paths = 0;
	  }

    // qc: pairs of <leaf1,leaf2> must be unique in dissims
    if (qc_on)
    {
	    Vector<Pair<const Leaf*>> pairs;  pairs. reserve (dissims. size ());
	    for (const Dissim& d : dissims)
	    	pairs << Pair<const Leaf*> (d. leaf1, d. leaf2);
	    const size_t size_init = pairs. size ();
	    pairs. sort ();
	    pairs. uniq ();
	    if (pairs. size () != size_init)
	    {
	    	cout << size_init << " -> " << pairs. size () << endl;
	    	throw runtime_error (fName + ": non-unique pairs");
	    }
	  }
  } 
  
  
  cleanTopology ();
  if (getDiscernibles (). size () <= 1)
  	throw runtime_error ("No discernible objects");
  	
     
  if (loadDissim)
  {
    setPaths ();        

    if (optimizeP)
    {  
      const Chronometer_OnePass cop ("Optimizing new leaves");  
      cout << "Optimizing new leaves ..." << endl;
      constexpr uint radius = 2 * areaRadius_std;  // PAR
      if (threads_max > 1)
      {
      	{
		      Progress prog (newLeaves. size ());
		      while (! newLeaves. empty ())
		      {
			      VectorOwn<OptimizeSmallSubgraph> osss;  osss. reserve (newLeaves. size ());
			      FFOR (size_t, i, newLeaves. size ())
			      {
			      	if (osss. size () >= threads_max)
			      		break;
			      	auto oss = new OptimizeSmallSubgraph (*this, newLeaves [i] -> getDiscernible (), radius);  
			      	for (const OptimizeSmallSubgraph* oss1 : osss)
			     			if (oss->close (*oss1))
			     			{
			     				delete oss;
			     				oss = nullptr;
			     				break;
			     			}
			     		if (! oss)
			     			continue;
			     		osss << oss;
			     		newLeaves [i] = nullptr;
			      }
			      ASSERT (! osss. empty ());
						{
							Unverbose unv;
							Threads th (osss. size () - 1);
							FFOR_START (size_t, i, 1, osss. size ())						
							  th << thread (processOptimizeSmallSubgraph, ref (* const_cast <OptimizeSmallSubgraph*> (osss [i])));
			      	const_cast <OptimizeSmallSubgraph*> (osss [0]) -> process ();
						}
						newLeaves. filterValue ([] (const Leaf* leaf) { return ! leaf; });
						// Use Threads for apply() for sibling subtrees ??
		      	for (const OptimizeSmallSubgraph* oss : osss)
		      	{
		      		EXEC_ASSERT (var_cast (oss) -> apply ());
			        prog (absCriterion2str ());
		      	}
			    }
			  }
			  setPredictionAbsCriterion ();
			  reportErrors (cerr);
		  }
      else
      {
	      Progress prog (newLeaves. size ());
	      for (const Leaf* leaf : newLeaves)
	      {
	        Unverbose unv;
	        optimizeSmallSubgraph (leaf->getDiscernible (), radius);  
	          // optimizeReinsert() ??
	        prog (absCriterion2str ());
	      #ifndef NDEBUG
	       	for (const DiGraph::Node* node : nodes)
	       	{
	       	  const DTNode* dtNode = static_cast <const DTNode*> (node);
	          if (dtNode->isTransient ())
	          {
	            cout << dtNode << endl;
	            dtNode->getParent () -> saveText (cout);
	            ERROR;
	          }
	        }
	      #endif
	      }
	    }
    }
  }  


  ASSERT (! dissimDs. get ());
  ASSERT (! dissimAttr);
}



namespace 
{

struct Newick
/* Grammar:
     tree ::= subtree ;
     subtree ::= item [: length]
     item ::= interior | leaf
     interior ::= ( list ) [interior_name]
     list ::= subtree | subtree , list
     leaf ::= name
     name ::= string | ' string_space '
     interior_name ::= bootstrap [: name]
     string ::=  // No ' ', ',', ';', '(', ')'
*/
{
private:
  istream& f;
  DistTree& tree;
  char c {' '};
    // Current character of f
  streamsize pos {0};
  static const string delimiters;
public:


  struct Error : runtime_error
  {
    Error (const string &what_arg,
           const Newick &newick)
      : runtime_error (what_arg + " at position " + toString (newick. pos))
      {}
  };


  Newick (istream &f_arg,
          DistTree &tree_arg)
    : f (f_arg)
    , tree (tree_arg)
    { 
      ASSERT (f. good ()); 
      ASSERT (tree. nodes. empty ());
      parseSubtree (nullptr);
      skipSpaces ();
      if (c != ';')
        throw Error ("No ending semicolon", *this);
    }
    
    
private:
  void readChar ()
  {
    f >> c;
    pos++;
  }
  
  void skipSpaces ()
  {
    while (c == ' ')
      readChar ();
    if (f. eof ())
      throw Error ("Unexpected end of file", *this);
  }

  void parseSubtree (Steiner* parent)
  {
    DTNode* node = parseItem (parent);
    skipSpaces ();
    if (c == ':')
    {
      Real len = NaN;
      const streamsize start = f. tellg ();
      f >> len;
      pos += f. tellg () - start;
      ASSERT (! isNan (len));
      node->len = max (0.0, len);
      readChar ();
    }
  }
  
  DTNode* parseItem (Steiner* parent)
  {
    DTNode* node = nullptr;
    skipSpaces ();
    if (c == '(')
    {
      auto st = new Steiner (tree, parent, NaN);
      parseInterior (st);
      node = st;
    }
    else
      node = new Leaf (tree, parent, NaN, parseName ());
    ASSERT (node);
    return node;
  }
  
  void parseInterior (Steiner* parent)
  {
    skipSpaces ();
    if (c != '(')
      throw Error ("'(' expected", *this);
    // List
    for (;;)
    {
      readChar ();
      parseSubtree (parent);
      skipSpaces ();
      if (c == ')')
        break;
      if (c != ',')
        throw Error ("Comma expected", *this);
    }
    readChar ();
    skipSpaces ();
    // Name
    if (! contains (delimiters, c))
    {
      // bootstrap is lost ??
      string name (parseName ());
      const size_t colonPos = name. find (':');      
      if (colonPos != string::npos)
        parent->name = name. substr (colonPos + 1);
    }
  }
  
  string parseName ()
  {
    skipSpaces ();
    string name;
    if (c == '\'')
    {
      for (;;)
      {
        readChar ();
        if (c == '\'')
          break;
        name += c;
      }
      readChar ();
    }
    else
      for (;;)
      {
        if (contains (delimiters, c))
          break;
        if (! printable (c))
          throw Error ("Non-printable character in name: " + string (1, c), *this);
        name += c;
        readChar ();
      }
    if (name. empty ())
      throw Error ("Empty name", *this);
    return name;
  }
};



const string Newick::delimiters (" ,;():");
  
}  // namespace



DistTree::DistTree (const string &newickFName)
: rand (seed_global)
{ 
  {
    ifstream f (newickFName);
    const Newick newick (f, *this);
  }

	ASSERT (root);
  ASSERT (nodes. front () == root);
  if (! static_cast <const DTNode*> (root) -> asSteiner ())
    throw runtime_error ("One-node tree");
  const_static_cast <DTNode*> (root) -> len = NaN;
    
  finishChanges ();

  setName2leaf ();        
}



DistTree::DistTree (Prob branchProb,
                    size_t leafNum_max)
: rand (seed_global)
{
	ASSERT (isProb (branchProb));
	ASSERT (branchProb < 1);
	ASSERT (leafNum_max);


  const Real len = 1;  // PAR

  Set<Steiner*> open;
  open << new Steiner (*this, nullptr, NaN);
  Bernoulli bernoulli;
  bernoulli. setParam (branchProb);
  size_t leafIndex = 0;
  size_t leaves = 1;  // isLeaf()
  while (! open. empty ())
  {
    auto it = open. begin ();
    ASSERT (it != open. end ());
    Steiner* st = *it;
    ASSERT (st);
    bernoulli. randVariable ();
    if (leaves < leafNum_max && bernoulli. variable)
    {
      if (! st->isLeaf ())
        leaves++;
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



DistTree::DistTree (Subgraph &subgraph,
                    Node2Node &newLeaves2boundary)
: subDepth (subgraph. tree. subDepth + 1)
, rand (seed_global)
{
  ASSERT (! subgraph. empty ());
  ASSERT (subgraph. subPaths. empty ());  // not finish()'ed
  ASSERT (newLeaves2boundary. empty ());


  // subgraph
  subgraph. finish ();
  subgraph. area2subPaths ();    
  subgraph. finishSubPaths ();
  subgraph. qc ();

  VectorPtr<TreeNode>& area     = subgraph. area;
  VectorPtr<TreeNode>& boundary = subgraph. boundary;
  const Steiner* &area_root     = subgraph. area_root;
  const DistTree& wholeTree     = subgraph. tree;
  IMPLY (subgraph. area_root, subgraph. area_root->graph == & wholeTree);


  // nodes
  Node2Node old2new;  // 1-1
  ASSERT (nodes. empty ());
  ASSERT (leafNum == 0);
  if (subgraph. unresolved ())
  {
    // Star topology
    auto st = new Steiner (*this, nullptr, NaN);
    root = st;
    if (area_root)
    {
      auto leaf = new Leaf (*this, st, 0, "L" + toString (leafNum));
      old2new [area_root] = leaf;  
    }
    for (const TreeNode* node_old : area)
    {
      VectorPtr<DiGraph::Node> children (node_old->getChildren ());
      children. sort ();
      if (! children. intersectsFast_merge (area))
      {
        auto leaf = new Leaf (*this, st, 0, "L" + toString (leafNum));
        old2new [node_old] = leaf;
      }
    }
    ASSERT (isStar ());
  }
  else
  {
    for (const TreeNode* node_old : area)
    {
      DTNode* node_new = nullptr;
      const Real len = static_cast <const DTNode*> (node_old) -> len;
      VectorPtr<DiGraph::Node> children (node_old->getChildren ());
      children. sort ();
      if (children. intersectsFast_merge (area))
        node_new = new Steiner (*this, nullptr, len);
      else  
        node_new = new Leaf (*this, nullptr, len, "L" + toString (leafNum));
      ASSERT (node_new);
      old2new [node_old] = node_new;
    }
    ASSERT (nodes. size () == area. size ());
    ASSERT (nodes. size () == old2new. size ());
    
    borrowArcs (old2new, true);  // Arc's are not parallel
    setRoot ();
    ASSERT (static_cast <const DTNode*> (root) -> asSteiner ());
    IMPLY (area_root, findPtr (old2new, area_root) == root);
    
    // subgraph.area_root --> Leaf
    {
      const VectorPtr<DiGraph::Node> children (root->getChildren ());
      ASSERT (! children. empty ());
      if (children. size () == 1)
      {
        ASSERT (area_root);
        // old2new[area_root] --> Leaf
        DTNode* child = const_static_cast <DTNode*> (children. front ());
        ASSERT (child);
        const TreeNode* root_ = root;                    // Old root in *this
        auto inter = new Steiner (*this, nullptr, NaN);  // New root in *this
        child->setParent (inter);
        child->len /= 2;
        ASSERT (root_->getChildren (). empty ());
        delete root_;
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
        ASSERT (! area_root);
      }
    }
  }
  ASSERT (isNan (static_cast <const DTNode*> (root) -> len));  

  
  setName2leaf ();
  ASSERT (boundary. size () == name2leaf. size ());

  // newLeaves2boundary
  for (const auto& it : old2new)
    if (static_cast <const DTNode*> (it. second) -> asLeaf ())
      newLeaves2boundary [it. second] = it. first;
  ASSERT (newLeaves2boundary. size () == boundary. size ());
    

  // dissims[] 
  // For some leaf pairs the dissimilarity may be missing
  dissims. reserve (getDissimSize_max ());
 	Vector<uint/*objNum*/> leaves2objNum (leafNum * leafNum, (uint) -1);
  const size_t reserve_size = getPathObjNums_size ();
 	for (DiGraph::Node* node : nodes)
 	{
 	  const DTNode* dtNode = static_cast <DTNode*> (node);
 	  ASSERT (dtNode->pathObjNums. empty ());
 	  var_cast (dtNode) -> pathObjNums. reserve (reserve_size); 
 	}
  for (const auto it2 : name2leaf)
    for (const auto it1 : name2leaf)
    {
      if (it1. first == it2. first)
        break;
      ASSERT (it1. first < it2. first);
      const Leaf* leaf1 = it1. second;
      const Leaf* leaf2 = it2. second;
      const uint objNum = leaves2dissims (var_cast (leaf1), var_cast (leaf2), 0, 0);
      ASSERT (leaf1->index < NO_INDEX);
      ASSERT (leaf2->index < NO_INDEX);
      leaves2objNum [leaf1->index * leafNum + leaf2->index] = objNum;
    }
  ASSERT (dissims. size () == getDissimSize_max ());


  // dissims[]: mult, target
  for (const SubPath& subPath : subgraph. subPaths)
  {
    const size_t wholeObjNum = subPath. objNum;
    
    const Real dist = wholeTree. dissims [wholeObjNum]. target - subPath. dist_hat_tails;
    if (isNan (dist))
      continue;
    // May be < 0
    const Real mult = wholeTree. dissims [wholeObjNum]. mult; 
    ASSERT (mult >= 0);
    if (mult == 0)  // wholeTree.dissims[wholeObjNum].target = INF
      continue;
    
    const DTNode* node1 = static_cast <const DTNode*> (findPtr (old2new, subPath. node1));
    const DTNode* node2 = static_cast <const DTNode*> (findPtr (old2new, subPath. node2));
    ASSERT (node1);
    ASSERT (node2);
    ASSERT (node1 != node2);
    const Leaf* leaf1 = node1->asLeaf ();
    const Leaf* leaf2 = node2->asLeaf ();
    ASSERT (leaf1);
    ASSERT (leaf2);    
    if (leaf2->name < leaf1->name)
      swap (leaf2, leaf1);

    const uint objNum = leaves2objNum [leaf1->index * leafNum + leaf2->index];
    ASSERT (objNum != (uint) -1);
    
    dissims [objNum]. mult   += mult;
    dissims [objNum]. target += mult * dist;
  }
  
  
  // Dissim::target, mult_sum, dissim2_sum
  ASSERT (mult_sum == 0);
  ASSERT (dissim2_sum == 0);
  for (Dissim& dissim : dissims)
    if (dissim. mult == 0)
      dissim. target = NaN;
    else
    {
      dissim. target /= dissim. mult;
      mult_sum    += dissim. mult;
      dissim2_sum += dissim. mult * sqr (dissim. target); 
    }
    

  setPaths ();
//ASSERT (setDiscernibles_ds () == 0);
}



void DistTree::loadTreeDir (const string &dir)
{
  ASSERT (! subDepth);
	ASSERT (! dir. empty ());
  ASSERT (isDirName (dir));

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
      new Leaf (*this, steiner, NaN, obj->name);
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
  auto s = new Steiner (*this, getName2steiner (prefix, name2steiner), NaN);
  name2steiner [name] = s;
  
  return s;
}



void DistTree::loadTreeFile (const string &fName)
{
  ASSERT (! subDepth);
	ASSERT (! fName. empty ());

  StringVector lines;
  try
  {
    LineInput in (fName, 10000);  // PAR
    lines = in. getVector ();
  }
  catch (const exception &e)
  {
    throw runtime_error ("Loading tree file " + strQuote (fName) + "\n" + e. what ());
  }
  ASSERT (! lines. empty ());

	size_t lineNum = 0; 
  Steiner* steiner = nullptr;  	
	EXEC_ASSERT (loadLines (lines, lineNum, steiner, 0));
}



namespace 
{
	string token2string (const string &s,
	                     const string &token)
	{
	  const string s1 (" " + s);
		const string token1 (" " + token + "=");
		const size_t pos = s1. find (token1);
		if (pos == string::npos)
			return string ();  
			
		string valueS (s1. substr (pos + token1. size ()));
		trim (valueS);
		if (valueS. empty ())
		  throw runtime_error ("Empty token '" + token + "' in: " + s);
		  
		if (valueS [0] == '\'')
		{
  		string name;
		  size_t i = 0;
		  for (;;)
		  {
		    i++;
		    if (! valueS [i])
		      throw runtime_error ("No closing single quote");
		    if (valueS [i] == '\'')
		      break;
		    else
		      name += valueS [i];
		  }
		  return name;
		}
		
	  return findSplit (valueS);
	}


	Real token2real (const string &s,
	                 const string &token)
	{
	  const string valueS (token2string (s, token));
	  if (valueS. empty ())
	    return NaN;
		return str2real (valueS);
	}
}



bool DistTree::loadLines (const StringVector &lines,
						              size_t &lineNum,
						              Steiner* parent,
						              size_t expectedOffset)
{ 
  ASSERT (! subDepth);
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
		throw runtime_error ("Tree file is damaged");
	}
	
	lineNum++;

  string idS (findSplit (s));
	ASSERT (isRight (idS, ":"));
	idS. erase (idS. size () - 1);
  const Real   len          = token2real (s, "len");
  const Real   paths        = token2real (s, "paths");
  const Real   errorDensity = token2real (s, "err_density");
  const Real   relCriterion = token2real (s, rel_errorS);
  const string name         = token2string (s, "name");
  const bool indiscernible = contains (s, Leaf::non_discernible);
  IMPLY (parent, len >= 0);
  DTNode* dtNode = nullptr;
	if (isLeft (idS, "0x"))
	{
	  ASSERT (! indiscernible);
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
		leaf->discernible = ! indiscernible;
		dtNode = leaf;
	}
	ASSERT (dtNode);
	
	if (! isNan (paths))
	{
  	dtNode->paths = (size_t) DM_sp::round (paths);
  	dtNode->errorDensity = errorDensity;
  }
  if (! isNan (relCriterion))
  	dtNode->relCriterion = relCriterion;
  if (! name. empty ())
  {
    ASSERT (! dtNode->asLeaf ());
    dtNode->name = name;
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
                             const string &dissimAttrName)
{
  ASSERT (! subDepth);
	ASSERT (! dissimDs. get ());
	ASSERT (! dissimAttr);
	ASSERT (! optimizable ());

	if (dissimFName. empty ())
		throw runtime_error ("Dataset " + dissimFName + " must exist");


  // dissimAttr, dissimDs
  {
    Unverbose unv;
    dissimDs. reset (new Dataset (dissimFName));
  }
  
  {
    const Attr* attr = dissimDs->name2attr (dissimAttrName);
    if (! attr)
    	throw runtime_error ("Attribute " + dissimAttrName + " must exist in the dataset " + dissimFName);
    dissimAttr = attr->asPositiveAttr2 ();
  }
  if (! dissimAttr)
   	throw runtime_error ("Attribute " + dissimAttrName + " must have type Positive2");
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



bool DistTree::getConnected () 
{
  ASSERT (! subDepth);
	ASSERT (dissimDs. get ());
	ASSERT (dissimAttr);
	ASSERT (! optimizable ());


  // Leaf::DisjointCluster
 	for (DiGraph::Node* node : nodes)
 	{
 	  const DTNode* dtNode = static_cast <DTNode*> (node);
 	  if (Leaf* leaf = const_cast <Leaf*> (dtNode->asLeaf ()))
 	    leaf->DisjointCluster::init ();
 	}
  FFOR (size_t, row, dissimDs->objs. size ())
    if (const Leaf* leaf1 = findPtr (name2leaf, dissimDs->objs [row] -> name))
      FOR (size_t, col, row)  // dissimAttr is symmetric
        if (dissimAttr->get (row, col) < INF)
          if (const Leaf* leaf2 = findPtr (name2leaf, dissimDs->objs [col] -> name))
            const_cast <Leaf*> (leaf1) -> merge (* const_cast <Leaf*> (leaf2));

  LeafCluster leafCluster;
 	for (DiGraph::Node* node : nodes)
 	{
 	  const DTNode* dtNode = static_cast <DTNode*> (node);
 	  if (const Leaf* leaf = dtNode->asLeaf ())
      leafCluster [const_cast <Leaf*> (leaf) -> getDisjointCluster ()] << leaf;
 	}
  ASSERT (! leafCluster. empty ());
  
  if (leafCluster. size () == 1)
    return true;
    

  const DisjointCluster* cluster_main = nullptr;
  size_t size_max = 0;
  for (const auto& it : leafCluster)
    if (maximize (size_max, it. second. size ()))
      cluster_main = it. first;
  ASSERT (cluster_main);
 
  for (const auto& it : leafCluster)
    if (it. first != cluster_main)
    {
      cout << endl;
      cout << "Cluster:" << endl;
      for (const Leaf* leaf : it. second)
        cout << leaf->name << endl;
    }

  return false;
}



LeafCluster DistTree::getIndiscernibles ()
{
  ASSERT (optimizable ());

  // Leaf::DisjointCluster
 	for (DiGraph::Node* node : nodes)
 	{
 	  const DTNode* dtNode = static_cast <DTNode*> (node);
 	  if (Leaf* leaf = var_cast (dtNode->asLeaf ()))
 	    leaf->DisjointCluster::init ();
 	}

 	for (const Dissim& dissim : dissims)
 	  if (   dissim. valid ()
 	      && dissim. target <= 0
 	     )
      var_cast (dissim. leaf1) -> DisjointCluster::merge (* var_cast (dissim. leaf2));

  LeafCluster leafCluster;
 	for (DiGraph::Node* node : nodes)
 	{
 	  const DTNode* dtNode = static_cast <DTNode*> (node);
 	  if (const Leaf* leaf = dtNode->asLeaf ())
      leafCluster [var_cast (leaf) -> DisjointCluster::getDisjointCluster ()] << leaf;
 	}
  
  if (leafCluster. size () == 1)
  	throw runtime_error ("No discernible objects");
 
  for (Iter<LeafCluster> iter (leafCluster); iter. next (); )
  {
    const auto& itMap = *iter;
    ASSERT (itMap. second. size () > 0);
    if (itMap. second. size () == 1)
      iter. erase ();
  }

  return leafCluster;
}



size_t DistTree::leafCluster2discernibles (const LeafCluster &leafCluster)
{
  ASSERT (! subDepth);

  size_t n = 0;   
  for (const auto& it : leafCluster)
  {
    const VectorPtr<Leaf>& clusterNodes = it. second;
    ASSERT (! clusterNodes. empty ());
    if (clusterNodes. size () == 1)
      continue;
    const Leaf* first = clusterNodes [0];
    ASSERT (first);
    Steiner* parent = var_cast (static_cast <const DTNode*> (first->getParent ()) -> asSteiner ());
    ASSERT (parent);
    auto steiner = new Steiner (*this, parent, 0.0); 
    for (const Leaf* leaf_ : clusterNodes)
    {
      Leaf* leaf = const_cast <Leaf*> (leaf_);
      leaf->setParent (steiner);
      leaf->discernible = false;  
      leaf->len = 0;  
      steiner->pathObjNums << leaf->pathObjNums;
      n++;
    }
    ASSERT (! steiner->isTransient ());
    ASSERT (n);
    
    if (optimizable ())
    {
      VectorPtr<Leaf> leaves (clusterNodes);
      leaves. sort ();
      ASSERT (leaves. isUniq ());
      steiner->pathObjNums. sort ();
      steiner->pathObjNums. uniq ();
      for (const uint objNum : steiner->pathObjNums)
      {
        Dissim& d = dissims [objNum]; 
        if (   leaves. containsFast (d. leaf1)
            && leaves. containsFast (d. leaf2)
           )
        d. lca = steiner;  
      }
      steiner->pathObjNums. filterValue ([this, &leaves] 
                                         (const uint objNum) 
                                         { const Dissim& d = this->dissims [objNum]; 
                                           return    leaves. containsFast (d. leaf1)
                                                  && leaves. containsFast (d. leaf2);
                                         }
                                        );
    }
  }
  
  if (n)
    cleanTopology ();
    
  return n;
}



size_t DistTree::setDiscernibles_ds ()
{
  ASSERT (! subDepth);
  ASSERT (dissimDs. get ());
  ASSERT (dissimAttr);
  ASSERT (! optimizable ());


  // Leaf::DisjointCluster
 	for (DiGraph::Node* node : nodes)
 	{
 	  const DTNode* dtNode = static_cast <DTNode*> (node);
 	  if (Leaf* leaf = const_cast <Leaf*> (dtNode->asLeaf ()))
 	  {
 	    leaf->discernible = true;
 	    leaf->DisjointCluster::init ();
 	  }
 	}
  FFOR (size_t, row, dissimDs->objs. size ())
    if (const Leaf* leaf1 = findPtr (name2leaf, dissimDs->objs [row] -> name))
      FOR (size_t, col, row)  // dissimAttr is symmetric
        if (dissimAttr->get (row, col) <= 0)  // => !isNan()
          if (const Leaf* leaf2 = findPtr (name2leaf, dissimDs->objs [col] -> name))
            const_cast <Leaf*> (leaf1) -> merge (* const_cast <Leaf*> (leaf2));

  LeafCluster leafCluster;
 	for (DiGraph::Node* node : nodes)
 	{
 	  const DTNode* dtNode = static_cast <DTNode*> (node);
 	  if (const Leaf* leaf = dtNode->asLeaf ())
      leafCluster [const_cast <Leaf*> (leaf) -> getDisjointCluster ()] << leaf;
 	}
  
  if (leafCluster. size () == 1)
  	throw runtime_error ("No discernible objects");
 
  return leafCluster2discernibles (leafCluster);
}



size_t DistTree::setDiscernibles ()
{ 
  ASSERT (! subDepth);
  ASSERT (optimizable ());

 	for (DiGraph::Node* node : nodes)
 	{
 	  const DTNode* dtNode = static_cast <DTNode*> (node);
 	  if (Leaf* leaf = const_cast <Leaf*> (dtNode->asLeaf ()))
 	    leaf->discernible = true;
 	}
  const LeafCluster leafCluster (getIndiscernibles ());
  const size_t n = leafCluster2discernibles (leafCluster);

  qcPaths ();
  
  return n;
}



void DistTree::cleanTopology ()
{
 	Vector<DiGraph::Node*> nodeVec;  nodeVec. reserve (nodes. size ());

 	// delete isLeaf() && Steiner
 	insertAll (nodeVec, nodes);
 	for (DiGraph::Node* node : nodeVec)  
	  if (const Steiner* s = static_cast <DTNode*> (node) -> asSteiner ())
	    while (   s 
	           && s->isLeaf () 
	           && s->graph
	          )
	    {
   		  ASSERT (! s->inDiscernible ());
   		  ASSERT (s->childrenDiscernible ());  // Actually no children
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
 		  ASSERT (s->childrenDiscernible ());
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
  Tree::LcaBuffer buf;
  FFOR (size_t, row, dissimDs->objs. size ())
    if (const Leaf* leaf1 = findPtr (name2leaf, dissimDs->objs [row] -> name))
      FOR (size_t, col, row)  // dissimAttr is symmetric
      {
        if (const Leaf* leaf2 = findPtr (name2leaf, dissimDs->objs [col] -> name))
        {
          const Real d = dissimAttr->get (row, col);
          if (isNan (d))
            continue;
          const TreeNode* ancestor = getLca (leaf1, leaf2, buf);
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
    if (dtNode->inDiscernible ())
      { ASSERT (dtNode->len == 0); }
    else
      if (const DTNode* parent = static_cast <const DTNode*> (dtNode->getParent ()))
        dtNode->len = max (0.0, parent->subtreeLen. getMean () - dtNode->subtreeLen. getMean ());
  }

  clearSubtreeLen ();

  finishChanges ();
}



namespace
{
  struct LeafPair
  {
    // !nullptr, discernible
    array <const DTNode*, 2> nodes;
      // nodes[0] < nodes[1]
    Real dissim {NaN};
      // !isNan()

    LeafPair ()
      { nodes [0] = nullptr;
        nodes [1] = nullptr;
      }
    LeafPair (const Leaf* leaf1,
              const Leaf* leaf2,
              Real dissim_arg)
      : dissim (dissim_arg)
      { ASSERT (leaf1);
        ASSERT (leaf2);
        nodes [0] = leaf1->getDiscernible ();
        nodes [1] = leaf2->getDiscernible ();
        orderNodes ();
      }
      
    void print (ostream &os) const
      { os        << nodes [0] -> getLcaName () 
           << ' ' << nodes [1] -> getLcaName () 
           << ' ' << dissim 
           << endl; 
      }
    bool operator== (const LeafPair& other) const
      { return    nodes [0] == other. nodes [0]
               && nodes [1] == other. nodes [1];
      }
    void orderNodes ()
      { swapGreater (nodes [0], nodes [1]); }
    bool same () const
      { return nodes [0] == nodes [1]; }
    bool merge (const LeafPair &from)
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
    static bool strictlyLess (const LeafPair& n1,
                              const LeafPair& n2)
      { LESS_PART (n1, n2, nodes [0]);
        LESS_PART (n1, n2, nodes [1]);
        return false;  
      }
  };
}



void DistTree::neighborJoin ()
{
  ASSERT (isStar ());
  ASSERT (name2leaf. size () >= 2);
    
  // DTNode::len: sum of dissimilarities from other objects (dissim_sum)

  size_t n = 0;
	for (const DiGraph::Arc* arc : root->arcs [false])
	{
	  const_static_cast <DTNode*> (arc->node [false]) -> len = 0;
	  n++;
	}
    
  size_t missing = 0;
  Vector<LeafPair> leafPairs;  leafPairs. reserve (getDissimSize_max ());
  if (optimizable ())
  {
    ASSERT (dissims. size () == getDissimSize_max ());
    for (const Dissim& dissim : dissims)
    {
      if (isNan (dissim. target))
      {
      //throw runtime_error ("No distance for " + dissim. leaf1->name + " - " + dissim. leaf2->name);
        missing++;
        continue;
      }
      const LeafPair leafPair (dissim. leaf1, dissim. leaf2, /*max (0.0,*/ dissim. target/*)*/);  
      if (leafPair. same ())
        continue;
      if (leafPair. dissim == INF)
      {
      //throw runtime_error ("Infinite distance for " + dissim. leaf1->name + " - " + dissim. leaf2->name);
        missing++;
        continue;  
      }
      leafPairs << leafPair;
    }
    // non-complete dissimilarity matrix ??
  }
  else
  {
    cout << "Neighbor joining ..." << endl;  
  	ASSERT (dissimDs. get ());
  	ASSERT (dissimAttr);
    FFOR (size_t, row, dissimDs->objs. size ())
    {
      const Leaf* leaf1 = findPtr (name2leaf, dissimDs->objs [row] -> name);
      ASSERT (leaf1);
      FOR (size_t, col, row)  // dissimAttr is symmetric
      {
        const Leaf* leaf2 = findPtr (name2leaf, dissimDs->objs [col] -> name);
        ASSERT (leaf2);
        const LeafPair leafPair (leaf1, leaf2, dissimAttr->get (row, col));
        if (leafPair. same ())
          continue;
        if (leafPair. dissim < 0)
          throw runtime_error ("Negative distance for " + leaf1->name + " - " + leaf2->name);
	      if (leafPair. dissim == INF || isNan (leafPair. dissim))
	      {
	      //throw runtime_error ("Infinite distance for " + dissim. leaf1->name + " - " + dissim. leaf2->name);
	        missing++;
	        continue;  
	      }
        leafPairs << leafPair;
      }
    }
  }
  if (leafPairs. empty ())
    return;
  

  Progress prog (n - 1);
  LeafPair leafPair_best;
  for (;;)
  {
    prog ();
    
    // Remove duplicate LeafPair 
    leafPairs. sort (LeafPair::strictlyLess);
    leafPairs. filterIndex ([&leafPairs, &leafPair_best] (size_t i) 
                                 { return    leafPairs [i] == leafPair_best 
                                          || (i && leafPairs [i - 1]. merge (leafPairs [i]));
                                 }
                              );    

    if (leafPairs. size () == 1)
      break;
    ASSERT (n > 2);
    
    if (! leafPair_best. nodes [0])  // First iteration
      for (LeafPair& leafPair : leafPairs)
        for (const bool first : {false, true})
          const_cast <DTNode*> (leafPair. nodes [first]) -> len += leafPair. dissim;
          
    if (verbose (-2))  
    {
      cout << endl << "Nodes and sum:" << endl;
    	for (const DiGraph::Arc* arc : root->arcs [false])
    	{
    	  const DTNode* node = static_cast <const DTNode*> (arc->node [false]);
    	  cout << node->getLcaName () << ": " << node->len << endl;
    	}    	
    	cout << endl << "Pairs:" << endl;
      for (const LeafPair& leafPair : leafPairs)
        leafPair. print (cout);
    }
      
      
    Steiner* newNode = nullptr;
    {
      // leafPair_best
      Real dissim_min = NaN;
      {
        Real criterion = INF;
        size_t i_best = NO_INDEX;
        FFOR (size_t, i, leafPairs. size ())
          // P (criterion1 < criterion2) ??
          if (minimize (criterion, leafPairs [i]. getCriterion (n)))            
            i_best = i;
        ASSERT (i_best != NO_INDEX);
        leafPair_best = leafPairs [i_best];
        if (verbose (-1))
        {
          cout << endl << "Best: ";
          leafPair_best. print (cout);
        }
        dissim_min = leafPair_best. getParentDissim (n);
      }
    //ASSERT (dissim_min >= 0);
      ASSERT (! isNan (dissim_min));
      ASSERT (dissim_min <= leafPair_best. dissim);
      {
        DTNode* a = const_cast <DTNode*> (leafPair_best. nodes [0]);
        DTNode* b = const_cast <DTNode*> (leafPair_best. nodes [1]);
        const Real dissim_sum_a = a->len;
        const Real dissim_sum_b = b->len;
        a->len = max (0.0, dissim_min);
        b->len = max (0.0, leafPair_best. dissim - dissim_min);
        // dissim_sum
        const Real dissim_sum = (  dissim_sum_a - leafPair_best. dissim - (Real) (n - 2) * a->len
                                 + dissim_sum_b - leafPair_best. dissim - (Real) (n - 2) * b->len
                                ) / 2;
        newNode = new Steiner (*this, const_static_cast <Steiner*> (leafPair_best. nodes [0] -> getParent ()), dissim_sum);
        a->setParent (newNode);
        b->setParent (newNode);
        if (verbose (-1))
          cout << "New: " << newNode->getLcaName () << " " << a->len << " " << b->len << " " << newNode->len << endl;
      }
    }
    ASSERT (newNode);
    
    // leafPairs
    for (LeafPair& leafPair : leafPairs)
      if (! (leafPair == leafPair_best))
        for (const bool first : {false, true})
        {
          bool found = false;
          for (const bool best_first : {false, true})
            if (leafPair. nodes [first] == leafPair_best. nodes [best_first])
            {
              leafPair. nodes [first] = newNode;
              const_cast <DTNode*> (leafPair. nodes [! first]) -> len -= leafPair. dissim;
              leafPair. dissim -= leafPair_best. nodes [best_first] -> len;
              maximize (leafPair. dissim, 0.0);
              const_cast <DTNode*> (leafPair. nodes [! first]) -> len += leafPair. dissim / 2;  
                // Done twice for leafPair.nodes[!first]
              found = true;
            }
          if (found)
            leafPair. orderNodes ();
          ASSERT (! leafPair. same ());
        }

    n--;
  }
  ASSERT (leafPairs. size () == 1);
  ASSERT (n >= 2); 
  IMPLY (! missing, n == 2);
  
  
  LeafPair& leafPair = leafPairs [0];
  ASSERT (! leafPair. same ());
  for (const bool first : {false, true})  
  {
    DTNode* node = const_cast <DTNode*> (leafPair. nodes [first]);
    ASSERT (node->getParent () == root);
    node->len = leafPair. dissim / 2;  
  }
  // If dissims is not complete
	for (const DiGraph::Arc* arc : root->arcs [false])
	{
	  const DTNode* node = static_cast <const DTNode*> (arc->node [false]);
	  ASSERT (node->len >= 0);
	}
  
  if (optimizable ())
    setPaths ();

  finishChanges ();
  
  reroot (true);  // Reduces dissims.size() if dissims are sparse
}



void DistTree::dissimDs2dissims (bool sparse)
{
	ASSERT (dissimDs. get ());
	ASSERT (dissimAttr);
  ASSERT (! optimizable ());
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
  FFOR (size_t, objNum, dissimDs->objs. size ())
  {
  	const Obj* obj = dissimDs->objs [objNum];
    const Leaf* leaf = findPtr (name2leaf, obj->name);
    IMPLY (! extraObjs (), leaf);
    if (leaf)
      const_cast <Leaf*> (leaf) -> comment = obj->comment;
  }

  // Sparsing
  Vector<Pair<const Leaf*>> selectedPairs;
  if (sparse)
  {
    ASSERT (dissims. empty ());  // => dissims do not affect selectedPairs
    selectedPairs = getMissingLeafPairs_ancestors (0, false);
  }

  // dissims[], mult_sum, dissim2_sum
  loadDissimPrepare (getDissimSize_max ());
  FFOR (size_t, row, dissimDs->objs. size ())
  {
    const string name1 = dissimDs->objs [row] -> name;
    const Leaf* leaf1 = findPtr (name2leaf, name1);
    if (! leaf1)
      continue;
    FOR (size_t, col, row)  // dissimAttr is symmetric
    {
      const string name2 = dissimDs->objs [col] -> name;
      const Leaf* leaf2 = findPtr (name2leaf, name2);
      if (! leaf2)
        continue;
      const Real dissim = dissimAttr->get (row, col);
      Pair<const Leaf*> p (leaf1, leaf2);
      ASSERT (! p. same ());
      if (p. first->name > p. second->name)
        p. swap ();
      if (sparse && ! selectedPairs. containsFast (p /*getObjName (name1, name2)*/))
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


  setPaths ();
}



void DistTree::loadDissimPrepare (size_t pairs_max)
{
  ASSERT (dissims. empty ());
  
  if (verbose ())
    cout << "Leaf pairs -> data objects ..." << endl;

  dissims. reserve (pairs_max);

  const size_t reserve_size = getPathObjNums_size ();
 	for (DiGraph::Node* node : nodes)
 	{
 	  const DTNode* dtNode = static_cast <DTNode*> (node);
 	  ASSERT (dtNode->pathObjNums. empty ());
 	  const_cast <DTNode*> (dtNode) -> pathObjNums. reserve (reserve_size);  
 	}
}
  


uint DistTree::leaves2dissims (Leaf* leaf1,
                               Leaf* leaf2,
                               Real target,
                               Real mult)
{ 
  ASSERT (leaf1);
  ASSERT (leaf2);
  
  const size_t objNum_ = dissims. size ();
  dissims << Dissim (leaf1, leaf2, target, mult);
  if (objNum_ > (size_t) dissims_max)
    throw runtime_error ("leaves2dissims: Too large objNum");
  const uint objNum = (uint) objNum_;
  leaf1->pathObjNums << objNum;
  leaf2->pathObjNums << objNum;      
  
  return objNum;
}



bool DistTree::addDissim (Leaf* leaf1,
                          Leaf* leaf2,
                          Real dissim)
{
  ASSERT (leaf1);
  ASSERT (leaf2);
  ASSERT (dissim_coeff > 0);
  ASSERT (detachedLeaves. empty ());
  IMPLY (dissim == 0, leaf1->getCollapsed (leaf2));

  
  if (   isNan (dissim)  // prediction must be large ??
    //|| ! DM_sp::finite (dissim)  // For getSparseLeafPairs()
     )
  {
  	if (verbose ())
      cout << leaf1->name << " - " << leaf2->name << ": " << dissim << endl;
    return false;  
  }
  
  ASSERT (dissim >= 0);
  dissim *= dissim_coeff;  	
    
  const Real mult = dissim2mult (dissim);  
  if (mult)
  {
    mult_sum    += mult;
    dissim2_sum += mult * sqr (dissim);  // max (0.0, dissim);
  }
  
  leaves2dissims (leaf1, leaf2, dissim, mult);
  
  return true;
}



#if 0
namespace
{
	
void setPaths_ (const Vector<size_t> &subTree,
                Vector<Dissim> &dissims,
                Real &absCriterion)
// Update: dissims, absCriterion
{
	if (subTree. empty ())
		return;
  Progress prog (subTree. size (), dissim_progress); 
  if (prog. active)
  	cerr << "One thread ..." << endl;
  Tree::LcaBuffer buf;
  for (const size_t objNum : subTree) 
  {
    prog ();
    absCriterion += dissims [objNum]. setPathObjNums (objNum, buf);
  }
}
	
}
#endif



void DistTree::setPaths ()
{
//ASSERT (optimizable ());
  if (dissims. size () > dissims_max)
    throw runtime_error ("setPaths: Too large objNum");
    
  setLca ();  
  
 	for (DiGraph::Node* node : nodes)
    if (Steiner* st = const_cast <Steiner*> (static_cast <const DTNode*> (node) -> asSteiner ()))
      st->pathObjNums. clear ();  

  absCriterion = 0;
  Tree::LcaBuffer buf;
#if 0
  if (   subDepth 
  	  || threads_max == 1 
  	  || dissims. size () < smallDissims
  	 )
#endif
  {
	  Progress prog (dissims. size (), dissim_progress); 
	  FFOR (size_t, objNum, dissims. size ()) 
	  {
	  	prog ();
	    absCriterion += dissims [objNum]. setPathObjNums (objNum, buf);
	  }
 	}
#if 0
	else
  {    
	  for (Dissim& dissim : dissims)
	  {
	  	ASSERT (dissim. lca);
		  const_cast <Steiner*> (dissim. lca) -> lcaNum++;
	  }
	
	  VectorPtr<Steiner> cut_best;  cut_best. reserve (100);  // PAR
	  {
		  const Steiner* root_ = static_cast <const Steiner*> (root);
		  const_cast <Steiner*> (root_) -> setSubTreeWeight ();
		  size_t weight_min = root_->subTreeWeight;
		  ASSERT (weight_min == dissims. size ());
		  size_t upWeight = 0;
		  List<const Steiner*> cut; 
		  cut << root_;
		  List<const Steiner*>::const_iterator it_max = cut. begin ();
			while (cut. size () <= threads_max)
			{
				for (const DiGraph::Arc* arc : (*it_max)->arcs [false])
				  if (const Steiner* child = static_cast <DTNode*> (arc->node [false]) -> asSteiner ())
				  	cut << child;
				upWeight += (*it_max)->lcaNum;
				cut. erase (it_max);
			  size_t downWeight_max = 0;
			  CONST_ITER (List<const Steiner*>, it, cut)
			  	if (maximize (downWeight_max, (*it)->subTreeWeight))
			  		it_max = it;
			  if (! minimize (weight_min, downWeight_max + upWeight))
			  	break;
			  cut_best = VectorPtr<Steiner> (cut);
			}
		}
		ASSERT (cut_best. size () <= threads_max);
		
		if (cut_best. size () == 1)
		{
		  Progress prog (dissims. size (), dissim_progress);
		  FFOR (size_t, objNum, dissims. size ()) 
		  {
		    prog ();
		    absCriterion += dissims [objNum]. setPathObjNums (objNum, buf);
		  }
		}
		else
		{
			size_t threadNum = 0;
			for (const Steiner* st : cut_best)
			{
				threadNum++;
			  ASSERT (threadNum);
			  const_cast <Steiner*> (st) -> threadNum2subTree (threadNum);
			  const_cast <Steiner*> (st) -> threadNum2ancestors (0);
			}
	
		  // Due to tree disbalance, time reduction is < 10% 
		  // ??
		  cout << "root->lcaNum:" << static_cast <const Steiner*> (root) -> lcaNum << endl; 
		  cout << "root->subTreeWeight:" << static_cast <const Steiner*> (root) -> subTreeWeight << endl; 
			for (const Steiner* st : cut_best)
				cout << st->subTreeWeight << endl;
		
		  Vector<Vector<size_t>> subTrees (cut_best. size () + 1);
		  for (Vector<size_t>& subTree : subTrees)
		  	subTree. reserve (dissims. size ());  // ??
		  FFOR (size_t, objNum, dissims. size ()) 
		  {
		    Dissim& dissim = dissims [objNum];
		    const Steiner* st = dissim. lca;
		    ASSERT (st);
		  	ASSERT (st->threadNum <= cut_best. size ());
		 		subTrees [st->threadNum] << objNum;
		  }
		  
		  // ??
		  cout << endl;
		  cout << "dissims:" << dissims. size () << endl;
		  for (const auto& subTree: subTrees)
		    cout << subTree. size () << endl;  
		  
		  vector<thread> threads;  threads. reserve (cut_best. size () - 1);
		  MVector absCriteria (cut_best. size () - 1, 0);
		  size_t processed = 0;
		  FFOR_START (size_t, i, 1, subTrees. size () - 1)
		    if (! subTrees [i]. empty ())
			  {
				  threads. push_back (thread (setPaths_, cref (subTrees [i]), ref (dissims), ref (absCriteria [i - 1])));	
				  processed += subTrees [i]. size ();
				}
			ASSERT (processed + subTrees [0]. size () + subTrees [cut_best. size ()]. size () == dissims. size ());
		  setPaths_ (subTrees [cut_best. size ()], dissims, absCriterion);
			for (auto& t : threads)  
				t. join ();
		  absCriterion += absCriteria. sum ();
		  	  
		  setPaths_ (subTrees [0], dissims, absCriterion);
    }
	}
#endif
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

 	for (const DiGraph::Node* node : nodes)
 	{
 	  const DTNode* dtNode = static_cast <const DTNode*> (node);
    ASSERT (! dtNode->isTransient ());
  }

  const size_t leavesSize = root->getLeavesSize ();
  ASSERT (leavesSize);
  ASSERT (name2leaf. size () == leavesSize);
  Vector<size_t> leafIndices;  leafIndices. reserve (leavesSize);
  for (const auto& it : name2leaf)
  {
  	const Leaf* leaf = it. second;
    ASSERT (leaf);
    ASSERT (leaf->graph == this);
    leafIndices << leaf->index;
  }
  leafIndices. sort ();
  leafIndices. uniq ();
  ASSERT (leafIndices. size () == leavesSize);
  
  ASSERT ((bool) dissimDs. get () == (bool) dissimAttr);
  if (dissimAttr)
  {
    dissimDs->qc ();
    ASSERT (& dissimAttr->ds == dissimDs. get ());
  }


  size_t leafDissims = 0;
	if (optimizable ())
	{
    Set<pair<const Leaf*,const Leaf*>> leafSet;
    for (const Dissim& dissim : dissims)
    {
      dissim. qc ();
      leafSet. addUnique (pair<const Leaf*,const Leaf*> (dissim. leaf1, dissim. leaf2));
    }
    
   	if (! isNan (absCriterion))
   	{
   	  ASSERT (absCriterion >= 0);
    //ASSERT (absCriterion_delta >= 0);
    }
   	
   	   	
   	for (const DiGraph::Node* node : nodes)
   	{
   	  const DTNode* dtNode = static_cast <const DTNode*> (node);
      ASSERT ((dtNode == root) == dtNode->pathObjNums. empty ());  
   		if (const Leaf* leaf = dtNode->asLeaf ())
   		  leafDissims += leaf->pathObjNums. size ();
  	}  
  }
  

  const size_t discernibles = getDiscernibles (). size ();
  ASSERT (discernibles <= leavesSize);

  for (const Leaf* leaf : detachedLeaves)
  {
    ASSERT (leaf);
    ASSERT (! leaf->graph);
    leaf->qc ();
	  leafDissims += leaf->pathObjNums. size ();
  }
  const size_t allLeaves = leavesSize + detachedLeaves. size ();
  ASSERT (dissims. size () <= (allLeaves * (allLeaves - 1)) / 2);

 	ASSERT (leafDissims == 2 * dissims. size ());

/*
  // May happen due to removeLeaf()
  for (const DTNode* dtNode : toDelete)
  	ASSERT (! dtNode->asLeaf ());
*/

  if (! subDepth && optimizable ())
  {
    const LeafCluster leafCluster (var_cast (this) -> getIndiscernibles ());
    for (const auto& it : leafCluster)
    {
      const VectorPtr<Leaf>& leaves = it. second;
      ASSERT (leaves. size () >= 2);
      for (const Leaf* leaf : leaves)
      {
        ASSERT (leaf->graph == this);
        ASSERT (leaf->getParent () == leaves [0] -> getParent ());
      }
    }
   	for (const DiGraph::Node* node : nodes)
   	{
   	  const DTNode* dtNode = static_cast <const DTNode*> (node);
   	  if (dtNode->childrenDiscernible ())
   	    continue;
   	  const DisjointCluster* dc = nullptr;
    	for (const DiGraph::Arc* arc : dtNode->arcs [false])
    	{
    	  Leaf* child = var_cast (static_cast <DTNode*> (arc->node [false]) -> asLeaf ());
    	  ASSERT (child);
    	  if (dc)
    	  { 
    	    if (child->DisjointCluster::getDisjointCluster () != dc)
    	    {
    	      cout << "Bad indiscernible " << child->name << endl;
    	      ERROR;
    	    }
    	  }
    	  else
    	    dc = child->DisjointCluster::getDisjointCluster ();
   	  }
   	}
  }
}



void DistTree::deleteLeaf (TreeNode* node,
                           bool deleteTransientAncestor)
{
  ASSERT (! subDepth);
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



const DTNode* DistTree::lcaName2node (const string &lcaName,   
	                                    Tree::LcaBuffer &buf) const
{
  ASSERT (! lcaName. empty ());

  string s (lcaName);
  const string name1 = findSplit (s, objNameSeparator);
  const Leaf* leaf1 = findPtr (name2leaf, name1);
  ASSERT (leaf1);
  const Leaf* leaf2 = s. empty () ? leaf1 : findPtr (name2leaf, s);
  ASSERT (leaf2);
  const DTNode* node = static_cast<const DTNode*> (getLca (leaf1, leaf2, buf));
  ASSERT (node);
  
  return node;
}



Set<const DTNode*> DistTree::getDiscernibles () const
{
  Set<const DTNode*> s;
  for (const DiGraph::Node* node : nodes)
    if (const Leaf* leaf = static_cast <const DTNode*> (node) -> asLeaf ())
      s << leaf->getDiscernible ();
  return s;
}

  

void DistTree::printInput (ostream &os) const
{
  os << "INPUT:" << endl;
  os << "# Leaves: " << name2leaf. size () << endl;
  const size_t discernibles = getDiscernibles (). size ();
  os << "# Discernible leaves: " << discernibles << endl;
  os << "# Nodes: " << nodes. size () << endl;
  if (! optimizable ())
    return;
  {
	  const ONumber on (os, 2, false);  // PAR
	  os << "# Dissimilarities: " << dissims. size () << " (" << (Real) dissims. size () / (Real) getDissimSize_max () * 100 << " %)" << endl; 
	  os << "Dissimilarities factor: " << (Real) dissims. size () / ((Real) discernibles * log ((Real) discernibles)) << endl;
	}
	{
	  const ONumber on (os, dissimDecimals, false);
	  os << "Ave. dissimilarity = " << getDissim_ave () << endl;
	}
  reportErrors (os);
}



void DistTree::saveFeatureTree (const string &fName) const
{
  if (fName. empty ())
    return;
  
  OFStream f (fName);
  static_cast <const DTNode*> (root) -> saveFeatureTree (f, 0);
}



void DistTree::qcPaths () 
{
  if (! qc_on)
    return;

  const size_t displayNodes = 1000;  // PAR
  if (nodes. size () >= displayNodes)
    cerr << "QC paths ..." << endl;

  size_t pathObjNums_checked = 0;
  size_t lcaObjNums_checked = 0;
  {
    Progress prog (nodes. size (), displayNodes);  
    for (DiGraph::Node* node : nodes)
    {
      prog ();
      DTNode* dtNode = const_static_cast <DTNode*> (node);
      Vector<uint>& pathObjNums = dtNode->pathObjNums;
      pathObjNums. sort ();
      ASSERT (pathObjNums. isUniq ());
      for (const uint objNum : pathObjNums)
        if (dissims [objNum]. valid ())
          pathObjNums_checked++;
      const Vector<uint> lcaObjNums (dtNode->getLcaObjNums ());
      for (const uint objNum : lcaObjNums)
        if (dissims [objNum]. valid ())
        {
          ASSERT (dissims [objNum]. lca == dtNode->asSteiner ());
          lcaObjNums_checked++;
        }
    }
  }
  
  size_t pathObjNums_all = 0;
  size_t lcaObjNums_all = 0;
  {
    Progress prog (dissims. size (), dissim_progress);
    LcaBuffer buf;
    FFOR (size_t, objNum, dissims. size ())
    {
      prog ();
      if (dissims [objNum]. valid ())
      {
        const VectorPtr<TreeNode>& path = dissims [objNum]. getPath (buf);
        for (const TreeNode* node : path)
        {
          ASSERT (static_cast <const DTNode*> (node) -> pathObjNums. containsFast ((uint) objNum));
          pathObjNums_all++;
        }
        lcaObjNums_all++;
      }
    }
  }

  ASSERT (pathObjNums_all == pathObjNums_checked);
  ASSERT (lcaObjNums_all == lcaObjNums_checked);
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



void DistTree::clearSubtreeLen ()
{
  for (DiGraph::Node* node : nodes)
    static_cast <DTNode*> (node) -> subtreeLen. clear ();
}



namespace
{

void setPredictionAbsCriterion_thread (size_t from,
								                       size_t to,
								                       Real &absCriterion,
								                       Vector<Dissim> &dissims)
{
  absCriterion = 0;
  Tree::LcaBuffer buf;
  Progress prog (dissims. size (), dissim_progress); 
	FOR_START (size_t, i, from, to)
	{
		Dissim& dissim = dissims [i];
    if (! dissim. valid ())
    	continue;
    prog ();
	  const VectorPtr<Tree::TreeNode>& path = dissim. getPath (buf);
	  dissim. prediction = DistTree::path2prediction (path);  
    absCriterion += dissim. getAbsCriterion ();
	}
}

}



void DistTree::setPredictionAbsCriterion ()
{
  ASSERT (! subDepth);

  absCriterion = 0;
#if 1
  vector<Real> absCriteria;
  arrayThreads (setPredictionAbsCriterion_thread, dissims. size (), absCriteria, ref (dissims)); 
  for (const Real x : absCriteria)
  	absCriterion += x;
#else
  LcaBuffer buf;
  Progress prog (dissims. size (), dissim_progress); 
  for (Dissim& dissim : dissims)
    if (dissim. valid ())
    {
    	prog ();
		  const VectorPtr<TreeNode>& path = dissim. getPath (buf);
		  dissim. prediction = path2prediction (path);  
      absCriterion += dissim. getAbsCriterion ();
    }
#endif
}



void DistTree::qcPredictionAbsCriterion () const
{
  if (! qc_on)
    return;

  Real absCriterion_ = 0;
  LcaBuffer buf;
  for (const Dissim& dissim : dissims)
    if (dissim. valid ())
    {
		  const VectorPtr<TreeNode>& path = dissim. getPath (buf);
		  const Real prediction_ = path2prediction (path);
      ASSERT_EQ (prediction_, dissim. prediction, 1e-3);  // PAR
      absCriterion_ += dissim. getAbsCriterion ();
    }
  if (fabs (absCriterion - absCriterion_) > 1e-3)  // PAR
  {
  	cout << absCriterion << " " << absCriterion_ << endl;
  	ERROR;
  }
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



void DistTree::printAbsCriterion_halves () const
{
  ASSERT (optimizable ());
  
  const Real dissim_ave = getDissim_ave ();
  size_t pairs = 0;
  size_t size_half       [2/*bool*/] = {0, 0};
  Real absCriterion_half [2/*bool*/] = {0, 0};
  Real dissim2_half      [2/*bool*/] = {0, 0};
  for (const Dissim& dissim : dissims)
    if (dissim. mult)
    {
      const Real d = dissim. target;
      const Real dHat = dissim. prediction;
      ASSERT (! isNan (dHat));
      const bool half2 = d > dissim_ave;
      pairs++;
      size_half [half2] ++;
      absCriterion_half [half2] += dissim. mult * sqr (dHat - d);
      dissim2_half      [half2] += dissim. mult * sqr (d);
    }
    
  for (const bool half2 : {false, true})
  {
    const ONumber on (cout, criterionDecimals, false);  // PAR
    cout << "absCriterion [" << half2 << "] = " << absCriterion_half [half2];
    const Prob unexplainedFrac = absCriterion_half [half2] / dissim2_half [half2];
    const Prob r = 1 - unexplainedFrac;
    const Real arcLenError = sqrt ((1 - r) / r); 
    cout << "  Error density [" << half2 << "] = " << arcLenError * 100 << " %";
    cout << "  Fraction = " << (Real) size_half [half2] / (Real) pairs * 100 << " %" << endl;
  }
}



#if 0
void DistTree::setLeafAbsCriterion ()  
// Output: Leaf::{absCriterion,absCriterion_ave}
// Time: O(p)
{
  ASSERT (optimizable ());
  
  for (DiGraph::Node* node : nodes)
    if (Leaf* leaf = const_cast <Leaf*> (static_cast <DTNode*> (node) -> asLeaf ()))
    {
      leaf->absCriterion = 0;
      leaf->absCriterion_ave = 0;  // tepmorary: = mult
    }

  for (const Dissim& dissim : dissims)
    if (dissim. mult)
    {
      const Real residual2 = sqr (dissim. getResidual ());      
      ASSERT (DM_sp::finite (residual2));

      Leaf* leaf1 = const_cast <Leaf*> (dissim. leaf1);
      Leaf* leaf2 = const_cast <Leaf*> (dissim. leaf2);

      leaf1->absCriterion += residual2 * dissim. mult;
      leaf2->absCriterion += residual2 * dissim. mult;

      leaf1->absCriterion_ave += dissim. mult;
      leaf2->absCriterion_ave += dissim. mult;
    }

//OFStream f ("leaf_criterion");  
  for (DiGraph::Node* node : nodes)
    if (const Leaf* leaf = static_cast <DTNode*> (node) -> asLeaf ())
    {
      const_cast <Leaf*> (leaf) -> absCriterion_ave = leaf->absCriterion == 0 ? 0 : leaf->absCriterion / leaf->absCriterion_ave;
      ASSERT (leaf->absCriterion_ave >= 0);
    }
}
#endif



void DistTree::setNodeAbsCriterion () 
{
  ASSERT (optimizable ());
  
  for (DiGraph::Node* node : nodes)
  {
  	DTNode* dtNode = static_cast <DTNode*> (node);
    dtNode->absCriterion = 0;
    dtNode->absCriterion_ave = 0;  // tepmorary: = count
  }

  Tree::LcaBuffer buf;
  for (const Dissim& dissim : dissims)
  	if (dissim. mult)
	  {
	  	const Real criterion = dissim. getAbsCriterion ();
	    ASSERT (DM_sp::finite (criterion));
	    ASSERT (criterion >= 0);
	    
	    // For interior nodes ??
	    ASSERT (dissim. lca);	
	  	VectorPtr<TreeNode>& path = dissim. getPath (buf);
	  	path << dissim. lca;
	  	
	  	for (const TreeNode* node : path)
	  	{
		  	DTNode* dtNode = const_static_cast <DTNode*> (node);
	      dtNode->absCriterion     += criterion;  
	      dtNode->absCriterion_ave += 1; 
	  	}
	  }

  for (DiGraph::Node* node : nodes)
  {
  	DTNode* dtNode = static_cast <DTNode*> (node);
    dtNode->absCriterion_ave = dtNode->absCriterion == 0 ? 0 : dtNode->absCriterion / dtNode->absCriterion_ave;
    ASSERT (dtNode->absCriterion_ave >= 0);
  }
}



#if 0
void DistTree::getSkipRetain (DTNode* &toSkip,
                              DTNode* &toRetain)
// Output: toSkip, toRetain; may be nullptr
// Update: {toSkip,toRetain}->len
// toSkip->attr is redundant
// toSkip->getparent() = toRetain->getParent() = root: the only children
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
  VectorPtr<DTNode> dtNodes;  dtNodes. reserve (2 * name2leaf. size ());
  size_t arcs = 0;
  for (DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <DTNode*> (node);
    if (   dtNode != root
        && dtNode != toSkip
        && ! dtNode->inDiscernible ()
       )
    {
      dtNodes << dtNode;
      const_cast <DTNode*> (dtNode) -> index = arcs;
      arcs++;
    }
    IMPLY (dtNode->inDiscernible (), dtNode->len == 0);
  }
  ASSERT (arcs == dtNodes. size ());

  // DTNode::dissimSum, matr
  Matrix matr (arcs, 0);
//matr. putAll (0);
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
        
  FFOR (size_t, i, dtNodes. size ())
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
    FFOR (size_t, i, path. size ())
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
    if (a->inDiscernible ())
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
  ASSERT (a->len >= 0);
  ASSERT (b->len >= 0);
  return a->len > b->len;  
}

}



size_t DistTree::optimizeLenArc ()
{
  if (verbose (1))
    cout << "Optimizing arc lengths at each arc ..." << endl;
    
  VectorPtr<DTNode> dtNodes;  dtNodes. reserve (2 * name2leaf. size ());
  for (DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <DTNode*> (node);
    if (   dtNode != root
        && ! dtNode->inDiscernible ()
       )
      dtNodes << dtNode;
    IMPLY (dtNode->inDiscernible (), dtNode->len == 0);
  }
  dtNodes. sort (DTNode_len_strictlyLess);
  

  const uint iters = 10;  // PAR
  FOR (uint, iter, iters)
  {
    Real absCriterion_prev = absCriterion;
    {
      Progress prog (dtNodes. size (), (dtNodes. size () >= 10000) * 100);  // PAR
      for (const DTNode* node : dtNodes)
      {
        prog (absCriterion2str ());
      #ifndef NDEBUG    
        const Real absCriterion_old = absCriterion;  
      #endif
      
        Real arcAbsCriterion_old = 0;
        WeightedMeanVar mv;
        for (const uint objNum : node->pathObjNums)
        {
          const Dissim& dissim = dissims [objNum];
          if (dissim. valid ())
          {
            const Real dist_hat_tails = max (0.0, dissim. prediction - node->len);
            const Real arcTarget = dissim. target - dist_hat_tails;
            mv. add (arcTarget, dissim. mult);
            arcAbsCriterion_old += dissim. getAbsCriterion ();
          }
        }
        const Real len_new = max (0.0, mv. getMean ());
        
        Real arcAbsCriterion_new = 0;
        for (const uint objNum : node->pathObjNums)
        {
          Dissim& dissim = dissims [objNum];
          if (dissim. valid ())
          {
            dissim. prediction = max (0.0, dissim. prediction - node->len + len_new);
            arcAbsCriterion_new += dissim. getAbsCriterion ();
          }
        }
        ASSERT (leReal (arcAbsCriterion_new, arcAbsCriterion_old));
        
        const_cast <DTNode*> (node) -> len = len_new;
        absCriterion -= arcAbsCriterion_old;
        absCriterion += arcAbsCriterion_new;
        ASSERT (leReal (absCriterion, absCriterion_old));
        maximize (absCriterion, 0.0);
      }
    }
  
    if (Progress::enabled ())
      cerr << '\r' << iter + 1 << " / " << (uint) iters << " " << absCriterion2str ();
    if (absCriterion_prev / absCriterion - 1 <= 0.01)   // PAR
      break;
    absCriterion_prev = absCriterion;  
  }
  if (Progress::enabled ())
    cerr << endl;


  return finishChanges ();
}



namespace
{
  
struct Star
{
  // !nullptr
  const Steiner* center;
  VectorPtr<DiGraph::Node> arcNodes;
    // Arc's make up a star
  Real lenSum;
  
  explicit Star (const Steiner* center_arg)
    : center   (center_arg)
    , arcNodes (center->getChildren ())
    , lenSum (0)
    {
      ASSERT (center);
      ASSERT (! center->isTransient ());
      arcNodes << center;
      ASSERT (arcNodes. size () >= 3);
      for (const DiGraph::Node* node : arcNodes)
        lenSum += static_cast <const DTNode*> (node) -> len;
      ASSERT (lenSum >= 0);
    }
   
  static bool strictlyLess (const Star &a,
                            const Star &b)
    {
      ASSERT (! isNan (a. lenSum));
      ASSERT (! isNan (b. lenSum));
      return a. lenSum > b. lenSum;  
    }
};
  
}



size_t DistTree::optimizeLenNode ()  
{
  if (verbose (1))
    cout << "Optimizing arc lengths at each node ..." << endl;

  Vector<Star> stars;  stars. reserve (2 * name2leaf. size ());
  for (const DiGraph::Node* node : nodes)
  {
    const DTNode* center = static_cast <const DTNode*> (node);
    if (   center == root
        || center->asLeaf () 
        || ! center->childrenDiscernible ()
       )
      continue;
    stars << Star (center->asSteiner ());
  }
  stars. sort (Star::strictlyLess); 


  Progress prog (stars. size ());  
#ifndef NDEBUG
  const Real absCriterion_old1 = absCriterion;
#endif
  Tree::LcaBuffer buf;
  for (const Star& star : stars)
  {
    const VectorPtr<DiGraph::Node>& arcNodes = star. arcNodes;

    Subgraph subgraph (*this);
    subgraph. reserve (3);  // PAR
    for (const DiGraph::Node* node : star. arcNodes)
      subgraph. area << static_cast <const TreeNode*> (node);
    subgraph. area << star. center->getParent ();
    subgraph. boundary = subgraph. area;
    ASSERT (subgraph. boundary. size () >= 4);
    const size_t centerPos = subgraph. boundary. size () - 2;
    ASSERT (subgraph. boundary [centerPos] == star. center);
    subgraph. boundary. eraseAt (centerPos);
    subgraph. finish ();
    subgraph. area2subPaths ();
    subgraph. finishSubPaths ();
    subgraph. qc ();    
    // subpaths --> arc pairs ??
  
    Dataset starDs;
    starDs. objs. reserve (subgraph. subPaths. size ());
    auto targetAttr = new RealAttr1 ("target", starDs);
    FFOR (size_t, objNum, subgraph. subPaths. size ())
      starDs. appendObj ();
    Space1<NumAttr1> sp (starDs, false);  sp. reserve (star. arcNodes. size ());
    FFOR (size_t, i, star. arcNodes. size ())
    {
      auto attr = new ExtBoolAttr1 ("X" + toString (i + 1), starDs);;
      sp << attr;
      attr->setAll (EFALSE);
    }
    FFOR (size_t, objNum, subgraph. subPaths. size ())
    {
      const SubPath& subPath = subgraph. subPaths [objNum];
      const VectorPtr<TreeNode>& path = subgraph. getPath (subPath, buf);
      FFOR (size_t, i, star. arcNodes. size ())
        if (path. contains (static_cast <const TreeNode*> (star. arcNodes [i])))
          (* const_static_cast <ExtBoolAttr1*> (sp [i])) [objNum] = ETRUE;        
      const size_t wholeObjNum = subPath. objNum;
      const_cast <Obj*> (starDs. objs [objNum]) -> mult = dissims [wholeObjNum]. mult; 
      (*targetAttr) [objNum] = dissims [wholeObjNum]. target - subPath. dist_hat_tails;
    }
    starDs. qc ();
    sp. qc ();        
    const Sample sample (starDs);
    
    L2LinearNumPrediction lr (sample, sp, *targetAttr);
    ASSERT (lr. beta. size () == arcNodes. size ());
    FFOR (size_t, attrNum, lr. beta. size ())
      lr. beta [attrNum] = static_cast <const DTNode*> (arcNodes [attrNum]) -> len;
    const bool solved = lr. solveUnconstrainedFast (nullptr, true, 10, 0.01);  // PAR
    lr. qc ();
  
    // DTNode::len
    if (solved)
    {
      FFOR (size_t, attrNum, arcNodes. size ())
        const_static_cast <DTNode*> (arcNodes [attrNum]) -> len = lr. beta [attrNum];
      const Real absCriterion_old = absCriterion;  
      subgraph. subPaths2tree ();
      ASSERT (leReal (absCriterion, absCriterion_old, 1e-4));  // PAR
      prog (absCriterion2str ()); 
    }
  }

  
#ifndef NDEBUG
  if (verbose ())
  {
    const ONumber on (cout, criterionDecimals, false);
    cout << absCriterion_old1 << " -> " << absCriterion << endl;
  }
  if (! leReal (absCriterion, absCriterion_old1))
  {
    const ONumber on (cout, 6, true);
  	cout << absCriterion_old1 << ' ' << absCriterion << endl;
  	ERROR;
  }
#endif

  return finishChanges ();
}



void DistTree::optimize2 () 
{
  ASSERT (dissims. size () == 1);
  Dissim& dissim = dissims [0];
  ASSERT (! isNan (dissim. target));  // otherwise tree is disconnected

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

  const Real t =  max (0.0, dissim. target);
  for (const Leaf* leaf : leaves)
    const_cast <Leaf*> (leaf) -> len = t / 2;
  
  dissim. prediction = t;
  absCriterion = dissim. getAbsCriterion ();

//setPaths ();
}



void DistTree::optimize3 () 
{
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
    Dissim& dissim = dissims [i];
    dissim. prediction = 0;
    const Real t = isNan (dissim. target) ? 0 : dissim. target;
    for (const Leaf* leaf_ : leaves)
    {
      Leaf* leaf = const_cast <Leaf*> (leaf_);
      if (dissim. hasLeaf (leaf))
        leaf->len += t;
      else
        leaf->len -= t;
    }
  }

  for (const Leaf* leaf_ : leaves)
  {
    Leaf* leaf = const_cast <Leaf*> (leaf_);
    leaf->len /= 2;  // >= 0 <= Dissim::target is a distance and triangle inequality 
    maximize (leaf->len, 0.0);
  }

  absCriterion = 0;
  FOR (size_t, i, 3)
  {
    Dissim& dissim = dissims [i];
    for (const Leaf* leaf : leaves)
      if (dissim. hasLeaf (leaf))
        dissim. prediction += leaf->len;
    absCriterion += dissim. getAbsCriterion ();
  }
  
  cleanTopology ();

  setPaths ();
}



namespace
{
  
void reinsert_thread (size_t from, 
                      size_t to, 
                      VectorPtr<Change> &changes, 
                      const DistTree &tree,
                      const VectorPtr<DTNode> &nodeVec
                     )
{ 
  
  ASSERT (from <= to);
  ASSERT (changes. empty ());
  ASSERT (to <= nodeVec. size ());
  
  const size_t q_max = 10 * tree. getSparseDissims_size ();  // PAR
  Tree::LcaBuffer buf;
  Progress prog (to - from);
  FOR_START (size_t, i, from, to)
  {
    prog ();
    const DTNode* fromNode = nodeVec [i];
    Real nodeAbsCriterion_old = NaN;
    const NewLeaf nl (fromNode, q_max, nodeAbsCriterion_old);
    ASSERT (nodeAbsCriterion_old >= 0);
    nl. qc ();
    const DTNode* toNode = nl. location. anchor;
    ASSERT (toNode);
    const Real improvement = nodeAbsCriterion_old - nl. location. absCriterion_leaf;
    if (   fromNode->getParent () == toNode
        || fromNode->getParent () == toNode->getParent ()
        || ! positive (improvement)
       )
      continue;
 	  if (! Change::valid (fromNode, toNode))
 	    continue;
 	  const DistTree::TreeNode* lca = nullptr;
 	  const size_t arcDist = Tree::getPath (fromNode, toNode, nullptr, lca, buf). size ();
 	  ASSERT (arcDist > 1);
 	  if (verbose ())
      cout << fromNode->getLcaName () << " -> " << toNode->getLcaName () << ' ' << improvement << ' ' << arcDist << endl;
    if (arcDist < areaDiameter_std)
      continue;
 	  auto change = new Change (fromNode, toNode);
 	  change->improvement = improvement;
 	  changes << change;
  }
}
            
}



bool DistTree::optimizeReinsert ()
{
  ASSERT (! subDepth);

  VectorPtr<DTNode> nodeVec;  nodeVec. reserve (nodes. size ());
  for (const DiGraph::Node* node_ : nodes)
  {
    const DTNode* node = static_cast <const DTNode*> (node_);
    if (   ! node->inDiscernible ()
        && node != root
       )
      nodeVec << node;
  }
  nodeVec. randomOrder ();

  vector<VectorPtr<Change>> results;
  arrayThreads (reinsert_thread, nodeVec. size (), results, cref (*this), cref (nodeVec));

	VectorOwn<Change> changes;  changes. reserve (256);  // PAR	
  for (const VectorPtr<Change>& threadChanges : results)
    changes << threadChanges;

  return applyChanges (changes, true);
}



void DistTree::optimizeWholeIter (uint iter_max,
                                  const string &output_tree)
{
  for (DiGraph::Node* node : nodes)
    static_cast <DTNode*> (node) -> stable = false;

  uint iter = 0;
  while (! iter_max || iter < iter_max)
  {
    iter++;
    if (verbose (1))
      cout << endl << "Topology optimization iter = " << iter << endl;
    if (! optimizeWhole ())
      break;	    
    saveFile (output_tree);
    if (verbose (1))
      cout. flush ();
  }
}



bool DistTree::optimizeWhole () 
{ 
  ASSERT (dissims. size () >= 2 * name2leaf. size () - 2);

  
	VectorOwn<Change> changes;  changes. reserve (256);  // PAR
	{
	 	Vector<DTNode*> nodeVec;  nodeVec. reserve (2 * name2leaf. size ());
    for (DiGraph::Node* node : nodes)
    {
      DTNode* dtNode = static_cast <DTNode*> (node);
      if (! dtNode->stable)
        nodeVec << dtNode;
    }
		Progress prog (nodeVec. size ());
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
    
  return applyChanges (changes, false);
}



const Change* DistTree::getBestChange (const DTNode* from) 
{
	ASSERT (from);
	
 	const Change* bestChange = nullptr;
 	
 	const size_t reserve_size = (size_t) min (pow (2, areaRadius_std + 1), 2 * (Real) name2leaf. size ());
  VectorPtr<TreeNode> area;      area.     reserve (reserve_size); 
  VectorPtr<TreeNode> boundary;  boundary. reserve (reserve_size); 
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



bool DistTree::applyChanges (VectorOwn<Change> &changes,
                             bool byNewLeaf)
{ 
	ASSERT (toDelete. empty ());
	ASSERT (absCriterion >= 0);
	
	
  const Real absCriterion_init = absCriterion;


  // DTNode::stable: init
  for (DiGraph::Node* node : nodes)
    static_cast <DTNode*> (node) -> stable = true;


  if (verbose (1))
    cout << "# Changes: " << changes. size () << endl;
  size_t commits = 0;
  {
  //Unverbose un;
    changes. sort (Change::strictlyLess);	
    size_t nChange = 0;
    Progress prog (changes. size ());
  	for (const Change* ch_ : changes)
  	{
  	  prog (absCriterion2str ());
  	  nChange++;
  	  
  		Change* ch = const_cast <Change*> (ch_);
  		ASSERT (ch);
  		ASSERT (ch->improvement > 0);

  	  if (! ch->valid ())
    		continue;
    		
      Unverbose un;  
      if (verbose ())
  		  ch->qc ();  
  
  	#ifndef NDEBUG
  	  const bool first = ch_ == changes. front ();  
  	#endif
  
      qcPaths (); 
  	    
   	  const Real absCriterion_old = absCriterion;
  
      if (verbose (1))
      {
  	    cout << "Apply " << nChange << "/" << changes. size () << ": ";
     	  ch->print (cout);  
     	}
  	  const bool success = ch->apply ();
      if (verbose (1))
  	    cout << "Success: " << success << "  improvement = " << ch->improvement << endl;
  	  IMPLY (first && ! byNewLeaf, success && ch->improvement > 0);
  	  if (! success || ch->improvement <= 0) 
  	  {
  	  //ASSERT (! first);
        if (verbose (1))
  	      cout << "Restore" << endl;
  	  	ch->restore ();
  	  }
  	  else
  	  {
       	ch->commit ();
       	commits++;
    	  if (verbose ())
    	  {
     	    cout << "absCriterion = " << absCriterion << endl;
    	    Unverbose unv;
    	    if (verbose ())
    	      qc ();
    	  }
        ASSERT (leReal (absCriterion, absCriterion_old));
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
    ONumber on (cout, criterionDecimals, false);
    cout << "# Commits = " << commits << endl;
    cout << "Improvement = " << improvement /*<< "  from: " << absCriterion_init << " to: " << absCriterion*/ << endl;
  }
  ASSERT ((! commits) == (improvement == 0));

  if (commits)
  {
    if (byNewLeaf)
      finishChanges ();
    else
    {
      optimizeLenArc (); 
      optimizeLenNode ();  
    }
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

  ch->apply ();
  ASSERT (! isNan (ch->improvement));
  ch->restore ();

  Unverbose unv;
  if (verbose ())
    ch->print (cout); 
  
	if (Change::strictlyLess (ch, bestChange))
	{
		delete bestChange;
  	bestChange = ch;
  }
  else
  	delete ch;
}



namespace
{

void processLargeImage_thread (Image &image,
				                       const Steiner* subTreeRoot,
					                     const VectorPtr<Tree::TreeNode> possibleBoundary)  // needs a copy
{ 
	image. processLarge (subTreeRoot, possibleBoundary); 
}

}
	                      


void DistTree::optimizeLargeSubgraphs ()
{
  ASSERT (! subDepth);
  ASSERT (threads_max);
  
  constexpr size_t large_min = 100;  // PAR
  ASSERT (large_min > 1);


	if (threads_max == 1)
	{
		optimizeSmallSubgraphs (areaRadius_std, false);
    return;
	}
	if (name2leaf. size () <= 2 * large_min)  // PAR
	{
	  cout << "Too small tree" << endl;
		optimizeSmallSubgraphs (areaRadius_std, false);
    return;
	}
	

  const Keep<size_t> threads_max_kp (threads_max);
	for (;;)
	{
		VectorPtr<Steiner> boundary;  boundary. reserve (threads_max);   // isTransient()
		{
			// Improves speed without thread's ??!
			VectorPtr<Steiner> cuts;  cuts. reserve (threads_max);
			{
				const_static_cast <Steiner*> (root) -> subtreeSize2leaves ();
				const size_t rootLeaves = static_cast <const Steiner*> (root) -> leaves;
				ASSERT (rootLeaves + 1 == nodes. size ());
				Normal normal;
				normal. setSeed ((ulong) clock ());
				constexpr Prob topSubgraphFrac = 0.2;  // PAR  // was: 0.5
				const Real goalSize_init = (Real) rootLeaves / ((Real) threads_max - 1 + topSubgraphFrac);  
				normal. setParam (0, max<Real> (1, goalSize_init * 0.02));  // PAR
				FFOR (size_t, threadNum, threads_max - 1)
				{
				  // cuts nodes should vary after tree addition
					const size_t goalSize = max (large_min, (size_t) max<Real> (0, goalSize_init + normal. rand ()));
					const Steiner* st_best = nullptr;
					size_t diff = numeric_limits<size_t>::max ();
			    for (const DiGraph::Node* node : nodes)
			    	if (const Steiner* st = static_cast <const DTNode*> (node) -> asSteiner ())
			  	    if (   st->leaves > max<size_t> (large_min, goalSize / 2)  // PAR
			  	    	  && st->childrenDiscernible ()
			  	    	  && minimize (diff, difference (st->leaves, goalSize))
			  	    	 )
			  	    	st_best = st;
			  	if (! st_best)
			  		break;
			  	if (st_best == root)
			  		break;
			    cuts << st_best;
			    const size_t leaves_best = st_best->leaves;
			    ASSERT (leaves_best);
			    const Steiner* st = st_best;
			    while (st)
			    {
			    	if (st->leaves)
			    	{
				    	ASSERT (st->leaves >= leaves_best);
				    	const_cast <Steiner*> (st) -> leaves -= leaves_best;
				    }
			    	st = static_cast <const Steiner*> (st->getParent ());
			    }
			    const_cast <Steiner*> (st_best) -> setLeaves (0);
			    if (TreeNode* parent = const_cast <TreeNode*> (st_best->getParent ()))
			    	parent->leaves = 0;  // To prevent neighboring boundary nodes
				}
			}
			ASSERT (boundary. empty ());
			for (const Steiner* cut : cuts)
			{
				ASSERT (cut);
		    auto inter = new Steiner (*this, const_static_cast <Steiner*> (cut->getParent ()), 0);
		    inter->pathObjNums = cut->pathObjNums;
		    const_cast <Steiner*> (cut) -> setParent (inter);
		    ASSERT (inter->getParent ());
		    boundary << inter;
			}
			ASSERT (boundary. size () == cuts. size ());
		}
		if (boundary. empty ()) 
		{
			cout << "No cuts found" << endl;
  		optimizeSmallSubgraphs (areaRadius_std, false);
      return;
		}
	//qc ();  // Breaks due to transients
	
	
		bool failed = false;
	  {
			VectorOwn<Image> images;  images. reserve (threads_max);
		  Image mainImage (*this);  
			{
				Threads th (boundary. size ());  
				VectorPtr<Steiner> possibleBoundary;  possibleBoundary. reserve (boundary. size ());
				for (const Steiner* cut : boundary)
				{
					ASSERT (cut);
			  	ASSERT (cut->isTransient ());
					auto image = new Image (*this);
					images << image;
					ASSERT (cut->isTransient ());
				#if 1
				  // Use functional ??
				  th << thread (processLargeImage_thread, ref (*image), cut, possibleBoundary);
				#else
				  image->processLarge (cut, possibleBoundary);
				#endif
				  possibleBoundary << cut;
			  }	  
			  // Top subgraph
			  mainImage. processLarge (nullptr, possibleBoundary);
			}
			// Image::apply() can be done by Threads if it is done in the order of cuts and for sibling subtrees ??!
			{
				Progress prog (images. size () + 1);
				Unverbose unv;
				if (! mainImage. apply ())
				  failed = true;
				prog (absCriterion2str ()); 
			  for (const Image* image : images)
			  {
			  	if (! failed && ! var_cast (image) -> apply ())
			  	  failed = true;
					prog (absCriterion2str () + " (approx.)"); 
			  }
			}
		}
	  
	  // DTNode::stable
	  for (DiGraph::Node* node : nodes)
	    static_cast <DTNode*> (node) -> stable = true;

	  // delete transients, DTNode::stable
	  for (const Steiner* st : boundary)
	  {
	  	ASSERT (st->graph == this);
	  	if (st->isTransient ())
	  	{
		  	if (const TreeNode* parent = st->getParent ())
		  		const_static_cast <DTNode*> (parent) -> stable = false;
		 	  delayDeleteRetainArcs (const_cast <Steiner*> (st));
		 	}
		 	else
		  	const_cast <Steiner*> (st) -> stable = false;
	 	}
	  toDelete. deleteData ();
	  
	  setPredictionAbsCriterion ();
	  qc ();
	  
	  if (failed)
	  {
      cerr << "*** OUT OF MEMORY ***" << endl;
      ASSERT (threads_max > 1);
      threads_max--;
      if (threads_max == 1)
      {
        optimizeSmallSubgraphs (areaRadius_std, false);
        break;
      }
	  }
	  else
	  {
  	  cerr << "Optimizing cut nodes ..." << endl;
  	  optimizeSmallSubgraphs (areaRadius_std, true);  // PAR
  	  qc ();
  	  qcPredictionAbsCriterion ();
  	  break;
  	}
  }
}



void DistTree::optimizeSmallSubgraphs (uint areaRadius,
	                                     bool unstableOnly)
{
  ASSERT (areaRadius >= 1);
  
  if (! unstableOnly)
	  for (DiGraph::Node* node : nodes)
	    static_cast <DTNode*> (node) -> stable = false;

  // sort nodes: cf. optimizeLenNode() ??

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
      optimizeSmallSubgraph (center, areaRadius);
    }
    {
      ostringstream oss;
      const ONumber on (oss, criterionDecimals, false);  // PAR
      oss << stables << '/' << steiners << ' ' << absCriterion;
      prog (oss. str ());  
    }
  }

#if 0
  // Almost no improvement
  optimizeLenArc (); 
  optimizeLenNode ();  
  if (verbose (1))
    reportErrors (cout);
#endif
}



void DistTree::optimizeSmallSubgraph (const DTNode* center,
                                      uint areaRadius)
{ 
#ifndef NDEBUG
  const Real absCriterion_old = absCriterion;  
#endif

	Image image (*this);
  image. processSmall (center, areaRadius);
  EXEC_ASSERT (image. apply ());

  qc ();
  qcPredictionAbsCriterion ();

#ifndef NDEBUG
  if (greaterReal (absCriterion, absCriterion_old, 1e-5))   // PAR
  {
    const ONumber on (cout, criterionDecimals + 3, true);
    cout << absCriterion << " " << absCriterion_old << endl;
    ERROR;
  }
#endif
}



void DistTree::delayDeleteRetainArcs (DTNode* node)
{
	ASSERT (node);
	ASSERT (node->graph == this);

  if (verbose ())
    cout << "To delete: " << node->getName () << endl;
    
  if (const Leaf* leaf = node->asLeaf ())
    name2leaf. erase (leaf->name);

  const VectorPtr<DiGraph::Node> children (node->getChildren ());
  for (const DiGraph::Node* child : children)
    const_static_cast <DTNode*> (child) -> len += (node == root ? NaN : node->len);

	if (const Steiner* st = node->asSteiner ())
    if (const TreeNode* parent_ = st->getParent ())
    {
      const Steiner* parent = static_cast <const DTNode*> (parent_) -> asSteiner ();
      ASSERT (parent);
      ASSERT (parent->graph == this);
      const Vector<uint> lcaObjNums (node->getLcaObjNums ());
      for (const size_t objNum : lcaObjNums)
      {
        Dissim& dissim = dissims [objNum];
        if (! dissim. valid ())
          continue;
        ASSERT (dissim. lca == st);
        dissim. lca = parent;
      }
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
 	Vector<DiGraph::Node*> nodeVec;  nodeVec. reserve (2 * name2leaf. size ());
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
  		  && s->childrenDiscernible ()
  		 )
  	{
  		delayDeleteRetainArcs (const_cast <Steiner*> (s));
  		return true;
  	}
 			
  return false;
}



void DistTree::removeLeaf (Leaf* leaf,
	                         bool optimizeP)
{
  ASSERT (! subDepth);
  ASSERT (leaf);
  
  if (! leaf->graph)
    return;
    
  ASSERT (! dissimDs. get ());
  ASSERT (! dissimAttr);
  
  const TreeNode* parent = leaf->getParent ();
  ASSERT (parent);
    
  leaf->detachChildrenUp ();  
  ASSERT (! leaf->graph);
	detachedLeaves << leaf;  

  name2leaf. erase (leaf->name);
  
  // Clean topology, parent, Leaf::discernible consistency
  ASSERT (! parent->isLeaf ());
  if (const TreeNode* child = parent->isTransient ())
	{
    ASSERT (parent->graph); 
    if (const Leaf* childLeaf = static_cast <const DTNode*> (child) -> asLeaf ())
      var_cast (childLeaf) -> discernible = true;
	  const Steiner* st = static_cast <const DTNode*> (parent) -> asSteiner ();
	  ASSERT (st);
	  parent = parent->getParent ();
    if (! parent)
      parent = root;
		delayDeleteRetainArcs (var_cast (st));
  }  
  else if (! static_cast <const DTNode*> (parent) -> childrenDiscernible () && optimizable ())
  {
    for (const DiGraph::Arc* arc : parent->arcs [false])
    {
      const Leaf* child_ = static_cast <const DTNode*> (arc->node [false]) -> asLeaf ();
      ASSERT (child_);
      var_cast (child_) -> discernible = true;
   	}
	  const Steiner* st = static_cast <const DTNode*> (parent) -> asSteiner ();
	  ASSERT (st);
    const LeafCluster leafCluster (var_cast (st) -> getIndiscernibles ());
    leafCluster2discernibles (leafCluster);
  }

      
  if (optimizable ())
  {
		for (const uint objNum : leaf->pathObjNums)
		{
			Dissim& dissim = dissims [objNum];
			ASSERT (dissim. hasLeaf (leaf));
			if (dissim. mult)
	  	{
	      absCriterion -= dissim. getAbsCriterion ();	      
        mult_sum     -= dissim. mult;
        dissim2_sum  -= dissim. mult * sqr (dissim. target); 
	      dissim. mult = 0;
	    }
    }
	  maximize (absCriterion, 0.0);
	  ASSERT (dissim2_sum >= absCriterion);
	  ASSERT (mult_sum > 0);
	
	  qcPaths (); 
	
	  if (optimizeP)
	  {
	    Unverbose unv;
	    ASSERT (parent);
	    optimizeSmallSubgraph (static_cast <const DTNode*> (parent) -> asSteiner (), 2 * areaRadius_std);  // PAR 
	  }
	}


  toDelete. deleteData (); 
}



void DistTree::reroot (DTNode* underRoot,
                       Real arcLen) 
{
  ASSERT (underRoot);
  ASSERT (& underRoot->getTree () == this);
  ASSERT (underRoot != root);  
  ASSERT (! underRoot->inDiscernible ());

  ASSERT (! isNan (arcLen));
  ASSERT (arcLen >= 0);
  ASSERT (arcLen <= underRoot->len);
  
  DTNode* root_ = const_static_cast<DTNode*> (root);
    
  const Steiner* newRootParent = static_cast <const DTNode*> (underRoot->getParent ()) -> asSteiner ();
  auto newRoot = new Steiner (*this, const_cast <Steiner*> (newRootParent), underRoot->len - arcLen);
  newRoot->pathObjNums = underRoot->pathObjNums;
  underRoot->setParent (newRoot); 
  underRoot->len = arcLen;
  
  newRoot->makeDTRoot ();
  ASSERT (newRoot == root);
  ASSERT (root_ != root);
  
  setLca ();  

  if (root_->isTransient ())
    delayDeleteRetainArcs (root_);

  finishChanges ();

  qcPaths ();  
}



Real DistTree::reroot (bool topological)
{
  DTNode* root_ = const_static_cast<DTNode*> (root);
  
  if (! root_->childrenDiscernible ())
    return 0;
  
  root_->setSubtreeLenUp (topological);
//cout << "Old radius: " << root_->getHeight_ave () << endl;  
  
  DTNode* bestDTNode = nullptr;
  Real bestDTNodeLen_new = NaN;
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
  for (const Dissim& dissim : dissims)
    if (dissim. mult)
      s += dissim. mult * dissim. getResidual ();
    ASSERT (! isNan (s));
  
  return s / mult_sum;
}



Real DistTree::getMinLeafLen () const
{
  Real len_min = INF;
  for (DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <DTNode*> (node);
    if (   (  dtNode->isLeaf () && ! dtNode->inDiscernible ())
        || (! dtNode->isLeaf () && ! dtNode->childrenDiscernible ())
       )
      minimize (len_min, dtNode->len);
  }
  return len_min;
}



Real DistTree::getSqrResidualCorr () const
{
  ASSERT (optimizable ());

  Correlation corr;
  for (const Dissim& dissim : dissims)
    if (dissim. mult)
      corr. add (dissim. target, sqr (dissim. getResidual ()));
  
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
  Tree::LcaBuffer buf;
  for (const Dissim& dissim : dissims)
    if (dissim. mult)
    {
      const Real d = dissim. target;
    //ASSERT (d >= 0);
      const Real dHat = dissim. prediction;
      ASSERT (dHat >= 0);
      if (nullReal (dHat))
      {
        epsilon2_0 += dissim. mult * sqr (d);
        continue;
      }
  
      ASSERT (positive (d));
      const Real a = dissim. mult * sqr (dHat - d) / dHat;
      ASSERT (DM_sp::finite (a));
      
  		const VectorPtr<TreeNode>& path = dissim. getPath (buf);
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
    if (dtNode->inDiscernible ())
      continue;
    dtNode->errorDensity = dtNode->paths ? sqrt (dtNode->errorDensity / (Real) dtNode->paths) : INF;
  }
  
  return epsilon2_0;
}



Dataset DistTree::getLeafErrorDataset () const
{
	ASSERT (optimizable ());
	
	Dataset ds;
  ds. objs. reserve (name2leaf. size () / 10 + 1);  // PAR
  auto criterionAttr = new PositiveAttr1 ("leaf_error", ds);  
  for (const auto& it : name2leaf)
  {
  	const Leaf* leaf = it. second;
  	ASSERT (leaf);
    if (leaf->graph)
    {
      const Real err = leaf->getRelCriterion ();
      ASSERT (err >= 0);
  	  const size_t index = ds. appendObj (it. first);
  	  (*criterionAttr) [index] = err;
    }
  }

  return ds;
}



namespace
{
  
bool leafRelCriterionStrictlyGreater (const Leaf* a,
                                      const Leaf* b)
{
  ASSERT (a);
  ASSERT (b);
  return a->getRelCriterion () > b->getRelCriterion ();
}

}



VectorPtr<Leaf> DistTree::findCriterionOutliers (Real outlier_EValue_max,
                                                 Real &outlier_min) const
{
  const Dataset ds (getLeafErrorDataset ());
  ASSERT (ds. attrs. size () == 1)

  const PositiveAttr1* criterionAttr = ds. attrs. front () -> asPositiveAttr1 ();
  ASSERT (criterionAttr);
  
  const Sample sample (ds);

  Normal distr;  // Beta1 ??
  outlier_min = criterionAttr->distr2outlier (sample, distr, true, outlier_EValue_max);

  VectorPtr<Leaf> res;
  if (! isNan (outlier_min))
    for (const auto& it : name2leaf)
    {
	  	const Leaf* leaf = it. second;
      if (leaf->graph)
        if (geReal (leaf->getRelCriterion (), outlier_min))
          res << leaf;
    }
    
  res. sort (leafRelCriterionStrictlyGreater);
  
  // Find a maximal "independent" subset ??
      
	return res;
}



namespace
{
	

struct RequestCandidate
{
	const Leaf* leaf1 {nullptr};
	const Leaf* leaf2 {nullptr};
	bool inDissims {false};
	
	RequestCandidate (const Leaf* leaf1_arg,
	                  const Leaf* leaf2_arg)
	  : leaf1 (leaf1_arg)
	  , leaf2 (leaf2_arg)
	  { ASSERT (leaf1);
	  	ASSERT (leaf2);
	  	ASSERT (leaf1 != leaf2);
			if (leaf1->name > leaf2->name)
		  	swap (leaf1, leaf2);	
	  }
	RequestCandidate ()
	  {}
	  
	bool operator< (const RequestCandidate &other) const
    { return Pair<const Leaf*> (leaf1, leaf2) < Pair<const Leaf*> (other. leaf1, other. leaf2); }
  bool operator== (const RequestCandidate &other) const
    { return Pair<const Leaf*> (leaf1, leaf2) == Pair<const Leaf*> (other. leaf1, other. leaf2); }
};



void hybrid2requests (size_t from, 
                      size_t to, 
                      Vector<RequestCandidate> &res, 
                      const VectorPtr<Leaf> &leaves)
{
  ASSERT (from <= to);
  ASSERT (to <= leaves. size ());
  ASSERT (res. empty ());
  
  Tree::LcaBuffer buf;
  FOR_START (size_t, i, from, to)
  {
  	const Leaf* child = leaves [i];
  	RequestCandidate req;
  	Real hybridness_tree = hybridness_min;
  	for (const Neighbor& badNeighbor1 : child->badNeighbors)
	  	for (const Neighbor& badNeighbor2 : child->badNeighbors)
	  	{
	  		if (& badNeighbor1 == & badNeighbor2)
	  			break;
		    const Tree::TreeNode* lca_ = nullptr;
		    const VectorPtr<Tree::TreeNode>& path = Tree::getPath ( badNeighbor1. leaf
            																						      , badNeighbor2. leaf
            																                  , nullptr
            																                  , lca_
            																                  , buf
            																                  );
		    ASSERT (! path. empty ());
		    const Real hybridness = DistTree::path2prediction (path) / (badNeighbor1. target + badNeighbor2. target);
		    if (maximize (hybridness_tree, hybridness))
		    {
		    	req. leaf1 = badNeighbor1. leaf;
		    	req. leaf2 = badNeighbor2. leaf;
		    }
	    }
	  if (req. leaf1)
	  	res << req;
	}
}



void dissim2hybridness (size_t from,
                        size_t to,
                        Notype&,
                        const Vector<Dissim> &dissims)
{
  ASSERT (from <= to);
  ASSERT (to <= dissims. size ());
  
  FOR_START (size_t, objNum, from, to)
	{
		const Dissim& dissim = dissims [objNum];
		if (! dissim. mult)
			continue;
		ASSERT (dissim. target > 0);
  	for (const Neighbor& neighbor : dissim. leaf1->badNeighbors)
  	{
    	Leaf* child = const_cast <Leaf*> (neighbor. leaf);
    	if (const Neighbor* otherNeighbor = child->findNeighbor (dissim. leaf2))
    	{
    	  const Real hybridness = dissim. target / (neighbor. target + otherNeighbor->target);
    	  lock_guard<mutex> lg {child->mtx};
	    	if (maximize (child->hybridness, hybridness))
	    		child->hybridParentsDissimObjNum = (uint) objNum;
	    }
    }
  }
}


}  // namespace



Vector<TriangleParentPair> DistTree::findHybrids (Real dissimOutlierEValue_max,
	                                                Vector<Pair<const Leaf*>> *dissimRequests) const
{
  ASSERT (! subDepth);
	ASSERT (optimizable ());
	ASSERT (dissimOutlierEValue_max > 0);
	ASSERT (hybridness_min > 1);


  Vector<TriangleParentPair> triangleParentPairs;


	constexpr Real hybridness_min_init = 1;  // PAR		
	ASSERT (hybridness_min > hybridness_min_init);
	

	Real outlier_min = NaN;
	{
	  // dissim.ssize() = 35M => 80 sec. 
		Dataset ds;
	  ds. objs. reserve (dissims. size ());  
	  auto criterionAttr = new PositiveAttr1 ("dissim_error", ds);  
	  for (const Dissim& dissim : dissims)  	
	  	if (dissim. mult)
		  {
	      const Real err = dissim. getAbsCriterion ();
	      ASSERT (err >= 0);
	  	  const size_t index = ds. appendObj ();
	  	  (*criterionAttr) [index] = err;
		  }
	  const Sample sample (ds);	
	  Normal distr;  // Beta1 ??
	  // Time: O(p log(p))
	  outlier_min = criterionAttr->distr2outlier (sample, distr, true, dissimOutlierEValue_max);
	}
	if (isNan (outlier_min))
		return triangleParentPairs;
	if (verbose ())
	  cout << "outlier_min = " << outlier_min << endl;  
	  
  for (const auto& it : name2leaf)
  {
  	Leaf* leaf = const_cast <Leaf*> (it. second);
  	ASSERT (leaf);
  	leaf->badNeighbors. clear ();
  	leaf->badNeighbors. reserve (name2leaf. size () / 100);  // PAR
  	leaf->hybridness = hybridness_min_init;  
  	leaf->hybridParentsDissimObjNum = dissims_max;
  	leaf->badCriterion = 0;
  }
	
	// Leaf::badNeighbors
  for (const Dissim& dissim : dissims)  	
  	if (   dissim. mult
  		  && dissim. getAbsCriterion () >= outlier_min
  		 )
    {
      var_cast (dissim. leaf1) -> badCriterion += dissim. getAbsCriterion ();
      var_cast (dissim. leaf2) -> badCriterion += dissim. getAbsCriterion ();
      if (dissim. target < dissim. prediction)  // <= hybrid
  		{
      	var_cast (dissim. leaf1) -> badNeighbors << Neighbor (dissim. leaf2, dissim. target);
      	var_cast (dissim. leaf2) -> badNeighbors << Neighbor (dissim. leaf1, dissim. target);
      }
    }
  for (const auto& it : name2leaf)
  {
  	Leaf* leaf = const_cast <Leaf*> (it. second);
  	leaf->badNeighbors. sort ();
  	ASSERT (leaf->badNeighbors. isUniq ());
  }

  // Leaf::{hybridness,hybridParentsDissimObjNum}
  {
	  vector<Notype> notypes;
    arrayThreads (dissim2hybridness, dissims. size (), notypes, cref (dissims));
  }

  {
	  Vector<Triangle> triangles;  
	  for (const auto& it : name2leaf)
	  {
	  	const Leaf* child = it. second;
	  	if (! child->graph)
	  		continue;
	  	ASSERT (child->hybridness >= hybridness_min_init);
	  	if (child->hybridness < hybridness_min)  
	  		continue;
			ASSERT (child->hybridParentsDissimObjNum < DistTree::dissims_max);
			const Dissim& dissim = dissims [child->hybridParentsDissimObjNum];
			const Real parent1_dissim = checkPtr<const Neighbor> (child->findNeighbor (dissim. leaf1)) -> target;
			const Real parent2_dissim = checkPtr<const Neighbor> (child->findNeighbor (dissim. leaf2)) -> target;
		  triangles << Triangle (child, child->hybridness, dissim. leaf1, dissim. leaf2, parent1_dissim, parent2_dissim);
	  }
	  if (triangles. empty ())
		  return triangleParentPairs;
	  	  	
	  ASSERT (triangleParentPairs. empty ());
	  for (const Triangle& tr : triangles)
	  {
	  	tr. qc ();
	  	triangleParentPairs << TriangleParentPair (tr. parents [0]. leaf, tr. parents [1]. leaf, tr. parentsDissim ());
	  }
	}
  triangleParentPairs. sort ();
  triangleParentPairs. uniq ();  

  for (TriangleParentPair& tpp : triangleParentPairs)
  	tpp. setTriangles (*this);
  triangleParentPairs. sort (TriangleParentPair::compareHybridness);  // PAR

  Set<const Leaf*> hybrids;
  {
	  Set<const Leaf*> nonHybrids;
	  for (TriangleParentPair& tpp : triangleParentPairs)
	  {
	  	tpp. finish (*this, hybrids);
	  	tpp. qc ();
		  const Set<const Leaf*> hybrids_    (tpp. getHybrids (true));
		  const Set<const Leaf*> nonHybrids_ (tpp. getHybrids (false));
		  IMPLY (hybrids_. empty (), nonHybrids_. empty ());
		  if (   hybrids_.    intersects (nonHybrids)
		  	  || nonHybrids_. intersects (hybrids)
		  	 )
				tpp. triangles. clear ();
			else
			{
				hybrids.    insertAll (hybrids_);
				nonHybrids. insertAll (nonHybrids_);
			}
	  }
	}
	
	triangleParentPairs. filterValue ([] (const TriangleParentPair& tpp) { return tpp. triangles. empty () || tpp. dissimError (); });
	if (triangleParentPairs. empty ())
		return triangleParentPairs;
  

  if (dissimRequests)  
  {  
	  Vector<RequestCandidate> requests;  requests. reserve (name2leaf. size () / 10); // PAR
	  {
      // dissims.size() = 35M => 89 sec.  
      VectorPtr<Leaf> leaves;  leaves. reserve (name2leaf. size ());
  	  for (const auto& it : name2leaf)
  	  {
  	  	const Leaf* child = it. second;
  	  	if (! child->graph)
  	  		continue;
  	  	if (child->badNeighbors. size () <= 1)
  	  		continue;
  	  	if (hybrids. contains (child))
  	  		continue;
  	  	leaves << child;
  	  }
  	  leaves. randomOrder ();

      vector<Vector<RequestCandidate>> requests_vec;  
      arrayThreads (hybrid2requests, leaves. size (), requests_vec, cref (leaves));
      for (const auto& requests_ : requests_vec)
        requests << requests_;
  	}
	  if (verbose ())
	    cout << "# Requests: " << requests. size () << endl;  
	  for (RequestCandidate& req : requests)
	  	if (req. leaf1->name > req. leaf2->name)
	  		swap (req. leaf1, req. leaf2);
	  requests. sort ();
	  requests. uniq ();
	  for (const Dissim& dissim : dissims)  	
	  {
	  	const RequestCandidate req (dissim. leaf1, dissim. leaf2);
	  	const size_t index = requests. binSearch (req);
	  	if (index != NO_INDEX)
	  		requests [index]. inDissims = true;
	  }
	  for (const RequestCandidate& req : requests)
	  	if (! req. inDissims)
	  		*dissimRequests << Pair<const Leaf*> (req. leaf1, req. leaf2);
	}
    
  
  for (const auto& it : name2leaf)
  {
  	Leaf* leaf = const_cast <Leaf*> (it. second);
  	ASSERT (leaf);
  	leaf->badNeighbors. wipe ();
  }


  return triangleParentPairs;
}



#if 0
Vector<Triangle> DistTree::findHybrids (Real parentLengthFrac_min,
	                                      Vector<Pair<const Leaf*>> &dissimRequests) const
{
	ASSERT (optimizable ());
	ASSERT (hybridness_min > 1);
	ASSERT (parentLengthFrac_min > 1);

#if 0
  const Real height_min = getAveArcLength () * 0.1;  // PAR
  ASSERT (height_min > 0);
#endif

  Vector<Triangle> triangles;
  Vector<RequestCandidate> requests;  
  for (const auto& it : name2leaf)
  {
  	Leaf* leaf = const_cast <Leaf*> (it. second);
  	ASSERT (leaf);
  	if (! leaf->graph)
  		continue;
  	const Leaf* parent1 = nullptr;
  	const Leaf* parent2 = nullptr;
  	Real parentDissimilarity = NaN;
  	const Real hybridness = leaf->getHybridness (/*height_min,*/ parentLengthFrac_min, parent1, parent2, parentDissimilarity);
  	if (hybridness < hybridness_min)
  		continue;  	
  	if (isNan (hybridness))
  	{
  		if (parent1)
  		  requests << RequestCandidate (parent1, parent2);
  	}
  	else
  	  triangles << Triangle (leaf, hybridness, parent1, parent2, parentDissimilarity);
  }
  
  // Instead: search in dissims ??
  if (! requests. empty ())
  {
	  requests. sort ();
	  requests. uniq ();
	  for (const Dissim& dissim : dissims)  	
	  {
	  	const RequestCandidate req (dissim. leaf1, dissim. leaf2);
	  	const size_t index = requests. binSearch (req);
	  	if (index != NO_INDEX)
	  		requests [index]. inDissims = true;
	  }
	  for (const RequestCandidate& req : requests)
	  	if (! req. inDissims)
	  		dissimRequests << Pair<const Leaf*> (req. leaf1, req. leaf2);
  }    
    
  return triangles;
}
#endif



#if 0
namespace
{

struct NodeDissim
{
  const DTNode* node {nullptr};
  Real dissim_min {INF};
  
  explicit NodeDissim (const DTNode* node_arg)
    : node (node_arg)
    {
      ASSERT (node);
      const DistTree& tree = node->getDistTree ();
      for (const size_t objNum : node->pathObjNums)
      {
        const Dissim& dissim = tree. dissims [objNum];
        if (dissim. valid ())
          minimize (dissim_min, dissim. target);
      }
    }
    
  void print (ostream &os) const
    { ONumber on (os, 6, true);
      os << node->getLcaName () << '\t' << dissim_min << endl; 
    }
};

}



VectorPtr<DTNode> DistTree::findOutlierArcs (Real outlier_EValue_max,
                                             Real &dissimOutlier_min) const
{
  VectorPtr<DTNode> res;  
  
  Vector<NodeDissim> nodeDissims;  nodeDissims. reserve (2 * name2leaf. size () / 10 + 1);  // PAR
    // Distinct NodeDissim::node's
  for (const DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <const DTNode*> (node);
    if (! dtNode->getParent ())
      continue;
    if (dtNode->inDiscernible ())
      continue;
    const NodeDissim nd (dtNode);
    if (nd. dissim_min < INF)
      nodeDissims << nd;
  }

	Dataset ds;
  auto lenAttr = new PositiveAttr1 ("ArcLength", ds);
  ds. objs. reserve (nodeDissims. capacity ());
  for (const NodeDissim& nd : nodeDissims)
  {
	  const size_t index = ds. appendObj ();
	  (*lenAttr) [index] = nd. dissim_min; 
  }
  ASSERT (nodeDissims. size () == ds. objs. size ());
  Sample sample (ds);
  
  Exponential distr;  

  dissimOutlier_min = lenAttr->distr2outlier (sample, distr, true, 1e-6 /*outlier_EValue_max*/);  // PAR

  if (isNan (dissimOutlier_min))
    return res;
    
  // Mixture of 2 types of dissimilarities (> 2 ??)
  sample. mult. setAll (0);
//OFStream f ("arc_len");  
  FFOR (size_t, objNum, nodeDissims. size ())
    if (geReal (nodeDissims [objNum]. dissim_min, dissimOutlier_min))  
    {
      sample. mult [objNum] = 1;
    //nodeDissims [objNum]. print (f);  
    }
  sample. finish ();
  if (sample. mult_sum >= (Real) name2leaf. size () * 0.001)  // PAR
    dissimOutlier_min = lenAttr->distr2outlier (sample, distr, true, outlier_EValue_max);  // PAR ??

  res. reserve (name2leaf. size () / 1000 + 1);  // PAR
  for (const NodeDissim& nd : nodeDissims)
    if (geReal (nd. dissim_min, dissimOutlier_min))  
      res << nd. node;
      
	return res;
}
#endif



VectorPtr<Leaf> DistTree::findDepthOutliers () const
{
  ASSERT (! subDepth);

  const Real outlier_EValue_max = 1e-3 / (Real) name2leaf. size ();  // PAR
  
  VectorPtr<Leaf> outliers;  outliers. reserve (name2leaf. size () / 100 + 1);  
 	VectorPtr<DTNode> descendants;  descendants. reserve (2 * powInt (2, (uint) areaDiameter_std));
  for (const DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <const DTNode*> (node);
    if (dtNode->asLeaf ())
      continue;
    if (! dtNode->childrenDiscernible ())
      continue;
  	descendants. clear ();
  	dtNode->getDescendants (descendants, areaDiameter_std, nullptr);  
  	if (descendants. size () <= areaDiameter_std)
  	  continue;
  	Dataset ds;
    auto lenAttr = new PositiveAttr1 ("PathLength", ds);
    ds. objs. reserve (descendants. capacity ());
  	for (const TreeNode* descendant : descendants)
  	{
  	  // Move the computation of pathLen into getDescendants() ??
  	  Real pathLen = 0;
  	  const TreeNode* node1 = descendant;
  	  while (node1 != node)
  	  {
  	    pathLen += static_cast <const DTNode*> (node1) -> len;
  	    node1 = node1->getParent ();
  	  }
  	  const size_t index = ds. appendObj ();
  	  (*lenAttr) [index] = pathLen; 
  	}
  	ASSERT (descendants. size () == ds. objs. size ());
    Sample sample (ds);    
    Normal distr;   
    const Real outlier_min = lenAttr->distr2outlier (sample, distr, true, outlier_EValue_max); 
    if (isNan (outlier_min))
      continue;
    FFOR (size_t, objNum, ds. objs. size ())
      if (geReal ((*lenAttr) [objNum], outlier_min))
        outliers << descendants [objNum] -> getReprLeaf ();
  }  
  outliers. sort ();
  outliers. uniq ();
  
  return outliers;
}
  


Vector<Pair<const Leaf*>> DistTree::getMissingLeafPairs_ancestors (size_t depth_max,
                                                                   bool refreshDissims) const
{
  ASSERT (! subDepth);

  Vector<Pair<const Leaf*>> pairs;  pairs. reserve (name2leaf. size () * getSparseDissims_size ());
  {
    Progress prog (name2leaf. size (), 1000);  // PAR
    for (const auto& it : name2leaf)
    {
      prog ();
      const Leaf* leaf = it. second;
      ASSERT (leaf);
      if (! leaf->graph)
        continue;
      const VectorPtr<Leaf> matches (leaf->getSparseLeafMatches (depth_max, true, refreshDissims));
      for (const Leaf* match : matches)
    	  if (leaf->getDiscernible () != match->getDiscernible ())
    	  {
      	  Pair<const Leaf*> p (leaf, match);
      	  ASSERT (! p. same ());
      	  if (p. first->name > p. second->name)
      	    p. swap ();
      	  pairs << p;
        }
    }
  }
  pairs. sort ();
  pairs. uniq ();
       
  return pairs;
}



Vector<Pair<const Leaf*>> DistTree::getMissingLeafPairs_subgraphs () const
{
  ASSERT (! subDepth);

  Vector<Pair<const Leaf*>> pairs;  pairs. reserve (name2leaf. size ()); 
  {
    Subgraph subgraph (*this);
    subgraph. reserve (areaRadius_std);
    Progress prog (nodes. size (), 100);  // PAR
    for (const DiGraph::Node* node : nodes)
    {
      prog ();
      const DTNode* center = static_cast <const DTNode*> (node);
      if (center->inDiscernible ())
        continue;
      subgraph. clear ();
      center->getArea (areaRadius_std, subgraph. area, subgraph. boundary);
      subgraph. removeIndiscernibles ();  
      subgraph. finish ();
      for (const TreeNode* node2 : subgraph. boundary)
      {      
        const DTNode* dtNode2 = static_cast <const DTNode*> (node2);
        const Vector<uint>& pathObjNums2 = subgraph. boundary2pathObjNums (dtNode2);
        const_cast <Vector<uint>&> (pathObjNums2). sort ();
        for (const TreeNode* node1 : subgraph. boundary)
        {
          if (node1 == node2)
            break;
          const DTNode* dtNode1 = static_cast <const DTNode*> (node1);
          const Vector<uint>& pathObjNums1 = subgraph. boundary2pathObjNums (dtNode1);
          if (pathObjNums1. intersectsFast_merge (pathObjNums2))
            continue;
     	    Pair<const Leaf*> leafPair (subgraph. getReprLeaf (dtNode1), subgraph. getReprLeaf (dtNode2));
     	    ASSERT (leafPair. first);
     	    ASSERT (leafPair. second);
     	    ASSERT (! leafPair. same ());
     	    if (leafPair. first->name > leafPair. second->name)
     	      leafPair. swap ();
     	    pairs << leafPair;
        }
      }
    }
  }

#if 0
  Real arcLen_outlier_min = NaN;
  const VectorPtr<DTNode> tooLongArcs (findOutlierArcs (0.1, arcLen_outlier_min));  // PAR
  {
    Progress prog (tooLongArcs. size (), 100);  // PAR
    for (const DTNode* node2 : tooLongArcs)
    {
      prog ();
      ASSERT (node2 != root);
      ASSERT (! node2->inDiscernible ());
      const_cast <Vector<uint>&> (node2->pathObjNums). sort ();
      for (const DTNode* node1 : tooLongArcs)
      {
        ASSERT (node1 != root);
        if (node1 == node2)
          break;
        if (node1->pathObjNums. intersectsFast_merge (node2->pathObjNums))
          continue;
        const TreeNode* lca = getLca (node1, node2);
        ASSERT (lca);
        const DTNode* reprNode1 = lca == node1 ? static_cast <const DTNode*> (node1->getParent () -> getDifferentChild (node1)) : node1;
        const DTNode* reprNode2 = lca == node2 ? static_cast <const DTNode*> (node2->getParent () -> getDifferentChild (node2)) : node2;
   	    Pair<const Leaf*> leafPair (reprNode1->getReprLeaf (), reprNode2->getReprLeaf ());
   	    ASSERT (leafPair. first);
   	    ASSERT (leafPair. second);
   	    ASSERT (! leafPair. same ());
   	    if (leafPair. first->name > leafPair. second->name)
   	      leafPair. swap ();
   	    pairs << leafPair;
      }
    }
  }
#endif

  pairs. sort ();
  pairs. uniq ();  
      
  return pairs;
}



Vector<Pair<const Leaf*>> DistTree::leaves2missingLeafPairs (const VectorPtr<Leaf> &leaves) const
{
  ASSERT (! subDepth);

  Vector<Pair<const Leaf*>> pairs;  pairs. reserve (leaves. size () * (leaves. size () - 1) / 2);  
  {
    Progress prog (leaves. size (), 100);  // PAR
    for (const Leaf* leaf2 : leaves)
    {
      prog ();
      for (const Leaf* leaf1 : leaves)
      {
        if (leaf1 == leaf2)
          break;
   	    Pair<const Leaf*> leafPair (leaf1, leaf2);
   	    ASSERT (! leafPair. same ());
   	    if (leafPair. first->name > leafPair. second->name)
   	      leafPair. swap ();
   	    pairs << leafPair;        
      }
    }
  }
  
  pairs. sort ();
  pairs. uniq ();
       
  if (! dissims. empty ())
    pairs. filterValue ([this] (Pair<const Leaf*> p) { const Dissim dissim (p. first, p. second); return dissims. containsFast (dissim); });

  return pairs;
}



#if 0
void DistTree::findTopologicalClusters () 
{
 	for (DiGraph::Node* node : nodes)
 	{
 	  const DTNode* dtNode = static_cast <DTNode*> (node);
 	  if (Leaf* leaf = const_cast <Leaf*> (dtNode->asLeaf ()))
 	    leaf->DisjointCluster::init ();
 	}
  {
    VectorPtr<Tree::TreeNode> area;      area.     reserve (name2leaf. size () * 2); 
    VectorPtr<Tree::TreeNode> boundary;  boundary. reserve (name2leaf. size ());
    for (DiGraph::Node* node : nodes)
    {
      const Leaf* leaf1 = static_cast <DTNode*> (node) -> asLeaf ();
      if (! leaf1)
        continue;
      if (! leaf1->graph)
        continue;
      area. clear ();
      boundary. clear ();
      leaf1->getArea (sparsingDepth, area, boundary);  // PAR
      for (const TreeNode* other : boundary)
        if (const Leaf* leaf2 = static_cast <const DTNode*> (other) -> asLeaf ())
          if (leaf2->graph)
            const_cast <Leaf*> (leaf1) -> merge (* const_cast <Leaf*> (leaf2));
    }
  }

  LeafCluster leafCluster;
 	for (const DiGraph::Node* node : nodes)
 	{
 	  const DTNode* dtNode = static_cast <const DTNode*> (node);
 	  if (const Leaf* leaf = dtNode->asLeaf ())
 	    if (leaf->graph)
        leafCluster [const_cast <Leaf*> (leaf) -> getDisjointCluster ()] << leaf;
 	}
  ASSERT (! leafCluster. empty ());
  
  // ??
  OFStream f ("leaf_clusters");
  for (const auto& it : leafCluster)
  {
    f << static_cast <const Leaf*> (it. first) -> name << '\t' <<  it. second. size () << '\t';
    for (const Leaf* leaf : it. second)
      f << ' ' << leaf->name;
    f << endl;
  }
}
#endif



VectorPtr<DTNode> DistTree::findDepthClusters (size_t clusters_min) const
// --> Tree ??
{
  ASSERT (clusters_min >= 2);

  Vector<TreeNode::NodeDist> nodeHeights;  nodeHeights. reserve (name2leaf. size ());
  root->getSubtreeHeights (nodeHeights);
  nodeHeights. sort (TreeNode::NodeDist::distLess);
  
  // Use outlier analysis ??
  size_t i_best = NO_INDEX;
  Real diff_max = 0;
  FFOR (size_t, i, nodeHeights. size ())
    if (i && maximize (diff_max, nodeHeights [i]. dist - nodeHeights [i - 1]. dist))
      i_best = i;
  ASSERT (i_best != NO_INDEX);
  
  if (nodeHeights. size () - i_best < clusters_min)
    i_best = nodeHeights. size () - clusters_min;
    
  VectorPtr<DTNode> clusters;  clusters. reserve (nodeHeights. size () - i_best + 1);
  FFOR_START (size_t, i, i_best, nodeHeights. size ())
    clusters << static_cast <const DTNode*> (nodeHeights [i]. node);
    
#ifndef NDEBUG
  clusters. sort ();
  for (const DTNode* node : clusters)
    if (const TreeNode* parent = node->getParent ())
      { ASSERT (clusters. containsFast (static_cast <const DTNode*> (parent))); }
#endif
  
  return clusters;
}



void DistTree::findGenogroups (Real genogroup_dist_max) 
{
  ASSERT (genogroup_dist_max > 0);
  ASSERT (genogroup_dist_max < INF);

  // DTNode::DisjointCluster
  VectorPtr<Leaf> leaves;  leaves. reserve (nodes. size ());
 	for (DiGraph::Node* node : nodes)
 	{
 	  DTNode* dtNode = const_static_cast <DTNode*> (node);
    dtNode->DisjointCluster::init ();
 	  if (const Leaf* leaf = dtNode->asLeaf ())
 	    if (leaf->graph)
 	      leaves << leaf;
 	}

  // Optimize using getDistanceHeight() ??!
  Progress prog (leaves. size (), 1000);  // PAR
  Tree::LcaBuffer buf;
  VectorPtr<TreeNode> area;
  VectorPtr<TreeNode> boundary;
  for (const Leaf* leaf : leaves)
  {
    prog ();
    area. clear ();
    boundary. clear ();
    leaf->getDistanceArea (genogroup_dist_max, area, boundary);
    for (const TreeNode* treeNode : boundary)
      if (const Leaf* other = static_cast <const DTNode*> (treeNode) -> asLeaf ())
      {
        ASSERT (other->graph);
        var_cast (leaf) -> DisjointCluster::merge (* var_cast (other));
      }
  }
}



#if 0
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
    const Real residual = fabs (dHat - d);      
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
#endif



void DistTree::saveDissim (ostream &os,
                           bool addPredictionAbsCriterion) const
{
  ASSERT (optimizable ());

  const ONumber on (os, dissimDecimals, true);
  for (const Dissim& dissim : dissims)
    if (dissim. valid ())
    {
      os         << dissim. leaf1->name
         << '\t' << dissim. leaf2->name
         << '\t' << dissim. target;
      if (addPredictionAbsCriterion)
        os << '\t' << dissim. prediction 
           << '\t' << dissim. getAbsCriterion ();
      os << endl;
    }
}




//////////////////////////////////////////////////////////////////////

// NewLeaf::Location

void NewLeaf::Location::qc () const
{
  if (! qc_on)
    return;

  ASSERT (anchor);
  ASSERT (! anchor->inDiscernible ());
  if (anchor == anchor->getTree (). root)
  {
    ASSERT (isNan (anchor->len));
    ASSERT (arcLen == 0);
  }
  else
  {
    IMPLY (absCriterion_leaf > 0 && anchor != anchor->getDistTree (). root, anchor->len > 0);
    ASSERT (anchor->len < INF);
    ASSERT (arcLen >= 0);
    ASSERT (arcLen <= anchor->len);
  }

  ASSERT (leafLen >= 0);
  ASSERT (leafLen < INF);
  
//IMPLY (leafLen == 0 && arcLen == 0, anchor->asLeaf ());
  
  ASSERT (absCriterion_leaf >= 0);  
}



void NewLeaf::Location::setAbsCriterion (const NewLeaf& nl)
{
  ASSERT (leafLen >= 0);
  ASSERT (arcLen  >= 0);

	absCriterion_leaf = 0;
	for (const Leaf2dissim& ld : nl. leaf2dissims)
	  if (ld. mult)
  	  absCriterion_leaf += ld. mult * sqr (ld. getEpsilon (*this));
	ASSERT (absCriterion_leaf >= 0);
}


    

// NewLeaf::Leaf2dissim

NewLeaf::Leaf2dissim::Leaf2dissim (const Leaf* leaf_arg,
                                   Real dissim_arg,
                                   Real mult_arg,
                                   const DTNode* anchor)
: leaf (leaf_arg)
, dissim (dissim_arg)
, mult (isNan (mult_arg) ? dissim2mult (dissim_arg) : mult_arg)
, dist_hat (0)
, leafIsBelow (false)
{ 
  if (! leaf)
  	throw runtime_error ("No other leaf found for a new leaf placement");
//ASSERT (dissim >= 0);
  ASSERT (! isNan (dissim));
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
                  const string &dataDir,
                  const string &name_arg,
                  bool init)
: Named (name_arg)
, tree (tree_arg)
{
  ASSERT (isDirName (dataDir));  
  const string nameDir (dataDir + "/" + name + "/");
  const string dissimFName  (nameDir + "dissim");
  const string leafFName    (nameDir + "leaf");
  const string requestFName (nameDir + "request");
  process (init, dissimFName, leafFName, requestFName);
}



NewLeaf::NewLeaf (const DTNode* dtNode,
                  size_t q_max,
                  Real &nodeAbsCriterion_old)
: Named (checkPtr (dtNode) -> getLcaName ())
, tree  (checkPtr (dtNode) -> getDistTree ())
, node_orig (dtNode)
{
  ASSERT (dtNode);
  ASSERT (dtNode != tree. root);
  ASSERT (q_max);
  ASSERT (node_orig);

  location. anchor = static_cast <const DTNode*> (tree. root);
  location. leafLen = 0;
  location. arcLen = 0;
  
  // Field of DTNode ??
  Vector<Tree::TreeNode::NodeDist> leafDepths;  leafDepths. reserve (tree. name2leaf. size () / 8 + 1);  // PAR
  dtNode->getLeafDepths (leafDepths);
  leafDepths. sort ();
    
  for (const uint objNum : dtNode->pathObjNums)
  {
    const Dissim& dissim = tree. dissims [objNum];
    if (! dissim. valid ())
      continue;
    size_t index = leafDepths. binSearch (Tree::TreeNode::NodeDist {dissim. leaf1, 0});
    const Leaf* leaf = dissim. leaf2;
    if (index == NO_INDEX)
    {
      index = leafDepths. binSearch (Tree::TreeNode::NodeDist {dissim. leaf2, 0});
      leaf = dissim. leaf1;
    }
    ASSERT (index != NO_INDEX);
    Leaf2dissim ld (leaf, dissim. target - leafDepths [index]. dist, dissim. mult, location. anchor);
    ld. treeAbsCriterion = dissim. getAbsCriterion ();
    leaf2dissims << ld;
  }

  if (leaf2dissims. size () > q_max)
  {
    leaf2dissims. sort (Leaf2dissim::multLess);
    leaf2dissims. resize (q_max);
    leaf2dissims. shrink_to_fit ();  
  }

  nodeAbsCriterion_old = 0;
  for (const Leaf2dissim& ld : leaf2dissims)
    nodeAbsCriterion_old += ld. treeAbsCriterion;  
  
  leaf2dissims. sort ();
  optimize ();
}



void NewLeaf::process (bool init,
                       const string &dissimFName,
                       const string &leafFName,
                       const string &requestFName)
{
  ASSERT (! name. empty ());
  ASSERT (! dissimFName. empty ());
  ASSERT (! leafFName. empty ());
  ASSERT (! requestFName. empty ());
  

  if (fileExists (requestFName))
  {
    cout << "File " << strQuote (requestFName) << " exists" << endl;
    return;
  }
    
  location. anchor = static_cast <const DTNode*> (tree. root);
  location. leafLen = 0;
  location. arcLen = 0;

  if (init)
  { 
    if (fileExists (leafFName))
      throw runtime_error ("File " + strQuote (leafFName) + " exists");
    if (fileExists (dissimFName))
      throw runtime_error ("File " + strQuote (dissimFName) + " exists");
    { OFStream of (dissimFName); }  // destructor creates an empty file
  }
  else
  {
    if (! fileExists (dissimFName))
      throw runtime_error ("File " + strQuote (dissimFName) + " does not exist");
    {
      LineInput f (dissimFName);
      string name1, name2, dissimS;
      while (f. nextLine ())
        try
        {
          istringstream iss (f. line);  // slow ??
          iss >> name1 >> name2 >> dissimS;
          ASSERT (iss. eof ());
          ASSERT (name1 < name2);
          if (name1 != name)
          	swap (name1, name2);
          if (name1 != name)
            throw runtime_error (dissimFName + " must contain " + name);
          const Real dissim = str2real (dissimS);
          if (dissim < 0)
            throw runtime_error ("Dissimilarity must be non-negative");
          leaf2dissims << Leaf2dissim (findPtr (tree. name2leaf, name2), dissim, NaN, location. anchor);
        }          
        catch (...)
        {
          cout << f. line << endl;
          throw;
        }
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

  saveLeaf (leafFName);
  saveRequest (requestFName);
}



void NewLeaf::saveRequest (const string &requestFName) const
{	
  VectorPtr<Leaf> requested (location. anchor->getSparseLeafMatches (sparsingDepth, false, false));
  requested. filterValue ([this] (const Leaf* leaf) { const Leaf2dissim ld (leaf); return leaf2dissims. containsFast (ld); });
  
  OFStream f (requestFName);
  for (const Leaf* leaf : requested)
  {
    ASSERT (name != leaf->name);
    const string* n1 = & name;
    const string* n2 = & leaf->name;
    if (*n1 > *n2)
      swap (n1, n2);
    f << *n1 << '\t' << *n2 << endl;
  }
}



void NewLeaf::optimize ()
{
  bool dissimExists = false;
	for (const Leaf2dissim& ld : leaf2dissims)
	  if (ld. dissim != INF)
  	{
  	  dissimExists = true;
  	  if (ld. dissim == 0)  // if all dissim's = INF ??
  	  {
  	    location. anchor = ld. leaf->getDiscernible ();
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
    if (location. anchor->len == 0)
      location. anchor = static_cast <const DTNode*> (location. anchor->getParent ());
  }
  else if (   node_orig  // greedy search
           && location. absCriterion_leaf / location_best. absCriterion_leaf - 1 > 0.01  // PAR
          )
    return;

  if (! location. anchor->childrenDiscernible ())
    return;
    
  const auto& arcs = location. anchor->arcs [false];
	for (const DiGraph::Arc* arc : arcs)
	{
	  const DTNode* child = static_cast <const DTNode*> (arc->node [false]);
    ASSERT (! child->inDiscernible ());
    if (child == node_orig)
      continue;
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
	//ASSERT (ld. dissim > 0);
	  if (ld. mult)
	  {
  	  mult_sum  += ld. mult;
  	  u_avg     += ld. mult * ld. getU ();
  	  delta_avg += ld. mult * ld. getDelta ();
	  }
	}	
	ASSERT (mult_sum > 0);  // <= pre-condition
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
	IMPLY (u_var > 0, mult_sum > 0);

  bool adjusted = ! (u_var > 0);

  // location.arcLen, adjusted
  location. arcLen = 0;
  if (location. anchor == tree. root)
    {	ASSERT (u_var == 0); }
  else
    if (u_var > 0)  // dissimilarity to all leaves may be INF
    {
    	location. arcLen = delta_u_cov / u_var;
    	if (maximize (location. arcLen, 0.0))
    	  adjusted = true;
    	if (minimize (location. arcLen, location. anchor->len))
    	  adjusted = true;
    }

  // location.leafLen, adjusted
	location. leafLen = mult_sum > 0 ? (delta_avg - u_avg * location. arcLen) : 0;
	if (maximize (location. leafLen, 0.0))
	  adjusted = true;
	
	
  // Case 1
  location. setAbsCriterion (*this);	
  if (adjusted)
  {
    Location loc (location);
    // Case 2
    loc. arcLen = 0;
    loc. leafLen = max (0.0, delta_avg);
    loc. setAbsCriterion (*this);	
    if (location. absCriterion_leaf > loc. absCriterion_leaf)
      location = loc;
    // Case 3
    if (location. anchor != tree. root)
    {
      loc. arcLen = location. anchor->len;
      loc. leafLen = max (0.0, delta_avg - loc. arcLen * u_avg);
      loc. setAbsCriterion (*this);	
      if (location. absCriterion_leaf > loc. absCriterion_leaf)
        location = loc;
    }
  }
}



bool NewLeaf::descend (const DTNode* anchor_new)
{
  ASSERT (anchor_new);
  ASSERT (anchor_new->getParent () == location. anchor);
  ASSERT (anchor_new != node_orig);
  
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



void NewLeaf::qc () const
{
  if (! qc_on)
    return;
  Named::qc ();
    
  location. qc ();    
  ASSERT (& location. anchor->getDistTree () == & tree);
}



}





/* TO DO ??

not solved LinearRegression
  
non-stability of results

*/
