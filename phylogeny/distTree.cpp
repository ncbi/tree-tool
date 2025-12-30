// distTree.cpp

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
*   Distance tree
*
*/


#undef NDEBUG

#include "distTree.hpp"

#include "../dm/prediction.hpp"
#include "../dm/optim.hpp"

#include "../common.inc"



namespace DistTree_sp
{


Chronometer chron_getBestChange ("getBestChange");
Chronometer chron_tree2subgraph ("tree2subgraph");
Chronometer chron_subgraphOptimize ("subgraphOptimize");
Chronometer chron_subgraph2tree ("subgraph2tree");




// VarianceType
const StringVector varianceTypeNames {"lin", "sqr", "pow", "exp", "linExp", "none"};
VarianceType varianceType = varianceType_none;
Real variancePower = NaN;
Real variance_min = 0.0;




// DissimParam

void DissimParam::qc () const 
{
  if (! qc_on)
    return;
    
  Root::qc ();

  QC_ASSERT (coeff > 0.0);
  QC_ASSERT (hybridness_min >= 1.0);
  QC_IMPLY (! isNan (boundary), boundary > 0.0);    
}



void DissimParam::transform (Real &dissim) const
{
  ASSERT (dissim >= 0.0);  
  if (power != 1.0)
    dissim = pow (dissim, power);
  dissim *= coeff;     
}




//

#define BAD_CRITERION(func)  \
  { \
    couterr << "!!" #func << endl;  \
    const ONumber on (cout, absCriterionDecimals, true); \
    PRINT (absCriterion); \
    PRINT (absCriterion_old); \
    PRINT (subDepth); \
  }




// Neighbor

struct Neighbor 
{ 
	const Leaf* leaf {nullptr}; 
	  // !nullptr
	size_t dissimType {no_index};
	Real target {NaN};
	  // > 0

	  
  Neighbor (const Leaf* leaf_arg,
            size_t dissimType_arg,
            Real target_arg)
    : leaf (leaf_arg)
    , dissimType (dissimType_arg)
    , target (target_arg)
    { 
      ASSERT (leaf);
      ASSERT (leaf->graph);
      ASSERT (target >= 0.0); 
    }
  Neighbor (const Leaf* leaf_arg,
            size_t dissimType_arg)
    : leaf (leaf_arg)
    , dissimType (dissimType_arg)
    { 
      ASSERT (leaf);
      ASSERT (leaf->graph);
    }

	  
  bool operator< (const Neighbor &other) const
    { 
      LESS_PART (*this, other, dissimType);
      return leaf < other. leaf; 
    }
	bool operator== (const Neighbor &other) const
	  { return    leaf       == other. leaf
	           && dissimType == other. dissimType; 
	  }
};




// Triangle

Triangle::Triangle (const Leaf* child_arg,
                    Real parentDissim_arg /*hybridness_arg*/,
                    const Leaf* parent1,
                    const Leaf* parent2,
                    Real parent1_dissim,
                    Real parent2_dissim,
				  	        size_t dissimType_arg)
: child (child_arg)
, hybridness (parentDissim_arg / (parent1_dissim + parent2_dissim) /*hybridness_arg*/)
, dissimType (dissimType_arg)
{
  ASSERT (parent1);
  ASSERT (parent2);
  parents [0]. leaf = parent1;
  parents [1]. leaf = parent2;
  parents [0]. dissim = parent1_dissim;
  parents [1]. dissim = parent2_dissim;
  if (parents [0]. leaf->getName () > parents [1]. leaf->getName ())
    swap (parents [0], parents [1]);
}



void Triangle::qc () const
{
  if (! qc_on)
    return;
    
  QC_ASSERT (child);
  QC_ASSERT (child->graph);
  QC_ASSERT (child->badCriterion >= 0.0);
  QC_IMPLY (child->good, ! child_hybrid);
  for (const bool i : {false, true})
  {
    const Parent& p = parents [i];
    QC_ASSERT (p. leaf); 
    QC_ASSERT (p. leaf->graph);
    QC_ASSERT (p. leaf->badCriterion >= 0.0);
    QC_ASSERT (child->getDiscernible () != p. leaf->getDiscernible ());
    QC_ASSERT (p. dissim > 0.0);
    QC_IMPLY (p. leaf->good, ! p. hybrid);
  }
  QC_ASSERT (! (   child_hybrid 
                && parents [false]. hybrid 
                && parents [true].  hybrid 
               )
            );
  QC_ASSERT (parents [0]. leaf->getDiscernible () != parents [1]. leaf->getDiscernible ());
  QC_ASSERT (parents [0]. leaf->getName () < parents [1]. leaf->getName ());
  QC_ASSERT (hybridness >= getHybridness_min ());
}



void Triangle::print (ostream &os) const
{ 
  const ONumber on (os, 3, false);  // PAR
  const Parent& p1 = parents [0];
  const Parent& p2 = parents [1];
  os         << child->getName ()
     << '\t' << hybridness 
     << '\t' << p1. leaf->getName ()
     << '\t' << p2. leaf->getName ()
     << '\t' << p1. dissim 
     << '\t' << p2. dissim 
     << '\t' << child_hybrid 
     << '\t' << p1. hybrid 
     << '\t' << p2. hybrid;
  if (dissimType != no_index)
    os << '\t' << dissimType;
  os << endl; 
}



bool Triangle::operator< (const Triangle &other) const
{ 
  LESS_PART (*this, other, child->getName ());
  LESS_PART (*this, other, parents [0]. leaf->getName ());
  LESS_PART (*this, other, parents [1]. leaf->getName ());
  LESS_PART (*this, other, dissimType);
  return false;
}



Real Triangle::getHybridness_min () const
{
  ASSERT (child);
  return child->getDistTree (). dissimParam. hybridness_min;
}



void Triangle::qcMatchHybrids (const VectorPtr<Leaf> &hybrids) const
{
  if (! qc_on)
    return;
  #define QC_MATCH_HYBRIDS(leaf,isHybrid)  if ((isHybrid) != hybrids. containsFast (leaf)) { cout << (leaf)->getName () << endl; print (cout); ERROR; } 
  QC_MATCH_HYBRIDS (child, child_hybrid);
  QC_MATCH_HYBRIDS (parents [0]. leaf, parents [0]. hybrid);
  QC_MATCH_HYBRIDS (parents [1]. leaf, parents [1]. hybrid);
  #undef QC_MATCH_HYBRIDS
}




//

void addTriangle (Vector<Triangle> &triangles,
                  const Leaf* leaf1,
                  const Leaf* leaf2,
                  const Leaf* leaf3,
                  Real target23,
                  Real target13,
                  Real target12,
                  size_t dissimType)
{
  ASSERT (leaf1);
  ASSERT (leaf2);
  ASSERT (leaf3);
  ASSERT (& leaf1->getDistTree () == & leaf2->getDistTree ());
  ASSERT (& leaf1->getDistTree () == & leaf3->getDistTree ());
  
  const Triangle tr1 ( leaf1
                     , target23
                     , leaf2
                     , leaf3
                     , target12
                     , target13
                     , dissimType
                     );
  const Triangle tr2 ( leaf2
                     , target13
                     , leaf1
                     , leaf3
                     , target12
                     , target23
                     , dissimType
                     );
  const Triangle tr3 ( leaf3
                     , target12
                     , leaf1
                     , leaf2
                     , target13
                     , target23
                     , dissimType
                     );
#ifndef NDEBUG
  const size_t size_old = triangles. size ();
#endif
  if (tr1. hybridness >= tr1. getHybridness_min ())  triangles << tr1;
  if (tr2. hybridness >= tr2. getHybridness_min ())  triangles << tr2;
  if (tr3. hybridness >= tr3. getHybridness_min ())  triangles << tr3;
  ASSERT (triangles. size () <= size_old + 1);
}




// TriangleParentPair

void TriangleParentPair::setTriangles (const DistTree &tree) 
{
  for (const bool i : {false, true})
    ASSERT (parents [i]. leaf);
  ASSERT (parents [0]. leaf->getName () < parents [1]. leaf->getName ());
  ASSERT (parentsDissim >= 0.0);
  ASSERT (getDissimParam (). hybridness_min > 1.0);
  ASSERT (triangles. empty ());

  Vector<Neighbor> hybridParents;  hybridParents. reserve (parents [0]. leaf->pathDissimNums. size ());
  for (const uint dissimNum : parents [0]. leaf->pathDissimNums)
  {
    const Dissim& dissim = tree. dissims [dissimNum];
    if (! dissim. validMult ())
      continue;
    if (dissim. type != dissimType)
      continue;
    const Leaf* other = dissim. getOtherLeaf (parents [0]. leaf);
    ASSERT (other->graph);
    hybridParents << Neighbor (other, dissim. type, dissim. target);
  }
  hybridParents. sort ();  // --> unordered_set ??
  ASSERT (hybridParents. isUniq ());

  for (const uint dissimNum : parents [1]. leaf->pathDissimNums)
  {
    const Dissim& dissim = tree. dissims [dissimNum];
    if (! dissim. validMult ())
      continue;
    if (dissim. type != dissimType)
      continue;
    const Leaf* other = dissim. getOtherLeaf (parents [1]. leaf);
    ASSERT (other->graph);
    const size_t i = hybridParents. binSearch (Neighbor (other, dissim. type));
    if (i == no_index)
      continue;
    addTriangle ( triangles
                , other
                , parents [0]. leaf
                , parents [1]. leaf
                , parentsDissim
                , dissim. target
                , hybridParents [i]. target
                , dissim. type
                );
  }
//ASSERT (! triangles. empty ());
}



void TriangleParentPair::triangles2hybridness_ave ()
{
  hybridness_ave = getDissimParam (). hybridness_min;
  if (triangles. empty ())
    return;

  MeanVar hybridness_mv;
  for (const Triangle& tr : triangles)
    hybridness_mv << tr. hybridness;
  hybridness_ave = hybridness_mv. getMean ();
  ASSERT (hybridness_ave >= getDissimParam (). hybridness_min);
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
  ASSERT (triangle_best_index == no_index);

  triangles. sort ();
  triangles. uniq ();

  triangles. filterValue ([&hybrids] (const Triangle &tr) { return hybrids. contains (tr. child); }); 
  if (triangles. empty ())
    return;

  Triangle* best = nullptr;
  Real hybridness = 0.0;
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
    if (   best->parents [! i]. dissim > getDissimParam (). boundary
        && best->parents [  i]. dissim < getDissimParam (). boundary
       ) 
    {
      const Parent& p = parents [i];    
      if (   ! childrenGood () 
          && child_classSize < all * classSizeFrac_max
         )
        setChildrenHybrid ();
      else if (   ! best->parents [i]. leaf->good 
               && (Real) p. classSize < all * classSizeFrac_max
              )
        best->parents [i]. hybrid = true;   
      else if (   ! childrenGood () 
               && ! best->parents [i]. leaf->good
               && child_classSize + (Real) p. classSize < all * classSizeFrac_max
              )
      {
        setChildrenHybrid ();
        best->parents [i]. hybrid = true;   
      }
      return;
    }
  
  {
    const Prob p = child_classSize / all;
    if (   ! childrenGood () 
        && p <= classSizeFrac_max
       )
    {
      // aa bb ab
      setChildrenHybrid ();
      return;
    }
    if (   ! best->parents [0]. leaf->good
        && ! best->parents [1]. leaf->good
        && p >= 1.0 - classSizeFrac_max
       )
    {
      // ab ac aa
      // common(ab,ac) = 1/4 aa <= random halves of aa
      best->parents [0]. hybrid = true;
      best->parents [1]. hybrid = true;
      return;       
    }
  }
  

#if 0
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
    size_t i_bad = no_index;
    Real badCriterion = triangles [0]. child->badCriterion;
    for (const bool i : {false, true})  
      if (maximize (badCriterion, parents [i]. leaf->badCriterion))
        i_bad = (size_t) i;
    if (i_bad == no_index)
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
#endif
}



void TriangleParentPair::qc () const
{
  if (! qc_on)
    return;

  for (const bool i : {false, true})  
    QC_ASSERT (parents [i]. leaf);
  QC_ASSERT (parents [0]. leaf->getName () < parents [1]. leaf->getName ());
  QC_ASSERT (parentsDissim > 0.0);
  
  if (triangles. empty ())
    return;
  
  bool child_hybrid;
  bool first = true;
  for (const Triangle& tr : triangles)
  {
    tr. qc ();
    for (const bool i : {false, true})  
      QC_ASSERT (tr. parents [i]. leaf == parents [i]. leaf);
    QC_ASSERT (eqReal (parentsDissim, tr. parentsDissim ()));
    if (first)
    {
      child_hybrid = tr. child_hybrid;
      first = false;
    }
    else
      QC_ASSERT (child_hybrid == tr. child_hybrid);
    QC_ASSERT (tr. dissimType == dissimType);
  }
  
  QC_ASSERT (& getBest ());

  for (const bool i : {false, true})  
    QC_ASSERT (parents [i]. classSize);

  QC_ASSERT (hybridness_ave >= getDissimParam (). hybridness_min);
}



void TriangleParentPair::print (ostream &os) const
{ 
  const ONumber on (os, 2, false);  // PAR
  const Parent& p1 = parents [0];
  const Parent& p2 = parents [1];
  os         << getBest (). child->getName ()
     << '\t' << p1. leaf->getName ()
     << '\t' << p2. leaf->getName ()
     << '\t' << triangles. size () 
     << '\t' << p1. classSize 
     << '\t' << p2. classSize
     << '\t' << hybridness_ave 
     << '\t' << getBest (). parents [0]. dissim 
     << '\t' << getBest (). parents [1]. dissim
     << '\t' << getBest (). child_hybrid 
     << '\t' << getBest (). parents [0]. hybrid 
     << '\t' << getBest (). parents [1]. hybrid
     << '\t' << (int) dissimType
     << endl;
}



const DissimParam& TriangleParentPair::getDissimParam () const
{
  ASSERT (parents [0]. leaf);
  return parents [0]. leaf -> getDistTree (). dissimParam;
}


bool TriangleParentPair::operator< (const TriangleParentPair &other) const
{ 
  LESS_PART (*this, other, parents [0]. leaf);
  LESS_PART (*this, other, parents [1]. leaf);
  LESS_PART (*this, other, dissimType);
  return false;
}



bool TriangleParentPair::compareHybridness (const TriangleParentPair &tpp1,
                                            const TriangleParentPair &tpp2)
{ 
  LESS_PART (tpp2, tpp1, hybridness_ave); 
  LESS_PART (tpp1, tpp2, parents [0]. leaf->getName ());
  LESS_PART (tpp1, tpp2, parents [1]. leaf->getName ());
  return false;
}



bool TriangleParentPair::dissimError () const  
{ 
  if (   getBest (). parent_dissim_ratio () < 0.25  // PAR
	    && hybridness_ave < 1.25  // PAR
	   )
	  return true;
	for (const bool i : {false, true})
	  if (getDissimParam (). at_boundary (getBest (). parents [i]. dissim))
		  return true; 
  return false;
}


size_t TriangleParentPair::child_parent2parents (const DistTree &tree,
                                                 const Leaf* child,
                                                 const Leaf* parent,
                                                 Real parentDissim) const
{
  ASSERT (child);
  ASSERT (parent);
  ASSERT (child != parent);
  ASSERT (parentDissim > 0.0);
  ASSERT (getDissimParam (). hybridness_min > 1.0);
    
  Vector<Neighbor> hybridParents;  hybridParents. reserve (child->pathDissimNums. size ());
  for (const uint dissimNum : child->pathDissimNums)
  {
    const Dissim& dissim = tree. dissims [dissimNum];
    if (! dissim. validMult ())
      continue;
    if (dissim. type != dissimType)
      continue;
    const Leaf* other = dissim. getOtherLeaf (child);
    ASSERT (other->graph);
    hybridParents << Neighbor (other, dissim. type, dissim. target);
  }
  hybridParents. sort ();
  ASSERT (hybridParents. isUniq ());

  size_t n = 0;
  for (const uint dissimNum : parent->pathDissimNums)
  {
    const Dissim& dissim = tree. dissims [dissimNum];
    if (! dissim. validMult ())
      continue;
    if (dissim. type != dissimType)
      continue;
    const Leaf* otherParent = dissim. getOtherLeaf (parent);
    ASSERT (otherParent->graph);
    const size_t i = hybridParents. binSearch (Neighbor (otherParent, dissim. type));
    if (i == no_index)
      continue;
    const Real hybridness = dissim. target / (parentDissim + hybridParents [i]. target);
    if (hybridness >= getDissimParam (). hybridness_min)  // otherParent's may be too diverse to be all hybrid
      n++;
  }

  return n;
}



bool TriangleParentPair::childrenGood () const
{ 
  for (const Triangle& tr : triangles)
	  if (tr. child->good)
	    return true;
	return false;
}



void TriangleParentPair::setChildrenHybrid ()
{ 
  ASSERT (! childrenGood ());
  for (Triangle& tr : triangles)
	  tr. child_hybrid = true;
}




//

static const string lenS               ("len");
static const string err_densityS       ("err_density");
static const string normCriterionS     ("norm_criterion");
static const string deformationS       ("deformation");
static const string deformation_criterionS ("deformation_criterion");



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
    QC_IMPLY (getParent (), ! isNan (len));
    QC_IMPLY (! getParent (), isNan (len));
    if (! childrenDiscernible ())
    {
    //QC_ASSERT (! DistTree_sp::variance_min);
      QC_ASSERT (! inDiscernible ());
      for (const DiGraph::Arc* arc : arcs [false])
      {
        const Leaf* leaf = static_cast <const DTNode*> (arc->node [false]) -> asLeaf ();
        QC_ASSERT (leaf);
        QC_ASSERT (leaf->inDiscernible ());
        QC_ASSERT (leaf->good == static_cast <const DTNode*> (arcs [false]. front () -> node [false]) -> asLeaf () -> good);
      }
    }
  }
  
  QC_IMPLY (! isNan (len), len >= 0.0);    
  QC_ASSERT (! contains (name, '\''));
}



void DTNode::saveContent (ostream& os) const
{ 
  {
    const ONumber oLen (os, dissimDecimals, true);
    os << lenS << "=" << len;
  }
  
  const ONumber oNum (os, relCriterionDecimals, true);  // PAR

  if (! isNan (errorDensity))
    os << "  " << err_densityS << "=" << errorDensity;
    
	if (maxDeformationDissimNum != dissims_max)
	  os << "  " << getDeformationS ();
  else if (const DistTree::DeformationPair* dp = findPtr (getDistTree (). node2deformationPair, this))
    os << "  " << deformationS << '=' << dp->leafName1 
                               << ':' << dp->leafName2 
       << "  " << deformation_criterionS << '=' << dp->deformation;
}



const DistTree& DTNode::getDistTree () const
{
  return static_cast <const DistTree&> (getTree ());
}



const Leaf* DTNode::inDiscernible () const
{ 
  if (const Leaf* g = asLeaf ())
  {
    if (g->discernible)
      return nullptr;
    else
    {
    //ASSERT (! DistTree_sp::variance_min);
      return g;
    }
  }
  return nullptr;
}



const DTNode* DTNode::getDiscernible () const
{ 
  if (const Leaf* leaf = inDiscernible ())
    return static_cast <const DTNode*> (leaf->getParent ());
  return this;
}



Prob DTNode::getArcExistence () const
{
  WeightedMeanVar mv;
  for (const uint dissimNum : pathDissimNums)
  {
    const Dissim& dissim = getDistTree (). dissims [dissimNum];
    if (! dissim. validMult ())
      continue;
    mv. add (dissim. target - (dissim. prediction - len) > 0.0 ? 1.0 : 0.0, dissim. mult);
  }
  return mv. getMean ();
}



Real DTNode::getDeformation () const
{ 
  if (maxDeformationDissimNum == dissims_max)
    return NaN;
  const Real deformation = getDistTree (). dissims [maxDeformationDissimNum]. getDeformation ();
  ASSERT (deformation >= 0.0);
  return deformation;
}



string DTNode::getDeformationS () const
{
  ASSERT (maxDeformationDissimNum != dissims_max);
  
  ostringstream os;
  const Dissim& dissim = getDistTree (). dissims [maxDeformationDissimNum];
  os << deformationS << "=" << dissim. leaf1->getName () 
                     << ':' << dissim. leaf2->getName ()
     << "  " << deformation_criterionS << "=" << dissim. getDeformation ();
  return os. str ();
}



void DTNode::saveFeatureTree (ostream &os,
                              bool withTime,
                              size_t offset) const
{
  FOR (size_t, i, offset)
    os << ' ';
  const ONumber on (os, dissimDecimals, true);
  string s = (asLeaf () ? "s" : "") + getName ();
  os << s << ": ";
  if (withTime)
    os << " t=" << (isNan (len) ? 0 : len);
  os << "  C=0  dC=+0-0" << endl;
  if (asLeaf ())
  {
    FOR (size_t, i, offset + 2)
      os << ' ';
    os << 'g' << getName () << ": " << endl;
  }
  else
    for (const DiGraph::Arc* arc : arcs [false])
      static_cast <const DTNode*> (arc->node [false]) -> saveFeatureTree (os, withTime, offset + 2);
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
        subtreeLen1. addValue (1.0);
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
      lenDelta = len / 2.0;
      globalLen. add (      subtreeLen);
      globalLen. add (parentSubtreeLen);
      globalLen. addValue (1);
    }
    else
    {
      // Optimal
      lenDelta =   len / 2.0 
                 + (  (parentSubtreeLen. getMean () + parentSubtreeLen. weights)
                    - (      subtreeLen. getMean () +       subtreeLen. weights)
                   ) / 4.0;
      maximize (lenDelta, 0.0);
      minimize (lenDelta, len);
      globalLen. add (      subtreeLen. getMean () + lenDelta,               subtreeLen. weights + lenDelta);
      globalLen. add (parentSubtreeLen. getMean () + (len - lenDelta), parentSubtreeLen. weights + (len - lenDelta));
    }
    ASSERT (! isNan (lenDelta));
    
  //constexpr Real zScore = 3.0;  // PAR
    if (! inDiscernible ())
      if (   ! bestDTNode 
          || bestGlobalLen. /*getOutlier_min (zScore)*/ getMean () > globalLen. /*getOutlier_min (zScore)*/ getMean ()
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
    var_cast (child) -> setSubtreeLeaves ();
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



Vector<uint> DTNode::getLcaDissimNums () 
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
    Vector<uint>& childPathObjNums = const_static_cast <DTNode*> (children [i]) -> pathDissimNums;
  #if 0
    // Faster for small n
    for (const size_t dissimNum : childPathObjNums)
      if (childDissims [dissimNum])
        lcaObjNums << dissimNum;
      else
        childDissims [dissimNum] = true;
  #else
    childPathObjNums. sort ();
    FOR (size_t, j, i)
      for (const uint dissimNum : static_cast <const DTNode*> (children [j]) -> pathDissimNums)
        if (childPathObjNums. containsFast (dissimNum))
          lcaObjNums << dissimNum;
  #endif
  }
        
  return lcaObjNums;
}



void DTNode::setErrorDensity (Real absCriterion_ave)
{
  ASSERT (absCriterion_ave >= 0.0);
  
  Real a = 0.0;
  Real b = 0.0;
  for (const uint dissimNum : pathDissimNums)
  {
    const Dissim& dissim = getDistTree (). dissims [dissimNum];
    if (! dissim. validMult ())
      continue;
    const Real criterion = dissim. getAbsCriterion ();
    ASSERT (criterion >= 0.0);
    ASSERT (criterion < inf);
    ASSERT (dissim. prediction >= 0.0);
    if (! dissim. prediction)
      continue;
    a += (criterion / absCriterion_ave - 1.0) / dissim. prediction;
    b += 1.0 / sqr (dissim. prediction);
  }
  
  errorDensity = a / sqrt (2.0 * b);
}



VectorPtr<Leaf> DTNode::getSparseLeafMatches (const string &targetName,
                                              size_t depth_max,
                                              bool subtractDissims,
                                              bool refreshDissims) const
{
  IMPLY (depth_max, depth_max >= sparsingDepth);
  IMPLY (subtractDissims, asLeaf ());
  IMPLY (subtractDissims, getDistTree (). optimizable ());
  
//hash<string> hash_str;
  const ulong seed = str_hash (targetName) + 1;

  unordered_map <const Steiner* /*lca*/, VectorPtr<Leaf>> lca2leaves;  lca2leaves. rehash (pathDissimNums. size ());  
  if (subtractDissims)
  {
    // Time: ~ O(p/n) ~ O(log(n))
    for (const uint dissimNum : pathDissimNums)
    {
      const Dissim& dissim = getDistTree (). dissims [dissimNum];
      if (! dissim. validMult ())
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
  size_t height = 0;  // From ancestor to *this
  const DTNode* ancestor = this;
  const DTNode* ancestor_prev = nullptr;
  constexpr size_t depth_min = 3;  // PAR
  constexpr Prob p_delta = 1.0 - 1e-2;  // PAR
  Prob p = 1.0;
  // Time: ~ O(log^2(n))
  while (ancestor) 
  {
    descendants. clear ();
    size_t searchDepth = sparsingDepth;
    {
      if (depth_max && height > depth_max)
      {
        const size_t dec = height - depth_max;
        if (searchDepth > dec)
          searchDepth -= dec;
        else
          searchDepth = 1;
      }
      maximize<size_t> (searchDepth, depth_min);  
      ASSERT (searchDepth);
      ancestor->getDescendants (descendants, searchDepth, ancestor_prev);  
    }
    ASSERT (! descendants. empty ());
    ASSERT (searchDepth >= depth_min);
    if (searchDepth == depth_min)
    {
      descendants. randomOrder ();
      p *= p_delta;
      Rand rand (seed_global);
      FFOR (size_t, i, descendants. size ())
        if (rand. getProb () >= p)
        {
          descendants. eraseMany (i, descendants. size ());
          break;
        }
    }
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
        const Leaf* repr = descendant->getReprLeaf (seed);
        ASSERT (repr);
        if (otherLeaves && otherLeaves->containsFast (repr))
          continue;
        matches << repr;
        added++;
      }
    }
    ancestor_prev = ancestor;
    ancestor = static_cast <const DTNode*> (ancestor->getParent ());
    height++;
  }
  matches. sort ();
  matches. uniq ();

  return matches;
}




// DTNode::ClosestLeaf

void DTNode::ClosestLeaf::qc () const
{
  if (! qc_on)
    return;
  QC_ASSERT (dc);
  QC_ASSERT (dist >= 0.0);
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
    
  QC_ASSERT (! isLeaf ()); 
}



void Steiner::saveContent (ostream& os) const 
{ 
  DTNode::saveContent (os);
  if (! name. empty ())
    os << "  name='" << name << '\'';
}



const Leaf* Steiner::getReprLeaf (ulong seed) const
{ 
  const size_t n = arcs [false]. size ();
  ASSERT (n);
  Rand& rand = getDistTree (). rand;
  if (seed)
    rand. setSeed (seed);
  const size_t index = rand. get (n);
  return static_cast <const DTNode*> (arcs [false]. at (index) -> node [false]) -> getReprLeaf (0);
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
#if 0
  cout << this << " setSubtreeLenUp: ";
  subtreeLen. saveText (cout); 
#endif
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

  for (DiGraph::Arc* arc : arcs [false])
  {
    DTNode* child = static_cast <DTNode*> (arc->node [false]);
    child->setLca ();  // DFS
    DisjointCluster::merge (*child);
    DTNode* cluster = static_cast <DTNode*> (getDisjointCluster ());
    ASSERT (cluster);
    cluster->tarjanLca = this;
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
    var_cast (oldParent -> asSteiner ()) -> reverseParent (target, this);
  }
  
  if (child)
  {
    setParent (child);
    // Arc-specific data
    len            = child->len;
    pathDissimNums = child->pathDissimNums;
    errorDensity   = child->errorDensity;
  }
}



void Steiner::makeRoot (Steiner* ancestor2descendant)
{
  ASSERT (ancestor2descendant);
  ASSERT (descendantOf (ancestor2descendant));

  const DTNode* parent_new = static_cast <const DTNode*> (ancestor2descendant->getParent ());
  // Arc-specific data
  const auto len_new          = ancestor2descendant->len;
  const auto pathObjNums_new  = ancestor2descendant->pathDissimNums;
  const auto errorDensity_new = ancestor2descendant->errorDensity;
  
  reverseParent (ancestor2descendant, nullptr);

  setParent (var_cast (parent_new));
  // Arc-specific data
  len             = len_new;
  pathDissimNums  = pathObjNums_new;
  errorDensity    = errorDensity_new;
}




const Steiner* Steiner::makeDTRoot ()
{ 
  const Steiner* root = static_cast <const DTNode*> (getTree (). root) -> asSteiner ();
  ASSERT (root);
  makeRoot (var_cast (root)); 
  ASSERT (this == getTree (). root);
  
  return root;
}



Cluster2Leaves Steiner::getIndiscernibles ()
{
  ASSERT (! getDistTree (). subDepth);
  ASSERT (getDistTree (). optimizable ());
  ASSERT (! childrenDiscernible ());
  
  Cluster2Leaves cluster2leaves;  cluster2leaves. rehash (arcs [false]. size ());

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

  for (const Leaf* leaf : children)
    for (const uint dissimNum : leaf->pathDissimNums)
    {
      const Dissim& dissim = getDistTree (). dissims [dissimNum];
      if (   dissim. valid ()
          && dissim. target <= 0.0
        //&& dissim. mult == inf
          && children. containsFast (dissim. getOtherLeaf (leaf))
         )
        var_cast (dissim. leaf1) -> DisjointCluster::merge (* var_cast (dissim. leaf2));
    }

  for (const Leaf* leaf : children)
    cluster2leaves [var_cast (leaf) -> DisjointCluster::getDisjointCluster ()] << leaf;
  
  return cluster2leaves;
}



Vector<DTNode::ClosestLeaf> Steiner::findGenogroups (Real genogroup_dist_max) 
{
  ASSERT (genogroup_dist_max > 0.0);
  ASSERT (genogroup_dist_max < inf);


  Vector<ClosestLeaf> res;
  
  Vector<ClosestLeaf> vec;  vec. reserve (16); // PAR
  for (const DiGraph::Arc* arc : arcs [false])
    vec << static_cast <DTNode*> (arc->node [false]) -> findGenogroups (genogroup_dist_max);
  if (vec. empty ())
    return res;

  {  
    ClosestLeaf* cl_min = nullptr;
    Real dist_min = inf;
    for (ClosestLeaf &cl : vec)
    {
      cl. qc ();
      ASSERT (cl. dist <= genogroup_dist_max);
      if (minimize (dist_min, cl. dist))
        cl_min = & cl;
    }
    ASSERT (cl_min);
    for (ClosestLeaf &cl : vec)
      if (   cl_min != & cl
          && cl_min->dist + cl. dist <= genogroup_dist_max
         )
        var_cast (cl_min->dc) -> DisjointCluster::merge (* var_cast (cl. dc));
  }
        
  unordered_map<const DisjointCluster*,Real/*dist*/> dc2dist;
  dc2dist. rehash (vec. size ());
  for (ClosestLeaf &cl : vec)
  {
    const DisjointCluster* dc = var_cast (cl. dc) -> getDisjointCluster ();
    const auto it = dc2dist. find (dc);
    if (it == dc2dist. cend ())
      dc2dist [dc] = cl. dist;
    else
      minimize (it->second, cl. dist);
  }
  
  ASSERT (res. empty ());
  res. reserve (dc2dist. size ()); 
  for (const auto& it : dc2dist)
  {
    const Real dist = it. second + len;
    if (dist <= genogroup_dist_max)
      res << ClosestLeaf {it. first, dist};
  }
   
  return res;
}



void Steiner::copySubtree (Steiner &to,
                           Real len_coeff) const
{
  ASSERT (len_coeff > 0.0);
  ASSERT (DM_sp::finite (len_coeff));


  DistTree& tree_to = var_cast (to. getDistTree ());
  ASSERT (& tree_to != & getDistTree ());
      
  for (const DiGraph::Arc* arc : arcs [false])
  {
    ASSERT (arc);
    const DTNode* child = static_cast <const DTNode*> (arc->node [false]);
    ASSERT (child);
    const Real len_to = child->len * len_coeff;
    // child_to
    DTNode* child_to = nullptr;
    if (const Leaf* leaf = child->asLeaf ())
    {
      auto leaf_to = new Leaf (tree_to, & to, len_to, leaf->getName ());
      leaf_to->normCriterion = leaf->normCriterion;
      child_to = leaf_to;
    }
    else
    {
      auto steiner_to = new Steiner (tree_to, & to, len_to);
      ASSERT (child->asSteiner ());
      child->asSteiner () -> copySubtree (*steiner_to, len_coeff);
      child_to = steiner_to;
    }      
    ASSERT (child_to);
    child_to->errorDensity = child->errorDensity;
    // child->maxDeformationDissimNum ??
  }  
}



void Steiner::replaceSubtree (const DistTree &from) 
{
  ASSERT (& from != & getDistTree ());  
  ASSERT (leaves);
  ASSERT (const_static_cast<DTNode*> (from. root) -> leaves == from. name2leaf. size ());
  
  
  const VectorPtr<DiGraph::Node> children (getChildren ());
 	for (const DiGraph::Node* child_ : children)
 	{
    const DTNode* child = static_cast <const DTNode*> (child_);  
    ASSERT (child);
    
    if (child->leaves <= 2)
      continue;      

    const Steiner* st = child->asSteiner ();
    ASSERT (st);

    const double subtree_len = child->getSubtreeLength ();
    ASSERT (subtree_len >= 0.0);
    if (! subtree_len)
      continue;
      
    VectorPtr<Tree::TreeNode> leafVec;  // Subset of Leaf's of *this
    child->getLeaves (leafVec);

    Set<const Tree::TreeNode*> leaves_from;  // Subset of Leaf's of tree_from
    for (const Tree::TreeNode* leaf : leafVec)
    {
      ASSERT (leaf);
      const Tree::TreeNode* leaf_from = findPtr (from. name2leaf, leaf->getName ());
      ASSERT (leaf_from);
      leaves_from << leaf_from;
    }
    ASSERT (leafVec. size () == leaves_from. size ());

    const VectorPtr<Tree::TreeNode> roots_from (from. leaves2lcas (leaves_from));
    ASSERT (! roots_from. empty ());

    const DTNode* root_from = static_cast <const DTNode*> (roots_from. front ());
    ASSERT (root_from);
    if (   roots_from. size () == 1
        && root_from->leaves == leaves_from. size ()
       )
    {
      // ratio
      const Real len_from = root_from->getSubtreeLength ();
      QC_ASSERT (len_from > 0.0);
      const Real ratio = subtree_len / len_from;
      ASSERT (ratio >= 0.0);

      ASSERT (root_from->asSteiner ());
      root_from->asSteiner () -> copySubtree (* var_cast (st), ratio);  

      for (const Tree::TreeNode* leaf : leafVec)
        var_cast (getDistTree ()). removeLeaf (const_static_cast<Leaf*> (leaf), false);
    }
    else
      var_cast (st) -> replaceSubtree (from);
  }
}



int Steiner::arcExistence_compare (const void* a, 
                                   const void* b)
{
  ASSERT (a);
  ASSERT (b);
  const Steiner* st_a = static_cast <const Steiner*> (a);
  const Steiner* st_b = static_cast <const Steiner*> (b);
  if (st_a->arcExistence < st_b->arcExistence)
    return 1;
  if (st_a->arcExistence > st_b->arcExistence)
    return -1;
  return 0;
}



void Steiner::arcExistence_index (Steiner &st, 
                                  size_t index)
{
  st. heapIndex = index;
}                         
  


#if 0
void Steiner::setSubTreeWeight ()
{
  subTreeWeight = lcaNum;
  for (const DiGraph::Arc* arc : arcs [false])
    if (const Steiner* child = static_cast <DTNode*> (arc->node [false]) -> asSteiner ())
    {
      var_cast (child) -> setSubTreeWeight ();
      subTreeWeight += child->subTreeWeight;
    }
}



void Steiner::threadNum2subTree (size_t threadNum_arg)
{
  ASSERT (threadNum_arg != no_index);
  
  threadNum = threadNum_arg;
  for (const DiGraph::Arc* arc : arcs [false])
    if (Steiner* child = var_cast (static_cast <DTNode*> (arc->node [false]) -> asSteiner ()))
      child->threadNum2subTree (threadNum_arg);
}



void Steiner::threadNum2ancestors (size_t threadNum_arg)
{
  ASSERT (threadNum_arg != no_index);

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
  var_cast (getDistTree ()). leafNum ++;  
//ASSERT (! isNan (len));
}



void Leaf::qc () const
{ 
  if (! qc_on)
    return;
  DTNode::qc ();
    
  if (graph)
  {
    QC_ASSERT (getParent ());
    QC_ASSERT (isLeaf ()); 
    if (getParent () && discernible != static_cast <const DTNode*> (getParent ()) -> childrenDiscernible ())
    {
      cout << getName () << " " << discernible << " " << getParent () -> getName () 
               << " " << static_cast <const DTNode*> (getParent ()) -> childrenDiscernible () << endl;
      ERROR;
    }
    for (const uint dissimNum : pathDissimNums)
      QC_ASSERT (getDistTree (). dissims [dissimNum]. hasLeaf (this));
    QC_IMPLY (getDistTree (). subDepth, discernible);
  }

  QC_ASSERT (! name. empty());
  QC_ASSERT (! isLeft (name, "0x"));
//QC_IMPLY (DistTree_sp::variance_min, discernible);
  if (! isNan (len) && ! discernible && len)
  {
    cout << getName () << " " << len << endl;
    ERROR;
  }  
  QC_ASSERT (getReprLeaf (0) == this);
  
  QC_ASSERT (index < no_index);
}



void Leaf::saveContent (ostream& os) const 
{ 
  DTNode::saveContent (os);
    
  if (! isNan (normCriterion))
    os << "  " << normCriterionS << '=' << normCriterion;    
	if (! discernible)
	  os << "  " << non_discernible;
}



void Leaf::setLca ()
{
  DTNode::setLca ();

  for (const uint dissimNum : pathDissimNums)
  {
    const Leaf* other = getDissimOther (dissimNum);
    ASSERT (other != this);
    if (! other->graph)
      continue;
    if (! other->tarjanLca)
      continue;
    DTNode* cluster = static_cast <DTNode*> (var_cast (other) -> getDisjointCluster ());
    ASSERT (cluster);
    const DTNode* lca = cluster->tarjanLca;
    ASSERT (lca);
    const Steiner* st = lca->asSteiner ();
    ASSERT (st);
    const Dissim& dissim = getDistTree (). dissims [dissimNum];
    ASSERT (! dissim. lca);
    var_cast (dissim). lca = st;
  }
}



const Leaf* Leaf::getDissimOther (size_t dissimNum) const
{ 
  const Dissim& dissim = getDistTree (). dissims [dissimNum];
  ASSERT (dissim. hasLeaf (this));
  return dissim. leaf1 == this ? dissim. leaf2 : dissim. leaf1; 
}



bool Leaf::isMainIndiscernible () const
{
  if (discernible)
    return true;
    
  const Tree::TreeNode* parent = getParent ();
  ASSERT (parent);
  for (const DiGraph::Arc* arc : parent->arcs [false])
  {
    const Leaf* sibling = static_cast <const DTNode*> (arc->node [false]) -> asLeaf ();
    ASSERT (sibling);
    ASSERT (! sibling->discernible);
    if (sibling->name < name)
      return false;
  }
  
  return true;
}



void Leaf::collapse (Leaf* other)
{
  ASSERT (! getDistTree (). subDepth);
  ASSERT (this != other);
//ASSERT (! DistTree_sp::variance_min);
//ASSERT (len == 0.0);
//ASSERT (discernible);

  if (! other)
    return;

  ASSERT (graph);
  ASSERT (graph == other->graph);
  
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
        indiscernibles << var_cast (child);
      }
      ASSERT (indiscernibles. size () >= 2);
    }
  }
  ASSERT (indiscernibles. contains (this));
    
  if (other->discernible)
  {
    auto st = new Steiner ( var_cast (getDistTree ())
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
    ASSERT (! other->len);
    if (getParent () != other->getParent ())
      for (Leaf* leaf : indiscernibles)
        leaf->setParent (const_static_cast <Steiner*> (other->getParent ()));
  }
  
  len = 0.0;
  discernible = false;

  const Tree::TreeNode* parent = getParent ();
  ASSERT (parent);
  ASSERT (parent == other->getParent ());
#ifndef NDEBUG
  for (const DiGraph::Arc* arc : parent -> arcs [false])
  {
    const Leaf* child = static_cast <const DTNode*> (arc->node [false]) -> asLeaf ();
    ASSERT (child);
    ASSERT (! child->discernible);
  }
#endif
}



void Leaf::addHybridTriangles (Vector<Triangle> &triangles) const
{
  ASSERT (getDistTree (). dissimParam. hybridness_min > 1.0);
  ASSERT (graph);
  
  const Vector<Dissim>& dissims = getDistTree (). dissims;

  Vector<Neighbor> neighbors;  neighbors. reserve (pathDissimNums. size ());
  for (const uint dissimNum1 : pathDissimNums)
  {
    const Dissim& dissim1 = dissims [dissimNum1];
    if (! dissim1. validMult ())
      continue;
    const Leaf* parent1 = dissim1. getOtherLeaf (this);
    ASSERT (parent1->graph);
    ASSERT (parent1 != this);
    neighbors << Neighbor (parent1, dissim1. type, dissim1. target);
  }
  neighbors. sort ();  // --> unordered_set ??
  ASSERT (neighbors. isUniq ());

  // Time: ~ O(log^4(n))  
  for (const Neighbor& neighbor1 : neighbors)
  {
    const Leaf* parent1 = neighbor1. leaf;
    for (const uint dissimNum2 : parent1->pathDissimNums)
    {
      const Dissim& dissim2 = dissims [dissimNum2];
      if (! dissim2. validMult ())
        continue;
      if (dissim2. type != neighbor1. dissimType)
        continue;
      const Leaf* parent2 = dissim2. getOtherLeaf (parent1);
      ASSERT (parent2->graph);
      ASSERT (parent2 != parent1);
      if (parent2 == this)
        continue;
      const size_t i = neighbors. binSearch (Neighbor (parent2, dissim2. type));
      if (i == no_index)
        continue;
      addTriangle ( triangles
                  , this
                  , parent1
                  , parent2
                  , dissim2. target
                  , neighbors [i]. target
                  , neighbor1. target
                  , neighbor1. dissimType
                  );
    }
  }
}




// SubPath

void SubPath::qc () const
{
  if (! qc_on)
    return;

  QC_ASSERT (dissimNum != dissims_max);
  
  QC_ASSERT (node1);
  QC_ASSERT (node2);
  QC_ASSERT (node1 != node2);
  QC_ASSERT (node1->graph);
  QC_ASSERT (node1->graph == node2->graph);
  
  QC_ASSERT (dist_hat_tails >= 0.0);
}




void SubPath::saveText (ostream &os) const
{
  os         << dissimNum
     << '\t' << node1
     << '\t' << node2
     << '\t' << dist_hat_tails
     << '\n';
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
  QC_ASSERT (area. ascending == etrue);
  QC_ASSERT (area. isUniq ());
  for (const Tree::TreeNode* node : area)
  {
    QC_ASSERT (node);
    QC_ASSERT (node->graph == & tree);
    QC_ASSERT (! static_cast <const DTNode*> (node) -> inDiscernible ());
  }
    
  // boundary
  QC_ASSERT (boundary. size () >= 2);
  QC_ASSERT (boundary. ascending == etrue);
  QC_ASSERT (boundary. isUniq ());
  QC_ASSERT (area. containsFastAll (boundary));
  QC_ASSERT ((area. size () == 2) == (area. size () == boundary. size ()));
  for (const Tree::TreeNode* node : boundary)
    QC_ASSERT ((area. size () == 2) == boundary. containsFast (node->getParent ()));

  QC_ASSERT ((bool) area_root == (bool) area_underRoot);

  if (area_root)
  {
    // area_root
  //QC_ASSERT (area_root->asSteiner ());
    QC_ASSERT (area_root->graph == & tree);
    QC_ASSERT (area. containsFast (area_root));
    QC_ASSERT (! area. containsFast (area_root->getParent ()));
    QC_ASSERT (boundary. containsFast (area_root));  
    // area_underRoot
    QC_ASSERT (area_underRoot->getParent () == area_root);
    QC_ASSERT (area. containsFast (area_underRoot));
    QC_ASSERT ((area. size () == 2) == boundary. containsFast (area_underRoot))
  }

  // subPaths
  QC_ASSERT (dissimNums. empty ());
  QC_ASSERT (! subPaths. empty ());
  {
    Vector<uint> dissimNums_;  dissimNums_. reserve (subPaths. size ());
    for (const SubPath& subPath : subPaths)
    {
      subPath. qc ();
      QC_ASSERT (boundary. containsFast (subPath. node1));
      QC_ASSERT (boundary. containsFast (subPath. node2));
      QC_ASSERT (tree. dissims [subPath. dissimNum]. valid ());
      dissimNums_ << subPath. dissimNum;
    }
    dissimNums_. sort ();
    QC_ASSERT (dissimNums_. isUniq ());
  }
}  



void Subgraph::reserve (uint radius)
{
  ASSERT (radius);
  
  const size_t boundary_size = Tree::radius2boundarySize (radius);  // PAR
  area.     reserve (boundary_size * 2); 
  boundary. reserve (boundary_size);
  
  const Real path_size_max = 2.0 * log2 ((Real) tree. name2leaf. size ());
  const Prob area_prob = min (1.0, (Real) boundary_size / (Real) tree. name2leaf. size ());
  const Prob p = 1.0 - pow (1.0 - area_prob, path_size_max);
  ASSERT (isProb (p));
  dissimNums. reserve ((size_t) ((Real) tree. dissims. size () * p)); 
}



void Subgraph::removeIndiscernibles ()
{
  ASSERT (! area. empty ());
  
//if (DistTree_sp::variance_min)
  //return;
  
  for (Iter<VectorPtr<Tree::TreeNode>> iter (area); iter. next (); )
    if (static_cast <const DTNode*> (*iter) -> inDiscernible ())
      iter. erase ();

  VectorPtr<Tree::TreeNode> newBoundary;  newBoundary. reserve (boundary. size ());
  for (Iter<VectorPtr<Tree::TreeNode>> iter (boundary); iter. next ();)
  {
    const DTNode* node = static_cast <const DTNode*> (*iter);
    if (node->inDiscernible ())
    {
      iter. erase ();
      newBoundary << node->getParent ();
    }
  }
  newBoundary. sort ();
  newBoundary. uniq ();
  if (qc_on)
    for (const Tree::TreeNode* node : boundary)
      { QC_ASSERT (! newBoundary. contains (node)); }
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



void Subgraph::dissimNums2subPaths ()
{
  ASSERT (! subPathsAbsCriterion);
  ASSERT (! dissimNums. empty ());
  ASSERT (subPaths. empty ());
  ASSERT (boundary. size () >= 2);
  

  unordered_map<uint/*dissimNum*/,size_t/*SubPath index*/> dissimNum2subPath;
  dissimNum2subPath. rehash (dissimNums. size ());
  for (const uint dissimNum : dissimNums)
    dissimNum2subPath [dissimNum] = no_index;

  dissimNums. wipe ();

  subPaths. reserve (dissimNum2subPath. size ());
  for (auto& it : dissimNum2subPath)
  {
    subPaths << SubPath (it. first);
    it. second = subPaths. size () - 1;
  }


  // SubPath::{node1,node2}
  for (const Tree::TreeNode* node : boundary)  
  {
    const DTNode* dtNode = static_cast <const DTNode*> (node);
    const Vector<uint>& pathDissimNums = boundary2pathDissimNums (dtNode);
    for (const uint dissimNum : pathDissimNums)
    {
      if (! tree. dissims [dissimNum]. valid ())  
        continue;
      size_t index = no_index;
      if (! find (dissimNum2subPath, dissimNum, index))
      {        
        ASSERT (! completeBoundary);
        continue;   
      }
      ASSERT (index != no_index); 
      SubPath& subPath = subPaths [index];
      ASSERT (subPath. dissimNum == dissimNum);
      ASSERT (! subPath. contains (dtNode));  
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
    const Dissim& dissim = tree. dissims [subPath. dissimNum];

    subPathsAbsCriterion += dissim. getAbsCriterion ();

    const VectorPtr<Tree::TreeNode>& path = getPath (subPath, buf);
    ASSERT (! path. empty ());
    const Real dist_hat_sub = DistTree::path2prediction (path);

    ASSERT (isNan (subPath. dist_hat_tails));
    subPath. dist_hat_tails = max (0.0, dissim. prediction - dist_hat_sub);

    subPath. qc ();
  }
}



#if 0
Real Subgraph::getImprovement (const DiGraph::Node2Node &boundary2new) const
{
  ASSERT (dissimNums. empty ());
  
  Real s = 0.0;
  Tree::LcaBuffer buf;
  for (const SubPath& subPath : subPaths)
  {
    const Dissim& dissim = tree. dissims [subPath. dissimNum];
    if (! dissim. validMult ())
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
#endif



namespace
{
  
void subPath2tree_dissim (Subgraph &subgraph,
                          const SubPath &subPath,
                          Tree::LcaBuffer &buf,
                          Real &absCriterion
                        #ifdef MUTEX
                         ,bool threadsUsed
                        #endif
                         )
// Update: absCriterion
{
  subPath. qc ();
  ASSERT (subPath. node1->graph == & subgraph. tree);
  ASSERT (subPath. node2->graph == & subgraph. tree);
  const uint dissimNum = subPath. dissimNum;
  const Tree::TreeNode* lca_ = nullptr;
  const VectorPtr<Tree::TreeNode>& path = Tree::getPath (subPath. node1, subPath. node2, subgraph. area_root, lca_, buf);
  ASSERT (lca_);
  Dissim& dissim = var_cast (subgraph. tree). dissims [dissimNum];
  if (! subgraph. viaRoot (subPath))
  {
    const Steiner* lca = static_cast <const DTNode*> (lca_) -> asSteiner ();
    ASSERT (lca);
    dissim. lca = lca;
  }
  for (const Tree::TreeNode* node : path)
  {
    ASSERT (node);
    const DTNode* dtNode = static_cast <const DTNode*> (node);
    ASSERT (dtNode != subgraph. area_root);
    if (qc_on)
      { QC_ASSERT (subgraph. boundary. containsFast (node) == subPath. contains (dtNode)); }
    if (! subPath. contains (dtNode))
    {
      Steiner* st = var_cast (dtNode->asSteiner ());
      ASSERT (st);
      if (qc_on && verbose ())
        { QC_ASSERT (! st->pathDissimNums. contains (dissimNum)); }
    #ifdef MUTEX
      const Lock (st->mtx, threadsUsed);
    //if (threadsUsed)
      //st->mtx. lock ();
    #endif
      st->pathDissimNums << dissimNum;
    #ifdef MUTEX
    //if (threadsUsed)
      //st->mtx. unlock ();
    #endif
    }
  }
  dissim. prediction = subPath. dist_hat_tails + DistTree::path2prediction (path);
  absCriterion += dissim. getAbsCriterion ();
}
  


#ifdef MUTEX
void subPath2tree_dissim_array (size_t from,
                                size_t to,
                                Real &absCriterion,
                                Subgraph &subgraph)
// Output: absCriterion
{
  absCriterion = 0.0;
  Tree::LcaBuffer buf;
  FOR_START (size_t, i, from, to)
    subPath2tree_dissim (subgraph, subgraph. subPaths [i], buf, absCriterion, true);
}
#endif



#if 0
void subPath2tree_subPathDissimsSets_array (size_t from,
                                            size_t to,
                                            unordered_set<uint>* &subPathDissimsSet,
                                            const Vector<SubPath> &subPaths)
{
  ASSERT (! subPathDissimsSet);
  subPathDissimsSet = new unordered_set<uint> ();
  subPathDissimsSet->rehash (to - from);
  FOR_START (size_t, i, from, to)
    subPathDissimsSet->insert (subPaths [i]. dissimNum);  
}



void subPath2tree_pathDissimNums_array (size_t from,
                                        size_t to,
                                        Notype /*&res*/,
                                        const VectorPtr<Tree::TreeNode> &area,
                                        const VectorPtr<Tree::TreeNode> &boundary,
                                        const Vector<unordered_set<uint>*> &subPathDissimsSets)
{
  FOR_START (size_t, i, from, to)
  {
    const Tree::TreeNode* node = area [i];
    if (! boundary. containsFast (node))
      const_static_cast <DTNode*> (node) -> pathDissimNums. filterValue 
        ([&subPathDissimsSets] (uint dissimNum) 
          { for (const unordered_set<uint>* subPathDissimsSet : subPathDissimsSets)
              if (subPathDissimsSet->find (dissimNum) != subPathDissimsSet->end ())
                return true;
            return false;
          }
        );
  }
}
#endif



void subPath2tree_subPathDissimsVec_array (size_t from,
                                           size_t to,
                                           Notype /*&res*/,
                                           const VectorPtr<Tree::TreeNode> &area,
                                           const VectorPtr<Tree::TreeNode> &boundary,
                                           const Vector<bool> &subPathDissimsVec)
{
  // Time: O(|area| (log |boundary| + p log(n) / n))
  FOR_START (size_t, i, from, to)
  {
    const Tree::TreeNode* node = area [i];
    if (! boundary. containsFast (node))
      const_static_cast <DTNode*> (node) -> pathDissimNums. filterValue 
        ([&subPathDissimsVec] (uint dissimNum) { return subPathDissimsVec [dissimNum]; });
  }
}



void subPath2tree_subPathDissimsSet_array (size_t from,
                                           size_t to,
                                           Notype /*&res*/,
                                           const VectorPtr<Tree::TreeNode> &area,
                                           const VectorPtr<Tree::TreeNode> &boundary,
                                           const unordered_set<uint> &subPathDissimsSet)
{
  FOR_START (size_t, i, from, to)
  {
    const Tree::TreeNode* node = area [i];
    if (! boundary. containsFast (node))
      const_static_cast <DTNode*> (node) -> pathDissimNums. filterValue 
        ([&subPathDissimsSet] (uint dissimNum) { return subPathDissimsSet. find (dissimNum) != subPathDissimsSet. end (); });
  }
}

}



void Subgraph::subPaths2tree ()
{
  ASSERT (dissimNums. empty ());

  chron_subgraph2tree. start ();  


  DistTree& tree_ = var_cast (tree);

//const size_t dissims_big = 5 * 1024 * 1024;  // PAR  
#ifdef MUTEX
  const bool useThreads = (threads_max > 1 && Threads::empty () && subPaths. size () >= dissims_big);  // PAR 
#endif
  
  // Delete subPaths from tree
  if (true /*tree. dissims. size () < dissims_big*/)  // PAR  
  {
  #ifdef MUTEX
    ASSERT (! useThreads);
  #endif
    // Time: const * p/8
    Vector<bool> subPathDissimsVec (tree. dissims. size (), false);
    // Time: O(|subPaths|)
    for (const SubPath& subPath : subPaths)
      subPathDissimsVec [subPath. dissimNum] = true;
    // Time: O(|area| (log(|boundary|) + p/n log(n)))  
    vector<Notype> notypes;
    arrayThreads (true, subPath2tree_subPathDissimsVec_array, area. size (), notypes, cref (area), cref (boundary), cref (subPathDissimsVec));
  }
  else
  {
  #if 0
    if (useThreads)
    {
      Vector<unordered_set<uint>*> subPathDissimsSets;  subPathDissimsSets. reserve (threads_max);
      arrayThreads (true, subPath2tree_subPathDissimsSets_array, subPaths. size (), subPathDissimsSets, cref (subPaths));
      vector<Notype> notypes;
      arrayThreads (true, subPath2tree_pathDissimNums_array, area. size (), notypes, cref (area), cref (boundary), cref (subPathDissimsSets));
      for (unordered_set<uint>* s : subPathDissimsSets)
        delete s;
    }
    else
  #endif
    {
      // bad_alloc
      unordered_set<uint> subPathDissimsSet;  
      // Time: O(|subPaths|)
      subPathDissimsSet. rehash (subPaths. size ());
      for (const SubPath& subPath : subPaths)
        subPathDissimsSet. insert (subPath. dissimNum);  
      // Time: O(|area| (log(|boundary|) + p/n log(n)))  
      vector<Notype> notypes;
      arrayThreads (true, subPath2tree_subPathDissimsSet_array, area. size (), notypes, cref (area), cref (boundary), cref (subPathDissimsSet));
    }
  }
  
  // Add subPaths in area to tree
  tree_. absCriterion -= subPathsAbsCriterion;
  ASSERT (tree. absCriterion < inf);
  // Time: O(|subPaths| log(|area|))
#ifdef MUTEX
  if (useThreads)  // slow 
  {
    vector<Real> absCriteria;  absCriteria. reserve (threads_max);
    arrayThreads (true, subPath2tree_dissim_array, subPaths. size (), absCriteria, ref (*this));
    for (const Real& absCriterion : absCriteria)
      tree_. absCriterion += absCriterion;
  }
  else
#endif
  {
    Tree::LcaBuffer buf;
    for (const SubPath& subPath : subPaths)
      subPath2tree_dissim (*this, subPath, buf, tree_. absCriterion /*, false*/);
  }
  ASSERT (tree. absCriterion < inf);
  maximize (tree_. absCriterion, 0.0);
  
  
  chron_subgraph2tree. stop ();
  
  tree_. qcPaths (); 
}



void Subgraph::node2dissimNums (const DTNode* node)
{ 
  ASSERT (node);
  ASSERT (node->graph == & tree);
  for (const uint dissimNum : node->pathDissimNums)
    if (tree. dissims [dissimNum]. valid ())
      dissimNums << dissimNum;
}



void Subgraph::area2dissimNums ()
// Approximation: parameter sparse: Use O(boundary.size()) dissimilarities, use getReprLeaf()'s like in sparsing ??
{
  ASSERT (dissimNums. empty ());
  ASSERT (! completeBoundary);
  
  bool first = true;
  for (const Tree::TreeNode* node : boundary)
  {
    const DTNode* dtNode = static_cast <const DTNode*> (node);
    if (dtNode == area_root)
      dtNode = area_underRoot;
    if (first)
      first = false;  // One of dtNode's is redundant given the others
    else
      node2dissimNums (dtNode);
  }
  
  completeBoundary = true;
}




// Change

void Change::qc () const
{
  if (! qc_on)
    return;
  Root::qc ();
    
  QC_ASSERT (valid ());
  QC_IMPLY (! isNan (improvement), improvement > 0.0);
  QC_ASSERT (! targets. empty ());
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
  improvement = 0.0;

  
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
    {
      Tree::LcaBuffer buf;
      area = std::move (Tree::getPath (from, to, tried ? nullptr : from->getAncestor (areaRadius_std), lca_, buf)); 
    }
    ASSERT (area. size () >= 1);
    ASSERT (lca_);
    const Steiner* lca = static_cast <const DTNode*> (lca_) -> asSteiner ();
    ASSERT (lca);
    if (   area. front () == from
        || area. front () == to
       )
      area. eraseAt (0);
    if (! area. empty ())
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
    ASSERT (area. size () >= 2);
    subgraph. finish ();
    ASSERT (subgraph. boundary. containsFast (from));
    ASSERT (subgraph. area. containsFast (to));
    IMPLY (arcEnd, subgraph. area. containsFast (arcEnd));

    subgraph. node2dissimNums (from);
    subgraph. node2dissimNums (to);
  //cout << "Change::apply: " << from << " -> " << to << endl;  
    subgraph. dissimNums2subPaths ();

    subgraph. qc ();
    if (qc_on)
      for (const SubPath& subPath : subgraph. subPaths)
        { QC_ASSERT (   subPath. contains (from)
                     || subPath. contains (to)
                     || subPath. contains (subgraph. area_root)
                    );
        }
  }
  

  // Topology
  inter = new Steiner (var_cast (tree), arcEnd, 0.0);
  var_cast (from) -> setParent (inter); 
  var_cast (to)   -> setParent (inter);
  
  // DTNode::len
  var_cast (from) -> len = 0.0;
  var_cast (to)   -> len = 0.0;
  ASSERT (! inter->len);
  

  Dataset ds;
  ds. objs. reserve (subgraph. subPaths. size ());
  auto target = new RealAttr1 ("target", ds);
  FFOR (size_t, objNum, subgraph. subPaths. size ())
    ds. appendObj ();
  auto fromAttr  = new ExtBoolAttr1 ("from",  ds);
  auto toAttr    = new ExtBoolAttr1 ("to",    ds);
  auto interAttr = new ExtBoolAttr1 ("inter", ds);
  fromAttr ->setAll (efalse);
  toAttr   ->setAll (efalse);
  interAttr->setAll (efalse);
  Tree::LcaBuffer buf;
  FFOR (size_t, objNum, subgraph. subPaths. size ())
  {
    const SubPath& subPath = subgraph. subPaths [objNum];
    const VectorPtr<Tree::TreeNode>& path = subgraph. getPath (subPath, buf);
    if (path. contains (from))
    {
      (*fromAttr) [objNum] = etrue;
      const bool toVia    = path. contains (to);
      const bool interVia = path. contains (inter);
      ASSERT (toVia != interVia);
      if (toVia)
        (*toAttr) [objNum] = etrue;
      if (interVia)
        (*interAttr) [objNum] = etrue;
    }
    else
    {
      const bool toUsed = path. contains (to);
      ASSERT (toUsed);
      ASSERT (toUsed == path. contains (inter));
      if (toUsed)
      {
        (*toAttr)    [objNum] = etrue;
        (*interAttr) [objNum] = etrue;
      }
    }
    const Dissim& dissim = tree. dissims [subPath. dissimNum];
    ASSERT (dissim. valid ());
    ASSERT (dissim. mult < inf);
    var_cast (ds. objs [objNum]) -> mult = dissim. mult; 
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
    var_cast (from) -> len = lr. beta [0];
    var_cast (to)   -> len = lr. beta [1];
    inter           -> len = lr. beta [2];
  }
  else
  {
    var_cast (from) -> len = lr. beta [0] / 2.0;
    var_cast (to)   -> len = from->len;
    inter           -> len = NaN;
    ASSERT (inter == tree. root);
  }

  Real subPathsAbsCriterion = 0.0;
  FFOR (size_t, objNum, subgraph. subPaths. size ())
  {
    const SubPath& subPath = subgraph. subPaths [objNum];
    const Dissim& dissim = tree. dissims [subPath. dissimNum];
    const Real prediction = max (0.0, dissim. target - lr. getResidual (objNum));
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
  var_cast (from) -> len = fromLen;
  var_cast (to)   -> len = toLen;

  // Topology
  var_cast (from) -> setParent (oldParent);  
  var_cast (to)   -> setParent (arcEnd);
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
    inter->pathDissimNums = to->pathDissimNums;
  else
  {
    ASSERT (to->pathDissimNums. empty ());
    var_cast (to) -> pathDissimNums = from->pathDissimNums;
  }

  subgraph. area << inter;
  subgraph. area. sort ();  
  subgraph. subPaths2tree (); 
  
  if (oldParent->isTransient ())
    var_cast (tree). delayDeleteRetainArcs (oldParent);

  subgraph. clear ();
}



bool Change::strictlyBetter (const Change* a, 
                             const Change* b)
{ 
  ASSERT (a);
  ASSERT (! isNan (a->improvement));
  
  if (a->improvement <= 0.0) 
    return false;
  
  if (a == b)  
    return false;
  if (! b)  
    return true;
  ASSERT (b->improvement >= 0.0);
    
  if (a->improvement > b->improvement)  return true;  
  if (a->improvement < b->improvement)  return false;
  if (a < b)  return true;  // non-stable result ??
    
  return false;
}



bool Change::longer (const Change* a, 
                     const Change* b)
{ 
  ASSERT (a);
  ASSERT (b);
  LESS_PART (*b, *a, arcDist);  
  return false;
}




// DissimType

DissimType::DissimType (const PositiveAttr2* dissimAttr_arg)
: Named (checkPtr (dissimAttr_arg) -> name)
, dissimAttr (dissimAttr_arg)
{
  ASSERT (dissimAttr);
}



void DissimType::qc () const
{
  if (! qc_on)
    return;
    
  Named::qc ();
  QC_ASSERT (goodName (name));
  QC_ASSERT (scaleCoeff >= 0.0);
  QC_ASSERT (scaleCoeff < inf);
}




// Dissim

Dissim::Dissim (const Leaf* leaf1_arg,
                const Leaf* leaf2_arg,
                Real target_arg,
                Real mult_arg,
                size_t type_arg)
: leaf1 (leaf1_arg)
, leaf2 (leaf2_arg)
, target (target_arg)
, type (type_arg)
, prediction (target_arg)
, mult (mult_arg)
{
  ASSERT (leaf1);
  ASSERT (leaf2);
  ASSERT (leaf1->graph);
  ASSERT (leaf1->graph == leaf2->graph);
  if (leaf1->name > leaf2->name)
    swap (leaf1, leaf2);
  ASSERT (mult >= 0.0);  
}



void Dissim::qc () const
{
  if (! qc_on)
    return;
  
  QC_ASSERT (leaf1);
  QC_ASSERT (leaf2);
  QC_ASSERT (leaf1 != leaf2);
  QC_ASSERT (leaf1->name < leaf2->name);
  
  if (isNan (mult))
    return;

  QC_ASSERT (mult >= 0.0);
        
  if (! mult)
    return;

  QC_ASSERT (valid ());
  QC_ASSERT (& leaf1->getDistTree () == & leaf2->getDistTree ());
  QC_ASSERT (target < inf);
  QC_ASSERT (prediction >= 0.0);
  QC_ASSERT (prediction < inf);
  QC_IMPLY (indiscernible (), ! prediction);
  QC_ASSERT (indiscernible () == (mult == inf));
  QC_ASSERT (lca);
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
  ASSERT (mult >= 0.0);
  if (! mult)
    return 0.0;
  const Real epsilon = prediction_arg - target;
  ASSERT (! isNan (epsilon));
  ASSERT (fabs (epsilon) < inf);
  if (! epsilon)
    return 0.0;
  return mult * sqr (epsilon);
}



void Dissim::setPathDissimNums (size_t dissimNum,
                                Tree::LcaBuffer &buf)
{
  ASSERT (valid ());

  const VectorPtr<Tree::TreeNode>& path = getPath (buf);
  prediction = DistTree::path2prediction (path);  
  for (const Tree::TreeNode* node : path)
    if (const Steiner* st = static_cast <const DTNode*> (node) -> asSteiner ())
      var_cast (st) -> pathDissimNums << (uint) dissimNum;  
}
  
  

bool Dissim::operator< (const Dissim &other) const
{ 
  LESS_PART (*this, other, leaf1);
  LESS_PART (*this, other, leaf2);
  return type < other. type;
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
    ASSERT (area. size () >= 2);
    if (areaRadius == 1 || boundary. size () <= boundary_size_max_std)
      break;
    area.     clear ();
    boundary. clear ();    
    ASSERT (areaRadius > 1);
    areaRadius--;
  }

  tree = new DistTree (subgraph, new2old, false);
  tree->qc ();
  rootInArea = (! subgraph. area_root || subgraph. area_root == wholeTree. root);

  ASSERT (area. containsFast (center));
  if (boundary. size () > boundary_size_max_std)
  //throw runtime_error (FUNC "Boundary size " + to_string (boundary. size ()) + " > " + to_string (boundary_size_max_std));
    couterr << "!!Boundary size = " << boundary. size () << endl;
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
      tree->optimizeLenWhole ();
      tree->optimizeLenArc ();
      tree->optimizeLenNode ();  
      if (subgraph. large () && areaRadius > 1)
        tree->optimizeSmallSubgraphs (areaRadius / 2);  // PAR
      else
        tree->optimizeWholeIter (20, noString);  // PAR
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
                          const VectorPtr<Tree::TreeNode> &possibleBoundary,
                          const VectorOwn<Change>* changes)
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
    ASSERT (area. size () >= 2);
    tree = new DistTree (subgraph, new2old, true);
    tree->qc ();
    rootInArea = (! subgraph. area_root || subgraph. area_root == wholeTree. root);

    // Optimization
    if (changes)
    {
      VectorOwn<Change> newChanges;  newChanges. reserve (changes->size ());
      const DiGraph::Node2Node old2new (DiGraph::reverse (new2old));
      Tree::LcaBuffer buf;
      for (const Change* change : *changes)
        if (const DTNode* from = getOld2new (change->from, old2new, buf))
        {
          ASSERT (from->graph == tree);
          const DTNode* to = getOld2new (change->to, old2new, buf);
          if (! to)
          {
            const Tree::TreeNode* lca_ = nullptr;
            const VectorPtr<Tree::TreeNode>& path = Tree::getPath (change->from, change->to, nullptr, lca_, buf);
            ASSERT (lca_);
            for (const Tree::TreeNode* node : path)
              if (node != change->from && boundary. containsFast (node))
              {
                EXEC_ASSERT (to = static_cast <const DistTree_sp::DTNode*> (findPtr (old2new, node)));
                break;
              }
          }
          const Tree::TreeNode* fromParent = from->getParent ();
          if (to)
          {
            const Tree::TreeNode* toParent = to->getParent ();
            if ( (   fromParent == to
                  || fromParent == toParent 
                 ) 
                && fromParent->arcs [false]. size () <= 2
               )  
            to = nullptr;
          }
          if (to)
          {
            auto newChange = new Change (from, to);
            ASSERT (newChange->valid ());
            newChange->improvement = change->improvement;
            newChange->arcDist     = change->arcDist;  // May be not correct due to *to change
            newChanges << newChange;
          }
        }
    //PRINT (newChanges. size ());  
      tree->applyChanges (newChanges, true); 
    }
    else
      tree->optimizeLargeSubgraphs (nullptr);
    
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
  
  
  DistTree& wholeTree = var_cast (subgraph. tree);
  const bool unstableCutP = ! wholeTree. unstableCut. empty ();
  

  DiGraph::Node2Node boundary2new (DiGraph::reverse (new2old));
  ASSERT (boundary2new. size () == new2old. size ());

  // Optimization of tree is not needed any more
  for (DiGraph::Node* node_new : tree->nodes)
    static_cast <DTNode*> (node_new) -> pathDissimNums. clear ();  

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
      section ("Re-rooting subgraph", true);
    // Replace leaf by st
    auto st = new Steiner (*tree, const_static_cast <Steiner*> (leaf->getParent ()), leaf->len);
    EXEC_ASSERT (new2old. erase (leaf) == 1);
    tree->delayDeleteRetainArcs (var_cast (leaf));
    ASSERT (isNan (static_cast <const DTNode*> (tree->root) -> len));
    const Steiner* root_ = st->makeDTRoot ();  // Old root
    boundary2new [subgraph. area_root] = st;
    new2old [st] = subgraph. area_root;
    ASSERT (root_ != st);
    ASSERT (root_->len >= 0);
    if (root_->isTransient ())
      tree->delayDeleteRetainArcs (var_cast (root_));
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
      if (unstableCutP)
        if (! dtNode->stable)
          if (const Steiner* st = dtNode->asSteiner ())
          {
            wholeTree. unstableCut. erase (st);
            wholeTree. unstableProcessed++;
          }
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


  if (center && contains (boundary2new, center))
  {
    ASSERT (center->graph == & wholeTree);
    if (unstableCutP)
      if (! center->stable)
        if (const Steiner* st = center->asSteiner ())
        {
          wholeTree. unstableCut. erase (st);
          wholeTree. unstableProcessed++;
        }
    var_cast (center) -> stable = true;
  }
  
  // wholeTree.unstableCut
  if (unstableCutP)
    for (const Tree::TreeNode* node_ : subgraph. boundary)
    {
      const DTNode* node = static_cast <const DTNode*> (node_);
      ASSERT (node->graph == & wholeTree);
      if (! node->stable)
      {
        ASSERT (node->getParent ());
        ASSERT (static_cast <const DTNode*> (node->getParent ()) -> stable);
        if (const Steiner* st = node->asSteiner ())
          wholeTree. unstableCut. insert (st);
      }
    }


  if (verbose (1))
    wholeTree. reportErrors (cout);

    
  return true;
}



const DTNode* Image::getOld2new (const DTNode* old,
                                 const DiGraph::Node2Node &old2new,
                                 Tree::LcaBuffer &buf) const
{
  ASSERT (old);
  ASSERT (tree);
  ASSERT (old->graph == & subgraph. tree);
  
  if (! subgraph. area. containsFast (old))
    return nullptr;
  
  const DTNode* node1 = old;
  while (! subgraph. boundary. containsFast (node1))
  {
    ASSERT (node1->arcs [false]. size () >= 2);
    node1 =  static_cast <const DTNode*> (node1->arcs [false]. front () -> node [false]);
  }
   
  const DTNode* node2 = old;
  while (! subgraph. boundary. containsFast (node2))
  {
    ASSERT (node2->arcs [false]. size () >= 2);
    node2 =  static_cast <const DTNode*> (node2->arcs [false]. back () -> node [false]);
  }

  const Tree::TreeNode* lca = Tree::getLca ( static_cast <const Tree::TreeNode*> (findPtr (old2new, node1))
                                           , static_cast <const Tree::TreeNode*> (findPtr (old2new, node2))
                                           , buf
                                           );
  ASSERT (lca);
   
  return static_cast <const DTNode*> (lca);
}




// DistTree

DistTree::DistTree (const DissimParam &dissimParam_arg,
	                  const string &treeDirFName,
                    const string &dissimFName,
                    const string &dissimAttrName,
	                  const string &multAttrName)
: dissimParam (dissimParam_arg)
, rand (seed_global)
{
  // Initial tree topology
  if (isDirName (treeDirFName))
    loadTreeDir (treeDirFName);  // Directory name
      // no DTNode::len
  else  
    loadTreeFile (treeDirFName);  
  ASSERT (root);
  ASSERT (nodes. front () == root);
  ASSERT (static_cast <const DTNode*> (root) -> asSteiner ());

  setName2leaf ();
        
  if (! isDirName (treeDirFName) && dissimFName. empty ())
    return;

  loadDissimDs (dissimFName, dissimAttrName, multAttrName);
  
  if (! getConnected ())
    throw runtime_error (FUNC "Disconnected objects");
  if (! setDiscernibles_ds ())
    return;

  if (isDirName (treeDirFName))
    setGlobalLen ();  // --> after dissim2Ds(), use dissims ??

  dissimDs2dissims ();      
}



DistTree::DistTree (const DissimParam &dissimParam_arg,
	                  const string &dissimFName,
                    const string &dissimAttrName,
	                  const string &multAttrName)
: dissimParam (dissimParam_arg)
, rand (seed_global)
{
  loadDissimDs (dissimFName, dissimAttrName, multAttrName);

  // Initial tree topology: star topology 
  ASSERT (! root);
  auto root_ = new Steiner (*this, nullptr, NaN);
  ASSERT (dissimDs. get ());
  for (const Obj* obj : dissimDs->objs)
    new Leaf (*this, root_, 0, obj->name);
  ASSERT (root);
  ASSERT (nodes. front () == root);
  ASSERT (static_cast <const DTNode*> (root) -> asSteiner ());

  setName2leaf ();
        
  if (! getConnected ())
    throw runtime_error (FUNC "Disconnected objects");
  if (! setDiscernibles_ds ())
    return;

  neighborJoin ();
  
  dissimDs2dissims ();
}



namespace
{

void processDissimLine (size_t from,
                        size_t to,
                        Notype /*&res*/,
                        Vector<DissimLine> &dls,
                        const DistTree::Name2leaf &name2leaf)
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
  const DTNode* center {nullptr};
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



DistTree::DistTree (const DissimParam &dissimParam_arg,
	                  const string &dataDirName,
	                  const string &treeFName,
                    bool loadNewLeaves,
                    bool loadDissim,
                    bool optimizeP)
: dissimParam (dissimParam_arg)
, rand (seed_global)
{
  ASSERT (! dataDirName. empty ());
  ASSERT (dataDirName. back () == '/');
  IMPLY (optimizeP, loadDissim);
  
  
  // Initial tree topology
  loadTreeFile (nvl (treeFName, dataDirName + "tree"));  
  ASSERT (root);
  ASSERT (nodes. front () == root);
  ASSERT (static_cast <const DTNode*> (root) -> asSteiner ());
  
  setName2leaf ();  
  
  
  VectorPtr<Leaf> newLeaves;  newLeaves. reserve (name2leaf. size () / 10 + 1);  // PAR
  if (loadNewLeaves)
  {
    section ("Loading new leaves", true);
    LineInput f (dataDirName + "leaf", 1);  // PAR
    string leafName, anchorName;
    Real leafLen, arcLen;
    Tree::LcaBuffer buf;
    Istringstream iss;
    while (f. nextLine ())
    {
      iss. reset (f. line);
      iss >> leafName >> anchorName >> leafLen >> arcLen;
      QC_ASSERT (leafLen >= 0.0);
      QC_ASSERT (arcLen >= 0.0);
    //ASSERT (iss. eof ());  // Extra fields

      const DTNode* anchor = lcaName2node (anchorName, buf);
      QC_ASSERT (anchor);
      if (const Leaf* anchorLeaf = anchor->asLeaf ())
        if (! anchorLeaf->discernible)
          anchor = static_cast <const DTNode*> (anchor->getParent ());
      ASSERT (anchor);
      while (   anchor != root 
             && arcLen > anchor->len
            )
      {
        arcLen -= anchor->len;
        anchor = static_cast <const DTNode*> (anchor->getParent ());
      }
      maximize (arcLen, 0.0);

      Leaf* leaf = nullptr;
      if (anchor == root)  
        leaf = new Leaf (*this, const_static_cast <Steiner*> (root), arcLen + leafLen, leafName);
      else 
      {
        auto st = new Steiner ( *this
                              , const_static_cast <Steiner*> (anchor->getParent ())
                              , anchor->len - arcLen
                              );
        ASSERT (st);
        var_cast (anchor) -> setParent (st);
        var_cast (anchor) -> len = arcLen;
        leaf = new Leaf (*this, st, leafLen, leafName);
        if (! arcLen && anchor->asSteiner ())  // anchor will be delete'd
        {          
          st->errorDensity = anchor->errorDensity;
          ASSERT (anchor->maxDeformationDissimNum == dissims_max);
          if (const DeformationPair* dp = findPtr (node2deformationPair, anchor))
            node2deformationPair [st] = std::move (DeformationPair {dp->leafName1, dp->leafName2, dp->deformation});
        }
      }
      ASSERT (leaf);

      name2leaf [leaf->name] = leaf;
      newLeaves << leaf;
    }
  }


  if (loadDissim)
  {
    loadDissimPrepare (name2leaf. size () * getSparseDissims_size ()); 
    const string fName (dataDirName + "dissim");
    {
      const Vector<DissimLine> dissimLines (getDissimLines (fName, dissims. capacity (), false));
      {
        Progress prog (dissimLines. size (), dissim_progress);
        for (const DissimLine &dl : dissimLines)
        {
          prog ();
          dl. apply (*this);
        }
      }
    }
    if (optimizeP)
      // qc
      for (const auto& it : name2leaf)
      {
        const Leaf* leaf = it. second;
        if (isNan (leaf->badCriterion))
          throw runtime_error (FUNC "No dissimilarities for object " + leaf->name);
        var_cast (leaf) -> badCriterion = NaN;
      }

    // qc: pairs of <leaf1,leaf2> must be unique in dissims
    if (qc_on)
    {
      Vector<LeafPair> pairs;  pairs. reserve (dissims. size ());
      for (const Dissim& d : dissims)
        pairs << LeafPair (d. leaf1, d. leaf2);
      const size_t size_init = pairs. size ();
      pairs. sort ();
      pairs. uniq ();
      if (pairs. size () != size_init)
      {
        cout << size_init << " -> " << pairs. size () << endl;
        throw runtime_error (FUNC + fName + ": non-unique pairs");
      }
    }
  } 
  
  
  cleanTopology ();
  if (getDiscernibles (). size () <= 1)
    throw runtime_error (FUNC "No discernible objects");
    
     
  if (loadDissim)
  {
    setPaths (true);  

    if (optimizeP)
    {  
      const Chronometer_OnePass cop ("Optimizing new leaves");  
      section ("Optimizing new leaves", true);
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
                th << thread (processOptimizeSmallSubgraph, ref (* var_cast (osss [i])));
              var_cast (osss [0]) -> process ();
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
  ASSERT (! multAttr);
  ASSERT (dissimTypes. empty ());
}



namespace 
{

struct Newick
/* Grammar:
     tree ::= [comment] subtree ;
     subtree ::= item [: length]
     item ::= interior | leaf
     interior ::= ( list ) [interior_name]
     list ::= subtree | subtree , list
     leaf ::= name
     name ::= string | ' string_space '
     interior_name ::= bootstrap [: name]
     string ::=  // No ' ', ',', ';', '(', ')'
     comment ::= \[ string_space \]
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
      : runtime_error (what_arg + " at position " + toString (newick. pos + 1))
      {}
  };


  Newick (istream &f_arg,
          DistTree &tree_arg)
    : f (f_arg)
    , tree (tree_arg)
    { 
      ASSERT (f. good ()); 
      ASSERT (tree. nodes. empty ());

      // Comment
      skipSpaces ();
      if (c == '[')
      {
        for (;;)
        {
          readChar ();
          if (c == ']')
            break;
        }
        readChar ();
      }

      parseSubtree (nullptr);

      skipSpaces ();
      if (c != ';')
        throw Error ("No ending semicolon", *this);
    }
    
    
private:
  void readChar ()
  {
    for (;;)
    {
      f. get (c);
      pos++;
      if (verbose ())
        cout << pos << ' ' << c << endl;
      if (   c != '\n'
          && c != '\r'
         )
        break;
    }
  }
  
  void skipSpaces ()
  {
    while (c == ' ')
      readChar ();
    if (f. eof ())
      throw Error ("Unexpected end of file", *this);
  }
  
  string parseWord (bool isName)
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
          throw Error (string ("Non-printable character in ") + (isName ? "name" : "number") + ": " + to_string ((int) c), *this);
        name += c;
        readChar ();
      }
  #if 0
    if (name. empty ())
      throw Error (string ("Empty ") + (isName ? "name" : "number"), *this);
  #else
    if (name. empty () && ! isName)
      throw Error ("Empty number", *this);
  #endif
    return unPercent (name);
  }

  void parseSubtree (Steiner* parent)
  {
    DTNode* node = parseItem (parent);
    skipSpaces ();
    Real len = 0.0;
    if (c == ':')
    {
      readChar ();
      string number (parseWord (false));
      strLower (number);
      try 
      {
        len = str2real (number);
      }
      catch (const exception &e)
      {
        throw Error (e. what (), *this);
      }
      if (isNan (len))
        len = inf;
    }
    node->len = max (0.0, len);
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
      node = new Leaf (tree, parent, NaN, parseWord (true));
    ASSERT (node);
    return node;
  }
  
  void parseInterior (Steiner* parent)
  {
    ASSERT (parent);
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
      const string s (parseWord (true));
      // bootstrap is lost ??
      const size_t colonPos = s. find (':');      
      if (colonPos != string::npos)
        parent->name = s. substr (colonPos + 1);
    }
  }
};



const string Newick::delimiters (" ,;():");
  
}  // namespace



DistTree::DistTree (DistTree::NewickFormat,
                    const string &newickFName)
: rand (seed_global)
{ 
  {
    IFStream f (newickFName);
    const Newick newick (f, *this);
  }
  ASSERT (root);
  ASSERT (nodes. front () == root);

  if (! static_cast <const DTNode*> (root) -> asSteiner ())
    throw runtime_error (FUNC "One-node tree");
  const_static_cast <DTNode*> (root) -> len = NaN;
  fixTransients ();
  finishChanges ();
  setName2leaf ();        
}



DistTree::DistTree (Prob branchProb,
                    size_t leafNum_max)
: rand (seed_global)
{
  ASSERT (isProb (branchProb));
  ASSERT (branchProb < 1.0);
  ASSERT (leafNum_max);


  const Real len = 1.0;  // PAR

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
    throw runtime_error (FUNC "One-node tree");

  setName2leaf ();
}



DistTree::DistTree (Subgraph &subgraph,
                    Node2Node &newLeaves2boundary,
                    bool sparse)
: dissimParam (subgraph. tree. dissimParam)
, subDepth (subgraph. tree. subDepth + 1)
, multFixed (true)
, rand (seed_global)
{
  ASSERT (! subgraph. empty ());  
  ASSERT (subgraph. subPaths. empty ());  // not finish()'ed  
  ASSERT (newLeaves2boundary. empty ());


  // subgraph
  subgraph. finish ();
  subgraph. area2dissimNums ();    
  subgraph. dissimNums2subPaths ();
  subgraph. qc ();

  VectorPtr<TreeNode>& area     = subgraph. area;
  VectorPtr<TreeNode>& boundary = subgraph. boundary;
  const Steiner* &area_root     = subgraph. area_root;
  const DistTree& wholeTree     = subgraph. tree;
  IMPLY (area_root, area_root->graph == & wholeTree);


  // nodes
  Node2Node old2new (area. size () + 1);  // 1-1
  ASSERT (nodes. empty ());
  ASSERT (leafNum == 0);
  if (subgraph. unresolved ())
  {
    // Star topology
    auto st = new Steiner (*this, nullptr, NaN);
    root = st;
    if (area_root)
    {
      auto leaf = new Leaf (*this, st, 0.0, "L" + toString (leafNum));
      old2new [area_root] = leaf;  
    }
    for (const TreeNode* node_old : area)
    {
      VectorPtr<DiGraph::Node> children (node_old->getChildren ());
      children. sort ();
      if (! children. intersectsFast_merge (area))
      {
        auto leaf = new Leaf (*this, st, 0.0, "L" + toString (leafNum));
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
      {
        ASSERT (! isNan (len));
        node_new = new Leaf (*this, nullptr, len, "L" + toString (leafNum));
      }
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
        ASSERT (! isNan (child->len));
        const TreeNode* root_ = root;                    // Old root in *this
        auto inter = new Steiner (*this, nullptr, NaN);  // New root in *this
        child->setParent (inter);
        child->len /= 2.0;
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
  // qc
  for (const Node* n : nodes)
  {
    ASSERT (n);
    if (n != root)
      ASSERT (! isNan (static_cast <const DTNode*> (n) -> len));
  }

  
  setName2leaf ();
  ASSERT (boundary. size () == name2leaf. size ());

  // newLeaves2boundary
  newLeaves2boundary. reserve (boundary. size ());
  for (const auto& it : old2new)
    if (static_cast <const DTNode*> (it. second) -> asLeaf ())
      newLeaves2boundary [it. second] = it. first;
  ASSERT (newLeaves2boundary. size () == boundary. size ());
    

  // dissims[]
  // For some leaf pairs the dissimilarity may be missing
  {
    unordered_map <size_t, Dissim> dissimMap;  // sparse
    Vector<uint/*dissimNum*/> leaves2dissimNum (leafNum * leafNum, dissims_max);  // !sparse  // PAR
    if (sparse)
      dissimMap. reserve (name2leaf. size () * getSparseDissims_size ());
    else
    {
      dissims. resize (getDissimSize_max ());
      ASSERT (leafNum == name2leaf. size ());
      {
        size_t dissimNum = 0;
        for (const auto& it2 : name2leaf)
          for (const auto& it1 : name2leaf)
          {
            const Leaf* leaf1 = it1. second;
            const Leaf* leaf2 = it2. second;
            if (leaf1 == leaf2)
              break;
            dissims [dissimNum] = std::move (Dissim (leaf1, leaf2, 0.0, 0.0, no_index));
            ASSERT (leaf1->index != leaf2->index);
            ASSERT (leaf1->index < no_index);
            ASSERT (leaf2->index < no_index);
            if (leaf1->index > leaf2->index)
              swap (leaf1, leaf2);
            leaves2dissimNum [leaf1->index * leafNum + leaf2->index] = (uint) dissimNum;
            dissimNum++;
            if (dissimNum >= (size_t) dissims_max)
              throw runtime_error (FUNC "dissimNum is out of uint limit");
          }
        ASSERT (dissims. size () == dissimNum);
      }
      ASSERT (dissims. size () == getDissimSize_max ());
    }


    // dissims[]: mult, target
    for (const SubPath& subPath : subgraph. subPaths)
    {
      const Dissim& wholeDissim = wholeTree. dissims [subPath. dissimNum];
      ASSERT (wholeDissim. valid ());
      
      const Real dist = wholeDissim. target - subPath. dist_hat_tails;
        // May be < 0
      if (isNan (dist))
        continue;
    
      const Real mult = wholeDissim. mult; 
      ASSERT (mult >= 0.0);
      ASSERT (mult < inf);
      if (! mult)  // wholeDissim.target = inf
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
      if (leaf1->index > leaf2->index)
        swap (leaf1, leaf2);

      const size_t index = leaf1->index * leafNum + leaf2->index;

      Dissim* dissim = nullptr;
      if (sparse)
      {
        dissim = & dissimMap [index];
        if (! dissim->leaf1)
          *dissim = std::move (Dissim (leaf1, leaf2, 0.0, 0.0, no_index));
      }
      else
      {
        const uint dissimNum = leaves2dissimNum [index];
        ASSERT (dissimNum < dissims_max);
        dissim = & dissims [dissimNum];
      }
      ASSERT (dissim);
      ASSERT (dissim->leaf1);
      dissim->target += mult * dist;
      dissim->mult   += mult;
    }
    
    
    if (sparse)
    {
      dissims. reserve (dissimMap. size ());
      for (auto& it : dissimMap) 
        dissims << std::move (it. second);
    }
  }
  
  
  // Dissim::target, mult_sum, target2_sum
  mult_sum    = 0.0;
  target2_sum = 0.0;
  for (Dissim& dissim : dissims)
    if (dissim. validMult ())
    {
      dissim. target /= dissim. mult;
      mult_sum    += dissim. mult;
      target2_sum += dissim. mult * sqr (dissim. target); 
    }
    else
      dissim. target = NaN;
    

  // DTNode::pathDissimNums[]
  const size_t reserve_size = getPathDissimNums_size ();
  for (DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <DTNode*> (node);
    ASSERT (dtNode->pathDissimNums. empty ());
    var_cast (dtNode) -> pathDissimNums. reserve (reserve_size); 
  }
  FFOR (uint, dissimNum, (uint) dissims. size ())
  {
    Dissim& dissim = dissims [dissimNum];
    var_cast (dissim. leaf1) -> pathDissimNums << dissimNum;
    var_cast (dissim. leaf2) -> pathDissimNums << dissimNum;      
  }
  setPaths (false);


  qcPaths ();
}



void DistTree::loadTreeDir (const string &dir)
{
  ASSERT (! subDepth);
  ASSERT (! dir. empty ());
  ASSERT (isDirName (dir));

#ifndef _MSC_VER
  Name2steiner name2steiner;
  DirItemGenerator dig (0, dir, false);
  string name;
  while (dig. next (name))
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
#else
  NOT_IMPLEMENTED;
#endif
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

  const StringVector lines (fName, (size_t) 10000, false);  // PAR
  QC_ASSERT (! lines. empty ());

  size_t lineNum = 0; 
  Steiner* steiner = nullptr;   
  EXEC_ASSERT (loadLines (lines, lineNum, steiner, 0));
  ASSERT (lineNum <= lines. size ());
  if (lineNum < lines. size ())
    throw runtime_error ("Only " + to_string (lineNum) + " line(s) of " + strQuote (fName) + " have been loaded");
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
      return noString;  
      
    string valueS (s1. substr (pos + token1. size ()));
    trim (valueS);
    if (valueS. empty ())
      throw runtime_error (FUNC "Empty token '" + token + "' in: " + s);
      
    if (valueS [0] == '\'')
    {
      string name;
      size_t i = 0;
      for (;;)
      {
        i++;
        if (! valueS [i])
          throw runtime_error (FUNC "No closing single quote");
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
    throw runtime_error (FUNC "Tree file is damaged");
  }
  
  lineNum++;

  const size_t colonPos = s. find (": ");
  if (colonPos == string::npos)
    throw runtime_error (FUNC "No colon in line " + toString (lineNum));
  const string idS = s. substr (0, colonPos);
  s. erase (0, colonPos + 2);
  if (idS. empty () || s. empty ())
    throw runtime_error (FUNC "Bad format of line " + toString (lineNum));
  const Real   len               = token2real (s, lenS);
  const Real   errorDensity      = token2real (s, err_densityS);
  const Real   normCriterion     = token2real (s, normCriterionS);  
        string deformationObj    = token2string (s, deformationS);
  const Real   deformation       = token2real (s, deformation_criterionS);
  const string name              = token2string (s, "name");
  const bool indiscernible = contains (s, Leaf::non_discernible);
  QC_IMPLY (parent, len >= 0.0);
  DTNode* dtNode = nullptr;
  if (isLeft (idS, "0x"))
  {
    QC_ASSERT (! indiscernible);
    Steiner* steinerParent = nullptr;
    if (parent)
    {
      steinerParent = var_cast (parent->asSteiner ());
      ASSERT (steinerParent);
    }
    auto steiner = new Steiner (*this, steinerParent, len);
    while (loadLines (lines, lineNum, steiner, expectedOffset + Offset::delta))
      ;
    dtNode = steiner;
    QC_ASSERT (isNan (normCriterion));
  }
  else
  {
    QC_ASSERT (parent);
    auto leaf = new Leaf (*this, parent, len, idS);
    leaf->discernible = /*DistTree_sp::variance_min ||*/ ! indiscernible;
    leaf->normCriterion = normCriterion;
    dtNode = leaf;
  }
  ASSERT (dtNode);
  
  dtNode->errorDensity = errorDensity;
  if (! name. empty ())
  {
    QC_ASSERT (! dtNode->asLeaf ());
    dtNode->name = name;
  }
  
  QC_ASSERT (deformationObj. empty () == isNan (deformation));
  if (! deformationObj. empty ())
  {
    const string leafName1 (findSplit (deformationObj, ':'));
    QC_ASSERT (! deformationObj. empty ());
    QC_ASSERT (deformation >= 0.0);
    ASSERT (! contains (node2deformationPair, dtNode));
    node2deformationPair [dtNode] = std::move (DeformationPair {leafName1, deformationObj, deformation});
  }
  
  return true;
}



void DistTree::setName2leaf ()
{
  name2leaf. clear ();
  name2leaf. rehash (nodes. size ());
  for (const DiGraph::Node* node : nodes)
    if (const Leaf* leaf = static_cast <const DTNode*> (node) -> asLeaf ())
      name2leaf [leaf->getName ()] = leaf;
}



namespace
{

void checkPositiveAttr2 (PositiveAttr2& attr)
{
  attr. qc ();

  {
    Real maxCorrection;
    size_t row_bad, col_bad;
    attr. matr. symmetrize (maxCorrection, row_bad, col_bad);
    if (maxCorrection > 2.0 * pow (10.0, - (Real) attr. decimals))  // PAR
      cout << "maxCorrection = " << maxCorrection 
           << " at " << attr. ds. objs [row_bad] -> name 
           << ", "   << attr. ds. objs [col_bad] -> name 
           << endl;
  }
}

}



void DistTree::loadDissimDs (const string &dissimFName,
                             const string &dissimAttrName,
                             const string &multAttrName)
{
  ASSERT (! subDepth);
  ASSERT (! dissimDs. get ());
  ASSERT (! dissimAttr);
  ASSERT (! multAttr);
  ASSERT (dissimTypes. empty ());
  ASSERT (! optimizable ());
  IMPLY (dissimAttrName. empty (), multAttrName. empty ());

  if (dissimFName. empty ())
    throw runtime_error (FUNC "Dataset " + dissimFName + " must exist");


  // dissimDs, dissimAttr, multAttr
  {
    Unverbose unv;
    dissimDs. reset (new Dataset (dissimFName));
  }

  if (dissimAttrName. empty ())
  {
    for (const Attr* attr_ : dissimDs->attrs)
      if (const PositiveAttr2* attr = attr_->asPositiveAttr2 ())
      {
        checkPositiveAttr2 (* var_cast (attr));
        dissimTypes << std::move (DissimType (attr));
      }
    if (dissimTypes. empty ())
      throw runtime_error (FUNC "No dissimilarities in " + strQuote (dissimFName));
    if (dissimTypes. size () == 1)
    {
      dissimAttr = dissimTypes [0]. dissimAttr;
      dissimTypes. clear ();
    }
    else
      mergeDissimAttrs ();
    ASSERT (! multAttr);
  }
  else
  {
    {
      const Attr* attr = dissimDs->name2attr (dissimAttrName);
      if (! attr)
        throw runtime_error (FUNC "Attribute " + strQuote (dissimAttrName) + " must exist in the dataset " + dissimFName);
      dissimAttr = attr->asPositiveAttr2 ();
    }
    if (! dissimAttr)
      throw runtime_error (FUNC "Attribute " + strQuote (dissimAttrName) + " must have type Positive2");      
    checkPositiveAttr2 (* var_cast (dissimAttr));

    if (! multAttrName. empty ())
    {
      {
        const Attr* attr = dissimDs->name2attr (multAttrName);
        if (! attr)
          throw runtime_error (FUNC "Attribute " + strQuote (multAttrName) + " must exist in the dataset " + dissimFName);
        multAttr = attr->asPositiveAttr2 ();
      }
      if (! multAttr)
        throw runtime_error (FUNC "Attribute " + strQuote (multAttrName) + " must have type Positive2");      
      checkPositiveAttr2 (* var_cast (multAttr));
      multFixed = true;
    }
  }
  
  ASSERT (dissimAttr);
}



void DistTree::mergeDissimAttrs ()
{
  ASSERT (! subDepth);
  ASSERT (dissimDs. get ());
  ASSERT (! dissimAttr);
  ASSERT (! multAttr);  
  ASSERT (dissimTypes. size () >= 2);
  ASSERT (! optimizable ());
  

  Vector<Real> vec;  vec. reserve (sqr (dissimDs->objs. size ()));
  streamsize decimals = 0;
  for (DissimType& dt : dissimTypes)
  {
    ASSERT (dt. dissimAttr);
    Real s = 0.0;
    size_t n = 0; 
    {
      Real x = NaN;   
      FFOR (size_t, i, dissimDs->objs. size ())
        FFOR (size_t, j, dissimDs->objs. size ())
          if (   dt. dissimAttr->matr. get (false, i, j, x)
              && x < inf
             )
          {
            s += x;
            n++;
          }
    }
    if (n == 0)
      throw runtime_error (FUNC "No data in dissimilarity " + strQuote (dt. dissimAttr->name));
    if (n == 1)
      throw runtime_error (FUNC "Only one value in dissimilarity " + strQuote (dt. dissimAttr->name));
    const Real average = s / (Real) n;
    ASSERT (average >= 0.0);
    if (average == 0.0)
      throw runtime_error (FUNC "Zero dissimilarity " + strQuote (dt. dissimAttr->name));
    dt. scaleCoeff = 1.0 / average;
    dt. qc ();
    ASSERT (dt. scaleCoeff > 0.0);
    maximize (decimals, dt. dissimAttr->decimals);
  }


  normalizeDissimCoeffs ();  
  
      
  dissimAttr = new PositiveAttr2 (dissimDs->findNewAttrName ("merged"), * var_cast (dissimDs. get ()), decimals + 1);
  FFOR (size_t, i, dissimDs->objs. size ())
    FFOR (size_t, j, dissimDs->objs. size ())
    {
      Real dissim_sum = 0.0;
      Real mult_sum_  = 0.0;
      Real x = NaN;
      ebool allZero = enull;
      for (const DissimType& dt : dissimTypes)
        if (   dt. dissimAttr->matr. get (false, i, j, x) 
            && x < inf
           )
        {
          ASSERT (x >= 0.0);
          if (x)
            allZero = efalse;
          else if (allZero == enull)
            allZero = etrue;            
          ASSERT (dt. scaleCoeff > 0.0);
          const Real mult = 1.0/*temporary*/ / sqr (dt. scaleCoeff);  
          dissim_sum += mult * x * dt. scaleCoeff;
          mult_sum_  += mult;
        }
      ASSERT (mult_sum_ >= 0.0);
      if (allZero == etrue)  
        var_cast (dissimAttr) -> matr. put (false, i, j, 0.0);  // To collapse()
      else if (mult_sum_)
      {
        ASSERT (dissim_sum >= 0.0);
        var_cast (dissimAttr) -> matr. put (false, i, j, dissim_sum / mult_sum_);
      }
    }
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
    if (Leaf* leaf = var_cast (dtNode->asLeaf ()))
      leaf->DisjointCluster::init ();
  }
  FFOR (size_t, row, dissimDs->objs. size ())
    if (const Leaf* leaf1 = findPtr (name2leaf, dissimDs->objs [row] -> name))
      FOR (size_t, col, row)  // dissimAttr is symmetric
        if (dissimAttr->get (row, col) < inf)
          if (const Leaf* leaf2 = findPtr (name2leaf, dissimDs->objs [col] -> name))
            var_cast (leaf1) -> merge (* var_cast (leaf2));

  Cluster2Leaves cluster2leaves;  cluster2leaves. rehash (nodes. size ());
  for (DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <DTNode*> (node);
    if (const Leaf* leaf = dtNode->asLeaf ())
      cluster2leaves [var_cast (leaf) -> getDisjointCluster ()] << leaf;
  }
  ASSERT (! cluster2leaves. empty ());
  
  if (cluster2leaves. size () == 1)
    return true;
    

  const DisjointCluster* cluster_main = nullptr;
  size_t size_max = 0;
  for (const auto& it : cluster2leaves)
    if (maximize (size_max, it. second. size ()))
      cluster_main = it. first;
  ASSERT (cluster_main);
 
  for (const auto& it : cluster2leaves)
    if (it. first != cluster_main)
    {
      cout << endl;
      cout << "Non-main cluster:" << endl;
      for (const Leaf* leaf : it. second)
        cout << leaf->name << '\n';
    }

  return false;
}



bool DistTree::getDissimConnected () 
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
        && dissim. mult
       )
      var_cast (dissim. leaf1) -> merge (* var_cast (dissim. leaf2));

  Cluster2Leaves cluster2leaves;  cluster2leaves. rehash (nodes. size ());
  for (DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <DTNode*> (node);
    if (const Leaf* leaf = dtNode->asLeaf ())
      cluster2leaves [var_cast (leaf) -> getDisjointCluster ()] << leaf;
  }
  ASSERT (! cluster2leaves. empty ());
  
  return cluster2leaves. size () == 1;
}



Cluster2Leaves DistTree::getIndiscernibles ()
{
  ASSERT (! subDepth);
  ASSERT (optimizable ());
  
//if (DistTree_sp::variance_min)
  //return Cluster2Leaves ();

  // Leaf::DisjointCluster
  for (DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <DTNode*> (node);
    if (Leaf* leaf = var_cast (dtNode->asLeaf ()))
      leaf->DisjointCluster::init ();
  }

  for (const Dissim& dissim : dissims)
    if (   dissim. valid ()
        && dissim. target <= 0.0
        && dissim. mult == inf
       )
      var_cast (dissim. leaf1) -> DisjointCluster::merge (* var_cast (dissim. leaf2));

  Cluster2Leaves cluster2leaves;  cluster2leaves. rehash (nodes. size ());
  for (DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <DTNode*> (node);
    if (const Leaf* leaf = dtNode->asLeaf ())
      cluster2leaves [var_cast (leaf) -> DisjointCluster::getDisjointCluster ()] << leaf;
  }
  
  if (cluster2leaves. size () == 1)
    throw runtime_error (FUNC "No discernible objects");
 
  for (Iter<Cluster2Leaves> iter (cluster2leaves); iter. next (); )
  {
    const auto& itMap = *iter;
    ASSERT (itMap. second. size () >= 1);
    if (itMap. second. size () == 1)
      iter. erase ();
  }

  return cluster2leaves;
}



void DistTree::leafCluster2discernibles (const Cluster2Leaves &cluster2leaves)
{
  ASSERT (! subDepth);

  size_t n = 0;   
  for (const auto& it : cluster2leaves)
  {
    const VectorPtr<Leaf>& clusterNodes = it. second;
    ASSERT (! clusterNodes. empty ());
    if (clusterNodes. size () == 1)
      continue;
      
  //ASSERT (! DistTree_sp::variance_min);
    const Leaf* first = clusterNodes [0];
    ASSERT (first);
    const DTNode* parent_ = static_cast <const DTNode*> (first->getParent ());
    ASSERT (parent_);
    Steiner* parent = var_cast (parent_->asSteiner ());
    ASSERT (parent);
    auto steiner = new Steiner (*this, parent, 0.0); 
    for (const Leaf* leaf_ : clusterNodes)
    {
      Leaf* leaf = var_cast (leaf_);
      leaf->setParent (steiner);
      leaf->discernible = false;  
      leaf->len = 0.0;  
      steiner->pathDissimNums << leaf->pathDissimNums;
      for (const uint dissimNum : leaf->pathDissimNums)
      {
        Dissim& dissim = dissims [dissimNum];
        if (dissim. indiscernible ())
          dissim. mult = inf;
      }
      n++;
    }
    ASSERT (! steiner->isTransient ());
    ASSERT (n);
    
    if (optimizable ())
    {
      VectorPtr<Leaf> leaves (clusterNodes);
      leaves. sort ();
      ASSERT (leaves. isUniq ());
      steiner->pathDissimNums. sort ();
      steiner->pathDissimNums. uniq ();
      for (const uint dissimNum : steiner->pathDissimNums)
      {
        Dissim& d = dissims [dissimNum]; 
        if (   leaves. containsFast (d. leaf1)
            && leaves. containsFast (d. leaf2)
           )
          d. lca = steiner;  
      }
      steiner->pathDissimNums. filterValue ([this, &leaves] 
                                            (uint dissimNum) 
                                            { const Dissim& d = this->dissims [dissimNum]; 
                                              return    leaves. containsFast (d. leaf1)
                                                     && leaves. containsFast (d. leaf2);
                                            }
                                           );
    }
  }
  
  if (n)
    cleanTopology ();
}



bool DistTree::setDiscernibles_ds ()
{
  ASSERT (! subDepth);
  ASSERT (dissimDs. get ());
  ASSERT (dissimAttr);
  ASSERT (! optimizable ());


  // Leaf::DisjointCluster
  for (DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <DTNode*> (node);
    if (Leaf* leaf = var_cast (dtNode->asLeaf ()))
    {
      leaf->discernible = true;
      leaf->DisjointCluster::init ();
    }
  }
  
  FFOR (size_t, row, dissimDs->objs. size ())
    if (const Leaf* leaf1 = findPtr (name2leaf, dissimDs->objs [row] -> name))
      FOR (size_t, col, row)  // dissimAttr is symmetric
        if (dissimAttr->get (row, col) <= 0.0)  // => !isNan()
          if (const Leaf* leaf2 = findPtr (name2leaf, dissimDs->objs [col] -> name))
            var_cast (leaf1) -> merge (* var_cast (leaf2));

  Cluster2Leaves cluster2leaves;  cluster2leaves. rehash (nodes. size ());
  for (DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <DTNode*> (node);
    if (const Leaf* leaf = dtNode->asLeaf ())
      cluster2leaves [var_cast (leaf) -> getDisjointCluster ()] << leaf;
  }
  ASSERT (! cluster2leaves. empty ());
  
  leafCluster2discernibles (cluster2leaves);
  
  return cluster2leaves. size () > 1;
}



void DistTree::setDiscernibles ()
{ 
  ASSERT (! subDepth);
  ASSERT (optimizable ());

  for (DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <DTNode*> (node);
    if (Leaf* leaf = var_cast (dtNode->asLeaf ()))
      leaf->discernible = true;
  }  
  
//if (DistTree_sp::variance_min)
  //return 0;
  
  const Cluster2Leaves cluster2leaves (getIndiscernibles ());
  leafCluster2discernibles (cluster2leaves);

  qcPaths ();
}



size_t DistTree::fixTransients ()
{ 
  ASSERT (! subDepth);

  size_t n = 0;
  for (;;)
  {
    bool found = false;
    Vector<DiGraph::Node*> nodeVec;  nodeVec. reserve (nodes. size ());
    insertAll (nodeVec, nodes);
    for (DiGraph::Node* node_ : nodeVec)
    {
      DTNode* node = static_cast <DTNode*> (node_);
      if (   node->isTransient ()
          || (node->isLeaf () && node->asSteiner ())
         )
      {
        delayDeleteRetainArcs (node);
        found = true;
        n++;
      }
    }
    if (! found)
      break;
  }
  toDelete. deleteData ();
  
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
        delayDeleteRetainArcs (var_cast (s));
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
      delayDeleteRetainArcs (var_cast (s));
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
    if (Leaf* leaf = var_cast (dtNode->asLeaf ()))
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
          Steiner* s = var_cast (static_cast <const DTNode*> (ancestor) -> asSteiner ());
          ASSERT (s);
          s->subtreeLen. add (max (0.0, d) / 2);  
        }
      }

  // DTNode::len
  for (DiGraph::Node* node : nodes)
  {
    DTNode* dtNode = static_cast <DTNode*> (node);
    if (dtNode->inDiscernible ())
      { ASSERT (! dtNode->len); }
    else
      if (const DTNode* parent = static_cast <const DTNode*> (dtNode->getParent ()))
        dtNode->len = max (0.0, parent->subtreeLen. getMean () - dtNode->subtreeLen. getMean ());
  }

  clearSubtreeLen ();

  finishChanges ();
}



namespace
{
  struct NodePair
  {
    // !nullptr, discernible
    array <const DTNode*, 2> nodes;
      // nodes[0] < nodes[1]
    Real dissim {NaN};
      // !isNan()

    NodePair ()
      { nodes [0] = nullptr;
        nodes [1] = nullptr;
      }
    NodePair (const Leaf* leaf1,
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
           << '\n'; 
      }
    bool operator== (const NodePair& other) const
      { return    nodes [0] == other. nodes [0]
               && nodes [1] == other. nodes [1];
      }
    void orderNodes ()
      { swapGreater (nodes [0], nodes [1]); }
    bool same () const
      { return nodes [0] == nodes [1]; }
    bool merge (const NodePair &from)
      { if (! (*this == from))
          return false;
        dissim = (dissim + from. dissim) / 2.0;
        return true;
      }
    Real getCriterion (size_t n) const
      { ASSERT (n > 2);
        return dissim - (nodes [0] -> len + nodes [1] -> len) / (Real) (n - 2); 
      }
      // James A. Studier, Karl J. Keppler, A Note on the Neighbor-Joining Algorithm of Saitou and Nei
    Real getParentDissim (size_t n) const
      { ASSERT (n > 2);
        return min (dissim, max (0.0, 0.5 * (dissim + (nodes [0] -> len - nodes [1] -> len) / (Real) (n - 2)))); 
      }
    static bool strictlyLess (const NodePair& n1,
                              const NodePair& n2)
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
  for (DiGraph::Arc* arc : root->arcs [false])
  {
    static_cast <DTNode*> (arc->node [false]) -> len = 0.0;
    n++;
  }
    
  size_t missing = 0;
  Vector<NodePair> leafPairs;  leafPairs. reserve (getOneDissimSize_max ());  
  if (optimizable ())
  {
    ASSERT (dissims. size () == getOneDissimSize_max ());
    ASSERT (dissimTypes. empty ());
    for (const Dissim& dissim : dissims)
    {
      if (isNan (dissim. target))
      {
        missing++;
        continue;
      }
      const NodePair leafPair (dissim. leaf1, dissim. leaf2, dissim. target);  
      if (leafPair. same ())
        continue;
      if (leafPair. dissim == inf)
      {
        missing++;
        continue;  
      }
      leafPairs << leafPair;
    }
  }
  else
  {
    section ("Neighbor joining", true);  
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
        NodePair leafPair (leaf1, leaf2, dissimAttr->get (row, col));
        if (leafPair. same ())
          continue;
        if (leafPair. dissim < 0.0)
          throw runtime_error (FUNC "Negative distance for " + leaf1->name + " - " + leaf2->name);
        if (leafPair. dissim == inf || isNan (leafPair. dissim))
        {
          missing++;
          continue;  
        }
        dissimParam. transform (leafPair. dissim);   
        leafPairs << leafPair;
      }
    }
  }
  if (leafPairs. empty ())
    return;
  

  Progress prog (n - 1);
  NodePair leafPair_best;
  for (;;)
  {
    prog ();
    
    // Remove duplicate NodePair 
    leafPairs. sort (NodePair::strictlyLess);
    leafPairs. filterIndex ([&leafPairs, &leafPair_best] (size_t i) 
                                 { return    leafPairs [i] == leafPair_best 
                                          || (i && leafPairs [i - 1]. merge (leafPairs [i]));
                                 }
                              );    

    if (leafPairs. size () == 1)
      break;
    ASSERT (n > 2);
    
    if (! leafPair_best. nodes [0])  // First iteration
      for (NodePair& leafPair : leafPairs)
        for (const bool first : {false, true})
          var_cast (leafPair. nodes [first]) -> len += leafPair. dissim;
          
    if (verbose (-2))  
    {
      cout << endl << "Nodes and sum:" << endl;
      for (const DiGraph::Arc* arc : root->arcs [false])
      {
        const DTNode* node = static_cast <const DTNode*> (arc->node [false]);
        cout << node->getLcaName () << ": " << node->len << '\n';
      }     
      cout << endl << "Pairs:" << endl;
      for (const NodePair& leafPair : leafPairs)
        leafPair. print (cout);
    }
      
      
    Steiner* newNode = nullptr;
    {
      // leafPair_best
      Real dissim_min = NaN;
      {
        Real criterion = inf;
        size_t i_best = no_index;
        FFOR (size_t, i, leafPairs. size ())
          // P (criterion1 < criterion2) ??
          if (minimize (criterion, leafPairs [i]. getCriterion (n)))            
            i_best = i;
        if (i_best == no_index)
          throw runtime_error (FUNC "Bad dissimilarity");
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
        DTNode* a = var_cast (leafPair_best. nodes [0]);
        DTNode* b = var_cast (leafPair_best. nodes [1]);
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
    for (NodePair& leafPair : leafPairs)
      if (! (leafPair == leafPair_best))
        for (const bool first : {false, true})
        {
          bool found = false;
          for (const bool best_first : {false, true})
            if (leafPair. nodes [first] == leafPair_best. nodes [best_first])
            {
              leafPair. nodes [first] = newNode;
              var_cast (leafPair. nodes [! first]) -> len -= leafPair. dissim;
              leafPair. dissim -= leafPair_best. nodes [best_first] -> len;
              maximize (leafPair. dissim, 0.0);
              var_cast (leafPair. nodes [! first]) -> len += leafPair. dissim / 2.0;  
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
  
  
  NodePair& leafPair = leafPairs [0];
  ASSERT (! leafPair. same ());
  for (const bool first : {false, true})  
  {
    DTNode* node = var_cast (leafPair. nodes [first]);
    ASSERT (node->getParent () == root);
    node->len = leafPair. dissim / 2.0;  
  }
  // If dissims is not complete
  for (const DiGraph::Arc* arc : root->arcs [false])
  {
    const DTNode* node = static_cast <const DTNode*> (arc->node [false]);
    ASSERT (node->len >= 0.0);
  }
  
  if (optimizable ())
    setPaths (true);

  finishChanges ();
  
  reroot (true);  // Reduces dissims.size() if dissims are sparse
}



void DistTree::dissimDs2dissims ()
{
  ASSERT (dissimDs. get ());
  ASSERT (dissimAttr);
  ASSERT (! optimizable ());
  ASSERT (detachedLeaves. empty ());


  if (qc_on)
  {
    const StringVector objNames (dissimDs->getObjNames ());
    if (! objNames. containsFastAll (name2leaf))
      throw runtime_error (FUNC "Tree has more objects than the dataset");
  }

  // Leaf::comment
  FFOR (size_t, objNum, dissimDs->objs. size ())
  {
    const Obj* obj = dissimDs->objs [objNum];
    const Leaf* leaf = findPtr (name2leaf, obj->name);
    if (! leaf)
      throw runtime_error (FUNC "Object " + strQuote (obj->name) + " is not in the tree");
    var_cast (leaf) -> comment = obj->comment;
  }

  if (name2leaf. size () != dissimDs->objs. size ())
    throw runtime_error (FUNC "Mismatch of the tree objects and the dataset objects");


  // dissims[]
  loadDissimPrepare (getDissimSize_max ());
  FFOR (size_t, row, dissimDs->objs. size ())
  {
    const string name1 = dissimDs->objs [row] -> name;
    const Leaf* leaf1 = findPtr (name2leaf, name1);
    if (! leaf1)
      continue;
    FOR (size_t, col, row)  // dissimAttr, multAttr are symmetric
    {
      const string name2 = dissimDs->objs [col] -> name;
      const Leaf* leaf2 = findPtr (name2leaf, name2);
      if (! leaf2)
        continue;
      if (dissimTypes. empty ())
      {
        const Real dissim = dissimAttr->get (row, col);
        const Real mult = multAttr ? multAttr->get (row, col) : 1.0/*temporary*/;
        addDissim (name1, name2, dissim, mult, no_index);
      }
      else
        FFOR (size_t, type, dissimTypes. size ())
        {
          const DissimType& dt = dissimTypes [type];
          const Real dissim = dt. dissimAttr->get (row, col);
          addDissim (name1, name2, dissim, 1.0/*temporary*/, type);
        }
    }
  }
  
  dissimDs. reset (nullptr);
  dissimAttr = nullptr;
  multAttr = nullptr;
  for (DissimType& dt : dissimTypes)
    dt. dissimAttr = nullptr;


  setPaths (true);
}



void DistTree::loadDissimPrepare (size_t pairs_max)
{
  ASSERT (dissims. empty ());
  
  if (verbose ())
    section ("Leaf pairs -> data objects", true);

  dissims. reserve (pairs_max);

  const size_t reserve_size = getPathDissimNums_size ();
  for (DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <DTNode*> (node);
    ASSERT (dtNode->pathDissimNums. empty ());
    var_cast (dtNode) -> pathDissimNums. reserve (reserve_size);  
  }
}



Vector<DissimLine> DistTree::getDissimLines (const string& fName,
                                             size_t reserveSize,
                                             bool mayBeEmpty) const
{
  Vector<DissimLine> dissimLines;  dissimLines. reserve (reserveSize);
  {
    {
      section ("Loading " + fName, true);
      LineInput f (fName, dissim_progress);  
      while (f. nextLine ())
      {
        DissimLine dl (f. line, f. lineNum);
        if (DM_sp::finite (dl. dissim))
        {
          dissimLines << std::move (dl);
          ASSERT (dl. name1. empty ());
        }
      }
      if (! mayBeEmpty && ! f. lineNum)
        throw runtime_error (FUNC "Empty " + fName);
    }
    section ("Sorting dissimilarities", true);
    dissimLines. sort ();
    dissimLines. uniq ();
  }
  vector<Notype> notypes;
  arrayThreads (true, processDissimLine, dissimLines. size (), notypes, ref (dissimLines), cref (name2leaf));
  
  return dissimLines;
}
  


bool DistTree::addDissim (Leaf* leaf1,
                          Leaf* leaf2,
                          Real target,
                          Real mult,
                          size_t type)
{
  ASSERT (dissimParam. power > 0.0);
  ASSERT (dissimParam. coeff > 0.0);

  ASSERT (detachedLeaves. empty ());

  ASSERT (leaf1);
  ASSERT (leaf2);
  IMPLY (! multFixed, mult == 1.0);  // mult will be set later
  IMPLY (/*! DistTree_sp::variance_min &&*/ ! target && dissimTypes. empty (), leaf1->getCollapsed (leaf2));  

  
  if (   isNan (target)  // prediction must be large ??
      || target == inf
     )
  {
    if (verbose ())
      cout << leaf1->name << " - " << leaf2->name << ": " << target << endl;
    return false;  
  }
  
  ASSERT (mult >= 0.0);
//ASSERT (mult < inf);

  dissimParam. transform (target);  
    
  if (type != no_index)
    target *= dissimTypes [type]. scaleCoeff;
  
  Dissim d (leaf1, leaf2, target, mult, type);  
  
  const size_t dissimNum_ = dissims. size ();
  if (dissimNum_ > (size_t) dissims_max)
    throw runtime_error (FUNC "Too large dissimNum");
  dissims << std::move (d);
  const uint dissimNum = (uint) dissimNum_;
  leaf1->pathDissimNums << dissimNum;
  leaf2->pathDissimNums << dissimNum;      
  
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
    section ("One thread", false);
  Tree::LcaBuffer buf;
  for (const size_t dissimNum : subTree) 
  {
    prog ();
    dissims [dissimNum]. setPathDissimNums (dissimNum, buf);
    absCriterion += dissims [dissimNum]. getAbsCriterion ();
  }
  ASSERT (absCriterion < inf);
}
  
}
#endif



void DistTree::setPaths (bool setDissimMultP)
{
//ASSERT (optimizable ());
  if (dissims. size () > (size_t) dissims_max)
    throw runtime_error (FUNC "Too large dissimNum");
    
  setLca ();  
  
  for (DiGraph::Node* node : nodes)
    if (Steiner* st = var_cast (static_cast <const DTNode*> (node) -> asSteiner ()))
      st->pathDissimNums. clear ();  

  absCriterion = 0.0;
  Tree::LcaBuffer buf;
#if 0
  if (   subDepth 
      || threads_max == 1 
      || dissims. size () < smallDissims
     )
#endif
  {
    Progress prog (dissims. size (), dissim_progress); 
    FFOR (size_t, dissimNum, dissims. size ()) 
    {
      prog ();
      dissims [dissimNum]. setPathDissimNums (dissimNum, buf);
      if (! setDissimMultP)
        absCriterion += dissims [dissimNum]. getAbsCriterion ();
    }
  ASSERT (absCriterion < inf);
  }
#if 0
  else
  {    
    for (Dissim& dissim : dissims)
    {
      ASSERT (dissim. lca);
      var_cast (dissim. lca) -> lcaNum++;
    }
  
    VectorPtr<Steiner> cut_best;  cut_best. reserve (100);  // PAR
    {
      const Steiner* root_ = static_cast <const Steiner*> (root);
      var_cast (root_) -> setSubTreeWeight ();
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
      FFOR (size_t, dissimNum, dissims. size ()) 
      {
        prog ();
        dissims [dissimNum]. setPathDissimNums (dissimNum, buf);
        absCriterion += dissims [dissimNum]. getAbsCriterion ();
      }
      ASSERT (absCriterion < inf);
    }
    else
    {
      size_t threadNum = 0;
      for (const Steiner* st : cut_best)
      {
        threadNum++;
        ASSERT (threadNum);
        var_cast (st) -> threadNum2subTree (threadNum);
        var_cast (st) -> threadNum2ancestors (0);
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
      FFOR (size_t, dissimNum, dissims. size ()) 
      {
        Dissim& dissim = dissims [dissimNum];
        const Steiner* st = dissim. lca;
        ASSERT (st);
        ASSERT (st->threadNum <= cut_best. size ());
        subTrees [st->threadNum] << dissimNum;
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

  if (setDissimMultP)
    setDissimMult (true);
}



Json* DistTree::toJson (JsonContainer* parent_arg,
                        const string& name_arg) const
{
  ASSERT (parent_arg);

  map<const DTNode*, size_t> node2index;  
  {
    size_t index = 1;
    for (const DiGraph::Node* node : nodes)
    {
      const DTNode* dtNode = static_cast <const DTNode*> (node);
      node2index [dtNode] = index;
      index++;
    }
  }

  auto arr = new JsonMap (parent_arg, name_arg);
  size_t index = 1;
  for (const DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <const DTNode*> (node);
    auto j = new JsonMap (arr, to_string (index));
    if (const DTNode* parent = static_cast <const DTNode*> (dtNode->getParent ()))
      new JsonInt ((long long) node2index [parent], j, "parent");
    else
      new JsonNull (j, "parent");
    dtNode->toJson (j, noString);
    index++;
  }
  
  return arr;
}



void DistTree::qc () const
{ 
  if (! qc_on)
    return;
  Tree::qc ();


  // Global ??
  QC_ASSERT (radius2boundarySize (areaRadius_std) <= boundary_size_max_std);
  QC_ASSERT (DistTree_sp::variance_min >= 0.0);
  QC_ASSERT (! isNan (DistTree_sp::variancePower) == (DistTree_sp::varianceType == varianceType_pow));
  QC_IMPLY (! isNan (DistTree_sp::variancePower), DistTree_sp::variancePower >= 0.0);
    
  dissimParam. qc ();


  QC_ASSERT (nodes. size () >= 2);
    
  QC_ASSERT (root);
  QC_ASSERT (root->graph == this);
  const Steiner* root_ = static_cast <const DTNode*> (root) -> asSteiner ();
  QC_ASSERT (root_);

  for (const DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <const DTNode*> (node);
    QC_ASSERT (! dtNode->isTransient ());
  }

  const size_t leavesSize = root->getLeavesSize ();
  QC_ASSERT (leavesSize);
  QC_ASSERT (name2leaf. size () == leavesSize);
  Vector<size_t> leafIndices;  leafIndices. reserve (leavesSize);
  for (const auto& it : name2leaf)
  {
    const Leaf* leaf = it. second;
    QC_ASSERT (leaf);
    QC_ASSERT (leaf->graph == this);
    leafIndices << leaf->index;
  }
  leafIndices. sort ();
  leafIndices. uniq ();
  QC_ASSERT (leafIndices. size () == leavesSize);
  
  QC_ASSERT ((bool) dissimDs. get () == (bool) dissimAttr);
  if (dissimAttr)
  {
    dissimDs->qc ();
    QC_ASSERT (& dissimAttr->ds == dissimDs. get ());
  }
  QC_IMPLY (multAttr, dissimAttr);
  QC_IMPLY (multAttr, & multAttr->ds == dissimDs. get ());
  
  if (! dissimTypes. empty ())
  {
    QC_ASSERT (dissimTypes. size () >= 2);
    QC_ASSERT (! multAttr);
    QC_ASSERT (! subDepth);
    QC_ASSERT (eqReal (getDissimCoeffProd (), 1.0, dissimCoeffProd_delta)); 
    for (const DissimType& dt : dissimTypes)
    {
      dt. qc ();
      QC_ASSERT ((bool) dt. dissimAttr == (bool) dissimAttr);
      QC_IMPLY (dt. dissimAttr, & dissimAttr->ds == & dt. dissimAttr->ds);
    }
  }

  QC_IMPLY (subDepth, multFixed);


  size_t leafDissims = 0;
  if (optimizable ())
  {
    Set<LeafPair> leafSet;
    for (const Dissim& dissim : dissims)
      if (dissim. mult)
      {
        dissim. qc ();
        QC_IMPLY (! subDepth, dissim. target >= 0.0);
        QC_IMPLY (/*! DistTree_sp::variance_min &&*/ ! subDepth && ! dissim. target, dissim. indiscernible ());
        leafSet. addUnique (LeafPair (dissim. leaf1, dissim. leaf2));
        if (subDepth)
          { QC_ASSERT (dissim. mult < inf); }
        else
          { QC_IMPLY (! dissim. target && ! DistTree_sp::variance_min, dissim. mult == inf); }
        if (dissimTypes. empty ())
          { QC_ASSERT (dissim. type == no_index); }
        else
          { QC_ASSERT (dissim. type < dissimTypes. size ()); } 
      }
    
    QC_ASSERT (absCriterion >= 0.0);
    QC_ASSERT (absCriterion < inf);
    QC_ASSERT (target2_sum >= 0.0);    
    QC_ASSERT (mult_sum >= 0.0);    
        
    for (const DiGraph::Node* node : nodes)
    {
      const DTNode* dtNode = static_cast <const DTNode*> (node);
      QC_ASSERT ((dtNode == root) == dtNode->pathDissimNums. empty ());  
      if (const Leaf* leaf = dtNode->asLeaf ())
        leafDissims += leaf->pathDissimNums. size ();
    }  
  }
  else
  {
    QC_ASSERT (! dissimDs. get ());
    QC_ASSERT (! dissimAttr);
    QC_ASSERT (! multAttr);
    QC_ASSERT (dissimTypes. empty ());
  }
  

  const size_t discernibles = getDiscernibles (). size ();
  QC_ASSERT (discernibles <= leavesSize);

  for (const Leaf* leaf : detachedLeaves)
  {
    QC_ASSERT (leaf);
    QC_ASSERT (! leaf->graph);
    leaf->qc ();
    leafDissims += leaf->pathDissimNums. size ();
  }
  const size_t allLeaves = leavesSize + detachedLeaves. size ();
  QC_ASSERT (dissims. size () <= (allLeaves * (allLeaves - 1)) / 2 * dissimTypesNum ());

  QC_ASSERT (leafDissims == 2 * dissims. size ());

/*
  // May happen due to removeLeaf()
  for (const DTNode* dtNode : toDelete)
    QC_ASSERT (! dtNode->asLeaf ());
*/

  if (! subDepth && optimizable ())
  {
    {
      const Cluster2Leaves cluster2leaves (var_cast (this) -> getIndiscernibles ());
      for (const auto& it : cluster2leaves)
      {
        const VectorPtr<Leaf>& leaves = it. second;
        QC_ASSERT (leaves. size () >= 2);
      //QC_ASSERT (! DistTree_sp::variance_min);
        for (const Leaf* leaf : leaves)
        {
          QC_ASSERT (leaf->graph == this);
          if (leaf->getParent () != leaves [0] -> getParent ())
            throw runtime_error (FUNC "Indiscernibles have different parents: " + leaves [0] -> name + " and " + leaf->name);
        }
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
        QC_ASSERT (child);
        if (dc)
        { 
          if (child->DisjointCluster::getDisjointCluster () != dc)
            throw runtime_error (FUNC "Bad indiscernible: " + child->name);
        }
        else
          dc = child->DisjointCluster::getDisjointCluster ();
      }
    }
  }
  
  for (const auto& it : node2deformationPair)
  {
    const DeformationPair& cp = it. second;
    QC_ASSERT (cp. leafName1 < cp. leafName2);
    QC_ASSERT (contains (name2leaf, cp. leafName1));
    QC_ASSERT (contains (name2leaf, cp. leafName2));
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
  const DTNode* parent_ = static_cast <const DTNode*> (leaf->getParent ());
  ASSERT (parent_);
  const Steiner* parent = parent_->asSteiner ();
  ASSERT (parent);
  delayDeleteRetainArcs (var_cast (leaf));
  if (deleteTransientAncestor && parent->isTransient ())
    delayDeleteRetainArcs (var_cast (parent));
  
  toDelete. deleteData ();
  
  node2deformationPair. clear ();
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
  if (! leaf1)
    throw runtime_error (FUNC "Object " + strQuote (name1) + " is not in the tree");
  const Leaf* leaf2 = s. empty () ? leaf1 : findPtr (name2leaf, s);
  if (! leaf2)
    throw runtime_error (FUNC "Object " + strQuote (s) + " is not in the tree");
  const DTNode* node = static_cast<const DTNode*> (getLca (leaf1, leaf2, buf));
  ASSERT (node);
  
  return node;
}



VectorPtr<DTNode> DistTree::getDiscernibles () const
{
  VectorPtr<DTNode> s;  s. reserve (name2leaf. size ());
  for (const DiGraph::Node* node : nodes)
    if (const Leaf* leaf = static_cast <const DTNode*> (node) -> asLeaf ())
      s << leaf->getDiscernible ();
  s. sort ();
  s. uniq ();
  
  return s;
}

  
  
void DistTree::setGoodLeaves (const string &goodFName)
{
  LineInput f (goodFName);
  while (f. nextLine ()) 
  {
    trim (f. line);
    if (strBlank (f. line))
      continue;
    if (const Leaf* leaf = findPtr (name2leaf, f. line))
    {
      const DTNode* discernible = leaf->getDiscernible ();
      ASSERT (discernible);
      for (const DiGraph::Arc* arc : discernible->arcs [false])
      {
        const Leaf* child = static_cast <const DTNode*> (arc->node [false]) -> asLeaf ();
        ASSERT (child);
        if (child->good)
          break;  // All children are good
        var_cast (child) -> good = true;
      }
      var_cast (leaf) -> good = true;
    }
  }
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
    os << "Dissimilarities factor: " << (Real) dissims. size () / ((Real) dissimTypesNum () * (Real) discernibles * log ((Real) discernibles)) << endl;
  }
  {
    const ONumber on (os, dissimDecimals, false);
    os << "Ave. dissimilarity = " << getDissim_ave () << endl;
  }
  if (! dissimTypes. empty ())
    os << "# Dissimilarity types = " << dissimTypes. size () << endl;
  if (multFixed)
    os << "Variance is fixed" << endl;
  reportErrors (os);
}



void DistTree::saveDissimCoeffs (const string &fName) const
{
  if (fName. empty ())
    return;  
        
  OFStream f (fName);
  for (const DissimType& dt : dissimTypes)
    f << dt. name << '\t' << dt. scaleCoeff << endl;
}



void DistTree::saveFeatureTree (const string &fName,
                                bool withTime) const
{
  if (fName. empty ())
    return;  
  OFStream f (fName);
  static_cast <const DTNode*> (root) -> saveFeatureTree (f, withTime, 0);
}



void DistTree::qcPaths () 
{
  if (! qc_on)
    return;

  QC_ASSERT (const_static_cast <const DTNode*> (root) -> pathDissimNums. empty ());

  size_t pathObjNums_checked = 0;
  size_t lcaObjNums_checked = 0;
  {
    const size_t displayNodes = 1000;  // PAR
    Progress prog (nodes. size (), displayNodes);  
    if (prog. active && nodes. size () >= displayNodes)
      section ("QC paths", false);
    for (DiGraph::Node* node : nodes)
    {
      prog ();
      DTNode* dtNode = static_cast <DTNode*> (node);
      Vector<uint>& pathDissimNums = dtNode->pathDissimNums;
      pathDissimNums. sort ();
      QC_ASSERT (pathDissimNums. isUniq ());
      for (const uint dissimNum : pathDissimNums)
        if (dissims [dissimNum]. valid ())
          pathObjNums_checked++;
      const Vector<uint> lcaObjNums (dtNode->getLcaDissimNums ());
      for (const uint dissimNum : lcaObjNums)
        if (dissims [dissimNum]. valid ())
        {
          QC_ASSERT (dissims [dissimNum]. lca == dtNode->asSteiner ());
          lcaObjNums_checked++;
        }
    }
  }
  
  size_t pathObjNums_all = 0;
  size_t lcaObjNums_all = 0;
  {
    Progress prog (dissims. size (), dissim_progress);
    LcaBuffer buf;
    FFOR (size_t, dissimNum, dissims. size ())
    {
      prog ();
      if (dissims [dissimNum]. valid ())
      {
        const VectorPtr<TreeNode>& path = dissims [dissimNum]. getPath (buf);
        for (const TreeNode* node : path)
        {
          QC_ASSERT (static_cast <const DTNode*> (node) -> pathDissimNums. containsFast ((uint) dissimNum));
          pathObjNums_all++;
        }
        lcaObjNums_all++;
      }
    }
  }

  QC_ASSERT (pathObjNums_all == pathObjNums_checked);
  QC_ASSERT (lcaObjNums_all == lcaObjNums_checked);
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
  absCriterion = 0.0;
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
    if (dissim. mult < inf)
      absCriterion += dissim. getAbsCriterion ();
  }
  ASSERT (absCriterion < inf);
}

}



void DistTree::setPredictionAbsCriterion ()
{
  if (subDepth)
    setPredictionAbsCriterion_thread (0, dissims. size (), absCriterion, dissims);
  else
  {
    absCriterion = 0.0;
    vector<Real> absCriteria;
    arrayThreads (true, setPredictionAbsCriterion_thread, dissims. size (), absCriteria, ref (dissims)); 
    for (const Real x : absCriteria)
      absCriterion += x;
    ASSERT (absCriterion < inf);
  }
}



void DistTree::qcPredictionAbsCriterion () const
{
  if (! qc_on)
    return;

  Real absCriterion_ = 0.0;
  LcaBuffer buf;
  for (const Dissim& dissim : dissims)
    if (dissim. validMult ())
    {
      const VectorPtr<TreeNode>& path = dissim. getPath (buf);
      const Real prediction_ = path2prediction (path);
      QC_ASSERT_EQ (prediction_, dissim. prediction, 1e-3);  // PAR
      absCriterion_ += dissim. getAbsCriterion ();
    }
  if (   fabs (absCriterion - absCriterion_) > 1e-3        // PAR
      && fabs (log (absCriterion / absCriterion_)) > 1e-3  // PAR
     )  
  {
    cout << absCriterion << " " << absCriterion_ << endl;
    ERROR;
  }
}



Real DistTree::path2prediction (const VectorPtr<TreeNode> &path) 
{
  Real dHat = 0.0;
  for (const TreeNode* node : path)
  {
    ASSERT (node);
    dHat += static_cast <const DTNode*> (node) -> len;
  }
  if (isNan (dHat))
  {
    for (const TreeNode* node : path)
      cout << node 
           << ' ' << static_cast <const DTNode*> (node) -> name
           << ' ' << static_cast <const DTNode*> (node) -> len 
           << endl;
    static_cast <const DTNode*> (path. front ()) -> getDistTree (). saveFile ("aa.tree");  // ??
    ERROR;
  }
  ASSERT (dHat >= 0.0);
  return dHat;
}



void DistTree::setDissimMult (bool usePrediction)
{ 
  ASSERT (optimizable ());
  ASSERT (absCriterion < inf);  


  // To keep absCriterion < inf
  if (! multFixed /*&& ! variance_min*/ && usePrediction)
  {  
    unordered_map<const Leaf*,Real/*dissim.target*/> leaf2target_min;  
    leaf2target_min. rehash (name2leaf. size () / 100 + 1);  // PAR
    for (Dissim& dissim : dissims)
      if (   dissim. valid ()
          && ! dissim. prediction
          && dissim. target > 0.0       
         )
      {
        const array<const Leaf*,2> leaves (dissim. getLeaves ());
        const Real target_half = 0.5 * dissim. target;
        ASSERT (target_half > 0.0);
        for (const Leaf* leaf : leaves)
          if (leaf->discernible)
          {
            const auto it = leaf2target_min. find (leaf);          
            if (it == leaf2target_min. cend ())
              leaf2target_min [leaf] = target_half;
            else
              minimize (it->second, target_half);
          }
      }
    for (const auto& it : leaf2target_min)
    {
      const Real inc = it. second;
      ASSERT (inc > 0.0);
      Leaf* leaf = var_cast (it. first);
      ASSERT (! leaf->len);
      leaf->len = inc;
    }
    if (! leaf2target_min. empty ())
    {
      Real inc = NaN;
      for (Dissim& dissim : dissims)
        if (dissim. valid ())
        {
          if (find (leaf2target_min, dissim. leaf1, inc))
            dissim. prediction += inc;
          if (find (leaf2target_min, dissim. leaf2, inc))
            dissim. prediction += inc;
        }
    }
  }
  
  
  mult_sum = 0.0;
  target2_sum = 0.0;
  absCriterion = 0.0;
  // Use Threads ??
  for (Dissim& dissim : dissims)
    setDissimMult (dissim, usePrediction);
  ASSERT (absCriterion < inf);
}



void DistTree::setDissimMult (Dissim& dissim,
                              bool usePrediction) 
{
  ASSERT (optimizable ());
  ASSERT (absCriterion < inf);  

  if (! dissim. valid ())
    return;

  if (! multFixed)
  {
    // dissim.mult
    if (dissim. indiscernible ())
      dissim. mult = inf;
    else
    { 
      const Real scale = (dissim. type == no_index ? 1.0 : dissimTypes [dissim. type]. scaleCoeff);
      ASSERT (scale > 0.0);
      dissim. mult = dist2mult ((usePrediction ? dissim. prediction : dissim. target) / scale) / sqr (scale);
      if (dissim. mult == inf)  
        dissim. mult = dist2mult (epsilon);  // PAR
      ASSERT (dissim. mult < inf);
    }
  }

  if (dissim. mult < inf)
  {
    absCriterion += dissim. getAbsCriterion ();
    mult_sum     += dissim. mult;
    target2_sum  += dissim. mult * sqr (dissim. target);
  }
  ASSERT (absCriterion < inf);
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
    section ("Optimizing arc lengths", true);

  DTNode* toSkip = nullptr;  
  DTNode* toRetain = nullptr; 
  getSkipRetain (toSkip, toRetain);

  // DTNode::index, dtNodes, arcs
  for (DiGraph::Node* node : nodes)
    static_cast <DTNode*> (node) -> index = no_index;
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
      var_cast (dtNode) -> index = arcs;
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
      if (dtNode1->index == no_index)
        continue;
      dtNode1->dissimSum += dissim * it. mult;
      for (const TreeNode* node2 : path)
      {
        const DTNode* dtNode2 = static_cast <const DTNode*> (node2);
        if (dtNode2->index != no_index)
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
    if (dtNode->index != no_index)
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
    var_cast (dtNodes [i]) -> len = max (0.0, beta [i]);
      
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
    section ("Optimizing arc lengths by quartets", true);


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
    ASSERT (belowSum >= 0.0);
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
    ASSERT (aboveSum >= 0.0);
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
      var_cast (b) -> len = a->len;
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



bool DistTree::optimizeLenWhole ()
{
  node2deformationPair. clear ();
  
  // beta
  Real covar = 0.0;
  Real predict2 = 0.0;
  for (const Dissim& dissim : dissims)  
    if (dissim. validMult ())
    {
      covar    += dissim. mult * dissim. target * dissim. prediction;
      predict2 += dissim. mult * sqr (dissim. prediction);
    }
  if (covar <= 0.0)
  {
  #if 0
    if (subDepth)
      throw runtime_error (FUNC "Strange dissimilarities");
  #endif
    return false;
  }
  if (! predict2)
    throw runtime_error (FUNC "No arcs");
  if (   covar    > 1.0 / epsilon
      && predict2 > 1.0 / epsilon
     )
    return false; 
  const Real beta = covar / predict2;
  if (isNan (beta))
    return false;
  if (! beta)
    return false;
  ASSERT (beta > 0.0);
  ASSERT (beta < inf);
      
  for (DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <DTNode*> (node);
    if (dtNode->len > 0.0)
      var_cast (dtNode) -> len *= beta;
  }

  const Real absCriterion_old = absCriterion;
  absCriterion = 0.0;
  for (Dissim& dissim : dissims)  
  {
    if (dissim. validMult ())
      dissim. prediction *= beta;
    if (dissim. mult < inf)
      absCriterion += dissim. getAbsCriterion ();
  }
  ASSERT (absCriterion < inf);
  if (! leRealRel (absCriterion, absCriterion_old, 1e-3))  // PAR
    BAD_CRITERION (optimizeLenWhole);
  
  return true;
}



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
  node2deformationPair. clear ();

  VectorPtr<DTNode> dtNodes;  dtNodes. reserve (2 * name2leaf. size ());
  for (DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <DTNode*> (node);
    if (   dtNode != root
        && ! dtNode->inDiscernible ()
       )
      dtNodes << dtNode;
    IMPLY (dtNode->inDiscernible (), ! dtNode->len);
  }
  dtNodes. sort (DTNode_len_strictlyLess);
  

  {
    const uint iters = 10;  // PAR
    Progress progIter (iters);
    FOR (uint, iter, iters)
    {
      Real absCriterion_prev = absCriterion;
      {
        Progress prog (dtNodes. size (), (dtNodes. size () >= 10000) * 100);  // PAR
      //ASSERT (! prog. active);  ??
        for (const DTNode* node : dtNodes)
        {
          prog (absCriterion2str ());
        #ifndef NDEBUG    
          const Real absCriterion_old = absCriterion;  
        #endif
        
          Real arcAbsCriterion_old = 0.0;
          WeightedMeanVar mv;
          for (const uint dissimNum : node->pathDissimNums)
          {
            const Dissim& dissim = dissims [dissimNum];
            if (dissim. validMult ())
            {
              const Real dist_hat_tails = max (0.0, dissim. prediction - node->len);
              const Real arcTarget = dissim. target - dist_hat_tails;
              mv. add (arcTarget, dissim. mult);
              arcAbsCriterion_old += dissim. getAbsCriterion ();
            }
          }
          const Real len_new = max (0.0, mv. getMean ());
          
          Real arcAbsCriterion_new = 0.0;
          for (const uint dissimNum : node->pathDissimNums)
          {
            Dissim& dissim = dissims [dissimNum];
            if (dissim. validMult ())
            {
              dissim. prediction = max (0.0, dissim. prediction - node->len + len_new);
              arcAbsCriterion_new += dissim. getAbsCriterion ();
            }
          }
          minimize (arcAbsCriterion_new, arcAbsCriterion_old);
          
          var_cast (node) -> len = len_new;
          absCriterion -= arcAbsCriterion_old;
          absCriterion += arcAbsCriterion_new;
          ASSERT (leReal (absCriterion, absCriterion_old));
          maximize (absCriterion, 0.0);
        }
      }
    
    #if 1
      progIter (absCriterion2str ());
    #else
      if (Progress::enabled ())
        cerr << '\r' << iter + 1 << " / " << (uint) iters << " " << absCriterion2str ();
    #endif
      if (absCriterion_prev / absCriterion - 1.0 <= 0.01)   // PAR
        break;
      absCriterion_prev = absCriterion;  
    }
  #if 0
    if (Progress::enabled ())
      cerr << endl;
  #endif
  }


  return finishChanges ();
}



namespace
{
  
struct Star
{
  // !nullptr
  const Steiner* center {nullptr};
  VectorPtr<DiGraph::Node> arcNodes;
    // Arc's make up a star
  Real lenSum {0.0};
  
  explicit Star (const Steiner* center_arg)
    : center   (center_arg)
    , arcNodes (center->getChildren ())
    {
      ASSERT (center);
      ASSERT (! center->isTransient ());
      arcNodes << center;
      ASSERT (arcNodes. size () >= 3);
      for (const DiGraph::Node* node : arcNodes)
        lenSum += static_cast <const DTNode*> (node) -> len;
      ASSERT (lenSum >= 0.0);
    }
   
  size_t liveNodes () const
    { size_t n = 0;
      for (const DiGraph::Node* node : arcNodes)
        if (static_cast <const DTNode*> (node) -> len)
          n++;
      return n;
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
  node2deformationPair. clear ();

  Vector<Star> stars;  stars. reserve (2 * name2leaf. size ());
  for (const DiGraph::Node* node : nodes)
  {
    const DTNode* center = static_cast <const DTNode*> (node);
    if (   center == root
        || center->asLeaf () 
        || ! center->childrenDiscernible ()
       )
      continue;
    stars << std::move (Star (center->asSteiner ()));
  }
  stars. sort (Star::strictlyLess); 


  Progress prog (stars. size ());  
#ifndef NDEBUG
  const Real absCriterion_old1 = absCriterion;
#endif
  Tree::LcaBuffer buf;
  for (const Star& star : stars)
  {
    if (star. liveNodes () <= 2)  
      continue;
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
    ASSERT (subgraph. area. size () >= 2);
    subgraph. finish ();
    subgraph. area2dissimNums ();
    subgraph. dissimNums2subPaths ();
    subgraph. qc ();    
  
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
      attr->setAll (efalse);
    }
    FFOR (size_t, objNum, subgraph. subPaths. size ())
    {
      const SubPath& subPath = subgraph. subPaths [objNum];
      const size_t wholeObjNum = subPath. dissimNum;
      const VectorPtr<TreeNode>& path = subgraph. getPath (subPath, buf);
      FFOR (size_t, i, star. arcNodes. size ())
        if (path. contains (static_cast <const TreeNode*> (star. arcNodes [i])))
      //if (star. arcNodes [i] -> pathDissimNums. containsFast (wholeObjNum))  // needs sorting ??
          (* const_static_cast <ExtBoolAttr1*> (sp [i])) [objNum] = etrue;        
      var_cast (starDs. objs [objNum]) -> mult = dissims [wholeObjNum]. mult; 
      (*targetAttr) [objNum] = dissims [wholeObjNum]. target - subPath. dist_hat_tails;
    }
    starDs. qc ();
    sp. qc ();        
    const Sample sample (starDs);
    
    L2LinearNumPrediction lr (sample, sp, *targetAttr);
    ASSERT (lr. beta. size () == arcNodes. size ());
    FFOR (size_t, attrNum, lr. beta. size ())
      lr. beta [attrNum] = static_cast <const DTNode*> (arcNodes [attrNum]) -> len;
  #if 0
    // For "!!optimizeLenNode1"
    const auto beta_init (lr. beta);
    lr. setAbsCriterion ();
    const auto absCriterion_lr_init = lr. absCriterion;
  #endif
    const bool solved = lr. solveUnconstrainedFast (nullptr, true, 10, 0.01);  // PAR
    lr. qc ();
  
    // DTNode::len
    if (solved)
    {
      FFOR (size_t, attrNum, arcNodes. size ())
        const_static_cast <DTNode*> (arcNodes [attrNum]) -> len = lr. beta [attrNum];
      const Real absCriterion_old = absCriterion;  
      subgraph. subPaths2tree ();
      if (! leRealRel (absCriterion, absCriterion_old, 1e-3))  // PAR 
        BAD_CRITERION (optimizeLenNode);
    #if 0
      { 
        cout << "!!optimizeLenNode1" << endl;  
        PRINT (absCriterion_lr_init);
        PRINT (lr. absCriterion);
        PRINT (absCriterion_old);
        PRINT (absCriterion);
        PRINT (subDepth);
        PRINT (star. center);
        PRINT (star. center->getParent ());
        cout << "beta_init:" << endl;
        FFOR (size_t, attrNum, arcNodes. size ())
          PRINT (beta_init [attrNum]);
        cout << "beta result:" << endl;
        FFOR (size_t, attrNum, arcNodes. size ())
          PRINT (lr. beta [attrNum]);
        OFStream ofs ("optimizeLenNode.dm");
        starDs. saveText (ofs);
        // 2 arcs of an equal length, and all the other arcs have a positive length close to 0: why ??
        ERROR;
      }
    #endif
      prog (absCriterion2str ()); 
    }
  }

  
  if (! leRealRel (absCriterion, absCriterion_old1, 1e-4))  // PAR
  {
    const Real absCriterion_old = absCriterion_old1;
    BAD_CRITERION (optimizeLenNode2);
  }

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
    var_cast (leaf) -> len = t / 2.0;
  
  dissim. prediction = t;

  setPaths (true);
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
    Leaf* leaf = var_cast (leaf_);
    leaf->len = 0.0;
    leaf->setParent (const_static_cast <DTNode*> (root));
  }

  FOR (size_t, i, 3)
  {
    Dissim& dissim = dissims [i];
    dissim. prediction = 0.0;
    const Real t = isNan (dissim. target) ? 0 : dissim. target;
    for (const Leaf* leaf_ : leaves)
    {
      Leaf* leaf = var_cast (leaf_);
      if (dissim. hasLeaf (leaf))
        leaf->len += t;
      else
        leaf->len -= t;
    }
  }

  for (const Leaf* leaf_ : leaves)
  {
    Leaf* leaf = var_cast (leaf_);
    leaf->len /= 2.0;  // >= 0 <= Dissim::target is a distance and triangle inequality 
    maximize (leaf->len, 0.0);
  }

  absCriterion = 0.0;
  FOR (size_t, i, 3)
  {
    Dissim& dissim = dissims [i];
    for (const Leaf* leaf : leaves)
      if (dissim. hasLeaf (leaf))
        dissim. prediction += leaf->len;
    absCriterion += dissim. getAbsCriterion ();
  }
  ASSERT (absCriterion < inf);
  
  cleanTopology ();

  setPaths (true);
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
  size_t arcDist_prev = 0;
  FOR_START (size_t, i, from, to)
  {
    prog (to_string (arcDist_prev));
    const DTNode* fromNode = nodeVec [i];
    Real nodeAbsCriterion_old = NaN;
    const NewLeaf nl (fromNode, q_max, nodeAbsCriterion_old);
    ASSERT (nodeAbsCriterion_old >= 0.0);
    nl. qc ();
    const DTNode* toNode = nl. location. anchor;
    ASSERT (toNode);
    const Real improvement = nodeAbsCriterion_old - nl. location. absCriterion_leaf;
    if (   improvement <= 0.0
        || fromNode->getParent () == toNode
        || fromNode->getParent () == toNode->getParent ()
       )
      continue;
    if (! Change::valid (fromNode, toNode))
      continue;
    const DistTree::TreeNode* lca = nullptr;
    const size_t arcDist = Tree::getPath (fromNode, toNode, nullptr, lca, buf). size ();
    ASSERT (arcDist > 1);
    if (verbose ())
      cout << fromNode->getLcaName () << " -> " << toNode->getLcaName () << ' ' << improvement << ' ' << arcDist << endl;
    if (arcDist < areaRadius_std)  // PAR 
      continue;
    arcDist_prev = arcDist;
    auto change = new Change (fromNode, toNode);
    change->improvement = improvement;
    change->arcDist = arcDist;
  #if 0
    if (const Leaf* leaf = change->from->asLeaf ())  
      if (leaf->name == "AY662658.1")
      {
        const ONumber on (cout, 12, true);
        PRINT (change->improvement);
      }
  #endif
    changes << change;
  }
}
            
}



void DistTree::optimizeReinsert ()
{
  ASSERT (! subDepth);


  node2deformationPair. clear ();

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
  arrayThreads (false, reinsert_thread, nodeVec. size (), results, cref (*this), cref (nodeVec));

  VectorOwn<Change> changes;  changes. reserve (256);  // PAR 
  for (const VectorPtr<Change>& threadChanges : results)
    changes << threadChanges;

  changes. sort (Change::longer); 
//PRINT (changes. size ());  

  optimizeLargeSubgraphs (& changes);
    // Invokes: applyChanges (changes, true); 
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

  
  node2deformationPair. clear ();

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
        ASSERT (bestChange->improvement > 0.0);
        changes << bestChange;
      }
      chron_getBestChange. stop ();
    }
  }

  changes. sort (Change::strictlyBetter); 
    
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
    ASSERT (to->graph);
    if (Change::valid (from, to))
      tryChange (new Change (from, to), bestChange);
  }
  
  if (bestChange)
  {
    if (verbose (1))
      cerr << "found ";
    ASSERT (bestChange->improvement > 0.0);
    return bestChange;  
  }

  return nullptr;
}



bool DistTree::applyChanges (const VectorOwn<Change> &changes,
                             bool byNewLeaf)
{ 
  ASSERT (toDelete. empty ());
  ASSERT (absCriterion >= 0.0);
  
  
  if (verbose (1))
    cout << "# Changes: " << changes. size () << endl;


  const Real absCriterion_init = absCriterion;


  // DTNode::stable: init
  for (DiGraph::Node* node : nodes)
    static_cast <DTNode*> (node) -> stable = true;


  size_t commits = 0;
  {
  //Unverbose un;
    size_t nChange = 0;
    Progress prog (changes. size ());
    for (const Change* ch_ : changes)
    {
      Change* ch = var_cast (ch_);
      ASSERT (ch);
      ASSERT (ch->improvement > 0.0);

      prog (absCriterion2str () + ifS (byNewLeaf, " " + to_string (ch->arcDist)));
      nChange++;
      
      if (! ch->valid ())
        continue;
        
      Unverbose un;  
      if (verbose ())
        ch->qc ();  
  
    #ifndef NDEBUG
      const bool first = ch_ == changes. front ();  
    #endif
  
      qcPaths (); 
        
    //const Real absCriterion_old = absCriterion;
  
      if (verbose (1))
      {
        cout << "Apply " << nChange << "/" << changes. size () << ": ";
        ch->saveText (cout);  
        cout << endl;
      }
      const bool success = ch->apply (); 
      if (verbose (1))
        cout << "Success: " << success << "  improvement = " << ch->improvement << endl;
      IMPLY (first && ! byNewLeaf, success && ch->improvement > 0.0);
        // byNewLeaf => not all dissimilarities have been used to create Change

    #if 0
      if (const Leaf* leaf = ch->from->asLeaf ())  
        if (leaf->name == "AY662658.1")
        {
          const ONumber on (cout, 12, true);
          PRINT (ch->improvement);
          PRINT (success);
        }
    #endif

      if (! success || ch->improvement <= 0.0)
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
      //ASSERT (leReal (absCriterion, absCriterion_old));  ??
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
  

  ASSERT (absCriterion < inf);
  const Real improvement = max (0.0, absCriterion_init - absCriterion);  // max ??
  if (verbose (1))
  {
    const ONumber on (cout, absCriterionDecimals, true);
    cout << "# Commits = " << commits << endl;
    cout << "Improvement = " << improvement /*<< "  from: " << absCriterion_init << " to: " << absCriterion*/ << endl;
  }
  IMPLY (! commits, ! improvement);

  if (commits)
  {
    if (byNewLeaf)
      finishChanges ();
    else
    {
      optimizeLenWhole ();
      optimizeLenArc (); 
      optimizeLenNode ();  
    }
    if (verbose (1))
      reportErrors (cout);
  }
  
  qc ();
  
  return improvement > 0.0;
}



void DistTree::tryChange (Change* ch,
                          const Change* &bestChange)
{ 
  ASSERT (absCriterion < inf);
  ASSERT (ch);
  ASSERT (ch->from->graph == this);
  ASSERT (isNan (ch->improvement));
  
  if (verbose ())
  {
    ch->saveText (cout); 
    cout << endl;
    ch->qc ();
  }

  ch->apply ();
  ASSERT (! isNan (ch->improvement));
  ch->restore ();

  Unverbose unv;
  if (verbose ())
  {
    ch->saveText (cout); 
    cout << endl;
  }
  
  if (Change::strictlyBetter (ch, bestChange))
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
                               const VectorPtr<Tree::TreeNode> possibleBoundary,  // needs a copy
                               const VectorOwn<Change>* changes)
{ 
  image. processLarge (subTreeRoot, possibleBoundary, changes); 
}

}
                        


void DistTree::optimizeLargeSubgraphs (const VectorOwn<Change>* changes)
{
  ASSERT (threads_max);
  
  
  node2deformationPair. clear ();

  const bool threadsUsed = (threads_max > 1 && ! subDepth);
  
//unique_ptr<const Chronometer_OnePass> cop (subDepth == 1 ? new Chronometer_OnePass ("optimizeLargeSubgraphs for " + toString (name2leaf. size ())) : nullptr);   
  
  constexpr size_t large_min = 1000;  // PAR  
  ASSERT (large_min > 1);

  const size_t largeParts = name2leaf. size () / large_min;

  if (! largeParts)  
  {
    if (changes)
      applyChanges (*changes, true); 
    else
      optimizeSmallSubgraphs (areaRadius_std);
    return;
  }
    
  size_t parts_max = threads_max;  
  if (! threadsUsed)
    parts_max = (size_t) ceil (6.0 * log ((Real) largeParts + 0.5)) + 1;  // PAR  
  minimize (parts_max, largeParts + 1);
  ASSERT (parts_max >= 2);

  
  VectorPtr<Steiner> boundary;  boundary. reserve (parts_max);   // isTransient()
  {
    VectorPtr<Steiner> cuts;  cuts. reserve (parts_max);
    {
      size_t goalSize = 0;
      const size_t rootNodes = nodes. size () - 1;
      {
        size_t indiscernibleChldren_max = 0;
        for (const DiGraph::Node* node : nodes)
          if (const Steiner* st = static_cast <const DTNode*> (node) -> asSteiner ())
            if (! st->childrenDiscernible ())
              maximize (indiscernibleChldren_max, st->arcs [false]. size ());
        ASSERT (rootNodes > large_min);
        constexpr Prob topSubgraphFrac = 0.2;  // PAR  
          // Top subtree is slower
        const Real goalSize_init = (Real) rootNodes / ((Real) parts_max - 1 + topSubgraphFrac);  
        goalSize = max ({large_min, indiscernibleChldren_max, (size_t) goalSize_init});
      }
      ASSERT (goalSize);
      for (;;)
      {
        const_static_cast <Steiner*> (root) -> subtreeSize2leaves ();  
          // Node::leaves = number of Arc's in the subtree
        ASSERT (root->leaves <= rootNodes);
        IMPLY (cuts. empty (), root->leaves == rootNodes);
        // Minimize maximum subtree size
        while (root->leaves > goalSize + 2)  // root->leaves - 2 > goalSize
        {
          const Steiner* st_best = nullptr;
          {
            size_t diff = numeric_limits<size_t>::max ();
            for (const DiGraph::Node* node : nodes)
              if (const Steiner* st = static_cast <const DTNode*> (node) -> asSteiner ())
                if (   st->leaves
                    && st != root
                    && st->leaves <= goalSize  // => st->leaves < root->leaves - 2
                   )
                {
                  ASSERT (root->leaves > st->leaves + 2);  // Each subgraph must be smaller than *this, otherwise there can be a recursion cycling
                  if (minimize (diff, goalSize - st->leaves))
                    st_best = st;
                }
          }
          if (! st_best)
            throw runtime_error (FUNC "Too few interior nodes");  // Try optimize anyway ??
          cuts << st_best;
          // Node::leaves
          const size_t leaves_best = st_best->leaves;
          ASSERT (leaves_best);
          const Steiner* st = st_best;
          while (st)
          {
            if (st->leaves)
            {
              ASSERT (st->leaves >= leaves_best);
              var_cast (st) -> leaves -= leaves_best;
            }
            st = static_cast <const Steiner*> (st->getParent ());
          }
          var_cast (st_best) -> setLeaves (0);
        }
        if (! threadsUsed || cuts. size () <= parts_max - 1)
          break;
        goalSize += (size_t) ((Real) goalSize * 0.01) + 1;  // PAR  
        cuts. clear ();
      }
    }
    if (cuts. empty ()) 
    {
      if (changes)
        applyChanges (*changes, true); 
      else
        optimizeSmallSubgraphs (areaRadius_std);
      return;
    }
    ASSERT (boundary. empty ());
    for (const Steiner* cut : cuts)
    {
      ASSERT (cut);
      ASSERT (cut != root);
      auto inter = new Steiner (*this, const_static_cast <Steiner*> (cut->getParent ()), 0.0);
      inter->pathDissimNums = cut->pathDissimNums;
      var_cast (cut) -> setParent (inter);
      ASSERT (inter->getParent ());
      boundary << inter;
    }
    ASSERT (boundary. size () == cuts. size ());
  }
//qc ();  // Breaks due to transients


  bool failed = false;
  {
    VectorOwn<Image> images;  images. reserve (boundary. size ());
    Image mainImage (*this);  
    {
      unique_ptr<Threads> th;
      if (threadsUsed)
        th. reset (new Threads (boundary. size ()));
      VectorPtr<Tree::TreeNode> possibleBoundary;  possibleBoundary. reserve (boundary. size ());
      Progress prog (boundary. size () + 1, ! threadsUsed);
      for (const Steiner* cut : boundary)
      {
        prog ();
        ASSERT (cut);
        ASSERT (cut->isTransient ());
        auto image = new Image (*this);
        images << image;
        ASSERT (cut->isTransient ());
        if (th. get ())
          *th << thread (processLargeImage_thread, ref (*image), cut, possibleBoundary, changes);
            // Use functional ??
        else
          image->processLarge (cut, possibleBoundary, changes);
        possibleBoundary << cut;
      }   
      // Top subgraph
      prog ();
      {
        Unverbose unv;
        mainImage. processLarge (nullptr, possibleBoundary, changes);
      }
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
        prog (absCriterion2str () + " (approx.)" /*+ " " + to_string (image->subgraph. subPaths. size ())*/); 
      }
    }
  }
  if (failed)
    throw runtime_error (FUNC "OUT OF MEMORY");   // if threads_num > 1 then try thread-free execution ??
    
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
      delayDeleteRetainArcs (var_cast (st));
    }
    else
      var_cast (st) -> stable = false;
  }
  toDelete. deleteData ();
  
  setPredictionAbsCriterion ();
  qc ();
  qcPaths ();
  
  if (! subDepth)
    section ("Optimizing cut nodes", false);
  optimizeSmallSubgraphsUnstable (areaRadius_std);  // PAR
  qc ();
  qcPredictionAbsCriterion ();
}



void DistTree::optimizeSmallSubgraphs (uint areaRadius)
{
  ASSERT (areaRadius >= 1);
  ASSERT (unstableCut. empty ());
  
  node2deformationPair. clear ();

  for (DiGraph::Node* node : nodes)
    static_cast <DTNode*> (node) -> stable = false;
    
  unstableCut. reset (name2leaf. size ());
  unstableCut. insert (static_cast <const Steiner*> (root));

  Progress prog;
  const size_t unstables_max = nodes. size () - name2leaf. size ();
  unstableProcessed = 0;
  size_t unstables_prev = numeric_limits<size_t>::max ();  // valid if qc_on
  while (! unstableCut. empty ())
  {
    size_t unstables = 0;  // valid if qc_on
    if (qc_on)
    {
      unstableCut. qc ();
      size_t inCut = 0;
      for (const DiGraph::Node* node_ : nodes)
      {
        const DTNode* node = static_cast <const DTNode*> (node_);
        if (! node->stable && node->asSteiner ())
        {
          unstables++;
          if (! node->getParent () || static_cast <const DTNode*> (node->getParent ()) -> stable)
            inCut++;
        }
        if (node != root && node->stable)
        {
          const DTNode* parent = static_cast <const DTNode*> (node->getParent ());
          ASSERT (parent);
          QC_ASSERT (parent->stable);
        }
      }
      for (const Steiner* st : unstableCut. getVec ())
      {
        QC_ASSERT (! st->stable);
        QC_ASSERT (! st->getParent () || static_cast <const DTNode*> (st->getParent ()) -> stable);
      }
      QC_ASSERT (inCut == unstableCut. size ());
      QC_ASSERT (unstableCut. size () <= unstables);
    }
    const Steiner* center = unstableCut. getVec (). getRandom (rand);
    {
      Unverbose unv;
      optimizeSmallSubgraph (center, areaRadius);
    }    
    ASSERT (unstableProcessed <= unstables_max);
    prog (to_string (unstableProcessed) + " / " + to_string (unstables_max) + " " + absCriterion2str ());
    // The number of un-stable DTNode's decreases at least by 1
    if (qc_on)
    {
      QC_ASSERT (unstables_prev > unstables);
      unstables_prev = unstables;
    }
  }
  ASSERT (unstableProcessed == unstables_max);
  ASSERT (unstableCut. empty ());
}



void DistTree::optimizeSmallSubgraphsUnstable (uint areaRadius)
{
  ASSERT (areaRadius >= 1);
  
  node2deformationPair. clear ();

  Progress prog;
  for (;;)
  {
    const Steiner* center = nullptr;
    for (const DiGraph::Node* node : nodes)
      if (const Steiner* st = static_cast <const DTNode*> (node) -> asSteiner ())
        if (   ! st->stable
            && (! st->getParent () || static_cast <const DTNode*> (st->getParent ()) -> stable)
           )
        {
          center = st;
          break;
        }
    if (! center)
      break;
    {
      Unverbose unv;
      optimizeSmallSubgraph (center, areaRadius);
    }
    prog (absCriterion2str ());
    // The number of un-stable DTNode's decreases at least by 1
  }
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

  if (! subDepth && ! leRealRel (absCriterion, absCriterion_old, 1e-4))  // PAR
    BAD_CRITERION (optimizeSmallSubgraph);
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
      const Vector<uint> lcaObjNums (node->getLcaDissimNums ());
      for (const size_t dissimNum : lcaObjNums)
      {
        Dissim& dissim = dissims [dissimNum];
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
        && s->childrenDiscernible ()
        && ! s->len
       )
    {
      delayDeleteRetainArcs (var_cast (s));
      return true;
    }
      
  return false;
}



size_t DistTree::deleteQuestionableArcs (Prob arcExistence_min)
{
  ASSERT (isProb (arcExistence_min));
  ASSERT (optimizable ());
  
  if (! arcExistence_min)
    return 0;
  
  Heap<Steiner> heap (Steiner::arcExistence_compare, Steiner::arcExistence_index, nodes. size ());
  for (const DiGraph::Node* node : nodes)
    if (const Steiner* st = static_cast <const DTNode*> (node) -> asSteiner ())
      if (   st->getParent () 
          && st->childrenDiscernible ()
         )
      {
        var_cast (st) -> arcExistence = st->getArcExistence ();
        var_cast (st) -> heapIndex = no_index;
        if (st->arcExistence >= arcExistence_min)
          continue;
        heap << var_cast (st);
      }
  
  size_t n = 0;
  {
    Progress prog (heap. size ());
    while (! heap. empty ())
    {
      const Steiner* st = heap. getMaximum ();
      ASSERT (st);
      if (st->arcExistence >= arcExistence_min)
        break;
      heap. deleteMaximum ();
      prog ();
      n++;
      const VectorPtr<DiGraph::Node> children (st->getChildren ());
      delayDeleteRetainArcs (var_cast (st));
      for (const DiGraph::Node* child : children)
        if (const Steiner* childSt = static_cast <const DTNode*> (child) -> asSteiner ())
          if (childSt->graph && childSt->heapIndex != no_index)
          {
          #ifndef NDEBUG
            const Prob arcExistence_old = childSt->arcExistence;
          #endif
            var_cast (childSt) -> arcExistence = childSt->getArcExistence ();
            ASSERT (arcExistence_old <= childSt->arcExistence);
            heap. decreaseKey (childSt->heapIndex);
          }
    }
  }
      
  toDelete. deleteData ();

  return n;  
}



namespace
{

  
struct DissimCoeffFunc final : Func1
{
  // Input
  const Vector<Real>& covar;
  const Vector<Real>& predict2;
  // Output
  Vector<Real> beta;
  
  
  DissimCoeffFunc (const Vector<Real>& covar_arg,
                   const Vector<Real>& predict2_arg)
    : covar (covar_arg)
    , predict2 (predict2_arg)
    , beta (covar_arg. size (), 0.0)
    {
      ASSERT (covar. size () == predict2. size ());
    }  
    
  Real f (Real lambda) final
    {
      Real sum = 0.0;
      FFOR (size_t, i, covar. size ())
        if (covar [i] > 0.0)
        {
          ASSERT (predict2 [i] > 0.0);
          beta [i] = (predict2 [i] + lambda) / covar [i];
          ASSERT (beta [i] >= 0.0);
          ASSERT (beta [i] < inf);
          sum += log (beta [i]);
        }
        else
          beta [i] = 0.0;
      return  sum;
    }
};

  
}



void DistTree::optimizeDissimCoeffs ()
{
  ASSERT (optimizable ());


  if (dissimTypes. empty ())
    return;
  ASSERT_EQ (getDissimCoeffProd (), 1.0, dissimCoeffProd_delta); 

    
  // Linear regression
  Vector<Real> covar    (dissimTypes. size (), 0.0);  
  Vector<Real> predict2 (dissimTypes. size (), 0.0);  
  for (const Dissim& dissim : dissims)  
    if (dissim. validMult ())
    {
      const Real mult = dissim. mult * sqr (dissimTypes [dissim. type]. scaleCoeff); 
      covar    [dissim. type] += mult * dissim. target * dissim. prediction;
      predict2 [dissim. type] += mult * sqr (dissim. prediction);
    }
    
  bool removed = false;
  FFOR (size_t, type, dissimTypes. size ())
    if (covar [type] <= 0.0)
    {
      removeDissimType (type);
      removed = true;
    }
  if (removed && ! getDissimConnected ())
    throw runtime_error (FUNC "Disconnected objects"); 


  DissimCoeffFunc func (covar, predict2);
  {
    Real lambda_min = -inf;
    FFOR (size_t, type, dissimTypes. size ())  
      if (dissimTypes [type]. scaleCoeff)
        maximize (lambda_min, - predict2 [type]);
    ASSERT (lambda_min < 0.0);
    
    Real lambda_max = 1.0; 
    while (func. f (lambda_max) <= 0.0)
      lambda_max *= 2.0;
      
    const Real lambda_best = func. findZero (lambda_min, lambda_max, dissimCoeffProd_delta);  
    if (lambda_best == lambda_max)
      throw runtime_error (FUNC "Cannot optimize");    
    ASSERT (lambda_best > lambda_min);
  }


  FFOR (size_t, type, func. beta. size ())
  {
    const Real fix = func. beta [type];
    ASSERT (fix >= 0.0);
    dissimTypes [type]. scaleCoeff *= fix;  
    dissimTypes [type]. qc ();
  }


  ASSERT_EQ (getDissimCoeffProd (), 1.0, 1e-3);  // PAR
  normalizeDissimCoeffs ();


#ifndef NDEBUG
  const Real absCriterion_old = absCriterion;
#endif
  mult_sum = 0.0;
  target2_sum = 0.0;
  absCriterion = 0.0;
  for (Dissim& dissim : dissims)
    if (dissim. validMult ())
    {
      if (const Real fix = func. beta [dissim. type])
      {
        dissim. target *= fix;
        dissim. mult /= sqr (fix);
      }
      else
        dissim. mult = 0.0;
      dissim. qc ();
      if (dissim. mult)
      {
        mult_sum     += dissim. mult;
        target2_sum  += dissim. mult * sqr (dissim. target);  
        absCriterion += dissim. getAbsCriterion ();
      }
    }    


  cerr << absCriterion2str () << endl;
  ASSERT (leRealRel (absCriterion / absCriterion_old, 1.0, 1e-3));  // PAR
  qcPredictionAbsCriterion ();
}



Real DistTree::normalizeDissimCoeffs ()
{
  ASSERT (! dissimTypes. empty ());
  
  Real coeffProd = 1.0;
  size_t n = 0;
  for (const DissimType& dt : dissimTypes)
    if (dt. scaleCoeff > 0.0)
    {
      coeffProd *= dt. scaleCoeff;
      n++;
    }
  if (! n)
    throw runtime_error (FUNC "all coefficients are 0");
  const Real fix = pow (coeffProd, 1.0 / (Real) n);
  ASSERT (fix > 0.0);
  
  const Real multiplier = 1.0 / fix;
  ASSERT (multiplier > 0.0);
  for (DissimType& dt : dissimTypes)
    dt. scaleCoeff *= multiplier;

  ASSERT_EQ (getDissimCoeffProd (), 1.0, dissimCoeffProd_delta); 
  
  return multiplier;
}



void DistTree::removeDissimType (size_t type)  
{ 
  ASSERT (dissimTypes [type]. scaleCoeff > 0.0);
  
  dissimTypes [type]. scaleCoeff = 0.0;  
  const Real multiplier = normalizeDissimCoeffs ();
  ASSERT (multiplier > 0.0);
  if (multiplier == 1.0)
    return;
  
  mult_sum = 0.0;
  target2_sum = 0.0;
  absCriterion = 0.0;
  for (Dissim& dissim : dissims)
    if (dissim. validMult ())
    {
      if (dissim. type == type)
        dissim. mult = 0.0;
      if (dissim. mult)
      {
        dissim. target     *= multiplier;
        dissim. prediction *= multiplier;
        dissim. mult       /= sqr (multiplier);
        mult_sum     += dissim. mult;
        target2_sum  += dissim. mult * sqr (dissim. target);  
        absCriterion += dissim. getAbsCriterion ();
      }
    }    
  ASSERT (absCriterion < inf);

  for (DiGraph::Node* node : nodes)
    static_cast <DTNode*> (node) -> len *= multiplier;
}



namespace
{
  
Real setDissimWeightAttrs (const Leaf* leaf1,
                           const Leaf* leaf2,
                           Real s,
                           Real s2,
                           Real w,
                           PositiveAttr2* dissimAttr,
                           PositiveAttr2* weightAttr)
// Return: dissimTypeError
{
  ASSERT ((bool) leaf1 == (bool) leaf2);

  if (! leaf1)
    return 0.0;

  ASSERT (leaf1);
  ASSERT (leaf2);
  ASSERT (leaf1 != leaf2);
  ASSERT (s >= 0.0);
  ASSERT (w >= 0.0);
//ASSERT (w < inf);
  ASSERT (dissimAttr);
  ASSERT (weightAttr);
  ASSERT (& dissimAttr->ds == & weightAttr->ds);
  
  if (! w)
    return 0.0;
    
  const Real ave = s / w;
  const size_t objNum1 = dissimAttr->ds. getName2objNum (leaf1->name);
  const size_t objNum2 = dissimAttr->ds. getName2objNum (leaf2->name);
  ASSERT (isNan (dissimAttr->get (objNum1, objNum2)));
  ASSERT (isNan (weightAttr->get (objNum1, objNum2)));
  dissimAttr->putSymm (objNum1, objNum2, ave);
  weightAttr->putSymm (objNum1, objNum2, w);
  
  if (! s2 && ! s)
    return 0.0;
  
  const Real var = s2 / w - sqr (ave);
  return var * w;
}                        
  
}



Dataset DistTree::getDissimWeightDataset (Real &dissimTypeError) const
{
  ASSERT (optimizable ());
  ASSERT (dissims. ascending == etrue);
  
  
  Dataset ds;
  ds. objs. reserve (name2leaf. size ());  
  for (const auto& it : name2leaf)
    ds. appendObj (it. first);
  ds. setName2objNum ();
  
  auto dissimAttr_ = new PositiveAttr2 ("dissim", ds, dissimDecimals);  
  auto weightAttr_ = new PositiveAttr2 ("weight", ds, dissimDecimals);  
  
  dissimAttr_->setDiag (0.0);
  weightAttr_->setDiag (0.0);
    
  dissimTypeError = 0.0;
  const Leaf* leaf1_old = nullptr;
  const Leaf* leaf2_old = nullptr;
  Real s  = 0.0;
  Real s2 = 0.0;
  Real w  = 0.0;
  for (const Dissim& dissim : dissims)    
    if (dissim. valid ())
    {
      if (   leaf1_old != dissim. leaf1
          || leaf2_old != dissim. leaf2
         )
      {
        dissimTypeError += setDissimWeightAttrs (leaf1_old, leaf2_old, s, s2, w, dissimAttr_, weightAttr_);
        s  = 0.0;
        s2 = 0.0;
        w  = 0.0;
        leaf1_old = dissim. leaf1;
        leaf2_old = dissim. leaf2;
      }
      if (dissim. mult < inf)
      {
        s  += dissim. mult *      dissim. target;
        s2 += dissim. mult * sqr (dissim. target);
      }
      w  += dissim. mult;
    }  
  dissimTypeError += setDissimWeightAttrs (leaf1_old, leaf2_old, s, s2, w, dissimAttr_, weightAttr_);
  
  ASSERT (leReal (dissimTypeError, absCriterion));
    

  return ds;  
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
  ASSERT (! multAttr);
  
  node2deformationPair. clear ();

  const TreeNode* parent = leaf->getParent ();
  if (! parent)
    throw runtime_error (FUNC "removeLeaf: Empty tree");
  ASSERT (parent->graph);
    
  leaf->detachChildrenUp ();  
  ASSERT (! leaf->graph);
  detachedLeaves << leaf;  

  name2leaf. erase (leaf->name);
  ASSERT (name2leaf. size () >= 1);
  if (name2leaf. size () == 1)
    throw runtime_error (FUNC "one leaf tree");
  
  // Clean topology, parent, Leaf::discernible consistency
  ASSERT (! parent->isLeaf ());
  if (const TreeNode* child = parent->isTransient ())
  {
    if (const Leaf* childLeaf = static_cast <const DTNode*> (child) -> asLeaf ())
      var_cast (childLeaf) -> discernible = true;
    const Steiner* st = static_cast <const DTNode*> (parent) -> asSteiner ();
    ASSERT (st);
    parent = parent->getParent ();
    delayDeleteRetainArcs (var_cast (st));
    if (! parent)
      parent = root;
    ASSERT (parent->graph);
  }  
  else if (! static_cast <const DTNode*> (parent) -> childrenDiscernible () && optimizable ())
  {
    const Steiner* st = static_cast <const DTNode*> (parent) -> asSteiner ();
    ASSERT (st);
    const Cluster2Leaves cluster2leaves (var_cast (st) -> getIndiscernibles ());
    ASSERT (! cluster2leaves. empty ());
    if (cluster2leaves. size () > 1)
    {
      VectorPtr<Leaf> cluster_old;  cluster_old. reserve (parent->arcs [false]. size ());
      for (const DiGraph::Arc* arc : parent->arcs [false])
      {
        const Leaf* child_ = static_cast <const DTNode*> (arc->node [false]) -> asLeaf ();
        ASSERT (child_);
        var_cast (child_) -> discernible = true;
        cluster_old << child_;
      }
      leafCluster2discernibles (cluster2leaves);
      ASSERT (! cluster_old. empty ());
      cluster_old. sort ();
      ASSERT (cluster_old. isUniq ());
      for (const Leaf* child_ : cluster_old)
        for (const uint dissimNum : child_->pathDissimNums)
        {
          Dissim& dissim = dissims [dissimNum];
          if (   dissim. valid ()
              && cluster_old. contains (dissim. getOtherLeaf (child_))
              && ! dissim. indiscernible ()
              && dissim. mult == inf
             )
          {
            setDissimMult (dissim, false);  // !dissim.prediction
            ASSERT (dissim. mult < inf);
          }
        }
      EXEC_ASSERT (parent = cluster_old [0] -> getParent ());
    }
    ASSERT (parent->graph);
  }

      
  if (optimizable ())
  {
    for (const uint dissimNum : leaf->pathDissimNums)
    {
      Dissim& dissim = dissims [dissimNum];
      ASSERT (dissim. hasLeaf (leaf));
      ASSERT ( ! dissim. valid ());
      ASSERT (dissim. mult >= 0.0);
      if (dissim. mult < inf)
      {
        absCriterion -= dissim. getAbsCriterion ();       
        mult_sum     -= dissim. mult;
        target2_sum  -= dissim. mult * sqr (dissim. target); 
      }
      dissim. mult = 0.0;
    }
    ASSERT (absCriterion < inf);
    maximize (absCriterion, 0.0);
  //ASSERT (target2_sum >= absCriterion);  // Can occur just after neighbor joining
    maximize (mult_sum, 0.0);
  
    qcPaths (); 
  
    if (optimizeP)
    {
      Unverbose unv;
      ASSERT (parent);
      ASSERT (parent->graph);
      optimizeSmallSubgraph (static_cast <const DTNode*> (parent) -> asSteiner (), /*2 * */ areaRadius_std);  // PAR 
    }
  }


  toDelete. deleteData (); 
}



void DistTree::reroot (DTNode* underRoot,
                       Real arcLen) 
{
  ASSERT (underRoot);
  ASSERT (& underRoot->getTree () == this);
  ASSERT (! underRoot->inDiscernible ());

  
//node2deformationPair. clear ();


  if (underRoot != root)
  {
    ASSERT (arcLen >= 0.0);
    ASSERT (arcLen <= underRoot->len);
    
    DTNode* root_ = const_static_cast<DTNode*> (root);
      
    const Steiner* newRootParent = static_cast <const DTNode*> (underRoot->getParent ()) -> asSteiner ();
    auto newRoot = new Steiner (*this, var_cast (newRootParent), underRoot->len - arcLen);
    // Arc-specific data
    newRoot->pathDissimNums = underRoot->pathDissimNums;
    newRoot->errorDensity   = underRoot->errorDensity;
    //
    underRoot->setParent (newRoot); 
    underRoot->len = arcLen;
    
    newRoot->makeDTRoot ();
    ASSERT (newRoot == root);
    ASSERT (root_ != root);
    
    setLca ();  

    if (root_->isTransient ())
      delayDeleteRetainArcs (root_);

    finishChanges ();
  }


  sort ();
  
  qcPaths ();  
}



Real DistTree::reroot (bool topological)
{
  DTNode* root_ = const_static_cast<DTNode*> (root);
  
  if (! root_->childrenDiscernible ())
    return 0.0;
  
  
  DTNode* bestDTNode = nullptr;
  Real bestDTNodeLen_new = NaN;

  root_->setSubtreeLenUp (topological);
//cout << "Old radius: " << root_->getHeight_ave () << endl;  
  
  WeightedMeanVar bestGlobalLen;
  root_->setGlobalLenDown (topological, bestDTNode, bestDTNodeLen_new, bestGlobalLen);
  ASSERT (bestDTNode);
  ASSERT (bestDTNode != root);
  
  const Real height = bestGlobalLen. getMean ();
  
  if (verbose ())
  {
    PRINT (bestDTNode);
    PRINT (bestDTNodeLen_new);
    PRINT (bestGlobalLen. getMean ());
    saveText (cout);
  }
  
  clearSubtreeLen ();


  reroot (bestDTNode, bestDTNodeLen_new);
    
  
  return height;
}



Real DistTree::getMeanResidual () const
{
  ASSERT (optimizable ());

  Real s = 0.0;
  for (const Dissim& dissim : dissims)
    if (dissim. validMult ())
      s += dissim. mult * dissim. getResidual ();
  ASSERT (! isNan (s));
  
  return s / mult_sum;
}



Real DistTree::getMinLeafLen () const
{
  Real len_min = inf;
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
    if (dissim. validMult ())
      corr. add (dissim. target, sqr (dissim. getResidual ()));
  
  return corr. getCorrelation ();
}
  
  

Real DistTree::getUnoptimizable () const
{
  ASSERT (optimizable ());

  Real epsilon2_0 = 0.0;
  for (const Dissim& dissim : dissims)
    if (   dissim. validMult ()
        && dissim. indiscernible ()
       )
    {
      ASSERT (! dissim. prediction);
      epsilon2_0 += dissim. mult * sqr (dissim. target);
    }

  return epsilon2_0;
}



namespace
{

void setErrorDensity_array (size_t from,
                            size_t to,
                            Notype /*&res*/,
                            const VectorPtr<DiGraph::Node> &nodeVec,
                            Real c)
{
  FOR_START (size_t, i, from, to)
  {
    DTNode* dtNode = const_static_cast <DTNode*> (nodeVec [i]);
    dtNode->setErrorDensity (c);
  }
}

}



void DistTree::setErrorDensities () 
{
  ASSERT (optimizable ());

  size_t n = 0;
  for (const Dissim& dissim : dissims)
    if (dissim. validMult ())
      n++;
  ASSERT (n);

  const Real c = absCriterion / (Real) n;
  ASSERT (c >= 0.0);

  VectorPtr<DiGraph::Node> nodeVec;  nodeVec. reserve (nodes. size ());
  for (DiGraph::Node* node : nodes)  
    nodeVec << node;
  nodeVec. randomOrder ();

  vector<Notype> notypes;
  arrayThreads (false, setErrorDensity_array, nodeVec. size (), notypes, cref (nodeVec), c);
}



void DistTree::setLeafNormCriterion ()  
{
  ASSERT (optimizable ());
  
  for (DiGraph::Node* node : nodes)
    if (Leaf* leaf = var_cast (static_cast <DTNode*> (node) -> asLeaf ()))
      leaf->normCriterion = 0.0;  // temporary

  size_t n = 0;
  for (const Dissim& dissim : dissims)
    if (dissim. validMult ())
    {
      const Real criterion = dissim. getAbsCriterion ();
      ASSERT (criterion < inf);
      ASSERT (criterion >= 0.0);
      
      Leaf* leaf1 = var_cast (dissim. leaf1);
      Leaf* leaf2 = var_cast (dissim. leaf2);

      leaf1->normCriterion += criterion; 
      leaf2->normCriterion += criterion; 
      
      n++;
    }
  ASSERT (n);

  const Real c = absCriterion / (Real) n;
  ASSERT (c >= 0.0);
  for (DiGraph::Node* node : nodes)
    if (const Leaf* leaf = static_cast <DTNode*> (node) -> asLeaf ())
      var_cast (leaf) -> normCriterion = (leaf->normCriterion / c - (Real) leaf->pathDissimNums. size ()) / sqrt (2.0 * (Real) leaf->pathDissimNums. size ());
}



namespace
{

void setNodeMaxDeformationDissimNum_array (size_t from,
                                           size_t to,
                                           Notype /*&res*/,
                                           const VectorPtr<DiGraph::Node> &nodeVec,
                                           const DistTree &tree)
{
  FOR_START (size_t, i, from, to)
  {
    DTNode* dtNode = const_static_cast <DTNode*> (nodeVec [i]);
    dtNode->maxDeformationDissimNum = dissims_max;
    Real target = 0.0;
    for (const uint dissimNum : dtNode->pathDissimNums)
    {
      const Dissim& dissim = tree. dissims [dissimNum];
      if (dissim. validMult ())
        if (maximize (target, dissim. getDeformation ()))
          dtNode->maxDeformationDissimNum = dissimNum;
    }
  }
}

}



void DistTree::setNodeMaxDeformationDissimNum ()
{
  ASSERT (optimizable ());
  
  VectorPtr<DiGraph::Node> nodeVec;  nodeVec. reserve (nodes. size ());
  for (DiGraph::Node* node : nodes)  
    nodeVec << node;
  nodeVec. randomOrder ();

  vector<Notype> notypes;
  arrayThreads (false, setNodeMaxDeformationDissimNum_array, nodeVec. size (), notypes, cref (nodeVec), cref (*this));
}



Real DistTree::getDeformation_mean () const
{
  ASSERT (optimizable ());

  Real deformation_mean = NaN;
  
  MeanVar mv;
  for (const Dissim& dissim : dissims)
    if (dissim. validMult ())
      mv << dissim. getDeformation ();
  deformation_mean = mv. getMean ();
  ASSERT (deformation_mean >= 0.0);
  
  return deformation_mean;
}



Dataset DistTree::getLeafErrorDataset (bool criterionAttrP,
                                       Real deformation_mean) const
{
  ASSERT (optimizable ());
  
  Dataset ds;
  ds. objs. reserve (name2leaf. size ());  
  PositiveAttr1* criterionAttr   = nullptr;
  PositiveAttr1* deformationAttr = nullptr;
  if (criterionAttrP)
    criterionAttr   = new PositiveAttr1 ("leaf_error",  ds);  
  if (! isNan (deformation_mean))
    deformationAttr = new PositiveAttr1 ("deformation", ds);  
  for (const auto& it : name2leaf)
  {
    const Leaf* leaf = it. second;
    ASSERT (leaf);
    if (leaf->graph)
    {
      const size_t index = ds. appendObj (it. first);
      if (criterionAttr)
        (*criterionAttr)   [index] = leaf->normCriterion;
      if (deformationAttr)
        (*deformationAttr) [index] = leaf->getDeformation () / deformation_mean;
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
  return a->normCriterion > b->normCriterion;
}

}



VectorPtr<Leaf> DistTree::findCriterionOutliers (const Dataset &leafErrorDs,
                                                 Real outlier_EValue_max,
                                                 Real &outlier_min_excl) const
{
  const Attr* attr = leafErrorDs. name2attr ("leaf_error");
  ASSERT (attr);
  const PositiveAttr1* criterionAttr = attr->asPositiveAttr1 ();
  ASSERT (criterionAttr);
  
  Normal distr;  
  const Sample sample (leafErrorDs);
  outlier_min_excl = criterionAttr->locScaleDistr2outlier (sample, distr, true, outlier_EValue_max);

  VectorPtr<Leaf> res;
  if (! isNan (outlier_min_excl))
    for (const auto& it : name2leaf)
    {
      const Leaf* leaf = it. second;
      if (   leaf->graph
          && ! leaf->good
          && leaf->normCriterion > outlier_min_excl
         )
        res << leaf;
    }
    
  res. sort (leafRelCriterionStrictlyGreater);
  
  // Find a maximal "independent" subset ??
      
  return res;
}



namespace
{
  
bool leafDeformationStrictlyGreater (const Leaf* a,
                                     const Leaf* b)
{
  ASSERT (a);
  ASSERT (b);
  return a->getDeformation () > b->getDeformation ();
}

}



VectorPtr<Leaf> DistTree::findDeformationOutliers (Real deformation_mean,
                                                   Real outlier_EValue_max,
                                                   Real &outlier_min_excl) const
{
  ASSERT (deformation_mean >= 0.0);
  ASSERT (outlier_EValue_max >= 0.0);

  VectorPtr<Leaf> res;

  if (! deformation_mean)
    return res;

  // outlier_min_excl
  {
    Chi2 chi2;
    chi2. setParam (1.0);
    MaxDistribution maxD;
    const size_t n = max<size_t> (1, 2 * dissims. size () / name2leaf. size ());
    maxD. setParam (& chi2, n);
    maxD. qc ();
    outlier_min_excl = maxD. getQuantileComp (1.0 - outlier_EValue_max / (Real) name2leaf. size (), 0.0, 1e6);  // PAR
  }
  ASSERT (outlier_min_excl >= 0.0)
  
  for (const auto& it : name2leaf)
  {
    const Leaf* leaf = it. second;
    if (   leaf->graph
        && ! leaf->good
        && leaf->getDeformation () / deformation_mean > outlier_min_excl
       )
      res << leaf;
  }
  res. sort (leafDeformationStrictlyGreater);  

  // Find a maximal "independent" subset ??
      
  return res;
}



namespace
{
  

struct TriangleType
{
  const Leaf* parent1;
  const Leaf* parent2;
  size_t dissimType;
  
  bool operator== (const TriangleType &other) const
    { return    parent1    == other. parent1
             && parent2    == other. parent2
             && dissimType == other. dissimType;
    }

  struct Hash
  {
    size_t operator() (const TriangleType &tt) const
      { static hash<const Leaf*> leafHash;
        static hash<size_t>      sizeHash;
        return leafHash (tt. parent1) ^ leafHash (tt. parent2) ^ sizeHash (tt. dissimType);
      }
  };
};



#if 0
??
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
  RequestCandidate () = default;
    
  bool operator< (const RequestCandidate &other) const
    { return LeafPair (leaf1, leaf2) < LeafPair (other. leaf1, other. leaf2); }
  bool operator== (const RequestCandidate &other) const
    { return LeafPair (leaf1, leaf2) == LeafPair (other. leaf1, other. leaf2); }
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
    const Leaf* leaf = leaves [i];
    RequestCandidate req;
    Real hybridness_tree = DistTree_sp::hybridness_min;
    for (const uint dissimNum1 : leaf->pathDissimNums)
    {
      const Dissim& dissim1 = leaf->getDistTree (). dissims [dissimNum1];
      if (! dissim1. validMult ())
        continue;
      FOR (uint, dissimNum2, dissimNum1)
      {
        const Dissim& dissim2 = leaf->getDistTree (). dissims [dissimNum2];
        if (! dissim2. validMult ())
          continue;
        ??
        const Tree::TreeNode* lca_ = nullptr;
        const VectorPtr<Tree::TreeNode>& path = Tree::getPath ( dissim1. getOtherLeaf (leaf)
                                                              , dissim2. getOtherLeaf (leaf)
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
    }
    if (req. leaf1)
      res << req;
  }
}
#endif



void addHybridTriangles_thread (size_t from,
                                size_t to,
                                Vector<Triangle>& res,
                                const VectorPtr<Leaf> &badLeaves)
{
  ASSERT (from <= to);
  ASSERT (to <= badLeaves. size ());
  ASSERT (res. empty ());
    
  FOR_START (size_t, i, from, to)
    badLeaves [i] -> addHybridTriangles (res);
}



void setTriangles_thread (size_t from,
                          size_t to,
                          Notype&,
                          Vector<TriangleParentPair> &triangleParentPairs_init,
                          const DistTree &tree)
{
  ASSERT (from <= to);
  ASSERT (to <= triangleParentPairs_init. size ());
    
  FOR_START (size_t, i, from, to)
    triangleParentPairs_init [i]. setTriangles (tree);
}


}  // namespace



Vector<TriangleParentPair> DistTree::findHybrids (Real dissimOutlierEValue_max,
                                                  Vector<LeafPair>* dissimRequests) const
{
  ASSERT (! subDepth);
  ASSERT (optimizable ());
  ASSERT (dissimOutlierEValue_max > 0.0);
  ASSERT (dissimParam. hybridness_min > 1.0);


//const Chronometer_OnePass cop ("findHybrids"); 

  Vector<TriangleParentPair> triangleParentPairs;

  constexpr Real hybridness_min_init = 1.0;  // PAR   
  ASSERT (dissimParam. hybridness_min > hybridness_min_init);

  for (auto& it : name2leaf)
  {
    const Leaf* leaf = it. second;
    if (! leaf->graph)
      continue;
    var_cast (leaf) -> badCriterion = 0.0;
  }

  Vector<TriangleParentPair> triangleParentPairs_init;  triangleParentPairs_init. reserve (dissims. size ());    

  // Time: O(p log(p))
  {
    // Bad dissims
    Dataset ds;
    ds. objs. reserve (dissims. size ());  
    auto criterionAttr = new PositiveAttr1 ("dissim_error", ds);  
    for (const Dissim& dissim : dissims)    
      if (dissim. validMult ())
      {
        const size_t index = ds. appendObj ();
        const Real err = dissim. getDeformation ();
        ASSERT (err >= 0.0);
        (*criterionAttr) [index] = err;
      }
    ds. qc ();
    const Sample sample (ds); 
    Chi2 chi2;
    chi2. setParam (1.0);
    chi2. qc ();
    const Real outlier_min_excl = criterionAttr->contDistr2outlier (sample, chi2, true, dissimOutlierEValue_max * (Real) dissimTypesNum () * 1e1);  // PAR
  #if 0
    Normal normal; 
    outlier_min_excl [criterionType] = criterionAttrs [criterionType] -> locScaleDistr2outlier (sample, normal, true, dissimOutlierEValue_max * (Real) dissimTypesNum () * 1e-2);  // PAR
  #endif
    for (const Dissim& dissim : dissims)    
      if (dissim. validMult ())
      {
        const Real err = dissim. getDeformation ();
        if (err <= outlier_min_excl)
          continue;
        var_cast (dissim. leaf1) -> badCriterion += dissim. getAbsCriterion ();
        var_cast (dissim. leaf2) -> badCriterion += dissim. getAbsCriterion ();
        triangleParentPairs_init << TriangleParentPair ( dissim. leaf1
                                                       , dissim. leaf2
                                                       , dissim. target
                                                       , dissim. type
                                                       );
          // Size: O(p)
      }
  }
  
  {
    // Bad leaves
  //const Chronometer_OnePass cop1 ("badLeaves"); 
    VectorPtr<Leaf> badLeaves;  badLeaves. reserve (name2leaf. size ());
    for (const auto& it : name2leaf)
    {
      const Leaf* leaf = it. second;
      if (! leaf->graph)
        continue;
      const DTNode* discernible = leaf->getDiscernible ();
      ASSERT (discernible);
      if (! discernible->len)  
        badLeaves << leaf;
    }
    Real outlier_min_excl = NaN;
  #if 1
    const Dataset leafErrorDs (getLeafErrorDataset (true, NaN));
    badLeaves << findCriterionOutliers (leafErrorDs, dissimOutlierEValue_max * 1e0, outlier_min_excl);  // PAR  
    badLeaves. sort (leafRelCriterionStrictlyGreater);
  #else  
    // Too few hybrids
    badLeaves << findDeformationOutliers (getDeformation_mean (), dissimOutlierEValue_max, outlier_min_excl); 
    badLeaves. sort (leafDeformationStrictlyGreater);
  #endif
    badLeaves. uniq ();
    const Real nLeaves = (Real) name2leaf. size ();
    const size_t badLeaves_size = min (badLeaves. size (), (size_t) (nLeaves / sqr (log (nLeaves))) + 1);  
    badLeaves. resize (badLeaves_size);
      // Size: O(n/log^2(n))
    if (! qc_on)
      badLeaves. randomOrder ();
    vector<Vector<Triangle>> resVec;
    // Time: O(n/log^2(n) * p^2/n^2 log^2(n) / threads_max) = O(p^2/n / threads_max)
    arrayThreads (false, addHybridTriangles_thread, badLeaves_size, resVec, cref (badLeaves));
    for (const Vector<Triangle>& res : resVec)
      for (const Triangle& tr : res)  
        triangleParentPairs_init << TriangleParentPair ( tr. parents [0]. leaf
                                                       , tr. parents [1]. leaf
                                                       , tr. parentsDissim ()
                                                       , tr. dissimType
                                                       );
  }
  
  {
    triangleParentPairs_init. sort ();
    triangleParentPairs_init. uniq ();
    const Real dissims_size = (Real) dissims. size ();
    const size_t triangleParentPairs_init_size = min (triangleParentPairs_init. size (), (size_t) (dissims_size / log (dissims_size)) + 1);  
    triangleParentPairs_init. resize (triangleParentPairs_init_size);
      // Size: O(p/log(p))
  }  
   
  // Time: O (p/log(p) * p/n log(n) / threads_max) = O(p^2/n / threads_max)
  {
    triangleParentPairs_init. randomOrder ();
    vector<Notype> notypes;  
    arrayThreads (false, setTriangles_thread, triangleParentPairs_init. size (), notypes, ref (triangleParentPairs_init), cref (*this));
    triangleParentPairs_init. sort ();
  }

  // triangleParentPairs_init --> triangleParentPairs
  ASSERT (triangleParentPairs. empty ());
  {
    size_t trianglesSize = 0;
    for (const TriangleParentPair& tpp : triangleParentPairs_init)
      trianglesSize += tpp. triangles. size ();
    unordered_map <TriangleType,size_t/*index in triangleParentPairs*/,TriangleType::Hash> triangleType2tpp;  triangleType2tpp. rehash (trianglesSize);
    for (const TriangleParentPair& tpp : triangleParentPairs_init)
      for (const Triangle& tr : tpp. triangles)
      {
        const TriangleType triangleType { tr. parents [0]. leaf
                                        , tr. parents [1]. leaf
                                        , tr. dissimType
                                        };
        size_t index = no_index;
        if (! find (triangleType2tpp, triangleType, index))
        {
          triangleParentPairs << TriangleParentPair ( tr. parents [0]. leaf
                                                    , tr. parents [1]. leaf
                                                    , tr. parentsDissim ()
                                                    , tr. dissimType
                                                    );
          index = triangleParentPairs. size () - 1;
          triangleType2tpp [triangleType] = index;
        }
        ASSERT (index != no_index);
        triangleParentPairs [index]. triangles << tr;
      }
  }    
  for (TriangleParentPair& tpp : triangleParentPairs)
    tpp. triangles2hybridness_ave ();
  triangleParentPairs. sort (TriangleParentPair::compareHybridness);  
  {
    const Real dissims_size = (Real) dissims. size ();
    const size_t triangleParentPairs_size = min (triangleParentPairs. size (), (size_t) (dissims_size / log (dissims_size)) + 1);  
    triangleParentPairs. resize (triangleParentPairs_size);
      // Size: O(p/log(p))
  }
    
  // Time: O(p/log(p) * p/n log(n)) = O(p^2/n)
  Set<const Leaf*> hybrids;
  {
    Set<const Leaf*> nonHybrids;
    for (TriangleParentPair& tpp : triangleParentPairs)
    {
      tpp. finish (*this, hybrids);
      tpp. qc ();
      if (tpp. triangles. empty ())
        continue;
      if (tpp. undecided ())
        continue;  
      const Set<const Leaf*> hybrids_    (tpp. getHybrids (true));
      const Set<const Leaf*> nonHybrids_ (tpp. getHybrids (false));
      IMPLY (hybrids_. empty (), nonHybrids_. empty ());
      if (   hybrids_.    intersects (nonHybrids)  // Wrongly decided
          || nonHybrids_. intersects (hybrids)
         )
        tpp. triangles. clear ();  // Preserve ??
      else
      {
        hybrids.    insertAll (hybrids_);
        nonHybrids. insertAll (nonHybrids_);
      }
    }
  }
  
  triangleParentPairs. filterValue ([] (const TriangleParentPair& tpp) { return tpp. triangles. empty () || tpp. dissimError (); });
  

  if (dissimRequests)  
  {  
  #if 0
    ??
    Vector<RequestCandidate> requests;  requests. reserve (name2leaf. size () / 10); // PAR
    triangleParentPairs_init. filterValue ([] (const TriangleParentPair& tpp) { return ! tpp. triangles. empty (); });
    triangleParentPairs_init. randomOrder ();
    {
      vector<Vector<RequestCandidate>> requests_vec;  
      arrayThreads (true, hybrid2requests, leaves. size (), requests_vec, cref (triangleParentPairs_init));
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
      if (index != no_index)
        requests [index]. inDissims = true;
    }
    for (const RequestCandidate& req : requests)
      if (! req. inDissims)
        *dissimRequests << LeafPair (req. leaf1, req. leaf2);
  #endif
  }
    
  
  return triangleParentPairs;
}



#if 0
namespace
{

struct NodeDissim
{
  const DTNode* node {nullptr};
  Real dissim_min {inf};
  
  explicit NodeDissim (const DTNode* node_arg)
    : node (node_arg)
    {
      ASSERT (node);
      const DistTree& tree = node->getDistTree ();
      for (const size_t dissimNum : node->pathDissimNums)
      {
        const Dissim& dissim = tree. dissims [dissimNum];
        if (dissim. valid ())
          minimize (dissim_min, dissim. target);
      }
    }
    
  void print (ostream &os) const
    { ONumber on (os, 6, true);
      os << node->getLcaName () << '\t' << dissim_min << '\n'; 
    }
};

}



VectorPtr<DTNode> DistTree::findOutlierArcs (Real outlier_EValue_max,
                                             Real &dissimOutlier_min_excl) const
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
    if (nd. dissim_min < inf)
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

  dissimOutlier_min_excl = lenAttr->locScaleDistr2outlier (sample, distr, true, outlier_EValue_max);

  if (isNan (dissimOutlier_min_excl))
    return res;
    
  // Mixture of 2 types of dissimilarities (> 2 ??)
  sample. mult. setAll (0.0);
//OFStream f ("arc_len");  
  FFOR (size_t, objNum, nodeDissims. size ())
    if (nodeDissims [objNum]. dissim_min >= dissimOutlier_min)
    {
      sample. mult [objNum] = 1.0;
    //nodeDissims [objNum]. print (f);  
    }
  sample. finish ();
  if (sample. mult_sum >= (Real) name2leaf. size () * 0.001)  // PAR
    dissimOutlier_min_excl = lenAttr->locScaleDistr2outlier (sample, distr, true, outlier_EValue_max);

  res. reserve (name2leaf. size () / 1000 + 1);  // PAR
  for (const NodeDissim& nd : nodeDissims)
    if (nd. dissim_min > dissimOutlier_min_excl)
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
    const Real outlier_min_excl = lenAttr->locScaleDistr2outlier (sample, distr, true, outlier_EValue_max); 
    if (isNan (outlier_min_excl))
      continue;
    FFOR (size_t, objNum, ds. objs. size ())
      if ((*lenAttr) [objNum] > outlier_min_excl)
        outliers << descendants [objNum] -> getReprLeaf (0);
  }  
  outliers. sort ();
  outliers. uniq ();
  
  return outliers;
}
  


Vector<LeafPair> DistTree::getMissingLeafPairs_ancestors (size_t depth_max,
                                                          bool refreshDissims) const
{
  ASSERT (! subDepth);

  Vector<LeafPair> pairs;  pairs. reserve (name2leaf. size () * getSparseDissims_size ());
  {
    Progress prog (name2leaf. size (), 1000);  // PAR
    for (const auto& it : name2leaf)
    {
      prog ();
      const Leaf* leaf = it. second;
      ASSERT (leaf);
      if (! leaf->graph)
        continue;
      if (! leaf->isMainIndiscernible ())
        continue;
      const VectorPtr<Leaf> matches (leaf->getSparseLeafMatches (leaf->name, depth_max, ! refreshDissims, refreshDissims));
      for (const Leaf* match : matches)
        if (leaf->getDiscernible () != match->getDiscernible ())
        {
          LeafPair p (leaf, match);
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



Vector<LeafPair> DistTree::getMissingLeafPairs_subgraphs () const
{
  ASSERT (! subDepth);

  Vector<LeafPair> pairs;  pairs. reserve (name2leaf. size ()); 
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
        const Vector<uint>& pathObjNums2 = subgraph. boundary2pathDissimNums (dtNode2);
        var_cast (pathObjNums2). sort ();
        for (const TreeNode* node1 : subgraph. boundary)
        {
          if (node1 == node2)
            break;
          const DTNode* dtNode1 = static_cast <const DTNode*> (node1);
          const Vector<uint>& pathObjNums1 = subgraph. boundary2pathDissimNums (dtNode1);
          if (pathObjNums1. intersectsFast_merge (pathObjNums2))
            continue;
          LeafPair leafPair (subgraph. getReprLeaf (dtNode1), subgraph. getReprLeaf (dtNode2));
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
      var_cast (node2->pathDissimNums). sort ();
      for (const DTNode* node1 : tooLongArcs)
      {
        ASSERT (node1 != root);
        if (node1 == node2)
          break;
        if (node1->pathDissimNums. intersectsFast_merge (node2->pathDissimNums))
          continue;
        const TreeNode* lca = getLca (node1, node2);
        ASSERT (lca);
        const DTNode* reprNode1 = lca == node1 ? static_cast <const DTNode*> (node1->getParent () -> getDifferentChild (node1)) : node1;
        const DTNode* reprNode2 = lca == node2 ? static_cast <const DTNode*> (node2->getParent () -> getDifferentChild (node2)) : node2;
        LeafPair leafPair (reprNode1->getReprLeaf (), reprNode2->getReprLeaf ());
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



Vector<LeafPair> DistTree::leaves2missingLeafPairs (const VectorPtr<Leaf> &leaves) const
{
  ASSERT (! subDepth);

  Vector<LeafPair> pairs;  pairs. reserve (leaves. size () * (leaves. size () - 1) / 2);  
  {
    Progress prog (leaves. size (), 100);  // PAR
    for (const Leaf* leaf2 : leaves)
    {
      prog ();
      for (const Leaf* leaf1 : leaves)
      {
        if (leaf1 == leaf2)
          break;
        LeafPair leafPair (leaf1, leaf2);
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
    pairs. filterValue ([this] (const LeafPair& p) { const Dissim dissim (p. first, p. second, NaN, 0.0, no_index); return dissims. containsFast (dissim); });

  return pairs;
}



#if 0
void DistTree::findTopologicalClusters () 
{
  for (DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <DTNode*> (node);
    if (Leaf* leaf = var_cast (dtNode->asLeaf ()))
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
            var_cast (leaf1) -> merge (* var_cast (leaf2));
    }
  }

  Cluster2Leaves cluster2leaves;  cluster2leaves. rehash (nodes. size ());
  for (const DiGraph::Node* node : nodes)
  {
    const DTNode* dtNode = static_cast <const DTNode*> (node);
    if (const Leaf* leaf = dtNode->asLeaf ())
      if (leaf->graph)
        cluster2leaves [var_cast (leaf) -> getDisjointCluster ()] << leaf;
  }
  ASSERT (! cluster2leaves. empty ());
  
  // ??
  OFStream f ("leaf_clusters");
  for (const auto& it : cluster2leaves)
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
  size_t i_best = no_index;
  Real diff_max = 0.0;
  FFOR (size_t, i, nodeHeights. size ())
    if (i && maximize (diff_max, nodeHeights [i]. dist - nodeHeights [i - 1]. dist))
      i_best = i;
  ASSERT (i_best != no_index);
  
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



#if 0
RealAttr1* DistTree::getResiduals2 () 
{
  ASSERT (optimizable ());

  unique_ptr<RealAttr1> resid2Attr (new RealAttr1 ("resid2", ds, target->decimals + 1));
  for (Iterator it (dsSample); it ();)  
  {
    const Real d = (*target) [*it];
  //ASSERT (d >= 0.0);
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

  unique_ptr<RealAttr1> logDiffAttr (new RealAttr1 ("logDiff", ds, target->decimals + 1));
  for (Iterator it (dsSample); it ();)  
  {
    const Real d = (*target) [*it];
  //ASSERT (d >= 0.0);
    const Real dHat = (*prediction) [*it];
    ASSERT (! isNan (dHat));
    if (   d    > 0.0
        && dHat > 0.0
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
                           bool redundantIndiscernible,
                           bool addExtra) const
{
  ASSERT (optimizable ());

  const ONumber on (os, dissimDecimals, true);
  Progress prog (dissims. size (), 1000);  // PAR
  for (const Dissim& dissim : dissims)
  {
    prog ();
    if (! dissim. valid ())
      continue;
    if (! redundantIndiscernible && dissim. redundantIndiscernible ())
      continue;
    os         << dissim. leaf1->name
       << '\t' << dissim. leaf2->name
       << '\t' << dissim. target;
    if (addExtra)
      os << '\t' << dissim. prediction 
         << '\t' << dissim. getAbsCriterion ()
         << '\t' << sqr (dissim. target - dissim. prediction);
    os << '\n';
  }
}




//////////////////////////////////////////////////////////////////////

// DissimLine

DissimLine::DissimLine (string &line,
                        uint lineNum)
{ 
  replace (line, '\t', ' ');
  name1 = findSplit (line);
  name2 = findSplit (line);
  #define MSG  getErrorStr (lineNum) + 
  if (name2. empty ())
    throw runtime_error (MSG "empty name2");
  if (name1 == name2)
    throw runtime_error (MSG "name1 == name2");
  if (name1 > name2)
    swap (name1, name2);
  dissim = str2real (line);
  if (isNan (dissim))
    throw runtime_error (MSG "dissimilarity is NaN");
  if (dissim < 0.0)
    throw runtime_error (MSG "dissimilarity is negative");
//if (! DM_sp::finite (dissim))
  //throw runtime_error (MSG "dissimilarity is infinite");
  #undef MSG
}



void DissimLine::process (const DistTree::Name2leaf &name2leaf)
{
  ASSERT (! leaf1);
  ASSERT (! leaf2);
  
  leaf1 = var_cast (findPtr (name2leaf, name1));
  if (! leaf1)
    return;
  //throw runtime_error (FUNC "Tree has no object " + name1);
  leaf1->badCriterion = -1.0;  // temporary

  leaf2 = var_cast (findPtr (name2leaf, name2));
  if (! leaf2)
    return;
  //throw runtime_error (FUNC "Tree has no object " + name2);
  leaf2->badCriterion = -1.0;  // temporary
}



void DissimLine::apply (DistTree &tree) const
{ 
  ASSERT (! isNan (dissim));
  ASSERT (dissim >= 0.0);
  if (! leaf1)
    return;
  if (! leaf2)
    return;
  if (   /*! DistTree_sp::variance_min 
      &&*/ ! dissim 
      && ! leaf1->getCollapsed (leaf2)  // Only for new Leaf's
     )  
    leaf1->collapse (leaf2);
  if (! tree. addDissim (leaf1, leaf2, dissim, 1.0/*temporary*/, no_index))
    throw runtime_error (FUNC "Cannot add dissimilarity: " + name1 + " " + name2 + " " + toString (dissim));
}


bool DissimLine::operator< (const DissimLine &other) const
{ 
  LESS_PART (*this, other, name1);
  LESS_PART (*this, other, name2);
  return false;
}




////////////////////////////////////////////////////////////////////////

// NewLeaf::Location

void NewLeaf::Location::qc () const
{
  if (! qc_on)
    return;

  QC_ASSERT (anchor);
  QC_ASSERT (! anchor->inDiscernible ());
  if (anchor->getParent ())
  {
  //QC_IMPLY (absCriterion_leaf > 0.0, anchorLen > 0.0);
    QC_ASSERT (anchorLen < inf);
    QC_ASSERT (anchor->len <= anchorLen);
    QC_ASSERT (anchor->len >= 0.0);
    QC_ASSERT (arcLen >= 0.0);
    QC_ASSERT (arcLen <= anchorLen);
  }
  else
  {
    QC_ASSERT (isNan (anchor->len));
    QC_ASSERT (isNan (anchorLen));
    QC_ASSERT (! arcLen);
  }

  QC_ASSERT (leafLen >= 0.0);
  QC_ASSERT (leafLen < inf);
  QC_IMPLY (indiscernibleFound, ! leafLen && ! arcLen);
  
//QC_IMPLY (leafLen == 0 && arcLen == 0, anchor->asLeaf ());
  
  QC_ASSERT (absCriterion_leaf >= 0.0);  
}



void NewLeaf::Location::setAbsCriterion_leaf (const NewLeaf& nl)
{
  ASSERT (leafLen >= 0.0);
  ASSERT (arcLen  >= 0.0);

  absCriterion_leaf = 0.0;
  for (const Leaf2dissim& ld : nl. leaf2dissims)
    if (ld. mult)
      absCriterion_leaf += ld. mult * sqr (ld. getEpsilon (*this));
  ASSERT (absCriterion_leaf >= 0.0);
}


    

// NewLeaf::Leaf2dissim

NewLeaf::Leaf2dissim::Leaf2dissim (const Leaf* leaf_arg,
                                   Real dissim_arg,
                                   Real mult_arg)
: leaf (leaf_arg)
, dissim (dissim_arg)
, mult (isNan (mult_arg) ? dist2mult (dissim_arg) : mult_arg)
{ 
  if (! leaf)
    throw runtime_error (FUNC "The other leaf is not found in the tree for a new leaf placement");
//ASSERT (dissim >= 0.0);
  ASSERT (! isNan (dissim));
  ASSERT (mult >= 0.0);
//ASSERT (mult < inf);
  
  // dist_hat, leafIsBelow
  const DTNode* node = leaf;
  while (node != leaf->getTree (). root)
  { 
    dist_hat += node->len;
    node = static_cast <const DTNode*> (node->getParent ());
  }
  ASSERT (dist_hat >= 0.0);
  ASSERT (dist_hat < inf);
}



void NewLeaf::Leaf2dissim::qc () const
{
  if (! qc_on)
    return;
    
  QC_ASSERT (leaf);
  
  QC_ASSERT (dissim >= 0.0);
  QC_ASSERT (DM_sp::finite (dissim));
  
  QC_ASSERT (mult >= 0.0);  
//QC_ASSERT (DM_sp::finite (mult));  // distTree_inc_place.sh can be used for an indiscernible object
}




// NewLeaf

NewLeaf::NewLeaf (const DistTree &tree_arg,
                  const string &dataDir,
                  const string &name_arg,
                  bool init)
: Named (name_arg)
, tree (tree_arg)
, location (tree_arg)
{
  ASSERT (isDirName (dataDir));  
  const string nameDir (dataDir + name + "/");
  const string dissimFName  (nameDir + "dissim");
  const string leafFName    (nameDir + "leaf");
  const string requestFName (nameDir + "request");
  process (init, dissimFName, leafFName, requestFName);
}



namespace
{
  
struct DissimMult
{
  Real dissim;
  Real mult;
  Real absCriterion;
};
  
}



NewLeaf::NewLeaf (const DTNode* dtNode,
                  size_t q_max,  
                  Real &nodeAbsCriterion_old)
: Named (checkPtr (dtNode) -> getLcaName ())
, tree  (checkPtr (dtNode) -> getDistTree ())
, node_orig (dtNode)
, location (checkPtr<const DTNode> (dtNode) -> getDistTree ())
{
  ASSERT (dtNode);
  ASSERT (dtNode->graph);
  ASSERT (dtNode != tree. root);
  ASSERT (q_max);
  ASSERT (node_orig);

#if 0
  location. anchor = static_cast <const DTNode*> (tree. root);
  location. leafLen = 0.0;
  location. arcLen = 0.0;
#endif

  // leaf2dissims[]
  {  
    unordered_map <const Leaf*, Vector<DissimMult>> leaf2dissimMults;
    leaf2dissimMults. rehash (dtNode->pathDissimNums. size ());
    {
      Vector<Tree::TreeNode::NodeDist> leafDepths;  leafDepths. reserve (tree. name2leaf. size () / 8 + 1);  // PAR
      dtNode->getLeafDepths (leafDepths);
      leafDepths. sort ();
      for (const uint dissimNum : dtNode->pathDissimNums)
      {
        const Dissim& dissim = tree. dissims [dissimNum];
        if (! dissim. validMult ())
          continue;
        size_t index = leafDepths. binSearch (Tree::TreeNode::NodeDist {dissim. leaf1, 0.0});
        const Leaf* leaf = dissim. leaf2;
        if (index == no_index)
        {
          index = leafDepths. binSearch (Tree::TreeNode::NodeDist {dissim. leaf2, 0.0});
          leaf = dissim. leaf1;
        }
        ASSERT (index != no_index);
        ASSERT (leaf);
        leaf2dissimMults [leaf] << DissimMult {max (0.0, dissim. target - leafDepths [index]. dist), dissim. mult, dissim. getAbsCriterion ()};
      }
    }
    for (const auto& it : leaf2dissimMults)
    {
      Real dissim_sum = 0.0;
      Real mult_sum   = 0.0;
      for (const DissimMult& dm : it. second)
      {
        dissim_sum += dm. mult * dm. dissim;
        mult_sum   += dm. mult;
      }
      if (mult_sum)
      {
        const Real dissim_ave = dissim_sum / mult_sum;
        Leaf2dissim ld (it. first, dissim_ave, mult_sum);
        {
          Real absCriterion = 0.0;
          Real dissim2 = 0.0;
          for (const DissimMult& dm : it. second)
          {
            absCriterion += dm. absCriterion;
            dissim2      += dm. mult * sqr (dm. dissim - dissim_ave);          
          }
          if (! geReal (absCriterion, dissim2, 1e-6))  // PAR
          {
            const Real absCriterion_old = dissim2; 
            const uint subDepth = 0;
            BAD_CRITERION (NewLeaf);
          }
          ld. absCriterion_sub = max (0.0, absCriterion - dissim2);
        }
        leaf2dissims << ld;
      }
    }
  }

  if (leaf2dissims. size () > q_max)  
  {
    leaf2dissims. sort (Leaf2dissim::dissimLess);
    leaf2dissims. resize (q_max);
    leaf2dissims. shrink_to_fit ();  
  }

  nodeAbsCriterion_old = 0.0;
  for (const Leaf2dissim& ld : leaf2dissims)
    nodeAbsCriterion_old += ld. absCriterion_sub;  
  ASSERT (nodeAbsCriterion_old >= 0.0);
  
  leaf2dissims. sort ();
  ASSERT (leaf2dissims. isUniq ());
  
  optimize ();
}



NewLeaf::NewLeaf (const DistTree &tree_arg,
                  Vector<NewLeaf::Leaf2dissim> &&leaf2dissims_arg)
: Named ("root")
, tree (tree_arg)
, location (tree_arg)
, leaf2dissims (std::move (leaf2dissims_arg))
{
  leaf2dissims. sort ();
  ASSERT (leaf2dissims. isUniq ());
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
  
#if 0
  VectorPtr<Leaf> subset;
  if (! subsetFName. empty ())
  {
    LineInput f (subsetFName);
    while (f. nextLine ())
    {
      trim (f. line);
      const Leaf* leaf = findPtr (name2leaf, f. line);
      if (! leaf)
        throw runtime_error (FUNC "Object " + strQuote (f. line) + " is not in the tree");
      subset << leaf;
    }
    subset. sort ();
    const size_t dup = subset. findDuplicate ();
    if (dup != no_index)
      throw runtime_error (FUNC "Duplicate object " + strQuote (subset [dup] -> name));    
  }
#endif
    
#if 0
  location. anchor = static_cast <const DTNode*> (tree. root);
  location. leafLen = 0.0;
  location. arcLen = 0.0;
#endif

  if (init)
  { 
    if (fileExists (leafFName))
      throw runtime_error (FUNC "File " + strQuote (leafFName) + " exists");
    if (fileExists (dissimFName))
      throw runtime_error (FUNC "File " + strQuote (dissimFName) + " exists");
    { OFStream of (dissimFName); }  // destructor creates an empty file
  }
  else
  {
    if (! fileExists (dissimFName))
      throw runtime_error (FUNC "File " + strQuote (dissimFName) + " does not exist");
    {
      LineInput f (dissimFName);
      string name1, name2, dissimS;
      Istringstream iss;  
      while (f. nextLine ())
        try
        {
          iss. reset (f. line);  
          name1. clear ();
          name2. clear ();
          dissimS. clear ();
          iss >> name1 >> name2 >> dissimS;
          QC_ASSERT (iss. eof ());
          QC_ASSERT (name1 != name2);
          QC_ASSERT (! dissimS. empty ());
          if (name1 != name)
            swap (name1, name2);
          if (name1 != name)
            throw runtime_error (FUNC + dissimFName + " must contain " + name);
          Real dissim = str2real (dissimS);
          tree. dissimParam. transform (dissim);  
          if (isNan (dissim))
            dissim = inf;  // To process "incomparable" objects by distTree_inc_add.sh
          if (dissim < 0.0)
            throw runtime_error (FUNC "Dissimilarity must be non-negative");
          leaf2dissims << Leaf2dissim (findPtr (tree. name2leaf, name2), dissim, NaN);
        }          
        catch (const exception &e)
        {          
          throw runtime_error (FUNC + f. lineStr () + " of " + dissimFName + "\n" + f. line + "\n" + e. what ());
        }
    }
    leaf2dissims. sort ();
    ASSERT (leaf2dissims. isUniq ());
    optimize ();
  }

  saveLeaf (leafFName);
  saveRequest (requestFName);
}



void NewLeaf::saveLeaf (const string &leafFName) const
{ 
  OFStream f (leafFName); 
  f << name << '\t';
  location. saveText (f);
  f << endl;
}



void NewLeaf::saveRequest (const string &requestFName) const
{ 
  OFStream f (requestFName);

  if (location. indiscernibleFound)
    return;

  VectorPtr<Leaf> requested (location. anchor->getSparseLeafMatches (name, sparsingDepth, false, false));
  requested. filterValue ([this] (const Leaf* leaf) { const Leaf2dissim ld (leaf); return leaf2dissims. containsFast (ld); });
  
  for (const Leaf* leaf : requested)
  {
    if (name == leaf->name)
      throw runtime_error (FUNC "Object " + name + " already exists in the tree");
    const string* n1 = & name;
    const string* n2 = & leaf->name;
    if (*n1 > *n2)
      swap (n1, n2);
    f << *n1 << '\t' << *n2 << '\n';
  }
}



void NewLeaf::optimize ()
{
  ASSERT (! location. indiscernibleFound);

  bool dissimExists = false;
  for (const Leaf2dissim& ld : leaf2dissims)
  {
    ASSERT (ld. mult >= 0.0);
    if (ld. mult > 0.0 /*ld. dissim != inf*/)
    {
      ASSERT (ld. dissim != inf);
      dissimExists = true;
      if (! ld. dissim)  
      {
        location. anchor = ld. leaf->getDiscernible ();
        location. leafLen = 0.0;
        location. arcLen = 0.0;
        location. absCriterion_leaf = 0.0;
        ASSERT (location. anchor);
        location. anchorLen = location. anchor->len;
        location. indiscernibleFound = true;
        {
          Tree::LcaBuffer buf;
          const Tree::TreeNode* lca = nullptr;
          for (Leaf2dissim& ld_ : leaf2dissims)
            if (ld_. leaf == ld. leaf)
              ld_. dist_hat = 0.0;
            else
            {
              const VectorPtr<Tree::TreeNode>& path = Tree::getPath (ld. leaf, ld_. leaf, nullptr, lca, buf);
              ld_. dist_hat = tree. path2prediction (path);
              ld_. leafIsBelow = false;
            }
        }
        return;
      }
    }
  }
    
  if (! dissimExists)
  {
    location. anchor = static_cast <const DTNode*> (tree. root);
    location. leafLen = inf;
    location. arcLen = inf;
    location. absCriterion_leaf = 0.0;
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
  if (location. anchorLen > 0.0)
  {
    anchor2location ();
    if (minimize (location_best. absCriterion_leaf, location. absCriterion_leaf))
    {
      location_best = location;
      leaf2dissims_best = leaf2dissims;
    }
  /*
    else if (   node_orig  // greedy search
             && location. absCriterion_leaf / location_best. absCriterion_leaf > 1.01  // PAR
            )
      return;
  */
  }

  if (! location. anchor->childrenDiscernible ())
    return;
    
  for (const DiGraph::Arc* arc : location. anchor->arcs [false])
  {
    const DTNode* child = static_cast <const DTNode*> (arc->node [false]);
    ASSERT (! child->inDiscernible ());
  //if (child == node_orig)
    //continue;
    const Keep<Location> location_old (location);
    const Keep<Vector<Leaf2dissim>> leaf2dissims_old (leaf2dissims);
    if (descend (child))  
      optimizeAnchor (location_best, leaf2dissims_best);
      
  #if 0
    if (name == "111838")
    {
      location. saveText (cout);
      cout << endl;
    }
  #endif
  }
}



void NewLeaf::anchor2location ()
{
  ASSERT (location. anchorLen > 0.0);

  Real mult_sum = 0.0;
  Real u_avg = 0.0; 
  Real delta_avg = 0.0;
  Real delta_u_sum = 0.0;
  for (const Leaf2dissim& ld : leaf2dissims)
  {
  //ASSERT (ld. dissim > 0);
    if (ld. mult)
    {
      mult_sum    += ld. mult;
      u_avg       += ld. mult * ld. getU ();
      delta_avg   += ld. mult * ld. getDelta ();
      delta_u_sum += ld. mult * ld. getDelta () * ld. getU ();
    }
  } 
  ASSERT (mult_sum > 0.0);  // <= pre-condition
  u_avg     /= mult_sum;
  delta_avg /= mult_sum;
  
  Real u_var = 0.0;
  Real delta_u_cov = 0.0;
  for (const Leaf2dissim& ld : leaf2dissims)
    if (ld. mult)
    {
      delta_u_cov += ld. mult * (ld. getDelta () - delta_avg) * (ld. getU () - u_avg);
      u_var       += ld. mult * sqr (ld. getU () - u_avg);
    }
  u_var       /= mult_sum;
  delta_u_cov /= mult_sum;  
  ASSERT (u_var >= 0.0);  
//IMPLY (u_var > 0.0, mult_sum > 0.0);


  IMPLY (location. anchor == tree. root, ! u_var);
  

  // location.{arcLen,leafLen}
  bool bad = false;
  if (u_var)  
  {
    location. arcLen = delta_u_cov / u_var;
    if (   location. arcLen >= 0.0 
        && location. arcLen <= location. anchorLen
       )
    {
      location. leafLen = delta_avg - u_avg * location. arcLen;
      if (location. leafLen < 0.0)
        bad = true;
    }
    else
      bad = true;
    if (bad)
    {
      location. leafLen = 0.0;
      location. arcLen = min (max (0.0, delta_u_sum / mult_sum), location. anchorLen);
      location. setAbsCriterion_leaf (*this); 
      const Location loc1 (location);
      //
      location. leafLen = max (0.0, delta_avg);
      location. arcLen = 0.0;
      location. setAbsCriterion_leaf (*this); 
      const Location loc2 (location);
      //
      location. leafLen = max (0.0, delta_avg - location. anchorLen * u_avg);
      location. arcLen = location. anchorLen;
      location. setAbsCriterion_leaf (*this); 
      const Location loc3 (location);
      //
      if (   loc1. absCriterion_leaf < loc2. absCriterion_leaf
          && loc1. absCriterion_leaf < loc3. absCriterion_leaf
         )
        location = loc1;
      else
        if (loc2. absCriterion_leaf < loc3. absCriterion_leaf)
          location = loc2;
        else
          location = loc3;
    }
  }
  else
  {
    location. leafLen = max (0.0, delta_avg);
    location. arcLen = 0.0;
  //ASSERT (location. anchor);
  //location. anchorLen = location. anchor->len;
  }
  if (! bad)
    location. setAbsCriterion_leaf (*this); 

  location. qc ();
}



bool NewLeaf::descend (const DTNode* anchorChild)
{
  ASSERT (anchorChild);
  ASSERT (anchorChild->getParent () == location. anchor);
  ASSERT (location. anchor);
  
  
  VectorPtr<Tree::TreeNode> nodeVec;  nodeVec. reserve (leaf2dissims. size ());
#ifndef NDEBUG
	size_t below = 0;
#endif
  for (Leaf2dissim& ld : leaf2dissims)
  {
    ASSERT (ld. leaf);
    if (! ld. leafIsBelow)
      continue;
  #ifndef NDEBUG
    below++;
  #endif
    if (ld. leaf->descendantOf (anchorChild))
      nodeVec << ld. leaf;
    else
      ld. leafIsBelow = false;
  }
  if (nodeVec. empty ())
    return false;
    
  IMPLY (location. anchor->getParent (), nodeVec. size () < below);  // <= location.anchor is the LCA of leaf2dissim::leaf's which are below location.anchor
    
	Tree::LcaBuffer buf;
  const Tree::TreeNode* lca_ = Tree::getLca (nodeVec, buf);  
  ASSERT (lca_);
  
  const DTNode* lca = static_cast <const DTNode*> (lca_);
  lca = lca->getDiscernible ();
  ASSERT (lca);
  ASSERT (lca != location. anchor);
  ASSERT (lca->descendantOf (anchorChild));
  
  IMPLY (node_orig, ! lca->descendantOf (node_orig));  // <= node_orig has no neighbors under itself
  
  location. anchorLen = lca->getPathLength (location. anchor);
  location. anchor = lca;
  ASSERT (location. anchorLen >= location. anchor->len);
  
  for (Leaf2dissim& ld : leaf2dissims)
    if (ld. mult)
      ld. dist_hat -= location. anchorLen * ld. getU ();
  
  return true;
}



void NewLeaf::qc () const
{
  if (! qc_on)
    return;
  Named::qc ();
    
  location. qc ();    
  QC_ASSERT (& location. anchor->getDistTree () == & tree);
  
  for (const Leaf2dissim& ld : leaf2dissims)
  {
    ld. qc ();
    ASSERT (ld. leaf);
    QC_ASSERT (& ld. leaf ->getDistTree () == & tree);
  }
}



void NewLeaf::saveResult (const string &fName,
                          size_t closest_num) const
{ 
  if (fName. empty ())
    return;
  if (! closest_num)
    return;
    
  OFStream f (fName);  

  Vector<Tree::TreeNode::NodeDist> neighbors;  neighbors. reserve (closest_num + 1);
  location. anchor->getClosestLeaves (closest_num, neighbors);
  if (const Leaf* leaf = static_cast <const DTNode*> (location. anchor) -> asLeaf ())
    neighbors << Tree::TreeNode::NodeDist {leaf, 0.0};
  for (Tree::TreeNode::NodeDist& nd : neighbors)
  {
    ASSERT (nd. node);
    nd. dist += location. leafLen + location. arcLen;
  }
  if (const DTNode* parent = static_cast <const DTNode*> (location. anchor->getParent ()))
  {
    Vector<Tree::TreeNode::NodeDist> neighbors2;  neighbors2. reserve (closest_num + 1);
    parent->getClosestLeaves (closest_num, neighbors2);
    location. anchor->getClosestLeaves (closest_num, neighbors2);
    for (Tree::TreeNode::NodeDist& nd : neighbors2)
    {
      ASSERT (nd. node);
      ASSERT (location. anchorLen >= location. arcLen);
      nd. dist += location. leafLen + (location. anchorLen - location. arcLen);
    }
    neighbors << neighbors2;
  }
  neighbors. sort (Tree::TreeNode::NodeDist::distLess);
  Set<const Tree::TreeNode*> nodes;
  for (const Tree::TreeNode::NodeDist& nd : neighbors)
  {
    ASSERT (nd. node);
    if (nodes. contains (nd. node))
      continue;
    f << nd. node->getName () << '\t' << nd. dist << '\n';
    nodes << nd. node;
    if (nodes. size () == closest_num)
      break;
  }
}



}



/* TO DO ??

non-stability of results

*/
