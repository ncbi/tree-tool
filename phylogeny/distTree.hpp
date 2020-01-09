// distTree.hpp

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


#ifndef DISTTREE_HPP
#define DISTTREE_HPP

#include "../common.hpp"
#include "../graph.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
#include "../dm/dataset.hpp"
using namespace DM_sp;



namespace DistTree_sp
{


extern Chronometer chron_getBestChange;
extern Chronometer chron_tree2subgraph;
extern Chronometer chron_subgraphOptimize;
extern Chronometer chron_subgraph2tree;



// PAR
constexpr streamsize dissimDecimals = 6;  
constexpr streamsize absCriterionDecimals = 4;  // Small for stability
constexpr streamsize relCriterionDecimals = 3;
constexpr uint areaRadius_std = 5;  
constexpr uint areaDiameter_std = 2 * areaRadius_std; 
constexpr uint subgraphDepth = areaRadius_std;  
constexpr uint boundary_size_max_std = 500;  
constexpr uint sparsingDepth = areaDiameter_std;  // must be: >= areaRadius_std
constexpr Prob rareProb = 0.01; 
constexpr size_t dissim_progress = 100000;
constexpr Real dissimCoeffProd_delta = 1e-6; 



// For Time: 
//   n = # Tree leaves
//   p = # distances = DistTree::dissims.size()
//   p >= n
//   O(): log log = 1



// --> DistTree ??
// Dissimilarity variance
enum VarianceType { varianceType_lin     // Dissimilarity ~ Poisson
                  , varianceType_sqr   
                  , varianceType_pow  
                  , varianceType_exp     // Dissimilarity = -ln(P), var P = const
                  , varianceType_linExp  // Dissimilarity = -ln(P), var P = p*(1-p)
                  , varianceType_none
                  };
extern const StringVector varianceTypeNames;
extern VarianceType varianceType;
extern Real variancePower;
extern Real variance_min;

extern Real dissim_power;
extern Real dissim_coeff;  // Irrelevant if varianceType = varianceType_lin
extern Real hybridness_min;
extern Real dissim_boundary;


inline VarianceType str2varianceType (const string &s)
  { size_t index = 0;
    if (varianceTypeNames. find (s, index))
      return (VarianceType) index;
    throw logic_error ("Unknown dissimilarity variance " + s);
  }      

// Input: varianceType
inline Real dist2mult (Real dist)
  { if (dist < 0.0)
      throw runtime_error ("Negative dist");
    Real var = NaN;  // Variance function
    switch (varianceType)
    { case varianceType_lin:    var = dist; break;
      case varianceType_sqr:    var = sqr (dist); break;
      case varianceType_pow:    var = pow (dist, variancePower); break;
      case varianceType_exp:    var = exp (2.0 * dist); break;
      case varianceType_linExp: var = exp (dist) - 1.0; break;
      case varianceType_none:   throw runtime_error ("Variance function is not specified");
      default:                  throw logic_error ("Unknown variance function");
    }  
    return 1.0 / (variance_min + var);
  }
  // Return: >= 0
  //         0 <=> dist = INF

inline Real dist_max ()
  { Real x = NaN;
    switch (varianceType)
    { case varianceType_lin:    x = 1.0 / epsilon; break;
      case varianceType_sqr:    x = 1.0 / sqrt (epsilon); break;
      case varianceType_pow:    x = pow (epsilon, - 1.0 / variancePower); break;
      case varianceType_exp:    x = - 0.5 * log (epsilon); break;
      case varianceType_linExp: x = log (1.0 / epsilon + 1.0); break;
      case varianceType_none:   x = INF; break;
      default:                  throw logic_error ("Unknown variance function");
    }  
    return pow (x / dissim_coeff, 1.0 / dissim_power);
  }
  // Solution of: dist2mult(dist_max) = epsilon
  // dist < dist_max() <=> !nullReal(dist2mult(dist))


void dissimTransform (Real &target);
  // Update: target


inline bool at_dissim_boundary (Real dissim)
  { return     dissim_boundary >= dissim 
  	       && (dissim_boundary -  dissim) / dissim_boundary <= 0.05;  // PAR 
  }
  // Has a small probability <= choice of dissim_boundary


constexpr static uint dissims_max {numeric_limits<uint>::max ()};



struct DistTree;
struct Image;

//struct DTNode;
  struct Steiner;
  struct Leaf;

struct NewLeaf;

struct DissimLine;



struct Neighbor 
{ 
	const Leaf* leaf {nullptr}; 
	  // !nullptr
	size_t dissimType {NO_INDEX};
	Real target {NaN};
	  // > 0

	  
	Neighbor (const Leaf* leaf_arg,
	          size_t dissimType_arg,
	          Real target_arg);
  Neighbor (const Leaf* leaf_arg,
            size_t dissimType_arg);

	  
	bool operator< (const Neighbor &other) const;
	bool operator== (const Neighbor &other) const
	  { return    leaf       == other. leaf
	           && dissimType == other. dissimType; 
	  }
};



struct Triangle  
// Triple of Leaf's with a triangle inequality violation
{ 
	struct Parent
	{ 
	  const Leaf* leaf {nullptr};
		Real dissim {NaN};
		  // = d(Triangle::child,leaf)
		  // > 0
		bool hybrid {false};
			// Cause of the triangle inequality violation
	};

	// !nullptr
	const Leaf* child {nullptr};
	Real hybridness {NaN};
	  // = d(parent1,parent2) / (d(child,parent1) + d(child,parent2))
	  // > 1
	array<Parent, 2> parents;	
	bool child_hybrid {false};
		// Cause of the triangle inequality violation
	size_t dissimType {NO_INDEX};

	  
	Triangle (const Leaf* child_arg,
		        Real parentDissim_arg,
				  	const Leaf* parent1,
				  	const Leaf* parent2,
				  	Real parent1_dissim,
				  	Real parent2_dissim,
				  	size_t dissimType_arg);
	Triangle () = default;
	void qc () const;
	void print (ostream &os) const;
	  // Matches PositiveAttr2::hybrid_format
	
	
	bool operator== (const Triangle &other) const
    { return    child             == other. child
             && parents [0]. leaf == other. parents [0]. leaf
             && parents [1]. leaf == other. parents [1]. leaf
             && dissimType        == other. dissimType;
    }
	bool operator< (const Triangle &other) const;
	Real parentsDissim () const
	  { return (parents [0]. dissim + parents [1]. dissim) * hybridness; }
	  // Return: d(parent1,parent2)
	Prob parent_dissim_ratio () const
	  { return min ( parents [0]. dissim / parents [1]. dissim
                 , parents [1]. dissim / parents [0]. dissim
                 );
    }
	bool hasHybrid () const
	  { return    child_hybrid 
	  	       || parents [0]. hybrid
	  	       || parents [1]. hybrid;
	  }
	  // Return: false <= undecided
	VectorPtr<Leaf> getHybrids (bool hybrid) const
	  { VectorPtr<Leaf> vec;  vec. reserve (3);
	  	if (child_hybrid        == hybrid)  vec << child;
	  	if (parents [0]. hybrid == hybrid)  vec << parents [0]. leaf;
	  	if (parents [1]. hybrid == hybrid)  vec << parents [1]. leaf;
	  	return vec;
	  }
	void qcMatchHybrids (const VectorPtr<Leaf> &hybrids) const;
	  // Requires: hybrids: sort()'ed
};



void addTriangle (Vector<Triangle> &triangles,
                  const Leaf* leaf1,
                  const Leaf* leaf2,
                  const Leaf* leaf3,
                  Real target23,
                  Real target13,
                  Real target12,
                  size_t dissimType);
  // Update: triangles: append by <= 1 element



struct TriangleParentPair
{
	// Input
	struct Parent
	{ 
	  const Leaf* leaf {nullptr};
			// !nullptr
		size_t classSize {0};
		  // Output
	};
	array<Parent,2> parents;	
	Real parentsDissim {NaN};
	  // = f(parents[0].leaf,parents[1].leaf)
	size_t dissimType {NO_INDEX};
	
	// Output
	Vector<Triangle> triangles;
	  // Triangle::parents[i].leaf = parents[i].leaf
	  // Clusterize Triangle::child's ??
	  // May be empty()
    // Average size: O(p/n log(n))  
private:
	size_t triangle_best_index {NO_INDEX};
	  // Index in triangles
public:
	Real hybridness_ave {NaN};

	
	TriangleParentPair (const Leaf* parent1,
		                  const Leaf* parent2,
		                  Real parentsDissim_arg,
		                  size_t dissimType_arg)
		: parentsDissim (parentsDissim_arg)
		, dissimType (dissimType_arg)
		{ parents [0]. leaf = parent1;
			parents [1]. leaf = parent2;
	  }
	TriangleParentPair () = default;
	void setTriangles (const DistTree &tree);
	  // Output: triangles
	  // Average time: O(p/n log(n))
  void triangles2hybridness_ave ();
	  // Output: hybridness_ave
  void finish (const DistTree &tree,
               const Set<const Leaf*> &hybrids);
    // Input: triangles, Leaf::badCriterion
    // Output: Triangle::*hybrid
    // Invokes: child_parent2parents()
	void qc () const;
  void print (ostream &os) const;
  static constexpr const char* format {"<child> <parent1> <parent2> <# children> <# parents 1> <# parents 2> <hybridness> <d(child,parent1)> <d(child,parent2)> <child is hybrid> <parent1 is hybrid> <parent2 is hybrid> <dissimilarity type>"};


  bool operator== (const TriangleParentPair &other) const
    { return    parents [0]. leaf == other. parents [0]. leaf
             && parents [1]. leaf == other. parents [1]. leaf
             && dissimType        == other. dissimType;
    }
  bool operator< (const TriangleParentPair &other) const;
  static bool compareHybridness (const TriangleParentPair &hpp1,
                                 const TriangleParentPair &hpp2);
  const Triangle& getBest () const
    { if (triangle_best_index < triangles. size ())
    	  return triangles [triangle_best_index]; 
    	throw logic_error ("TriangleParentPair::getBest()");
    }
  bool dissimError () const  
    { if (   getBest (). parent_dissim_ratio () < 0.25  // PAR
    	    && hybridness_ave < 1.25  // PAR
    	   )
    	  return true;
			for (const bool i : {false, true})
			  if (at_dissim_boundary (getBest (). parents [i]. dissim))
				  return true; 
      return false;
    }
  Vector<Triangle> getHybridTriangles () const
    { Vector<Triangle> vec;  vec. reserve (triangles. size ());
    	for (const Triangle& tr : triangles)
    		if (tr. hasHybrid ())
    			vec << tr;
    	return vec;
    }
  VectorPtr<Leaf> getHybrids (bool hybrid) const
    { VectorPtr<Leaf> vec;  vec. reserve (triangles. size ());
    	for (const Triangle& tr : triangles)
    		if (tr. hasHybrid ())
    			vec << move (tr. getHybrids (hybrid));
    	return vec;
    }
  void qcMatchHybrids (const VectorPtr<Leaf> &hybrids) const
    { if (! qc_on)
        return;
      for (const Triangle& tr : triangles)
    		if (tr. hasHybrid ())
	    		tr. qcMatchHybrids (hybrids);
    }
	  // Requires: hybrids: sort()'ed
  bool undecided () const
    { return    ! triangles. empty ()
             && ! getBest (). hasHybrid ();
    }
private:
	size_t child_parent2parents (const DistTree &tree,
                               const Leaf* child,
                               const Leaf* parent,
                               Real parentDissim) const;
	  // Average time: O(p/n log(n))
	void setChildrenHybrid ()
	  { for (Triangle& tr : triangles)
	  	  tr. child_hybrid = true;
	  }
};



typedef  unordered_map <const DisjointCluster*, VectorPtr<Leaf>>  LeafCluster;



struct DTNode : Tree::TreeNode 
{
  friend DistTree;
  friend Image;
  friend Steiner;
  friend Leaf;
  friend NewLeaf;

  string name;  
    // !empty() => from Newick
	Real len;
	  // Arc length between *this and *getParent()
	  // *this is root => NaN
  Vector<uint/*dissimNum*/> pathDissimNums; 
    // Paths: function of getDistTree().dissims
    // Dissimilarity paths passing through *this arc
    // asLeaf() => getDistTree().dissims[dissimNum].hasLeaf(this)  
    //             aggregate size = 2 p
    // !asLeaf(): aggregate size = O(p log(n))
    // Distribution of size ??
    // Max. size() is at about the topological center
private:
  bool stable {false};
    // Init: false
public:

  // Temporary
private:
  const DTNode* tarjanLca {nullptr};
public:
protected:
  WeightedMeanVar subtreeLen; 
    // Average subtree height 
    // weights = topological ? # leaves : sum of DTNode::len in the subtree excluding *this
public:
  Real errorDensity {NaN};
    // ~ Normal(0,1)
  uint maxDeformationDissimNum {dissims_max};
    // Index of DistTree::dissims[]


protected:
	DTNode (DistTree &tree,
          Steiner* parent_arg,
	        Real len_arg);
public:
  void qc () const override;
  void saveContent (ostream& os) const override;


  virtual const Steiner* asSteiner () const
    { return nullptr; }
  virtual const Leaf* asLeaf () const
    { return nullptr; }


	double getParentDistance () const final
	  { return isNan (len) ? 0.0 : len; }

  const DistTree& getDistTree () const;

  const Leaf* inDiscernible () const;
    // Return: this or nullptr
  bool childrenDiscernible () const
    { return arcs [false]. empty () || ! static_cast <DTNode*> ((*arcs [false]. begin ()) -> node [false]) -> inDiscernible (); }
  const DTNode* getDiscernible () const;
    // Return: this or getParent(); !nullptr
  Real getHeight_ave () const
    { return subtreeLen. getMean (); }    
    // After: DistTree::setHeight()
  Real getDeformation () const;
    // Input: maxDeformationDissimNum
  virtual const Leaf* getReprLeaf (ulong seed) const = 0;
    // Return: !nullptr, in subtree
    // For sparse *getDistTree().dissimAttr
    // Deterministic <=> (bool)seed
    // Invokes: getDistTree().rand
    // Time: O(log(n))
private:
  void saveFeatureTree (ostream &os,
                        size_t offset) const;
  virtual void setSubtreeLenUp (bool topological) = 0;
    // Output: subtreeLen
  void setGlobalLenDown (bool topological,
                         DTNode* &bestDTNode,
                         Real &bestDTNodeLen_new,
                         WeightedMeanVar &bestGlobalLen);
    // Output: subtreeLen: Global len = average path length from *this to all leaves
  virtual void getDescendants (VectorPtr<DTNode> &descendants,
                               size_t depth,
                               const DTNode* exclude) const = 0;
    // Update: descendants (append)
  virtual void setLca ();
    // Off-line LCA algorithm by Tarjan
    // Requires: !tarjanLca 
  Vector<uint/*dissimNum*/> getLcaDissimNums ();
    // Return: dissimNum's s.t. getDistTree().dissims[dissimNum].lca = this
    // Invokes: DTNode::pathDissimNums.sort()
  void setErrorDensity (Real absCriterion_ave);
    // Output: errorDensity
    // Time: O(|pathDissimNums|)
  VectorPtr<Leaf> getSparseLeafMatches (const string &targetName,
                                        size_t depth_max,
                                        bool subtractDissims,
                                        bool refreshDissims) const;
    // Return: size = O(log(n)); sort()'ed, uniq()'ed
    //         getDistTree().reroot(true) reduces size()
    // Input: targetName: for Rand::setSeed()
    //        depth_max: 0 <=> no restriction
    //        refreshDissims => improves criterion and quality; number of new dissims = ~10% of dissims
    // Time: O(log^2(n)) 

  struct ClosestLeaf
  {
    const Leaf* leaf;
    Real dist;
    void qc () const;
  };
  virtual Vector<ClosestLeaf> findGenogroups (Real genogroup_dist_max) = 0;
    // Return: dist <= genogroup_dist_max
};



struct Steiner : DTNode
// Steiner node
{
	Steiner (DistTree &tree,
	         Steiner* parent_arg,
	         Real len_arg);
	void qc () const override;
  void saveContent (ostream& os) const final;


  const Steiner* asSteiner () const final
    { return this; }

  bool isInteriorType () const final
    { return childrenDiscernible (); }

private:
  const Leaf* getReprLeaf (ulong seed) const final;
  void setSubtreeLenUp (bool topological) final;
  void getDescendants (VectorPtr<DTNode> &descendants,
                       size_t depth,
                       const DTNode* exclude) const final;
  void setLca () final;

  void reverseParent (const Steiner* target, 
                      Steiner* child);
    // Until target
    // Input: target: !nullptr
    //        child: nullptr <=> *this becomes getTree().root
    // Requires: descendantOf(target)
    // Invokes: setParent(child)
public:
  void makeRoot (Steiner* ancestor2descendant);
    // Opposite: ancestor2descendant->makeRoot(this);
    // Invokes: setParent(ancestor2descendant->getParent()); contents = ancestor2descendant->contents
  const Steiner* makeDTRoot ();
    // Return: Old root, !nullptr
    // Invokes: makeRoot(getTree().root)
  LeafCluster getIndiscernibles ();
    // Invokes: Leaf->DisjointCluster
  Vector<ClosestLeaf> findGenogroups (Real genogroup_dist_max) final;
    // Time: O(n log(n))+

#if 0
private:
	friend DistTree;
  void setSubTreeWeight ();
    // Input: lcaNum
    // Update: subTreeWeight
  // Update: threadNum
  void threadNum2subTree (size_t threadNum_arg);
  void threadNum2ancestors (size_t threadNum_arg);
#endif
};



struct Leaf : DTNode
// name: !empty()
{
	friend DistTree;
	
  string comment;
  static const string non_discernible;
  bool discernible {true};  // May be not used: parameter ??
    // false => getParent()->getChildren() is an equivalence class of indiscernibles
  Real normCriterion {NaN};
    // ~ Normal(0,1)
    
  // Temporary
private:
  size_t index {NO_INDEX};
public:
  // For DistTree::findHybrids()
  Real badCriterion {NaN};
  

	Leaf (DistTree &tree,
	      Steiner* parent_arg,
	      Real len_arg,
	      const string &name_arg);
	void qc () const final;
  void saveContent (ostream& os) const final;


  const Leaf* asLeaf () const final
    { return this; }


  string getName () const final
    { return name; }
  string getNewickName (bool minimal) const final
    { if (minimal)
        return name;
      string s = name + prepend (" ", comment); 
      if (! isNan (normCriterion))
        s += " " + real2str (normCriterion, 1, false);  // PAR
      return s;
    }
  bool isLeafType () const final
    { return true; }

private:
  const Leaf* getReprLeaf (ulong /*seed*/) const final
    { return this; }
  void setSubtreeLenUp (bool topological) final
    { subtreeLen. clear ();
      if (topological)
    	  subtreeLen. add (0.0, 1.0);
    }
  void getDescendants (VectorPtr<DTNode> &descendants,
                       size_t /*depth*/,
                       const DTNode* exclude) const final
    { if (this != exclude)
    	  descendants << this; 
    }
  void setLca () final;
public:

  const Leaf* getDissimOther (size_t dissimNum) const;
    // Return: !nullptr; != this
  bool getCollapsed (const Leaf* other) const
    { return    other
             && getParent () == other->getParent ()
             && ! discernible
             && ! other->discernible;
    }
private:
  friend DissimLine;
  void collapse (Leaf* other);
    // Output: discernible = false
    // Requires: !DistTree_sp::variance_min
    // Invokes: setParent()
    // To be followed by: DistTree::cleanTopology()
public:	
  void addHybridTriangles (Vector<Triangle> &triangles) const;
    // Invokes: addTriangle()
    // Average time: O(p^2/n^2 log^2(n))  
  Vector<ClosestLeaf> findGenogroups (Real genogroup_dist_max) final
    { if (len <= genogroup_dist_max)
      { Vector<ClosestLeaf> res {{this, len}};
        return res; 
      }
      return {};
    }
};



struct SubPath
// Path going through a connected subgraph
{
  uint dissimNum {dissims_max};    
    // Index of DistTree::dissims
  const DTNode* node1 {nullptr};
  const DTNode* node2 {nullptr};
    // !nullptr, different
  Real dist_hat_tails {NaN};

    
  SubPath () = default;
  explicit SubPath (uint dissimNum_arg)
    : dissimNum (dissimNum_arg)
    {}
  void qc () const;

  
  bool contains (const DTNode* node) const
    { return    node1 == node
             || node2 == node;
    }
};



struct Subgraph : Root
{
  const DistTree& tree;
  // !nullptr
  VectorPtr<Tree::TreeNode> area;  
    // Connected area
    // sort()'ed
  VectorPtr<Tree::TreeNode> boundary;
    // Of area
    // Size: O(|area|)
    // sort()'ed
  const Steiner* area_root {nullptr};
    // May be nullptr
    // boundary.contains(area_root)
private:  
  const DTNode* area_underRoot {nullptr};
    // May be nullptr
    // area.contains(area_underRoot)
  // (bool)area_underRoot = (bool)area_root
  Vector<uint> dissimNums;
    // tree.dissims passing through area which can be changed
    // Size: O(|bounadry| p/n log(n))
  bool completeBoundary {false};
public:  
  Vector<SubPath> subPaths;
    // Size: O(|bounadry| p/n log(n))
  Real subPathsAbsCriterion {0.0};
  
  
  explicit Subgraph (const DistTree &tree_arg);
  void qc () const override;
  bool empty () const override
    { return    area. empty ()
             && boundary. empty ()
             && ! area_root
             && ! area_underRoot
             && dissimNums. empty ()
             && ! completeBoundary
             && subPaths. empty () 
             && ! subPathsAbsCriterion;
    }
  void clear () override
    { area. wipe ();
      boundary. wipe ();
      area_root = nullptr;
      area_underRoot = nullptr;
      dissimNums. wipe ();
      completeBoundary = false;
      subPaths. wipe ();
      subPathsAbsCriterion = 0.0;
    }

  
  // Usage:
//set area, boundary
  void reserve (uint radius);
  void removeIndiscernibles ();
    // Update: area, boundary
  void finish ();
    // Output: area_root, area_underRoot
    // Time: O(|area| log(|area|))
//set dissims
  void dissimNums2subPaths ();
    // Output: subPaths, subPathsAbsCriterion
    // Time: O(|dissimNums| log(|area|)) 
//change topology of tree within area
#if 0
  // Not used
  Real getImprovement (const DiGraph::Node2Node &boundary2new) const;
    // Time: O(|subPaths| (log(|boundary|) + log(|area|)))
#endif
  void subPaths2tree ();
    // Update: tree: Paths, absCriterion, Dissim::prediction
    // Time: O(|subPaths| + |area| (log(|boundary| + p/n log(n)) + |subPaths| log(|area|)

  bool large () const
    { return boundary. size () > 64; } // PAR
  bool unresolved () const
    { const Real resolution = (Real) area. size () / (Real) boundary. size ();
        // 1..2; 2 <=> completely resolved
      return resolution <= 1.2;   // PAR 
    }
  bool viaRoot (const SubPath &subPath) const
    { return    subPath. node1 == area_root 
             || subPath. node2 == area_root;
    }
  void node2dissimNums (const DTNode* node);
    // Time: O(p/n log(n))
  void area2dissimNums ();
    // Output: dissimNums, completeBoundary
    // Invokes: node2subPaths()
    // Time: O(|boundary| p/n log(n))
  VectorPtr<Tree::TreeNode>& getPath (const SubPath &subPath,
  	                                  Tree::LcaBuffer &buf) const
    // Return: reference to buf
    { const Tree::TreeNode* lca_ = nullptr;
      // tree.dissims[subPath.dissimNum].lca can be used instead of area_root if viaRoot(subPath) and tree topology has not been changed ??
      return Tree::getPath (subPath. node1, subPath. node2, area_root, lca_, buf);
    }
    // Requires: subPath in subPaths
  const Vector<uint>& boundary2pathDissimNums (const DTNode* dtNode) const
    { return dtNode == area_root /* && area_underRoot */
               ? area_underRoot->pathDissimNums
               : dtNode        ->pathDissimNums;
    }
  const Leaf* getReprLeaf (const DTNode* dtNode) const
    { return dtNode == area_root 
               ? static_cast <const DTNode*> (dtNode->getDifferentChild (area_underRoot)) -> getReprLeaf (0)
               : dtNode->getReprLeaf (0);
    }
};



struct Change : Root
// Of topology
// *to becomes a sibling of *from
// Enough to transform any topology to any topology. Proof: by induction by node depth descending
{
private:
	const DistTree& tree;
public:
	const DTNode* from;
	  // !nullptr
	const DTNode* to;
	  // !nullptr
	// Output of apply_()
	VectorPtr<DTNode> targets;  
	  // DTNode's whose len may be changed 
	Real improvement {NaN};
	  // isNan() or positive()
    // Too small values are noise => not stable in tree sampling
private:
	Real fromLen {NaN};
	Real toLen {NaN};
	// !nullptr
	Steiner* oldParent {nullptr};
	  // Old from->getParent()
	Steiner* arcEnd {nullptr};
	  // Old to->getParent()
	Steiner* inter {nullptr};
	  // Between *to and *arcEnd
  Subgraph subgraph;
  enum Status {eInit, eApplied, eDone};
  Status status {eInit};
public:

	
	Change (const DTNode* from_arg,
				  const DTNode* to_arg)
		: tree (var_cast (from_arg->getDistTree ()))
		, from (from_arg)
		, to (to_arg)
		, targets {from, to}
		, subgraph (tree)
		{}
    // Requires: valid()
	static bool valid (const DTNode* from_arg,
	                   const DTNode* to_arg)
	  { return    from_arg
             && from_arg->graph
             && ! from_arg->inDiscernible ()
	           && to_arg
	  	       && to_arg->graph == from_arg->graph
	  	       && to_arg != from_arg
             && ! to_arg->inDiscernible ()
    	       && from_arg->getParent ()
	  	       && ! to_arg->descendantOf (from_arg)
	  	       && ! (from_arg->getParent() == to_arg->getParent() && from_arg->getParent() -> arcs [false]. size () == 2)  
	  	       && ! (from_arg->getParent() == to_arg              && from_arg->getParent() -> arcs [false]. size () == 2); 
	  }
 ~Change ();
	void qc () const override;
	  // Invokes: valid()
protected:
	void saveText (ostream& os) const override
	  { os << from->getName () << " (parent = " << (from->getParent () ? from->getParent () -> getName () : "null") << ") -> " << to->getName () 
         << "  improvement = " << improvement; 
	  }
public:
	void print (ostream& os) const override
	  { saveText (os);
	    os << endl; 
	  }


  bool valid () const
    { return valid (from, to); }
  // Update: tree topology, DTNode::len, tree.dissims[].prediction
	bool apply ();
	  // Return: success
	  // Minimum change to compute tree.absCriterion
	  // status: eInit --> eApplied|eFail
	  // Time: O(log^4(n))
	void restore ();
	  // Output: tree.dissims[].prediction
	  // status: eApplied --> eInit
	void commit ();
	  // status: eApplied --> eDone
    // May invoke: tree.delayDeleteRetainArcs()
    // Time: O(log^2(n))
	static bool strictlyBetter (const Change* a, 
	                           const Change* b);
    // Requires: (bool)a
};



struct DissimType : Named
{
  const PositiveAttr2* dissimAttr {nullptr};
    // In *DistTree::dissimDs
  Real scaleCoeff {NaN}; 
    // >= 0, < INF
    
  explicit DissimType (const PositiveAttr2* dissimAttr_arg);
  void qc () const override;
  void saveText (ostream &os) const override
    { const ONumber on (os, dissimDecimals, true);
      os << name << ' ' << scaleCoeff << endl; 
    }
}; 



struct Dissim
{
	// Input
  // !nullptr
  // leaf1->name < leaf2->name
  const Leaf* leaf1 {nullptr};
  const Leaf* leaf2 {nullptr};
  //
  // Use 4-byte variables ??
  Real target {NaN};
    // Dissimilarity between leaf1 and leaf2; !isNan()
    // < INF
    // Update: = original target * DissimType::scaleCoeff
  size_t type {NO_INDEX};
    // < DistTree::dissimTypes.size()
  
  // Output
  Real prediction {NaN};
    // Tree distance
    // >= 0
  Real mult {NaN};
    // >= 0
    // INF <=> leaf1 and leaf2 must be collapse()'ed
  const Steiner* lca {nullptr};
    // Paths
  

  Dissim (const Leaf* leaf1_arg,
          const Leaf* leaf2_arg,
          Real target_arg,
          Real mult_arg,
          size_t type_arg);
  Dissim () = default;
  void saveText (ostream &os) const
    { os <<         leaf1->name 
         << '\t' << leaf2->name 
         << '\t' << target 
         << '\t' << type
         << '\t' << prediction
         << '\t' << mult
         << endl;
    }
  void qc () const;

          
  bool valid () const
    { return    leaf1->graph
             && leaf2->graph;
    }
    // For topology
  bool validMult () const
    { return    valid ()
             && mult
             && mult < INF;
    }
  bool hasLeaf (const Leaf* leaf) const
    { return    leaf == leaf1
             || leaf == leaf2;
    }
  const Leaf* getOtherLeaf (const Leaf* leaf) const
    { if (leaf == leaf1) return leaf2;
    	if (leaf == leaf2) return leaf1;
    	throw logic_error ("getOtherLeaf");
    }
  bool indiscernible () const
    { return    ! leaf1->discernible
             && ! leaf2->discernible
             && leaf1->getParent () == leaf2->getParent ();
    }
  string getObjName () const;
  VectorPtr<Tree::TreeNode>& getPath (Tree::LcaBuffer &buf) const;
  	// Return: reference to buf
  Real getResidual () const
    { return prediction - target; }
  Real getAbsCriterion (Real prediction_arg) const;
  Real getAbsCriterion () const
    { return getAbsCriterion (prediction); }
  Real getDeformation () const
    { const Real residual = sqr (target - prediction);
      if (! residual)
        return 0.0;
      return residual / min (prediction, target);
    }
    // Return: distribution is Chi^2_1 if mean = 1
    
  void setPathDissimNums (size_t dissimNum,
                          Tree::LcaBuffer &buf);
    // Output: prediction, Steiner::pathDissimNums
  array<const Leaf*,2> getLeaves () const
    { array<const Leaf*, 2> leaves;
      leaves [0] = leaf1;
      leaves [1] = leaf2;
      return leaves;
    }                
    
  bool operator< (const Dissim &other) const;
  bool operator== (const Dissim &other) const
   { return    Pair<const Leaf*> (leaf1, leaf2) == Pair<const Leaf*> (other. leaf1, other. leaf2)
            && type == other. type; 
   }
};




struct Image : Nocopy
// Tree subgraph replica
{
  Subgraph subgraph;
  const DTNode* center {nullptr};
    // In subgraph.tree
    // May be delete'd
  DistTree* tree {nullptr};
    // = &subgraph.tree
    // nullptr <=> bad_alloc
  DiGraph::Node2Node new2old;  
    // Initially: newLeaves2boundary
  bool rootInArea {false};
  bool neighborJoinP {false};

  
  explicit Image (const DistTree &mainTree);
 ~Image ();


  void processSmall (const DTNode* center_arg,
			               uint areaRadius);
	  // Invokes: tree->{optimizeLenArc(),optimizeLenNode(),optimizeWholeIter() or optimizeSmallSubgraphs()}
	  // Time: ~ O(|area| (log(|area|) log^2(n) + |area|) + Time(optimizeWholeIter(|area|)))
	void processLarge (const Steiner* subTreeRoot,
	                   const VectorPtr<Tree::TreeNode> &possibleBoundary);
	  // Time: ~ O(|area| log(|area|) log^2(n) + Time(optimizeSmallSubgraphs(|area|)))
  bool apply ();
    // Return: false <=> bad_alloc
	  // Output: DTNode::stable = true
    // Time: ~ O(|area| log(|area|) log^2(n))
};




///////////////////////////////////////////////////////////////////////////

struct SteinerHash;



struct DistTree : Tree
// Of DTNode*
// Least-squares distance tree
// Steiner tree
// nodes.size() >= 2
{
  friend DTNode;
  friend Steiner;
  friend Leaf;
  friend Change;
  friend Subgraph;
  friend Image;
  friend DissimLine;

  const uint subDepth {0};
    // > 0 => *this is a subgraph of a tree with subDepth - 1
  typedef  unordered_map<string/*Leaf::name*/,const Leaf*>  Name2leaf;
  Name2leaf name2leaf;  // subDepth => replace by leavesSize ??
    // 1-1

private:
  // Temporary
  // Dissimilarity
  // May be nullptr
  unique_ptr<Dataset> dissimDs;
    // Original data
  const PositiveAttr2* dissimAttr {nullptr};
    // In *dissimDs
  const PositiveAttr2* multAttr {nullptr};
    // In *dissimDs   
public:
    
  Vector<Dissim> dissims;
  Vector<DissimType> dissimTypes;
    // Product(DissimType::scaleCoeff) = 1.0
  bool multFixed {false};
  Real mult_sum {NaN};
  Real target2_sum {NaN};
    // = sum_{dissim in dissims} dissim.target^2 * dissim.mult        
  Real absCriterion {NaN};
    // = L2LinearNumPrediction::absCriterion  

private:
  size_t leafNum {0};
    // For Leaf::index
	VectorOwn<DTNode> toDelete;
	VectorOwn<Leaf> detachedLeaves;
	  // !Leaf::graph
	mutable Rand rand;
public:
  
  struct DeformationPair
  {
    string leafName1;
    string leafName2;
    Real deformation;
  };
  unordered_map<const DTNode*,DeformationPair> node2deformationPair;
    // Requires: topology is unchanged

private:
  RandomSet<const Steiner*> unstableCut;
    // !Steiner::stable, but getParent()->stable
  size_t unstableProcessed {0};
    // For Progress
public:


  // Input: dissimFName: <dmSuff>-file without <dmSuf>, contains attributes dissimAttrName and multAttrName
  //                     may contain more objects than *this contains leaves
  //        dissimFName, dissimAttrName, multAttrName: all may be empty	  
	DistTree (const string &treeDirFName,
	          const string &dissimFName,
	          const string &dissimAttrName,
	          const string &multAttrName);
	  // Input: treeDirName: if directory anme then contains the result of mdsTree.sh; ends with '/'
	  // Invokes: loadTreeFile() or loadTreeDir(), loadDissimDs(), dissimDs2dissims(), setGlobalLen()
	DistTree (const string &dissimFName,
	          const string &dissimAttrName,
	          const string &multAttrName);
	  // Invokes: loadDissimDs(), dissimDs2dissims(), neighborJoin()
	DistTree (const string &dataDirName,
	          const string &treeFName,
            bool loadNewLeaves,
	          bool loadDissim,
	          bool optimizeP);
	  // Input: dataDirName: ends with '/': incremental distance tree directory:
	  //          temporary file name       line/file format                              meaning
	  //          -------------------       -----------------------------                 ----------------------------
    //          leaf                      <obj_new> <obj1>-<obj2> <leaf_len> <arc_len>
    //         [dissim.add[-req]]
    //          search/<obj_new>/                                                       Initialization of search for <obj_new> location
    //        [ search/<obj_new>/dissim   <obj_new> <obj> <dissimilarity>        
    //          search/<obj_new>/leaf     = as in leaf 
    //         [search/<obj_new>/request  <obj_new> <obj>]                              Request to compute dissimilarity
    //        ]
    //         [dissim.bad]               <obj1> <obj2> nan
	  //         [dissim_request]           <obj1> <obj2>                                 Request to compute dissimilarity
	  //       <dissimilarity>: >= 0, < INF
	  // Invokes: optimizeSmallSubgraph() for each added Leaf; Threads
	  // Time: if loadDissim then O(p log(n) + Time(optimizeSmallSubgraph) * new_leaves)
	  //       if !loadDissim then O(n log(n) + new_leaves)
  //  
  explicit DistTree (const string &newickFName);
  DistTree (Prob branchProb,
            size_t leafNum_max);
    // Random tree: DTNode::len = 1
    // Time: O(n)
  DistTree (Subgraph &subgraph,
            Node2Node &newLeaves2boundary,
            bool sparse);
    // Connected subgraph of subgraph.tree: boundary of subgraph.area are Leaf's of *this
    // If subgraph.unresolved() then the topology of *this is changed to a star
    // Input: subgraph: !empty(), not finish()'ed
    // Output: subgraph: area: contains newLeaves2boundary.values(); Leaf::discernible
    //         newLeaves2boundary
	  // Time: ~ O(|area| (log(|area|) log^2(subgraph.tree.n) + (sparse ? log(|area|) : |area|)))
  DistTree () = default;
private:
  void loadTreeDir (const string &dir);
	  // Input: dir: Directory with a tree of <dmSuff>-files
	  // Uses: temporary file "<dirFile>/.list"
	  // Invokes: getName2steiner()
  typedef  map<string,Steiner*>  Name2steiner;  
  Steiner* getName2steiner (const string &name,
                            Name2steiner &name2steiner);
    // Update: name2steiner
  void loadTreeFile (const string &fName);
    // InvokesL loadLines()
  bool loadLines (const StringVector &lines,
		              size_t &lineNum,
		              Steiner* parent,
		              size_t expectedOffset);
    // Return: a child of parent has been loaded
    // Output: topology, DTNode::len, Leaf::discernible
    // Update: lineNum
  void setName2leaf ();
  size_t getPathDissimNums_size () const
    { return 10 * (size_t) log (name2leaf. size () + 1); }  // PAR
  void loadDissimDs (const string &dissimFName,
                     const string &dissimAttrName,
                     const string &multAttrName);
    // Output: dissimDs
    // invokes: dissimDs->setName2objNum()
  // Input: dissimDs
  void mergeDissimAttrs ();
  bool getConnected ();
    // Find connected components of leaves where pairs have dissimilarities with positive multiplicity
    // Return: true <=> 1 connected component
    // Input: dissimDs
    // Output: DisjointCluster
    //         cout: print other connected components
  bool getDissimConnected ();
    // Find connected components of leaves where pairs have dissimilarities with positive multiplicity
    // Return: true <=> 1 connected component
    // Input: dissims
    // Output: DisjointCluster
  LeafCluster getIndiscernibles ();
    // Return: VectorPtr::size() >= 2
    // Invokes: Leaf->DisjointCluster
    // Time: ~ O(p)
  size_t leafCluster2discernibles (const LeafCluster &leafCluster);
    // Return: Number of indiscernible leaves
    // Output: Leaf::len = 0, Leaf::discernible = false, topology
    // Invokes: cleanTopology()
    // Time: O(n)
  size_t setDiscernibles_ds ();
    // Invokes: leafCluster2discernibles()
  void cleanTopology ();
    // Time: O(n)
  void setGlobalLen ();
    // Molecular clock 
    // Output: DTNode::len
    // Temporary: DTNode::subtreeLen
	  // Time: O(p log(n))
  void neighborJoin ();
    // Greedy
    // Assumes: Obj::mult = 1
    // Requires: isStar()
    // Invokes: reroot(true)
    // Time: O(n^3)
  //
  void dissimDs2dissims ();
    // Update: dissimDs: delete
    // Output: dissims etc.
    //         if an object is absent in dissimDs then it is deleted from the Tree
    // Invokes: getSelectedPairs(), setPaths()
  void loadDissimPrepare (size_t pairs_max);
    // Output: Dissim::target
  bool addDissim (Leaf* leaf1,
                  Leaf* leaf2,
                  Real target,
                  Real mult,
                  size_t type);
	  // Return: Dissim is added
    // Append: dissims[], Leaf::pathDissimNums
  bool addDissim (const string &name1,
                  const string &name2,
                  Real target,
                  Real mult,
                  size_t type)
    { return addDissim ( var_cast (findPtr (name2leaf, name1))
    	                 , var_cast (findPtr (name2leaf, name2))
    	                 , target
    	                 , mult
    	                 , type
    	                 );
    }
  void setPaths (bool setDissimMultP);
    // Output: dissims::Dissim, DTNode::pathDissimNums, absCriterion
    // Invokes: setLca(), setDissimMult()
    // Time: O(p log(n))
public:
	void qc () const override;
	  // Invokes: getIndiscernibles()


  size_t dissimTypesNum () const
    { return dissimTypes. empty () ? 1 : dissimTypes. size (); }
  void deleteLeaf (TreeNode* leaf,
                   bool deleteTransientAncestor) final;
    // Requires: !optimizable()
    
  static string getObjName (const string &name1,
                            const string &name2);
  const DTNode* lcaName2node (const string &lcaName,   
                              Tree::LcaBuffer &buf) const;
    // Return: !nullptr
    // Input: lcaName: <leaf1 name> <objNameSeparator> <leaf2 name>
  size_t getOneDissimSize_max () const
    { return name2leaf. size () * (name2leaf. size () - 1) / 2; }	
  size_t getDissimSize_max () const
    { return getOneDissimSize_max () * dissimTypesNum (); }	
  size_t getSparseDissims_size () const
    { return 7 * getPathDissimNums_size (); }  // PAR
  VectorPtr<DTNode> getDiscernibles () const;
    // Logical leaves
  static void printParam (ostream &os) 
    { os << "PARAMETERS:" << endl;
      os << "# Threads: " << threads_max << endl;
      os << "Variance function: " << varianceTypeNames [varianceType] << endl;
      if (! isNan (variancePower))
        os << "Variance power: " << variancePower << endl;
      if (DistTree_sp::variance_min)
        os << "Min. variance: " << variance_min << endl;
      os << "Max. possible distance: " << dist_max () << endl;
      if (dissim_power != 1.0)
        os << "Dissimilarity power: " << dissim_power << endl;
      if (dissim_coeff != 1.0)
        os << "Dissimilarity coefficient: " << dissim_coeff << endl;
    /*if (hybridness_min != 1.0)  // always > 1.0
        os << "Min. hybridness: " << hybridness_min << endl;
      if (! isNan (dissim_boundary))
        os << "Dissimilarity boundary (for hybrids): " << dissim_boundary << endl;*/
      os << "Subgraph radius: " << areaRadius_std << endl;
    }
	void printInput (ostream &os) const;
	bool optimizable () const  
	  { return ! dissims. empty (); }
	Real getDissim_ave () const
	  { WeightedMeanVar mv;
	    for (const Dissim& dissim : dissims)
	      if (dissim. validMult ())
	        mv. add (dissim. target, dissim. mult);
	    return mv. getMean ();
	  }
  Real getAbsCriterion_ave () const
    { return absCriterion / (Real) dissims. size (); }
    // Approximate: includes !Dissim::validMult() ?? 
  Prob getUnexplainedFrac (Real unoptimizable) const
    { return (absCriterion - unoptimizable) / (target2_sum - unoptimizable); }
  Real getRelCriterion (Real unoptimizable) const
    { return sqrt (getUnexplainedFrac (unoptimizable)); }
  string absCriterion2str (Real unoptimizable = 0.0) const
    { return real2str (absCriterion - unoptimizable, absCriterionDecimals); }
  void reportErrors (ostream &os,
                     Real unoptimizable = 0.0) const
    { const ONumber on (os, relCriterionDecimals, false);  
      os << "Abs. criterion = " << absCriterion2str (unoptimizable)
         << "  Rel. criterion = " << getRelCriterion (unoptimizable) * 100.0 << " %"
         << endl;
    }    
  void saveDissimCoeffs (const string &fName) const;
  void saveFeatureTree (const string &fName) const;

private:
  void qcPaths ();
    // Sort: DTNode::pathDissimNums 
    // Time: ~ O(p log(n))
  void setLca ();
    // Input: Leaf::pathDissimNums
    // Output: Dissim::lca
  void clearSubtreeLen ();
    // Invokes: DTNode::subtreeLen.clear()
  void setPredictionAbsCriterion ();
    // Output: Dissim::prediction, absCriterion
    // Invokes: Threads
    // Time: O(p log(n) / threads_max)
  void qcPredictionAbsCriterion () const;
public:
  size_t setDiscernibles ();
    // Invokes: getIndiscernibles(), leafCluster2discernibles()
  size_t fixTransients ();
    // Return: number of transient nodes deleted
  static Real path2prediction (const VectorPtr<TreeNode> &path);
    // Return: >= 0
	  // Input: DTNode::len
	  // Time: O(|path|)
	void setDissimMult (bool usePrediction);
	  // Input: multFixed
	  // Output: Dissim::mult, absCriterion, mult_sum, target2_sum
private:
  void setDissimMult (Dissim& dissim,
                      bool usePrediction);
	  // Input: multFixed
	  // Output: dissim::mult
public:
	  
  // Optimization	  
  bool optimizeLenWhole ();
    // Rerturn: success
	size_t optimizeLenArc ();
	  // Return: # nodes delete'd
	  // Update: DTNode::len
	  // Output: Dissim::prediction, absCriterion
	  // Time: O(p log(n))
  size_t optimizeLenNode ();
	  // Return: # nodes delete'd
	  // Update: DTNode::len
	  // Output: Dissim::prediction, absCriterion
    // After: deleteLenZero()
    // Postcondition: Dissim: prediction = 0 => target = 0 
    // Not idempotent
    // Time: O(n log^4(n))
  // Topology
	void optimize2 ();
	  // Optimal solution
	  // Requires: 2 leaves
	void optimize3 ();
	  // Optimal solution, does not depend on Obj::mult
	  // Requires: 3 leaves
  bool optimizeReinsert ();
    // Re-inserts subtrees with small DTNode::pathDissimNums.size()
	  // Return: false <=> finished
    // Invokes: NewLeaf(DTNode*), Change, applyChanges(), Threads
	void optimizeWholeIter (uint iter_max,
	                        const string &output_tree);
	  // Input: iter_max: 0 <=> infinity
	  // Update: cout
	  // Invokes: optimizeWhole(), saveFile(output_tree)
private:
	bool optimizeWhole ();
	  // Update: DTNode::stable
	  // Return: false <=> finished
	  // Requries: getConnected()
	  // Invokes: getBestChange(), applyChanges()
	  // Time of 1 iteration: O(n Time(getBestChange))  
  const Change* getBestChange (const DTNode* from);
    // Return: May be nullptr
    // Invokes: tryChange()
    // Time: O(min(n,2^areaRadius_std) log^4(n))
  bool applyChanges (VectorOwn<Change> &changes,
                     bool byNewLeaf);
	  // Return: false <=> finished
	  // Input: changes: !byNewLeaf <=> Change::apply()/restore() was done
    // Update: topology, changes (sort by Change::improvement descending)
    // Output: DTNode::stable
    // Invokes: once: finishChanges(), optimizeLen(), optimizeLenLocal(), reportErrors(cout)
	void tryChange (Change* ch,
	                const Change* &bestChange);
    // Update: bestChange: positive(improvement)
    // Invokes: Change::{apply(),restore()}
public:
  void optimizeLargeSubgraphs ();
    // Invokes: optimizeSmallSubgraphs(), Threads
	  // Time: ~ O(threads_max n log^3(n))

private:
	void optimizeSmallSubgraphs (uint areaRadius);
	  // Invokes: optimizeSmallSubgraph()
	  // Time: O(p log^2(n) * Time(optimizeSmallSubgraph))
  void optimizeSmallSubgraphsUnstable (uint areaRadius);
    // Input:: DTNode::stable
	  // Invokes: optimizeSmallSubgraph()
	  // Time: (number of !DTNode::stable node's) * Time(optimizeSmallSubgraph))
	void optimizeSmallSubgraph (const DTNode* center,
	                            uint areaRadius);
  void delayDeleteRetainArcs (DTNode* node);
    // Invokes: s->detachChildrenUp()
  size_t finishChanges ();
    // Return: deleteLenZero()
  size_t deleteLenZero ();
    // Delete arcs where len = 0
    // Does not delete root
    // Invokes: deleteLenZero(node), delayDeleteRetainArcs()
  bool deleteLenZero (DTNode* node);
    // Return: success
public:
  Real getDissimCoeffProd () const
    { Real prod = 1.0;
      for (const DissimType& dt : dissimTypes)
        if (dt. scaleCoeff)
          prod *= dt. scaleCoeff;
      return prod;
    }  
  void optimizeDissimCoeffs ();
    // Update: DissimType::scaleCoeff, Dissim::{target,mult}
private:
  Real normalizeDissimCoeffs ();
    // Return: multiplier
  void removeDissimType (size_t type);
public:
  Dataset getDissimWeightDataset (Real &dissimTypeError) const;
    // Return: attributes: "dissim", "weight"
    // Output: dissimTypeError - part of absCriterion
    // Requires: dissims.searchSorted
  void removeLeaf (Leaf* leaf,
                   bool optimizeP);
    // Invokes: leaf->detachChildrenUp(), optimizeSmallSubgraph(), toDelete.deleteData()
    // Update: detachedLeaves
    // !leaf->getParent()->childrenDiscernible(), number of children > 1 and !optimizable() => may produce incorrect tree
	  // Time: Time(optimizeSmallSubgraph)    
        
  // After optimization
  void setHeight ()
    { const_static_cast<DTNode*> (root) -> setSubtreeLenUp (false); }
    // Input: DTNode::len
    // Output: DTNode::subtreeLen
    // Time: O(n)
  void reroot (DTNode* underRoot,
               Real arcLen);
    // Invokes: sort()
  Real reroot (bool topological);
    // Center of the tree w.r.t. DTNode::setGlobalLenDown(); !topological => molecular clock
    // Return: root->getHeight()
    // Invokes: setGlobalLenDown(), reroot(,)
    // Time: O(n)
    
  // Quality
  Real getMeanResidual () const;
    // Input: Dissim::prediction
	  // Time: O(p)
  Real getMinLeafLen () const;
    // Return: min. length of discernible leaf arcs 
  Real getSqrResidualCorr () const;
    // Return: correlation between squared residual and Dissim::target
    // Input: Dissim::prediction
	  // Time: O(p)
  Real getUnoptimizable () const;
    // Return: epsilon2_0
  void setErrorDensities ();
    // Invokes: DTNode::setErrorDensity()
	  // Time: O(p log(n))
	void setLeafNormCriterion ();
    // Output: Leaf::{absCriterion,absCriterion_ave}
    // Time: O(p)
	void setNodeMaxDeformationDissimNum ();
    // Output: DTNode::maxDeformationDissimNum
    // Time: O(p log(n))
  Real getDeformation_mean () const;
    // Return: >= 0
    // Time: O(p)
  Dataset getLeafErrorDataset (bool criterionAttrP,
                               Real deformation_mean) const;
    // Input: deformation_mean: may be NaN
    // Return: attrs: PositiveAttr1 "leaf_error", "deformation" if deformation_mean is not NaN
    // Invokes: Leaf::getDeformation()
    // Requires: setLeafNormCriterion(), setNodeMaxDeformationDissimNum()
    // Time: O(n)

  // Outliers
  // Return: distinct
  VectorPtr<Leaf> findCriterionOutliers (const Dataset &leafErrorDs,
                                         Real outlier_EValue_max,
                                         Real &outlier_min_excl) const;
    // Relative average absolute criterion
    // Idempotent
    // Return: sort()'ed by Leaf::normCriterion descending
    // Output: outlier_min_excl
    // Invokes: RealAttr2::locScaleDistr2outlier()
    // Requires: after setLeafNormCriterion()
    // Time: O(n log(n))
  VectorPtr<Leaf> findDeformationOutliers (Real deformation_mean,
                                           Real outlier_EValue_max,
                                           Real &outlier_min_excl) const;
    // Return: sort()'ed by Dissim::getDeformation() descending
    // Output: outlier_min_excl
    // Invokes: Leaf::getDeformation(), MaxDistribution::getQuantileComp()
    // Time: O(n log(n))  // sorting of result
  Vector<TriangleParentPair> findHybrids (Real dissimOutlierEValue_max,
	                                        Vector<Pair<const Leaf*>> *dissimRequests) const;
    // ~Idempotent w.r.t. restoring hybrids in the tree
    // Update (append): *dissimRequests if !nullptr  // Not implemented ??
    // After: setLeafNormCriterion() 
    // Invokes: RealAttr2::normal2outlier(), findCriterionOutliers()
    // Time: O(p^2/n)
  VectorPtr<Leaf> findDepthOutliers () const;
    // Invokes: DTNode::getReprLeaf()
#if 0
  VectorPtr<DTNode> findOutlierArcs (Real outlier_EValue_max,
                                     Real &dissimOutlier_min_excl) const;
    // Output: dissimOutlier_min
#endif
    
  // Missing dissimilarities
  // Return: not in dissims; sort()'ed, uniq()'ed
  Vector<Pair<const Leaf*>> getMissingLeafPairs_ancestors (size_t depth_max,
                                                           bool refreshDissims) const;
    // Return: almost a superset of getMissingLeafPairs_subgraphs()
    // Invokes: DTNode::getSparseLeafMatches()
    // Time: ~ O(n log^2(n))
  Vector<Pair<const Leaf*>> getMissingLeafPairs_subgraphs () const;
  Vector<Pair<const Leaf*>> leaves2missingLeafPairs (const VectorPtr<Leaf> &leaves) const;
    // After: dissims.sort()

  // Clustering
//void findTopologicalClusters ();
    // Output: DisjointCluster::<Leaf>
  VectorPtr<DTNode> findDepthClusters (size_t clusters_min) const;
    // Return: connected subgraph including root
  void findGenogroups (Real genogroup_dist_max);  
    // Output: DTNode->DisjointCluster
    // For different genogroups their interior nodes do not intersect
    // Time: O(n log(n)) 

  // Statistics
#if 0
  ??
  RealAttr1* getResiduals2 ();
    // Non-weighted squared residuals
    // Return: !nullptr
  RealAttr1* getLogPredictionDiff ();
    // log(target) - log(predict);
    // Return: !nullptr
  void pairResiduals2dm (const RealAttr1* resid2Attr,
                         const RealAttr1* logDiffAttr,
                         ostream &os) const;
    // Output: os: <dmSuff>-file with attributes: dissim, distHat, resid2, logDiff
#endif

  static constexpr const char* dissimExtra {"<prediction>, <absCriterion>, <squared difference>"};
  void saveDissim (ostream &os,
                   bool addExtra) const;
    // Input: addExtra: add dissimExtra
};




///////////////////////////////////////////////////////////////////////////

struct NewLeaf : Named
// To become Leaf
// name = Leaf::name
// For Time: q = leaf2dissims.size()
{
private:
  const DistTree& tree;
  const DTNode* node_orig {nullptr};
public:
  

  struct Location : Root
  {
    const DTNode* anchor {nullptr};
      // Current best position of NewLeaf
      // = LCA of Leaf2dissim::leaf's or NewLeaf::tree.root
    Real anchorLen {NaN};
      // >= anchor->len
      // Distance to the previous anchor
    Real leafLen {NaN};
    Real arcLen {NaN};
      // From anchor to leaf->getParent()
    Real absCriterion_leaf {INF};
      // To be minimized, >= 0
      
    void qc () const override;
    void saveText (ostream &os) const override
      { const ONumber on (os, dissimDecimals, true);
        os         << anchor->getLcaName ()
           << '\t' << leafLen 
           << '\t' << arcLen  
           // Not used 
           << '\t' << anchorLen
           << '\t' << anchor->len
           << '\t' << absCriterion_leaf;
      }
      
    void setAbsCriterion_leaf (const NewLeaf& nl);
  };
  Location location;


  struct Leaf2dissim
  {
    // Input
    const Leaf* leaf {nullptr};
    // Below are functions of leaf
    Real dissim {NaN};
      // Between NewLeaf and leaf
    Real mult {NaN};
    Real absCriterion_sub {NaN};  
      // For DistTree::optimizeReinsert()
      
    // Output
    // Function of NewLeaf::Location::anchor
    Real dist_hat {0.0};
      // From leaf to NewLeaf::Location::anchor
    bool leafIsBelow {true};
    
    Leaf2dissim (const Leaf* leaf_arg,
                 Real dissim_arg,
                 Real mult_arg);
      // Input: anchor = DistTree::root
    explicit Leaf2dissim (const Leaf* leaf_arg)
      : leaf (leaf_arg)
      {}
    Leaf2dissim () = default;
      
    Real getDelta () const
      { return dissim - dist_hat; }
    Real getU () const
      { return leafIsBelow ? 1.0 : -1.0; }
    Real getEpsilon (const Location& loc) const
      { const Real leaf_dist_hat = dist_hat + loc. arcLen * getU () + loc. leafLen;
        return dissim - leaf_dist_hat;
      }
      
    bool operator< (const Leaf2dissim &other) const
      { return leaf < other. leaf; }
    bool operator== (const Leaf2dissim &other) const
      { return leaf == other. leaf; }

    static bool dissimLess (const Leaf2dissim &ld1,
                            const Leaf2dissim &ld2)
      { return ld1. dissim < ld2. dissim; }
  };
  Vector<Leaf2dissim> leaf2dissims;
    // Leaf2dissim::leaf: distinct, sort()'ed


  // Find best location, greedy
  // Time: O(q^2 log(n))
  NewLeaf (const DistTree &tree_arg,
           const string &dataDir_arg,
           const string &name_arg,
           bool init);
    // Invokes: process()
  NewLeaf (const DistTree &tree_arg,
           const string &name_arg,
           const string &dissimFName,
           const string &leafFName,
           const string &requestFName,
           bool init)
    : Named (name_arg)
    , tree (tree_arg)
    { process (init, dissimFName, leafFName, requestFName); }
  NewLeaf (const DTNode* dtNode,
           size_t q_max,
           Real &nodeAbsCriterion_old);
    // q = q_max
    // Output: nodeAbsCriterion_old: in subgraph, restricted by q_max
    // Invokes: optimize()
    // Time: O(n + p log(n) / n)
    // Cumulative time for all DTNode's: O(n log(n) + p log(n)) = O(p log(n))
private:
  void process (bool init,
                const string &dissimFName,
                const string &leafFName,
                const string &requestFName);
    // Invokes: saveLeaf(), saveRequest()
  void saveLeaf (const string &leafFName) const
    { OFStream f (leafFName);
      f << name << '\t';
      location. saveText (f);
      f << endl;
    }
  void saveRequest (const string &requestFName) const;
    // Input: location.anchor
    // Output: file requestFName
    // Invokes: DTNode::getSparseLeafMatches()
    // Time: O(log(n) (log(n) + log(q)))
  void optimize ();
    // Output: location
    // Update: leaf2dissims.{dist_hat,leafIsBelow}
    // Invokes: optimizeAnchor()
  void optimizeAnchor (Location &location_best,
                       Vector<Leaf2dissim> &leaf2dissims_best);
    // Depth-first search, greedy
    // Update: location, leaf2dissims, location_best
    // Output: leaf2dissims_best
    // Invokes: anchor2location(), descend()
  void anchor2location ();
    // Input: leaf2dissims
    // Output: location.{leafLen,arcLen,absCriterion_leaf}
    // Time: O(q)
  bool descend (const DTNode* anchorChild);
    // Return: true <=> location.anchor has leaves in leaf2dissims
    // Update: leaf2dissims, location.anchor
    // Time: O(q log(n))
public:
  void qc () const override;
};



}



#endif


