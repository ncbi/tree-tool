// distTree.hpp

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
constexpr streamsize criterionDecimals = 3;
constexpr uint areaRadius_std = 5;  
  // The greater the better DistTree::absCriterion
constexpr size_t areaDiameter_std = 2 * areaRadius_std; 
constexpr uint subgraphDepth = areaRadius_std;  
constexpr uint boundary_size_max_std = 500;  // was: 100
constexpr size_t sparsingDepth = areaDiameter_std;  // must be: >= areaRadius_std
constexpr Prob rareProb = 0.01; 
constexpr size_t dissim_progress = 1e5;



// For Time: 
//   n = # Tree leaves, p = # distances = DistTree::dissims.size()
//   p >= n
//   ~ O(): p/n = log(n); log log = 1



// --> DistTree ??
// Dissimilarity variance
enum VarianceType { varianceType_lin     // Dissimilarity ~ Poisson
                  , varianceType_exp     // Dissimilarity = -ln(P), var P = const
                  , varianceType_linExp  // Dissimilarity = -ln(P), var P = p*(1-p)
                  };
extern const StringVector varianceTypeNames;
extern VarianceType varianceType;

inline VarianceType str2varianceType (const string &s)
  { size_t index = 0;
    if (varianceTypeNames. find (s, index))
      return (VarianceType) index;
    throw logic_error ("Unknown dissimilarity variance " + s);
  }      

// Input: varianceType
inline Real dissim2mult (Real dissim)
  { switch (varianceType)
    { case varianceType_lin:    return positive (dissim) ? 1 / dissim : 0;
      case varianceType_exp:    return exp (- 2 * max (dissim, 0.0)); 
      case varianceType_linExp: return positive (dissim) ? 1 / (exp (dissim) - 1) : 0; 
    }  
    throw logic_error ("Unknown dissimilarity variance");
  }
  // Return: >= 0
  //         0 <=> dissim = INF
inline Real dissim_max ()
  { switch (varianceType)
    { case varianceType_lin:    return 1 / epsilon;
      case varianceType_exp:    return - 0.5 * log (epsilon);
      case varianceType_linExp: return log (1 / epsilon + 1);
    }  
    throw logic_error ("Unknown dissimilarity variance");
  }
  // Solution of: dissim2mult(dissim_max) = epsilon
  // dissim < dissim_max() <=> !nullReal(dissim2mult(dissim))


extern Real dissim_coeff;  // Irrelevant if varianceType = varianceType_lin
extern Real hybridness_min;
extern Real dissim_boundary;

inline bool at_dissim_boundary (Real dissim)
  { return     dissim_boundary >= dissim 
  	       && (dissim_boundary -  dissim) / dissim_boundary <= 0.05;  // PAR 
  }
  // Has a small probability <= choice of dissim_boundary



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
	Real target {NAN};
	  // > 0

	  
	Neighbor (const Leaf* leaf_arg,
	          Real target_arg);
	explicit Neighbor (const Leaf* leaf_arg);

	  
	bool operator< (const Neighbor &other) const
	  { return leaf < other. leaf; }
	bool operator== (const Neighbor &other) const
	  { return leaf == other. leaf; }
};



struct Triangle  
// Triple of Leaf's with a triangle inequality violation
{ 
	struct Parent
		{ const Leaf* leaf {nullptr};
			Real dissim {NAN};
			  // = d(Triangle::child,leaf)
			  // > 0
			bool hybrid {false};
				// Cause of the triangle inequality violation
		};

	// !nullptr
	const Leaf* child {nullptr};
	Real hybridness {NAN};
	  // = d(parent1,parent2) / (d(child,parent1) + d(child,parent2))
	  // > 1
	array<Parent, 2> parents;	
	bool child_hybrid {false};
		// Cause of the triangle inequality violation

	  
	Triangle (const Leaf* child_arg,
		        Real hybridness_arg,
				  	const Leaf* parent1,
				  	const Leaf* parent2,
				  	Real parent1_dissim,
				  	Real parent2_dissim);
	void qc () const;
	void print (ostream &os) const;
	static constexpr const char* format {"<child> <hybridness> <parent1> <parent2> <d(child,parent1)> <d(child,parent2)> <child is hybrid> <parent1 is hybrid> <parent2 is hybrid>"};
	
	
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



struct TriangleParentPair
{
	// Input
	struct Parent
		{ const Leaf* leaf {nullptr};
				// !nullptr
			size_t classSize {0};
			  // Output
		};
	array<Parent,2> parents;	
	Real parentsDissim {NAN};
	  // = f(parents[0].leaf,parents[1].leaf)
	
	// Output
	Vector<Triangle> triangles;
	  // Triangle::parents[i].leaf = parents[i].leaf
	  // Clusterize Triangle::child's ??
	  // May be empty()
private:
	size_t triangle_best_index {NO_INDEX};
	  // Index in triangles
public:
	Real hybridness_ave {NAN};

	
	TriangleParentPair (const Leaf* parent1,
		                  const Leaf* parent2,
		                  Real parentsDissim_arg)
		: parentsDissim (parentsDissim_arg)
		{ parents [0]. leaf = parent1;
			parents [1]. leaf = parent2;
	  }
	TriangleParentPair ()
	  {}
	void setTriangles (const DistTree &tree);
	  // Output: triangles, hybridness_ave
	  // Time: ~ O(log(n))
  void finish (const DistTree &tree,
               const Set<const Leaf*> &hybrids);
    // Input: triangles
    // Output: Triangle::*hybrid
	void qc () const;
  void print (ostream &os) const;
  static constexpr const char* format {"<child> <parent1> <parent2> <# children> <# parents 1> <# parents 2> <hybridness> <d(child,parent1)> <d(child,parent2)> <avg. child size> <parent1 size> <parent2 size> <child is hybrid> <parent1 is hybrid> <parent2 is hybrid>"};


  Pair<const Leaf*> getPair () const
    { return Pair<const Leaf*> (parents [0]. leaf, parents [1]. leaf); }
	bool operator< (const TriangleParentPair &other) const
    { return getPair () < other. getPair (); }
  bool operator== (const TriangleParentPair &other) const
    { return getPair () == other. getPair (); }
  static bool compareHybridness (const TriangleParentPair &hpp1,
                                 const TriangleParentPair &hpp2)
    { return hpp1. hybridness_ave > hpp2. hybridness_ave; }
    

  const Triangle& getBest () const
    { if (triangle_best_index < triangles. size ())
    	  return triangles [triangle_best_index]; 
    	throw logic_error ("TriangleParentPair::getBest()");
    }
  bool dissimError () const  
    { if (    getBest (). parent_dissim_ratio () < 0.25  // PAR
    	    && hybridness_ave < 1.25  // PAR
    	   )
    	  return true;
			for (const bool i : {false, true})
			  if (at_dissim_boundary (getBest (). parents [i]. dissim))
				  return true; 
      return false;
    }
  Vector<Triangle> getHybridTriangles () const
    { Vector<Triangle> vec;  
    	for (const Triangle &tr : triangles)
    		if (tr. hasHybrid ())
    			vec << tr;
    	return vec;
    }
  VectorPtr<Leaf> getHybrids (bool hybrid) const
    { VectorPtr<Leaf> vec;  
    	for (const Triangle &tr : triangles)
    		if (tr. hasHybrid ())
    			vec << tr. getHybrids (hybrid);
    	return vec;
    }
  void qcMatchHybrids (const VectorPtr<Leaf> &hybrids) const
    { for (const Triangle &tr : triangles)
    		if (tr. hasHybrid ())
	    		tr. qcMatchHybrids (hybrids);
    }
	  // Requires: hybrids: sort()'ed
private:
	size_t child_parent2parents (const DistTree &tree,
                               const Leaf* child,
                               const Leaf* parent,
                               Real parentDissim) const;
	  // Time: ~ O(log(n))
	void setChildrenHybrid ()
	  { for (Triangle &tr : triangles)
	  	  tr. child_hybrid = true;
	  }
};



struct DTNode : Tree::TreeNode 
{
  friend DistTree;
  friend Image;
  friend Steiner;
  friend Leaf;
  friend NewLeaf;

	Real len;
	  // Arc length between *this and *getParent()
	  // *this is root => NAN
  Vector<uint/*objNum*/> pathObjNums; 
    // Paths: function of getDistTree().dissims
    // Dissimilarity paths passing through *this arc
    // asLeaf() => getDistTree().dissims[objNum].hasLeaf(this)  
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
  size_t paths {0};  // = pathObjNums.size() ??
    // Number of paths going through *this arc
  Real errorDensity {NAN};
  Real absCriterion {NAN};
    // = sum(dissim.getAbsCriterion())
    // For Leaf's: sum = 2*getTree().absCriterion
  Real absCriterion_ave {NAN};
    // = absCriterion/count
  Real relCriterion {NAN};
    // !isNan() => = getRelCriterion()


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
	  { return isNan (len) ? 0 : len; }

  const DistTree& getDistTree () const;

  const Leaf* inDiscernible () const;
    // Return: this or nullptr
  bool childrenDiscernible () const
    { return arcs [false]. empty () || ! static_cast <DTNode*> ((*arcs [false]. begin ()) -> node [false]) -> inDiscernible (); }
  Real getHeight_ave () const
    { return subtreeLen. getMean (); }    
    // After: DistTree::setHeight()
  Real getEpsilon2 () const
    { return (Real) paths * sqr (errorDensity) * len; }
  Real getRelCriterion () const;
  virtual const Leaf* getReprLeaf () const = 0;
    // Return: !nullptr, in subtree
    // For sparse *getDistTree().dissimAttr
    // Invokes: getDistTree().rand
    // Time: ~ O(log(n))
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
    // Update: descendants
  virtual void setLca ();
    // Off-line LCA algorithm by Tarjan
    // Requires: !tarjanLca 
  Vector<uint/*objNum*/> getLcaObjNums ();
    // Return: objNum's s.t. getDistTree().dissims[objNum].lca = this
    // Invokes: DTNode::pathObjNums.sort()
  VectorPtr<Leaf> getSparseLeafMatches (size_t depth_max,
                                        bool subtractDissims) const;
    // Return: size = O(log(n)); sort()'ed, uniq()'ed
    //         getDistTree().reroot(true) reduces size()
    // Input: depth_max: 0 <=> no restriction
    // Time: O(log^2(n)) 
};



struct Steiner : DTNode
// Steiner node
{
#if 0
private:
  size_t threadNum {NO_INDEX};
  // Temporary
  size_t lcaNum {0};  
  size_t subTreeWeight {0};  // union ??
public:
#endif


	Steiner (DistTree &tree,
	         Steiner* parent_arg,
	         Real len_arg);
	void qc () const override;
#if 0
  void saveContent (ostream& os) const final
    { DTNode::saveContent (os);
    	if (frequentChild)
    	  os << " " << "frequent";
    	os << " fd=" << frequentDegree;
    }
#endif


  const Steiner* asSteiner () const final
    { return this; }

  bool isInteriorType () const final
    { return childrenDiscernible (); }

private:
  const Leaf* getReprLeaf () const final;
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
{
	friend DistTree;
	
  string name;  
    // !empty()
  string comment;
  static const string non_discernible;
  bool discernible {true}; 
    // false => getParent()->getChildren() is an equivalence class of indiscernibles
    
private:
  // Temporary
  size_t index {NO_INDEX};
public:

  // Hybrid data
private:
	friend Triangle;
  Vector<Neighbor> badNeighbors;
    // Neighbor::leaf is unique
  Real hybridness {1};
    // >= 1    
    // = max d(parent1,parent2) / (d(this,parent1) + d(this,parent2))
    // > 1 => triangle inequality violation
  uint hybridParentsDissimObjNum;  
    // hybridness > 1 => in DistTree::dissims
public:
  

	Leaf (DistTree &tree,
	      Steiner* parent_arg,
	      Real len_arg,
	      const string &name_arg);
	Leaf (DistTree &tree,
	      Leaf* other,
	      const string &name_arg);
	  // Invokes: collapse(other)
	void qc () const final;
  void saveContent (ostream& os) const final
    { DTNode::saveContent (os);
    	if (! discernible)
    	  os << "  " << non_discernible;
    }


  const Leaf* asLeaf () const final
    { return this; }


  string getName () const final
    { return name; }
  string getNewickName (bool minimal) const final
    { if (minimal)
        return name;
      string s = name + prepend (" ", comment); 
      if (! isNan (relCriterion))
        s += " " + real2str (relCriterion, 1);  // PAR
      else if (! isNan (absCriterion))
        s += " " + real2str (getRelCriterion (), 1);  // PAR
      return s;
    }
  bool isLeafType () const final
    { return true; }

private:
  const Leaf* getReprLeaf () const final
    { return this; }
  void setSubtreeLenUp (bool topological) final
    { subtreeLen. clear ();
      if (topological)
    	  subtreeLen. add (0, 1);
    }
  void getDescendants (VectorPtr<DTNode> &descendants,
                       size_t /*depth*/,
                       const DTNode* exclude) const final
    { if (this != exclude)
    	  descendants << this; 
    }
  void setLca () final;
public:

  const Leaf* getDissimOther (size_t objNum) const;
    // Return: !nullptr; != this
  const DTNode* getDiscernible () const;
    // Return: this or getParent(); !nullptr
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
    // Invokes: setParent()
    // To be followed by: DistTree::cleanTopology()
public:	
#if 0
	bool possibleHybrid () const
	  { return getDiscernible () -> len == 0; }
#endif
private:
	const Neighbor* findNeighbor (const Leaf* leaf) const
	  { const size_t i = badNeighbors. binSearch (Neighbor (leaf));
    	if (i == NO_INDEX)    		
    		return nullptr;
    	return & badNeighbors [i];
    }
public:
  Real getHybridness (/*Real height_min,*/
                    //Real parentLengthFrac_min,
	                    const Leaf* &parent1,
	                    const Leaf* &parent2,
	                    Real &parent1_dissim,
	                    Real &parent2_dissim) const;
		// Return: 0: not a hybrid
		//         > 1: hybrid
		//                = max d(parent1,parent2) / (d(this,parent1) + d(this,parent2))
		//                triangle inequality violation
		//         NAN: request of the dissimilarity between parent1 and parent2, if (bool)parent1
	  // Input: parentLengthFrac_min > 1
	  // Time: ~ O(log^3(n))
	void getHybridTriangles (Vector<Triangle> &hybrids,
		                       Set<const Leaf*> &tried) const;
		// Update: hybrids, tried; hybrids is a subset of tried
		// Invokes: getHybridness()
		// Recursive
};



struct SubPath
// Path going through a connected subgraph
{
  uint objNum {(uint) -1};    
    // Index of DistTree::dissims
  const DTNode* node1 {nullptr};
  const DTNode* node2 {nullptr};
    // !nullptr, different
  Real dist_hat_tails {NAN};

    
  SubPath ()
    {}
  explicit SubPath (uint objNum_arg)
    : objNum (objNum_arg)
    {}
  void qc () const;

  
  bool operator== (const SubPath &other) const
    { return objNum == other. objNum;
    }
  bool operator< (const SubPath &other) const
    { return objNum < other. objNum; }
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
  VectorPtr<Tree::TreeNode> boundary;
    // Make boundary and area disjoint ??
  const Steiner* area_root {nullptr};
    // May be nullptr
    // boundary.contains(area_root)
  const DTNode* area_underRoot {nullptr};
    // May be nullptr
    // area.contains(area_underRoot)
  // (bool)area_underRoot = (bool)area_root
  Vector<SubPath> subPaths;
    // Some paths of tree.dissims passing through area
    // Size: O(|bounadry| p/n log(n)) ~ O(|area| log^2(n))
  Real subPathsAbsCriterion {0};
  
  
  explicit Subgraph (const DistTree &tree_arg);
  void qc () const override;
  bool empty () const override
    { return    area. empty ()
             && boundary. empty ()
             && ! area_root
             && ! area_underRoot
             && subPaths. empty () 
             && ! subPathsAbsCriterion;
    }
  void clear () override
    { area. wipe ();
      boundary. wipe ();
      area_root = nullptr;
      area_underRoot = nullptr;
      subPaths. wipe ();
      subPathsAbsCriterion = 0;
    }

  
  // Usage:
//set area, boundary
  void reserve ();
  void removeIndiscernibles ();
    // Update: area, boundary
  void finish ();
    // Time: O(|area| log(|area|))
//set subPaths
  void finishSubPaths ();
    // Time: O(|boundary| p/n log(n) log(|subPaths|) + |subPaths| log(|area|)) 
    //       ~ O(|area| log(|area| log^2(n)))
//change topology of tree within area
  Real getImprovement (const DiGraph::Node2Node &boundary2new) const;
    // Time: O(|subPaths| (log(|boundary|) + log(|area|)))
  void subPaths2tree ();
    // Update: tree: Paths, absCriterion, Dissim::prediction
    // Time: O(p/*fast*/ + |area| (log(|boundary|) + p/n log(n)) + |subPaths| log(|area|) log(|boundary|)) 
    //       ~ O(|area| log^2(|area|) log^2(n))

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
  void node2subPaths (const DTNode* node);
    // Time: O(p/n log(n))
  void area2subPaths ();
    // Invokes: node2subPaths()
    // Time: O(|boundary| p/n log(n))
  VectorPtr<Tree::TreeNode>& getPath (const SubPath &subPath,
  	                                  Tree::LcaBuffer &buf) const
    // Return: reference to buf
    { const Tree::TreeNode* lca_ = nullptr;
      // tree.dissims[subPath.objNum].lca can be used instead of area_root if viaRoot(subPath) and tree topology has not been changed ??
      return Tree::getPath (subPath. node1, subPath. node2, area_root, lca_, buf);
    }
    // Requires: subPath in subPaths
  const Vector<uint>& boundary2pathObjNums (const DTNode* dtNode) const
    { return dtNode == area_root /* && area_underRoot */
               ? area_underRoot->pathObjNums
               : dtNode        ->pathObjNums;
    }
  const Leaf* getReprLeaf (const DTNode* dtNode) const
    { return dtNode == area_root 
               ? static_cast <const DTNode*> (dtNode->getDifferentChild (area_underRoot)) -> getReprLeaf ()
               : dtNode->getReprLeaf ();
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
	Real improvement {NAN};
	  // isNan() or positive()
    // Too small values are noise => not stable in tree sampling
private:
	Real fromLen {NAN};
	Real toLen {NAN};
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
		: tree (const_cast <DistTree&> (from_arg->getDistTree ()))
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
	static bool strictlyLess (const Change* a, 
	                          const Change* b);
    // Requires: (bool)a
};



struct Dissim
{
	// Input
  // !nullptr
  // leaf1->name < leaf2->name
  const Leaf* leaf1 {nullptr};
  const Leaf* leaf2 {nullptr};
  //
  Real target {NAN};
    // Dissimilarity between leaf1 and leaf2; !isNan()
  Real mult {NAN};
  
  // Output
  Real prediction {NAN};
    // Tree distance
  const Steiner* lca {nullptr};
    // Paths
  

  Dissim (const Leaf* leaf1_arg,
          const Leaf* leaf2_arg,
          Real target_arg = NAN,
          Real mult_arg = 0);
  void qc () const;

          
  bool valid () const
    { return    leaf1->graph
             && leaf2->graph;
    }
    // For topology
  bool hasLeaf (const Leaf* leaf) const
    { return    leaf == leaf1
             || leaf == leaf2;
    }
  const Leaf* getOtherLeaf (const Leaf* leaf) const
    { if (leaf == leaf1) return leaf2;
    	if (leaf == leaf2) return leaf1;
    	throw logic_error ("getOtherLeaf");
    }
  string getObjName () const;
  void print () const
    { cout << leaf1 << ' ' << leaf2 << ' ' << mult << ' ' << lca << endl; }
  VectorPtr<Tree::TreeNode>& getPath (Tree::LcaBuffer &buf) const;
  	// Return: reference to buf
  Real getResidual () const
    { return prediction - target; }
  Real getAbsCriterion (Real prediction_arg) const;
  Real getAbsCriterion () const
    { return getAbsCriterion (prediction); }
    // LSE is MLE <= sqrt(mult)*(prediction-target) ~ N(0,sigma^2)
  Real process (size_t objNum,
                Tree::LcaBuffer &buf);
    // Return: getAbsCriterion()
    // Output: prediction, Steiner::pathObjNums
    
  bool operator< (const Dissim &other) const
    { return Pair<const Leaf*> (leaf1, leaf2) < Pair<const Leaf*> (other. leaf1, other. leaf2); }
  bool operator== (const Dissim &other) const
    { return Pair<const Leaf*> (leaf1, leaf2) == Pair<const Leaf*> (other. leaf1, other. leaf2); }
};




struct Image : Nocopy
// Tree subgraph replica
{
  Subgraph subgraph;
  const DTNode* center {nullptr};
    // In subgraph.tree
    // May be delete'd
  DistTree* tree {nullptr};
    // Subgraph tree
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
    // Time: ~ O(n^2 / threads_max^2)
  void apply ();
	  // Output: DTNode::stable = true
    // Time: ~ O(|area| log^2(|area|) log^2(n))
};




///////////////////////////////////////////////////////////////////////////

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
  map<string/*Leaf::name*/,const Leaf*> name2leaf;
    // 1-1

private:
  // Temporary
  // Dissimilarity
  // May be nullptr
  Common_sp::AutoPtr<Dataset> dissimDs;
    // Original data
  const PositiveAttr2* dissimAttr {nullptr};
    // In *dissimDs
public:
    
  constexpr static uint dissims_max {numeric_limits<uint>::max ()};
  Vector<Dissim> dissims;
  Real mult_sum {0};
  Real dissim2_sum {0};
    // = sum_{dissim in dissims} dissim.target^2 * dissim.mult        
  Real absCriterion {NAN};
    // = L2LinearNumPrediction::absCriterion  
private:
  size_t leafNum {0};
    // For Leaf::index
	VectorOwn<DTNode> toDelete;
	VectorOwn<Leaf> detachedLeaves;
	  // !Leaf::graph
	mutable Rand rand;
public:


  // Input: dissimFName: <dmSuff>-file without <dmSuf>, contains attribute dissimAttrName
  //                     may contain more objects than *this contains leaves
  //        dissimFName and dissimAttrName: both may be empty	  
	DistTree (const string &treeFName,
	          const string &dissimFName,
	          const string &dissimAttrName,
	          bool sparse);
	  // Invokes: loadTreeFile(), loadDissimDs(), dissimDs2dissims()
	DistTree (const string &dirName,
	          const string &dissimFName,
	          const string &dissimAttrName);
	  // Input: dirName: contains the result of mdsTree.sh; ends with '/'
	  // Invokes: loadTreeDir(), loadDissimDs(), dissimDs2dissims(), setGlobalLen()
	DistTree (const string &dissimFName,
	          const string &dissimAttrName,
	          bool sparse);
	  // Invokes: loadDissimDs(), dissimDs2dissims(), neighborJoin()
	DistTree (const string &dataDirName,
            bool loadNewLeaves,
	          bool loadDissim,
	          bool optimizeP);
	  // Input: dataDirName: ends with '/': directory for incremental data structure:
	  //          file name                 line/file format                              meaning
	  //          ---------                 -----------------------------                 ----------------------------
	  //          tree                                           
	  //          dissim                    <obj1> <obj2> <dissimilarity>                 <ob1>, <ob2> are tree leaves
    //          leaf                      <obj_new> <obj1>-<obj2> <leaf_len> <arc_len>
    //         [dissim.add[-req]]
    //          new/                      <obj>                                         New objects to add to the tree
    //          search/<obj_new>/                                                       Initialization of search for <obj_new> location
    //        [ search/<obj_new>/dissim   <obj_new> <obj> <dissimilarity>        
    //          search/<obj_new>/leaf     = as in leaf 
    //         [search/<obj_new>/request  <obj_new> <obj>]                              Request to compute dissimilarity
    //        ]
    //         [alien]                                                                  Objects similar to no other objects
	  //          outlier/                  <obj>                                         Tree outlier objects
	  //         [dissim_request]           <obj1> <obj2>                                 Request to compute dissimilarity
	  //          version                   <natural number>
	  //          hist/                     {tree,makeDistTree,leaf,outlier,
	  //                                     makeFeatureTree
	  //                                    }.<version>                                   Historic versions of data
    //          grid_min                  <number>                                      Min. number of dissimilarity requests to be processed on a grid (e.g., 2000 for fast dissimilarities)
	  //         [phen/]                                                                  Link to a directory with phenotypes for makeFeatureTree
	  //          runlog                                                                  Invocations of distTree_inc_new.sh
	  //         [stop]                     /dev/null                                     Stop distTree_inc_new.sh
    //          request2dissim.sh         executable with parameters: request, dissim.add, log
    //          objects_in_tree.sh        executable with parameters: list of objects, 0/1/null
    //          request_closest.sh        executable with parameter: object; output: pairs of objects to request dissimilarities for
    //         [hybrid identification]
	  //          hybridness_min            <number>                                      Min. hybridness, >1
	  //          dissim_boundary           <number>                                      Boundary between two merged dissmilarity measures causing discontinuity
    //          hybrid2db.sh              executable with parameter: tab-delimited hybrid file (outpuf of makeDistTree -find_hybrids)
    //         -db2hybrid.sh              executable with parameter: <Obj>; output: empty or <inter_parent_dissim>\n<parent1>\n<parent2>
    //          db2unhybrid.sh            executable: output: <Obj> list to move to new/
	  //       <dissimilarity>: >= 0, < INF
	  // Invokes: optimizeSmallSubgraph() for each added Leaf; Threads
	  // Time: if loadDissim then O(p log(n) + Time(optimizeSmallSubgraph) * new_leaves)
	  //       if !loadDissim then O(n + new_leaves)
  //  
  explicit DistTree (const string &newickFName);
  DistTree (Prob branchProb,
            size_t leafNum_max);
    // Random tree: DTNode::len = 1
    // Time: O(n)
  DistTree (Subgraph &subgraph,
            Node2Node &newLeaves2boundary);
    // Connected subgraph of subgraph.tree: boundary of subgraph.area are Leaf's of *this
    // If subgraph.unresolved() then the topology of *this is changed to a star
    // Input: subgraph: !empty(), not finish()'ed
    // Output: subgraph: area: contains newLeaves2boundary.values(); discernible
    //         newLeaves2boundary
	  // Time: ~ O(|area| (log(|area|) log^2(subgraph.tree.n) + |area| /*fast*/))
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
  void newick2node (ifstream &f,
                    Steiner* parent);
  void setName2leaf ();
  void loadDissimDs (const string &dissimFName,
                     const string &dissimAttrName);
    // Output: dissimDs
    // invokes: dissimDs->setName2objNum()
  // Input: dissimDs
  bool getConnected ();
    // Find connected components of leaves where pairs have dissimilarities with positive multiplicity
    // Return: true <=> 1 connected component
    // Output: DisjointCluster
    //         cout: print other connected components
  size_t setDiscernible ();
    // Return: Number of indiscernible leaves
    // Output: Leaf::len = 0, Leaf::discernible = false, topology
    // Invokes: cleanTopology()
  void cleanTopology ();
  void setGlobalLen ();
    // Molecular clock 
    // Output: DTNode::len
    // Temporary: DTNode::subtreeLen
	  // Time: O(p log(n))
  void neighborJoin ();
    // Greedy
    // Assumes: Obj::mult = 1
    // Requires: complete *dissimAttr, isStar()
    // Invokes: reroot(true)
    // Time: O(n^3)
  //
  void dissimDs2dissims (bool sparse);
    // Update: dissimDs: delete
    // Output: dissims etc.
    //         if an object is absent in dissimDs then it is deleted from the Tree
    // Invokes: getSelectedPairs(), setPaths()
  void loadDissimPrepare (size_t pairs_max);
    // Output: Dissim::target
  uint leaves2dissims (Leaf* leaf1,
                       Leaf* leaf2,
                       Real target,
                       Real mult);
    // Return: dissims.size() - 1
    // Append: dissims[], Leaf::pathObjNums
  bool addDissim (Leaf* leaf1,
                  Leaf* leaf2,
                  Real dissim);
	  // Return: dissim is added
	  // Update: Dissim, mult_sum, dissim2_sum
	  // Invokes: leaves2dissims()
  bool addDissim (const string &name1,
                  const string &name2,
                  Real dissim)
    { return addDissim ( const_cast <Leaf*> (findPtr (name2leaf, name1))
    	                 , const_cast <Leaf*> (findPtr (name2leaf, name2))
    	                 , dissim
    	                 );
    }
  void setPaths ();
    // Output: dissims::Dissim, DTNode::pathObjNums, absCriterion
    // Invokes: setLca()
    // Time: O(p log(n))
public:
	void qc () const override;


  void deleteLeaf (TreeNode* leaf,
                   bool deleteTransientAncestor) final;
    // Requires: !optimizable()
    
  static string getObjName (const string &name1,
                            const string &name2);
  const DTNode* lcaName2node (const string &lcaName,   
                              Tree::LcaBuffer &buf) const;
    // Return: !nullptr
    // Input: lcaName: <leaf1 name> <objNameSeparator> <leaf2 name>
  size_t extraObjs () const
    { return dissimDs. get () ? dissimDs->objs. size () - name2leaf. size () : 0; }
  size_t getDissimSize_max () const
    { return name2leaf. size () * (name2leaf. size () - 1) / 2; }	
  size_t getSparseDissims_size () const
    { return ((size_t) log (name2leaf. size ()) + 1) * 70; }  // PAR
  Set<const DTNode*> getDiscernibles () const;
    // Logical leaves
  static void printParam (ostream &os) 
    { os << "PARAMETERS:" << endl;
      os << "# Threads: " << threads_max << endl;
      os << "Dissimilarity variance: " << varianceTypeNames [varianceType] << endl;
      os << "Max. possible dissimilarity: " << dissim_max () << endl;
      os << "Subgraph radius: " << areaRadius_std << endl;
    }
	void printInput (ostream &os) const;
	bool optimizable () const  
	  { return ! dissims. empty (); }
	Real getDissim_ave () const
	  { WeightedMeanVar mv;
	    for (const Dissim& dissim : dissims)
	      mv. add (dissim. target, dissim. mult);
	    return mv. getMean ();
	  }
  Real getAbsCriterion_ave () const
    { return absCriterion / (Real) dissims. size (); }
    // Approximate: includes !Dissim::valid() ?? 
  Prob getUnexplainedFrac () const
    { return absCriterion / dissim2_sum; }
  Real getErrorDensity () const
    { const Prob r = 1 - getUnexplainedFrac ();   
      return sqrt ((1 - r) / r);   
    }
    // Requires: interpretation <= linear variance of dissimilarities  ??
  string absCriterion2str () const
    { return real2str (absCriterion, criterionDecimals); }
  void reportErrors (ostream &os) const
    { const ONumber on (os, criterionDecimals, false);  // PAR
      os << "absCriterion = " << absCriterion 
         << "  Error density = " << getErrorDensity () * 100 << " %"
         << endl;
    }    
  void saveFeatureTree (const string &fName) const;

private:
  void qcPaths ();
    // Sort: DTNode::pathObjNums 
  void setLca ();
    // Input: Leaf::pathObjNums
    // Output: Dissim::lca
  void clearSubtreeLen ();
    // Invokes: DTNode::subtreeLen.clear()
  void setPredictionAbsCriterion ();
    // Output: Dissim::prediction, absCriterion
    // Invokes: Threads
    // Time: O(p log(n) / threads_max)
  void qcPredictionAbsCriterion () const;
public:
  static Real path2prediction (const VectorPtr<TreeNode> &path);
    // Return: >= 0
	  // Input: DTNode::len
	  // Time: O(|path|)
  void printAbsCriterion_halves () const;
	void setNodeAbsCriterion ();
    // Output: DTNode::{absCriterion,absCriterion_ave}
    // Time: O(p log(n))
	  
  // Optimization	  
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
    // Re-inserts subtrees with small DTNode::pathObjNums.size()
	  // Return: false <=> finished
    // Invokes: NewLeaf(DTNode*), Change, applyChanges()
    // Time: O(n log^3(n))
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
    // Time: ~ O(n^2 / threads_max^2 + threads_max n log^4(n))
private:
	void optimizeSmallSubgraphs (uint areaRadius,
	                             bool unstableOnly);
	  // Invokes: optimizeSmallSubgraph()
	  // Time: ~ O(n * Time(optimizeSmallSubgraph))
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
  void removeLeaf (Leaf* leaf,
                   bool optimizeP);
    // Invokes: leaf->detachChildrenUp(), optimizeSmallSubgraph(), toDelete.deleteData()
    // Update: detachedLeaves
	  // Time: Time(optimizeSmallSubgraph)    
        
  // After optimization
  void setHeight ()
    { const_static_cast<DTNode*> (root) -> setSubtreeLenUp (false); }
    // Input: DTNode::len
    // Output: DTNode::subtreeLen
    // Time: O(n)
  void reroot (DTNode* underRoot,
               Real arcLen);
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
  Real setErrorDensities ();
    // Requires: linear variance of dissimilarities
    // Return: epsilon2_0
    // Input: Dissim::prediction
    // Output: DTNode::{paths,errorDensity}
    // Requires: Leaf::discernible is set
	  // Time: O(p log(n))

  // Outliers
  // Return: distinct
  Dataset getLeafErrorDataset () const;
    // Return: attrs = {PositiveAttr1}
    // After: setNodeAbsCriterion()    
  VectorPtr<Leaf> findCriterionOutliers (Real outlier_EValue_max,
                                         Real &outlier_min) const;
    // Idempotent
    // Output: outlier_min
    // Invokes: getLeafErrorDataset(), RealAttr2::normal2outlier() 
    // Time: O(n log(n))
  Vector<TriangleParentPair> findHybrids (Real dissimOutlierEValue_max,
	                                        Vector<Pair<const Leaf*>> *dissimRequests) const;
    // ~Idempotent w.r.t. restoring hybrids in the tree
    // Update (append): *dissimRequests if !nullptr
    // Invokes: RealAttr2::normal2outlier() 
    // Time: ~ O(n log^2(n)) without dissimRequests
    // Time: ~ O(n log^3(n)) with    dissimRequests
#if 0
  Vector<Triangle> findHybrids (Real parentLengthFrac_min,
                                Vector<Pair<const Leaf*>> &dissimRequests) const;
	  // Input: parentLengthFrac_min > 1
    // Update (append): dissimRequests
    // Invokes: Leaf::getHybridness()
    // Time: ~ O(n log^3(n))
#endif
  VectorPtr<Leaf> findDepthOutliers () const;
    // Invokes: DTNode;:getReprLeaf()
#if 0
  VectorPtr<DTNode> findOutlierArcs (Real outlier_EValue_max,
                                     Real &dissimOutlier_min) const;
    // Output: dissimOutlier_min
#endif
    
  // Missing dissimilarities
  // Return: not in dissims; sort()'ed, uniq()'ed
  Vector<Pair<const Leaf*>> getMissingLeafPairs_ancestors (size_t depth_max) const;
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
  void findSpecies (Real species_dist_max);
    // Output: DTNode::DisjointCluster

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

  void saveDissim (ostream &os) const;
};




///////////////////////////////////////////////////////////////////////////

struct NewLeaf : Named
// To become Leaf
// name = Leaf::name
// For time: q = leaf2dissims.size()
{
private:
  const DistTree& tree;
    // !nullptr
  const string dataDir;
  const DTNode* node_orig {nullptr};
public:
  

  struct Location : Root
  {
    const DTNode* anchor {nullptr};
      // Current best position of NewLeaf
      // New Leaf* leaf: leaf->getParent() is on the anchor arc
    Real leafLen {NAN};
    Real arcLen {NAN};
      // From anchor to leaf->getParent()
    Real absCriterion_leaf {INF};
      // To be minimized, >= 0
      
    void qc () const;
    void saveText (ostream &os) const
      { const ONumber on (os, dissimDecimals, true);
        os         << anchor->getLcaName ()
           << '\t' << leafLen 
           << '\t' << arcLen  
           // Not used 
           << '\t' << anchor->len
           << '\t' << absCriterion_leaf;
      }
      
    void setAbsCriterion (const NewLeaf& nl);
  };
  Location location;


  struct Leaf2dissim
  {
    const Leaf* leaf {nullptr};
    Real dissim {NAN};
      // Between NewLeaf and leaf
    Real mult {NAN};
    Real treeAbsCriterion {NAN};
      
    // Function of NewLeaf::Location::anchor
    Real dist_hat {NAN};
      // From leaf to NewLeaf::Location::anchor
    bool leafIsBelow {false};
    
    Leaf2dissim (const Leaf* leaf_arg,
                 Real dissim_arg,
                 Real mult_arg,
                 const DTNode* anchor);
      // Input: anchor = NewLeaf::Location::anchor
    explicit Leaf2dissim (const Leaf* leaf_arg)
      : leaf (leaf_arg)
      {}
    Leaf2dissim ()
      {}
      
    Real getDelta () const
      { return dissim - dist_hat; }
    Real getU () const
      { return leafIsBelow ? 1 : -1; }
    Real getEpsilon (const Location& loc) const
      { const Real leaf_dist_hat = dist_hat + loc. arcLen * getU () + loc. leafLen;
        return dissim - leaf_dist_hat;
      }
      
    bool operator< (const Leaf2dissim &other) const
      { return leaf < other. leaf; }
    bool operator== (const Leaf2dissim &other) const
      { return leaf == other. leaf; }

    static bool multLess (const Leaf2dissim &ld1,
                          const Leaf2dissim &ld2)
      { return ld1. mult > ld2. mult; }
  };
  Vector<Leaf2dissim> leaf2dissims;
    // Leaf2dissim::leaf: distinct, sort()'ed


  // Find best location, greedy
  NewLeaf (const DistTree &tree_arg,
           const string &dataDir_arg,
           const string &name_arg,
           bool init);
    // Invokes: saveLeaf(), saveRequest()
  NewLeaf (const DTNode* dtNode,
           size_t q_max,
           Real &nodeAbsCriterion_old);
    // q = q_max
private:
  string getNameDir      () const  { return dataDir + "/" + name + "/"; }
  string getDissimFName  () const  { return getNameDir () + "dissim"; }
  string getLeafFName    () const  { return getNameDir () + "leaf"; }
  string getRequestFName () const  { return getNameDir () + "request"; }
  void saveLeaf () const
    { OFStream f (getLeafFName ());
      f << name << '\t';
      location. saveText (f);
      f << endl;
    }
  void saveRequest () const;
    // Input: location.anchor
    // Output: file getRequestFName()
    // Invokes: DTNode::getSparseLeafMatches()
    // Time: O(log(n) log(q))
  void optimize ();
    // Output: location
    // Update: leaf2dissims.{dist_hat,leafIsBelow}
    // Invokes: optimizeAnchor()
  void optimizeAnchor (Location &location_best,
                       Vector<Leaf2dissim> &leaf2dissims_best);
    // Depth-first search, greedy
    // Update: location, leaf2dissims, location_best
    // Output: leaf2dissims_best
    // Invokes: setLocation(), descend()
    // Time: O(q log^2(n))
  void setLocation ();
    // Output: location.{leafLen,arcLen,absCriterion_leaf}
    // Time: O(q)
  bool descend (const DTNode* anchor_new);
    // Return: true <=> location.anchor has leaves in leaf2dissims
    // Update: leaf2dissims
    // Output: location.anchor
    // Time: O(q log(n))
public:
  void qc () const override;
};



}



#endif


