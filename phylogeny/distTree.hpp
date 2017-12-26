// distTree.hpp

#ifndef DISTTREE_HPP
#define DISTTREE_HPP

#include "../common.hpp"
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
      case varianceType_exp:    return exp (- 2 * max (0.0, dissim)); 
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



// PAR
constexpr streamsize dissimDecimals = 6;
constexpr streamsize criterionDecimals = 3;
constexpr uint areaRadius_std = 5;  
  // The greater the better DistTree::absCriterion
constexpr uint subgraphDepth = areaRadius_std;  
constexpr uint boundary_size_max_std = 500;  // was: 100
constexpr size_t sparsingDepth = 2 * areaRadius_std;  
constexpr Prob rareProb = 0.01; 



struct DistTree;

//struct DTNode;
  struct Steiner;
  struct Leaf;

struct NewLeaf;



// For Time: n = # Tree leaves, p = # distances = DistTree::dissims.size()
//           p >= n



struct DTNode : Tree::TreeNode 
{
  friend DistTree;
  friend Steiner;
  friend Leaf;
  friend NewLeaf;

	Real len;
	  // Arc length between *this and *getParent()
  Vector<size_t/*objNum*/> pathObjNums; 
    // Paths: function of getDistTree().dissims
    // Dissimilarity paths passing through *this arc
    // asLeaf() => getDistTree().dissims[objNum].hasLeaf(this)  
    //             aggregate size = 2 p
    // !asLeaf(): aggregate size = O(p log(n))
    // Distribution of size ??
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
  const Leaf* reprLeaf {nullptr};
    // In subtree
    // For sparse *getDistTree().dissimAttr


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

  const Leaf* inDiscernable () const;
    // Return: this or nullptr
  bool childrenDiscernable () const
    { return arcs [false]. empty () || ! static_cast <DTNode*> ((*arcs [false]. begin ()) -> node [false]) -> inDiscernable (); }
  Real getHeight_ave () const
    { return subtreeLen. getMean (); }    
    // Requires: after DistTree::setHeight()
  Real getEpsilon2 () const
    { return (Real) paths * sqr (errorDensity) * len; }
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
  virtual void setRepresentative () = 0;
    // Output: reprLeaf
    // Requires: after getDistTree().sort()
  virtual void getDescendants (VectorPtr<DTNode> &descendants,
                               size_t depth) const = 0;
    // Update: descendants
  virtual void setLca ();
    // Off-line LCA algorithm by Tarjan
    // Requires: !tarjanLca 
  Vector<size_t/*objNum*/> getLcaObjNums ();
    // Return: objNum's s.t. getDistTree().dissims[objNum].lca = this
    // Invokes: DTNode::pathObjNums.sort()
  VectorPtr<Leaf> getSparseLeafMatches (size_t depth_max) const;
    // Return: size = O(log(n)); !contains(this->reprLeaf); sort()'ed, uniq()'ed
    //         getDistTree().reroot(true) reduces size()
    // Input: depth_max: 0 <=> no restriction
    // Requires: after getDistTree().setReprLeaves()
};



struct Steiner : DTNode
// Steiner node
{
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
    { return childrenDiscernable (); }

private:
  void setSubtreeLenUp (bool topological) final;
  void setRepresentative () final
    { reprLeaf = nullptr;
      for (const DiGraph::Arc* arc : arcs [false])
      { DTNode* node = static_cast <DTNode*> (arc->node [false]);
        node->setRepresentative ();
        if (! reprLeaf)
          reprLeaf = node->reprLeaf;
      }
    }
  void getDescendants (VectorPtr<DTNode> &descendants,
                       size_t depth) const final;
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
};



struct Leaf : DTNode
{
  string name;  
    // !empty()
  string comment;
  
  static const string non_discernable;
  bool discernable {true}; 
    // false => getParent()->getChildren() is an equivalence class of indiscernables
    
  // Temporary
  Real absCriterion {NAN};
    // sum = contribution to 2*getTree().absCriterion
  Real absCriterion_ave {NAN};
  Real relCriterion {NAN};
    // !isNan() => = getRelCriterion()
  size_t index {NO_INDEX};
  

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
      if (! isNan (absCriterion))
      { const ONumber oNum (os, criterionDecimals, true);  // PAR
        os << "  leaf_error=" << getRelCriterion ();
      }
    	if (! discernable)
    	  os << "  " << non_discernable;
    #if 0
      if (frequent)
    	  os << " " << "frequent";
     	os << " fd=" << frequentDegree;
    #endif
    }


  const Leaf* asLeaf () const final
    { return this; }


  string getName () const final
    { return name; }
  string getNewickName (bool minimal) const final
    { if (minimal)
        return name;
      string s = name + prepend (" ", comment); 
      if (! isNan (absCriterion))
        s += " " + real2str (getRelCriterion (), 1);  // PAR
      return s;
    }
  bool isLeafType () const final
    { return true; }

private:
  void setSubtreeLenUp (bool topological) final
    { subtreeLen. clear ();
      if (topological)
    	  subtreeLen. add (0, 1);
    }
  void setRepresentative () final
    { reprLeaf = this; }
  void getDescendants (VectorPtr<DTNode> &descendants,
                       size_t /*depth*/) const final
    { descendants << this; }
  void setLca () final;
public:

  const Leaf* getDissimOther (size_t objNum) const;
    // Return: !nullptr; != this
  const DTNode* getDiscernable () const;
    // Return: this or getParent(); !nullptr
  Real getRelCriterion () const;
  bool getCollapsed (const Leaf* other) const
    { return    other
             && getParent () == other->getParent ()
             && ! discernable
             && ! other->discernable;
    }
private:
  friend DistTree;
  void collapse (Leaf* other);
    // Output: discernable = false
    // Invokes: setParent()
    // To be followed by: DistTree::cleanTopology()
public:
};



struct SubPath
// Path going through a connected subgraph
{
  size_t objNum {NO_INDEX};    
    // Index of DistTree::dissims
  const DTNode* node1 {nullptr};
  const DTNode* node2 {nullptr};
    // !nullptr, different
  Real dist_hat_tails {NAN};

    
  SubPath ()
    {}
  explicit SubPath (size_t objNum_arg)
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
    // Size: O(|bounadry| p/n log(n))
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
    { area. clear ();
      boundary. clear ();
      area_root = nullptr;
      area_underRoot = nullptr;
      subPaths. clear ();
      subPathsAbsCriterion = 0;
    }
  
  // Usage:
//set asea, boundary
  void removeIndiscernables ();
  void finish ();
    // Time: O(|area| log(|area|))
//set subPaths
  void finishSubPaths ();
    // Time: O(|boundary| p/n log(n) log(|subPaths|) + |subPaths| log(|area|))
//change topology of tree within area
  Real getImprovement (const DiGraph::Node2Node &boundary2new) const;
    // Time: O(|subPaths| (log(|boundary2new|) + log(|area|)))
  void subPaths2tree ();
    // Update: tree: Paths, absCriterion, Dissim::prediction
    // Time: O(|area| log(|subPaths|) + |subPaths| |area| log(|area|))

  bool large () const
    { return boundary. size () > 64; } // PAR
  bool dense () const
    { const Real density = (Real) area. size () / (Real) boundary. size ();
        // 1..2; 2 <=> sparse
      return density <= 1.2;   // PAR 
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
  VectorPtr<Tree::TreeNode> getPath (const SubPath &subPath) const
    { const Tree::TreeNode* lca_ = nullptr;
      // tree.dissims[subPath.objNum].lca can be used instead of area_root if viaRoot(subPath) and tree topology has not been changed ??
      return Tree::getPath (subPath. node1, subPath. node2, area_root, lca_);
    }
    // Requires: subPath in subPaths
  const Vector<size_t>& boundary2pathObjNums (const DTNode* dtNode) const
    { return dtNode == area_root /* && area_underRoot */
               ? area_underRoot->pathObjNums
               : dtNode        ->pathObjNums;
    }
  const Leaf* getReprLeaf (const DTNode* dtNode) const
    { return dtNode == area_root /* && area_underRoot */
               ? static_cast <const DTNode*> (dtNode->getDifferentChild (area_underRoot)) -> reprLeaf
               : dtNode->reprLeaf;
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
    // Too small values are noise => not stable in bootstrap
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
             && ! from_arg->inDiscernable ()
	           && to_arg
	  	       && to_arg->graph == from_arg->graph
	  	       && to_arg != from_arg
             && ! to_arg->inDiscernable ()
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
	static bool strictlyLess (const Change* a, 
	                          const Change* b);
    // Requires: (bool)a
};



struct Dissim
{
  const Leaf* leaf1 {nullptr};
  const Leaf* leaf2 {nullptr};
    // !nullptr
    // leaf1->name < leaf2->name
  Real target {NAN};
    // Dissimilarity between leaf1 and leaf2; !isNan()
  Real mult {NAN};
  
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
  string getObjName () const;
  void print () const
    { cout << leaf1 << ' ' << leaf2 << ' ' << mult << ' ' << lca << endl; }
  VectorPtr<Tree::TreeNode> getPath () const;
  Real getResidual () const
    { return prediction - target; }
  Real getAbsCriterion (Real prediction_arg) const;
  Real getAbsCriterion () const
    { return getAbsCriterion (prediction); }

  bool operator< (const Dissim &other) const
    { return Pair<const Leaf*> (leaf1, leaf2) < Pair<const Leaf*> (other. leaf1, other. leaf2); }
  bool operator== (const Dissim &other) const
    { return Pair<const Leaf*> (leaf1, leaf2) == Pair<const Leaf*> (other. leaf1, other. leaf2); }
};




///////////////////////////////////////////////////////////////////////////

struct DistTree : Tree
// Of DTNode*
// Least-squares distance tree
// Steiner tree
// nodes.size() >= 2
{
  friend DTNode;
  friend Leaf;
  friend Change;
  friend Subgraph;

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
    
  Vector<Dissim> dissims;
  Real mult_sum {0};
public:
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
public:


  // Input: dissimFName: <dmSuff>-file without <dmSuf>, contains attribute attrName,
  //                     may contain more objects than *this contains leaves
	DistTree (const string &treeFName,
	          const string &dissimFName,
	          const string &attrName,
	          bool sparse);
	  // Input: dissimFName and attrName: may be both empty
	  // Invokes: loadTreeFile(), loadDissimDs(), dissimDs2dissims()
	DistTree (const string &dirName,
	          const string &dissimFName,
	          const string &attrName);
	  // Input: dirName: contains the result of mdsTree.sh; ends with '/'
	  // Invokes: loadTreeDir(), loadDissimDs(), dissimDs2dissims(), setGlobalLen()
	DistTree (const string &dissimFName,
	          const string &attrName,
	          bool sparse);
	  // Invokes: loadDissimDs(), dissimDs2dissims(), neighborJoin()
	DistTree (const string &dataDirName,
            bool loadNewLeaves,
	          bool loadDissim);
	  // Input: dataDirName: ends with '/'
	  //          files: tree, leaf, dissim
	  //        <dataDirName/> contains:
	  //          file name                 line format                                   meaning
	  //          ---------                 -----------------------------                 ----------------------------
	  //          tree                                           
	  //          outlier/<obj>             <obj>                                         Outliers. Their dissimilarities are either in dissim or in dissim.outlier
	  //          dissim                    <obj1> <obj2> <dissimilarity>                 <ob1>, <ob2> are tree leaves
    //          leaf                      <obj_new> <obj1>-<obj2> <leaf_len> <arc_len>
    //         [dissim.add]
    //          new/<obj>                                                               Ne wobjects to add to the tree
    //          search/<obj_new>/                                                       Initialization of search for <obj_new> location
    //        [ search/<obj_new>/dissim   <obj_new> <obj> <dissimilarity>        
    //          search/<obj_new>/leaf     = as in leaf 
    //         [search/<obj_new>/request  <obj_new> <obj>]                              Request to compute dissimilarity
    //        ]
    //          request2dissim.sh         executable with parameters: search/<obj_new>/request search/<obj_new>/dissim.add
    //          strong_outliers           ""|-strong_outliers                           makeDistTree parameter
	  //       ?? request                   <obj1> <obj2>                                 Request to compute dissimilarity
	  //       ?? delete                    <obj>
	  //          version                   <natural number>
	  //          old/{tree,makeDistTree,leaf}.<version>                                  Old versions of data
	  //       <dissimilarity>: >= 0, < INF
	  // Invokes: optimizeSubgraph() for each added Leaf
	  // Time: if loadDissim then O(p log(n) + Time(optimizeSubgraph) * new_leaves)
	  //       if !loadDissim then O(n + new_leaves)
  //  
  explicit DistTree (const string &newickFName);
  DistTree (Prob branchProb,
            size_t leafNum_max);
    // Random tree: DTNode::len = 1
    // Time: O(n)
  DistTree (const DTNode* center,
            uint &areaRadius,
            Subgraph &subgraph,
            Node2Node &newLeaves2boundary);
    // Connected subgraph of center->getTree(); boundary of area are Leaf's of *this
    // If subgraph.dense() then *this has a star topology
    // Update: areaRadius: >= 1; may be decreased
    // Output: subgraph: tree = center->getTree()
    //                   area: contains center, newLeaves2boundary.values(); discernable
    //         newLeaves2boundary
	  // Time: O(log^2(wholeTree.n)) + f(|area|), where wholeTree = center->getDistTree()
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
    // Output: topology, DTNode::len, Leaf::discernable
    // Update: lineNum
  void newick2node (ifstream &f,
                    Steiner* parent);
  void setName2leaf ();
  void loadDissimDs (const string &dissimFName,
                     const string &attrName);
    // Output: dissimDs
    // invokes: dissimDs->setName2objNum()
  // Input: dissimDs
  bool getConnected ();
    // Find connected components of leaves where pairs have dissimilarities with positive multiplicity
    // Return: true <=> 1 connected component
    // Output: DisjointCluster
    //         cout: print other connected components
  size_t setDiscernable ();
    // Return: Number of indiscernable leaves
    // Output: Leaf::len = 0, Leaf::discernable = false, topology
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
  size_t leaves2dissims (Leaf* leaf1,
                         Leaf* leaf2,
                         Real target,
                         Real mult);
    // Append: dissims[]
    // Return: dissims.size() - 1
  bool addDissim (const string &name1,
                  const string &name2,
                  Real dissim);
	  // Return: dissim is added
	  // Update: Dissim, mult_sum, dissim2_sum
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
  const DTNode* lcaName2node (const string &lcaName) const;
    // Return: !nullptr
    // Input: lcaName: <leaf1 name> <objNameSeparator> <leaf2 name>
  size_t extraObjs () const
    { return dissimDs. get () ? dissimDs->objs. size () - name2leaf. size () : 0; }
  size_t getDissimSize_max () const
    { return name2leaf. size () * (name2leaf. size () - 1) / 2; }	
  size_t getSparseDissims_size () const
    { return (size_t) log ((Real) name2leaf. size ()) * 70; }  // PAR
  Set<const DTNode*> getDiscernables () const;
    // Logical leaves
  static void printParam (ostream &os) 
    { os << "PARAMETERS:" << endl;
      os << "Dissimilarity variance: " << varianceTypeNames [varianceType] << endl;
      os << "Max. possible dissimilarity = " << dissim_max () << endl;
      os << "Subgraph radius = " << areaRadius_std << endl;
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
  Real getLeafAbsCriterion () const
  //{ return 2 * absCriterion / (Real) name2leaf. size (); }
    { return absCriterion / mult_sum; }
  Prob getUnexplainedFrac () const
    { return absCriterion / dissim2_sum; }
  Real getErrorDensity () const
    { const Prob r = 1 - getUnexplainedFrac ();   
      return sqrt ((1 - r) / r);   
    }
    // Requires: linear variance of dissimilarities  ??
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
public:
  static Real path2prediction (const VectorPtr<TreeNode> &path);
    // Return: >= 0
	  // Input: DTNode::len
	  // Time: O(|path|)
  void printAbsCriterion_halves () const;
  void setLeafAbsCriterion ();
    // Output: Leaf::absCriterion
    // Time: O(p)
	  
  // Optimization	  
	void optimizeLenArc ();
	  // Return: success
	  // Update: DTNode::len
	  // Output: Dissim::prediction, absCriterion
	  // To be followed by: finishChanges()
	  // Time: O(p log(n)); return = 1.04 * optimal
  void optimizeLenNode ();
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
	bool optimize ();
	  // Update: DTNode::stable
	  // Return: false <=> finished
	  // Requries: getConnected()
	  // Invokes: getBestChange(), applyChanges()
	  // Time of 1 iteration: O(n Time(getBestChange))  
	void optimizeIter (uint iter_max,
	                   const string &output_tree);
	  // Input: iter_max: 0 <=> infinity
	  // Update: cout
	  // Invokes: optimize(), saveFile(output_tree)
	void optimizeSubgraphs (uint areaRadius);
	  // Invokes: optimizeSubgraph()
	  // Time: O(n * Time(optimizeSubgraph))
private:
	uint optimizeSubgraph (const DTNode* center,
	                       uint areaRadius);
	  // Return: adjusted areaRadius
	  // Input: center: may be delete'd
	  // Output: DTNode::stable = true
	  // Invokes: DistTree(center,areaRadius,).optimizeIter()
	  // Time: O(log^2(n) + Time(optimizeIter(),n = min(this->n,2^areaRadius))
  const Change* getBestChange (const DTNode* from);
    // Return: May be nullptr
    // Invokes: tryChange()
    // Time: O(min(n,2^areaRadius_std) log^4(n))
  bool applyChanges (VectorOwn<Change> &changes);
	  // Return: false <=> finished
    // Update: topology, changes (sort by Change::improvement descending)
    // Output: DTNode::stable
    // Invokes: once: finishChanges(), optimizeLen(), optimizeLenLocal()
	void tryChange (Change* ch,
	                const Change* &bestChange);
    // Update: bestChange: positive(improvement)
    // Invokes: Change::{apply(),restore()}
public:	
  // Auxiliary
  void delayDeleteRetainArcs (DTNode* node);
    // Invokes: s->detachChildrenUp()
  size_t finishChanges ();
    // Return: deleteLenZero()
  size_t deleteLenZero ();
    // Delete arcs where len = 0
    // Does not delete root
    // Invokes: deleteLenZero(node), delayDeleteRetainArcs()
private:
  bool deleteLenZero (DTNode* node);
    // Return: success
public:
  void removeLeaf (Leaf* leaf);
    // Invokes: leaf->detachChildrenUp(), optimizeSubgraph(), toDelete.deleteData()
    // Update: detachedLeaves
	  // Time: Time(optimizeSubgraph)
    
  void setReprLeaves ()
    { sort ();
      const_static_cast<DTNode*> (root) -> setRepresentative ();
    }  
    // Output: DTNode::reprLeaf
  Vector<Pair<const Leaf*>> getMissingLeafPairs_ancestors ();
    // Return: not in dissims
    // Invokes: setReprLeaves(), dissims.sort(), DTNode::getSparseLeafPairs()
  Vector<Pair<const Leaf*>> getMissingLeafPairs_subgraphs ();
    // Invokes: setReprLeaves(), findTooLongArcs()
        
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
    // Return: min. length of discernable leaf arcs 
  Real getSqrResidualCorr () const;
    // Return: correlation between squared residual and Dissim::target
    // Input: Dissim::prediction
	  // Time: O(p)
  Real setErrorDensities ();
    // Requires: linear variance of dissimilarities
    // Return: epsilon2_0
    // Input: Dissim::prediction
    // Output: DTNode::{paths,errorDensity}
    // Requires: Leaf::discernable is set
	  // Time: O(p log(n))
  VectorPtr<Leaf> findOutliers (Real &outlier_min) const;
    // Idempotent
    // Output: outlier_min
    // Requires: after setLeafAbsCriterion()    
    // Invokes: Leaf::getRelCriterion(), RealAttr2::normal2outlier()
    // Time: O(n log(n))
  VectorPtr<DTNode> findTooLongArcs (Real &arcLen_outlier_min) const;
    // Output: arcLen_outlier_min

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
{
private:
  const DistTree& tree;
    // !nullptr
  const string dataDir;
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
  };
  Location location;


  struct Leaf2dissim
  {
    const Leaf* leaf {nullptr};
    Real dissim {NAN};
      // Between NewLeaf and leaf
    Real mult {NAN};
      // Function of dissim
      
    // Function of NewLeaf::Location::anchor
    Real dist_hat {NAN};
      // From leaf to NewLeaf::Location::anchor
    bool leafIsBelow {false};
    
    Leaf2dissim (const Leaf* leaf_arg,
                 Real dissim_arg,
                 const DTNode* anchor);
      // Input: anchor = NewLeaf::Location::anchor
    explicit Leaf2dissim (const Leaf* leaf_arg)
      : leaf (leaf_arg)
      {}
      
    Real getDelta () const
      { return dissim - dist_hat; }
    Real getU () const
      { return leafIsBelow ? 1 : -1; }
    Real getEpsilon (const NewLeaf& nl) const
      { const Real leaf_dist_hat = dist_hat + nl. location. arcLen * getU () + nl. location. leafLen;
        return dissim - leaf_dist_hat;
      }
      
    bool operator< (const Leaf2dissim &other) const
      { return leaf < other. leaf; }
    bool operator== (const Leaf2dissim &other) const
      { return leaf == other. leaf; }
  };
  Vector<Leaf2dissim> leaf2dissims;
    // Leaf2dissim::leaf: distinct, sort()'ed
    // Assumption: size() = O(log(n))


  NewLeaf (const DistTree &tree_arg,
           const string &dataDir_arg,
           const string &name_arg,
           bool init);
private:
  string getNameDir      () const  { return dataDir + "/" + name + "/"; }
  string getDissimFName  () const  { return getNameDir () + "dissim"; }
  string getLeafFName    () const  { return getNameDir () + "leaf"; }
  string getRequestFName () const  { return getNameDir () + "request"; }
  void saveLeaf () const
    { OFStream of (getLeafFName ());
      of << name << '\t';
      location. saveText (of);      
      of << endl;
    }
  void saveRequest () const;
    // Input: location.anchor
    // Output: file getRequestFName()
    // Time: O(log^2(n))
  void optimize ();
    // Output: location
    // Update: leaf2dissims.{dist_hat,leafIsBelow}
    // Time: average: O(log^3(n))
    // Invokes: optimizeAnchor()
  void optimizeAnchor (Location &location_best,
                       Vector<Leaf2dissim> &leaf2dissims_best);
    // Update: location, leaf2dissims, location_best
    // Output: leaf2dissims_best
    // Invokes: setLocation(), descend()
  void setLocation ();
    // Output: location.{leafLen,arcLen,absCriterion_leaf}
    // Time: O(log(n))
  bool descend (const DTNode* anchor_new);
    // Return: true <=> location.anchor has leaves in leaf2dissims
    // Update: leaf2dissims
    // Output: location.anchor
    // Time: O(log^2(n))
public:
  void qc () const;
};



}



#endif


