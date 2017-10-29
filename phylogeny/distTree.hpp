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
constexpr uint areaRadius_std = 5;  
  // The greater then better DistTree::absCriterion
  // >= 4 <= ChangeToCousin can be applied  
constexpr size_t sparsingDepth = 2 * areaRadius_std;  
constexpr Prob rareProb = 0.01; 




struct DistTree;

//struct DTNode;
  struct Steiner;
  struct Leaf;



struct DTNode : Tree::TreeNode 
{
  friend DistTree;
  friend Steiner;
  friend Leaf;

  // For optimization
	Real len;
	  // Arc length between *this and *getParent()
	const CompactBoolAttr1* attr {nullptr};
    // If !nullptr then in getDistTree().ds
    // ~DTNode() does not delete
    // ExtBoolAttr1: makes faster by 5 % 

private:
  bool inPath {false};
  bool stable {false};
    // Init: false
public:
protected:
  WeightedMeanVar subtreeLen; 
    // Temporary
    // Average subtree height 
    // weights = topological ? # leaves : sum of DTNode::len in the subtree excluding *this
public:
  size_t paths {0};
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

  void addAttr ();
    // Output: attr
    // Time: O(p)
  const Leaf* inDiscernable () const;
    // Return: this or nullptr
  bool childrenDiscernable () const
    { return arcs [false]. empty () || ! static_cast <DTNode*> ((*arcs [false]. begin ()) -> node [false]) -> inDiscernable (); }
  Real getHeight_ave () const
    { return subtreeLen. getMean (); }    
    // Requires: after DistTree::setHeight()
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
//void setSubtreeLeaves ();
    // Output: subtreeLeaves
    // Time: O(n^2)
public:
  Real getEpsilon2 () const
    { return (Real) paths * sqr (errorDensity) * len; }
  virtual void setRepresentative () = 0;
    // Output: reprLeaf
    // Requires: after getDistTree().sort()
#if 0
  virtual const Leaf* selectRepresentative (const Leaf2dist &leaf2dist) const = 0;
    // Return: nullptr or !isNan(leaf2dist [return->index])
    // Depth-first search
#endif
  virtual void getDescendants (VectorPtr<DTNode> &descendants,
                               size_t depth) const = 0;
    // Update: descendants
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
#if 0
  const Leaf* selectRepresentative (const Leaf2dist &leaf2dist) const final
    { for (const DiGraph::Arc* arc : arcs [false])
        if (const Leaf* leaf = static_cast <const DTNode*> (arc->node [false]) -> selectRepresentative (leaf2dist))
          return leaf;
      return nullptr;
    }
#endif
  void getDescendants (VectorPtr<DTNode> &descendants,
                       size_t depth) const final;

private:
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
  
  Real absCriterion {NAN};
    // sum = contribution to 2*getTree().absCriterion
  Real relCriterion {NAN};
    // !isNan() => = getRelCriterion()
  

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
      { const ONumber oNum (os, 6, true);  // PAR
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

  void setSubtreeLenUp (bool topological) final
    { subtreeLen. clear ();
      if (topological)
    	  subtreeLen. add (0, 1);
    }
  void setRepresentative () final
    { reprLeaf = this; }
#if 0
  const Leaf* selectRepresentative (const Leaf2dist &leaf2dist) const final
    { return isNan (leaf2dist [index]) ? nullptr : this; }
#endif
  void getDescendants (VectorPtr<DTNode> &descendants,
                       size_t /*depth*/) const final
    { descendants << this; }

  const DTNode* getDiscernable () const;
    // Return: !nullptr
  Real getRelCriterion () const;
    // Return: average over all Leaf's = 1
    //         If getDistTree().ds.objs has no all pairs of Leaf's then approximate
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




//struct Change
  //struct Move
    //struct ChangeToSibling
    //struct ChangeToChild
  //struct Swap
    //struct ChangeToUncle
    //struct ChangeToCousin



struct Change : Root  
// Of topology
// Input: tree.dsSample
{
protected:
	const DistTree& tree;
public:
	const DTNode* from;
	  // !nullptr
	Real improvement {NAN};
	  // isNan() or positive()
    // Too small values are noise => not stable in bootstrap
	// Output of apply_()
protected:
	Real fromLen {NAN};
	const DTNode* to;
	  // !nullptr
	Real toLen {NAN};
public:
	VectorPtr<Tree::TreeNode> targets;  
	  // DTNode's whose len may be changed 
	//
private:
  enum Status {eInit, eApplied, eDone};
  Status status;
public:
	
	
protected:
	Change (const DTNode* from_arg,
				  const DTNode* to_arg)
		: tree (const_cast <DistTree&> (from_arg->getDistTree ()))
		, from (from_arg)
		, to (to_arg)
		, targets (4, nullptr)  
		, status (eInit)
		{ targets. clear ();
		  targets << from << to; 
		}
    // Requires: valid_()
	static bool valid_ (const DTNode* from_arg,
	                    const DTNode* to_arg)
	  { return    from_arg
             && from_arg->graph
             && ! from_arg->inDiscernable ()
	           && to_arg
	  	       && to_arg->graph == from_arg->graph
	  	       && to_arg != from_arg
             && ! to_arg->inDiscernable ();
	  }
    // Requires: parameters are the same as in the constructor
public:
 ~Change ();
	void qc () const override;
	  // Invokes: valid()
protected:
	void saveText (ostream& os) const override
	  { os << from->getName () << " (parent = " << (from->getParent () ? from->getParent () -> getName () : "null") << ") -> " << to->getName () 
	       << "  " << type () 
         << "  improvement = " << improvement; 
	  }
public:
	void print (ostream& os) const override
	  { saveText (os);
	    os << endl; 
	  }


	virtual const char* type () const = 0;
	virtual bool valid () const = 0;
  // Update: Topology, DTNode::len, tree.prediction
	bool apply ();
	  // Return: success
	  // Minimum change to compute tree.absCriterion
	  // Output: *tree.prediction_old
	  // status: eInit --> eApplied|eFail
	  // Invokes: apply_()
	void restore ();
	  // Output: *tree.prediction
	  // status: eApplied --> eInit
	  // Invokes: restore_()
	void commit ();
	  // status: eApplied --> eDone
	  // Invokes: commit_()
private:
  // Update: Tree::tree
  // Output: DTNode::len, *tree.prediction
	virtual bool apply_ () = 0;
	  // Return: success  
	  // Output: *tree.prediction
	virtual void restore_ () = 0;
	virtual void commit_ () = 0;
    // May invoke: tree.delayDeleteRetainArcs()
public:
	static bool strictlyLess (const Change* a, 
	                          const Change* b);
    // Requires: (bool)a
};



struct Move : Change
{
protected:
	// !nullptr
	Steiner* oldParent {nullptr};
	  // Old from->getParent()
	Steiner* arcEnd {nullptr};
	  // Old to->getParent()
	Steiner* inter {nullptr};
	  // Between *to and *arcEnd
public:


protected:	
	Move (const DTNode* from_arg,
			  const DTNode* to_arg)
		: Change (from_arg, to_arg)
		{}
};



struct ChangeToSibling : Move
// *to becomes a sibling of *from
// This Change is enough to transform any topology to any topology. Proof: by induction by node depth descending
{
	ChangeToSibling (const DTNode* from_arg,
					         const DTNode* to_arg)
		: Move (from_arg, to_arg)
		{}
	static bool valid_ (const DTNode* from_arg,
	                    const DTNode* to_arg)
	  { return    Move::valid_ (from_arg, to_arg)
    	       && from_arg->getParent ()
	  	       && ! to_arg->descendantOf (from_arg)
	  	       && ! (from_arg->getParent() == to_arg->getParent() && from_arg->getParent() -> arcs [false]. size () == 2)  
	  	       && ! (from_arg->getParent() == to_arg              && from_arg->getParent() -> arcs [false]. size () == 2); 
	  }


	const char* type () const final
	  { return "sibling"; }
	bool valid () const final
	  { return valid_ (from, to); }
private:
	bool apply_ () final;
	  // Time: O(p log(n))
	void restore_ () final;
	void commit_ () final;
};



#if 0   // Error density: increase by 0.004; time: decrease by 30 %
struct ChangeToChild : Move
// New graph: *to -> *inter -> *oldParent
// Similarity: ChangeToSibling: from ->setParent(inter)
//             ChangetoChild:   inter->setParent(from->getParent())
{
	ChangeToChild (const DTNode* from_arg,
					       const DTNode* to_arg)
		: Move (from_arg, to_arg)
		{}
	static bool valid_ (const DTNode* from_arg,
	                    const DTNode* to_arg)
	  { return    Move::valid_ (from_arg, to_arg)
	           && from_arg->asSteiner ()
	  	       && to_arg->descendantOf (from_arg)
	  	     //&& ! (to_arg->getParent() == from_arg && from_arg -> arcs [false]. size () == 2)
	  	       && to_arg->getParent () != from_arg;  // *this is redundant given ChildToSibling: from_arg->getParent()=to_arg
	  }


	const char* type () const final
	  { return "child"; }
	bool valid () const final
	  { return valid_ (from, to); }
private:
	bool apply_ () final;
	void restore_ () final;
	void commit_ () final;
};




struct Swap : Change
// Swap *from and *to
{
protected:
	Swap (const DTNode* from_arg,
			  const DTNode* to_arg)
		: Change (from_arg, to_arg)
		{}
	static bool valid_ (const DTNode* from_arg,
	                    const DTNode* to_arg)
	  { return    Change::valid_ (from_arg, to_arg)
    	       && from_arg->getParent ()
	  	       && ! to_arg->descendantOf (from_arg)
	  	       && from_arg->getParent () != to_arg->getParent ();
	  }


private:
	bool apply_ () final;
	void restore_ () final;
	void commit_ () final
	  {}
};



struct ChangeToUncle : Swap
// *to = uncle
// Includes the Felsenstein "swap" change
{
	ChangeToUncle (const DTNode* from_arg,
					       const DTNode* to_arg)
		: Swap (from_arg, to_arg)
		{}
	static bool valid_ (const DTNode* from_arg,
	                    const DTNode* to_arg)
	  { return    Swap::valid_ (from_arg, to_arg)
	  	       && from_arg->getParent () -> descendantOf (to_arg->getParent ())
	  	       && ! from_arg->getParent () -> descendantOf (to_arg);
	  }


	const char* type () const final
	  { return "uncle"; }
	bool valid () const final
	  { return valid_ (from, to); }
};



struct ChangeToCousin : Swap
// *to = first cousin
{
	ChangeToCousin (const DTNode* from_arg,
  					      const DTNode* to_arg)
		: Swap (from_arg, to_arg)
		{}
	static bool valid_ (const DTNode* from_arg,
	                    const DTNode* to_arg)
	  { return    Swap::valid_ (from_arg, to_arg)
	  	       && to_arg->getParent ()  
	  	       && from_arg->getParent () -> getParent () == to_arg->getParent () -> getParent ()
	  	       && from_arg < to_arg;
	  }


	const char* type () const final
	  { return "cousin"; }
	bool valid () const final
	  { return valid_ (from, to); }
};
#endif




///////////////////////////////////////////////////////////////////////////

struct DistTree : Tree
// Of DTNode*
// Least-squares distance tree
// Steiner tree
// nodes.size() >= 2
// For Time: n = # leaves, p = # distances = ds.objs.size()
//           p >= n
{
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

  friend Change;
  friend ChangeToSibling;
  Dataset ds;
    // objs: pairs of Leaf's
    //       objs[i].name = obj2leaf1[i]->name + objNameSeparator + obj2leaf2[i]->name
	  //       objs[i]->mult = dissim2mult(target[i]); positive()
    // attrs: DTNode::attr: 0/1: 1 <=> on the path between the pair of Leaf's    
    // attrs.size() <= 2*n
    // objs.size() <= n*(n-1)/2
    // size = O(p)
  Sample dsSample;  
    // = Sample(ds)
  // !nullptr
  const RealAttr1* target {nullptr};
    // target[i] = dissimilarity between obj2leaf1[i] and obj2leaf2[i]; !isMissing()
public:
  Real dissim2_sum {0};
    // = sum target[i]^2 * mult
private:
  // For optimization
  const RealAttr1* prediction {nullptr};
    // Tree distances
  // For Change
  const CompactBoolAttr1* fromAttr_new {nullptr};
  const CompactBoolAttr1* toAttr_new {nullptr};
  const CompactBoolAttr1* interAttr {nullptr};
  const RealAttr1* target_new {nullptr};
  const RealAttr1* prediction_old {nullptr};
  bool nodeAttrExist {false};
  
  // ds-tree relations
  // !nullptr
  VectorPtr<Leaf> obj2leaf1, obj2leaf2;
    // size() = ds.objs.size()
    // obj2leaf1[i]->name < obj2leaf2[i]->name
public:
    
  Real absCriterion {NAN};
    // = L2LinearNumPrediction::absCriterion  
private:
  friend DTNode;
  friend Leaf;
  size_t attrNum_max {0};
  Real absCriterion_delta {NAN};
	VectorOwn<DTNode> toDelete;
	VectorOwn<Leaf> detachedLeaves;
public:


  // Input: dissimFName: <dmSuff>-file without <dmSuf>, contains attribute attrName,
  //                     may contain more objects than *this contains leaves
	DistTree (const string &treeFName,
	          const string &dissimFName,
	          const string &attrName,
	          bool sparse);
	  // Input: dissimFName and attrName: may be both empty
	  // Invokes: loadTreeFile(), loadDissimDs(), dissimDs2ds(), topology2attrs()
	DistTree (const string &dirName,
	          const string &dissimFName,
	          const string &attrName);
	  // Input: dirName: contains the result of mdsTree.sh; ends with '/'
	  // Invokes: loadTreeDir(), loadDissimDs(), dissimDs2ds(), setGlobalLen(), topology2attrs()
	DistTree (const string &dissimFName,
	          const string &attrName,
	          bool sparse);
	  // Invokes: loadDissimDs(), dissimDs2ds(), neighborJoin(), topology2attrs()
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
	  //       ?? request                   <obj1> <obj2>                                 Request to compute dissimilarity
	  //       ?? deleted                   <obj>
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
            uint areaRadius,
            VectorPtr<TreeNode> &area,
            const DTNode* &area_root,
            Node2Node &newLeaves2boundary);
    // Connected subgraph of center->getTree(); boundary of area are Leaf's of *this
    // Input: areaRadius: >= 1
    // Output: area: contains center, area_root, newLeaves2boundary.values(); getTree() = center->getTree()
    //               discernable
    //         area_root: !nullptr
    //         newLeaves2boundary
	  // Time: O(wholeDs.p log(wholeDs.n)) + f(|area|), where wholeDs = center->getDistTree().ds
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
  void dissimDs2ds (bool sparse);
    // Update: dissimDs: delete
    // Output: ds etc.
    //         Tree: if an object is absent in ds then it is deleted from the Tree
    // Invokes: getSelectedPairs(), loadDissimFinish()
  void loadDissimPrepare (size_t pairs_max,
                          streamsize target_decimals);
    // Output: target
  bool addDissim (const string &name1,
                  const string &name2,
                  Real dissim);
	  // Return: dissim is added
	  // Update: ds.objs, dissim2_sum, *target, objLeaf1, objLeaf2
  void loadDissimFinish ();
    // Output: dsSample, absCriterion_delta
    // Invokes: addAttr()
    // Time: O(p n)
public:
	void qc () const override;
	  // Invokes: ASSERT (eqReal (absCriterion, getAbsCriterion ()))


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
	  { return ! ds. objs. empty (); }
	Real getDissim_ave () const
	  { Real average, scatter;
	    target->getAverageScatter (dsSample, average, scatter);
	    return average;
	  }
  Real getLeafAbsCriterion () const
    { return 2 * absCriterion / (Real) name2leaf. size (); }
  Prob getUnexplainedFrac () const
    { return absCriterion / dissim2_sum; }
  Real getErrorDensity () const
    { const Prob r = 1 - getUnexplainedFrac ();   
      return sqrt ((1 - r) / r);   
    }
    // Requires: linear variance of dissimilarities  ??
  void reportErrors (ostream &os) const
    { const ONumber on (os, 6, false);  // PAR
      os << "absCriterion = " << absCriterion 
         << "  Error density = " << getErrorDensity () * 100 << " %"
         << endl;
    }    
  void saveFeatureTree (const string &fName) const;

private:
#if 0
  friend ChangeToChild;
  friend ChangeToSibling;
  friend Swap;
#endif
  void resetAttrs ();
    // Time: O(p n)
  void topology2attrs ();
    // Output: DTNode::attr
	  // Time: O(p log(n))
  void topology2attrs (const List<DiGraph::Node*>& nodes_arg);
    // Input: nodes_arg: subset of nodes
    // Output: DTNode::attr
    // Temporary: DTNode::inPath
	  // Time: O(p (log(n) + |nodes_arg|))
  void clearSubtreeLen ();
    // Invokes: DTNode::subtreeLen.clear()
public:
  static Real path2prediction (const VectorPtr<TreeNode> &path);
    // Return: >= 0
	  // Input: DTNode::len
	  // Time: O(path.size()) = O(log(n))
  void setPrediction ();
    // Output: *prediction
    // Invokes: path2prediction()
	  // Time: O(p log(n))
  Real getAbsCriterion () const;
    // Input: prediction
	  // Time: O(p)
  void setAbsCriterion ()
    { absCriterion = getAbsCriterion (); }
    // More precise than L2LinearNumPrediction::absCriterion and includes !discernable nodes where dissimilarity != 0  
  void checkAbsCriterion (const string &title) const;
    // Invokes: qc(), getAbsCriterion()
  void printAbsCriterion_halves () const;
  void setLeafAbsCriterion ();
    // Output: Leaf::absCriterion
	  
  // Optimization	  
  // Input: DTNode::attr
private:
  void getSkipRetain (DTNode* &toSkip,
                      DTNode* &toRetain);
    // Output: toSkip, toRetain; may be nullptr
    // Update: {toSkip,toRetain}->len
    // toSkip->attr is redundant
    // toSkip->getparent() = toRetain->getParent() = root: the only children
public:
	bool optimizeLenArc ();
	  // Return: success
	  // Update: DTNode::len
	  // Output: prediction, absCriterion
	  // Invokes: setPrediction(), setAbsCriterion()
	  // To be followed by: finishChanges()
	  // Time: O(p n); return = 1.04 * optimal
  void optimizeLenNode ();
	  // Update: DTNode::len
	  // Output: prediction, absCriterion
    // After: deleteLenZero()
    // Postcondition: (*prediction)[] = 0 => (*target)[] = 0 
    // Not idempotent
    // Time: O(n p log(n))
  void checkAttrPrediction () const;
    // Input: ds, DTNode::attr
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
	  // Time of 1 iteration: O(n min(n,2^areaRadius_std) p log(n))  
	void optimizeIter (const string &output_tree);
	  // Update: cout
	  // Invokes: optimize(), saveFile(output_tree)
	void optimizeSubgraphs ();
	  // Invokes: optimizeSubgraph()
	  // Time: O(n * Time(optimizeSubgraph))
private:
//void setSubtreeLeaves ();
    // Output: DTNode::subtreeLeaves, Leaf::index
    // Invokes: DTNode::setSubtreeLeaves()
    // Time: O(n^2)
#if 0
  void addSubtreeLeaf (Leaf* leaf);
    // Update: DTNode::subtreeLeaves
    // Output: leaf->index
    // Time: O(n)
#endif
	Real optimizeSubgraph (const Steiner* center);
	  // Return: min. distance to boundary
	  // Input: center: may be delete'd
	  // Update: DTNode::attr
	  // Output: DTNode::stable = true
	  // Invokes: DistTree(center,areaRadius_std,).optimizeIter(), setAbsCriterion()
	  // Time: O(p (log(n) + min(n,2^areaRadius_std)) + Time(optimizeIter,n = min(this->n,2^areaRadius_std)))
  const Change* getBestChange (const DTNode* from);
    // Return: May be nullptr
    // Invokes: tryChange()
    // Time: O(min(n,2^areaRadius_std) p log(n))
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
  void delayDeleteRetainArcs (DTNode* s);
    // Invokes: s->isolateChildrenUp()
  size_t finishChanges ();
    // Return: deleteLenZero()
  size_t deleteLenZero ();
    // Delete arcs where len = 0
    // Does not delete root
    // Invokes: delayDeleteRetainArcs()
  void removeLeaf (Leaf* leaf);
    // Invokes: leaf->isolateChildrenUp(), optimizeSubgraph(), toDelete.deleteData()
    // Update: detachedLeaves
	  // Time: Time(optimizeSubgraph)
    
  void setReprLeaves ()
    { sort ();
      const_static_cast<DTNode*> (root) -> setRepresentative ();
    }  
    // Output: DTNode::reprLeaf
  Set<string> selectPairs ();  
    // Return: pairs of Leaf->name's 
    //         O(log(n)) pairs for each Leaf
    //         not in ds.objs
    //         reroot(true) reduces size()
    // Invokes: setReprLeaves(), ds.setName2objNum(), getObjName()
        
  // After optimization
  void setHeight ()
    { const_static_cast<DTNode*> (root) -> setSubtreeLenUp (false); }
    // Input: DTNode::len
    // Output: DTNode::subtreeLen
  void reroot (DTNode* underRoot,
               Real arcLen);
  Real reroot (bool topological);
    // Center of the tree w.r.t. DTNode::setGlobalLenDown()
    // Return: root->getHeight()
    // Invokes: setGloballenDown(), reroot(,)
    
  // Quality
  Real getMeanResidual () const;
    // Input: prediction
	  // Time: O(p)
  Real getMinLeafLen () const;
    // Return: min. length of discernable leaf arcs 
  Real getSqrResidualCorr () const;
    // Return: correlation between squared residual and target
    // Input: prediction
	  // Time: O(p)
  Real setErrorDensities ();
    // Requires: linear variance of dissimilarities
    // Return: epsilon2_0
    // Input: prediction
    // Output: DTNode::{paths,errorDensity}
    // Requires: Leaf::discernable is set
	  // Time: O(p log(n))
  VectorPtr<Leaf> findOutliers (bool strong,
                                Real &outlier_min) const;
    // Idempotent
    // Input: strong: use log(Leaf::getRelCriterion())
    // Output: outlier_min
    // Requires: after setLeafAbsCriterion()    
    // Invokes: Leaf::getRelCriterion(), RealAttr2::normal2outlier()

  // Statistics
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
    // Leaf2dissim::leaf: distinct, ordered
    // Assumption: size() = O(log(n))


  NewLeaf (const DistTree &tree_arg,
           const string &dataDir_arg,
           const string &name_arg,
           bool init);
private:
  string getNameDir () const       { return dataDir + "/" + name + "/"; }
  string getDissimFName () const   { return getNameDir () + "dissim"; }
  string getLeafFName ( )const     { return getNameDir () + "leaf"; }
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


