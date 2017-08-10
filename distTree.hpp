// distTree.hpp

#ifndef DISTTREE_HPP
#define DISTTREE_HPP

#include "../cpp/common.hpp"
using namespace Common_sp;
#include "numeric.hpp"
#include "dataset.hpp"
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



constexpr uint areaRadius_std = 5;   // PAR
  // The greater then better DistTree::absCriterion
  // >= 4 <= ChangeToCousin can be applied  



struct DistTree;

//struct DTNode;
  struct Steiner;
  struct Leaf;



typedef  Vector<Real>  Leaf2dist;  // Leaf2dissim ??
  // Index: Leaf::index



struct DTNode : Tree::TreeNode 
{
  friend struct DistTree;
  friend struct Leaf;

  // For optimization
	Real len;
	  // Arc length between *this and *getParent()
	const CompactBoolAttr1* attr {nullptr};
    // In getDistTree().ds
    // ~DTNode() does not delete
    // ExtBoolAttr1: makes faster by 5 % 
  WeightedMeanVar subtreeLen;
    // Global len = len from *this to the leaves
    // weights = sum of DTNode::len in the subtree excluding *this
private:
  bool stable {false};
    // Init: false
public:
  size_t paths {0};
  Real errorDensity {NAN};
protected:
  Vector<bool> subtreeLeaves;  
    // Indexed by Leaf::index
public:
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
  const Leaf* inDiscernable () const;
    // Return: this or nullptr
  bool childrenDiscernable () const
    { return arcs [false]. empty () || ! static_cast <DTNode*> ((*arcs [false]. begin ()) -> node [false]) -> inDiscernable (); }
  Real getHeight () const
    { return subtreeLen. getMean (); }    
private:
  void saveFeatureTree (ostream &os,
                        size_t offset) const;
  void setSubtreeLenUp ();
    // Output: subtreeLen
  void setGlobalLenDown (DTNode* &bestDTNode,
                         Real &bestDTNodeLen_new,
                         WeightedMeanVar &bestGlobalLen);
  void setSubtreeLeaves ();
    // Output: subtreeLeaves
    // Time: O(n^2)
public:
  Real getEpsilon2 () const
    { return (Real) paths * sqr (errorDensity) * len; }
  virtual void setRepresentative () = 0;
    // Output: reprLeaf
    // Requires: after getDistTree().sort()
  virtual const Leaf* selectRepresentative (const Leaf2dist &leaf2dist) const = 0;
    // Return: nullptr or !isNan(leaf2dist [return->index])
    // Depth-first search

  struct Closest : Root
  {
    const DTNode* node;
    Real absCriterion_delta;
    Real leafLen;      
    Real arcLen;
      // From node to node->getParent()
      
    Closest ()
      : node (nullptr)
      , absCriterion_delta (INF)
      , leafLen (NAN)
      , arcLen (NAN)
      {}
    Closest (const DTNode* node_arg,
             Real absCriterion_delta_arg,
             Real leafLen_arg,
             Real arcLen_arg)
      : node (node_arg)
      , absCriterion_delta (absCriterion_delta_arg)
      , leafLen (leafLen_arg)
      , arcLen (arcLen_arg)
      {}
    void qc () const;
      // Invokes: ASSERT(node)
    void saveText (ostream &os) const
      { os <<         (node ? node->getName () : "<none>")
           << "  " << "absCriterion_delta = " << absCriterion_delta 
           << "  " << "leafLen = " << leafLen 
           << "  " << "arcLen = " << arcLen;
      }

    Steiner* insert ();
      // Return: new
      // Update: *node
  };
private:
  void findClosestNode (const Leaf2dist &leaf2dist,
                        Leaf2dist &leaf2hat_dist,
                        Closest &closest) const;
    // Input: subtreeLeaves
    // Update: lead2hat_dist: matches *getParent()
    //         closest
    // Time: O(n^2)
public:
};



struct Steiner : DTNode
// Steiner node
{
  size_t bootstrap [2/*bool*/];  // not used ??
    // 0 - # mismatches, 1 - # matches
  
  
	Steiner (DistTree &tree,
	         Steiner* parent_arg,
	         Real len_arg);
	void qc () const override;


  const Steiner* asSteiner () const final
    { return this; }

  bool isInteriorType () const final
    { return true; }

  void setRepresentative () final
    { reprLeaf = nullptr;
      for (const DiGraph::Arc* arc : arcs [false])
      { DTNode* node = static_cast <DTNode*> (arc->node [false]);
        node->setRepresentative ();
        if (! reprLeaf)
          reprLeaf = node->reprLeaf;
      }
    }
  const Leaf* selectRepresentative (const Leaf2dist &leaf2dist) const final
    { for (const DiGraph::Arc* arc : arcs [false])
        if (const Leaf* leaf = static_cast <const DTNode*> (arc->node [false]) -> selectRepresentative (leaf2dist))
          return leaf;
      return nullptr;
    }

private:
  void reverseParent (const Steiner* target, 
                      Steiner* child);
    // Until target
    // Input: target: !nullptr
    //        child: nullptr <=> *this becomes getTree().root
    // Requires: descendentOf(target)
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
  bool discernable {true}; 
    // false => getParent()->getChildren() is an equivalence class of indiscernables
  Real relLenError {NAN};
private:
  friend struct DistTree;
  friend struct DTNode;
  size_t index {NO_INDEX};
public:
  
  WeightedMeanVar absCriterion;
    // sum = contribution to 2*getTree().absCriterion
  

	Leaf (DistTree &tree,
	      Steiner* parent_arg,
	      Real len_arg,
	      const string &name_arg);
	void qc () const final;
  void saveContent (ostream& os) const final
    { DTNode::saveContent (os);
      if (! isNan (relLenError))
      { ONumber oNum (os, 6, true);  // PAR
        os << "  leaf_error=" << relLenError;
      }
    	if (! discernable)
    	  os << "  non-discernable";
    }


  const Leaf* asLeaf () const final
    { return this; }


  string getName () const final
    { return name; }
  string getNewickName (bool minimal) const final
    { if (minimal)
        return name;
      string s = name + prepend (" ", comment); 
      if (absCriterion. weights)
        s += " " + real2str (getRelLenError (), 1);  // PAR
      return s;
    }
  bool isLeafType () const final
    { return true; }

  void setRepresentative () final
    { reprLeaf = this; }
  const Leaf* selectRepresentative (const Leaf2dist &leaf2dist) const final
    { return isNan (leaf2dist [index]) ? nullptr : this; }

  const DTNode* getDiscernable () const;
    // Return: !nullptr
  Real getLenError () const
    { return sqrt (absCriterion. getMean ()); }
  Real getRelLenError () const;
    // Invokes: getLenError()
  Steiner* collapse (Leaf* other);
    // Return: new; may be nullptr
    // Output: discernable = false
    // Requires: !positive(distance(this,other))
    // Invokes: setParent()
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
  // Update: Topology, DTNode::{attr,len}, tree.prediction
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
	static bool compare (const Change* a, 
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
	  	       && ! to_arg->descendentOf (from_arg)
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
	  	       && to_arg->descendentOf (from_arg)
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
	  	       && ! to_arg->descendentOf (from_arg)
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
	  	       && from_arg->getParent () -> descendentOf (to_arg->getParent ())
	  	       && ! from_arg->getParent () -> descendentOf (to_arg);
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




// DistTree

struct DistTree : Tree
// Of DTNode*
// Least-squares distance tree
// Steiner tree
// nodes.size() >= 2
// For Time: n = # leaves, p = # distances = ds.objs.size()
{
  map<string/*Leaf::name*/,const Leaf*> name2leaf;

  // Dissimilarity
  // May be nullptr
  Common_sp::AutoPtr<Dataset> dissimDs;
    // Original data
  const PositiveAttr2* dissimAttr {nullptr};
    // In *dissimDs

  Dataset ds;
    // objs: pairs of Leaf's
    //       objs[i].name = obj2leaf1[i]->name + "-" + obj2leaf2[i]->name
	  //       objs[i]->mult = dissim2mult(target[i]); positive()
    // attrs: DTNode::attr: 0/1: 1 <=> on the path between the pair of Leaf's    
    // attrs.size() <= 2*n
    // objs.size() <= n*(n-1)/2
    // size = O(n^3)
  Sample dsSample;  
  // !nullptr
  const RealAttr1* target {nullptr};
    // target[i] = dissimilarity between obj2leaf1[i] and obj2leaf2[i]; !isMissing()
  Real dissim2_sum {0};
    // = sum target[i]^2 * mult
  // For optimization
  const RealAttr1* prediction {nullptr};
    // Tree distances
  // For Change
  const CompactBoolAttr1* fromAttr_new {nullptr};
  const CompactBoolAttr1* toAttr_new {nullptr};
  const CompactBoolAttr1* interAttr {nullptr};
  const RealAttr1* target_new {nullptr};
  const RealAttr1* prediction_old {nullptr};
  
  // Tree-ds relations
  size_t attrNum_max {0};
  // !nullptr
  VectorPtr<Leaf> obj2leaf1, obj2leaf2;
    // size() = ds.objs.size()
    // obj2leaf1[i]->name < obj2leaf2[i]->name
    
  Real absCriterion {NAN};
    // = L2LinearNumPrediction::absCriterion  
private:
  Real absCriterion_delta {NAN};
	VectorOwn<DTNode> toDelete;
public:


  // Input: dissimFName: <dmSuff>-file without <dmSuf>, contains attribute attrName
	DistTree (const string &treeFName,
	          const string &dissimFName,
	          const string &attrName,
	          bool sparse);
	  // Input: dissimFName and attrName: may be both empty
	  // Invokes: loadTreeFile(), loadDissimDs(), dissimDs2ds()
	DistTree (const string &dirName,
	          const string &dissimFName,
	          const string &attrName);
	  // Input: dirName: contains the result of mdsTree.sh; ends with '/'
	  // Invokes: loadTreeDir(), loadDissimDs(), dissimDs2ds(), setGlobalLen()
	DistTree (const string &dissimFName,
	          const string &attrName);
	  // Invokes: loadDissimDs(), dissimDs2ds(), neighborJoin()
  //
  explicit DistTree (const string &newickFName);
  DistTree (Prob branchProb,
            size_t leafNum_max);
    // Random tree: DTNode::len = 1
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
    // Output: topology, DTNode::len
    // Update: lineNum
  void setName2leaf ();
  void loadDissimDs (const string &dissimFName,
                     const string &attrName);
    // Output: dissimDs
  void dissimDs2ds (bool sparse);
    // Update: dissimDs: delete
    // Output: ds etc.
    //         Tree: if an object is absent in ds then it is deleted from the Tree
    // Invokes: setReprLeaves(), loadDissimFinish()
  bool addDissim (const string &name1,
                  const string &name2,
                  Real dissim);
	  // Return: dissim is added
	  // Update: ds.objs, dissim2_sum, *target, objLeaf1, objLeaf2
  bool getConnected ();
    // Find connected components of leaves where pairs have dissimilarities with positive multiplicity
    // Return: true <=> 1 connected component
    // Input: dissimDs
    // Output: DisjointCluster
    //         cout: print other connected components
  size_t setDiscernable ();
    // Return: Number of indiscernable leaves
    // Input: dissimDs
    // Output: Leaf::len = 0, Leaf::discernable = false
  void neighborJoin ();
    // Greedy
    // Requires: no missing dissimilarities, star topology
    // Time: O(p n)
  void newick2node (ifstream &f,
                    Steiner* parent);
  void loadDissimFinish ();
    // Output: dsSample, absCriterion_delta
public:
	void qc () const override;
	  // Invokes: ASSERT (eqReal (absCriterion, getAbsCriterion ()))
	void qcAttrs () const;


  void deleteLeaf (TreeNode* leaf,
                   bool deleteTransientAncestor) final;
    // Requires: !optimizable()
    
  static string getObjName (const string &name1,
                            const string &name2);
  const DTNode* lcaName2node (const string &lcaName) const;
    // Return: !nullptr
    // Input: lcaName: <leaf1 name>-<leaf2 name>
  size_t extraObjs () const
    { return dissimDs. get () ? dissimDs->objs. size () - name2leaf. size () : 0; }
  size_t dissimSize_max () const
    { return name2leaf. size () * (name2leaf. size () - 1) / 2; }	
  Set<const DTNode*> getDiscernables () const;
	void printInput (ostream &os) const;
	bool optimizable () const
	  { return ! ds. attrs. empty (); }
	Real getDissim_ave () const
	  { Real average, scatter;
	    target->getAverageScatter (dsSample, average, scatter);
	    return average;
	  }
  Real getLenError () const
    { return sqrt (absCriterion / dsSample. multSum); }
  Prob getUnexplainedFrac () const
    { return absCriterion / dissim2_sum; }
  Real getErrorDensity () const
    { const Prob r = 1 - getUnexplainedFrac ();   
      return sqrt ((1 - r) / r);   
    }
    // Requires: linear variance of dissimilarities  ??
  void reportErrors (ostream &os) const
    { ONumber on (os, 6, false);  // PAR
      os << "absCriterion = " << absCriterion 
         << "  Error density = " << getErrorDensity () * 100 << " %"
         << endl;
    }    
  void saveFeatureTree (const string &fName) const;

private:
#if 0
  friend struct ChangeToChild;
  friend struct ChangeToSibling;
  friend struct Swap;
#endif
  void topology2attrs (const List<DiGraph::Node*>& nodes_arg);
    // Output: DTNode::attrs
	  // Time: O(p (log(n) + |nodes_arg|))
  void clearSubtreeLen ();
    // Invokes: DTNode::subtreeLen.clear()
  void setGlobalLen ();
    // Output: DTNode::subtreeLen, DTNode::len
	  // Time: O(p log(n))
public:
  static Real path2prediction (const VectorPtr<TreeNode> &path);
    // Return: >= 0
	  // Input: DTNode::len
	  // Time: O(path.size()) = O(log(n))
  void setPrediction ();
    // Output: *prediction
    // Invokes: path2prediction()
	  // Time: O(p log(n))
  void checkPrediction () const;
    // Input: ds
  Real getAbsCriterion () const;
    // Input: prediction
	  // Time: O(p)
  void setAbsCriterion ()
    { absCriterion = getAbsCriterion (); }
    // More precise than L2LinearNumPrediction::absCriterion and includes !discernable nodes where dissimilarity != 0  
  void checkAbsCriterion (const string &title) const;
    // Invokes: qc(), checkPrediction(), getAbsCriterion()
  void printAbsCriterion_halves () const;
  void setLeafAbsCriterion ();
    // Output: Leaf::{absCriterion,relLenError}
	  
  // Optimization	  
	bool optimizeLen ();
	  // Return: success
	  // Input: DTNode::attr
	  // Update: DTNode::len
	  // Output: prediction, absCriterion
	  // Invokes: setPrediction(), setAbsCriterion()
	  // To be followed by: finishChanges()
	  // Time: O(p n); 3 min./903 leaves: return = 1.04 * optimal
  void optimizeLenLocal ();
	  // Input: DTNode::attr
	  // Update: DTNode::len
	  // Output: prediction, absCriterion
    // After: deleteLenZero()
    // Postcondition: (*prediction)[] = 0 => (*target)[] = 0 
    // Not idempotent
    // Time: O(n p log(n))
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
	  // Invokes: getBestChange(), applyChanges()
	  // Time of 1 iteration: O(n min(n,2^areaRadius_std) p log(n))  
	void optimizeIter (const string &output_tree);
	  // Update: cout
	  // Invokes: optimize(), saveFile(output_tree)
	void optimizeSubtrees ();
	  // Invokes: optimizeSubtree()
	  // Time: O(n * Time(optimizeSubtree))
	void optimizeAdd (bool sparse,
	                  const string &output_tree);
	  // Input: dissimDs, dissimAttr
	  // Requires: (bool)dissimAttr
	  // Invokes: addDissim(), optimizeSubtree() if leafRelCriterion is large, root->findClosestNode(), DTNode::selectRepresentative()
	  // Time: !sparse: ~30 sec./1 new leaf for 3000 leaves ??
	  //       sparse:    4 sec./1 new leaf for 3500 leaves  ??
private:
  void setSubtreeLeaves ();
    // Output: DTNode::subtreeLeaves, Leaf::index
    // Invokes: DTNode::setSubtreeLeaves()
    // Time: O(n^2)
  void addSubtreeLeaf (Leaf* leaf);
    // Update: DTNode::subtreeLeaves
    // Output: leaf->index
    // Time: O(n)
	Real optimizeSubtree (const Steiner* center);
	  // Return: min. distance to boundary
	  // Input: center: may be delete'd
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
  void finishChanges ();
    // Invokes: deleteLenZero()
  size_t deleteLenZero ();
    // Delete arcs where len = 0
    // Does not delete root
    // Invokes: delayDeleteRetainArcs()
    
  void setReprLeaves ()
    { sort ();
      const_static_cast<DTNode*> (root) -> setRepresentative ();
    }  
    // Output: DTNode::reprLeaf
        
  // After optimization
  void reroot (DTNode* underRoot,
               Real arcLen);
  void reroot ();
    // Molecular clock
    // Invokes: reroot(,)
  void setHeight ()
    { const_static_cast<DTNode*> (root) -> setSubtreeLenUp (); }
    // Input: DTNode::len
    
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
  size_t printLeafRelLenErros (ostream &os,
                               Real relErr_min) const;
    // Return: # outliers
    // Input: relErr_min >= 0
    // Requires: after setLeafAbsCriterion()

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
};



}



#endif


