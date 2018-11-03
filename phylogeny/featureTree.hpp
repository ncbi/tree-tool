// featureTree.hpp

#ifndef FEATURE_TREE_HPP
#define FEATURE_TREE_HPP

#include "../common.hpp"
#include "../graph.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
#include "../dm/dataset.hpp"
using namespace DM_sp;



namespace FeatureTree_sp
{



struct Phyl;
  struct Species;
    struct Fossil;
    struct Strain;
  struct Genome;



struct Feature : Named
// Boolean attribute of Genome
// Gene or phenotype
// name: key
{
  typedef string Id;
  // For FeatureTree::len
  // Valid if !allTimeZero
//Prob lambda_0;  ??
  // Stats
  size_t genomes {0};
  VectorPtr<Phyl> gains;
  VectorPtr<Phyl> losses;


  Feature ()
    {}
	explicit Feature (const Id &name_arg);
	void qc () const override;
	void saveText (ostream& os) const override
	  { os << name << " +" << gains. size () << " -" << losses. size () << " / " << genomes << endl; }

	
  bool operator== (const Feature &other) const
    { return name == other. name; }
  bool operator< (const Feature &other) const
    { return name < other. name; }
  static bool statEqual (const Feature& a,
                         const Feature& b)
    { return    a. gains  == b. gains
    	       && a. losses == b. losses;
    }
  static bool statLess (const Feature& a,
                        const Feature& b);
  size_t mutations () const
    { return gains. size () + losses. size (); }
  bool monophyletic () const
    { return gains. size () <= 1 && losses. empty (); }
  bool better (const Feature* other) const
    { return    ! other 
             || mutations () < other->mutations ()
             || (mutations () == other->mutations () && name.size () > other->name. size ());
    }
};



inline bool eqTreeLen (Real len1,
                       Real len2)
  { return eqReal (len1, len2, 1e-4 /*PAR*/); }  



struct FeatureTree;



struct Phyl : Tree::TreeNode 
{
  // For FeatureTree::len
	struct CoreEval
	{
		Real treeLen {0};
		  // >= 0
		ebool core {EFALSE};
		
    CoreEval ()
      {}
		CoreEval (Real treeLen_arg,
		          ebool core_arg)
		  : treeLen (treeLen_arg)
		  , core (core_arg)
		  {}
		  
		bool operator== (const CoreEval &other) const
		  { return    eqReal (treeLen, other. treeLen)
		  	       && core == other. core;
		  }
	};
	Vector<CoreEval> parent2core [2/*bool parentCore*/];
	  // CoreEval::core: optimal given parentCore
    // size() = getFeatureTree().features.size()
	Real weight [2/*thisCore*/] [2/*parentCore*/];
	  // = -log(prob); >= 0; may be INF
	Vector<bool> core;
    // size() = getFeatureTree().features.size()
	size_t index_init;
	  // Matches the orginial node number in DFS
private:
  bool stable {false};
  friend struct FeatureTree;
public:


protected:
	Phyl (FeatureTree &tree,
        Species* parent_arg);
	  // To be followed by setWeight()
public:
	void init ();
  void qc () const override;
protected:
  void saveContent (ostream& os) const override;
    // Input: core
public:
#if 0
private:
  bool getSaveSubtreeP () const;
public:
#endif


  virtual const Species* asSpecies () const
    { return nullptr; }
  virtual const Fossil* asFossil () const
    { return nullptr; }
  virtual const Strain* asStrain () const
    { return nullptr; }
  virtual const Genome* asGenome () const
    { return nullptr; }


  const FeatureTree& getFeatureTree () const;
  
  // weight[][]
  virtual void setWeight () = 0;
    // Output: weight[][]
	Real feature2weight (size_t featureIndex,
	                     bool thisCore,
	                     bool parentCore) const;
	  // Return: >= 0	  
	Real feature2weight (bool thisCore,
	                     bool parentCore) const
	  { return feature2weight (NO_INDEX/*??*/, thisCore, parentCore); }
	  
	virtual void setCore () = 0;
	  // Input: getParent()->core	
	  // Output: core
	// Input: core
	bool feature2parentCore (size_t featureIndex) const;
protected:
	bool feature2core (size_t featureIndex) const;
public:
	size_t getCoreSize () const;  
  size_t getCoreChange (bool gain) const;
    // !gain <=> loss
  size_t getCoreChange () const
    { return   getCoreChange (false)
    	       + getCoreChange (true);
    }
  Real getDistance () const;
    // Distance to getParent() for features in getFeatureTree().features
    // Return: >= 0
    // Input: core
    // Invokes: getPooledDistance()
	Real getSubtreeLength () const;  
	  // Return: sum_{n \in subtree} n->getDistance() 
#if 0
	Real getFeatureDistance (size_t featureIndex) const;
	Real getFeatureSubtreeLength (size_t featureIndex) const;  
	Real getFeatureTreeLen (size_t featureIndex) const
	  { return parent2core [feature2parentCore (featureIndex)] [featureIndex]. treeLen; }
#endif
  Real getNeighborDistance () const;
    // Invokes: getDistance()
  bool badNeighborDistance (Real &neighborDistance_stnd,
                            Real &depth_stnd) const;
    // Output: neighborDistance_stnd, depth_stnd
  void getBadNodes (VectorPtr<Phyl> &badNodes,
                    bool parentBad) const;
    // Update: badNodes (append)
    // DFS

  // parent2core[]
  // Sankoff algorithm
protected:
	virtual void setCoreEval (size_t featureIndex,
        	                  bool parentCore,
        	                  CoreEval ce) = 0;
	virtual void assignFeature (size_t featureIndex);
	  // Output: parent2core[]
	  // Input: children->parent2core[]
	virtual void assignFeatures ();
	  // Invokes: assignFeature()
public:
	void assignFeaturesDown ();
	  // Post-order DFS
	  // Invokes: assignFeatures()

  virtual void getParent2corePooled (size_t parent2corePooled [2/*thisCore*/] [2/*parentCore*/]) const;
    // Output: parent2corePooled[][]: core change for features not in getFeatureTree().features
	Real getPooledDistance () const;
    // Distance to getParent() for features not in getFeatureTree().features
    // Invokes: getParent2corePooled()
};



struct Species : Phyl
// core: Function of getFeatureTree().superRootCore[] and Phyl::parent2core[]
{
  string id;
    // !empty() => SPECIES.id
  // Valid if !allTimeZero
	Real time;
	  // isNan() || >= 0
    // Exponential distribution, if yes then use P_{prior} ??
//Real timeVariance;  // Read the MLE theory ??
	Real pooledSubtreeDistance {INF};  
	  // Part of the arc (this,getParent()) length due to FeatureTree::commonCore and Genome::singletons
	  // >= 0
private:
	Real weight_old [2/*thisCore*/] [2/*parentCore*/];
  Real time_old {NAN};
    // May be NAN
	Real pooledSubtreeDistance_old {NAN};  
  struct Movement
  {
  	bool parentCore;
  	size_t featureIndex;
    CoreEval from;  
  	Movement (bool parentCore_arg,
  	          size_t featureIndex_arg,
					    CoreEval from_arg)
			: parentCore (parentCore_arg)
			, featureIndex (featureIndex_arg)
			, from (from_arg)
			{}
	  void undo (Species* a) const;
  };
  bool movementsOn {false};
  Vector<Movement> movements;
    // Valid if movementsOn
public:


protected:
	Species (FeatureTree &tree,
	         Fossil* parent_arg,
	         const string &id_arg,
	         Real time_arg);
	  // To be followed by: init()
public:
  void qc () const override;
  void saveContent (ostream& os) const override;


  const Species* asSpecies () const final
    { return this; }


  string getName () const override
    { return id. empty () ? Phyl::getName () : id; }
	double getParentDistance () const final;
  void setWeight () final;
    // Input: time, getFeatureTree().lambda0
	void setCore () final;
private:
	void setCoreEval (size_t featureIndex,
	                  bool parentCore,
	                  CoreEval ce) final
	  {	CoreEval& ce_old = parent2core [parentCore] [featureIndex];
		  if (ce_old == ce)
		  	return;
		  if (movementsOn)
			  movements << Movement (parentCore, featureIndex, ce_old);
	    ce_old = ce;
	  }
public:
protected:
  void assignFeatures () override
    { Phyl::assignFeatures ();
      pooledSubtreeDistance = getPooledSubtreeDistance ();
    }
public:

private:
  Real getTime () const;
    // Input: core, getFeatureTree().{lambda0,timeOptimFrac}
  Real getTime_max () const;  // ??
    // Should be: weight[!c][c] > weight[c][c], otherwise the data are too few to create *this
    // qc(): time <= getTime_max()
public:
  void setTimeWeight ()
    { time = getTime (); 
    	setWeight ();
    }
	// Usage: assignTime() ... {restoreTime()|commitTime()}
	void assignTime ();
	  // Invokes: setTimeWeight()
	void restoreTime ();
	void commitTime ()
	  { time_old = NAN; }

	virtual Real getPooledSubtreeDistance () const = 0;
	  
	// Usage: assignFeaturesUp(a) ... {restoreFeaturesUp(a)|commitFeaturesUp(a))
  // Input: toParentExcluding: may be 0
	void assignFeaturesUp (const Fossil* toParentExcluding);
	  // Invokes: assignFeatures()
	void restoreFeaturesUp (const Fossil* toParentExcluding);
	void commitFeaturesUp (const Fossil* toParentExcluding);
	// Usage: {rememberFeatures()|rememberAllFeatures()} ... {restoreFeatures()|commitFeatures()}
	void rememberAllFeatures ();
protected:
	virtual void rememberFeatures ();
	virtual void restoreFeatures ();
	void commitFeatures ();
public:
};



struct Fossil : Species
{
	Fossil (FeatureTree &tree,
	        Fossil* parent_arg,
	        const string &id_arg,
	        Real time_arg)
		: Species (tree, parent_arg, id_arg, time_arg)
		{}
	void qc () const override;


  const Fossil* asFossil () const final
    { return this; }

  bool isInteriorType () const final
    { return true; }

  Real getPooledSubtreeDistance () const final;

  void setId (uint &id_arg);
    // Output: id
    // Update: id_arg
};



struct Strain : Species
// Distance from children to *this depends on Genome quality
{
  bool singletonsInCore {true};
    // getGenome()->singletons: true <=> in this->core, false <=> annotation errors in Genome
    // Init: true
    // Analog of parent2core[false] and core
private:
  bool singletonsInCore_old {true};
public:


	Strain (FeatureTree &tree,
	        Fossil* parent_arg,
	        const string &id_arg,
	        Real time_arg)
		: Species (tree, parent_arg, id_arg, time_arg)
		{}
  void qc () const override;
  void saveContent (ostream& os) const override;


  const Strain* asStrain () const final
    { return this; }


  string getName () const final
    { return "s" + id; }
  string getNewickName (bool /*minimal*/) const final
    { return string (); }

  void getParent2corePooled (size_t parent2corePooled [2/*thisCore*/] [2/*parentCore*/]) const final;
	void assignFeatures () final;
	Real getPooledSubtreeDistance () const final;
private:
	void rememberFeatures () final
	  { Species::rememberFeatures ();
	    singletonsInCore_old = singletonsInCore;
	  }
	void restoreFeatures () final
	  { singletonsInCore = singletonsInCore_old;
	    Species::restoreFeatures ();
	  }
public:

  const Genome* getGenome () const
    { return static_cast <const Phyl*> (static_cast <const DiGraph::Arc*> (arcs [false]. front ()) -> node [false]) -> asGenome (); }
    // Return: !nullptr
  size_t orphans () const;
};



struct Genome : Phyl
{
  string id;  
    // GENOME.id
    // !empty()
  
  Set<string> nominals;
    // Subset of getFeatureTree().nominals
	Vector<bool> optionalCore;
    // size() = getFeatureTree().features.size()
private:
  map<Feature::Id,bool/*optional*/> coreSet; 
  friend struct FeatureTree;
public:
  size_t coreNonSingletons {0};
    // Includes optionalCore[]  
  Set<Feature::Id> singletons;
//Set<Feature::Id> missedCore;  ??
    // Opposite to singletons


	Genome (FeatureTree &tree,
	        Strain* parent_arg,
	        const string &id_arg);
	  // To be followed by: initDir(), init()
  static string geneLineFormat ()
    { return "{{<gene> [<optional (0|1)>]} | {<nominal name>:<value>} \\n}*"; }
	void initDir (const string &geneDir);
	  // Input: file "geneDir/id" with the format: `geneLineFormat()`
	  // Output: coreSet, coreNonSingletons
	void coreSet2nominals ();
	  // Update: getFeatureTree().nominals, nominals
	void nominals2coreSet ();
	  // Update: coreSet
	void init (const map <Feature::Id, size_t/*index*/> &feature2index);
	  // Output: CoreEval::core, core
	void qc () const override;
	void saveContent (ostream& os) const override;


  const Genome* asGenome () const final
    { return this; }

  string getName () const final
    { return "g" + id; }
  string getLeafName () const final
    { return id; }
	double getParentDistance () const final
	  { return 0; }
  string getNewickName (bool /*minimal*/) const final
	  { return id /*+ " " + taxName + getNameExtra ()*/; 
	  }
  bool isLeafType () const final
    { return true; }

  void setWeight () final;
  void getParent2corePooled (size_t parent2corePooled [2/*thisCore*/] [2/*parentCore*/]) const final;
	void setCore () final;
private:
	void setCoreEval (size_t featureIndex,
	                  bool parentCore,
	                  CoreEval ce) final
	  {	parent2core [parentCore] [featureIndex] = ce; }
	void assignFeature (size_t featureIndex) final;
public:
  const Strain* getStrain () const
    { return static_cast <const Phyl*> (getParent ()) -> asStrain (); }
    // Return: !nullptr
  void saveFeatures (const string &dir) const;
    // Save Feature::name's
    // Input: dir: !empty()
  size_t getGenes () const
    { return coreNonSingletons + singletons. size (); }
	void getSingletons (Set<Feature::Id> &globalSingletons,
	                    Set<Feature::Id> &nonSingletons) const;
    // Update: globalSingletons, nonSingletons: !intersect()
	void getCommonCore (Set<Feature::Id> &commonCore) const
    { commonCore. intersect (coreSet); }
	size_t removeFromCoreSet (const Set<Feature::Id> &toRemove)
    { size_t n = 0;
      for (const Feature::Id& f : toRemove)
	      n += coreSet. erase (f);
	    return n;
    }
	void setSingletons (const Set<Feature::Id> &globalSingletons)
    { singletons = coreSet;
      singletons. intersect (globalSingletons);
      for (const Feature::Id& f : singletons)
	      coreSet. erase (f);
    }
};



struct ChangeTo;



struct Change : Root  
// Of topology
// An opposite Change exists
{
protected:
	const FeatureTree& tree;
public:
	const Species* from {nullptr};
	  // !nullptr
	Real improvement {INF};
	  // > 0
	VectorPtr<Tree::TreeNode> targets; 
	  // Lowest Phyl's whose parent2core may be changed

  // status fields
private:
  enum Status {eInit, eApplied, eCommitted};
  Status status;
public:
protected:
	Fossil* isolatedTransient;
	  // parent2core[] is not changed
	Species* isolatedTransientChild;
	// (bool)isolatedTransient = (bool)isolatedTransientChild
public:
	
	
protected:
  // To invoke: clear()
	explicit Change (const Species* from_arg)
		: tree (const_cast <FeatureTree&> (from_arg->getFeatureTree ()))
		, from (from_arg)
		, improvement (INF)
		{}
    // Requires: valid_()
	Change (const FeatureTree &tree_arg,
	        istream &is);
public:
	static bool valid_ (const Species* from_arg)
    { return    from_arg
             && from_arg->graph
    	       && from_arg->getParent ();   // Not needed for ChangeToSibling ??
    }
    // Requires: parameters are the same as in the constructor
 ~Change ();
	void qc () const override;
	  // Invokes: valid()
  void saveText (ostream& os) const override;
    // W/o endl
    // Requires: tree topology has not changed after loading the tree
	void print (ostream& os) const override
	  { os << from->getName () << " -> " << type () << "  improvement = ";
	    ONumber oNum (os, 1, false);  // PAR
	    os << improvementDeinflated () << endl; 
	  }
protected:
  void clear () override
    { status = eInit;
      isolatedTransient = nullptr;
		  isolatedTransientChild = nullptr;
    }
    // Output: status fields
public:


  virtual const ChangeTo* asChangeTo () const  
    { return nullptr; }

	  
	virtual const char* type () const = 0;
	virtual bool valid () const = 0;
	virtual Set<const Tree::TreeNode*> getCoreChanged () const
    { Tree::LcaBuffer buf;
    	return Tree::getParents (targets, buf); 
    }
	  // Return: Phyl whose core is supposed to be changed
  // Update: tree
	void apply ();
	  // status: eInit --> eApplied
	  // Output: improvement
	  // Invokes: apply_()
	void restore ();
	  // status: eApplied --> eInit
	  // Invokes: restore_(), clear()
	void commit ();
	  // status: eApplied --> eCommitted
	  // Invokes: commit_()
private:
	virtual void apply_ () = 0;
	virtual void restore_ () = 0;
	virtual void commit_ () = 0;
public:
protected:
  bool isolateTransient (Fossil* f);
//void restoreTransient (); ??
  Real improvementDeinflated () const;
public:
	static bool compare (const Change* a, 
	                     const Change* b)
    // Requires: (bool)a
    { if (a == b)  return false;
      if (! b)  return positive (a->improvement);
      if (a->improvement > b->improvement)  return true;
      if (a->improvement < b->improvement)  return false;
      const int cmp = strcmp (a->type (), b->type ());
      if (cmp < 0)  return true; 
      if (cmp > 0)  return false; 
      if (a->from->index_init < b->from->index_init)  return true;
      if (a->from->index_init > b->from->index_init)  return false;
      return a->better (b);
    }
  virtual bool better (const Change* /*other*/) const 
    { NEVER_CALL; return false; }
};



struct ChangeTo : Change  
{
	const Species* to {nullptr};
	  // !nullptr


protected:	
	ChangeTo (const Species* from_arg,
					  const Species* to_arg)
		: Change (from_arg)
		, to (to_arg)
		{}
	ChangeTo (const FeatureTree &tree_arg,
	          istream &is);
	static bool valid_ (const Species* from_arg,
	                    const Species* to_arg)
	  { return    Change::valid_ (from_arg)
	           && to_arg
	  	       && to_arg->graph == from_arg->graph
	  	       && ! to_arg->descendantOf (from_arg);
	  }
public:
  void saveText (ostream& os) const override;
	void print (ostream& os) const override
	  { os << from->getName () << " ->(" << type () << ") " << to->getName ()
	  	   << "  improvement = " << improvementDeinflated () << endl; 
	  }


  const ChangeTo* asChangeTo () const final
    { return this; }


  bool better (const Change* other) const final;
};



struct Move : ChangeTo
{
  // status fields
protected:
	// !nullptr
	const Fossil* lca {nullptr};
	Fossil* oldParent {nullptr};
	  // Old from->getParent()
	Species* oldParentRepr {nullptr};
	  // Representative of oldParent
public:


protected:	
	Move (const Species* from_arg,
			  const Species* to_arg)
		: ChangeTo (from_arg, to_arg)
		{}
	Move (const FeatureTree &tree_arg,
	      istream &is)
	  : ChangeTo (tree_arg, is)
	  {}
  void clear () override
    { ChangeTo::clear ();
      lca = nullptr;
		  oldParent = nullptr;
		  oldParentRepr = nullptr;
    }
public:
};



struct ChangeToSibling : Move
// *to becomes a sibling of *this
{
  // status fields
private:
	const Fossil* inter {nullptr};
	  // !nullptr <=> after apply_()
	Fossil* arcEnd {nullptr};
public:

	
	ChangeToSibling (const Species* from_arg,
					         const Species* to_arg)
		: Move (from_arg, to_arg)
		{ clear (); 
		  targets = getTargets (); 
		}
	ChangeToSibling (const FeatureTree &tree_arg,
				           istream &is)
    : Move (tree_arg, is)
		{ clear ();
		  targets = getTargets (); 
		}
	static bool valid_ (const Species* from_arg,
	                    const Species* to_arg)
	  { return    Move::valid_ (from_arg, to_arg)
	           && to_arg != to_arg->getTree (). root
	  	       && ! (from_arg->getParent() == to_arg->getParent() && to_arg->getParent() -> arcs [false]. size () == 2)
	  	       && ! (from_arg->getParent() == to_arg              && to_arg              -> arcs [false]. size () == 2);
	  }
private:
	void clear () final
	  { Move::clear ();
		  inter = nullptr;
		  arcEnd = nullptr;
	  }
	VectorPtr<Tree::TreeNode> getTargets () const
		{ return VectorPtr<Tree::TreeNode> {from->getParent (), to->getParent ()}; }
public:


  static const char* type_ ()
	  { return "sibling"; }
	const char* type () const final
	  { return type_ (); }
	bool valid () const final
	  { return valid_ (from, to); }
private:
	void apply_ () final;
	void restore_ () final;
	void commit_ () final;
public:
};



//struct ChangeToChild  ??
// Opposite to ChangeToSibling



struct ChangeToParent : Move
// *to becomes a parent of *this
// Redundant given ChangeToSibling 
{
	ChangeToParent (const Species* from_arg,
					        const Species* to_arg)
		: Move (from_arg, to_arg)
		{ clear ();
		  targets = getTargets (); 
		}
	ChangeToParent (const FeatureTree &tree_arg,
				          istream &is)
    : Move (tree_arg, is)
		{ clear ();
		  targets = getTargets (); 
		}
	static bool valid_ (const Species* from_arg,
	                    const Species* to_arg)
	  { return    Move::valid_ (from_arg, to_arg)
	           && to_arg->asFossil ()
	  	       && to_arg != from_arg->getParent ();
	  }
private:
	VectorPtr<Tree::TreeNode> getTargets () const
		{ return VectorPtr<Tree::TreeNode> {from->getParent (), to}; }
public:


  static const char* type_ ()
	  { return "parent"; }
	const char* type () const final
	  { return type_ (); }
	bool valid () const final
	  { return valid_ (from, to); }
private:
	void apply_ () final;
	void restore_ () final;
	void commit_ () final;
public:
};



//struct ChangeToRoot : ChangeTo  ??
  // from->arc middle -> parent = to 



struct ChangeToUncle : ChangeTo
// *to = uncle
// Swap to->parent and this->parent
// Includes the "swap" Change
{
	ChangeToUncle (const Species* from_arg,
					       const Species* to_arg)
		: ChangeTo (from_arg, to_arg)
		{ clear ();
		  targets = getTargets (); 
		}
	ChangeToUncle (const FeatureTree &tree_arg,
				         istream &is)
    : ChangeTo (tree_arg, is)
		{ clear ();
		  targets = getTargets (); 
		}
	static bool valid_ (const Species* from_arg,
	                    const Species* to_arg)
	  { return    ChangeTo::valid_ (from_arg, to_arg)
	  	       && to_arg->getParent ()
	  	       && from_arg->getParent () -> descendantOf (to_arg->getParent ())
	  	       && ! from_arg->getParent () -> descendantOf (to_arg)
	  	       && from_arg->getParent () != to_arg->getParent ();
	  }
private:
	VectorPtr<Tree::TreeNode> getTargets () const
		{ return VectorPtr<Tree::TreeNode> {from->getParent ()}; }
public:


  static const char* type_ ()
	  { return "uncle"; }
	const char* type () const final
	  { return type_ (); }
	bool valid () const final
	  { return valid_ (from, to); }
private:
	void apply_ () final;
	void restore_ () final;
	void commit_ () final;
public:
};



struct ChangeToCousin : ChangeTo
// *to = first cousin
// Swap to->parent and this->parent
{
	ChangeToCousin (const Species* from_arg,
  					      const Species* to_arg)
		: ChangeTo (from_arg, to_arg)
		{ clear ();
		  targets = getTargets (); 
		}
	ChangeToCousin (const FeatureTree &tree_arg,
				        istream &is)
    : ChangeTo (tree_arg, is)
		{ clear ();
		  targets = getTargets (); 
		}
	static bool valid_ (const Species* from_arg,
	                    const Species* to_arg)
	  { return    ChangeTo::valid_ (from_arg, to_arg)
	  	       && to_arg->getParent ()
	  	       && from_arg->getParent () -> getParent () == to_arg->getParent () -> getParent ()
	  	       && from_arg->getParent () != to_arg->getParent ()
	  	     //&& from_arg->index_init < to_arg->index_init
	  	       ;
	  }
private:
	VectorPtr<Tree::TreeNode> getTargets () const
		{ return VectorPtr<Tree::TreeNode> {from->getParent (), to->getParent ()}; }
public:


  static const char* type_ ()
	  { return "cousin"; }
	const char* type () const final
	  { return type_ (); }
	bool valid () const final
	  { return valid_ (from, to); }
private:
	void apply_ () final;
	void restore_ () final;
	void commit_ () final;
public:
};



struct ChangeRoot : Change
// tree.root = from
// Can be used once in FeatureTree::applyChanges()
// Requires: FeatureTree::allTimeZero
{
  // status fields
private:
	Species* root_old {nullptr};
	  // !nullptr
public:
	
	
	explicit ChangeRoot (const Species* from_arg)
		: Change (from_arg)
		{ clear ();
		  targets = getTargets (); 
		}
private:
	void clear () override
	  { Change::clear ();
	    root_old = nullptr;
	  }
	VectorPtr<Tree::TreeNode> getTargets () const
		{ return VectorPtr<Tree::TreeNode> {from, nullptr}; }
public:
	void qc () const override;


  static const char* type_ ()
	  { return "root"; }
	const char* type () const final
	  { return type_ (); }
	bool valid () const final
	  { return valid_ (from); }
	Set<const Tree::TreeNode*> getCoreChanged () const override
    { return Change::getCoreChanged () << nullptr; }
private:
	void apply_ () final;
	void restore_ () final;
	void commit_ () final;
public:
};



struct ChangeDel : Change
// from is delete'd
// Redundant if !tree.allTimeZero ??
{
  // status fields
private:
	VectorPtr<DiGraph::Node> fromChildren;
	  // = from.getChildren()
	Fossil* oldParent {nullptr};
	  // Of from
public:
	
	
	explicit ChangeDel (const Species* from_arg)
		: Change (from_arg)
		{ clear ();
		  targets << from->getParent (); 
		}
	ChangeDel (const FeatureTree &tree_arg,
		         istream &is)
		: Change (tree_arg, is)
		{ clear ();
		  targets << from->getParent (); 
		}
private:
	void clear () override
	  { Change::clear ();
		  fromChildren. clear ();
		  oldParent = nullptr;
	  }
public:
	static bool valid_ (const Species* from_arg)
	  { return    Change::valid_ (from_arg)
	  	       && from_arg->asFossil ();
	  }


  static const char* type_ ()
	  { return "del"; }
	const char* type () const final
	  { return type_ (); }
	bool valid () const final
	  { return valid_ (from); }
	Set<const Tree::TreeNode*> getCoreChanged () const override
    { return Set<const Tree::TreeNode*> () << targets [0]; }
private:
	void apply_ () final;
	void restore_ () final;
	void commit_ () final;
public:
};



typedef  Set<string>  Nominal;
typedef  map<string/*nominal attribute name*/, Nominal>  Nominals;



struct FeatureTree : Tree
// Of Phyl*
// !allTimeZero => effectively unrooted 
{
  static constexpr Real len_delta {1e-2};

	string inputTreeFName;
	  // !empty()

  // For FeatureTree::len
  static constexpr bool emptySuperRoot {false};  // PAR
  const bool preferGain;
	bool allTimeZero {false}; 
	  // true <=> parsimony method, Species::time is not used
	  // Init: isNan(Species::time) 
	Prob timeOptimFrac {1};
	  // Annealing
	  // allTimeZero => 1
	  // To be increased
  // Valid if !allTimeZero
  Vector<bool> superRootCore;
    // size() == features.size()
  Prob lambda0 {NAN};
		// > 0
  Real time_init {NAN};
    // > 0    
	
  Nominals nominals;
  
	Vector<Feature> features;
	  // Feature::name's are sorted, unique
	Set<Feature::Id> commonCore;
	size_t globalSingletonsSize {0};  
	size_t genomeGenes_ave {0};
	Real len {NAN};
	  // >= len_min
	  // !emptySuperRoot, allTimeZero => parsimony in an undirected tree
	Real len_min {NAN};
	  // >= 0
	Prob lenInflation {0};

  // Internal
	size_t nodeIndex_max {0};
	bool coreSynced {false};
	  // Phyl::core[] corresponds to Phyl::parent2core[]
private:
	Common_sp::AutoPtr <Progress> prog_;
	VectorOwn<Species> toDelete;
public:
  // Stats
	Common_sp::AutoPtr <const Normal> distDistr;
	Common_sp::AutoPtr <const Normal> depthDistr;
	  
	size_t reportFeature {NO_INDEX};

	  
	FeatureTree (const string &treeFName,
  	           const string &geneDir,
  	           const string &coreFeaturesFName,
  	           bool preferGain_arg);
	   // Invokes: loadPhylFile(), Genome::initDir(), finish()
private:
  bool loadPhylLines (const StringVector& lines,
		                  size_t &lineNum,
		                  Species* parent,
		                  size_t expectedOffset);
    // Return: a child of parent has been loaded
    // Output: topology, Species::time
    //         Genome::{human fields} are empty
    // Update: lineNum
  void loadPhylFile (/*int root_species_id,*/
	                   const string &treeFName);
	  // Invokes: loadPhylLines()
  void finish (/*Server* db,*/
               const string &coreFeaturesFName);
    // Input: db or coreFeaturesFName if !allTimeZero
    // Invokes: setLenGlobal(), setCore()
public:
	void qc () const override;
	void saveText (ostream& os) const override;
	  // os <<: if !allTimeZero
	  //          if !timeOptimWhole() then timeOptimFrac
	  //          lambda0
	  //          time_init
	  //        topology
	  // Requires: coreSynced
	void print (ostream& os) const override;

	
  void deleteLeaf (TreeNode* leaf,
                   bool deleteTransientAncestor) override;
    // Requires: features.empty()

	void printInput (ostream& os) const;
	void dump (const string &fName/*,
	           bool setIds*/);
	  // Invokes: setStats()
	void progInternal (const string& step = string ())
	  { if (prog_. get ()) 
	  	  (*prog_) (step);
	  }
	  // Requires: after Progress::Start
	bool featuresExist () const
	  { return ! features. empty (); }
	size_t getTotalGenes () const
	  { return commonCore. size () + features. size () + globalSingletonsSize; }  
  Real getLength_min ()  
    { if (! emptySuperRoot)
        return 0;
      // Needed for E. coli pan ??
      if (allTimeZero)  
    		return (Real) getTotalGenes ();
      // Simplistic model: all Genome::core's are the same, all Species::time = 0 except root
    	const Species* s = static_cast <const Species*> (root);
    	return (Real) getTotalGenes () * s->feature2weight (true,  false)
    	     /*+ commonMissings   * s->feature2weight (false, false)*/;
    }
  // Sankoff algorithm
private:
	Real feature2treeLength (size_t featureIndex) const;
public:
  Real getLength () const;

  // OPTIMIZATION 
  // Phyl::parent2core[]
  void setLenGlobal ()
    { const_static_cast <Species*> (root) -> assignFeaturesDown ();
    	len = getLength ();
    }
  // Phyl::core[]
  void setCore ()
    { coreSynced = true;
      const_static_cast <Species*> (root) -> setCore (); 
    }
  // Phyl::{time,weight[][]}
  Real getRootTime () const
    { return emptySuperRoot ? INF : 0; }
private:
  void setTimeWeight ();
    // Idempotent
    // Optimal if timeOptimWhole()
    // Invokes: setCore(), Species::setTimeWeight()
public:
  void optimizeTime ();
    // Invokes: setTimeWeight(),setLenGlobal()
    // timeOptimWhole() => local optimum is over-stable
	void useTime (const string &coreFeaturesFName);
	  // Output: allTimeZero, superRootCore
	  // Requires: allTimeZero
	  // Invokes: loadSuperRootCoreFile(coreFeaturesFName)
private:
  void getParent2core_sum (size_t parent2core [2/*thisCore*/] [2/*parentCore*/]) const;
    // Input: Species-nodes
    // Output: parent2core[][]
  static Real getLambda0_commonTime (const size_t parent2core [2/*thisCore*/] [2/*parentCore*/],
                                     Real commonTime);
  void optimizeLambda0 ();
    // Update: lambda0, PhyL::weight[][]
//void optimizeLambdaTime ();  
    // Invokes: optimizeTime()
  bool timeOptimWhole () const
    { return eqReal (timeOptimFrac, 1); }
public:
  // Topology
  const Change* getBestChange (const Species* from);
    // Return: May be nullptr
    // Invokes: tryChange()
  bool applyChanges (VectorOwn<Change> &changes);
	  // Return: false <=> finished
    // Update: topology, timeOptimFrac, changes (sort by Change::improvement descending), cout
    // Output: DTNode::stable
    // Print: len, timeOptimFrac
    // Invokes: optimizeTime(), finishChanges()
	bool optimize ();
	  // Return: false <=> finished
	  // Update: Phyl::stable
	  // Invokes: getBestChange(), applyChanges()
	string findRoot (size_t &bestCoreSize);
	  // arg min_node core size
    // Idempotent
    // Return: new root LCA name in the old tree; "" if new root = old root
	  // Output: bestCoreSize, topology, superRootCore
	  // Invokes: setCore(), sort()
	  // Requires: allTimeZero
  void resetSuperRootCore (size_t coreChange [2/*core2nonCore*/]);
    // Idempotent
    // Update: superRootCore, Phyl::core
    // Output: coreChange[]: !core2nonCore = nonCore2core
	  // Requires: !allTimeZero, coreSynced
private:
  void loadSuperRootCoreFile (const string &coreFeaturesFName);
    // Requires: !coreFeaturesFName.empty()
public:
  void saveSuperRootCore (const string &coreFeaturesFName) const;
  bool getSuperRootCore (size_t featureIndex) const
    { const auto& parent2core = static_cast <const Species*> (root) -> parent2core;
      return emptySuperRoot 
           	   ? false
           	   : allTimeZero 
           	     ? lessReal ( parent2core [true]  [featureIndex]. treeLen
           	                , parent2core [false] [featureIndex]. treeLen
           	                )  // superRoot should be minimal
  	    	       : superRootCore [featureIndex];
    }
  void setStats ();
    // Output: Stats
    // Invokes: Phyl::getNeighborDistance()
  void clearStats ();
    // Output: Stats
	void delayDelete (Species* s);
private:
	void tryChange (Change* ch,
	                const Change* &bestChange);
    // Update: bestChange: positive(improvement)
    // Invokes: Change::{apply(),restore()}
public:	
  void finishChanges ();
    // Invokes: deleteTimeZero()
private:
  size_t deleteTimeZero ();
    // Return: # Fossil's delete'd
public:
		
  const Genome* findGenome (const string &genomeId) const;
    // Return: May be nullptr
  size_t findFeature (const Feature::Id &featureName) const;
    // Return: may be NO_INDEX
};



}



#endif


