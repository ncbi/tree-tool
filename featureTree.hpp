// featureTree.hpp

#ifndef FEATURE_TREE_HPP
#define FEATURE_TREE_HPP

#include "common.hpp"
using namespace Common_sp;
#include "numeric.hpp"
#include "dataset.hpp"
using namespace DM_sp;



namespace FeatureTree_sp
{



struct Feature : Named
// name: key
{
  typedef string Id;
  bool isGene;
    // Otherwise a phenotype
  // For FeatureTree::len
  // Valid if !allTimeZero
//Prob lambda_0;  ??
  // Stats
  size_t genomes;
  size_t gains;
  size_t losses;


	Feature (const Id &name_arg,
	         bool isGene_arg);
	Feature ()
	  : isGene (true)
	  , genomes (0)
	  , gains (0)
	  , losses (0)
	  {}
	void qc () const;
	void saveText (ostream& os) const
	  { os << name << " +" << gains << " -" << losses << " / " << genomes << endl; }

	
  bool operator== (const Feature &other) const
    { return    name   == other. name
             && isGene == other. isGene; 
    }
  bool operator< (const Feature &other) const
    { if (isGene > other. isGene)
        return true;
      if (isGene < other. isGene)
        return false;      
      return name < other. name; 
    }
#if 0
  static int id2geneId (const Id &id)
    { return str2<int> (id); }  
  int geneId () const
    { return id2geneId (name); }  
    // For the database
#endif
  size_t mutations () const
    { return gains + losses; }
  bool better (const Feature* other) const
    { return    ! other 
             || mutations () < other->mutations ()
             || (mutations () == other->mutations () && name.size () > other->name. size ());
    }
};



#if 0
struct TargetFeature : Named  
{
  Feature::Id featureId;
    // !empty()
  size_t index;
    // May be NO_INDEX      
  size_t serial;

  TargetFeature ()
    : index (NO_INDEX)
    , serial (0)
    {}
  TargetFeature (const string &name_arg,
                 const Feature::Id &featureId_arg,
                 size_t index_arg)
    : Named (name_arg)
    , featureId (featureId_arg)
    , index (index_arg)
    , serial (0)
    {}
  bool empty () const
    { return Named::empty () && featureId. empty (); }
  void clear ()
    { Named::clear ();
      featureId. clear ();
      index = NO_INDEX;
      serial = 0;
    }
  void read (istream &is)
	  { is >> featureId;
	    Named::read (is); 
	  }
	  
	string getName () const
	  { return name + (serial ? toString (serial) : ""); }
};
#endif



inline bool eqTreeLen (Real len1,
                       Real len2)
  { return eqReal (len1, len2, 1e-4 /*PAR*/); }  



struct FeatureTree;

//struct Phyl;
  struct Species;
    struct Fossil;
    struct Strain;
  struct Genome;



struct Phyl : Tree::TreeNode 
{
  // For FeatureTree::len
	struct CoreEval
	{
		Real treeLen;
		  // >= 0
		ebool core;
		
		CoreEval ()
		  : treeLen (0)
		  , core (EFALSE)
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
	  // !getFeatureTree().inputTreeFName.empty() => matches the orginial node number in DFS
#if 0
	bool hasPhenChange;
	  // Phenotype is changed in the subtree inclusive
#endif
private:
  bool stable;
    // Init: false
  friend struct FeatureTree;
public:


protected:
	Phyl (FeatureTree &tree,
        Species* parent_arg);
	  // To be followed by setWeight()
public:
	void init ();
  void qc () const;
protected:
  void saveContent (ostream& os) const;
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
	bool feature2core (size_t featureIndex) const
	  { const bool parentCore = feature2parentCore (featureIndex);
			const ebool c = parent2core [parentCore] [featureIndex]. core;
		  return c == UBOOL ? parentCore /*PAR*/ : (bool) c;
	  }
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
#ifndef NDEBUG
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

#if 0
  string getPhenChange (bool skipSingleNonSystems) const;
    // Return: ';'-delimited
    // Input: Features:<Stats>
protected:
  virtual string getTargetFeatureChange (const TargetFeature &tf) const;
public:
  string getTargetFeaturesChange () const;
  void setHasPhenChange ();
    // Output: hasPhenChange
    // Invokes: getPhenChange(), getTargetFeaturesChange()
#endif

#if 0
  // Requires: getFeatureTree().coreSynced
  virtual void saveDatabaseTopology (Server &db,
                                     bool &parent_STRAIN) = 0;
    // Update: parent_STRAIN: valid if getParent()
  virtual void saveDatabaseNode (Server &db) const = 0;
    // Does not insert into *_PHEN tables
  virtual void saveDatabasePhen (Server &db) const = 0;
    // Insert into *_PHEN tables
#endif
};



struct Species : Phyl
// core: Function of getFeatureTree().rootCore[] and Phyl::parent2core[]
{
  string id;
    // !empty() => SPECIES.id
  // Valid if !allTimeZero
	Real time;
	  // isNan() || >= 0
    // Exponential distribution, if yes then use P_{prior} ??
//Real timeVariance;  // Read the MLE theory ??
	Real pooledSubtreeDistance;  
	  // Part of the arc (this,getParent()) length due to FeatureTree::commonCore and Genome::singletons
	  // >= 0
private:
	Real weight_old [2/*thisCore*/] [2/*parentCore*/];
  Real time_old;
    // May be NAN
	Real pooledSubtreeDistance_old;  
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
  bool movementsOn;
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
  void qc () const;
  void saveContent (ostream& os) const;


  const Species* asSpecies () const
    { return this; }


  string getName () const
    { return id. empty () ? Phyl::getName () : id; }
	double getParentDistance () const;
  string getNewickName (bool minimal) const;
  void setWeight ();
    // Input: time, getFeatureTree().lambda0
	void setCore ();
private:
	void setCoreEval (size_t featureIndex,
	                  bool parentCore,
	                  CoreEval ce)
	  {	CoreEval& ce_old = parent2core [parentCore] [featureIndex];
		  if (ce_old == ce)
		  	return;
		  if (movementsOn)
			  movements << Movement (parentCore, featureIndex, ce_old);
	    ce_old = ce;
	  }
public:
protected:
  void assignFeatures ()
    { Phyl::assignFeatures ();
      pooledSubtreeDistance = getPooledSubtreeDistance ();
    }
public:
#if 0
  void saveDatabaseTopology (Server &db,
                             bool &parent_STRAIN);
    // Output: id: incremented in DFS preorder
  void saveDatabaseNode (Server &db) const;
  void saveDatabasePhen (Server &db) const;
#endif

  // Sankoff algorithm
private:
	Real feature2treeLength (size_t featureIndex) const;
	  // Requires: this == getFeatureTree().root
public:
	Real root2treeLength () const;
	  // Tree length if *this is the root

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
	void qc () const;


  const Fossil* asFossil () const
    { return this; }

  bool isInteriorType () const
    { return true; }

  Real getPooledSubtreeDistance () const;

  void setId (uint &id_arg);
    // Output: id
    // Update: id_arg
};



struct Strain : Species
// Distance from children to *this depends on Genome quality
{
  bool singletonsInCore;
    // getGenome()->singletons: true <=> in this->core, false <=> annotation errors in Genome
    // Init: true
    // Analog of parent2core[false] and core
private:
  bool singletonsInCore_old;
public:


	Strain (FeatureTree &tree,
	        Fossil* parent_arg,
	        const string &id_arg,
	        Real time_arg)
		: Species (tree, parent_arg, id_arg, time_arg)
    , singletonsInCore (true)
    , singletonsInCore_old (true)
		{}
  void qc () const;
  void saveContent (ostream& os) const;


  const Strain* asStrain () const
    { return this; }


  string getName () const 
    { return "s" + id; }
  string getNewickName (bool /*minimal*/) const
    { return string (); }

  void getParent2corePooled (size_t parent2corePooled [2/*thisCore*/] [2/*parentCore*/]) const;
	void assignFeatures ();
	Real getPooledSubtreeDistance () const;
private:
	void rememberFeatures ()
	  { Species::rememberFeatures ();
	    singletonsInCore_old = singletonsInCore;
	  }
	void restoreFeatures ()
	  { singletonsInCore = singletonsInCore_old;
	    Species::restoreFeatures ();
	  }
#if 0
  string getTargetFeatureChange (const TargetFeature &tf) const;
#endif
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
  
	Vector<bool> optionalCore;
    // size() = getFeatureTree().features.size()
  bool hasPhens;
private:
  Set<Feature::Id> coreSet;
  map<Feature::Id,bool/*optional*/> phens; 
  friend struct FeatureTree;
public:
  size_t coreNonSingletons;
    // Includes optionalCore[]
  Set<Feature::Id> singletons;
//Set<Feature::Id> missedCore;  ??
    // Opposite to singletons

#if 0
  // Quality    
  // For annotProb
  uint L50;    
  // "System" phenotypes
  Feature::Id genus;
  Feature::Id taxSpecies;
  Feature::Id subspecies;
  Feature::Id serovar;
  Feature::Id pathovar;
  Feature::Id continent;

  // For human
  size_t oddCdss;  
  uint project_id;
  uint tax_id;
  string taxName;  
  string sequencer;
  uint pubmed;
  string outbreakName;
  uint outbreakYear;
  string phylogeneticClass;
  string isolate;
#endif
  
    
	Genome (FeatureTree &tree,
	        Strain* parent_arg,
	        const string &id_arg);
	  // To be followed by: {initDb()|initDir()}, init()
#if 0
	void initDb (Server &db);
	  // Output: coreSet, phens, coreNonSingletons
#endif
	void initDir (const string &geneDir,
                const string &phenDir);
	  // Input: file "geneDir/id" with the format: 
	  //          <gene>*
	  //        file "phenDir/id" with the format: 
	  //   [nophenotypes]
	  //   <phenotype> {0|1}*  // 1 <=> optional
	  // Output: coreSet, phens, coreNonSingletons
	  // Invokes: addPhen()
#if 0
private:
  void addSystemPhens ()
    { addPhen (subspecies, false);
      addPhen (taxSpecies, false);
      addPhen (genus,      false);
      addPhen (serovar,    false);
      addPhen (pathovar,   false);
      addPhen (continent,  false);
    }
public:
#endif
	void init (const map <Feature::Id, size_t/*index*/> &feature2index);
	  // Output: CoreEval::core, core
	void qc () const;
	void saveContent (ostream& os) const;


  const Genome* asGenome () const
    { return this; }


  string getName () const 
    { return "g" + id; }
	double getParentDistance () const
	  { return 0; }
  string getNewickName (bool /*minimal*/) const
	  { return id /*+ " " + taxName + getNameExtra ()*/; 
	  }
  bool isLeafType () const
    { return true; }

  void setWeight ();
  void getParent2corePooled (size_t parent2corePooled [2/*thisCore*/] [2/*parentCore*/]) const;
	void setCore ();
private:
	void setCoreEval (size_t featureIndex,
	                  bool parentCore,
	                  CoreEval ce)
	  {	parent2core [parentCore] [featureIndex] = ce; }
	void assignFeature (size_t featureIndex);
#if 0
  string getTargetFeatureChange (const TargetFeature &tf) const;
#endif
public:
#if 0
  void saveDatabaseTopology (Server &db,
                             bool &parent_STRAIN);
    // Invokes: SPECIES_redo()
  void saveDatabaseNode (Server &db) const;    
  void saveDatabasePhen (Server &db) const;
#endif

#if 0
  string getNameExtra () const
    { string s;
     	if (! outbreakName. empty ())
		    s += "  " + outbreakName + " (" + toString (outbreakYear) + ")";
      if (! isolate. empty ())
      	s += "  " + isolate;
      if (! pathovar. empty ())
      	s += "  " + pathovar;
      return s;
    }
#endif
  const Strain* getStrain () const
    { return static_cast <const Phyl*> (getParent ()) -> asStrain (); }
    // Return: !nullptr
#if 0
  bool isSystemPhen (const string &phenName) const
    { return    phenName == subspecies
             || phenName == taxSpecies
             || phenName == genus
             || phenName == serovar
             || phenName == pathovar
             || phenName == continent;
    }
#endif
  void saveFeatures (const string &dir) const;
    // Save Feature::name's
    // Input: dir: !empty()
private:
  void addPhen (const Feature::Id &phen,
                bool optional);
    // Update: phens
public:
  size_t getGenes () const
    { return coreNonSingletons + singletons. size (); }
  // Use optionalCore[] ??
	void getSingletons (Set<Feature::Id> &globalSingletons,
	                    Set<Feature::Id> &nonSingletons) const;
    // Update: globalSingletons, nonSingletons: !intersect()
	void getCommonCore (Set<Feature::Id> &commonCore) const
    { commonCore. intersect (coreSet); }
	size_t removeFromCoreSet (const Set<Feature::Id> &toRemove)
    { return coreSet. setMinus (toRemove); }
	void setSingletons (const Set<Feature::Id> &globalSingletons)
    { singletons = coreSet;
      singletons. intersect (globalSingletons);
      coreSet. setMinus (singletons); 
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
	const Species* from;
	  // !nullptr
	Real improvement;
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
	static bool valid_ (const Species* from_arg)
    { return    from_arg
             && from_arg->graph
    	       && from_arg->getParent ();   // Not needed for ChangeToSibling ??
    }
    // Requires: parameters are the same as in the constructor
public:
 ~Change ();
	void qc () const;
	  // Invokes: valid()
  void saveText (ostream& os) const;
    // W/o endl
    // Requires: tree topology has not changed after loading the tree
	void print (ostream& os) const
	  { os << from->getName () << " -> " << type () << "  improvement = ";
	    ONumber oNum (os, 1, false);  // PAR
	    os << improvementDeinflated () << endl; 
	  }
protected:
  void clear ()
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
    { return Tree::getParents (targets); }
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
	const Species* to;
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
	  	       && ! to_arg->descendentOf (from_arg);
	  }
public:
  void saveText (ostream& os) const;
	void print (ostream& os) const
	  { os << from->getName () << " ->(" << type () << ") " << to->getName ()
	  	   << "  improvement = " << improvementDeinflated () << endl; 
	  }


  const ChangeTo* asChangeTo () const  
    { return this; }


  bool better (const Change* other) const;
};



struct Move : ChangeTo
{
  // status fields
protected:
	// !nullptr
	const Fossil* lca;
	Fossil* oldParent;
	  // Old from->getParent()
	Species* oldParentRepr;
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
  void clear ()
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
	const Fossil* inter;
	  // !nullptr <=> after apply_()
	Fossil* arcEnd;
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
	void clear ()
	  { Move::clear ();
		  inter = nullptr;
		  arcEnd = nullptr;
	  }
	VectorPtr<Tree::TreeNode> getTargets () const
		{ return VectorPtr<Tree::TreeNode>::make (from->getParent (), to->getParent ()); }
public:


  static const char* type_ ()
	  { return "sibling"; }
	const char* type () const 
	  { return type_ (); }
	bool valid () const
	  { return valid_ (from, to); }
private:
	void apply_ ();
	void restore_ ();
	void commit_ ();
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
		{ return VectorPtr<Tree::TreeNode>::make (from->getParent (), to); }
public:


  static const char* type_ ()
	  { return "parent"; }
	const char* type () const 
	  { return type_ (); }
	bool valid () const
	  { return valid_ (from, to); }
private:
	void apply_ ();
	void restore_ ();
	void commit_ ();
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
	  	       && from_arg->getParent () -> descendentOf (to_arg->getParent ())
	  	       && ! from_arg->getParent () -> descendentOf (to_arg)
	  	       && from_arg->getParent () != to_arg->getParent ();
	  }
private:
	VectorPtr<Tree::TreeNode> getTargets () const
		{ return VectorPtr<Tree::TreeNode>::make (from->getParent ()); }
public:


  static const char* type_ ()
	  { return "uncle"; }
	const char* type () const 
	  { return type_ (); }
	bool valid () const
	  { return valid_ (from, to); }
private:
	void apply_ ();
	void restore_ ();
	void commit_ ();
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
		{ return VectorPtr<Tree::TreeNode>::make (from->getParent (), to->getParent ()); }
public:


  static const char* type_ ()
	  { return "cousin"; }
	const char* type () const 
	  { return type_ (); }
	bool valid () const
	  { return valid_ (from, to); }
private:
	void apply_ ();
	void restore_ ();
	void commit_ ();
public:
};



struct ChangeRoot : Change
// tree.root = from
// Can be used once in FeatureTree::applyChanges()
// Requires: FeatureTree::allTimeZero
{
  // status fields
private:
	Species* root_old;
	  // !nullptr
public:
	
	
	explicit ChangeRoot (const Species* from_arg)
		: Change (from_arg)
		{ clear ();
		  targets = getTargets (); 
		}
#if 0
	ChangeRoot (const FeatureTree &tree_arg,
		          istream &is)
		: Change (tree_arg, is)
	  , root_old (nullptr)
		{ clear ();
		  targets = getTargets (); 
		}
#endif
private:
	void clear ()
	  { Change::clear ();
	    root_old = nullptr;
	  }
	VectorPtr<Tree::TreeNode> getTargets () const
		{ return VectorPtr<Tree::TreeNode>::make (from, nullptr); }
public:
	void qc () const;


  static const char* type_ ()
	  { return "root"; }
	const char* type () const 
	  { return type_ (); }
	bool valid () const
	  { return valid_ (from); }
	Set<const Tree::TreeNode*> getCoreChanged () const
    { return Change::getCoreChanged () << nullptr; }
private:
	void apply_ ();
	void restore_ ();
	void commit_ ();
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
	Fossil* oldParent;
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
	void clear ()
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
	const char* type () const 
	  { return type_ (); }
	bool valid () const
	  { return valid_ (from); }
	Set<const Tree::TreeNode*> getCoreChanged () const
    { return Set<const Tree::TreeNode*> () << targets [0]; }
private:
	void apply_ ();
	void restore_ ();
	void commit_ ();
public:
};



struct FeatureTree : Tree
// Of Phyl*
{
  static const Real len_delta;

	string inputTreeFName;
	  // May be empty()
	bool genesExist;

  // For FeatureTree::len
  static const bool emptyRoot = false;
//Topology
	bool allTimeZero; 
	  // true <=> parsimony method, Species::time is not used
	  // Init: isNan(Species::time) 
	Prob timeOptimFrac;
	  // Annealing
	  // Init: 1
	  // allTimeZero => 1
	  // To be increased
  // Valid if !allTimeZero
  Vector<bool> rootCore;
    // size() == features.size()
  Prob lambda0;
		// > 0
  Real time_init;   
    // > 0    
	
	Vector<Feature> features;
	  // Genes and phenotypes in this order
	  // Feature::name's are unique
	size_t genes;
	  // <= features.size()
	Set<Feature::Id> commonCore;
	size_t globalSingletonsSize;  
	size_t genomeGenes_ave;
  List<string> taxNamePrefix;  // ??
	Real len;
	  // >= len_min
	Real len_min;
	  // >= 0
	Prob lenInflation;
	  // Init: 0

  // Internal
	size_t nodeIndex_max;
	bool coreSynced;
	  // Phyl::core[] corresponds to Phyl::parent2core[]
private:
	Common_sp::AutoPtr <Progress> prog_;
	VectorOwn<Species> toDelete;
public:
  // Stats
	Common_sp::AutoPtr <const Normal> distDistr;
	Common_sp::AutoPtr <const Normal> depthDistr;


#if 0
  bool savePhenChangesOnly;
    // Init: false
  Vector<TargetFeature> targetFeatures;
    // Same name => different serial
#endif
	  
	  
#if 0
	FeatureTree (int root_species_id,
            const string &treeFName,
            Server &db);
	   // Input: treeFName: may be empty()
	   // Invokes: if treeFName.empty() then loadPhylDb() else loadPhylLines(), Genome::initDb(), abbreviate(), finish()
#endif
	FeatureTree (const string &treeFName,
	          const string &geneDir,
	          const string &phenDir,
	          const string &coreFeaturesFName,
	          bool genesExist_arg);
	   // Invokes: loadPhylFile(), Genome::initDir(), finish()
private:
#if 0
  void loadPhylDb (Server &db,
                   int species_id,
                   Fossil* parent);
    // Output: topology, Species::time
#endif
  bool loadPhylLines (const Vector<string>& lines,
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
	void qc () const;
	void saveText (ostream& os) const;
	  // os <<: if !allTimeZero
	  //          if !timeOptimWhole() then timeOptimFrac
	  //          lambda0
	  //          time_init
	  //        topology
	  // Requires: coreSynced
	void print (ostream& os) const;

	
  void deleteLeaf (TreeNode* leaf,
                   bool deleteTransientAncestor);
    // Requires: features.empty()

	void printInput (ostream& os) const;
	void dump (const string &fName,
	           bool setIds);
	  // Invokes: setStats()
	void progInternal (const string& step = string ())
	  { if (prog_. get ()) 
	  	  (*prog_) (step);
	  }
	  // Requires: after Progress::Start
	size_t getTotalGenes () const
	  { return commonCore. size () + genes + globalSingletonsSize; }  
  Real getLength () const
    { return static_cast <const Species*> (root) -> root2treeLength (); }
  Real getLength_min ()  
    { if (! emptyRoot)
        return 0;
      // Needed for E. coli pan ??
      if (allTimeZero)  
    		return (Real) getTotalGenes ();
      // Simplistic model: all Genome::core's are the same, all Species::time = 0 except root
    	const Species* s = static_cast <const Species*> (root);
    	return   (Real) getTotalGenes () * s->feature2weight (true,  false)
    	     /*+ commonMissings   * s->feature2weight (false, false)*/;
    }

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
    { return emptyRoot ? INF : 0; }
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
	  // Output: allTimeZero, rootCore
	  // Requires: allTimeZero
	  // Invokes: loadRootCoreFile(coreFeaturesFName)
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
    // Invokes: optimizeTime(), finishChanges()
	bool optimize ();
	  // Return: false <=> finished
	  // Update: Phyl::stable
	  // Invokes: getBestChange(), applyChanges()
	void findRoot ();
    // Idempotent
	  // Output: topology, rootCore
	  // Invokes: setCore()
	  // Requires: allTimeZero
  void resetRootCore (size_t coreChange [2/*core2nonCore*/]);
    // Idempotent
    // Update: rootCore, Phyl::core
    // Output: coreChange[]: !core2nonCore = nonCore2core
	  // Requires: !allTimeZero, coreSynced
private:
#if 0
  void loadRootCoreDb (Server* db);
    // Requires: (bool)db
#endif
  void loadRootCoreFile (const string &coreFeaturesFName);
    // Requires: !coreFeaturesFName.empty()
public:
  void saveRootCore (const string &coreFeaturesFName) const;
  bool getRootCore (size_t featureIndex) const
    { const Species* root_ = static_cast <const Species*> (root);
      return emptyRoot 
           	   ? false
           	   : allTimeZero || ! genesExist
           	     ? root_->parent2core [false] [featureIndex]. treeLen >
  	    	         root_->parent2core [true]  [featureIndex]. treeLen
  	    	           // Lesser |core| is preferred: tie => false 
  	    	       : rootCore [featureIndex];
    }
  void setStats ();
    // Output: Stats
    // Invokes: Phyl::getNeighborDistance()
  void clearStats ();
    // Output: Stats
//size_t badNodesToRoot ();
    // Return: # bad nodes moved to root
    // Invokes: setDistr()
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
		
//size_t addStrains ();
	  // Return: # strains added

#if 0	
private:	
	void abbreviate ();
	  // Update: Genome::taxName's
	  // Output: taxNamePrefix
	  // Idempotent
public:
	string abbreviationLegend () const;
#endif
    
  const Genome* findGenome (const string& genomeId) const;
    // Return: May be nullptr
    
#if 0
  void loadTargetFeatures (const string &fName);
    // fName file line format: <Feature::name> <Print name>
#endif

#if 0
  // Update: prog_
  void saveDatabase (Server &db,
                     int root_species_id);
    // Update: Fossil::id
    // Invokes: saveDatabasePhen()
  void saveDatabasePhen (Server &db);
#endif
};



}



#endif


