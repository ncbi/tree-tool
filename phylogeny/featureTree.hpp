// featureTree.hpp

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
*   Feature tree
*
*/


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



struct FeatureTree;

struct Phyl;
  struct Species;
    struct Fossil;
    struct Strain;
  struct Genome;



struct Feature : Named
// Boolean attribute of Genome
// Id name: key
{
  typedef string Id;

  // Stats
  size_t genomes {0};
    // Non-optional
  size_t optionalGenomes {0}; 
  bool rootGain {false}; 
  VectorPtr<Phyl> gains;  
  VectorPtr<Phyl> losses; 
  array<Real,2/*bool*/> len;
    // Root-independent


	explicit Feature (const Id &name_arg)
    : Named (name_arg)
    {}
  Feature () = default;
	void qc () const override;
	void saveText (ostream& os) const override
	  { os         << name 
	       << "\t" << gains. size () 
	       << "\t" << losses. size () 
	       << "\t" << gains. size () + losses. size () - 1  // criterion
	       << "\t" << genomes 
	       << "\t" << optionalGenomes;
	    if (! isNan (getLambda (false)))
	    { os << "\t" << getLambda (false)
	         << "\t" << getLambda (true);
	    }
	    os << endl; 
	  }

	
	string getNominVar () const
	  { const size_t pos = name. find (':');
	    if (pos == string::npos)
	      return noString;
	    return name. substr (0, pos);
	  }
  bool operator== (const Feature &other) const
    { return name == other. name; }
  bool operator< (const Feature &other) const
    { return name < other. name; }    
  static bool statEqual (const Feature& a,
                         const Feature& b)
    { return    a. rootGain == b. rootGain
             && a. gains    == b. gains
    	       && a. losses   == b. losses;
    }
  size_t allGenomes () const
    { return genomes + optionalGenomes; }
  size_t realGains () const
    { return rootGain + gains. size (); }
  static bool statLess (const Feature& a,
                        const Feature& b);
  size_t mutations () const
    { return gains. size () + losses. size (); }
  Real getLambda (bool core) const
    { return (Real) mutations () / len [core]; }
    // Approximate estimation
  bool monophyletic () const
    { return realGains () == 1 && losses. empty (); }
  bool better (const Feature* other) const
    { return    ! other 
             || mutations () < other->mutations ()
             || (mutations () == other->mutations () && name.size () > other->name. size ());
    }
  void clearStats ()
    { genomes = 0;
      optionalGenomes = 0;
      rootGain = false;
      gains. clear ();
      losses. clear ();
      for (const bool b : {false, true})
        len [b] = 0.0;
    }
  size_t getLeaves (const Phyl* phyl,
                    bool gained) const;
    // Return: > 0
  
  static bool nominalSingleton (const Id &featureId);
};



typedef  unordered_map<Feature::Id, size_t>  Feature2index;



inline bool eqTreeLen (float len1,
                       float len2)
  { return eqReal (len1, len2, 1e-4 /*PAR*/); }  



struct Phyl : Tree::TreeNode 
{
  friend FeatureTree;

  // For FeatureTree::len
	struct CoreEval
	{
		float treeLen {0.0};
		  // >= 0
		ebool core {efalse};
		
		CoreEval (float treeLen_arg,
		          ebool core_arg)
		  : treeLen (treeLen_arg)
		  , core (core_arg)
		  {}
    CoreEval () = default;
		  
		bool operator== (const CoreEval &other) const
		  { return    eqReal (treeLen, other. treeLen)
		  	       && core == other. core;
		  }
	};
	Vector<CoreEval> parent2core [2/*bool parentCore*/];
	  // CoreEval::core: optimal given parentCore
    // size() = getFeatureTree().features.size()
	float weight [2/*thisCore*/] [2/*parentCore*/];
	  // = -log(prob); >= 0; may be inf
	Vector<bool> core;
    // size() = getFeatureTree().features.size()
	size_t index_init;
	  // Matches the orginial node number in DFS
private:
  bool stable {false};
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
	float feature2weight (size_t featureIndex,
	                      bool thisCore,
	                      bool parentCore) const;
	  // Return: >= 0	  
	float feature2weight (bool thisCore,
	                     bool parentCore) const
	  { return feature2weight (no_index/*??*/, thisCore, parentCore); }
	  
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
  float getDistance () const;
    // Distance to getParent() for features in getFeatureTree().features
    // Return: >= 0
    // Input: core
    // Invokes: getPooledDistance()
	float getSubtreeLength () const;  
	  // Return: sum_{n \in subtree} n->getDistance() 
#if 0
	float getFeatureDistance (size_t featureIndex) const;
	float getFeatureSubtreeLength (size_t featureIndex) const;  
	float getFeatureTreeLen (size_t featureIndex) const
	  { return parent2core [feature2parentCore (featureIndex)] [featureIndex]. treeLen; }
#endif
  float getNeighborDistance () const;
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
	float getPooledDistance () const;
    // Distance to getParent() for features not in getFeatureTree().features
    // Invokes: getParent2corePooled()
};



struct Species : Phyl
// core: Function of getFeatureTree().superRootCore[] and Phyl::parent2core[]
{
  string id;
    // !empty() => SPECIES.id
  // Valid if !allTimeZero
	Real time {NaN};
	  // isNan() || >= 0
    // Exponential distribution, if yes then use P_{prior} ??
//Real timeVariance;  // Read the MLE theory ??
	float pooledSubtreeDistance {numeric_limits<float>::infinity ()};  
	  // Part of the arc (this,getParent()) length due to FeatureTree::commonCore and Genome::singletons
	  // >= 0
private:
  friend FeatureTree;
	float weight_old [2/*thisCore*/] [2/*parentCore*/];
  Real time_old {NaN};
    // May be NaN
	float pooledSubtreeDistance_old {(float) NaN};  
  struct Movement
  {
  	bool parentCore {false};
  	size_t featureIndex {no_index};
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
  size_t middleCoreSize {0};
    // Size of core in the middle of the arc
    // Undefined for FeatureTree::root
public:


protected:
	Species (FeatureTree &tree,
	         Fossil* parent_arg,
	         const string &id_arg,
	         Real time_arg);
	  // To be followed by: init()
public:
  string getName () const override
    { return id. empty () ? Phyl::getName () : id; }
  void qc () const override;
  void saveContent (ostream& os) const override;


  const Species* asSpecies () const final
    { return this; }

  string getNewickName (bool /*minimal*/) const final
    { return noString; }

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
	  { time_old = NaN; }

	virtual float getPooledSubtreeDistance () const = 0;
	  
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
  
  bool getMiddleCore (size_t featureIndex) const
    { return getParent () ? feature2parentCore (featureIndex) && core [featureIndex] : false; }
  void setMiddleCoreSize ();
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

  float getPooledSubtreeDistance () const final;

  void setId (uint &id_arg);
    // Output: id
    // Update: id_arg
};



struct Strain : Species
// Distance from children to *this depends on Genome quality
// Not neded for maximum parsimony ??!
{
  bool singletonsInCore {true};
    // getGenome()->singletons: true <=> in this->core, false <=> annotation errors in Genome
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
  string getName () const final
    { return "s" + id; }
  void qc () const override;
  void saveContent (ostream& os) const override;


  const Strain* asStrain () const final
    { return this; }

	void assignFeatures () final;
  void getParent2corePooled (size_t parent2corePooled [2/*thisCore*/] [2/*parentCore*/]) const final;
	float getPooledSubtreeDistance () const final;
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
  friend FeatureTree;
  
  struct GenomeFeature  // = ObjFeature in evolution.hpp ??
  { 
    Feature::Id id;
    bool optional {false};
    
    GenomeFeature (const Feature::Id &id_arg,
                   bool optional_arg)
      : id (id_arg)
      , optional (optional_arg)
      {}
    GenomeFeature () = default;
    
    static bool isFeatureNominal (const string &s)
      { return contains (s, ':'); }
    bool isFeatureNominal () const
      { return isFeatureNominal (id); }
    bool operator< (const GenomeFeature& gf) const
      { return id < gf. id; }
    bool operator== (const GenomeFeature& gf) const
      { return id == gf. id; }
  };

  string id;  
    // !empty()
  
	Vector<bool> optionalCore;
    // size() = getFeatureTree().features.size()
    // Nominal attribute => true
private:
  // Temporary
  typedef  Vector<GenomeFeature>  CoreSet;
    // searchSorted
  CoreSet coreSet; 
  StringVector nominals;
    // Nominal attribute names
    // Subset of getFeatureTree().nominal2values
public:
  size_t coreNonSingletons {0};
    // Does not include optionalCore[]  
  Vector<Feature::Id> singletons;
    // searchSorted
//Set<Feature::Id> missedCore;  ??
    // Opposite to singletons


	Genome (FeatureTree &tree,
	        Strain* parent_arg,
	        const string &id_arg);
	  // To be followed by: initDir(), init()
  static string featureLineFormat ();
	void initDir (const string &featureDir,
	              bool large,
	              bool nominalSingletonIsOptional);
	  // Input: file "featureDir/id" with the format: `featureLineFormat()`
    //        large: files in featureDir are grouped into subdirectories named str2hash_class(<file name>)
	  // Output: coreSet, coreNonSingletons
private:
	void coreSet2nominals ();
	  // Update: getFeatureTree().nominal2values, nominals
	void nominals2coreSet ();
	  // Update: coreSet: add optional GenomeFeature's
	void getSingletons (Set<Feature::Id> &globalSingletons,
	                    Set<Feature::Id> &nonSingletons) const;
    // Update: globalSingletons, nonSingletons: !intersect()
	void getCommonCore (Set<Feature::Id> &commonCore) const
    { for (Iter<Set<Feature::Id>> iter (commonCore); iter. next (); )
				if (! coreSet. containsFast (GenomeFeature (*iter, false/*irrelevant*/)))
					iter. erase ();
    }
	void setSingletons (const Set<Feature::Id> &globalSingletons);
	  // Update: coreSet, singletons, coreNonSingletons
	void init (const Feature2index &feature2index);
	  // Output: CoreEval::core, core
	void init (const Vector<size_t> &featureBatch);
public:
  string getName () const final
    { return "g" + id; }
	void qc () const override;
	void saveContent (ostream& os) const override;


  const Genome* asGenome () const final
    { return this; }

  string getLeafName () const final
    { return id; }
	double getParentDistance () const final
	  { return 0.0; }
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
  size_t getFeatures () const
    { return coreNonSingletons + singletons. size (); }
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
	float improvement {numeric_limits<float>::infinity ()};
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
		: tree (var_cast (from_arg->getFeatureTree ()))
		, from (from_arg)
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
 ~Change ()
    {
    #ifndef NDEBUG
      if (status == eApplied)
        errorExit ("Change::status = eApplied");
    #endif
    }
	void qc () const override;
	  // Invokes: valid()
  void saveText (ostream& os) const override;
    // W/o endl
    // Requires: tree topology has not changed after loading the tree
	virtual void print (ostream& os) const
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
      if (! b)  return a->improvement > 0.0;
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
    { throwf ("NEVER_CALL"); }
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
	static bool valid_ (const Species* from_arg)
    { return    Change::valid_ (from_arg)
             && from_arg->getParent ();
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



struct FeatureTree : Tree
// Of Phyl*
// !allTimeZero => effectively unrooted 
{
  static constexpr Real len_delta {1e-2};

	string inputTreeFName;
	  // !empty()

  // For FeatureTree::len
  static constexpr bool emptySuperRoot {false};  // PAR
  const bool nominalSingletonIsOptional;
  const bool preferGain;  
	bool allTimeZero {false}; 
	  // true <=> parsimony method, Species::time is not used
	  // Init: = isNan(Species::time) 
	Prob timeOptimFrac {1.0};
	  // Annealing
	  // allTimeZero => 1
	  // To be increased
  // Valid if !allTimeZero
  Vector<bool> superRootCore;
    // size() == features.size()
  Prob lambda0 {NaN};
		// > 0
  Real time_init {NaN};
    // > 0    
	
private:
  friend Genome;
  map <string/*nominal attribute name*/, Set<string>/*nominal attribute values*/>  nominal2values;
public:  
	Vector<Feature> features;
	  // Feature::name's are sorted, unique
	  // !Feature::nominalSingleton(Feature::Id)
	Set<Feature::Id> commonCore;
	size_t globalSingletonsSize {0};  
	Real genomeFeatures_ave {0.0};
	float len {(float) NaN};
	  // >= len_min
	  // !emptySuperRoot, allTimeZero => parsimony in an undirected tree
	float len_min {(float) NaN};
	  // >= 0
	Prob lenInflation {0.0};
	bool oneFeatureInTree {false};
    // RAM is limited, no topology optimization, only maximum parsimony
  static constexpr size_t featureBatchSize {10000};  // PAR
    // 164777 Genome's, 114555 features: 26% of RAM of lmem21
    // "Assigning features": 15 min./10000 features

  // Internal
	size_t nodeIndex_max {0};
	bool coreSynced {false};
	  // Phyl::core[] corresponds to Phyl::parent2core[]
private:
	unique_ptr<Progress> prog_;
	VectorOwn<Species> toDelete;
public:
  // Stats
	unique_ptr<const Normal> distDistr;
	unique_ptr<const Normal> depthDistr;
	  
	  
	FeatureTree (const string &treeFName,
  	           const string &featureDir,
  	           bool large,
  	           const string &coreFeaturesFName,
  	           bool nominalSingletonIsOptional,
  	           bool preferGain_arg,
  	           bool oneFeatureInTree_arg);
    // Input: coreFeaturesFName if !allTimeZero
    //        large: files in featureDir are grouped into subdirectories named str2hash_class(<file name>)
    // Invokes: loadPhylFile(), Genome::initDir(), setLenGlobal(), setCore(), Threads
  FeatureTree (const string &treeFName,
      				 const string &genomesListFName,
  	           bool preferGain_arg);
  	// features.size() = 1
private:
  void processBatch (const VectorPtr<Phyl> &phyls,
                     const Vector<size_t> &featureBatch,
                     float &featureLen);
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
public:
	void qc () const override;
	void saveText (ostream& os) const override;
	  // os <<: if !allTimeZero
	  //          if !timeOptimWhole() then timeOptimFrac
	  //          lambda0
	  //          time_init
	  //        topology
	  // Requires: coreSynced

	
  void deleteLeaf (TreeNode* leaf,
                   bool deleteTransientAncestor) override;
    // Requires: features.empty()

	void printInput (ostream& os) const;
	void progInternal (const string& step = noString)
	  { if (prog_. get ()) 
	  	  (*prog_) (step);
	  }
	  // Requires: after Progress::Start
	bool featuresExist () const
	  { return ! features. empty (); }
	size_t getTotalFeatures () const
	  { return commonCore. size () + features. size () + globalSingletonsSize; }  
  float getLength_min ()  
    { if (! emptySuperRoot)
        return 0.0;
      // Needed for E. coli pan ??
      if (allTimeZero)  
    		return (float) getTotalFeatures ();
      // Simplistic model: all Genome::core's are the same, all Species::time = 0 except root
    	const Species* s = static_cast <const Species*> (root);
    	return (float) getTotalFeatures () * s->feature2weight (true,  false)
    	     /*+ commonMissings   * s->feature2weight (false, false)*/;
    }
  // Sankoff algorithm
private:
	float feature2treeLength (size_t featureIndex) const;
public:
  float getLength () const;

  // OPTIMIZATION 
  // Phyl::parent2core[]
  void setLenGlobal ();
  // Phyl::core[]
  void setCore ()
    { coreSynced = true;
      const_static_cast <Species*> (root) -> setCore (); 
    }
  // Phyl::{time,weight[][]}
  Real getRootTime () const
    { return emptySuperRoot ? inf : 0; }
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
  const Species* findRoot ();
	  // Return: arg min_{Species::middleCcoreSize()} or nullptr if current root is best
    // Idempotent (prove ??)
    // Input: Species::middleCoreSize
	  // Invokes: setLeaves();
	  // Requires: allTimeZero
	string changeRoot ();
    // Return: new root LCA name in the old tree; "" if new root = old root
	  // Output: topology, superRootCore
	  // Invokes: findRoot()
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
           	     ?   parent2core [true]  [featureIndex]. treeLen
           	       < parent2core [false] [featureIndex]. treeLen
           	           // superRoot should be minimal
  	    	       : superRootCore [featureIndex];
    }
  void setStats ();
    // Output: Stats
    // Invokes: Phyl::getNeighborDistance()
private:
  void clearStats ();
    // Output: Stats
  void setFeatureStats (Feature &f,
                        size_t coreIndex,
                        const Phyl* phyl);
	void tryChange (Change* ch,
	                const Change* &bestChange);
    // Update: bestChange: positive(improvement)
    // Invokes: Change::{apply(),restore()}
public:	
	void delayDelete (Species* s);
  void finishChanges ();
    // Invokes: deleteTimeZero()
private:
  size_t deleteTimeZero ();
    // Return: # Fossil's delete'd
public:
		
  const Genome* findGenome (const string &genomeId) const;
    // Return: May be nullptr
  size_t findFeature (const Feature::Id &featureName) const;
    // Return: may be no_index
};



}



#endif


