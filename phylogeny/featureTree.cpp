// featureTree.cpp

#undef NDEBUG
#include "../common.inc"

#include "featureTree.hpp"

#include "../dm/optim.hpp"



namespace FeatureTree_sp
{



#define NOMINAL_OTHER "_OTHER_"



// Feature

void Feature::qc () const
{
  if (! qc_on)
    return;

	ASSERT (! name. empty ());
	ASSERT (! contains (name, ";"));
//IMPLY (isGene, geneId ());
	ASSERT (realGains () <= genomes);
	IMPLY (genomes, realGains ());
	IMPLY (genomes <= 1, losses. empty ());
}



bool Feature::statLess (const Feature& a,
                        const Feature& b)
{ 
	LESS_PART (a, b, rootGain);
	LESS_PART (a, b, gains);
	LESS_PART (a, b, losses);
	return false;
}



bool Feature::singleton (const Id &featureId)
{ 
  return isRight (featureId, ":" NOMINAL_OTHER); 
}




// Phyl	

Phyl::Phyl (FeatureTree &tree,
	          Species* parent_arg)
: TreeNode (tree, parent_arg)  // FeatureTree must be declared
, index_init (tree. nodeIndex_max++)
{
	for (const bool i : {false, true})
  	for (const bool j : {false, true})
		  weight [i] [j] = NaN;
}



void Phyl::init () 
{ 
	const size_t n = getFeatureTree (). features. size ();
	for (const bool parentCore : {false, true})
	{
	  ASSERT (parent2core [parentCore]. empty ());
		parent2core [parentCore]. resize (n); 
	}
	
	ASSERT (core. empty ());
	core. resize (n, false); 
}



void Phyl::qc () const
{ 
  if (! qc_on)
    return;
	TreeNode::qc ();

	const size_t n = getFeatureTree (). features. size ();
	for (const bool parentCore : {false, true})
  	ASSERT (parent2core [parentCore]. size () == n);
	ASSERT (core. size () == n);
	for (const bool i : {false, true})
  	for (const bool j : {false, true})
	    ASSERT (weight [i] [j] >= 0);

	FFOR (size_t, i, core. size ())
	{
 	  ASSERT (parent2core [false] [i]. core <= parent2core [true] [i]. core); 
 	//IMPLY (getFeatureTree (). allTimeZero, fabs (parent2core [false] [i]. treeLen - parent2core [true] [i]. treeLen) <= 1.001); ??
 	  ASSERT (! (   parent2core [false] [i]. core == UBOOL
 	             && parent2core [true]  [i]. core == UBOOL
 	            )
 	         );
 	  if (getFeatureTree (). coreSynced)
 	  {
	    IMPLY (   parent2core [false] [i]. core  != UBOOL
	    	     && parent2core [false] [i]. core == parent2core [true] [i]. core,
	    	     core [i] == parent2core [false] [i]. core
	    	    );
	  }
 	}

	if (verbose ())
	{
		Real weight_ [2] [2];
  	for (const bool i : {false, true})
    	for (const bool j : {false, true})
		    weight_ [i] [j] = weight [i] [j];
		const_cast <Phyl*> (this) -> setWeight ();
  	for (const bool i : {false, true})
    	for (const bool j : {false, true})
		    if (! eqReal (weight_ [i] [j], weight [i] [j], 1e-3))
		    {
		    	ONumber oNum (cout, 5, false);
		      cout << getName () << ": w[" << (int) i << "][" << (int) j << "]=" << weight_ [i] [j] << " " << weight [i] [j] << endl;
		    	ERROR;
		    }
	}
}



void Phyl::saveContent (ostream& os) const
{ 
	os << "C=" << getCoreSize ();	

	if (! getFeatureTree (). allTimeZero)
	{
  	{ 
      const Real r = weight [true] [true] - weight [false] [true];
  		ONumber oNum (os, 2, false);  // PAR
      if (r >= 0)
        os << "  dW[][true]=" << r;  // May happen if getParent() == getFeatureTree().root and time >> 0
    }
	}
	os << "  dC="
	   << '+' << getCoreChange (true)
	   << '-' << getCoreChange (false);
	 
	const size_t reportFeature = getFeatureTree (). reportFeature;
	if (reportFeature != NO_INDEX)
	{
	  os << "  " << getFeatureTree (). features [reportFeature]. name 
	     << ": core=" << core [reportFeature];
	  for (const bool parentCore : {false, true})
	    os << " coreEval[" << parentCore << "]=" << parent2core [parentCore] [reportFeature]. core 
	       << " len["      << parentCore << "]=" << parent2core [parentCore] [reportFeature]. treeLen;	        
	}
}



const FeatureTree& Phyl::getFeatureTree () const
{
  return static_cast <const FeatureTree&> (getTree ());
}



float Phyl::feature2weight (size_t /*featureIndex ??*/,
	                          bool thisCore,
	                          bool parentCore) const
{ 
	const float w = /*getFeatureTree (). features [featureIndex].*/ weight [thisCore] [parentCore]; 
	ASSERT (w >= 0);
	return w;
}



bool Phyl::feature2parentCore (size_t featureIndex) const
{ 
  const Species* s = static_cast <const Species*> (getParent ());
  return s ? s->core [featureIndex] : getFeatureTree (). getSuperRootCore (featureIndex);
}



bool Phyl::feature2core (size_t featureIndex) const
{ 
  const bool parentCore = feature2parentCore (featureIndex);
	const ebool c = parent2core [parentCore] [featureIndex]. core;
  return c == UBOOL ? ! getFeatureTree (). preferGain : (bool) c; 
}



size_t Phyl::getCoreSize () const
{
  ASSERT (getFeatureTree (). coreSynced);   

  size_t n = getFeatureTree (). commonCore. size ();
  FFOR (size_t, i, getFeatureTree (). features. size ())
    if (core [i])
    	n++;
  return n;
}



size_t Phyl::getCoreChange (bool gain) const
{ 
	ASSERT (getFeatureTree (). coreSynced);

	size_t d = 0;
	FFOR (size_t, i, getFeatureTree (). features. size ())
	  if (   core [i]               == gain
	  	  && feature2parentCore (i) != gain
	  	 )
	  	d++;
	return d;
}



float Phyl::getDistance () const
{ 
	ASSERT (getFeatureTree (). coreSynced);

	float d = getPooledDistance ();
	FFOR (size_t, i, getFeatureTree (). features. size ())
  	d += feature2weight (i, core [i], feature2parentCore (i));
  ASSERT (d >= 0);
#ifndef NDEBUG
  if (d == INF)
  {
  	cout << getName () << "  " << endl;
  	saveContent (cout);
  	cout << endl;
  	for (const bool i : {false, true})
    	for (const bool j : {false, true})
	  	  cout << "weight [" << (int) i << "][" << (int) j << "] = " << weight [i] [j] << endl;
  	ERROR;
  }
#endif
  
	return d;
}



float Phyl::getSubtreeLength () const
{ 
	float d = getDistance ();
	for (const DiGraph::Arc* arc : arcs [false])
	  d += static_cast <Phyl*> (arc->node [false]) -> getSubtreeLength ();
	return d;
}	  



void Phyl::getParent2corePooled (size_t parent2corePooled [2/*thisCore*/] [2/*parentCore*/]) const
{
	for (const bool thisCore : {false, true})
  	for (const bool parentCore : {false, true})
      parent2corePooled [thisCore] [parentCore] = 0;
  parent2corePooled [true] [getParent () ? true : ! getFeatureTree (). emptySuperRoot] = getFeatureTree (). commonCore. size ();
  parent2corePooled [false] [false] = getFeatureTree (). globalSingletonsSize;
}



float Phyl::getPooledDistance () const
{ 
  size_t parent2corePooled [2/*thisCore*/] [2/*parentCore*/];
  getParent2corePooled (parent2corePooled);

  float d = 0;
	for (const bool thisCore : {false, true})
  	for (const bool parentCore : {false, true})
      d += (float) multiplyLog ((Real) parent2corePooled [thisCore] [parentCore], feature2weight (thisCore, parentCore));
  
  if (! (d >= 0))
  {
    cout << getName () << " " << d << endl;
  	for (const bool thisCore : {false, true})
    	for (const bool parentCore : {false, true})
        cout << "parent2corePooled [" << (int) thisCore << "] [" << (int) parentCore << "] = " << parent2corePooled [thisCore] [parentCore] << endl;
  	for (const bool thisCore : {false, true})
    	for (const bool parentCore : {false, true})
        cout << "feature2weight ["    << (int) thisCore << "] [" << (int) parentCore << "] = " << feature2weight (thisCore, parentCore) << endl;
    ERROR;
  }
  return d;
}



float Phyl::getNeighborDistance () const
{
	float d = getDistance ();
	for (const DiGraph::Arc* arc: arcs [false])
	  d += static_cast <Phyl*> (arc->node [false]) -> getDistance ();
	return d;
}



bool Phyl::badNeighborDistance (Real &neighborDistance_stnd,
                                Real &depth_stnd) const
{
  neighborDistance_stnd = NaN;
  depth_stnd = NaN;
  const Normal* distDistr = getFeatureTree (). distDistr. get ();
  if (! distDistr)
    return false;
  neighborDistance_stnd = distDistr->stnd (getNeighborDistance ());
  depth_stnd = getFeatureTree (). depthDistr->stnd ((Real) getTopologicalDepth ());
  // PAR
  return    neighborDistance_stnd >= 3.5  
         && depth_stnd > -2.0  
         && this != getFeatureTree (). root;
}



void Phyl::getBadNodes (VectorPtr<Phyl> &badNodes,
                        bool parentBad) const
{
	Real neighborDistance_stnd, depth_stnd;
	const bool bad = badNeighborDistance (neighborDistance_stnd, depth_stnd);
	if (bad || parentBad)
	  badNodes << this;
  
	for (const DiGraph::Arc* arc : arcs [false])
	  static_cast <Phyl*> (arc->node [false]) -> getBadNodes (badNodes, bad);
}



#if 0
Real Phyl::getFeatureDistance (size_t featureIndex) const
{ 
	ASSERT (getFeatureTree (). coreSynced);
	return feature2weight (featureIndex, core [featureIndex], feature2parentCore (featureIndex));
}



Real Phyl::getFeatureSubtreeLength (size_t featureIndex) const
{ 
	Real d = getFeatureDistance (featureIndex);
	for (const DiGraph::Arc* arc : arcs [false])
	  d += static_cast <Phyl*> (arc->node [false]) -> getFeatureSubtreeLength (featureIndex);
	return d;
}	  
#endif



void Phyl::assignFeature (size_t featureIndex)
{ 
  const_cast <FeatureTree&> (getFeatureTree ()). coreSynced = false;

	float childrenCoreDistance [2/*bool thisCore*/];
	for (const bool thisCore : {false, true})
	{
    childrenCoreDistance [thisCore] = 0;
		for (const DiGraph::Arc* arc : arcs [false])
			childrenCoreDistance [thisCore] += static_cast <Phyl*> (arc->node [false]) -> parent2core [thisCore] [featureIndex]. treeLen;
  }
  
	for (const bool parentCore : {false, true})
	{
		float distance [2/*bool thisCore*/];
  	for (const bool thisCore : {false, true})
		  distance [thisCore] =   childrenCoreDistance [thisCore] 
		                        + feature2weight (featureIndex, thisCore, parentCore);
		const ebool featureCore = eqReal (distance [false], distance [true])
		                            ? UBOOL   
		                            : (ebool) (distance [true] < distance [false]);
		setCoreEval (featureIndex, parentCore, CoreEval (min (distance [false], distance [true]), featureCore));
	}
}



void Phyl::assignFeatures ()
{ 
  FFOR (size_t, i, parent2core [false]. size ())
    assignFeature (i);
}



void Phyl::assignFeaturesDown ()
{ 
  // Threads ??
	for (DiGraph::Arc* arc : arcs [false])
		static_cast <Phyl*> (arc->node [false]) -> assignFeaturesDown ();
  assignFeatures ();
}




// Species

void Species::Movement::undo (Species* s) const
{ 
	ASSERT (s);
	s->parent2core [parentCore] [featureIndex] = from;
}



Species::Species (FeatureTree &tree,
 	                Fossil* parent_arg,
 	                const string &id_arg,
	                Real time_arg)
: Phyl (tree, parent_arg)  // FeatureTree must be declared
, id (id_arg)
, time (tree. allTimeZero ? NaN : time_arg)
{
	for (const bool i : {false, true})
  	for (const bool j : {false, true})
		  weight_old [i] [j] = NaN;
}



void Species::qc () const
{ 
  if (! qc_on)
    return;
	Phyl::qc ();
	
	IMPLY (! isNan (time), time >= 0);	
	IMPLY (getFeatureTree (). allTimeZero, isNan (time)); 
  ASSERT (pooledSubtreeDistance >= 0);	
	IMPLY (! movementsOn, movements. empty ());
}



void Species::saveContent (ostream& os) const
{ 
	if (! getFeatureTree (). allTimeZero)
  {
		ONumber oNum (os, 6, true);  // PAR
		os << "t=" << time << "  ";
	}
	
	Phyl::saveContent (os);
}



double Species::getParentDistance () const
{ 
  return getParent () 
           ? getFeatureTree (). allTimeZero 
             ? (double) getCoreChange ()
             : time 
           : -1; 
}



void Species::setWeight ()
{
  if (   getFeatureTree (). allTimeZero 
      || ! getFeatureTree (). featuresExist ()
      || isNan (getFeatureTree (). lambda0)
     )
  	for (const bool parentCore : {false, true})
	  {
	  	// Parsimony method
	    weight [  parentCore] [parentCore] = 0;
	    weight [! parentCore] [parentCore] = 1;  
		}
  else
  {
	  ASSERT (time >= 0);
	 	const Prob a = 1 - exp (- time);
	 	ASSERT (isProb (a));
  	for (const bool parentCore : {false, true})
	  {
	  	const Prob changeProb = negateProb (getFeatureTree (). lambda0, parentCore) * a;
			ASSERT (isProb (changeProb));
    	for (const bool thisCore : {false, true})
		  {
		  	const Prob p = negateProb (changeProb, parentCore == thisCore); 
		  	const Real w = - log (p);
		  	ASSERT (w >= 0);
		    weight [thisCore] [parentCore] = (float) w;
		  }
		}
	}
}



void Species::setCore ()
{
	FFOR (size_t, i, core. size ())
	  core [i] = feature2core (i);
	for (DiGraph::Arc* arc : arcs [false])
	  static_cast <Phyl*> (arc->node [false]) -> setCore ();
}



namespace {

struct TimeFunc : Func1
{
	const FeatureTree& tree;
    // For lambda0
  Real parent2core [2/*thisCore*/] [2/*parentCore*/];
  
//Prob b [2];


  TimeFunc (const FeatureTree& tree_arg,
            const size_t parent2core_arg [2/*thisCore*/] [2/*parentCore*/] /*,
            size_t nodes*/)
		: tree (tree_arg)
		{
      ASSERT (! tree. allTimeZero);
    //ASSERT (nodes);
    	for (const bool i : {false, true})
      	for (const bool j : {false, true})
   	      parent2core [i] [j] = (Real) parent2core_arg [i] [j];
   	//parent2core [false] [false] +=        tree->commonMissings      * (Real) nodes;
	  //parent2core [true] [true]   += (Real) tree->commonCore. size () /* * (Real) nodes*/;
    }

  
  Real f (Real x) 
  { 
    ASSERT (x >= 0);
    
  	const Prob a = 1 - exp (- x);
  	
  	Real lhs = 0;
  	for (const bool c : {false, true})
  		if (const Prob denomin = a)
  	    lhs += parent2core [! c] [c] / denomin;
  	  else
  	  	lhs = INF;
  	ASSERT (lhs >= 0);
  	  
  	Real rhs = 0;
  	for (const bool c : {false, true})
  		if (const Real denomin = 1 / negateProb (tree. lambda0, c) - a)
  	    rhs += parent2core [c] [c] / denomin;
  	  else
  	  	rhs = INF;
  	ASSERT (rhs >= 0);
  	  
  	const Real r = rhs - lhs;
  	ASSERT (! isNan (r));
  	return r;
  }  
};



Real coreChange2time (const FeatureTree& tree,
				              const size_t parent2core [2/*thisCore*/] [2/*parentCore*/] /*,
				              size_t nodes*/)
{
	ASSERT (! tree. allTimeZero);
//ASSERT (nodes);

	TimeFunc f (tree, parent2core /*, nodes*/);
	
	ASSERT (! positive (f. f (0)));
	if (nullReal (f. f (0)))
		return 0;

	if (! positive (f. f (INF)))
		return INF;
  
  Unverbose unv1;
  Unverbose unv2;
  const Real r = f. findZeroPositive (1.0, 1e-6);  // PAR
  ASSERT (r >= 0);

  return r;
}

}



Real Species::getTime () const
{
  if (getFeatureTree (). allTimeZero)
  	return NaN;

  if (! getParent ())
  	return getFeatureTree (). getRootTime ();;  
  
  if (! getFeatureTree (). timeOptimFrac)  	
    return getFeatureTree (). time_init;  

  size_t parent2core_ [2/*thisCore*/] [2/*parentCore*/];
  getParent2corePooled (parent2core_);	  
	FFOR (size_t, i, getFeatureTree (). features. size ())
	  parent2core_ [core [i]] [feature2parentCore (i)] ++;
	  
  const Real t = coreChange2time (getFeatureTree (), parent2core_ /*, 1*/);
  ASSERT (t >= 0);
  
  return exp (convexCombination (getFeatureTree (). timeOptimFrac, log (t), log (getFeatureTree (). time_init)));
}



void Species::assignTime ()
{ 
	IMPLY (! getFeatureTree (). allTimeZero, ! isNan (time));
	
	time_old = time;
	for (const bool i : {false, true})
  	for (const bool j : {false, true})
		  weight_old [i] [j] = weight [i] [j];

  setTimeWeight ();
}



void Species::restoreTime ()
{
	IMPLY (! getFeatureTree (). allTimeZero, ! isNan (time_old));
	
	time = time_old;
	for (const bool i : {false, true})
  	for (const bool j : {false, true})
		  weight [i] [j] = weight_old [i] [j];
	
	commitTime ();
}



void Species::assignFeaturesUp (const Fossil* toParentExcluding)
{
	Species* s = this;
	while (s != toParentExcluding)
	{
		ASSERT (s);
	  s->rememberFeatures ();
	  s->assignFeatures ();
	  s = const_static_cast <Species*> (s->getParent ());
	}
}



void Species::restoreFeaturesUp (const Fossil* toParentExcluding)
{
	Species* s = this;
	while (s != toParentExcluding)
	{
		ASSERT (s);
		s->restoreFeatures ();
	  s = const_static_cast <Species*> (s->getParent ());
	}
}



void Species::commitFeaturesUp (const Fossil* toParentExcluding)
{
	Species* s = this;
	while (s != toParentExcluding)
	{
		ASSERT (s);
		s->commitFeatures ();
	  s = const_static_cast <Species*> (s->getParent ());
	}
}



void Species::rememberAllFeatures ()
{
  rememberFeatures ();

  movements. reserve (2 * parent2core [false]. size ());
	for (const bool parentCore : {false, true})
    FFOR (size_t, featureIndex, parent2core [parentCore]. size ())
      movements << Movement (parentCore, featureIndex, parent2core [parentCore] [featureIndex]);
}



void Species::rememberFeatures ()
{
	ASSERT (! movementsOn);
	ASSERT (movements. empty ());
	movementsOn = true;
	movements. reserve (128);  // PAR
	
	ASSERT (isNan (pooledSubtreeDistance_old));
	pooledSubtreeDistance_old = pooledSubtreeDistance;
}



void Species::restoreFeatures ()
{
	/*C++11: CONST_*/ITER_REV (Vector<Movement>, it, movements)
	  it->undo (this);
	  
	ASSERT (! isNan (pooledSubtreeDistance_old));
	pooledSubtreeDistance = pooledSubtreeDistance_old;
	
	commitFeatures ();
}



void Species::commitFeatures ()
{
	ASSERT (movementsOn);
	movements. clear ();
	movementsOn = false;
	pooledSubtreeDistance_old = NaN;
}




// Fossil

void Fossil::qc () const
{ 
  if (! qc_on)
    return;
	Species::qc ();

	for (const DiGraph::Arc* arc : arcs [false])
	{
		const Phyl* sub = static_cast <Phyl*> (arc->node [false]);
		ASSERT (sub->asFossil () || sub->asStrain ());
	}	
	ASSERT (! isLeaf ());
}



float Fossil::getPooledSubtreeDistance () const 
{ 
  float d = getPooledDistance ();
	for (const DiGraph::Arc* arc : arcs [false])
	  d += static_cast <Species*> (arc->node [false]) -> pooledSubtreeDistance;
  return d; 
}



void Fossil::setId (uint &id_arg)
{
  id = toString (id_arg);
  id_arg++;
	for (const DiGraph::Arc* arc : arcs [false])
		if (const Fossil* f = static_cast <Phyl*> (arc->node [false]) -> asFossil ())
		  const_cast <Fossil*> (f) -> setId (id_arg);
}





// Strain

void Strain::saveContent (ostream& os) const
{
  Species::saveContent (os);

	os << "  S=" << (singletonsInCore ? getGenome () -> singletons. size () : 0);
  if (const size_t s = orphans ())
    os << "  orphans=" << s;
}



void Strain::qc () const
{
  if (! qc_on)
    return;
  Species::qc ();
  
  ASSERT (arcs [false]. size () == 1);  // Transient
  ASSERT (getGenome ());
	
	IMPLY (getFeatureTree (). allTimeZero, singletonsInCore);
	
	const Species* p = static_cast <const Species*> (getParent ());
	const Genome* g = getGenome ();
	FFOR (size_t, i, core. size ())
	  IMPLY (p->core [i] == g->core [i], p->core [i] == core [i]);  // <= (time <= time_max) ??
}



void Strain::getParent2corePooled (size_t parent2corePooled [2/*thisCore*/] [2/*parentCore*/]) const
{
  Species::getParent2corePooled (parent2corePooled);
  
  if (singletonsInCore)
  {
    const Genome* g = getGenome ();
    ASSERT (parent2corePooled [false] [false] >= g->singletons. size ());
    parent2corePooled [false] [false] -= g->singletons. size ();
    parent2corePooled [true] [false]  += g->singletons. size ();
  }
}



void Strain::assignFeatures ()
{ 
  // singletonsInCore
	float distance [2/*bool thisCore*/];
	const Genome* g = getGenome ();
	for (const bool thisCore : {false, true})
	  distance [thisCore] =   g->feature2weight (true,     thisCore)
	                        +    feature2weight (thisCore, false);
  singletonsInCore = leReal (distance [true], distance [false]);  

  Species::assignFeatures ();
}



float Strain::getPooledSubtreeDistance () const 
{ 
  return getPooledDistance () + getGenome () -> getPooledDistance (); 
}



size_t Strain::orphans () const
{
  ASSERT (getFeatureTree (). coreSynced);
  
  size_t n = 0;
  if (const Fossil* p = static_cast <const Fossil*> (getParent ()))
  {
    const Genome* g = getGenome ();
    FFOR (size_t, i, getFeatureTree (). features. size ())
      if (   p->core [i] == g->core [i] 
          && p->core [i] !=    core [i]
         )
        n++;
  }
  
  return n;
}




// Genome

Genome::Genome (FeatureTree &tree,
	              Strain* parent_arg,
				        const string &id_arg)
: Phyl (tree, parent_arg)
, id (id_arg)
{ 
	ASSERT (parent_arg);
	ASSERT (! id. empty ()); 
	ASSERT (! contains (id, ' '));
}



string Genome::geneLineFormat ()
{ 
  return "{{<Boolean> [<optional (0|1)>]} | {<nominal>:<value>} \\n}*, where <value> = " NOMINAL_OTHER " means singleton Boolean attributes for all values of the nominal attribute"; 
}



void Genome::initDir (const string &geneDir,
                      bool nominalSingletonIsOptional)
{
  ASSERT (! id. empty ());
  ASSERT (coreSet. empty ());
  ASSERT (! coreNonSingletons);
 
  if (! geneDir. empty ())
  {
    LineInput f (geneDir + "/" + id);
    // coreSet
    while (f. nextLine ())
    {
      Feature::Id geneName;
      bool optional = false;
      trim (f. line);
      replace (f. line, '\t', ' ');
      if (contains (f. line, " : "))
        throw runtime_error ("':' cannot have neighboring spaces");
      if (contains (f. line, ':'))
      {
        geneName = f. line;
        if (Feature::singleton (geneName) && nominalSingletonIsOptional)
          continue;
      }
      else if (   isRight (f. line, " 0")
               || isRight (f. line, " 1")
              )
      {
        optional = isRight (f. line, " 1");
        geneName = f. line. substr (0, f. line. size () - 2);
      }
      else
        geneName = f. line;
      trim (geneName);
      ASSERT (! geneName. empty ());
      if (contains (coreSet, geneName))
        throw runtime_error ("Gene " + strQuote (geneName) + " is duplicated in genome " + getName ());
      coreSet [geneName] = optional;
    }
  }
  
  coreNonSingletons = coreSet. size ();  // Includes singletons and common core
}



void Genome::coreSet2nominals () 
{
  ASSERT (nominals. empty ());
	for (const auto& it : coreSet)
	{
	  const string& geneName = it. first;
	  const size_t pos = geneName. find (':');
	  if (pos == string::npos)
	    continue;
	  if (it. second)
	    throw runtime_error ("Optional nominal attribute " + geneName);
	  const string attrName (geneName. substr (0, pos));
	  const string value    (geneName. substr (pos + 1));
	  if (attrName. empty ())
	    throw runtime_error ("Empty nominal attribute name: " + geneName);
	  if (value. empty ())
	    throw runtime_error ("Empty nominal attribute value: " + geneName);
	  if (! nominals. addUnique (attrName))
	    throw runtime_error ("Duplicate nominal attribute: " + attrName);
	  if (! Feature::singleton (geneName))
	    const_cast <FeatureTree&> (getFeatureTree ()). nominal2values [attrName] << value;
	}
}



void Genome::nominals2coreSet ()
{
  for (const auto& it : getFeatureTree (). nominal2values)
    if (! nominals. contains (it. first))
      for (const string& value : it. second)
      {
        const Feature::Id s (it. first + ":" + value);
        ASSERT (! contains (coreSet, s));
        coreSet [s] = true;
        coreNonSingletons++;
      }
        
  nominals. clear ();
}



void Genome::init (const map <Feature::Id, size_t/*index*/> &feature2index)
{
	IMPLY ( ! coreSet. empty (), getFeatureTree (). featuresExist ());

	Phyl::init ();

	ASSERT (optionalCore. empty ());
	optionalCore. resize (core. size (), false); 
	for (const auto& it : coreSet)
	{
		ASSERT (contains (feature2index, it. first));
		const size_t featureIndex = feature2index. at (it. first);
		for (const bool parentCore : {false, true})
	    parent2core [parentCore] [featureIndex]. core = ETRUE;
	  core [featureIndex] = true;
    optionalCore [featureIndex] = it. second;
	}

	coreSet. clear ();
}



void Genome::qc () const
{ 
  if (! qc_on)
    return;
	Phyl::qc ();
	  
	ASSERT (optionalCore. size () == core. size ());

  ASSERT (getParent ());
  ASSERT (getStrain ());
	ASSERT (isLeaf ());
	
//ASSERT (getName (). substr (1) == getStrain () -> getName (). substr (1));  

	ASSERT (lessReal (weight [true] [true], weight [false] [true]));
	ASSERT (lessReal (weight [false] [false], weight [true] [false]));
//ASSERT (leReal (weight [false] [false], weight [true] [true]));

  ASSERT (! id. empty ());
  ASSERT (! contains (id, ' '));
  ASSERT (coreNonSingletons >= getFeatureTree (). commonCore. size ());
	ASSERT (nominals. empty ());
	ASSERT (coreSet. empty ());
	ASSERT (singletons. searchSorted);
/*
	if (getFeatureTree (). coreSynced && coreNonSingletons < getCoreSize ())
	{
	  print (cout);
	  cout << endl;
	  cout << coreNonSingletons << ' ' << getCoreSize () << ' ' << singletons. size () << endl;
    FFOR (size_t, i, getFeatureTree (). features. size ())
      if (core [i])
      	cout << getFeatureTree (). features [i]. name << endl;
	  ERROR;
	}
*/
	
	ASSERT (! singletons. intersects (getFeatureTree (). commonCore));
	
	for (const string& attrName : nominals)
	  ASSERT (contains (getFeatureTree (). nominal2values, attrName));
}



void Genome::saveContent (ostream& os) const
{ 
  Phyl::saveContent (os);  

	os << "  S=" << (getStrain () -> singletonsInCore ? 0 : singletons. size ()) /*<< '+' << oddCdss*/;    
/*
	if (const size_t opt = coreNonSingletons - getCoreSize ())
	  os << "  OptC=" << opt;
*/
#if 0
	os << "  " << project_id  // << ':' << tax_id 
	   << "  " << taxName;
	os << "  L50=" << L50;
  if (! sequencer. empty ())
  	os << "  " << sequencer;
	if (pubmed)
		os << "  PMID:" << pubmed;
	const string extra (getNameExtra ());
	if (! extra. empty ())
	  os << extra;
  if (! phylogeneticClass. empty ())
  	os << "  " << phylogeneticClass;
#endif
}



void Genome::setWeight ()
{
  if (   getFeatureTree (). allTimeZero 
      || ! getFeatureTree (). featuresExist ()
      || isNan (getFeatureTree (). lambda0)
     )
  	for (const bool parentCore : {false, true})
	  {
	  	// Parsimony method
	    weight [  parentCore] [parentCore] = 0;
	    weight [! parentCore] [parentCore] = INF;  
		}
  else
  {
    const Real n_1 = (Real) getGenes ();                               // ~ # Present genes in getStrain()
    const Real n_0 = (Real) getFeatureTree (). getTotalGenes () - n_1;  // ~ # Absent  genes in getStrain()
    ASSERT (n_0 > 0);
  #if 0
    ASSERT (L50);
    const Real missedGenes = 30.5 + 0.63 * L50 + 0.48 * (Real) singletons. size ();  // PAR: estimated for Salmonella  ??
    const Prob annotError = missedGenes / n_1;  // was: 1 - 0.995
    const Prob misannotError = missedGenes * 0.3 / n_0;  // PAR: estimated for Salmonella ??
  #else
    const Prob annotError = 0.005;
    const Real missedGenes = annotError * n_1;
    const Prob misannotError = missedGenes * 0.3 / n_0;  // PAR: estimated for Salmonella ??
  #endif
    ASSERT (isProb (annotError));
    ASSERT (isProb (misannotError));
    weight [true]  [true]  = - (float) log (1 - annotError);
    weight [false] [true]  = - (float) log (annotError);
    weight [true]  [false] = - (float) log (misannotError);
    weight [false] [false] = - (float) log (1 - misannotError);
	}
}



void Genome::getParent2corePooled (size_t parent2corePooled [2/*thisCore*/] [2/*parentCore*/]) const
{
  Phyl::getParent2corePooled (parent2corePooled);

  ASSERT (parent2corePooled [false] [false] >= singletons. size ());
  parent2corePooled [false] [false] -= singletons. size ();
  if (getStrain () -> singletonsInCore)
    parent2corePooled [true] [true]  += singletons. size ();
  else
    parent2corePooled [true] [false] += singletons. size ();
}



void Genome::setCore ()
{
	FFOR (size_t, i, core. size ())
    if (optionalCore [i])
	    core [i] = feature2core (i);
}



void Genome::assignFeature (size_t featureIndex)
{
  if (optionalCore [featureIndex])
    Phyl::assignFeature (featureIndex);
  else
  	for (const bool parentCore : {false, true})
    {
    	const ebool core_old = parent2core [parentCore] [featureIndex]. core;
    	ASSERT (core_old != UBOOL);
    	setCoreEval (featureIndex, parentCore, CoreEval (feature2weight (featureIndex, core_old, parentCore), core_old));
    }
}



void Genome::getSingletons (Set<Feature::Id> &globalSingletons,
                            Set<Feature::Id> &nonSingletons) const
{ 
	if (getFeatureTree (). featuresExist () && coreSet. empty ())
	{
	  cout << getName () << endl;  
	  ERROR;
	}
	
	for (const auto& it : coreSet)
	{
	/*if (it. second)
	    continue; ?? */
	  const Feature::Id& featureId = it. first;
	  if (Feature::singleton (featureId))
	  {
	    ASSERT (! globalSingletons. contains (featureId));
	    ASSERT (! nonSingletons.    contains (featureId));
	  }
	  else
	  {
  	  if (nonSingletons. contains (featureId))
  	  	;
  	  else if (globalSingletons. contains (featureId))
  	  	setMove (& globalSingletons, & nonSingletons, featureId);
  	  else
  	  	EXEC_ASSERT (globalSingletons. insert (featureId). second);
  	}
	}
}



size_t Genome::setSingletons (const Set<Feature::Id> &globalSingletons)
{ 
  ASSERT (singletons. empty ());
  
  Vector<Feature::Id> featureSingletons;
  for (Iter<CoreSet> iter (coreSet); iter. next (); )
  {
	  const auto& mapIt = *iter;
		if (Feature::singleton (mapIt. first))
		{
		  if (! mapIt. second)
  		  featureSingletons << mapIt. first;
			iter. erase ();
		}
	}
		
	for (const auto& it : coreSet)
	  if (globalSingletons. contains (it. first))
	    singletons << it. first;
		
  for (const Feature::Id& f : singletons)
    coreSet. erase (f);
    
  singletons << featureSingletons;  // ??

  singletons. sort ();
  singletons. uniq ();
  
  ASSERT (coreNonSingletons >= singletons. size ());
  coreNonSingletons -= singletons. size ();

  return featureSingletons. size ();
}




// Change

Change::Change (const FeatureTree &tree_arg,
                istream &is)
: tree (tree_arg)
{ 
  size_t fromIndex; 
  is >> fromIndex >> improvement;
	from = static_cast <const Phyl*> (tree. nodes. at (fromIndex)) -> asSpecies ();
}



Change::~Change ()
{
  ASSERT (status != eApplied);
}



void Change::qc () const
{
  if (! qc_on)
    return;
	ASSERT (valid ());
	ASSERT (positive (improvement));
	ASSERT (! targets. empty ());
}



void Change::saveText (ostream& os) const
{ 
	os         << type () 
	   << '\t' << from->index_init
	   << '\t' << improvement;
}



void Change::apply ()
{
  ASSERT (status == eInit);
	ASSERT (! isolatedTransient);
	ASSERT (! isolatedTransientChild);

  apply_ ();
  improvement = tree. len * (1 + (float) tree. lenInflation) - tree. getLength ();
  status = eApplied;
}



void Change::restore ()
{
  ASSERT (status == eApplied);  
  status = eInit;
  
	ASSERT ((bool) isolatedTransient == (bool) isolatedTransientChild);

  restore_ ();
  clear ();
}



void Change::commit ()
{
  ASSERT (status == eApplied);
  status = eCommitted;
  
	ASSERT ((bool) isolatedTransient == (bool) isolatedTransientChild);

  commit_ ();
	const_cast <FeatureTree&> (tree). delayDelete (isolatedTransient);
}



bool Change::isolateTransient (Fossil* f)
{
  ASSERT (f);
  ASSERT (f->graph);
  ASSERT (! isolatedTransient);
  ASSERT (! isolatedTransientChild);
  
  isolatedTransientChild = static_cast <Species*> (f->isolateTransient ()); 
  if (isolatedTransientChild)
   	isolatedTransient = f;
	ASSERT ((bool) isolatedTransient == (bool) isolatedTransientChild);
   	
  return (bool) isolatedTransient;
}



Real Change::improvementDeinflated () const
{ 
  return improvement - tree. len * tree. lenInflation; 
}





// ChangeTo

ChangeTo::ChangeTo (const FeatureTree &tree_arg,
                    istream &is)
: Change (tree_arg, is)
{ 
  size_t toIndex; 
  is >> toIndex;
	to = static_cast <const Phyl*> (tree. nodes. at (toIndex)) -> asSpecies ();
}



void ChangeTo::saveText (ostream& os) const
{ 
	Change::saveText (os);
	os << '\t' << to->index_init;
}



bool ChangeTo::better (const Change* other) const
{
  ASSERT (other);
  ASSERT (from == other->from);
  const ChangeTo* c = other->asChangeTo ();
  ASSERT (c);
  if (to == c->to)
  {
 	  saveText (cout);
 	  cout << endl;
 	  other->saveText (cout);
 	  cout << endl;    
    ERROR;
  }
  return to->index_init < c->to->index_init;
}




// ChangeToSibling

void ChangeToSibling::apply_ ()
{
  ASSERT (! inter);
  ASSERT (! oldParent);
  ASSERT (! oldParentRepr);
  ASSERT (! arcEnd);
  ASSERT (! lca);


  oldParent = const_static_cast <Fossil*> (from->getParent ());
  ASSERT (oldParent);
  // May be: oldParent == to
  oldParentRepr = oldParent;
  arcEnd = const_static_cast <Fossil*> (to->getParent ());
    // May be nullptr
  ASSERT (from != arcEnd);

  // Topology
  auto inter_ = new Fossil (const_cast <FeatureTree&> (tree), arcEnd, string (), NaN);
  inter = inter_;
  inter_->init ();
  // PAR
  FFOR (size_t, i, inter_->core. size ())
    inter_->core [i] = (arcEnd ? arcEnd->core [i] : false) + from->core [i] + to->core [i] >= 2;
  inter_->setTimeWeight ();
  const_cast <Species*> (from) -> setParent (inter_); 
  const_cast <Species*> (to)   -> setParent (inter_);
  {
    Fossil* oldParent_parent = const_static_cast <Fossil*> (oldParent->getParent ());  // May be nullptr
    if (isolateTransient (oldParent))
    {
    	oldParent = oldParent_parent;
    	oldParentRepr = isolatedTransientChild;
      ASSERT (isolatedTransient != oldParent);
      ASSERT (oldParentRepr != oldParent);
      ASSERT (oldParentRepr != from);
    }
  }
  ASSERT (oldParentRepr);
  
  Species* fromSibling = const_cast <Species*> (static_cast <const Phyl*> (inter->getOtherChild (from)) -> asSpecies ());
  ASSERT (fromSibling);
  Tree::LcaBuffer buf;
	lca = static_cast <const Phyl*> (Tree::getLca (oldParentRepr, inter, buf)) -> asFossil ();
	ASSERT (lca);
  
  // Species::time
  const_cast <Species*> (from) -> assignTime (); 
  fromSibling->assignTime ();  
  if (isolatedTransientChild && oldParentRepr != fromSibling)
    oldParentRepr->assignTime (); 
  
  // Species::parent2core[]
  oldParentRepr->assignFeaturesUp (lca);  
  if (! oldParentRepr->descendantOf (fromSibling))
    fromSibling->assignFeaturesUp (inter);
  const_cast <Species*> (from) -> assignFeaturesUp (nullptr);
}



void ChangeToSibling::restore_ ()
{
  ASSERT (inter);
  ASSERT (oldParentRepr);
 	ASSERT (lca);
  if (isolatedTransientChild)
  { ASSERT (oldParentRepr == isolatedTransientChild); }
  else
  { ASSERT (oldParentRepr == oldParent); }
	

  Species* fromSibling = const_cast <Species*> (static_cast <const Phyl*> (inter->getOtherChild (from)) -> asSpecies ());
  ASSERT (fromSibling);

  // Species::parent2core[]
  const_cast <Species*> (from) -> restoreFeaturesUp (nullptr);
  if (! oldParentRepr->descendantOf (fromSibling))
    fromSibling->restoreFeaturesUp (inter);
  oldParentRepr->restoreFeaturesUp (lca);  
  
  // Species::time
  const_cast <Species*> (from) -> restoreTime (); 
  fromSibling->restoreTime ();  
  if (isolatedTransientChild && oldParentRepr != fromSibling)
    oldParentRepr->restoreTime (); 

  // Topology
  if (isolatedTransient)
  {
  	IMPLY (isolatedTransient == to,   isolatedTransientChild == inter->arcs [false]. back () -> node [false]);
  	IMPLY (isolatedTransient == arcEnd, isolatedTransientChild == inter);
  	// May be: isolatedTransientChild->getParent() == inter
  	isolatedTransient->attach (const_cast <FeatureTree&> (tree));
	  isolatedTransientChild->setParent (isolatedTransient);
	  isolatedTransient     ->setParent (oldParent);
	  oldParent = isolatedTransient;
  }
	ASSERT (oldParent);
  const_cast <Species*> (from) -> setParent (oldParent);  
  const_cast <Species*> (to)   -> setParent (arcEnd);
  delete inter;
}



void ChangeToSibling::commit_ ()
{
  ASSERT (inter);
  ASSERT (oldParentRepr);
 	ASSERT (lca);
  if (isolatedTransientChild)
  { ASSERT (oldParentRepr == isolatedTransientChild); }
  else
  { ASSERT (oldParentRepr == oldParent); }
  
  Species* fromSibling = const_cast <Species*> (static_cast <const Phyl*> (inter->getOtherChild (from)) -> asSpecies ());
  ASSERT (fromSibling);

  // Species::time
  const_cast <Species*> (from) -> commitTime (); 
  fromSibling->commitTime ();  
  if (isolatedTransientChild && oldParentRepr != fromSibling)
    oldParentRepr->commitTime (); 

  // Species::parent2core[]
  oldParentRepr->commitFeaturesUp (lca);
  if (! oldParentRepr->descendantOf (fromSibling))
    fromSibling->commitFeaturesUp (inter);
  const_cast <Species*> (from) -> commitFeaturesUp (nullptr);
}




// ChangeToParent

void ChangeToParent::apply_ ()
{
  ASSERT (! oldParent);
  ASSERT (! oldParentRepr);
  ASSERT (! lca);


  oldParent = const_static_cast <Fossil*> (from->getParent ());
  ASSERT (oldParent);
  oldParentRepr = oldParent;

  // Topology
  const_cast <Species*> (from) -> setParent (const_cast <Species*> (to)); 
  {
    Fossil* oldParent_parent = const_static_cast <Fossil*> (oldParent->getParent ());  // May be nullptr
    if (isolateTransient (oldParent))
    {
    	oldParent = oldParent_parent;
    	oldParentRepr = isolatedTransientChild;
      ASSERT (isolatedTransient != oldParent);
      ASSERT (oldParentRepr != oldParent);
    }
  }
  ASSERT (oldParentRepr);
  
  Tree::LcaBuffer buf;
	lca = static_cast <const Phyl*> (Tree::getLca (oldParentRepr, to, buf)) -> asFossil ();
	ASSERT (lca);

  // Species::time
  if (isolatedTransientChild)
    oldParentRepr->assignTime (); 
  const_cast <Species*> (from) -> assignTime (); 

  // Species::parent2core[]
  oldParentRepr->assignFeaturesUp (lca);  
  const_cast <Species*> (from) -> assignFeaturesUp (nullptr);
}



void ChangeToParent::restore_ ()
{
  ASSERT (oldParentRepr);
 	ASSERT (lca);
  if (isolatedTransientChild)
  { ASSERT (oldParentRepr == isolatedTransientChild); }
  else
  { ASSERT (oldParentRepr == oldParent); }
	

  // Species::parent2core[]
  const_cast <Species*> (from) -> restoreFeaturesUp (nullptr);
  oldParentRepr->restoreFeaturesUp (lca);  
  
  // Species::time
  if (isolatedTransientChild)
    oldParentRepr->restoreTime (); 
  const_cast <Species*> (from) -> restoreTime (); 

  // Topology
  if (isolatedTransient)
  {
  	ASSERT (isolatedTransient != to);
  	ASSERT (oldParent == isolatedTransientChild->getParent ());
  	isolatedTransient->attach (const_cast <FeatureTree&> (tree));
	  isolatedTransientChild->setParent (isolatedTransient);
	  isolatedTransient     ->setParent (oldParent);
	  oldParent = isolatedTransient;
  }
	ASSERT (oldParent);
  const_cast <Species*> (from) -> setParent (oldParent);  
}



void ChangeToParent::commit_ ()
{
  ASSERT (oldParentRepr);
 	ASSERT (lca);
  if (isolatedTransientChild)
  { ASSERT (oldParentRepr == isolatedTransientChild); }
  else
  { ASSERT (oldParentRepr == oldParent); }
  
  // Species::time
  if (isolatedTransientChild)
    oldParentRepr->commitTime (); 
  const_cast <Species*> (from) -> commitTime (); 

  // Species::parent2core[]
  oldParentRepr->commitFeaturesUp (lca);
  const_cast <Species*> (from) -> commitFeaturesUp (nullptr);  
}




// ChangeToUncle

void ChangeToUncle::apply_ ()
{
  Fossil* fromParent = const_static_cast <Fossil*> (from->getParent ());
  ASSERT (fromParent);
  Fossil* toParent = const_static_cast <Fossil*> (to->getParent ());
  ASSERT (toParent);

  // Topology
  const_cast <Species*> (from) -> setParent (toParent); 
  const_cast <Species*> (to)   -> setParent (fromParent); 

  // Species::time
  const_cast <Species*> (from) -> assignTime (); 
  const_cast <Species*> (to)   -> assignTime (); 

  // Species::parent2core[]
  const_cast <Species*> (from) -> assignFeaturesUp (toParent);
  const_cast <Species*> (to)   -> assignFeaturesUp (nullptr);
}



void ChangeToUncle::restore_ ()
{
  Fossil* fromParent = const_static_cast <Fossil*> (to->getParent ());
  ASSERT (fromParent);
  Fossil* toParent = const_static_cast <Fossil*> (from->getParent ());
  ASSERT (toParent);


  // Species::parent2core[]
  const_cast <Species*> (from) -> restoreFeaturesUp (toParent);
  const_cast <Species*> (to) -> restoreFeaturesUp (nullptr);
  
  // Species::time
  const_cast <Species*> (from) -> restoreTime (); 
  const_cast <Species*> (to)   -> restoreTime (); 

  // Topology
  const_cast <Species*> (from) -> setParent (fromParent); 
  const_cast <Species*> (to)   -> setParent (toParent); 
}



void ChangeToUncle::commit_ ()
{
  // Species::time
  const_cast <Species*> (from) -> commitTime (); 
  const_cast <Species*> (to)   -> commitTime (); 

  // Species::parent2core[]
  const_cast <Species*> (from) -> commitFeaturesUp (static_cast <const Fossil*> (from->getParent ()));
  const_cast <Species*> (to) -> commitFeaturesUp (nullptr);  
}




// ChangeToCousin

void ChangeToCousin::apply_ ()
{
  Fossil* fromParent = const_static_cast <Fossil*> (from->getParent ());
  ASSERT (fromParent);
  Fossil* toParent = const_static_cast <Fossil*> (to->getParent ());
  ASSERT (toParent);

  // Topology
  const_cast <Species*> (from) -> setParent (toParent); 
  const_cast <Species*> (to)   -> setParent (fromParent); 

  // Species::time
  const_cast <Species*> (from) -> assignTime (); 
  const_cast <Species*> (to)   -> assignTime (); 

  // Species::parent2core[]
  const_cast <Species*> (from) -> assignFeaturesUp (static_cast <const Fossil*> (toParent->getParent ()));
  const_cast <Species*> (to)   -> assignFeaturesUp (nullptr);
}



void ChangeToCousin::restore_ ()
{
  Fossil* fromParent = const_static_cast <Fossil*> (to->getParent ());
  ASSERT (fromParent);
  Fossil* toParent = const_static_cast <Fossil*> (from->getParent ());
  ASSERT (toParent);

  // Species::parent2core[]
  const_cast <Species*> (from) -> restoreFeaturesUp (static_cast <const Fossil*> (toParent->getParent ()));
  const_cast <Species*> (to)   -> restoreFeaturesUp (nullptr);
  
  // Species::time
  const_cast <Species*> (from) -> restoreTime (); 
  const_cast <Species*> (to)   -> restoreTime (); 

  // Topology
  const_cast <Species*> (from) -> setParent (fromParent); 
  const_cast <Species*> (to)   -> setParent (toParent); 
}



void ChangeToCousin::commit_ ()
{
  // Species::time
  const_cast <Species*> (from) -> commitTime (); 
  const_cast <Species*> (to)   -> commitTime (); 
  
  // Species::parent2core[]
  const_cast <Species*> (from) -> commitFeaturesUp (static_cast <const Fossil*> (from->getParent () -> getParent ()));
  const_cast <Species*> (to) -> commitFeaturesUp (nullptr);  
}




// ChangeRoot

void ChangeRoot::qc () const
{
  if (! qc_on)
    return;
	Change::qc ();
	
	ASSERT (tree. allTimeZero);
}



void ChangeRoot::apply_ ()
{
  ASSERT (! root_old);


  // Topology
  auto inter = new Fossil (const_cast <FeatureTree&> (tree), const_static_cast <Fossil*> (from->getParent ()), string (), NaN);
  inter->init ();
  FFOR (size_t, i, inter->core. size ())
    inter->core [i] = from->core [i];  // PAR
  inter->setTimeWeight ();
  const_cast <Species*> (from) -> setParent (inter); 
  //
	Fossil* root_old_ = const_cast <Fossil*> (static_cast <const Phyl*> (inter->makeRoot ()) -> asFossil ());
	ASSERT (root_old_);
	root_old = root_old_;
	//
  if (isolateTransient (root_old_))
  	root_old = isolatedTransientChild; 
	ASSERT (root_old);
	ASSERT (root_old != tree. root);
	  
  // Species::time
  root_old->assignTime (); 
 	const_static_cast <Species*> (tree. root) -> assignTime ();
	
  // Species::parent2core[]
	root_old->assignFeaturesUp (nullptr);
   
  ASSERT (eqTreeLen (tree. len, tree. getLength ()));  // <= tree.allTimeZero
}



void ChangeRoot::restore_ ()
{
	ASSERT (root_old);
	ASSERT (tree. root == from->getParent ());
	

  // Species::parent2core[]
	root_old->restoreFeaturesUp (nullptr);
	
  // Species::time
 	root_old->restoreTime (); 
 	const_static_cast <Species*> (tree. root) -> restoreTime ();
	
  // Topology
  if (isolatedTransient) 
  {
  	isolatedTransient->attach (const_cast <FeatureTree&> (tree));
  	const Tree::TreeNode* parent_old = isolatedTransientChild->getParent ();
	  isolatedTransientChild->setParent (isolatedTransient);
	  isolatedTransient     ->setParent (const_cast <Tree::TreeNode*> (parent_old));
	  root_old = isolatedTransient;
  }
	EXEC_ASSERT (root_old->makeRoot () == from->getParent ());
  const Fossil* inter = static_cast <const Fossil*> (from->getParent ());
  const_cast <Species*> (from) -> setParent (const_static_cast <Fossil*> (inter->getParent ()));  
  delete inter;
}



void ChangeRoot::commit_ ()
{
	ASSERT (root_old);

  // Species::time
 	root_old->commitTime (); 
 	const_static_cast <Species*> (tree. root) -> commitTime ();

  // Species::parent2core[]
	root_old->commitFeaturesUp (nullptr);
}




// ChangeDel

void ChangeDel::apply_ ()
{
	ASSERT (! oldParent);
  ASSERT (fromChildren. empty ());


	oldParent = const_static_cast <Fossil*> (from->getParent ());
	ASSERT (oldParent);
	fromChildren = from->getChildren ();
	ASSERT (fromChildren. size () >= 2);
	
  // Topology
	const_cast <Species*> (from) -> detachChildrenUp ();
	ASSERT (static_cast <const Species*> (fromChildren. front ()) -> getParent () == oldParent);

  // Species::time
  for (const DiGraph::Node* child : fromChildren)
  	const_static_cast <Species*> (child) -> assignTime ();
	
  // Species::parent2core[]
  for (const DiGraph::Node* child : fromChildren)
  	const_static_cast <Species*> (child) -> assignFeaturesUp (oldParent);
	oldParent->assignFeaturesUp (nullptr);
}



void ChangeDel::restore_ ()
{
	ASSERT (! isolatedTransient);
	ASSERT (! isolatedTransientChild);
	ASSERT (oldParent);
  ASSERT (! fromChildren. empty ());
	

  // Species::parent2core[]
  for (const DiGraph::Node* node : fromChildren)
  	const_static_cast <Species*> (node) -> restoreFeaturesUp (oldParent);
	oldParent->restoreFeaturesUp (nullptr);
	
  // Species::time
  for (const DiGraph::Node* node : fromChildren)
  	const_static_cast <Species*> (node) -> restoreTime ();
	
  // Topology
  Species* from_ = const_cast <Species*> (from);
	from_->attach (const_cast <FeatureTree&> (tree));
  from_->setParent (oldParent);
  for (const DiGraph::Node* node : fromChildren)
    const_static_cast <Species*> (node) -> setParent (from_);
}



void ChangeDel::commit_ ()
{
	ASSERT (! isolatedTransient);
	ASSERT (! isolatedTransientChild);
	ASSERT (oldParent);
  ASSERT (! fromChildren. empty ());

  // Species::time
  for (const DiGraph::Node* node : fromChildren)
  	const_static_cast <Species*> (node) -> commitTime ();
	
  // Species::parent2core[]
  for (const DiGraph::Node* node : fromChildren)
  	const_static_cast <Species*> (node) -> commitFeaturesUp (oldParent);
	oldParent->commitFeaturesUp (nullptr);

	const_cast <FeatureTree&> (tree). delayDelete (const_cast <Species*> (from));
}




// FeatureTree

namespace
{
  // PAR
  static const uint areaDistance_std = 5;  // >= 4 <= ChangeToCousin can be applied  was: 3

  const char* timeOptimFracS = "#time_optim_frac=";
  const char* lambda0S       = "#lambda0=";
  const char* time_initS     = "#time_init=";
}



FeatureTree::FeatureTree (const string &treeFName,
      						        const string &geneDir,
      						        const string &coreFeaturesFName,
      						        bool nominalSingletonIsOptional,
  	                      bool preferGain_arg)
: inputTreeFName (treeFName)
, preferGain (preferGain_arg)
{
	ASSERT (! treeFName. empty ());

  const Chronometer_OnePass cop ("Initialization");  

 	loadPhylFile (treeFName);
  ASSERT (root);
  ASSERT (nodes. front () == root);
  ASSERT (static_cast <const Phyl*> (root) -> asSpecies ());

  if (geneDir. empty ())
    return;
   
  size_t genomes = 0; 
  {
    cerr << "Genomes ..." << endl;
	  Progress prog;
	  // Threads ??
	 	for (const DiGraph::Node* node : nodes)
	 		if (const Genome* g = static_cast <const Phyl*> (node) -> asGenome ())
	 		{
	 			prog (g->getName ());
	      var_cast (g) -> initDir (geneDir, nominalSingletonIsOptional);
	      genomes++;
	    }
	}
	ASSERT (genomes);

  // nominal2values, Genome::nominals
  ASSERT (nominal2values. empty ());
 	for (const DiGraph::Node* node : nodes)
 		if (const Genome* g = static_cast <const Phyl*> (node) -> asGenome ())
 		  try { var_cast (g) -> coreSet2nominals (); }
 		    catch (const exception &e)
 		      { throw runtime_error ("In genome " + g->id + ": " + e. what ()); }

  // Genome::coreSet
  {
    Progress prog (genomes);
    // 86 sec./50K genomes
   	for (const DiGraph::Node* node : nodes)
   		if (const Genome* g = static_cast <const Phyl*> (node) -> asGenome ())
   		{
   		  prog ();
   		  var_cast (g) -> nominals2coreSet (); 
   		}
  }

  // optionalCore[i] in all Genome's => remove Feature i ??

  // globalSingletons, nonSingletons
  Set<Feature::Id> globalSingletons;
  Set<Feature::Id> nonSingletons;
  // 196 sec./50K genomes
  {
    Progress prog (genomes);
   	for (const DiGraph::Node* node : nodes)
   		if (const Genome* g = static_cast <const Phyl*> (node) -> asGenome ())
   		{
   		  prog ();
   		  g->getSingletons (globalSingletons, nonSingletons);
   		}
  }
  ASSERT (! globalSingletons. intersects (nonSingletons));
  globalSingletonsSize = globalSingletons. size ();
  {
    Progress prog (genomes);
   	for (const DiGraph::Node* node : nodes)
   		if (Genome* g = const_cast <Genome*> (static_cast <const Phyl*> (node) -> asGenome ()))
   		{
   		  prog ();
        globalSingletonsSize += g->setSingletons (globalSingletons);
      }
  }
  
  // commonCore
  // 64 sec./50K genomes
  ASSERT (commonCore. empty ());
  commonCore = nonSingletons;
  {
    Progress prog (genomes);
   	for (const DiGraph::Node* node : nodes)
   		if (const Genome* g = static_cast <const Phyl*> (node) -> asGenome ())
   		{
   		  prog ();
   		  g->getCommonCore (commonCore);
   		}
  }
  ASSERT (! globalSingletons. intersects (commonCore));
  {
    Progress prog (genomes);
   	for (const DiGraph::Node* node : nodes)
   		if (Genome* g = const_cast <Genome*> (static_cast <const Phyl*> (node) -> asGenome ()))
      { 
        prog ();
        EXEC_ASSERT (g->removeFromCoreSet (commonCore) == commonCore. size ()); 
      }
  }      
  for (const Feature::Id& fId : commonCore)
    nonSingletons. erase (fId);    
    
  if (featuresExist () && nonSingletons. empty ())
  	throw runtime_error ("All genes are singletons or common core");
  	
  // features, feature2index
  ASSERT (features. empty ());
  map<Feature::Id, size_t/*index*/> feature2index;
#ifndef NDEBUG
  Feature::Id prevFeature;
#endif
  for (const Feature::Id& fId : nonSingletons)
  {
  	ASSERT (prevFeature < fId);
  	feature2index [fId] = features. size ();
  	features << move (Feature (fId));
  #ifndef NDEBUG
  	prevFeature = fId;
  #endif
  }
  features. searchSorted = true;
//genes = features. size ();
  ASSERT (features. size () == nonSingletons. size ());
  ASSERT (feature2index. size () == features. size ());
  cerr << "# Features: " << features. size () << endl;

  // allTimeZero, Phyl::init()   
  // If no optimization then optimize each Feature separately !??
  cerr << "Memory ..." << endl;
  {   
    // 288 sec./50K genomes
    size_t timeNan = 0;
    size_t timeNonNan = 0;
    Progress prog (nodes. size ());
   	for (const DiGraph::Node* node : nodes)
   	{
   	  prog ();
   		const Phyl* p = static_cast <const Phyl*> (node);
   		if (const Species* s = p->asSpecies ())
   		{
     		if (isNan (s->time))
     			timeNan++;
     		else
     			timeNonNan++;
   			const_cast <Species*> (s) -> init ();
   		}
   		else if (const Genome* g = p->asGenome ())
   			const_cast <Genome*> (g) -> init (feature2index);
   		else
   			ERROR;
   	}
   	ASSERT (timeNan || timeNonNan);
   	ASSERT (! (timeNan && timeNonNan));
   	allTimeZero = timeNan;
  }

  genomeGenes_ave = 0;  
 	for (const DiGraph::Node* node : nodes)
 		if (const Genome* g = static_cast <const Phyl*> (node) -> asGenome ())
 			genomeGenes_ave += g->getGenes ();
  genomeGenes_ave /= genomes;

  if (! allTimeZero && featuresExist ())
    loadSuperRootCoreFile (coreFeaturesFName);
 	for (DiGraph::Node* node : nodes)  
 		static_cast <Phyl*> (node) -> setWeight ();
  cerr << "Assigning features ..." << endl;
  setLenGlobal ();  // 399 sec./50K genomes
  setCore ();  // 67 sec./50K genomes
  len_min = getLength_min ();

  ASSERT (nodes. front () == root);
#ifndef NDEBUG
  {
    size_t index = 0;
   	for (const DiGraph::Node* node : nodes)
   	{
   		ASSERT (static_cast <const Phyl*> (node) -> index_init == index);
   		index++;
   	}
  }
#endif
}



namespace 
{
	Real str2time (const string &s)
	{
	  const string s1 (" " + s);
		const string tStr (" t=");
		const size_t pos = s1. find (tStr);
		if (pos == string::npos)
			return NaN;  
		string timeS (s1. substr (pos + tStr. size ()));
		return str2real (findSplit (timeS));
	}
}



bool FeatureTree::loadPhylLines (const StringVector& lines,
  						                   size_t &lineNum,
  						                   Species* parent,
  						                   size_t expectedOffset)
{ 
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
		ERROR;
	}
	
	lineNum++;

  string idS (findSplit (s));
	ASSERT (isRight (idS, ":"));
	idS. erase (idS. size () - 1);
	if (isLeft (idS, "g"))
	{
	  idS. erase (0, 1);
	  ASSERT (parent);
	  Strain* st = const_cast <Strain*> (parent->asStrain ());
		new Genome (*this, st, idS);
	}
	else
	{
    const Real time = str2time (s);
    Fossil* fossilParent = nullptr;
    if (parent)
    {
      fossilParent = const_cast <Fossil*> (parent->asFossil ());
      ASSERT (fossilParent);
    }
	  Species* sp = nullptr;
	  if (isLeft (idS, "s"))
	  {
  	  idS. erase (0, 1);
  		sp = new Strain (*this, fossilParent, idS, time);
	  }
	  else
	  {
  		const string id (isLeft (idS, "0x") ? string () : idS);
  		sp = new Fossil (*this, fossilParent, id, time);
    }
    ASSERT (sp);
		while (loadPhylLines (lines, lineNum, sp, expectedOffset + Offset::delta))
			;
	}
	
	return true;
}



namespace
{

bool getParam (const string &line,
               const string &paramS,
               Real &param)
// Update: param
{
  if (! isLeft (line, paramS))
    return false;
  param = str2real (line. substr (paramS. size ()));
  return true;
}

}



void FeatureTree::loadPhylFile (/*int root_species_id,*/
	                              const string &treeFName)
{
	ASSERT (! treeFName. empty ());
	
		
  StringVector lines;
  {
    LineInput in (treeFName, 10000);  // PAR
    lines = in. getVector ();
  }
  ASSERT (! lines. empty ());
	size_t lineNum = 0; 

  if (getParam (lines [lineNum], timeOptimFracS, timeOptimFrac))
  	lineNum++;
  if (getParam (lines [lineNum], lambda0S, lambda0))
  	lineNum++;
  if (getParam (lines [lineNum], time_initS, time_init))
  	lineNum++;

  Species* sp = nullptr;  	
	EXEC_ASSERT (loadPhylLines (lines, lineNum, sp, 0));
#if 0
  // root_species_id
 	for (const DiGraph::Node* node : nodes)
 		if (const Species* s = static_cast <const Phyl*> (node) -> asSpecies ())
 			if (s->id == (uint) root_species_id)
 				const_cast <Species*> (s) -> id = 0;
	const_static_cast <Species*> (root) -> id = (uint) root_species_id;
#endif
}



void FeatureTree::qc () const
{ 
  if (! qc_on)
    return;
	Tree::qc ();
		
	ASSERT (! inputTreeFName. empty ());

	ASSERT (root);
	ASSERT (root->graph == this);
	const Species* root_ = static_cast <const Phyl*> (root) -> asSpecies ();
	ASSERT (root_);
//ASSERT (nodes. size () <= numeric_limits<unsigned short>::max ());  // To fit Phyl::CoreEval::treeLen
	if (! allTimeZero && featuresExist ())
  {
    ASSERT (root_->time == getRootTime ());
    if (! emptySuperRoot)
    	for (const bool i : {false, true})
    	{
    	  ASSERT (root_->weight [i] [  i] == 0);
    	  ASSERT (root_->weight [i] [! i] == INF);
    	}
  }

  Set<size_t> indices;
  size_t n = 0;
 	for (const DiGraph::Node* node : nodes)
 	{
 		const Phyl* p = static_cast <const Phyl*> (node);
 		ASSERT ((bool) p->asStrain () == (bool) p->isTransient ());
 		ASSERT (! indices. contains (p->index_init));
 		indices << p->index_init;
 		if (const Genome* g = p->asGenome ())
 		  n += g->singletons. size ();
 	}
 	ASSERT (n == globalSingletonsSize);

  if (allTimeZero || ! featuresExist ())
  { 
    ASSERT (timeOptimWhole ());
    ASSERT (isNan (lambda0)); 
    ASSERT (isNan (time_init)); 
  }
  else
  {
    ASSERT (superRootCore. size () == features. size ());
    ASSERT (isProb (timeOptimFrac));
	  ASSERT (isProb (lambda0));
	  ASSERT (lambda0 < 1 /*0.5*/);
    ASSERT (time_init > 0);
	}

	Feature prevFeature;
	FFOR (size_t, i, features. size ())
	{
	  const Feature f (features [i]);
	  f. qc ();
	  ASSERT (! Feature::singleton (f. name));
		ASSERT (/*prevFeature. isGene == f. isGene,*/ prevFeature. name < f. name);
	//ASSERT (prevFeature. isGene >= f. isGene);
	//ASSERT (f. isGene == (i < genes));
		prevFeature = f;
	}
	ASSERT (features. searchSorted);

//ASSERT (genes <= features. size ());
	
	ASSERT (geReal (len, len_min));
	ASSERT (len_min >= 0);
	
	IMPLY (! features. empty (), genomeGenes_ave > 0);
	ASSERT (genomeGenes_ave <= getTotalGenes ());
	
	ASSERT ((bool) distDistr. get () == (bool) depthDistr. get ());
	
	if (coreSynced) 
	{
		if (! emptySuperRoot)
			FFOR (size_t, i, root_->core. size ())
			{
			  ASSERT (root_->core [i] == root_->feature2parentCore (i));
			  IMPLY (! allTimeZero && featuresExist (), root_->core [i] == superRootCore [i]);
			}
		const float subtreeLen = static_cast <const Phyl*> (root) -> getSubtreeLength ();
		if (! eqReal (len, subtreeLen, len_delta)) 
		{
			cout << fixed; cout. precision (10);  // PAR
			cout        << subtreeLen 
			         << " " << len 
			       //<< " " << root_->pooledSubtreeDistance 
			         << " " << root_->getDistance ()
			         << endl;
			ERROR;
		}
	}
}



void FeatureTree::saveText (ostream& os) const
{ 
  ASSERT (coreSynced);
  
  if (! allTimeZero && featuresExist ())
  {
  	ONumber oN (os, 6, true);  // PAR
    if (! timeOptimWhole ())
  	  os << timeOptimFracS << timeOptimFrac << endl;
  	os << lambda0S       << lambda0       << endl;
  	os << time_initS     << time_init     << endl;
  }
  
	Tree::saveText (os); 
}



void FeatureTree::print (ostream& os) const
{ 
	os << "TREE" << endl;
	saveText (os);
#if 0
  if (! taxNamePrefix. empty ())
  { 
  	os << endl;
    os << "Legend:" << endl;
    os << abbreviationLegend () << endl;
  }
#endif
}



void FeatureTree::deleteLeaf (TreeNode* node,
                              bool deleteTransientAncestor)  
{
  ASSERT (features. empty ());
  
  ASSERT (node);
  ASSERT (& node->getTree () == this);
  ASSERT (node != root);
  
  const Genome* leaf = static_cast <Phyl*> (node) -> asGenome ();
  ASSERT (leaf);
  if (verbose ())
    cout << "Deleting: " << leaf->getName () << endl;

  const Strain* parent = static_cast <const Phyl*> (leaf->getParent ()) -> asStrain ();
  ASSERT (parent);
  ASSERT (parent->isTransient ());

  const Fossil* grandparent = static_cast <const Phyl*> (parent->getParent ()) -> asFossil ();
  ASSERT (grandparent);
  
  const_cast <Genome*> (leaf) -> detachChildrenUp ();
  
	const_cast <Strain*> (parent) -> detachChildrenUp ();
	delayDelete (const_cast <Strain*> (parent));

  if (deleteTransientAncestor && grandparent->isTransient ())
  {
  	const_cast <Fossil*> (grandparent) -> detachChildrenUp ();
  	delayDelete (const_cast <Fossil*> (grandparent));
  }
	
  delete leaf;
  toDelete. deleteData ();
}



void FeatureTree::printInput (ostream& os) const
{
  os << "Tree from file: " << inputTreeFName << endl;
	  
  const size_t genomes = root->getLeavesSize ();
  os << "# Genomes: " << genomes << endl;
  os << "# Species: " << nodes. size () - genomes << endl;
  os << "# Common core genes: " << commonCore. size () << endl;
  os << "# Singleton genes:   " << globalSingletonsSize << endl;
  os << "# Other genes:       " << features. size () << endl;
  os << "Genome size ave.:    " << genomeGenes_ave << endl;

  os << "Time: " << (allTimeZero ? "Not used" : "Used") << endl;
  if (allTimeZero)
  {
    if (preferGain)
      os << "Gains are preferred over losses" << endl;
  }
  else
  {
  	ONumber oNum (os, 3, true);  // PAR
    os << "Lambda_0          = " << lambda0 << endl;
    os << "Initial time      = " << time_init << endl;
    os << "timeOptimFrac     = " << timeOptimFrac << endl;
    if (! superRootCore. empty ())
    {
      size_t n = commonCore. size ();
      FFOR (size_t, i, features. size ())
        if (superRootCore [i])
          n++;
      os << "# Root core genes = " << n << endl;
    }
  }
  
 	ONumber oNum (os, 0, false);  // PAR
  os << "Tree length min.  = " << len_min << endl;
  os << "Tree length       = " << len << endl;
  os << endl;
}



void FeatureTree::dump (const string &fName)
{ 
  setCore ();
  setStats ();
#if 0
  if (setIds)
    if (const Fossil* f = static_cast <const Species*> (root) -> asFossil ())
    { 
      uint id = 1;  // 0 means "unknown"
      const_cast <Fossil*> (f) -> setId (id);
    }
#endif
  qc ();      
  saveFile (fName);
}



float FeatureTree::feature2treeLength (size_t featureIndex) const
{ 
  // Cf. Phyl:;feature2parentCore()
  const auto& parent2core = static_cast <const Species*> (root) -> parent2core;
  return emptySuperRoot 
           ? parent2core [false] [featureIndex]. treeLen
           : allTimeZero 
           	     ? min ( parent2core [false] [featureIndex]. treeLen 
      		             , parent2core [true]  [featureIndex]. treeLen
      		             )
  		           : parent2core [superRootCore [featureIndex]] [featureIndex]. treeLen;
}



float FeatureTree::getLength () const
{
  float s = static_cast <const Species*> (root) -> pooledSubtreeDistance;
  FFOR (size_t, i, features. size ())
    s += feature2treeLength (i);
  return s;
}



void FeatureTree::setTimeWeight ()
{
  setCore ();
 	for (DiGraph::Node* node : nodes)
 	  if (const Species* s = static_cast <Phyl*> (node) -> asSpecies ())
 		  const_cast <Species*> (s) -> setTimeWeight (); 
}



void FeatureTree::optimizeTime ()
{
  ASSERT (len < INF);
  
  Real len_old = INF;  
	while (len_old > len + len_delta) 
	{
		len_old = len;
	//cout << "Before: len = " << len << endl;  
		setTimeWeight ();
		setLenGlobal ();
	//cout << "After:  len = " << len << endl;  
		if (verbose ())
		{
	  	cout << "len (time optimized) = " << len << endl;
			qc ();
		}
		if (timeOptimWhole () && ! leReal (len, len_old + len_delta)) 
		{
			cout << len << " " << len_old << endl;
			ERROR;
		}
	}
}



#if 0
namespace 
{

struct LambdaFunc : Func1
{
	FeatureTree* tree;
	  // !nullptr

  Real f (Real x) 
  { 
  	ASSERT (tree);
  	  	
	 	for (DiGraph::Node* node : tree->nodes)  
	 		static_cast <Phyl*> (node) -> rememberAllFeatures ();

    tree->lambda0 = x;
  	Unverbose unv;
  	tree->optimizeTime ();
  	
   	for (DiGraph::Node* node : tree->nodes)
 		  static_cast <Phyl*> (node) -> restoreFeatures ();

  	return tree->len;
  }  
};

}



void FeatureTree::optimizeLambdaTime ()
{
	const Real len_old = len;

  if (verbose ())
    qc ();
	Real x = NaN;
	if (allTimeZero)
		x = 0.5;  // PAR
	else
	{
		LambdaFunc f;
		f. tree = this;
		Unverbose unv;
		EXEC_ASSERT (f. minimizeConvex (0, 0.5, lambda0, 1e-6, 100, x));  // PAR
	}
  if (verbose ())
    qc ();
  lambda0 = x;
  optimizeTime ();  
  if (verbose ())
    qc ();
  
 	ASSERT (leReal (len, len_old + len_delta));
  	
	if (verbose () && ! allTimeZero)
	{
	  ONumber oNum (cout, 3, true);  // PAR
	  cout << "lambda0 = " << lambda0
         << "  len = " << len
         << endl;
  }
}
#endif



void FeatureTree::useTime (const string &coreFeaturesFName)
{ 
	ASSERT (allTimeZero);
	ASSERT (featuresExist ());

	allTimeZero = false;
  loadSuperRootCoreFile (coreFeaturesFName);
  setCore ();
  timeOptimFrac = 0;


  size_t parent2core [2/*thisCore*/] [2/*parentCore*/];
  getParent2core_sum (parent2core);
	if (verbose ())
  	for (const bool i : {false, true})
    	for (const bool j : {false, true})
        cout << "parent2core[" << (int) i << "][" << (int) j << "]=" << parent2core [i] [j] << endl;

  float len_old;  
  float len1 = INF;;
	lambda0 = (Real) parent2core [true] [false] / ((Real) parent2core [true] [false] + (Real) parent2core [false] [true]);  // Optimal if time_init = 0
	do 
	{
		len_old = len1;

    time_init = coreChange2time (*this, parent2core);
    if (verbose ())
      cout << "time_init = " << time_init << endl;
    ASSERT (time_init > 0);
  
    lambda0 = getLambda0_commonTime (parent2core, time_init);  
    if (verbose ())
      cout << "lambda0 = " << lambda0 << endl;
  
   	for (DiGraph::Node* node : nodes)  
   	{
   	  Phyl* p = static_cast <Phyl*> (node);
   		if (const Species* s = p->asSpecies ())
  		  const_cast <Species*> (s) -> time = s->getParent () ? time_init : getRootTime ();
 		  p->setWeight ();
   	}
   	
   	len1 = static_cast <const Phyl*> (root) -> getSubtreeLength ();
	  if (verbose ())  
	  	cout << "len1 = " << len1 << endl;
		if (! leReal (len1, len_old + len_delta, 1e-3))  // PAR
		{
			cout << len1 << " " << len_old << endl;
			ERROR;
		}
	}
	while (len_old > len1 + len_delta);


  {
    Unverbose unv;
    if (verbose ())
    {
      const Prob lambda0_old = lambda0;
  	  optimizeLambda0 ();  
  	  ASSERT (eqReal (lambda0, lambda0_old, 1e-4));
    }
  }


  // Cf. FeatureTree::FeatureTree()
	setLenGlobal ();  
	ASSERT (leReal (len, len1));
	setCore ();  
  len_min = getLength_min ();

  if (verbose ())
    cout << "Use time: on" << endl;
  qc ();
}



void FeatureTree::getParent2core_sum (size_t parent2core [2/*thisCore*/] [2/*parentCore*/]) const
{
  ASSERT (coreSynced);

	for (const bool i : {false, true})
  	for (const bool j : {false, true})
      parent2core [i] [j] = 0;
 	for (const DiGraph::Node* node : nodes)  
 		if (const Species* s = static_cast <const Phyl*> (node) -> asSpecies ())
 		{
   		if (! s->getParent ())
   			continue;
  
      size_t parent2core_ [2/*thisCore*/] [2/*parentCore*/];
      s->getParent2corePooled (parent2core_);	  
  		FFOR (size_t, i, features. size ())
  		  parent2core_ [s->core [i]] [s->feature2parentCore (i)] ++;
  
    	for (const bool i : {false, true})
      	for (const bool j : {false, true})
          parent2core [i] [j] += parent2core_ [i] [j];
  	}
}



namespace 
{

struct LambdaCommonTimeFunc : Func1
{
  Real parent2core [2/*thisCore*/] [2/*parentCore*/];
  Real a;
  

  LambdaCommonTimeFunc (const size_t parent2core_arg [2/*thisCore*/] [2/*parentCore*/],
                        Real time)
		: a (1 - exp (- time))
		{
      ASSERT (time > 0);
      ASSERT (time < INF);
    	for (const bool i : {false, true})
      	for (const bool j : {false, true})
   	      parent2core [i] [j] = (Real) parent2core_arg [i] [j];
    }

  
  Real f (Real x) 
  { 
    ASSERT (x > 0);
    ASSERT (x < 1);
    const Real x1 = 1 - x;
  	return   parent2core [false] [true]  / x1
  	       + parent2core [false] [false] * a / (1 - x * a)
  	       - parent2core [true]  [false] / x
  	       - parent2core [true]  [true]  * a / (1 - x1 * a);
  }  
};

}



Real FeatureTree::getLambda0_commonTime (const size_t parent2core [2/*thisCore*/] [2/*parentCore*/],
                                         Real commonTime) 
{
	LambdaCommonTimeFunc f (parent2core, commonTime);
	const Real delta = 1e-6;  // PAR
	return f. findZero (delta, 1 - delta, delta);
}



namespace 
{
  
struct LambdaFunc : Func1
{
  FeatureTree& tree;
  
  LambdaFunc (FeatureTree &tree_arg)
  : tree (tree_arg)
  {
    ASSERT (! tree. allTimeZero);
    ASSERT (tree. coreSynced);
  }
  
  
  Real f (Real x) 
  { 
    ASSERT (x > 0);
    ASSERT (x < 1);

    tree. lambda0 = x;
   	for (DiGraph::Node* node : tree. nodes)  
   		static_cast <Phyl*> (node) -> setWeight ();
    return static_cast <const Species*> (tree. root) -> getSubtreeLength ();
  }  
};

}



void FeatureTree::optimizeLambda0 () 
{
  setCore ();
  
	LambdaFunc f (*this);
	const Real delta = 1e-6;  // PAR
	f. minimizeConvex (delta, 1 - delta, lambda0, delta, 100, lambda0);  // PAR
	f. f (lambda0);  // Output: Phyl::weight[][]
}



const Change* FeatureTree::getBestChange (const Species* from) 
{
	ASSERT (from);
	ASSERT (from != root);
	
 	const Change* bestChange = nullptr;

 	#define TRY_CHANGE(T,P)  if (T::valid_ P) tryChange (new T P, bestChange)

  VectorPtr<TreeNode> area, boundary;
  from->getArea (areaDistance_std, area, boundary);  
  if (verbose (1))
    cerr << " area=" << area. size () << " ";

 	for (const Tree::TreeNode* node : area)  
 		if (Species* to = const_cast <Species*> (static_cast <const Phyl*> (node) -> asSpecies ()))
 		{
   	  TRY_CHANGE (ChangeToSibling, (from, to));  
   	  TRY_CHANGE (ChangeToUncle,   (from, to));   
    //TRY_CHANGE (ChangeToCousin,  (from, to));  
   	//TRY_CHANGE (ChangeToParent,  (from, to));  
  	  if (verbose ())
  	  {
    	  ASSERT (to->graph);
  	    to->qc ();
  	  }
    }
  
//TRY_CHANGE (ChangeRoot, (from));
  TRY_CHANGE (ChangeDel,  (from));  

  #undef TRY_CHANGE
  
  if (bestChange)
  {
  	ASSERT (positive (bestChange->improvement));
	  if (bestChange->improvement / len >= (allTimeZero ? 1e-4 : 1e-6))  // PAR  --> 0 ??
	  {
      if (verbose (1))
        cerr << "found ";
      return bestChange;	
    }
  	delete bestChange;
	}

  return nullptr;
}



bool FeatureTree::applyChanges (VectorOwn<Change> &changes)
{ 
	ASSERT (toDelete. empty ());
	
	
	clearStats ();

  // Phyl::stable: init
  for (DiGraph::Node* node : nodes)
    static_cast <Phyl*> (node) -> stable = true;
	
	
  const float len_init = len;


  changes. sort (Change::compare);	
  Set<const TreeNode*> changedNodes;  	
	for (const Change* ch_ : changes)
	{
		Change* ch = const_cast <Change*> (ch_);
		ASSERT (ch);
  	  
	  if (! ch->valid ())
  		continue;

    // Redundant filtering ??
	  bool bad = false;
    for (const TreeNode* target : ch->targets)
	    if (   ! target->graph 
	        || changedNodes. contains (target)
	       )
	    {
	    	bad = true;
	    	break;
	    }
	  if (bad)
	  	continue;

    if (verbose ())
    {
  	  cout << "Apply: ";
  	  ch->print (cout);
  	}
		ch->qc ();  
	#ifndef NDEBUG
	  const bool first = ch == changes. front ();  
	#endif
	  changedNodes << ch->getCoreChanged ();
  	float len_old = len;
    ch->apply ();
	  len = getLength ();
	  if (verbose ())
	  {
  	  ONumber on (cout, 0, false);
  	  cout << "len = " << len << endl;
  	}
	  IMPLY (first, eqTreeLen (len, len_old - ch->improvement));  // Should be: old improvement ??
	  if (verbose ())
	  {
	  	qc ();
	  	const Real len_old1 = len;
	  	setLenGlobal ();
	  	if (! eqReal (len, len_old1, len_delta)) 
	  	{
	  		cout << len << " " << len_old1 << endl;
	  		ch->restore ();
	  	//setCore ();
	  	//print (cout);
	  		ERROR;
	  	}
	  }
	  
	  if (geReal (len, len_old))
	  {
	  	ch->restore ();
	  	ASSERT (! first);
	    len = getLength ();
	    if (verbose ())
	    {
  	    cout << "restore" << endl;
  		  cout << "len = " << len << endl;
  		}
	  	break;  // ??
	  }
	  else
	  {
	  	ch->commit ();

      // Phyl::stable
      for (const TreeNode* target : ch->targets)
  	    if (target->graph)  
  	    {
          VectorPtr<TreeNode> area, boundary;
          target->getArea (areaDistance_std, area, boundary);  
          for (const TreeNode* areaNode : area)
            const_static_cast <Phyl*> (areaNode) -> stable = false;
        }
	  }

  	len_old = len;
    optimizeTime (); 
    if (verbose ())
	    cout << "len = " << len << endl;
    if (! leReal (len, len_old + len_delta)) 
    {
    	if (timeOptimWhole ())
    	{
	    	cout << len << " " << len_old << endl;
	    	ERROR;
	    }
	    break;
    }
	}
	

  finishChanges ();
  

  const float improvement = max<float> (0.0, len_init - len);
  cout << "Improvement = " << improvement << endl;
  cout << "len = " << len << endl;
  IMPLY (timeOptimWhole (), geReal (improvement, - len_delta));

  const Prob improvementRel = improvement / (len_init - len_min);
  ASSERT (isProb (improvementRel));
  if (improvementRel <= 1e-4)  // PAR
  {
  	if (timeOptimWhole ())
  	{
      if (leReal (improvement, len_delta))
      	return false;
    }
  	else
  	{
  	  ASSERT (! allTimeZero);
    	optimizeLambda0 ();  
  		timeOptimFrac += 0.1;  // PAR 
  		minimize (timeOptimFrac, 1.0);
	  	cout << endl << "timeOptimFrac = " << timeOptimFrac << endl;
	  	cout. flush ();
    	optimizeTime ();
    	finishChanges ();
     	for (DiGraph::Node* node : nodes)
     	  static_cast <Phyl*> (node) -> stable = false;
  	}
  }

  qc ();
  
  return true;
}



bool FeatureTree::optimize () 
{ 
	VectorOwn<Change> changes;  changes. reserve (256);  // PAR
	{
	 	Vector<DiGraph::Node*> nodeVec;  nodeVec. reserve (nodes. size ());
    for (DiGraph::Node* node : nodes)
    {
      Phyl* phyl = static_cast <Phyl*> (node);
      if (! phyl->stable)
        nodeVec << phyl;
    }
		Progress prog (nodeVec. size ());
	 	for (const DiGraph::Node* node : nodeVec)  
	 	{
	   	prog ();
	 		if (const Species* from = static_cast <const Phyl*> (node) -> asSpecies ())
  	 	 	if (from != root)
    		 	if (const Change* bestChange = getBestChange (from))
    		  { 
    		  	ASSERT (positive (bestChange->improvement));
    		  	changes << bestChange;
    		  }
  	}
  }
    
  return applyChanges (changes);
}



string FeatureTree::findRoot (size_t &bestCoreSize) 
{
  ASSERT (allTimeZero);

  string rootLcaName;

  setCore ();

  const Species* bestFrom = nullptr;
 	bestCoreSize = 0;
  FFOR (size_t, i, features. size ())
    if (static_cast <const Species*> (root) -> core [i])
      bestCoreSize++;
//cout << "Init. bestCoreSize = " << bestCoreSize << endl;  
  setLeaves ();
  const size_t halfLeaves = root->leaves / 2;
 	for (const DiGraph::Node* node : nodes)
 		if (const Species* from = static_cast <const Phyl*> (node) -> asSpecies ())
 		{
 		#if 0
   		const Species* parent = static_cast <const Species*> (from->getParent ());
   		if (! parent)
   		  continue;
   	 	ASSERT (from != root);
   	#endif
      size_t coreSize = 0;
      FFOR (size_t, i, features. size ())
        if (   from  ->core [i] 
          //&& parent->core [i]
           )
        coreSize++;
      if (   minimize (bestCoreSize, coreSize)
          || (   bestCoreSize == coreSize 
              && bestFrom 
              &&   difference (    from->leaves, halfLeaves) 
                 < difference (bestFrom->leaves, halfLeaves)
             )
         )
        bestFrom = from;
   	}

  if (bestFrom)
  {
  	rootLcaName = bestFrom->getLcaName ();
    {
      ChangeRoot ch (bestFrom);    
      ch. apply ();
    	ch. commit ();
    }
    ASSERT (eqTreeLen (len, getLength ()));
    clearStats ();
    finishChanges ();
  }    
  
  // superRootCore
  bool changed = bestFrom;
  ASSERT (superRootCore. empty ());
  superRootCore. resize (features. size (), false);
  FFOR (size_t, i, features. size ())
    if (getSuperRootCore (i))
  	{
      superRootCore [i] = true;
      changed = true;
    }
      
  if (changed)
    setCore ();
    
  sort ();

  return rootLcaName;
}



void FeatureTree::resetSuperRootCore (size_t coreChange [2/*core2nonCore*/])
{
  ASSERT (! allTimeZero);
  ASSERT (coreSynced);
  
	for (const bool b : {false, true})
 	  coreChange [b] = 0;

	const Species* root_ = static_cast <const Species*> (root);
	if (root_->arcs [false]. size () < 2)
	  return;
	  
  FFOR (size_t, i, superRootCore. size ())
  {
    size_t n = 0;
  	for (const DiGraph::Arc* arc : root_->arcs [false])
    {
      const Phyl* s = static_cast <Phyl*> (arc->node [false]);
  	  if (s->core [i])
  	    n++;
  	}
  	if (n <= 1)
  	  if (superRootCore [i])
  	  {
  	    superRootCore [i] = false;
  	  //if (i < genes)
  	      coreChange [true] ++;
  	  }
  	if (n == root_->arcs [false]. size ())
  	  if (! superRootCore [i])
  	  {
  	    superRootCore [i] = true;
  	  //if (i < genes)
  	      coreChange [false] ++;
  	  }
  }

  len = getLength ();
  setCore ();
  clearStats ();
}



void FeatureTree::loadSuperRootCoreFile (const string &coreFeaturesFName)
{ 
	ASSERT (! allTimeZero);
	ASSERT (superRootCore. empty ());
	  
  superRootCore. resize (features. size (), false);

	if (coreFeaturesFName. empty ())
	  return;

  typedef  map<Feature::Id, size_t>  Feature2index;
  Feature2index feature2index;
  FFOR (size_t, i, features. size ())
    feature2index [features [i]. name] = i;
  ASSERT (feature2index. size () == features. size ());

  {
  	LineInput f (coreFeaturesFName, 10 * 1024);  // PAR
  	size_t miss = 0;
  	while (f. nextLine ())
  	{
  	  trim (f. line);
  	  Feature2index::const_iterator it = feature2index. find (f. line);
  	  if (it == feature2index. end ())
  	    miss++;
  	  else
  	    superRootCore [it->second] = true;
  	}
  	ASSERT (miss == commonCore. size ());
  }
}



void FeatureTree::saveSuperRootCore (const string &coreFeaturesFName) const
{
  ASSERT (! coreFeaturesFName. empty ());
//ASSERT (! superRootCore. empty ());

 	Set<Feature::Id> s (commonCore);
  FFOR (size_t, i, features. size ())
    if (superRootCore [i])
      s << features [i]. name;
      
  OFStream f (coreFeaturesFName);
  for (const Feature::Id& fId : s)
    f << fId << endl;
}



void FeatureTree::setStats () 
{
  ASSERT (coreSynced);

  clearStats ();  
  
  Dataset ds;
  auto dist  = new RealAttr1 ("Distance", ds);  
  auto depth = new RealAttr1 ("Depth",    ds);  
 	for (const DiGraph::Node* node : nodes)
 	{
 		const Phyl* phyl = static_cast <const Phyl*> (node);
  	const size_t n = ds. appendObj (phyl->getName ());
    (*dist)  [n] = phyl->getNeighborDistance ();
    (*depth) [n] = (Real) phyl->getTopologicalDepth ();
    const Genome* g = phyl->asGenome ();
    FFOR (size_t, i, features. size ())
    {
      Feature& f = features [i];
      if (g && g->core [i])
      {
        if (g->optionalCore [i])
          f. optionalGenomes++;
        else
          f. genomes++;
      }
    #if 0
      if (phyl->core [i])
      {
        f. coreNodes++;
        if (f. name == "13-00:Aquificales")   // ??
          cout << phyl->getLcaName () << endl;
      }
    #endif
      if (phyl == root)
      {
        if (phyl->core [i])
          f. rootGain = true;
      }
      else
      {
        if (! phyl->feature2parentCore (i) && phyl->core [i])  
          f. gains << phyl;
        if (phyl->feature2parentCore (i) && ! phyl->core [i])
          f. losses << phyl;
      }
    }
  }
  
  for (Feature& f : features)
  {
  	f. gains. sort ();
  	ASSERT (f. gains. isUniq ());
   	f. losses. sort ();
  	ASSERT (f. losses. isUniq ());
  }
    
  const Sample sample (ds);

  {
    const UniVariate<NumAttr1> an (sample, *dist);
    auto distDistr_ = new Normal;
    distDistr_->analysis = & an;
    distDistr_->estimate ();
    distDistr. reset (distDistr_);
  }

  {
    const UniVariate<NumAttr1> an (sample, *depth);
    auto depthDistr_ = new Normal;
    depthDistr_->analysis = & an;
    depthDistr_->estimate ();
    depthDistr. reset (depthDistr_);
  }
}
  


void FeatureTree::clearStats () 
{
  for (Feature& f : features)
  {
    f. genomes = 0;
    f. optionalGenomes = 0;
    f. rootGain = false;
    f. gains. clear ();
    f. losses. clear ();
  }
  
  distDistr. reset (nullptr);
  depthDistr. reset (nullptr);
}
  


void FeatureTree::delayDelete (Species* s)
{
	if (! s)
		return;
	ASSERT (! s->graph);
	toDelete << s;
}



void FeatureTree::tryChange (Change* ch,
	                           const Change* &bestChange)
{ 
	ASSERT (ch);
	ASSERT (ch->from->graph == this);

  if (verbose ())
  {
    ch->print (cout); 
    ch->qc ();
  }

  ch->apply ();
  ch->restore ();

  Unverbose unv;
  if (verbose ())
  {
    ch->print (cout); 
    qc ();  
    ASSERT (eqTreeLen (len, getLength ()));
    const Real len_old = len;
    setLenGlobal ();
    if (! eqReal (len, len_old, 0.1))  // PAR
    {
      cout << "len = " << len << "  len_old = " << len_old << endl;
      ERROR;
    }
  }
  
	if (Change::compare (ch, bestChange))
	{
		delete bestChange;
  	bestChange = ch;
  }
  else
  	delete ch;
}



void FeatureTree::finishChanges ()
{
  if (const size_t n = deleteTimeZero ())
    if (verbose ())
      cout << "# Nodes deleted = " << n << endl;
  toDelete. deleteData ();
  qc ();
}



size_t FeatureTree::deleteTimeZero ()
{
  const Real len_old = len;
  
  setCore ();

  size_t n = 0;
 	Vector<DiGraph::Node*> nodeVec;  nodeVec. reserve (nodes. size ());
 	insertAll (nodeVec, nodes);
 	for (const DiGraph::Node* node : nodeVec)  
 		if (const Fossil* f = static_cast <const Phyl*> (node) -> asFossil ())
 			if (   f->getParent () 
 				  && (allTimeZero || leReal (f->time, 1e-7))  // PAR
 				  && ! f->getCoreChange ()
 				 )
 			{
 			  if (verbose ())
 			    cout << "To delete:" << f->getName () << endl;
 				const_cast <Fossil*> (f) -> detachChildrenUp ();
 				delayDelete (const_cast <Fossil*> (f));
 				n++;
 			}
 			
  setLenGlobal ();
  
 	if (qc_on)
 	  if (! eqReal (len, len_old, 1e-3))
 	  {
 	  	cout << len << " " << len_old << endl;
 	  	ERROR;
 	  }
 				
  return n;  
}



#if 0
void FeatureTree::abbreviate ()
{
	ASSERT (taxNamePrefix. empty ());
	
	
  List<List<string> > names;
 	for (const DiGraph::Node* node : nodes)
 		if (const Genome* g = static_cast <const Phyl*> (node) -> asGenome ())
 		{
 			const List<string> words (str2list (g->taxName));
 		//ASSERT (! words. empty ());
 		  names. push_back (words);
 		}
 	
 		
  while (! names. empty ())
  {
    const size_t nNames = names. size ();

    typedef  map<string, size_t>  Freq;
	  Freq freq;
	 	for (List<string>& lst : names)
	 	  freq [lst. front ()]++;
	 	  
	 	string bestWord;
	 	size_t freq_max = 0;
	  for (const auto& it : freq)
	 	  if (maximize (freq_max, it. second))
	 	  	bestWord = it. first;
	 	  	
	 	if (freq_max < nNames / 2)
	 		break;
	 		
	 	taxNamePrefix. push_back (bestWord);

  #if 1
    for (Iter <List<List<string> > > iter (names); iter. next (); )
	 	  if (iter->front () == bestWord)
	 	  {
	 	  	iter->pop_front ();
	 	  	if (iter->empty ())
 	 	  	  iter. erase ();
	 	  }
	 	  else
	 	  	iter. erase ();
  #else	  
	 	for (List<List<string> >::iterator it = names. begin (); it != names. end ();)
	 	{
	 		List<List<string> >::iterator itCur = it;
	 		it++;
	 	  if (itCur->front () == bestWord)
	 	  {
	 	  	itCur->pop_front ();
	 	  	if (itCur->empty ())
 	 	  	  names. erase (itCur);
	 	  }
	 	  else
	 	  	names. erase (itCur);
	 	}
	#endif
	}


 	for (const DiGraph::Node* node : nodes)
 		if (const Genome* g = static_cast <const Phyl*> (node) -> asGenome ())
 		{
 			List<string> words (str2list (g->taxName));
      List<string>::iterator wordsIt = words. begin ();
    	for (const string& prefix : taxNamePrefix)
    	{ 
	    	if (   wordsIt == words. end ()
	    		  || prefix != *wordsIt
	    		 )
	      	break;
	      *wordsIt = string (1, prefix. at (0)) + ".";
	      wordsIt++;
 		  }
 		  const_cast <Genome*> (g) -> taxName = list2str (words);
 		}
}



string FeatureTree::abbreviationLegend () const
{
  string abbreviation;
  for (const string& pref : taxNamePrefix)
  {
  	if (! abbreviation. empty ())
  		abbreviation += " ";
  	abbreviation += string (1, pref. at (0)) + ".";
  }
  return abbreviation + " - " + list2str (taxNamePrefix);
}
#endif



const Genome* FeatureTree::findGenome (const string &genomeId) const
{
  ASSERT (! genomeId. empty ());

 	for (const DiGraph::Node* node : nodes)
 		if (const Genome* g = static_cast <const Phyl*> (node) -> asGenome ())
 			if (g->id == genomeId)
 				return g;
  return nullptr;
}



size_t FeatureTree::findFeature (const  Feature::Id &featureName) const
{
  FFOR (size_t, i, features. size ())
    if (features [i]. name == featureName)
      return i;
  return NO_INDEX;
}



}



/*
TODO: ??

MLE should sum over internal states

Borrow Change's from distTree
Local Change's involving more than 1 node
  
*/
