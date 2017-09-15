// featureTree.cpp

#undef NDEBUG
#include "../common.inc"

#include "featureTree.hpp"

#include "../dm/optim.hpp"



namespace FeatureTree_sp
{



// Feature

Feature::Feature (const Feature::Id &name_arg,
	                bool isGene_arg)
: Named (name_arg)
, isGene (isGene_arg)
, genomes (0)
, gains (0)
, losses (0)
{
  qc ();
}



void Feature::qc () const
{
  if (! qc_on)
    return;

	ASSERT (! name. empty ());
	ASSERT (! contains (name, ";"));
//IMPLY (isGene, geneId ());
	ASSERT (gains <= genomes);
	IMPLY (genomes, gains);
	IMPLY (genomes <= 1, ! losses);
}




// Phyl	

Phyl::Phyl (FeatureTree &tree,
	          Species* parent_arg)
: TreeNode (tree, parent_arg)  // FeatureTree must be declared
, index_init (tree. nodeIndex_max++)
//, hasPhenChange (false)
, stable (false)
{
	FOR (unsigned char, i, 2)
		FOR (unsigned char, j, 2)
		  weight [i] [j] = NAN;
}



void Phyl::init () 
{ 
	const size_t n = getFeatureTree (). features. size ();
	FOR (unsigned char, parentCore, 2)
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
	FOR (unsigned char, parentCore, 2)
  	ASSERT (parent2core [parentCore]. size () == n);
	ASSERT (core. size () == n);
	FOR (unsigned char, i, 2)
	  FOR (unsigned char, j, 2)
	    ASSERT (weight [i] [j] >= 0);
#if 0
	if (const Phyl* p = static_cast <const Phyl*> (getParent ()))
	  IMPLY (hasPhenChange, p->hasPhenChange);
#endif

	FOR (size_t, i, core. size ())
	{
 	  ASSERT (parent2core [false] [i]. core <= parent2core [true] [i]. core); 
 	//IMPLY (getFeatureTree (). allTimeZero, abs (parent2core [false] [i]. treeLen - parent2core [true] [i]. treeLen) <= 1.001); ??
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
		FOR (unsigned char, i, 2)
		  FOR (unsigned char, j, 2)
		    weight_ [i] [j] = weight [i] [j];
		const_cast <Phyl*> (this) -> setWeight ();
		FOR (unsigned char, i, 2)
		  FOR (unsigned char, j, 2)
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

#if 0	   
	const string phenChange (getPhenChange (true));
	if (! phenChange. empty ())
	  os << "  (" << phenChange + ")";
	   
	Real neighborDistance_stnd;
	Real depth_stnd;
	if (badNeighborDistance (neighborDistance_stnd, depth_stnd))
	{
		ONumber oNum (os, 1, false);  // PAR
	  os << "  distNghb=" << neighborDistance_stnd << "(depth=" << depth_stnd << ")";
	}

  const string targetFeaturesChange (getTargetFeaturesChange ());
  if (! targetFeaturesChange. empty ())
    os << "  " << targetFeaturesChange;

	if (verbose ())
	{
		ONumber oNum (os, 1, false);  // PAR
	  os << "  dist=" << getDistance ();
	}
#endif
	
#if 0
	os << fixed << setprecision (3);  // PAR
	FOR (unsigned char, j, 2)
	  FOR (unsigned char, i, 2)
	    os << " w[" << (int) i << "][" << (int) j << "]=" << weight [i] [j];
#endif
}



#if 0
bool Phyl::getSaveSubtreeP () const 
{ 
  return ! getFeatureTree (). savePhenChangesOnly || hasPhenChange;
}
#endif



const FeatureTree& Phyl::getFeatureTree () const
{
  return static_cast <const FeatureTree&> (getTree ());
}



Real Phyl::feature2weight (size_t /*featureIndex ??*/,
	                         bool thisCore,
	                         bool parentCore) const
{ 
	const Real w = /*getFeatureTree (). features [featureIndex].*/ weight [thisCore] [parentCore]; 
	ASSERT (w >= 0);
	return w;
}



bool Phyl::feature2parentCore (size_t featureIndex) const
{ 
  const Species* s = static_cast <const Species*> (getParent ());
  return s ? s->core [featureIndex] : getFeatureTree (). getRootCore (featureIndex);
}



size_t Phyl::getCoreSize () const
{
  ASSERT (getFeatureTree (). coreSynced);   

  size_t n = getFeatureTree (). commonCore. size ();
  FOR (size_t, i, getFeatureTree (). genes)
    if (core [i])
    	n++;
  return n;
}



size_t Phyl::getCoreChange (bool gain) const
{ 
	ASSERT (getFeatureTree (). coreSynced);

	size_t d = 0;
	FOR (size_t, i, getFeatureTree (). genes)
	  if (   core [i]               == gain
	  	  && feature2parentCore (i) != gain
	  	 )
	  	d++;
	return d;
}



Real Phyl::getDistance () const
{ 
	ASSERT (getFeatureTree (). coreSynced);

	Real d = getPooledDistance ();
	FOR (size_t, i, getFeatureTree (). genes)
  	d += feature2weight (i, core [i], feature2parentCore (i));
  ASSERT (d >= 0);
#ifndef NDEBUG
  if (d == INF)
  {
  	cout << getName () << "  " << endl;
  	saveContent (cout);
  	cout << endl;
  	FOR (unsigned char, i, 2)
	  	FOR (unsigned char, j, 2)
	  	  cout << "weight [" << (int) i << "][" << (int) j << "] = " << weight [i] [j] << endl;
  	ERROR;
  }
#endif
  
	return d;
}



Real Phyl::getSubtreeLength () const
{ 
	Real d = getDistance ();
	for (const DiGraph::Arc* arc : arcs [false])
	  d += static_cast <Phyl*> (arc->node [false]) -> getSubtreeLength ();
	return d;
}	  



void Phyl::getParent2corePooled (size_t parent2corePooled [2/*thisCore*/] [2/*parentCore*/]) const
{
  FOR (unsigned char, thisCore, 2)
    FOR (unsigned char, parentCore, 2)
      parent2corePooled [thisCore] [parentCore] = 0;
  parent2corePooled [true] [getParent () ? true : ! FeatureTree::emptyRoot] = getFeatureTree (). commonCore. size ();
  parent2corePooled [false] [false] = getFeatureTree (). globalSingletonsSize;
}



Real Phyl::getPooledDistance () const
{ 
  size_t parent2corePooled [2/*thisCore*/] [2/*parentCore*/];
  getParent2corePooled (parent2corePooled);

  Real d = 0;
  FOR (unsigned char, thisCore, 2)
    FOR (unsigned char, parentCore, 2)
      d += multiplyLog ((Real) parent2corePooled [thisCore] [parentCore], feature2weight (thisCore, parentCore));
  
  if (! (d >= 0))
  {
    cout << getName () << " " << d << endl;
    FOR (unsigned char, thisCore, 2)
      FOR (unsigned char, parentCore, 2)
        cout << "parent2corePooled [" << (int) thisCore << "] [" << (int) parentCore << "] = " << parent2corePooled [thisCore] [parentCore] << endl;
    FOR (unsigned char, thisCore, 2)
      FOR (unsigned char, parentCore, 2)
        cout << "feature2weight ["    << (int) thisCore << "] [" << (int) parentCore << "] = " << feature2weight (thisCore, parentCore) << endl;
    ERROR;
  }
  return d;
}



Real Phyl::getNeighborDistance () const
{
	Real d = getDistance ();
	for (const DiGraph::Arc* arc: arcs [false])
	  d += static_cast <Phyl*> (arc->node [false]) -> getDistance ();
	return d;
}



bool Phyl::badNeighborDistance (Real &neighborDistance_stnd,
                                Real &depth_stnd) const
{
  neighborDistance_stnd = NAN;
  depth_stnd = NAN;
  const Normal* distDistr = getFeatureTree (). distDistr. get ();
  if (! distDistr)
    return false;
  neighborDistance_stnd = distDistr->stnd (getNeighborDistance ());
  depth_stnd = getFeatureTree (). depthDistr->stnd ((Real) getDepth ());
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



#ifndef NDEBUG
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

	Real childrenCoreDistance [2/*bool thisCore*/];
	FOR (unsigned char, thisCore, 2)
	{
    childrenCoreDistance [thisCore] = 0;
		for (const DiGraph::Arc* arc : arcs [false])
			childrenCoreDistance [thisCore] += static_cast <Phyl*> (arc->node [false]) -> parent2core [thisCore] [featureIndex]. treeLen;
  }
  
	FOR (unsigned char, parentCore, 2)
	{
		Real distance [2/*bool thisCore*/];
		FOR (unsigned char, thisCore, 2)
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
  FOR (size_t, i, parent2core [false]. size ())
    assignFeature (i);
}



void Phyl::assignFeaturesDown ()
{ 
	for (DiGraph::Arc* arc : arcs [false])
		static_cast <Phyl*> (arc->node [false]) -> assignFeaturesDown ();
  assignFeatures ();
}



#if 0
string Phyl::getPhenChange (bool skipSingleNonSystems) const
{ 
	ASSERT (getFeatureTree (). coreSynced);
	ASSERT (getFeatureTree (). distDistr. get ());

	string s;
	FOR_START (size_t, i, getFeatureTree (). genes, core. size ())
	  if (core [i] != feature2parentCore (i))
	  {
	    const Feature& f (getFeatureTree (). features [i]);
	    f. qc ();
	    if (f. genomes == 1)
	    {
        const Strain* st = asStrain ();
	      const Genome* g = asGenome ();
	      if (st)
	        g = st->getGenome ();
	      if (getFeatureTree (). allTimeZero)
	      {
  	      ASSERT (   st 
  	              && core [i] 
  	              &&   st->getGenome () -> core [i] 
  	              && ! st->getGenome () -> optionalCore [i]
  	             );
  	    }
  	    else 
	      {
 	        ASSERT (g);
  	      if (! (core [i] && ! g->optionalCore [i]))
  	      {
  	        cout << getName () << " " << (int) core [i] << endl;
  	        cout << " " << (int) g->optionalCore [i] << endl;
  	        ERROR;
  	      }
  	    }
        ASSERT (g);
	      if (skipSingleNonSystems && ! g->isSystemPhen (f. name))
	        continue;
	    }
	    if (! s. empty ())
	      s += ";";
	  	s += f. name;
	  	if (! core [i])
	  	  s += "-";
	  }
	return s;
}



string Phyl::getTargetFeatureChange (const TargetFeature &tf) const
{ 
  if (   tf. index == NO_INDEX
      || core [tf. index] == feature2parentCore (tf. index)
     )
    return string ();
  return tf. getName () + (core [tf. index] ? "" : "-");
}



string Phyl::getTargetFeaturesChange () const
{ 
  string s;
  for (const TargetFeature& tf : getFeatureTree (). targetFeatures)
  {
    const string c (getTargetFeatureChange (tf));
    if (c. empty ())
      continue;
    if (! s. empty ())
      s += ",";
    s += c;
  }
  return s;
}



void Phyl::setHasPhenChange ()
{
  const string phenChange (! getFeatureTree (). targetFeatures. empty ()
                             ? getTargetFeaturesChange ()
                             : getPhenChange (true)
                          );
  hasPhenChange = ! phenChange. empty ();
	for (DiGraph::Arc* arc : arcs [false])
	{
	  Phyl* s = static_cast <Phyl*> (arc->node [false]);
	  s->setHasPhenChange ();
	  maximize (hasPhenChange, s->hasPhenChange);
	}
}
#endif




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
, time (tree. allTimeZero ? NAN : time_arg)
, pooledSubtreeDistance (INF)
, time_old (NAN)
, pooledSubtreeDistance_old (NAN)
, movementsOn (false)
{
	FOR (unsigned char, i, 2)
		FOR (unsigned char, j, 2)
		  weight_old [i] [j] = NAN;
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



string Species::getNewickName (bool /*minimal*/) const
{
	ASSERT (getFeatureTree (). coreSynced);
	ASSERT (getFeatureTree (). distDistr. get ());

  if (! getParent ())
    return string ();

	const Feature* bestF = nullptr;
	FOR_START (size_t, i, getFeatureTree (). genes, core. size ())
	  if (core [i] && ! feature2parentCore (i))
	  {
	    const Feature& f (getFeatureTree (). features [i]);
	    f. qc ();
	    if (f. better (bestF))
	      bestF = & f;
	  }

	return bestF ? bestF->name : string ();  
}



void Species::setWeight ()
{
  if (getFeatureTree (). allTimeZero || ! getFeatureTree (). genesExist)
		FOR (unsigned char, parentCore, 2)
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
		FOR (unsigned char, parentCore, 2)
	  {
	  	const Prob changeProb = negateProb (getFeatureTree (). lambda0, parentCore) * a;
			ASSERT (isProb (changeProb));
		  FOR (unsigned char, thisCore, 2)
		  {
		  	const Prob p = negateProb (changeProb, parentCore == thisCore); 
		  	const Real w = - log (p);
		  	ASSERT (w >= 0);
		    weight [thisCore] [parentCore] = w;
		  }
		}
	}
}



void Species::setCore ()
{
	FOR (size_t, i, core. size ())
	  core [i] = feature2core (i);
	for (DiGraph::Arc* arc : arcs [false])
	  static_cast <Phyl*> (arc->node [false]) -> setCore ();
}



#if 0
void Species::saveDatabaseTopology (Server &db,
	                                  bool &parent_STRAIN)
{
	ASSERT (getFeatureTree (). coreSynced);
	
	const_cast <FeatureTree&> (getFeatureTree ()). progInternal ();

	
  bool this_STRAIN = true;
  if (const Fossil* f = static_cast <const Fossil*> (getParent ()))
  {
  	// SPECIES, STRAIN
  	ASSERT (f->id);
	  CQuery query (db.NewQuery());
    query.SetParameter("@id", 0, eSDB_Int4, eSP_InOut);
	  if (parent_STRAIN)
	  {
		  query.SetParameter("@strain_id", (int) f->id);
		  query.ExecuteSP("STRAIN_split");
		}
	  else
	  {
		  query.SetParameter("@species_id", (int) f->id);
		  query.ExecuteSP("SPECIES_split");
		}
	  id = (uint) query.GetParameter("@id").AsInt4();
	  parent_STRAIN = false;
  }
  else
  {
    ASSERT (id);
	  CQuery query (db.NewQuery());
	  query.SetSql("select null from STRAIN where id = @id");
	  query.SetParameter("@id", (int) id);
 	  this_STRAIN = exists (query);
	}
  ASSERT (id);


#ifndef NDEBUG
  uint prevId = id;
#endif  
	for (const DiGraph::Arc* arc : arcs [false])
	{
		Phyl* s = static_cast <Phyl*> (arc->node [false]);
	  s->saveDatabaseTopology (db, this_STRAIN);
	#ifndef NDEBUG
	  if (const Fossil* f = s->asFossil ())
	  {
		  ASSERT (f->id > prevId);
		  prevId = f->id;
		}
	#endif  
	}
}



void Species::saveDatabaseNode (Server &db) const
{
	ASSERT (getFeatureTree (). coreSynced);
	ASSERT (id);

	
	FeatureTree& tree = const_cast <FeatureTree&> (getFeatureTree ());
	
	tree. progInternal (getName ());
  

  // SPECIES  
  {
	  CQuery speciesQuery (db.NewQuery());
	  speciesQuery.SetSql("\
update SPECIES \
  set core_size = @core_size \
    , [time] = @time \
  where id = @id");
	  speciesQuery.SetParameter("@id", (int) id);
	  speciesQuery.SetParameter("@core_size", (int) getCoreSize ());
	  speciesQuery.SetParameter("@time", time);
	  speciesQuery.Execute();
	  ASSERT (speciesQuery.GetRowCount() == 1);
	}


  // SPECIES_PROT
  {
    CQuery insQuery (db.NewQuery());
    insQuery.SetSql("\
insert into SPECIES_PROT \
         ( species,  prot, gain) \
  values (@species, @prot,    1)");
  	insQuery.SetParameter("@species", (int) id);
  	//
    CQuery delQuery (db.NewQuery());
    delQuery.SetSql("\
update SPECIES_PROT \
  set gain = 0 \
  where     species = @species \
        and prot    = @prot");
  	delQuery.SetParameter("@species", (int) id);
  	if (const Fossil* parent_ = static_cast <const Fossil*> (getParent ()))
  	{
  		ASSERT (parent_->id);
  		{
  		  CQuery borrowQuery (db.NewQuery());
  		  borrowQuery.SetSql("\
insert into SPECIES_PROT ( species, prot) \
  select                           @species, prot \
    from SPECIES_PROT \
    where     species = @parent_id \
          and isnull(gain,1) = 1");
  			borrowQuery.SetParameter("@species", (int) id);
  			borrowQuery.SetParameter("@parent_id", (int) parent_->id);
  	    borrowQuery.Execute();
  		}
  	  FOR (size_t, i, tree. genes)
  	    if (core [i] && ! parent_->core [i])
  	    {
  		    insQuery.SetParameter("@prot", tree. features [i]. geneId ());
  	      insQuery.Execute();
  	    }
  	    else if (! core [i] && parent_->core [i])
  	    {
  		    delQuery.SetParameter("@prot", tree. features [i]. geneId ());
  	      delQuery.Execute();
  	    }
  	}
  	else
  	{
  		for (const Feature::Id& fId : tree. commonCore)
  	  {
  	    insQuery.SetParameter("@prot", str2<int> (fId));
  	    insQuery.Execute();
  	  }
  	  FOR (size_t, i, tree. genes)
  	    if (core [i])
  	    {
  		    insQuery.SetParameter("@prot", tree. features [i]. geneId ());
  	      insQuery.Execute();
  	    }
  	}
  }


	for (const DiGraph::Arc* arc : arcs [false])
		static_cast <Phyl*> (arc->node [false]) -> saveDatabaseNode (db);
}



void Species::saveDatabasePhen (Server &db) const
{
	ASSERT (getFeatureTree (). coreSynced);
	ASSERT (id);

	
	FeatureTree& tree = const_cast <FeatureTree&> (getFeatureTree ());
	
	tree. progInternal (getName ());
  

  // SPECIES_PHEN
  {
    CQuery query (db.NewQuery());
    query.SetSql("delete from SPECIES_PHEN where species = @species");
  	query.SetParameter("@species", (int) id);
  	query.Execute();
  }
  // Cf. SPECIES_PROT
  {
    CQuery insQuery (db.NewQuery());
    insQuery.SetSql("\
insert into SPECIES_PHEN \
         ( species,  phen, gain) \
  values (@species, @phen,    1)");
  	insQuery.SetParameter("@species", (int) id);
  	//
    CQuery delQuery (db.NewQuery());
    delQuery.SetSql("\
update SPECIES_PHEN \
  set gain = 0 \
  where     species = @species \
        and phen    = @phen");
  	delQuery.SetParameter("@species", (int) id);
  	if (const Fossil* parent_ = static_cast <const Fossil*> (getParent ()))
  	{
  		ASSERT (parent_->id);
  		{
  		  CQuery borrowQuery (db.NewQuery());
  		  borrowQuery.SetSql("\
insert into SPECIES_PHEN ( species, phen) \
  select                  @species, phen \
    from SPECIES_PHEN \
    where     species = @parent_id \
          and isnull(gain,1) = 1");
  			borrowQuery.SetParameter("@species", (int) id);
  			borrowQuery.SetParameter("@parent_id", (int) parent_->id);
  	    borrowQuery.Execute();
  		}
  	  FOR_START (size_t, i, tree. genes, core. size ())
  	    if (core [i] && ! parent_->core [i])
  	    {
  		    insQuery.SetParameter("@phen", tree. features [i]. name);
  	      insQuery.Execute();
  	    }
  	    else if (! core [i] && parent_->core [i])
  	    {
  		    delQuery.SetParameter("@phen", tree. features [i]. name);
  	      delQuery.Execute();
  	    }
  	}
  	else
  	  FOR_START (size_t, i, tree. genes, core. size ())
  	    if (core [i])
  	    {
  		    insQuery.SetParameter("@phen", tree. features [i]. name);
  	      insQuery.Execute();
  	    }
  }


	for (const DiGraph::Arc* arc : arcs [false])
		static_cast <Phyl*> (arc->node [false]) -> saveDatabasePhen (db);
}
#endif



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
   	  FOR (unsigned char, i, 2)
   	    FOR (unsigned char, j, 2)
   	      parent2core [i] [j] = (Real) parent2core_arg [i] [j];
   	//parent2core [false] [false] +=        tree->commonMissings      * (Real) nodes;
	  //parent2core [true] [true]   += (Real) tree->commonCore. size () /* * (Real) nodes*/;
    }

  
  Real f (Real x) 
  { 
    ASSERT (x >= 0);
    
  	const Prob a = 1 - exp (- x);
  	
  	Real lhs = 0;
  	FOR (unsigned char, c, 2)
  		if (const Prob denomin = a)
  	    lhs += parent2core [! c] [c] / denomin;
  	  else
  	  	lhs = INF;
  	ASSERT (lhs >= 0);
  	  
  	Real rhs = 0;
  	FOR (unsigned char, c, 2)
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
  	return NAN;

  if (! getParent ())
  	return getFeatureTree (). getRootTime ();;  
  
  if (! getFeatureTree (). timeOptimFrac)  	
    return getFeatureTree (). time_init;  

  size_t parent2core_ [2/*thisCore*/] [2/*parentCore*/];
  getParent2corePooled (parent2core_);	  
	FOR (size_t, i, getFeatureTree (). genes)
	  parent2core_ [core [i]] [feature2parentCore (i)] ++;
	  
  const Real t = coreChange2time (getFeatureTree (), parent2core_ /*, 1*/);
  ASSERT (t >= 0);
  
  return exp (convexCombination (getFeatureTree (). timeOptimFrac, log (t), log (getFeatureTree (). time_init)));
}



void Species::assignTime ()
{ 
	IMPLY (! getFeatureTree (). allTimeZero, ! isNan (time));
	
	time_old = time;
	FOR (unsigned char, i, 2)
		FOR (unsigned char, j, 2)
		  weight_old [i] [j] = weight [i] [j];

  setTimeWeight ();
}



void Species::restoreTime ()
{
	IMPLY (! getFeatureTree (). allTimeZero, ! isNan (time_old));
	
	time = time_old;
	FOR (unsigned char, i, 2)
		FOR (unsigned char, j, 2)
		  weight [i] [j] = weight_old [i] [j];
	
	commitTime ();
}



Real Species::feature2treeLength (size_t featureIndex) const
{ 
  // Cf. feature2parentCore()
  return FeatureTree::emptyRoot 
           ? parent2core [false] [featureIndex]. treeLen
           : getFeatureTree (). allTimeZero || ! getFeatureTree (). genesExist
           	     ? min ( parent2core [false] [featureIndex]. treeLen 
      		             , parent2core [true]  [featureIndex]. treeLen
      		             )
  		           : parent2core [getFeatureTree (). rootCore [featureIndex]] [featureIndex]. treeLen;
}



Real Species::root2treeLength () const
{
  Real s = pooledSubtreeDistance;
  FOR (size_t, i, getFeatureTree (). genes)
    s += feature2treeLength (i);
  return s;
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
  FOR (unsigned char, parentCore, 2)
    FOR (size_t, featureIndex, parent2core [parentCore]. size ())
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
	pooledSubtreeDistance_old = NAN;
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



Real Fossil::getPooledSubtreeDistance () const 
{ 
  Real d = getPooledDistance ();
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
	FOR (size_t, i, core. size ())
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
	Real distance [2/*bool thisCore*/];
	const Genome* g = getGenome ();
	FOR (unsigned char, thisCore, 2)
	  distance [thisCore] =   g->feature2weight (true,     thisCore)
	                        +    feature2weight (thisCore, false);
  singletonsInCore = leReal (distance [true], distance [false]);  

  Species::assignFeatures ();
}



Real Strain::getPooledSubtreeDistance () const 
{ 
  return getPooledDistance () + getGenome () -> getPooledDistance (); 
}



#if 0
string Strain::getTargetFeatureChange (const TargetFeature &tf) const
{ 
  if (tf. index == NO_INDEX)
    if (singletonsInCore && getGenome () -> singletons. contains (tf. featureId))
      return tf. getName ();
    else
      return string ();
  else
    return Species::getTargetFeatureChange (tf);
}
#endif



size_t Strain::orphans () const
{
  ASSERT (getFeatureTree (). coreSynced);
  
  size_t n = 0;
  if (const Fossil* p = static_cast <const Fossil*> (getParent ()))
  {
    const Genome* g = getGenome ();
    FOR (size_t, i, getFeatureTree (). genes)
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
, hasPhens (true)
, coreNonSingletons (0)
#if 0
, L50 (0)
, oddCdss (0)
, project_id (0)
, tax_id (0)
, pubmed (0)
, outbreakYear (0)
#endif
{ 
	ASSERT (parent_arg);
	ASSERT (! id. empty ()); 
	ASSERT (! contains (id, ' '));
}



#if 0
void Genome::initDb (Server &db)
{
  ASSERT (coreSet. empty ());
  ASSERT (id);
  ASSERT (! coreNonSingletons);

	{
	  CQuery query (db.NewQuery());
		query.SetSql ("\
select /* 1*/ GENOME.serovar \
     , /* 2*/ GENOME.pathovar \
     , /* 3*/ GENOME.odd_cdss \
     , /* 4*/ GENOME.project_id \
     , /* 5*/ GENOME.tax \
     , /* 6*/ TAX.name \
     , /* 7*/ SUBSP.name \
     , /* 8*/ SP.name \
     , /* 9*/ GENUS.name \
     , /*10*/ COUNTRY.continent \
     , /*11*/ GENOME.L50 \
     , /*12*/ GENOME.sequencer \
     , /*13*/ GENOME.pubmed \
     , /*14*/ OUTBREAK.name \
     , /*15*/ OUTBREAK.year \
     , /*16*/ GENOME.class \
  from           GENOME \
            join TAX on TAX.id = GENOME.tax \
            join GENOMETAX on GENOMETAX.id = GENOME.tax \
       left join TAXSPECIES on TAXSPECIES.id = GENOMETAX.taxspecies \
       left join TAX SUBSP on SUBSP.id = GENOMETAX.subspecies \
       left join TAX SP on SP.id = GENOMETAX.taxspecies \
       left join TAX GENUS on GENUS.id = TAXSPECIES.genus \
       left join COUNTRY on COUNTRY.id = GENOME.country \
       left join OUTBREAK on OUTBREAK.id = GENOME.outbreak \
  where GENOME.id = @id");
	  query.SetParameter("@id", (int) id);
	  query.Execute();
	  ITERATE(CQuery, row, query) 
	  {
	    ASSERT (! project_id);
      serovar           =          row[ 1].AsString();
      pathovar          =          row[ 2].AsString();
      oddCdss           = (size_t) row[ 3].AsInt4();  
      project_id        = (uint)   row[ 4].AsInt4();
      tax_id            = (uint)   row[ 5].AsInt4();
      taxName           =          row[ 6].AsString();
      subspecies        =          row[ 7].AsString();
      taxSpecies        =          row[ 8].AsString();
      genus             =          row[ 9].AsString();
      continent         =          row[10].AsString();
      L50               = (uint)   row[11].AsInt4();
      sequencer         =          row[12].AsString();
      pubmed            = (uint)   row[13].AsInt4();
      outbreakName      =          row[14].AsString();
      outbreakYear      = (uint)   row[15].AsInt4();
      phylogeneticClass =          row[16].AsString();
    }
  }
  
  // subspecies
  if (! taxSpecies. empty () && isLeft (subspecies, taxSpecies + " "))
    subspecies. erase (0, taxSpecies. size () + 1);

  // taxSpecies    
  if (! genus. empty () && isLeft (taxSpecies, genus + " "))
    taxSpecies. erase (0, genus. size () + 1);

#if 0
  // subspecies
  const string subsp ("subsp. ");
  if (isLeft (subspecies, subsp))
    subspecies. erase (0, subsp. size ());
  if (! subspecies. empty () && subspecies == taxSpecies)
    subspecies = subsp + subspecies;
#endif
    
  // isolate
  if (project_id)
	{
	  CQuery query (db.NewQuery());
		query.SetSql ("\
select distinct isolate \
  from BIOPROJECT \
  where     gb = @gb \
        and isolate is not null");
	  query.SetParameter("@gb", (int) project_id);
	  query.Execute();
	  ITERATE(CQuery, row, query) 
	  {
	    ASSERT (isolate. empty ());
	    isolate = row[1].AsString();
	  }
  }

  // coreSet
  // Load partial CDSs as optionalCore[] = true ??
  {
	  CQuery query (db.NewQuery());
		query.SetSql ("\
select distinct ORF.prot \
  from      ORF \
       join GDNA on GDNA.nucl_gi = ORF.gdna \
  where     GDNA.genome = @genome \
        and ORF.consistent = 1 \
  option (loop join)");
	  query.SetParameter("@genome", (int) id);
	  query.Execute();
	  ITERATE(CQuery, row, query) 
	    coreSet << toString (row[1].AsInt4());
	}

  coreNonSingletons = coreSet. size ();  // Includes singletons and common core
  
  // phens
  {
	  CQuery query (db.NewQuery());
		query.SetSql ("\
select GENOME_PHEN.phen, GENOME_PHEN.optional \
  from           GENOME_PHEN \
            join GENOME on GENOME.id = GENOME_PHEN.genome \
       left join PAN_PHEN on     PAN_PHEN.pan = GENOME.pan \
                             and PAN_PHEN.phen = GENOME_PHEN.phen \
  where     GENOME_PHEN.genome = @genome \
        and (   PAN_PHEN.pan is null \
             or GENOME_PHEN.optional = 0 \
            )");
	  query.SetParameter("@genome", (int) id);
	  query.Execute();
	  ITERATE(CQuery, row, query) 
	    addPhen (row[1].AsString(), 
	             row[2].AsInt4());
	}
  {
	  CQuery query (db.NewQuery());
		query.SetSql ("\
select PAN_PHEN.phen \
  from           GENOME \
            join PAN_PHEN on PAN_PHEN.pan = GENOME.pan \
       left join GENOME_PHEN on     GENOME_PHEN.genome = GENOME.id \
                                and GENOME_PHEN.phen = PAN_PHEN.phen \
  where     GENOME.id = @genome \
        and GENOME_PHEN.genome is null");
	  query.SetParameter("@genome", (int) id);
	  query.Execute();
	  ITERATE(CQuery, row, query) 
	    addPhen (row[1].AsString(), true);
	}
	addSystemPhens ();
}
#endif



void Genome::initDir (const string &geneDir,
                      const string &phenDir)
{
  ASSERT (coreSet. empty ());
  ASSERT (! id. empty ());
  ASSERT (! coreNonSingletons);

 
  // geneDir
  {
    ASSERT (! geneDir. empty ());
    ifstream f (geneDir + "/" + id);
    ASSERT (f. good ());
    // coreSet
    for (;;)
    {
    	Feature::Id geneName;
    	f >> geneName;
  	  if (f. eof ())
  	  	break;
    	ASSERT (! geneName. empty ());
      coreSet << geneName; 
    }
  }
  if (verbose ())
    cout << id << ": core=" << coreSet. size () << endl;


  // phenDir
  if (! phenDir. empty ())
  {
    ifstream f (phenDir + "/" + id);
    ASSERT (f. good ());

    const size_t len = 1024;  // PAR
    char s [len];
  
  #if 0
    #define READ(P)   { f. getline (s, len); stringstream ss (string (s). substr (strlen (#P "="))); ss >> P; }
    #define READS(P)  { f. getline (s, len); stringstream ss (string (s). substr (strlen (#P "="))); P = ss. str (); trim (P); }
    READS (subspecies);
    READS (taxSpecies);
    READS (genus);
    READS (serovar);
    READS (pathovar);
    READS (continent);
    READ  (oddCdss);
    READ  (project_id);
    READ  (tax_id);
    READS (taxName);  
    READ  (L50);
    READS (sequencer);
    READ  (pubmed);
    READS (outbreakName);
    READ  (outbreakYear);
    READS (phylogeneticClass);
    READS (isolate);
    #undef READ
    #undef READS
    f. getline (s, len);
    ASSERT (strlen (s) == 0);
  #endif
    
    // phens
    while (! f. eof ())
    {
      f. getline (s, len);
      stringstream ss (s);
      Feature::Id phen;
      bool optional;
      ss >> phen;
      if (phen == "nophenotypes")
      {
        hasPhens = false;
        continue;
      }
      ss >> optional;
      addPhen (phen, optional);
    }
  #if 0
    addSystemPhens ();
  #endif
  }

  
  coreNonSingletons = coreSet. size ();  // Includes singletons and common core
}



void Genome::init (const map <Feature::Id, size_t/*index*/> &feature2index)
{
	IMPLY ( ! coreSet. empty (), getFeatureTree (). genesExist);

	Phyl::init ();

	ASSERT (optionalCore. empty ());
	optionalCore. resize (core. size (), false); 

  // coreSet
  // optionalCore[] ??
	for (const Feature::Id& fId : coreSet)
	{
		ASSERT (feature2index. find (fId) != feature2index. end ());
		const size_t featureIndex = feature2index. at (fId);
		FOR (unsigned char, parentCore, 2)
	    parent2core [parentCore] [featureIndex]. core = ETRUE;
	  core [featureIndex] = true;
	}

  // phens
  if (hasPhens)
  	for (const auto it : phens)
  	{
  		ASSERT (feature2index. find (it. first) != feature2index. end ());
  		const size_t featureIndex = feature2index. at (it. first);
  		FOR (unsigned char, parentCore, 2)
  	    parent2core [parentCore] [featureIndex]. core = ETRUE;
  	  core [featureIndex] = true;
  	  if (it. second)
  	    optionalCore [featureIndex] = true;
  	}
  else
  {
    ASSERT (phens. empty ());
    FOR_START (size_t, i, getFeatureTree (). genes, optionalCore.size ())
      optionalCore [i] = true;
  }

	coreSet. clear ();
	phens. clear ();
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
  ASSERT (coreNonSingletons >= getFeatureTree (). commonCore. size ());
	ASSERT (coreSet. empty ());
	IMPLY (getFeatureTree (). coreSynced, coreNonSingletons >= getCoreSize ());
	ASSERT (phens. empty ());
	
	ASSERT (! singletons. intersects (getFeatureTree (). commonCore));
	
	IMPLY (! hasPhens, phens. empty ());

#if 0
//ASSERT (project_id);
//ASSERT (tax_id);
//ASSERT (! taxName. empty ());
	ASSERT (outbreakName. empty () == ! outbreakYear);
//ASSERT (L50);
	
//ASSERT (! contains (serovar, ":"));
	ASSERT (! contains (pathovar, ":"));
#endif
}



void Genome::saveContent (ostream& os) const
{ 
  Phyl::saveContent (os);  

	os << "  S=" << (getStrain () -> singletonsInCore ? 0 : singletons. size ()) /*<< '+' << oddCdss*/;    
	if (const size_t opt = coreNonSingletons - getCoreSize ())
	  os << "  OptC=" << opt;
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
  if (getFeatureTree (). allTimeZero || ! getFeatureTree (). genesExist)
		FOR (unsigned char, parentCore, 2)
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
    weight [true]  [true]  = - log (1 - annotError);
    weight [false] [true]  = - log (annotError);
    weight [true]  [false] = - log (misannotError);
    weight [false] [false] = - log (1 - misannotError);
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
	FOR (size_t, i, core. size ())
    if (optionalCore [i])
	    core [i] = feature2core (i);
}



void Genome::assignFeature (size_t featureIndex)
{
  if (optionalCore [featureIndex])
    Phyl::assignFeature (featureIndex);
  else
  	FOR (unsigned char, parentCore, 2)
    {
    	const ebool core_old = parent2core [parentCore] [featureIndex]. core;
    	ASSERT (core_old != UBOOL);
    	setCoreEval (featureIndex, parentCore, CoreEval (feature2weight (featureIndex, core_old, parentCore), core_old));
    }
}



#if 0
string Genome::getTargetFeatureChange (const TargetFeature &tf) const
{ 
  if (tf. index == NO_INDEX)
    if (! getStrain () -> singletonsInCore && singletons. contains (tf. featureId))
      return tf. getName ();
    else
      return string ();
  else
    return Phyl::getTargetFeatureChange (tf);
}
#endif



#if 0
void Genome::saveDatabaseTopology (Server &db,
	                                 bool &parent_STRAIN)
{
	ASSERT (getFeatureTree (). coreSynced);
	ASSERT (parent_STRAIN);
  const Strain* f = getStrain ();
 	ASSERT (f->id);

	const_cast <FeatureTree&> (getFeatureTree ()). progInternal ();

  CQuery query (db.NewQuery());
  query.SetParameter("@genome", (int) id);
  query.SetParameter("@strain", (int) f->id);
  query.ExecuteSP("GENOME_strain");
}



void Genome::saveDatabaseNode (Server &db) const
{
	ASSERT (getFeatureTree (). coreSynced);
	ASSERT (id);
	
	FeatureTree& tree = const_cast <FeatureTree&> (getFeatureTree ());	
	tree. progInternal (getName ());  

	const Strain* st = getStrain ();
	ASSERT (st->id);

  // GENOME
  {
    CQuery query (db.NewQuery());
    query.SetSql("\
update GENOME \
  set missed_prots = @missed_prots \
    , extra_prots = @extra_prots \
    , singleton_prots = @singleton_prots \
    , singleton_prots_in_core = @singleton_prots_in_core \
  where id = @id");
    query.SetParameter("@id", (int) id);
    query.SetParameter("@missed_prots", (int) getCoreChange (false));
    query.SetParameter("@extra_prots",  (int) getCoreChange (true));
    query.SetParameter("@singleton_prots",  (int) singletons. size ());
    query.SetParameter("@singleton_prots_in_core",  (int) st->singletonsInCore);
    query.Execute();
    ASSERT (query.GetRowCount() == 1);
  }
  
  // STRAIN_PROT 
  {
    CQuery query (db.NewQuery());
    query.SetSql("\
insert into STRAIN_PROT \
         ( strain,  prot,  gain) \
  values (@strain, @prot, @gain)");
  	query.SetParameter("@strain", (int) st->id);
    FOR (size_t, i, tree. genes)
      if (core [i] != st->core [i])
      {
  	    query.SetParameter("@prot", tree. features [i]. geneId ());
  	    query.SetParameter("@gain", core [i]);
        query.Execute();
      }
    query.SetParameter("@gain", 1);
    for (const Feature::Id& fId : singletons)
    {
	    query.SetParameter("@prot", Feature::id2geneId (fId));
      query.Execute();
    }
  }
}



void Genome::saveDatabasePhen (Server &db) const
{
	ASSERT (getFeatureTree (). coreSynced);
	ASSERT (id);
	
	FeatureTree& tree = const_cast <FeatureTree&> (getFeatureTree ());	
	tree. progInternal (getName ());  

	const Strain* st = getStrain ();
	ASSERT (st->id);

  // STRAIN_PHEN
  {
    CQuery query (db.NewQuery());
    query.SetSql("delete from STRAIN_PHEN where strain = @strain");
  	query.SetParameter("@strain", (int) st->id);
  	query.Execute();
  }
  {
    CQuery query (db.NewQuery());
    query.SetSql("\
insert into STRAIN_PHEN \
         ( strain,  phen,  gain) \
  values (@strain, @phen, @gain)");
  	query.SetParameter("@strain", (int) st->id);
	  FOR_START (size_t, i, tree. genes, core. size ())
      if (core [i] != st->core [i])
      {
  	    query.SetParameter("@phen", tree. features [i]. name);
  	    query.SetParameter("@gain", core [i]);
        query.Execute();
      }
  }
}



void Genome::saveFeatures (const string &dir) const
{
	ASSERT (! dir. empty ());
	
	OFStream f (dir, toString (id), "");

  #define WRITE(P)  f << #P "=" << P << endl
  WRITE (subspecies);
  WRITE (taxSpecies);
  WRITE (genus);
  WRITE (serovar);
  WRITE (pathovar);
  WRITE (continent);
  WRITE (oddCdss);
  WRITE (project_id);
  WRITE (tax_id);
  WRITE (taxName);  
  WRITE (L50);
  WRITE (sequencer);
  WRITE (pubmed);
  WRITE (outbreakName);
  WRITE (outbreakYear);
  WRITE (phylogeneticClass);
  WRITE (isolate);
  #undef WRITE
  f << endl;
	
	// Phenotypes
	FOR_START (size_t, i, getFeatureTree (). genes, core. size ())  
	{
	  const Feature::Id phenName (getFeatureTree (). features [i]. name);
	  ASSERT (! singletons. contains (phenName));
	  if (isSystemPhen (phenName))
	    continue;
	  if (optionalCore [i])
	  	f << phenName << " 1" << endl;
	  else if (core [i])
	  	f << phenName << " 0" << endl;
	}
  f << endl;

  // Genes
  // optionalCore[] ??
	FOR (size_t, i, getFeatureTree (). genes)
 	{
	  const Feature::Id geneName (getFeatureTree (). features [i]. name);
	  ASSERT (! singletons. contains (geneName));
	  if (core [i])
	  	f << geneName << endl;
	}
	for (const Feature::Id& fId : getFeatureTree (). commonCore)
		f << fId << endl;
	for (const Feature::Id& fId : singletons)
		f << fId << endl;
}
#endif



void Genome::addPhen (const Feature::Id &phen,
                      bool optional)
{
  if (phen. empty ())
    return;
  ASSERT (phen [0]);
  ASSERT (phen [phen. size () - 1]);
  if (phens. find (phen) != phens. end ())
    ERROR_MSG ("Phenotype \"" + phen + "\" already exists in genome " + getName ());
  phens [phen] = optional;
}



void Genome::getSingletons (Set<Feature::Id> &globalSingletons,
                            Set<Feature::Id> &nonSingletons) const
{ 
	if (getFeatureTree (). genesExist && coreSet. empty ())
	{
	  cout << getName () << endl;  
	  ERROR;
	}
	
	for (const Feature::Id& fId : coreSet)
	  if (nonSingletons. contains (fId))
	  	;
	  else if (globalSingletons. contains (fId))
	  	setMove (& globalSingletons, & nonSingletons, fId);
	  else
	  	EXEC_ASSERT (globalSingletons. insert (fId). second);
}




// Change

Change::Change (const FeatureTree &tree_arg,
                istream &is)
: tree (tree_arg)
, from (nullptr)
, improvement (INF)
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
//const FeatureTree& tree = * from->getFeatureTree ();
  ASSERT (! tree. inputTreeFName. empty ());
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
  improvement = tree. len * (1 + tree. lenInflation) - tree. getLength ();
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
, to (nullptr)
{ 
  size_t toIndex; 
  is >> toIndex;
	to = static_cast <const Phyl*> (tree. nodes. at (toIndex)) -> asSpecies ();
}



void ChangeTo::saveText (ostream& os) const
{ 
  ASSERT (! tree /*to->getFeatureTree ().*/. inputTreeFName. empty ());
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
  auto inter_ = new Fossil (const_cast <FeatureTree&> (tree), arcEnd, string (), NAN);
  inter = inter_;
  inter_->init ();
  // PAR
  FOR (size_t, i, inter_->core. size ())
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
	lca = static_cast <const Phyl*> (Tree::getLowestCommonAncestor (oldParentRepr, inter)) -> asFossil ();
	ASSERT (lca);
  
  // Species::time
  const_cast <Species*> (from) -> assignTime (); 
  fromSibling->assignTime ();  
  if (isolatedTransientChild && oldParentRepr != fromSibling)
    oldParentRepr->assignTime (); 
  
  // Species::parent2core[]
  oldParentRepr->assignFeaturesUp (lca);  
  if (! oldParentRepr->descendentOf (fromSibling))
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
  if (! oldParentRepr->descendentOf (fromSibling))
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
  if (! oldParentRepr->descendentOf (fromSibling))
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
  
	lca = static_cast <const Phyl*> (Tree::getLowestCommonAncestor (oldParentRepr, to)) -> asFossil ();
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
  auto inter = new Fossil (const_cast <FeatureTree&> (tree), const_static_cast <Fossil*> (from->getParent ()), string (), NAN);
  inter->init ();
  FOR (size_t, i, inter->core. size ())
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
	const_cast <Species*> (from) -> isolateChildrenUp ();
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



// PAR
const Real FeatureTree::len_delta = 1e-2;  



#if 0
FeatureTree::FeatureTree (int root_species_id,
	                  const string &treeFName,
	                  Server &db)
: inputTreeFName (treeFName)
, genesExist (true)
, allTimeZero (false)
, timeOptimFrac (1)  
, lambda0 (NAN)
, time_init (NAN)
, genes (0)
, globalSingletonsSize (0)
, genomeGenes_ave (0)
, len (NAN)
, len_min (NAN)
, lenInflation (0)
, nodeIndex_max (0)
, coreSynced (false)
, savePhenChangesOnly (false)
{
	ASSERT (root_species_id > 0);
	
	
	if (treeFName. empty ())
	{
		{
  	  CQuery query (db.NewQuery());
  		query.SetSql ("select lambda0, time_init from PAN where id = (select T.pan from SPECIES T where T.id = @species)");
  	  query.SetParameter("@species", root_species_id);
  	  query.Execute();
  	  ITERATE(CQuery, row, query) 
  	  {
  	    ASSERT (isNan (lambda0));
  	    lambda0   = AsDoubleNan (row[1]);
  	    time_init = AsDoubleNan (row[2]);
		  }
		}
		loadPhylDb (db, root_species_id, 0);
  	}
  else
  	loadPhylFile (root_species_id, treeFName);
  ASSERT (root);
  ASSERT (nodes. front () == root);
	ASSERT (static_cast <const Species*> (root) -> id == (uint) root_species_id);	

  cerr << "Genomes ..." << endl;
  {
	  Progress prog;
	 	for (const DiGraph::Node* node : nodes)
	 		if (const Genome* g = static_cast <const Phyl*> (node) -> asGenome ())
	 		{
	 			prog (g->getName ());
	      const_cast <Genome*> (g) -> initDb (db);
	    }
	}

  abbreviate ();

  finish (& db, string ());  
}
#endif



FeatureTree::FeatureTree (const string &treeFName,
      						        const string &geneDir,
      						        const string &phenDir,
      						        const string &coreFeaturesFName,
      	                  bool genesExist_arg)
: inputTreeFName (treeFName)
, genesExist (genesExist_arg)
, allTimeZero (false)
, timeOptimFrac (1)
, lambda0 (NAN)
, time_init (NAN)
, genes (0)
, globalSingletonsSize (0)
, genomeGenes_ave (0)
, len (NAN)
, len_min (NAN)
, lenInflation (0)
, nodeIndex_max (0)
, coreSynced (false)
{
	ASSERT (! treeFName. empty ());
	IMPLY (genesExist, ! geneDir. empty ());

 	loadPhylFile (treeFName);
  ASSERT (root);
  ASSERT (nodes. front () == root);
  ASSERT (static_cast <const Phyl*> (root) -> asSpecies ());

  if (! genesExist)
    return;
    
  {
    cerr << "Genomes ..." << endl;
	  Progress prog;
	 	for (const DiGraph::Node* node : nodes)
	 		if (const Genome* g = static_cast <const Phyl*> (node) -> asGenome ())
	 		{
	 			prog (g->getName ());
	      const_cast <Genome*> (g) -> initDir (geneDir, phenDir);
	    }
	}

  finish (coreFeaturesFName);  
}



#if 0
void FeatureTree::loadPhylDb (Server &db,
	                         int species_id,
                           Fossil* parent)
{
	ASSERT (species_id > 0);

  Real time = NAN;	
	{
	  CQuery query (db.NewQuery());
		query.SetSql ("select [time] from SPECIES where id = @species_id");
	  query.SetParameter("@species_id", species_id);
	  query.Execute();
	  ITERATE(CQuery, row, query) 
	  {
	    ASSERT (isNan (time));
	    time = AsDoubleNan (row[1]);
	  }
	}

	bool isStrain = false;	
	{
	  CQuery query (db.NewQuery());
		query.SetSql ("select null from STRAIN where id = @species_id");
	  query.SetParameter("@species_id", species_id);
    isStrain = exists (query);;
	}
	
	if (isStrain)
	{
	  Vector<uint> ids;  
	  {
  	  CQuery query (db.NewQuery());
  		query.SetSql ("\
select id \
  from GENOME \
  where strain = @species_id");
  	  query.SetParameter("@species_id", species_id);
  	  query.Execute();
  	  ITERATE(CQuery, row, query) 
  	  	ids << (uint) row[1].AsInt4();
  	}
  	ASSERT (! ids. empty ());
  	if (ids. size () == 1)
  	{
      auto s = new Strain (*this, parent, toString (isNan (time) ? ids [0] : (uint) species_id), time);  
        // isNan(time) => there may be GENOME.id = species_id => duplicate Strain::getName()
	  	new Genome (*this, s, toString (ids [0]));
  	}
  	else
  	{
	    auto phyl = new Fossil (*this, parent, toString (species_id), time);
  	  for (const uint id : ids)
  	  {
    	  auto s = new Strain (*this, phyl, toString (id), isNan (phyl->time) ? NAN : 0.0);
  	  	new Genome (*this, s, toString (it));
  	  }
  	}
	}
	else
	{
    auto f = new Fossil (*this, parent, toString (species_id), time);
		Vector<int> subspecies;
	  {
		  CQuery query (db.NewQuery());
			query.SetSql ("select id from SPECIES where parent = @species_id");
		  query.SetParameter("@species_id", species_id);
		  query.Execute();
		  ITERATE(CQuery, row, query) 
		    subspecies << row[1].AsInt4();
	  }
	  ASSERT (! subspecies. empty ());
	  for (const int i : subspecies)
	  	loadPhylDb (db, i, f);
	}
}
#endif



namespace 
{
	Real str2time (const string &s)
	{
	  const string s1 (" " + s);
		const string tStr (" t=");
		const size_t pos = s1. find (tStr);
		if (pos == string::npos)
			return NAN;  
		string timeS (s1. substr (pos + tStr. size ()));
		return str2real (findSplit (timeS));
	}
}



bool FeatureTree::loadPhylLines (const Vector<string>& lines,
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
//ASSERT (root_species_id > 0);
	ASSERT (! treeFName. empty ());
	
		
  Vector<string> lines;
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



void FeatureTree::finish (const string &coreFeaturesFName)
{
  // optionalCore[i] in all Genome's => remove Feature i ??

  // globalSingletons, nonSingletons
  Set<Feature::Id> globalSingletons;
  Set<Feature::Id> nonSingletons;
 	for (const DiGraph::Node* node : nodes)
 		if (const Genome* g = static_cast <const Phyl*> (node) -> asGenome ())
 		  g->getSingletons (globalSingletons, nonSingletons);
  ASSERT (! globalSingletons. intersects (nonSingletons));
  globalSingletonsSize = globalSingletons. size ();
 	for (const DiGraph::Node* node : nodes)
 		if (Genome* g = const_cast <Genome*> (static_cast <const Phyl*> (node) -> asGenome ()))
 		{
      g->setSingletons (globalSingletons);
      ASSERT (g->coreNonSingletons >= g->singletons. size ());
      g->coreNonSingletons -= g->singletons. size ();
    }
  
  // commonCore
  ASSERT (commonCore. empty ());
  commonCore = nonSingletons;
 	for (const DiGraph::Node* node : nodes)
 		if (const Genome* g = static_cast <const Phyl*> (node) -> asGenome ())
 		  g->getCommonCore (commonCore);
  ASSERT (! globalSingletons. intersects (commonCore));
 	for (const DiGraph::Node* node : nodes)
 		if (Genome* g = const_cast <Genome*> (static_cast <const Phyl*> (node) -> asGenome ()))
      { EXEC_ASSERT (g->removeFromCoreSet (commonCore) == commonCore. size ()); }
      
  for (const Feature::Id& fId : commonCore)
    nonSingletons. erase (fId);    
  if (genesExist && nonSingletons. empty ())
  	throw runtime_error ("All genes are singletons or common core");
  	
  // features, feature2index
  ASSERT (features. empty ());
  map<Feature::Id, size_t/*index*/> feature2index;
  // Feature::isGene
#ifndef NDEBUG
  Feature::Id prevFeature;
#endif
  for (const Feature::Id& fId : nonSingletons)
  {
  	ASSERT (prevFeature < fId);
  	feature2index [fId] = features. size ();
  	features << Feature (fId, true);
  #ifndef NDEBUG
  	prevFeature = fId;
  #endif
  }
  genes = features. size ();
  ASSERT (genes == nonSingletons. size ());
  // !Feature::isGene
  {
    Set<Feature::Id> phens;
   	for (const DiGraph::Node* node : nodes)
   		if (const Genome* g = static_cast <const Phyl*> (node) -> asGenome ())
   	    for (const auto phenIt : g->phens)
   	      phens << phenIt. first; 
   	for (const Feature::Id& fId : phens)
   	{
  	  feature2index [fId] = features. size ();
   	  features << Feature (fId, false);
   	}
  }
  ASSERT (feature2index. size () == features. size ());

  // allTimeZero, Phyl::init()      
  size_t timeNan = 0;
  size_t timeNonNan = 0;
 	for (const DiGraph::Node* node : nodes)
 	{
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

  genomeGenes_ave = 0;  
  size_t genomes = 0;
 	for (const DiGraph::Node* node : nodes)
 		if (const Genome* g = static_cast <const Phyl*> (node) -> asGenome ())
 		{
 			genomeGenes_ave += g->getGenes ();
 			genomes++;
 	  }
  genomeGenes_ave /= genomes;

  if (! allTimeZero && genesExist)
    loadRootCoreFile (coreFeaturesFName);
 	for (DiGraph::Node* node : nodes)  
 		static_cast <Phyl*> (node) -> setWeight ();
  setLenGlobal ();
  setCore ();
  len_min = getLength_min ();

  ASSERT (nodes. front () == root);
#ifndef NDEBUG
  if (! inputTreeFName. empty ())
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



void FeatureTree::qc () const
{ 
  if (! qc_on)
    return;
	Tree::qc ();
		
	ASSERT (root);
	ASSERT (root->graph == this);
	const Species* root_ = static_cast <const Phyl*> (root) -> asSpecies ();
	ASSERT (root_);
//ASSERT (nodes. size () <= numeric_limits<unsigned short>::max ());  // To fit Phyl::CoreEval::treeLen
	if (! allTimeZero && genesExist)
  {
    ASSERT (root_->time == getRootTime ());
    if (! emptyRoot)
    	FOR (unsigned char, i, 2)
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

  if (allTimeZero || ! genesExist)
  { 
    ASSERT (timeOptimWhole ());
    ASSERT (isNan (lambda0)); 
    ASSERT (isNan (time_init)); 
  }
  else
  {
    ASSERT (rootCore. size () == features. size ());
    ASSERT (isProb (timeOptimFrac));
	  ASSERT (isProb (lambda0));
	  ASSERT (lambda0 < 1 /*0.5*/);
    ASSERT (time_init > 0);
	}

	Feature prevFeature;
	FOR (size_t, i, features. size ())
	{
	  const Feature f (features [i]);
	  f. qc ();
		IMPLY (prevFeature. isGene == f. isGene, prevFeature. name < f. name);
		ASSERT (prevFeature. isGene >= f. isGene);
		ASSERT (f. isGene == (i < genes));
		prevFeature = f;
	}

  ASSERT ((bool) genes == genesExist);	
	ASSERT (genes <= features. size ());
	
	ASSERT (geReal (len, len_min));
	ASSERT (len_min >= 0);
	
	IMPLY (genesExist, genomeGenes_ave > 0);
	ASSERT (genomeGenes_ave <= getTotalGenes ());
	
	ASSERT ((bool) distDistr. get () == (bool) depthDistr. get ());
	
	if (coreSynced) 
	{
		if (! emptyRoot)
			FOR (size_t, i, root_->core. size ())
			{
			  ASSERT (root_->core [i] == root_->feature2parentCore (i));
			  IMPLY (! allTimeZero && genesExist, root_->core [i] == rootCore [i]);
			}
		const Real subtreeLen = static_cast <const Phyl*> (root) -> getSubtreeLength ();
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
	
//IMPLY (! targetFeatures. empty (), savePhenChangesOnly);	
}



void FeatureTree::saveText (ostream& os) const
{ 
  ASSERT (coreSynced);
  
  if (! allTimeZero && genesExist)
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
  
  const_cast <Genome*> (leaf) -> isolateChildrenUp ();
  
	const_cast <Strain*> (parent) -> isolateChildrenUp ();
	delayDelete (const_cast <Strain*> (parent));

  if (deleteTransientAncestor && grandparent->isTransient ())
  {
  	const_cast <Fossil*> (grandparent) -> isolateChildrenUp ();
  	delayDelete (const_cast <Fossil*> (grandparent));
  }
	
  delete leaf;
  toDelete. deleteData ();
}



void FeatureTree::printInput (ostream& os) const
{
//os << "SPECIES.id: " << static_cast <const Species*> (root) -> id << endl;

  if (inputTreeFName. empty ())
	  os << "Tree from database" << endl;
	else
	  os << "Tree from file: " << inputTreeFName << endl;
	  
  const size_t genomes = root->getLeavesSize ();
  os << "# Genomes: " << genomes << endl;
  os << "# Species: " << nodes. size () - genomes << endl;
  os << "# Common core genes: " << commonCore. size () << endl;
  os << "# Singleton genes:   " << globalSingletonsSize << endl;
  os << "# Other genes:       " << genes << endl;
  os << "# Phenotypes:        " << features. size () - genes << endl;
  os << "Genome size ave.:    " << genomeGenes_ave << endl;

  os << "Time: " << (allTimeZero ? "Not used" : "Used") << endl;
  if (! allTimeZero)
  {
  	ONumber oNum (os, 3, true);  // PAR
    os << "Lambda_0          = " << lambda0 << endl;
    os << "Initial time      = " << time_init << endl;
    os << "timeOptimFrac     = " << timeOptimFrac << endl;
    if (! rootCore. empty ())
    {
      size_t n = commonCore. size ();
      FOR (size_t, i, genes)
        if (rootCore [i])
          n++;
      os << "# Root core genes = " << n << endl;
    }
  }
  
 	ONumber oNum (os, 0, false);  // PAR
  os << "Tree length min.  = " << len_min << endl;
  os << "Tree length       = " << len << endl;
  os << endl;
}



void FeatureTree::dump (const string &fName,
                        bool setIds)
{ 
 	ONumber oNum (cout, 0, false);  // PAR
  cout << "Tree length = " << len << endl;  
  cout. flush ();
  setCore ();
  sort ();
  setStats ();
  if (setIds)
    if (const Fossil* f = static_cast <const Species*> (root) -> asFossil ())
    { 
      uint id = 1;  // 0 means "unknown"
      const_cast <Fossil*> (f) -> setId (id);
    }
  qc ();      
  saveFile (fName);
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
	Real x = NAN;
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
  	
	if (! allTimeZero)
	{
		{
		  ONumber oNum (cout, 3, true);  // PAR
		  cout << "lambda0 = " << lambda0;
		}
		{
		  ONumber oNum (cout, 3, false);  // PAR
	    cout << "  len = " << len
	         << endl;
	  }
  }
}
#endif



void FeatureTree::useTime (const string &coreFeaturesFName)
{ 
	ASSERT (allTimeZero);
	ASSERT (genesExist);

	allTimeZero = false;
  loadRootCoreFile (coreFeaturesFName);
  setCore ();
  timeOptimFrac = 0;


  size_t parent2core [2/*thisCore*/] [2/*parentCore*/];
  getParent2core_sum (parent2core);
	if (verbose ())
    FOR (unsigned char, i, 2)
      FOR (unsigned char, j, 2)
        cout << "parent2core[" << (int) i << "][" << (int) j << "]=" << parent2core [i] [j] << endl;

  Real len_old;  
  Real len1 = INF;;
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
		if (! leReal (len1, len_old + len_delta)) 
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


  // Cf. finish()
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

  FOR (unsigned char, i, 2)
    FOR (unsigned char, j, 2)
      parent2core [i] [j] = 0;
 	for (const DiGraph::Node* node : nodes)  
 		if (const Species* s = static_cast <const Phyl*> (node) -> asSpecies ())
 		{
   		if (! s->getParent ())
   			continue;
  
      size_t parent2core_ [2/*thisCore*/] [2/*parentCore*/];
      s->getParent2corePooled (parent2core_);	  
  		FOR (size_t, i, genes)
  		  parent2core_ [s->core [i]] [s->feature2parentCore (i)] ++;
  
      FOR (unsigned char, i, 2)
        FOR (unsigned char, j, 2)
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
   	  FOR (unsigned char, i, 2)
   	    FOR (unsigned char, j, 2)
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
      TRY_CHANGE (ChangeToCousin,  (from, to));  
   	//TRY_CHANGE (ChangeToParent,  (from, to));  
  	  if (verbose ())
  	  {
    	  ASSERT (to->graph);
  	    to->qc ();
  	  }
    }
  
//TRY_CHANGE (ChangeRoot, (from));
  TRY_CHANGE (ChangeDel,  (from));  

#if 0  
  // From DistTree:
  if (   ! bestChange 
      && area. size () < 100  // PAR
     )
  {
    if (verbose (1))
      cerr << "more ";
   	for (const Tree::TreeNode* node : area)  
   	{
   		DTNode* to = const_static_cast <DTNode*> (node);
 	    // Applied rarely
 	    TRY_CHANGE (ChangeToChild,  (from, to));  
   	  TRY_CHANGE (ChangeToUncle,  (from, to));  
      TRY_CHANGE (ChangeToCousin, (from, to));  
  	  if (verbose ())
  	  {
    	  ASSERT (to->graph);
  	    to->qc ();
  	  }
    }
  }
#endif

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
	
	
//cout << endl;	    

  const Real len_init = len;


  Common_sp::sort (changes, Change::compare);	
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

	  cout << "Apply: ";
	  ch->print (cout);
		ch->qc ();  
	#ifndef NDEBUG
	  const bool first = ch == changes. front ();  
	#endif
	  changedNodes << ch->getCoreChanged ();
  	Real len_old = len;
    ch->apply ();
	  len = getLength ();
	  ONumber on (cout, 0, false);
	  cout << "len = " << len << endl;
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
	    cout << "restore" << endl;
		  cout << "len = " << len << endl;
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
  

  const Real improvement = max (0.0, len_init - len);
  cout << "Improvement = " << improvement << endl;
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
	  	cout << "timeOptimFrac = " << timeOptimFrac << endl;
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
		Progress prog ((uint) nodeVec. size ());
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



void FeatureTree::findRoot () 
{
  ASSERT (allTimeZero);

  setCore ();

  const Species* bestFrom = nullptr;
 	size_t bestCoreSize = 0;
  FOR (size_t, i, genes)
    if (static_cast <const Species*> (root) -> core [i])
      bestCoreSize++;
 	for (const DiGraph::Node* node : nodes)
 		if (const Species* from = static_cast <const Phyl*> (node) -> asSpecies ())
 		{
   		const Species* parent = static_cast <const Species*> (from->getParent ());
   		if (! parent)
   		  continue;
   	 	ASSERT (from != root);
      size_t coreSize = 0;
      FOR (size_t, i, genes)
        if (   from  ->core [i] 
            && parent->core [i]
           )
        coreSize++;
      if (minimize (bestCoreSize, coreSize))
        bestFrom = from;
   	}

  if (bestFrom)
  {
    {
      ChangeRoot ch (bestFrom);    
      ch. apply ();
    	ch. commit ();
    }
    ASSERT (eqTreeLen (len, getLength ()));
    clearStats ();
    finishChanges ();
  }    
  
  // rootCore
  ASSERT (rootCore. empty ());
  rootCore. resize (features. size (), false);
  FOR (size_t, i, features. size ())
    if (getRootCore (i))
      rootCore [i] = true;
      
  setCore ();
}



void FeatureTree::resetRootCore (size_t coreChange [2/*core2nonCore*/])
{
  ASSERT (! allTimeZero);
  ASSERT (coreSynced);
  
 	FOR (unsigned char, b, 2)
 	  coreChange [b] = 0;

	const Species* root_ = static_cast <const Species*> (root);
	if (root_->arcs [false]. size () < 2)
	  return;
	  
  FOR (size_t, i, rootCore. size ())
  {
    size_t n = 0;
  	for (const DiGraph::Arc* arc : root_->arcs [false])
    {
      const Phyl* s = static_cast <Phyl*> (arc->node [false]);
  	  if (s->core [i])
  	    n++;
  	}
  	if (n <= 1)
  	  if (rootCore [i])
  	  {
  	    rootCore [i] = false;
  	    if (i < genes)
  	      coreChange [true] ++;
  	  }
  	if (n == root_->arcs [false]. size ())
  	  if (! rootCore [i])
  	  {
  	    rootCore [i] = true;
  	    if (i < genes)
  	      coreChange [false] ++;
  	  }
  }

  len = getLength ();
  setCore ();
  clearStats ();
}



#if 0
void FeatureTree::loadRootCoreDb (Server* db)
{ 
  ASSERT (db);
	ASSERT (! allTimeZero);
	ASSERT (rootCore. empty ());
	    
  typedef  map<Feature::Id, size_t>  Feature2index;
  Feature2index feature2index;
  FOR (size_t, i, features. size ())
    feature2index [features [i]. name] = i;
  ASSERT (feature2index. size () == features. size ());

  rootCore. resize (features. size (), false);
	size_t miss = 0;
  CQuery query (db->NewQuery());
  query.SetSql("\
select prot \
  from SPECIES_PROT \
  where     species = @species \
        and isnull(gain,1) = 1");
	query.SetParameter("@species", (int) static_cast <const Species*> (root) -> id);
  query.Execute();
  ITERATE(CQuery, row, query) 
	{
	  Feature2index::const_iterator it = feature2index. find (toString (row [1].AsInt4()));
	  if (it == feature2index. end ())
	    miss++;
	  else
	    rootCore [it->second] = true;
	}
	ASSERT (miss == commonCore. size ());
}
#endif



void FeatureTree::loadRootCoreFile (const string &coreFeaturesFName)
{ 
	ASSERT (! coreFeaturesFName. empty ());
	ASSERT (! allTimeZero);
	ASSERT (rootCore. empty ());
	  
  typedef  map<Feature::Id, size_t>  Feature2index;
  Feature2index feature2index;
  FOR (size_t, i, features. size ())
    feature2index [features [i]. name] = i;
  ASSERT (feature2index. size () == features. size ());

  rootCore. resize (features. size (), false);
	LineInput f (coreFeaturesFName, 10 * 1024);  // PAR
	size_t miss = 0;
	while (f. nextLine ())
	{
	  trim (f. line);
	  Feature2index::const_iterator it = feature2index. find (f. line);
	  if (it == feature2index. end ())
	    miss++;
	  else
	    rootCore [it->second] = true;
	}
	ASSERT (miss == commonCore. size ());
}



void FeatureTree::saveRootCore (const string &coreFeaturesFName) const
{
  ASSERT (! coreFeaturesFName. empty ());
  ASSERT (! rootCore. empty ());

 	Set<Feature::Id> s (commonCore);
  FOR (size_t, i, features. size ())
    if (rootCore [i])
      s << features [i]. name;
      
  OFStream f ("", coreFeaturesFName, "");
  for (const Feature::Id& fId : s)
    f << fId << endl;
}



void FeatureTree::setStats () 
{
  ASSERT (coreSynced);
  
  for (Feature& f : features)
  {
    f. genomes = 0;
    f. gains = 0;
    f. losses = 0;
  }
  
  Dataset ds;
  auto dist  = new RealAttr1 ("Distance", ds);  
  auto depth = new RealAttr1 ("Depth",    ds);  
 	for (const DiGraph::Node* node : nodes)
 	{
 		const Phyl* phyl = static_cast <const Phyl*> (node);
  	const size_t n = ds. appendObj (phyl->getName ());
    (*dist)  [n] = phyl->getNeighborDistance ();
    (*depth) [n] = (Real) phyl->getDepth ();
    const Genome* g = phyl->asGenome ();
    FOR (size_t, i, features. size ())
    {
      if (g)
        if (! g->optionalCore [i])
          if (g->core [i])
            features [i]. genomes++;
      if ((! phyl->feature2parentCore (i) || phyl == root) && phyl->core [i])  // use unrooted tree ??
        features [i]. gains++;
      if (phyl->feature2parentCore (i) && ! phyl->core [i])
        features [i]. losses++;
    }
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
  
#if 0
  const_static_cast <Species*> (root) -> setHasPhenChange ();
#endif
}
  


void FeatureTree::clearStats () 
{
  for (Feature& f : features)
  {
    f. genomes = 0;
    f. gains = 0;
    f. losses = 0;
  }
  
  distDistr. reset (nullptr);
  depthDistr. reset (nullptr);

#if 0  
 	for (DiGraph::Node* node : nodes)
 		static_cast <Phyl*> (node) -> hasPhenChange = false;
#endif
}
  


#if 0
size_t FeatureTree::badNodesToRoot ()  
{
  setDistr ();
  
  VectorPtr<Phyl> badNodes;  
  const_static_cast <Species*> (root) -> getBadNodes (badNodes, false);
  
  for (const Phyl* phyl : badNodes)
  {
    ASSERT (phyl->getParent ());
    if (phyl->getParent () == root)
      continue;
    ChangeToParent ch (const_cast <Phyl*> (phyl), const_static_cast <Species*> (root));
    ch. qc ();
    ch. apply ();
	  len = getLength ();
    ch. commit ();
  }
  
  optimizeTime ();

	VectorOwn<Change> changes;  changes. reserve (badNodes. size ());  
  {
		Progress prog ((uint) badNodes. size ());
	 	for (const Phyl* phyl : badNodes)  
	 	  if (phyl->graph)  
  	 	{
  	   	prog ();
  	  //cout << endl;  
  		 	if (const Change* bestChange = getBestChange (phyl))
  		  { 
  		  	ASSERT (positive (bestChange->improvement));
  		  	changes << bestChange;
  		  //bestChange->print (cout); 
  		  }
  		//else
  		  //cout << (*it)->getName () << ": no change" << endl;  
  	  }
  	}

  finishChanges ();
  setCore ();
  
  qc ();
  
  return badNodes. size ();
}
#endif



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
 				const_cast <Fossil*> (f) -> isolateChildrenUp ();
 				delayDelete (const_cast <Fossil*> (f));
 				n++;
 			}
 			
  setLenGlobal ();
  
 	if (verbose ())
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
	  for (const auto it : freq)
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



#if 0
void FeatureTree::loadTargetFeatures (const string &fName)
{
  ASSERT (targetFeatures. empty ());
  
  if (fName. empty ())
    return;
  
  ObjectInput in (fName, 1024);  // PAR
  TargetFeature tf;
  while (in. next (tf))
  {
    tf. index = features. binSearch (Feature (tf. featureId, true));
    FOR_REV (size_t, i, targetFeatures. size ())
      if (targetFeatures [i]. name == tf. name)
      {
        tf. serial = targetFeatures [i]. serial + 1;
        break;
      }
    targetFeatures << tf;
  }

  if (! targetFeatures. empty ())
    savePhenChangesOnly = true;
}
#endif



#if 0
void FeatureTree::saveDatabase (Server &db,
                             int root_species_id)
{
	ASSERT (root_species_id > 0);
	ASSERT (coreSynced);
	ASSERT (! allTimeZero);
	ASSERT (timeOptimWhole ());


  {	
   	Transaction tr (db);
  
    {
    	cerr << "SPECIES,STRAIN ..." << endl;
  		bool parent_STRAIN = true;  // Ignored for root
  		const_static_cast <Species*> (root) -> id = (uint) root_species_id;
  		Progress::Start progStart (prog_, (uint) nodes. size ());
  	  const_static_cast <Species*> (root) -> saveDatabaseTopology (db, parent_STRAIN);
  	    // Genome::saveDatabaseTopology() invokes SPECIES_redo() which delete's from SPECIES_PROT
  	}
  
    {
  	  CQuery query (db.NewQuery());
    	query.SetParameter("@id", root_species_id);
  	  query.ExecuteSP("SPECIES_redo");
  	}
  }


  {	
   	Transaction tr (db);
  
    {
      CQuery query (db.NewQuery());
      query.SetSql("\
  update PAN \
    set dat = sysdatetime() \
      , lambda0 = @lambda0 \
      , time_init = @time_init \
      , genes_ave = @genes_ave \
    where id = (select T.pan from SPECIES T where T.id = @id)");
    	query.SetParameter("@id", (int) root_species_id);
    	query.SetParameter("@lambda0", lambda0);
    	query.SetParameter("@time_init", time_init);
    	query.SetParameter("@genes_ave", (int) genomeGenes_ave);
    	query.Execute();
  	  ASSERT (query.GetRowCount() == 1);
    }
   	
    {
    	cerr << "{SPECIES|STRAIN}_PROT ..." << endl;
  		Progress::Start progStart (prog_, (uint) (nodes. size ()));
  	  const_static_cast <Species*> (root) -> saveDatabaseNode (db);
  	}
  
    {
    	cerr << "SPECIES_finish ..." << endl;
  	  CQuery query (db.NewQuery());
  	  query.ExecuteSP("SPECIES_finish");
  	    // Deletes the old topology
  	}
  }

  	
  {
  	cerr << "SPECIES_length ..." << endl;
	  CQuery query (db.NewQuery());
  	query.SetParameter("@root_species", root_species_id);
	  query.ExecuteSP("SPECIES_length");
	}

    
  saveDatabasePhen (db);
}



void FeatureTree::saveDatabasePhen (Server &db) 
{
	ASSERT (coreSynced);
	ASSERT (timeOptimWhole ());

 	Transaction tr (db);

  {
  	cerr << "{SPECIES|STRAIN}_PHEN ..." << endl;
  	Progress::Start progStart (prog_, (uint) (nodes. size ()));
    const_static_cast <Species*> (root) -> saveDatabasePhen (db);
  }

  {
  	cerr << "SPECIES_PHEN_length ..." << endl;
	  CQuery query (db.NewQuery());
	  query.ExecuteSP("SPECIES_PHEN_length");
	}
}
#endif



}



/*
TODO: ??

MLE should sum over internal states

Borrow Change's from distTree
Local Change's involving more than 1 node
  
*/
