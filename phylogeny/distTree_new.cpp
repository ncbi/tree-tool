// distTree_new.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "distTree.hpp"
using namespace DistTree_sp;



namespace 
{



#if 0
typedef  Vector<Real>  Leaf2dist;  
  // Index: Leaf::index ??



Steiner* DTNode::Closest::insert () 
// Return: new
// Update: *node
{
  auto st = new Steiner ( const_cast <DistTree&> (node->getDistTree ())
                        , const_static_cast <Steiner*> (node->getParent ())
                        , node->len - arcLen
                        );
  const_cast <DTNode*> (node) -> setParent (st);
  const_cast <DTNode*> (node) -> len = arcLen;
  st->subtreeLeaves = node->subtreeLeaves;
  return st;
}




void DTNode::findClosestNode (const Leaf2dist &leaf2dist,
                              Leaf2dist &leaf2hat_dist,
                              Closest &closest) const
// Input: subtreeLeaves
// Update: lead2hat_dist: matches *getParent()
//         closest
// Time: O(n^2)
{
  ASSERT (! leaf2dist. empty ());
  ASSERT (leaf2dist. size () == subtreeLeaves. size ());
  ASSERT (leaf2hat_dist. empty () || leaf2hat_dist. size () == leaf2dist.size ());


  bool foundSubLeaf   = false;
  bool foundSuperLeaf = false;
	FOR (size_t, i, leaf2dist. size ())
	  if (! isNan (leaf2dist [i]))
	  {
	    if (subtreeLeaves [i])
	      foundSubLeaf = true;
	    if (! subtreeLeaves [i])
	      foundSuperLeaf = true;
	  }
	if (! foundSubLeaf)
	  return;


  // leaf2hat_dist
  // Depth-first order
  if (leaf2hat_dist. empty ())
  {
    if (arcs [false]. empty ())  // can it be done at root ??
    {
      leaf2hat_dist. resize (leaf2dist. size (), NAN);
      const Leaf* leaf = asLeaf();
      ASSERT (leaf);
      ASSERT (leaf->index != NO_INDEX);
      ASSERT (! isNan (leaf2dist [leaf->index]));
      for (const DiGraph::Node* node : getTree (). nodes)
        if (const Leaf* g = static_cast <const DTNode*> (node) -> asLeaf ())
          if (   g->index != NO_INDEX
              && ! isNan (leaf2dist [g->index])
             )
            leaf2hat_dist [g->index] = DistTree::path2prediction (DistTree::getPath (leaf, g));
    }
  }
  else
    // *getParent() --> *this
  	FOR (size_t, i, leaf2dist. size ())
  	{
  	  const Real u = subtreeLeaves [i] ? 1 : -1;
  	  leaf2hat_dist [i] -= len * u;
  	}

#ifndef NDEBUG
  Leaf2dist leaf2hat_dist_ (leaf2hat_dist);
#endif
  for (const DiGraph::Arc* arc : arcs [false])
  {
    const DTNode* child = static_cast <const DTNode*> (arc->node [false]);
    if (const Leaf* g = child->asLeaf ())
      if (g->index == NO_INDEX)
        continue;
    const_cast <DTNode*> (child) -> findClosestNode (leaf2dist, leaf2hat_dist, closest);
  #ifndef NDEBUG
    if (leaf2hat_dist_. empty ())
      leaf2hat_dist_ = leaf2hat_dist;
    else
    {
      ASSERT (leaf2hat_dist. size () == leaf2hat_dist_. size ());
    	FOR (size_t, i, leaf2hat_dist. size ())
      {
        ASSERT (isNan (leaf2hat_dist [i]) == isNan (leaf2hat_dist_ [i]));
        if (! isNan (leaf2hat_dist [i]))
          { ASSERT_EQ (leaf2hat_dist [i], leaf2hat_dist_ [i], 1e-3); } // PAR
      }
    }
  #endif
  }
  ASSERT (! leaf2hat_dist. empty ());
	  

	if (this == getTree (). root)
	  return;
	  

	if (len > 0 && foundSuperLeaf)
	{
  	Real mult_sum = 0;
  	Real u_avg = 0;	
  	Real delta_avg = 0;
  	FOR (size_t, i, leaf2dist. size ())
  	{
  	  const Real dist = leaf2dist [i];
	    if (isNan (dist))
	      continue;
  	  ASSERT (positive (dist));
  	  const Real mult = dissim2mult (dist);
  	  mult_sum += mult;
  	  const Real u = subtreeLeaves [i] ? 1 : -1;
  	  const Real delta = dist - leaf2hat_dist [i];
  	  u_avg     += mult * u;
  	  delta_avg += mult * delta;
  	}
  	ASSERT (positive (mult_sum));
  	u_avg     /= mult_sum;
  	delta_avg /= mult_sum;
  	
  	Real u_var = 0;
  	Real delta_u_cov = 0;
  	FOR (size_t, i, leaf2dist. size ())
  	{
  	  const Real dist = leaf2dist [i];
	    if (isNan (dist))
	      continue;
  	  const Real mult = dissim2mult (dist);
  	  const Real u = subtreeLeaves [i] ? 1 : -1;
  	  const Real delta = dist - leaf2hat_dist [i];
  	  delta_u_cov += mult * (delta - delta_avg) * (u - u_avg);
  	  u_var       += mult * sqr (u - u_avg);
  	}
  	u_var       /= mult_sum;
  	delta_u_cov /= mult_sum;  	
  	ASSERT (u_var > 0);

  	const Real arcLen = min (max (0.0, delta_u_cov / u_var), len);
  	const Real leafLen = max (0.0, delta_avg - u_avg * arcLen);
  	
  	Real absCriterion_delta = 0;
  	FOR (size_t, i, leaf2dist. size ())
  	{
  	  const Real dist = leaf2dist [i];
	    if (isNan (dist))
	      continue;
  	  const Real mult = dissim2mult (dist);
  	  const Real u = subtreeLeaves [i] ? 1 : -1;
  	  const Real delta = dist - leaf2hat_dist [i];
  	  absCriterion_delta += mult * sqr (delta - arcLen * u - leafLen);
  	}
  	ASSERT (absCriterion_delta >= 0);
  	
  	if (absCriterion_delta < closest. absCriterion_delta)
  	{
  	  closest = Closest (this, absCriterion_delta, leafLen, arcLen);
  	  closest. qc ();
  	}
  }


  // *this --> *getParent()
	FOR (size_t, i, leaf2dist. size ())
	{
	  const Real u = subtreeLeaves [i] ? 1 : -1;
	  leaf2hat_dist [i] += len * u;
	}
}
#endif




struct ThisApplication : Application
{
	ThisApplication ()
	: Application ("Find location of new objects in a distance tree")
	{
	  addPositional ("data", "Directory with data");
	  addFlag ("init", "Initialize search");
	  addKey ("variance", "Dissimilarity variance: " + varianceTypeNames. toString (" | "), varianceTypeNames [varianceType]);
	}
	
	
	
	void body () const
  {
	  const string dataDir      = getArg ("data");
	  const bool init           = getFlag ("init");
	               varianceType = str2varianceType (getArg ("variance"));  // Global    
    if (! isRight (dataDir, "/"))
      throw runtime_error ("\"" + dataDir + "\" must end with '/'");


    DistTree::printParam (cout);
    cout << endl;

    DistTree tree (dataDir, false);
    tree. setReprLeaves ();  
    if (verbose ())
      tree. qc ();     

    tree. printInput (cout);
    cout << endl;
    
    cout << "Processing new objects ..." << endl;
    const string newDir (dataDir + "search/");
    FileItemGenerator fig (1, true, newDir);  // PAR
	  string item;
	  while (fig. next (item))
    {
      const NewLeaf nl (tree, newDir, item, init);
      nl. qc ();
    }
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


