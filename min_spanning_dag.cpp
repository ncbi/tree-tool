// min_spanning_dag.cpp  

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
*   Minimum spanning DAG
*
*/


#undef NDEBUG

#include "common.hpp"
#include "graph.hpp"
using namespace Common_sp;
#include "version.inc"

#include "common.inc"



namespace 
{
  
  
struct Item final : DiGraph::Node
{
  const string name;


  Item (DiGraph &graph_arg,
        const string &name_arg)
    : DiGraph::Node (graph_arg)
    , name (name_arg)
    {}
};



struct Link final : DiGraph::Arc
{
  const double weight;
  
  
  Link (Item* start,
        Item* end,
        double weight_arg)
    : DiGraph::Arc (start, end)
    , weight (weight_arg)
    {
      ASSERT (! isNan (weight));
    }
  void saveText (ostream &os) const final
    { os         << static_cast <const Item*> (node [false]) -> name
         << '\t' << weight
         << '\t' << static_cast <const Item*> (node [true]) -> name
         << '\n';
    }
    
    
  bool operator< (const Link &other) const
    { return weight > other. weight; }
};



struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Print minimum spanning DAG")
    {
      version = VERSION;
  	  addPositional ("in", "File with lines: <item1> <item2> <weight>. Sorted by <weight> ascending, or descending if -max");
  	  addFlag ("max", "<weight> must be maximized in DAG, otherwise minimized");
  	  addKey ("rename", "File with pairs: <item> <tab> <new name>");
  	}


	void body () const final
	{
		const string inFName    = getArg ("in");
		const bool   maxP       = getFlag ("max");
		const string matchFName = getArg ("rename");
		
		
		KeyValue old2new;
		if (! matchFName. empty ())
		{
      LineInput fIn (matchFName);  
      Set<string> newNames;
      while (fIn. nextLine ())
      {
        string& newName = fIn. line;
        const string oldName (findSplit (newName, '\t'));
        QC_ASSERT (! newName. empty ());
        if (contains (old2new, oldName))
          throw runtime_error ("Duplicate old name: " + strQuote (oldName));
        if (newNames. contains (newName))
          throw runtime_error ("Duplicate new name: " + strQuote (newName));
        old2new [oldName] = newName;
      }
		}


    unordered_map<string/*Item::name*/,Item*> items;  items. rehash (100000);  // PAR
    DiGraph gr;  // Of Item, DAG
    {
      double weight_prev = - numeric_limits<double>::max ();
      if (maxP)
        weight_prev *= -1;
      LineInput fIn (inFName);  
      Istringstream iss;
      while (fIn. nextLine ())
      {
        string childName;
        string parentName;
        double weight = NaN;
        iss. reset (fIn. line);
        iss >> childName >> parentName >> weight;
        QC_ASSERT (! parentName. empty ());
        QC_ASSERT (! isNan (weight));        
        if (     (maxP && weight > weight_prev)
            || (! maxP && weight < weight_prev)
           )
          throw runtime_error ("Weight disorder: " + fIn. lineStr ());
        trim (childName);
        trim (parentName);
        if (childName == parentName)
          continue;
          
        childName  = find (old2new, childName,  true);
        parentName = find (old2new, parentName, true);
          
        // child, parent
        if (! contains (items, childName))
          items [childName] = new Item (gr, childName);
        if (! contains (items, parentName))
          items [parentName] = new Item (gr, parentName);
        Item* child = items [childName];
        Item* parent = items [parentName];
        ASSERT (child);
        ASSERT (parent);
        ASSERT (child != parent);
        
        gr. clearReachable ();
        child->setReachable (true);
        if (parent->reachable)
          continue;
        
        gr. clearReachable ();
        parent->setReachable (true);
        if (child->reachable)
          continue;

        new Link (child, parent, weight);
        child->DisjointCluster::merge (* parent);

        weight_prev = weight;
      }
    }
    ASSERT (gr. nodes. size () == items. size ());
    gr. qc ();


    map<const DisjointCluster*,VectorPtr<Item>> clusters;
    for (auto& it : items)
      clusters [it. second->getDisjointCluster ()] << it. second;


    bool first = false;
    for (const auto& it : clusters)
    {
      const VectorPtr<Item>& cluster = it. second;
      ASSERT (! cluster . empty ());
      VectorPtr<Link> links;  links. reserve (cluster. size () * 10);  // PAR      
      for (const Item* item : cluster)
        for (const DiGraph::Arc* arc : item->arcs [true])
        {
          ASSERT (arc);
          links << static_cast <const Link*> (arc);
        }
      links. sortPtr ();
      if (first)
        first = false;
      else
        cout << endl;
      for (const Link* link : links)
        link->saveText (cout);
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
