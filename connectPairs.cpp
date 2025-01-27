// connectPairs.cpp

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
*   Create sets of connected items
*
*/


#undef NDEBUG

#include "common.hpp"
using namespace Common_sp;
#include "version.inc"

#include "common.inc"


namespace 
{
  
  
static const string ext ("component");



struct Item : Named, DisjointCluster
{
  size_t n {1};
  
  
  explicit Item (const string &name_arg)
    : Named (name_arg)
    {}


	static bool lessPtr (const Item* x,
	                     const Item* y)
	  { ASSERT (x);
	    ASSERT (y);
	    LESS_PART (*y, *x, n);
	    return x->name < y->name; 
	  }
};

  
  

struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Create disjoint sets of items")
    {
      version = VERSION;
  	  addPositional ("in", "List of pairs of connected items: <item1> <item2>\nwhere <item1>, <item2> are strings with no spaces");
  	  addPositional ("out", "If -pairs then output file with pairs <item> <tab> <item_min>, otherwise output directory with the sets of connected items, where each set is a file named by its lexicographycally smallest item with the added extension " + strQuote ("." + ext));
  	  addKey ("subset", "List of items. Item pairs are restricted to this list");
  	  addFlag ("pairs", "<out> is a list of pairs: <item> <item_min>, where <item_min> is the lexicographycally smallest item of the cluster");
  	  addFlag ("center", "<item_min> is the item with the maximum number of connected items");
  	}


	void body () const final
	{
		const string inFName     = getArg ("in");
		const string outFName    = getArg ("out");
		const string subsetFName = getArg ("subset");
		const bool   pairsP      = getFlag ("pairs");
		const bool   centerP     = getFlag ("center");


    unique_ptr<StringVector> subset;
    if (! subsetFName. empty ())
    {
      subset. reset (new StringVector (subsetFName, (size_t) 10000, true));  // PAR 
      subset->sort ();
      QC_ASSERT (subset->isUniq ());
    }
      
    unordered_map<string/*Item::name*/,Item*> items;  items. rehash (100000);  // PAR
      // Not delete'd
    {
      LineInput fIn (inFName);  
      string s1, s2;
      Istringstream iss;
      while (fIn. nextLine ())
      {
        iss. reset (fIn. line);
        s2. clear ();
        iss >> s1 >> s2;
        QC_ASSERT (! s2. empty ());
        trim (s1);
        trim (s2);
        if (subset)
        {
          if (! subset->containsFast (s1))
            continue;
          if (! subset->containsFast (s2))
            continue;
        }
        Item*& it1 = items [s1];
        Item*& it2 = items [s2];
        if (it1) it1->n++; else it1 = new Item (s1);
        if (it2) it2->n++; else it2 = new Item (s2);
        ASSERT (it1);
        ASSERT (it2);
        it1->DisjointCluster::merge (*it2 );
      }
    }
    
    map<const DisjointCluster*,VectorPtr<Item>> clusters;
    for (auto& it : items)
      clusters [it. second->getDisjointCluster ()] << it. second;
      
    unique_ptr<OFStream> fOut;
    if (pairsP)
      fOut. reset (new OFStream (outFName));
    for (auto& it : clusters)
    {
      VectorPtr<Item>& cluster = it. second;
      if (centerP)
        cluster. sort (Item::lessPtr);
      else
        cluster. sort (Named::lessPtr);  
      ASSERT (! cluster . empty ());
      const string repr (cluster [0] -> name);
      if (! pairsP)  
        fOut. reset (new OFStream (outFName, repr, ext));
      ASSERT (fOut);
      for (const Item* item : cluster)
      {
        *fOut << item->name;
        if (pairsP)
          *fOut << '\t' << repr;
        *fOut  << '\n';
      }
      if (! pairsP)
        fOut. reset ();
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
