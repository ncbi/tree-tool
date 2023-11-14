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
  explicit Item (const string &name_arg)
    : Named (name_arg)
    {}
};

  
  

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Create sets of connected items")
  {
    version = VERSION;
	  addPositional ("in", "List of pairs of connected items: <item1> <item2>\nwhere <item1>, <item2> are strings with no spaces");
	  addPositional ("out", "Output directory with the sets of connected items. Each set is named by its lexicographycally smaller item with the added extension " + strQuote ("." + ext));
	  addKey ("subset", "List of items. Item pairs are restricted to this list");
	  addFlag ("pairs", "<out> is a list of pairs: <item1> <item_min>, where <item_min> is lexicographycally smallest item of the cluster");
	}


	void body () const final
	{
		const string inFName     = getArg ("in");
		const string outFName    = getArg ("out");
		const string subsetFName = getArg ("subset");
		const bool   pairs  = getFlag ("pairs");


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
        if (subset. get ())
        {
          if (! subset->containsFast (s1))
            continue;
          if (! subset->containsFast (s2))
            continue;
        }
        if (! contains (items, s1))
          items [s1] = new Item (s1);
        if (! contains (items, s2))
          items [s2] = new Item (s2);
        items [s1] -> DisjointCluster::merge (* items [s2]);
      }
    }
    
    map<const DisjointCluster*,VectorPtr<Item>> clusters;
    for (auto& it : items)
      clusters [it. second->getDisjointCluster ()] << it. second;
      
    unique_ptr<OFStream> fOut;
    if (pairs)
      fOut. reset (new OFStream (outFName));
    for (auto& it : clusters)
    {
      VectorPtr<Item>& cluster = it. second;
      cluster. sort (Named::lessPtr);
      ASSERT (! cluster . empty ());
      if (! pairs)
        fOut. reset (new OFStream (outFName, cluster [0] -> name, ext));
      for (const Item* item : cluster)
        if (pairs)
          *fOut << item->name << '\t' << cluster [0] -> name << endl;
        else
          *fOut << item->name << endl;
      if (! pairs)
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
