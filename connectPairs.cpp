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
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;
#include "version.inc"



namespace 
{
  
  
static const string ext = "component";



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
	}


	void body () const final
	{
		const string in  = getArg ("in");
		const string out = getArg ("out");


    unordered_map<string/*Item::name*/,Item*> items;  items. rehash (100000);  // PAR
      // Not delete'd
    {
      LineInput fIn (in, 1024 * 1024);  // PAR
      string s1, s2;
      Istringstream iss;
      while (fIn. nextLine ())
      {
        iss. reset (fIn. line);
        s2. clear ();
        iss >> s1 >> s2;
        QC_ASSERT (! s2. empty ());
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
      
    for (auto& it : clusters)
    {
      VectorPtr<Item>& cluster = it. second;
      cluster. sort (Named::lessPtr);
      ASSERT (! cluster . empty ());
      OFStream fOut (out, cluster [0] -> name, ext);
      for (const Item* item : cluster)
        fOut << item->name << endl;
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
