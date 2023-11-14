// graph_test.cpp

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
*   Test of graph.{hpp,cpp}
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


struct City : DiGraph::Node
{
  string name;

  City (DiGraph &g,
        const string &name_arg) 
    : DiGraph::Node (g)
    , name (name_arg)
    {}

  string getName () const final
    { return name; }
};



struct Distance : DiGraph::Arc
{
  int miles;

  Distance (City* start_arg,
            City* end_arg,
            int miles_arg) 
    : DiGraph::Arc (start_arg, end_arg)
    , miles (miles_arg)
    {}


  void saveContent (ostream &os) const final
    { os << miles; }
};



struct ThisApplication : Application
{
	ThisApplication ()
	: Application ("DiGraph test")
	{
    version = VERSION;
	  addPositional ("go", "Go");
	}
	
	
	
	void body () const final
  {
  	{
	    DiGraph g;
	    g. qc ();
	  
		  // Node's
		  auto london = new City (g, "London");
		  auto paris  = new City (g, "Paris");
		  auto berlin = new City (g, "Berlin");
		  auto madrid = new City (g, "Madrid");
		  // Arc's
		  auto london_paris  =   new Distance (london, paris, 100);
		/*auto london_berlin =*/ new Distance (london, berlin, 200);
		/*auto london_madrid =*/ new Distance (london, madrid, 50);
		  auto berlin_paris  =   new Distance (berlin, paris, 30);
		  auto berlin_paris2 =   new Distance (berlin, paris, 35);
		  auto berlin_berlin =   new Distance (berlin, berlin, 0);
		/*auto berlin_madrid =*/ new Distance (berlin, madrid, 200);  
	  
	    g. qc ();
	    g. saveText (cout);  
	    cout << endl;
	    
	    delete berlin_paris;
	    delete berlin_paris2;
	    delete london_paris;
	    delete berlin_berlin;  
	    g. qc ();
	    g. saveText (cout);
	    cout << endl;
	    
	    delete paris;
	    g. qc ();
	    g. saveText (cout);
	  }
	  
    
    // TopologicalSort
  	{
	    DiGraph g;
		  // Node's
		  auto london = new City (g, "London");
		  auto paris  = new City (g, "Paris");
		  auto berlin = new City (g, "Berlin");
		  auto madrid = new City (g, "Madrid");
		  // Arc's
      new Distance (london, paris, 100);
      new Distance (london, berlin, 200);
      new Distance (london, madrid, 50);
      new Distance (berlin, paris, 30);
      new Distance (berlin, paris, 35);
      new Distance (berlin, berlin, 0);
	    g. qc ();

      const size_t s = g. nodes. size ();
	    TopologicalSort ts (g, true);
	    VectorPtr<DiGraph::Node> nodes;
	    while (const DiGraph::Node* n = ts.  getFront ()) 
	    	nodes << n;
	    g. qc ();
	    QC_ASSERT (s == g. nodes. size ());
	    QC_ASSERT (nodes. size () == s);
	    nodes. sort ();
	    nodes. uniq ();
	    QC_ASSERT (nodes. size () == s);
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

