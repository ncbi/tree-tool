// graph_test.cpp

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;




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

  void saveContent (ostream &os) const final
    { os << ' ' << name; }
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
	  addPositional ("go", "Go");
	}
	
	
	
	void body () const final
  {
    DiGraph g;
  
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
    g. print (cout);  
    cout << endl;
    
    delete berlin_paris;
    delete berlin_paris2;
    delete london_paris;
    delete berlin_berlin;  
    g. print (cout);
    cout << endl;
    
    delete paris;
    g. print (cout);
  }
};


}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}

