// test.cpp

#include <iostream>
using namespace std;
#include <common.hpp>
#include <common.inc>
using namespace common;


class City: public Graph::Node
{
  string name;

public:
  City (const string &name_arg) :
    Graph::Node (),
    name (name_arg)
    {}


  void printContent () const
    { cout << name; }
};



class Distance: public Graph::Edge
{
  int miles;

public:
  Distance (City* start_arg,
            City* end_arg,
            int miles_arg) :
    Graph::Edge (start_arg, end_arg),
    miles (miles_arg)
    {}


  void printContent () const
    { cout << miles; }
};



int main (int argc,
          char* argv[])
{
  Graph g;

  City* london = new City ("London");
  City* paris = new City ("Paris");
  City* berlin = new City ("Berlin");
  City* madrid = new City ("Madrid");
  g. nodes. push_back (london);
  g. nodes. push_back (paris);
  g. nodes. push_back (berlin);
  g. nodes. push_back (madrid);

  Distance* london_paris = new Distance (london, paris, 100);
  Distance* london_berlin = new Distance (london, berlin, 200);
  Distance* london_madrid = new Distance (london, madrid, 50);
  Distance* berlin_paris = new Distance (berlin, paris, 30);
  Distance* berlin_paris2 = new Distance (berlin, paris, 35);
  Distance* berlin_berlin = new Distance (berlin, berlin, 0);
  Distance* berlin_madrid = new Distance (berlin, madrid, 200);

  g. print ();

  cout << endl;
  delete berlin_paris;
  delete berlin_paris2;
  delete london_paris;
  delete berlin_berlin;

  g. print ();

  
  return 0;
}

