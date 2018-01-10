// trav.cpp

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;



namespace
{
  
  
struct S : Root
{
  vector<int> a;
  vector<int> b;
  
  
  S ()
    : a (2)
    , b (3)
    { cout << "Default constructor" << endl; }
 ~S () 
    { cout << "Destructor" << endl; }
  S (const S& other)
    : a (other. a)
    , b (other. b)
    { cout << "Copy constructor" << endl; }
  S (S&& other)
    : a (move (other. a))
    , b (move (other. b))
    { cout << "Move constructor" << endl; }
}; 



S f ()
{
  S s;
  return s;
}


inline S g ()
{
  return f ();
}



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Test")
  	{
  	  addPositional ("go", "Go");
  	}
  


	void body () const
	{
	  Vector<double> vec;  vec. reserve (1000000);
	  cout << vec. capacity () << endl;
	  Rand rand;
	  FOR (size_t, i, vec. capacity ())
	    vec << rand. getProb ();
	  cout << vec. capacity () << endl;
	  cout << vec. size () << endl;
	
	  {    
      const Chronometer_OnePass chron ("sort");
  	  for (size_t i = 0; i < 100000; i++)
  	    vec. sort ();
  	}

	  {    
      const Chronometer_OnePass chron ("nth_element");
  	  for (size_t i = 0; i < 1000; i++)
  	    std::nth_element (vec. begin (), vec. begin () + 7000, vec. end ());
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


