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
	#if 0
	//const S s3 (g ());
	  vector<int> v;
	  cout << v [5] << endl;
	#endif
	  
	  vector<bool> vec;
	  vec. resize (5, false);
	  cout << vec. size () << endl;
	  FOR (size_t, i, 5)
	    cout << vec [i];
	  cout << endl;

	  vec. resize (10, true);
	  cout << vec. size () << endl;
	  FOR (size_t, i, 10)
	    cout << vec [i];
	  cout << endl;
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


