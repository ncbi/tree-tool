// cpp_test.cpp

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
	   const double d = 1.5;
	   const double r = abs (d);  // non-standard
	   cout << r << endl;
  #if 0
    int i = d;
    cout << i << endl;
  #endif  
  }
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


