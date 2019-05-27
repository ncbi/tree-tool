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
  
  
  explicit S (size_t n)
    : a (n)
    {}
  S ()
    { cout << "Default constructor" << endl; }
 ~S () 
    { cout << "Destructor: " << a. size () << endl; }
  S (const S& other)
    : a (other. a)
    { cout << "Copy constructor" << endl; }
  S& operator= (const S &other)
    { a = other. a;
      cout << "Copy assignment" << endl;
      return *this;
    }
  S (S&&) = default;
  S& operator= (S&&) = default;
#if 0
  S (S&& other)
    : a (move (other. a))
    { cout << "Move constructor" << endl; }
  S& operator= (S&& other)
    { a = move (other. a);
      cout << "Move assignment" << endl;
      return *this;
    }
#endif
  
  void f () const
    { throw runtime_error (FUNC  + ": error"); }
}; 



#if 0
void get (S &s)
{
  cout << "Move get" << endl;
  S s1;
  swap (s1, s);
}
#endif



#if 0
void get (const S &s)
{
  cout << "Const get" << endl;
  S s1 (s);
  S s2;
  swap (s1, s2);
}



void get (S &s)
{
  cout << "Update get" << endl;
  S s1;
  swap (s1, s);
}
#endif



#if 0
void testStack1 ()
{
  ERROR;
}



void testStack ()
{
  testStack1 ();
}
#endif



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Test")
  	{
  	  addPositional ("go", "Go");
  	}
  


 	void body () const final
	{
	/*
	  S s (3);
	  get (s);
	*/
	
	#if 0
	  cout << max (1.0, 2.0) << ' ' << max (2.0, 0.0/0.0) << ' ' << max (0.0/0.0, 2.0) << ' ' << max (0.0/0.0, 0.0/0.0) << endl;
	  
	  Vector<int> vec;
	  vec << 3 << 1 << 5;
	  vec. sort ();
	  
	  const Vector<int> vec1 (vec);
	  cout << vec1. searchSorted << endl;

	  Vector<int> vec2;
	  vec2 = vec;
	  cout << vec2. searchSorted << endl;
	  
	  Vector<int> vec3 (move (vec));
	  cout << vec3. searchSorted << ' ' << vec3. size () << endl;
	  cout << vec.  searchSorted << ' ' << vec.  size () << endl;
	#endif

  #if 0	
	  cout << to_string (5) << endl;
	  cout << to_string (5.5) << endl;
	  cout << to_string (-5) << endl;
	  cout << to_string (5.5e5) << endl;
	  cout << endl;
	  
	  cout << toString (5) << endl;
	  cout << toString (5.5) << endl;
	  cout << toString (-5) << endl;
	  cout << toString (5.5e5) << endl;
	#endif
	  
	//testStack ();

  #if 0	
	  double x = 0.0/0.0;
	  cout << x << endl;
	  if (x)
	    cout << "x" << endl;
	  if (! x)
	    cout << "!x" << endl;
	#endif
	
	#if 0
	  Vector<bool> vec (10, false);
	  cout << vec [3] << endl;
	  vec [3] = true;
	  cout << vec [3] << endl;
	  vec [3] = false;
	  cout << vec [3] << endl;
	#endif

  #if 0	
	  const double x = 0.0/0.0;
	  cout << x << endl;
	  cout << (bool) x << endl;
	  cout << (! (bool) x) << endl;
	#endif

  #if 0	
	  double x = 1.0;
	  double y = 2.0;
	  size_t i = 0;
    for (; i < 1000000000000; i++)
    {
      double c = x / y;
  	  if (maximize (x, c))
  	    y = x + 1.0;
  	}
  	cout << i << endl;
  #endif
  
    PRINT (sizeof (int));
    PRINT (sizeof (size_t));
    PRINT (sizeof (float));
    PRINT (sizeof (double));
  }
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


