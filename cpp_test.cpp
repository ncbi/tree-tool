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



struct Struct
{
  string s;
  Struct ()
    : s ("abcdefgihjkhjklhjkhjlhjhjkhjklhjklhjkljklhjhhhhhhhhhhhhhhhhhhhhhhhhh")
    {}
};



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Test")
  	{
  	  addPositional ("go", "Go");
  	}
  


 	void body () const
	{
	  constexpr size_t len = 10000000;  // PAR
	  Vector<ulong> vec;  vec. reserve (len);
	  Rand rand;
	  FOR (size_t, i, len)
	    vec << rand. get ((ulong) len);
	  {
      const Chronometer_OnePass cop ("map");  
  	  map<ulong,size_t> m;
  	  FOR (size_t, i, len)
  	    m [vec [i]] = i;
  	  FOR (size_t, i, len)
  	    m [vec [i]] ++;
  	}
	  {
      const Chronometer_OnePass cop ("unordered_map");  
  	  unordered_map<ulong,size_t> m (len);
  	  FOR (size_t, i, len)
  	    m [vec [i]] = i;
  	  FOR (size_t, i, len)
  	    m [vec [i]] ++;
  	}
	  {
      const Chronometer_OnePass cop ("unordered_map");  
  	  unordered_map<ulong,size_t> m ;
  	  m. reserve (len);
  	  FOR (size_t, i, len)
  	    m [vec [i]] = i;
  	  FOR (size_t, i, len)
  	    m [vec [i]] ++;
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


