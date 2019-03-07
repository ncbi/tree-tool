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
  
  void f () const
    { throw runtime_error (FUNC  + ": error"); }
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



template <typename K, typename V>
void um_stat (const unordered_map<K,V> &um)
{
  #define STAT(f)  cout << #f " = " << um. f () << endl;
//STAT (bucket_count);
  STAT (load_factor);
//STAT (max_load_factor);
  #undef STAT
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
	  double d = str2<double> ("1e0");
	  cout << d << endl;
	  
	  {
  	  istringstream iss ("1e");
  	  iss >> d;
  	  cout << d << ' ' << iss. eof () << ' ' << iss. fail () << ' ' << iss. rdstate () << endl;
  	}

	  {
  	  istringstream iss ("1e1");
  	  iss >> d;
  	  cout << d << ' ' << iss. eof () << ' ' << iss. fail () << ' ' << iss. rdstate () << endl;
  	}
	#else
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
	    cout << endl;
      const Chronometer_OnePass cop ("unordered_map: constructor");  
  	  unordered_map<ulong,size_t> m (len);
  	  um_stat (m);
  	  FOR (size_t, i, len)
  	    m [vec [i]] = i;
  	  FOR (size_t, i, len)
  	    m [vec [i]] ++;
  	  cout << endl;
  	  um_stat (m);
  	}
	  {
	    cout << endl;
      const Chronometer_OnePass cop ("unordered_map: reserve");  
  	  unordered_map<ulong,size_t> m;
  	  m. reserve (len);
  	  um_stat (m);
  	  FOR (size_t, i, len)
  	    m [vec [i]] = i;
  	  FOR (size_t, i, len)
  	    m [vec [i]] ++;
  	  cout << endl;
  	  um_stat (m);
  	}
	  {
	    cout << endl;
      const Chronometer_OnePass cop ("unordered_map: rehash");  
  	  unordered_map<ulong,size_t> m;
  	  m. rehash (len);
  	  um_stat (m);
  	  FOR (size_t, i, len)
  	    m [vec [i]] = i;
  	  FOR (size_t, i, len)
  	    m [vec [i]] ++;
  	  cout << endl;
  	  um_stat (m);
  	}
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


