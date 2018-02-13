// str2hash.cpp

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;




namespace 
{
  
  
  
struct ThisApplication : Application
{
  ThisApplication ()
  : Application ("Print sorted unique hash codes for a list of strings form stdin", false)
  {}



	void body () const final
	{
    hash<string> str_hash;
    static_assert (sizeof (size_t) == 8, "Size of size_t must be 8 bytes");
    
    Vector<size_t> hashes;  hashes. reserve (10000);   // PAR
	  {
	  	string s;
	    while (cin >> s)
	      hashes << str_hash (s);;
	  }

    hashes. sort ();
    hashes. uniq ();
    
  	cout << hashes;
  }  
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



