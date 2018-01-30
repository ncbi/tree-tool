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
  : Application ("Print hash codes for a list of strings")
  { 
	  addPositional ("in", "List of strings");
	  addPositional ("out", "Output file of sorted unique hashes");
	}



	void body () const final
	{
    const string in  = getArg ("in");  
    const string out = getArg ("out"); 
    ASSERT (! out. empty ()); 
        
    
    hash<string> str_hash;
    static_assert (sizeof (size_t) == 8, "Size of size_t must be 8 bytes");
    
    Vector<size_t> hashes;  hashes. reserve (10000);   // PAR
    {
	    LineInput f (in);
	    while (f. nextLine ())
	      hashes << str_hash (f. line);
	  }

    hashes. sort ();
    hashes. uniq ();
    
    {
      OFStream f (out);
      f << hashes;
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



