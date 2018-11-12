// setRandOrd.cpp

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;



namespace
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print a random reordering of the list\nRequires: Line length < 1024")
  	{
  	  addPositional ("items", "File with items (end-of-line separated)");
  	}



	void body () const final
	{
		const string items = getArg ("items");


    StringVector vec;
    {
	    LineInput in (items, 100 * 1024);  // PAR
	    vec = in. getVector ();
	  }
	  vec. randomOrder ();	  
	  for (const string& s : vec)
	    cout << s << endl;
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


