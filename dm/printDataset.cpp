// printDatset.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
using namespace DM_sp;



namespace
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Read a dataset and print it")
  	{
  	  addPositional ("file", dmSuff + "-file with Number attributes");
  	}
	
	
	
	void body () const final
	{
		const string inFName = getArg ("file");

    
    Dataset ds (inFName);
    ds. qc ();
    
    ds. saveText (cout);
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


