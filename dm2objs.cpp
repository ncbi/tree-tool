// dm2objs.cpp

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
using namespace DM_sp;



namespace
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print the object names")
    {
      addPositional ("file", dmSuff + "-file");
      addFlag ("comments", "Include comments");
    }
	
	
	
	void body () const
	{
		const string inFName = getArg ("file");
		const bool comments  = getFlag ("comments");


    Dataset ds (inFName);
    
    for (const Obj* obj : ds. objs)
    {
    	cout << obj->name;
    	if (comments)
    	  cout << '\t' << obj->comment;
    	cout << endl;
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


