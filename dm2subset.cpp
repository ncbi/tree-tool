// dm2subset.cpp

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
    : Application ("Print a subset of a " + dmSuff + "-file")
  {
	  addFlag ("exclude", "Exclude the list of objects, otherwise include");
	  addPositional ("file", dmSuff + "-file");
	  addPositional ("objNameFName", "File with object names defining the subset");
	}



	void body () const
	{
		const bool   exclude      = getFlag ("exclude");
		const string fName        = getArg ("file");
		const string objNameFName = getArg ("objNameFName");
		
		
    Dataset ds (fName);
    
    Set<string> objNames;
    {
      LineInput li (objNameFName);
      while (li. nextLine ())
        objNames << li. line;
    }
    
    Sample sm (ds);
    FOR (size_t, row, ds. objs. size ())
      if (objNames. contains (ds. objs [row] -> name) == exclude)
        sm. mult [row] = 0;
    sm. finish ();    

    sm. save (VectorPtr<Attr> (ds. attrs), cout);  
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



