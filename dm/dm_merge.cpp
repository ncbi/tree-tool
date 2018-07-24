// dm_merge.cpp

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
    : Application ("Merge 2 " + dmSuff + "-files and print the resulting " + dmSuff + "-file. Files must have same objects in the same order, and different attributes")
  {
	  addPositional ("file1", "First "  + dmSuff + "-file");
	  addPositional ("file2", "Second " + dmSuff + "-file");
	}



	void body () const final
	{
		const string fName1 = getArg ("file1");
		const string fName2 = getArg ("file2");
		
		
    Dataset ds1 (fName1);
    const Dataset ds2 (fName2);
    if (ds1. objs. size () != ds2. objs. size ())
    	throw runtime_error ("Datasets have different number of objects");
    FFOR (size_t, objNum, ds1. objs. size ())
      if (ds1. objs [objNum] -> name != ds2. objs [objNum] -> name)
	    	throw runtime_error ("Dataset 1 has object " + ds1. objs [objNum] -> name + ", but dataset 2 has object " + ds2. objs [objNum] -> name);

    for (const Attr* attr : ds2. attrs)
    	attr->copyToDataset (ds1);
    
    ds1. saveText (cout);
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



