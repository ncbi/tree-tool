// dm2space.cpp

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
    : Application ("Print a space of a " + dmSuff + "-file")
  {
	  addPositional ("file", dmSuff + "-file");
	  addPositional ("attrNameFName", "File with attribute names defining the space");
	  addFlag ("exclude", "Exclude the list of attributes, otherwise include");
	}



	void body () const final
	{
		const string fName         = getArg ("file");
		const string attrNameFName = getArg ("attrNameFName");
		const bool   exclude       = getFlag ("exclude");
		
		
    Dataset ds (fName);
    
    Set<string> attrNames;
    {
      LineInput li (attrNameFName);
      while (li. nextLine ())
      {
      	trim (li. line);
      	if (! li. line. empty ())
          attrNames << li. line;
      }
    }
    
    VectorPtr<Attr> attrs;  attrs. reserve (ds. attrs. size ());
    for (const string& name : attrNames)
    	if (const Attr* attr = ds. name2attr (name))
    		attrs << attr;
    	else
    	  throw runtime_error ("Attribute " + name + " is not in the dataset");
    	  
    if (exclude)
    {
    	VectorPtr<Attr> attrs1 (ds. attrs);
    	attrs. sort ();
    	attrs1. setMinus (attrs);
    	attrs = attrs1;
    }
    
    const Sample sm (ds);
    sm. save (attrs, cout);  
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



