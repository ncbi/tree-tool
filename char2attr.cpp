// char2attr.cpp

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
    : Application ("Add a Boolean attribute for the object names containing a specified substring")
  	{
  	  addFlag("equal_multiplicity", "the objects with the substring should have teh same weight as the other objects");
  	  addPositional ("in_file", "Input " + dmSuff + "-file");
  	  addPositional ("substr", "Object name substring");	
  	  addPositional ("attr_name", "Name of the new Boolean attribute");
  	  addPositional ("out_file", "Output " + dmSuff + "-file");
  	}
	
	
	
	void body () const final
	{
		const bool equal_multiplicity = getFlag ("equal_multiplicity");
		const string in_file          = getArg  ("in_file");
		const string substr           = getArg ("substr");
		const string attr_name        = getArg ("attr_name");
		const string out_file         = getArg ("out_file");
		ASSERT (! substr. empty ());
		ASSERT (! attr_name. empty ());


    Dataset ds (in_file);
    ExtBoolAttr1* attr = new ExtBoolAttr1 (attr_name, ds);   
    FOR (size_t, i, ds. objs. size ())
      (*attr) [i] = (ebool) contains (ds. objs. at (i) -> name, substr);
      
    if (equal_multiplicity)
    {
      Real classMult = 0;
      Real otherMult = 0;
      FOR (size_t, i, ds. objs. size ())
        if ((*attr) [i])
          classMult += ds. objs. at (i) -> mult;
        else
          otherMult += ds. objs. at (i) -> mult;
      ASSERT (classMult > 0);
      const Real ratio = otherMult / classMult;

      FOR (size_t, i, ds. objs. size ())
        if ((*attr) [i])
          const_cast <Obj*> (ds. objs. at (i)) -> mult *= ratio;
    }
      
    {
      OFStream of ("", out_file, dmExt);
      ds. saveText (of);
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


