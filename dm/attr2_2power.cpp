// attr2_2power.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "matrix.hpp"
#include "dataset.hpp"
using namespace DM_sp;



namespace
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print a " + dmSuff + "-file adding <attr>^power named <attr>_<power>")
    {
      addPositional ("file", dmSuff + "-file");
      addPositional ("attr", "2-way attribute name");
      addPositional ("power", "Power");
    }
	
	
	
	void body () const final
	{
		const string inFName  = getArg ("file");
		const string attrName = getArg ("attr");
		const Real power      = str2real (getArg ("power"));
		ASSERT (power > 0);


    Dataset ds (inFName);    
    const Attr* attr_ = ds. name2attr (attrName);
    ASSERT (attr_);
    const RealAttr2* attr = attr_->asRealAttr2 ();
    ASSERT (attr);
    
    RealAttr2* attr_new = attr->copyAttr (attrName + "_" + toString (power));
    FOR (size_t, i, ds. objs. size ())
	    FOR (size_t, j, ds. objs. size ())
	    {
	      const Real r = attr->matr. get (false, i, j);
	      if (isNan (r))
	        continue;
	      attr_new->matr. put (false, i, j, pow (r, power));
	    }
    
    ds. print (cout);
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


