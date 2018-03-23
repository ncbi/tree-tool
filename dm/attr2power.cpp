// attr2power.cpp

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
    : Application ("Print a " + dmSuff + "-file adding <attr>^power named <attr>_power")
    {
      addPositional ("file", dmSuff + "-file");
      addPositional ("attr", "Attribute name");
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
    const RealAttr1* attr = attr_->asRealAttr1 ();
    ASSERT (attr);
    
    auto attr1 = new RealAttr1 (attrName + toString (power), ds, attr->decimals);
    FOR (size_t, objNum, ds. objs. size ())
    {
      const Real r = (*attr) [objNum];
      if (isNan (r))
        continue;
      (*attr1) [objNum] = pow (r, power);
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


