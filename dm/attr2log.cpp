// attr2log.cpp

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
    : Application ("Print a " + dmSuff + "-file adding log(attr) named \"log_\" + attr")
    {
      addPositional ("file", dmSuff + "-file");
      addPositional ("attr", "Attribute name");
      addFlag ("inf2nan", "Convert -INF to NaN");
    }
	
	
	
	void body () const
	{
		const string inFName  = getArg ("file");
		const string attrName = getArg ("attr");
		const bool inf2nan    = getFlag ("inf2nan");


    Dataset ds (inFName);    
    const Attr* attr_ = ds. name2attr (attrName);
    ASSERT (attr_);
    const RealAttr1* attr = attr_->asRealAttr1 ();
    ASSERT (attr);
    
    auto logAttr = new RealAttr1 ("log_" + attrName, ds, attr->decimals + 1);
    FOR (size_t, objNum, ds. objs. size ())
    {
      const Real r = (*attr) [objNum];
      if (isNan (r))
        continue;
      ASSERT (r >= 0);
      Real res = log (r);
      if (inf2nan && res == -INF)
      	res = NaN;
      (*logAttr) [objNum] = res;
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


