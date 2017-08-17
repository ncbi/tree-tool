// attr2_multiply.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "numeric.hpp"
#include "matrix.hpp"
#include "dataset.hpp"
using namespace DM_sp;



namespace 
{


const string clusterAttrName = "Cluster";



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Multiply a ReaAttr2 and print the resulting data set")
    {
  	  addPositional ("file", dmSuff + "-file");
  	  addPositional ("attr2Name", "Attribute name of an object-object table in the " + dmSuff + "-file");
  	  addPositional ("coefficient", "Coeffciient to multiply <attr2Name> by");
  	}



	void body () const
	{
		const string fName     = getArg ("file");
		const string attr2Name = getArg ("attr2Name");
		const Real coefficient = str2real (getArg ("coefficient"));
		ASSERT (! isNan (coefficient));
		
		
    Dataset ds (fName);
    
    // dist
    RealAttr2* dist = nullptr;
    {
      const Attr* attr = ds. name2attr (attr2Name);
      ASSERT (attr);
      dist = const_cast <RealAttr2*> (attr->asRealAttr2 ());
    }
    ASSERT (dist);
    
    dist->decimals += (uint) max<long> (0, - DM_sp::round (log10 (coefficient)));
    
    Matrix& matr = dist->matr;
    FOR (size_t, row, ds. objs. size ())
    FOR (size_t, col, ds. objs. size ())
      matr. putProd (false, row, col, coefficient);
      
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



