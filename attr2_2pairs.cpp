// attr2_2pairs.cpp

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
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
    : Application ("Print pairs for an Attr2: <objName1>\t<objName2>\t<value>")
    {
  	  addPositional ("file", dmSuff + "-file");
  	  addPositional ("attr2Name", "Attribute name of an object-object table in the " + dmSuff + "-file");
  	  addFlag ("symmetric", "Attribute is symmetric, print only lines where <objName1> <= <objName2>");
  	  addFlag ("diagonal", "Print the diagonal elements");
  	}



	void body () const
	{
		const string fName     = getArg ("file");
		const string attr2Name = getArg ("attr2Name");
		const bool symmetric   = getFlag ("symmetric");
		const bool diagonal    = getFlag ("diagonal");
		
		
    Dataset ds (fName);
    
    // dist
    const RealAttr2* dist = nullptr;
    {
      const Attr* attr = ds. name2attr (attr2Name);
      ASSERT (attr);
      dist = attr->asRealAttr2 ();
    }
    ASSERT (dist);

    // Similarity -> distance ??
    
    Matrix& matr = const_cast <RealAttr2*> (dist) -> matr;
  //matr. similarity2sqrDistance ();
    
    FOR (size_t, row, ds. objs. size ())
    FOR (size_t, col, ds. objs. size ())
      if (! symmetric || ds. objs [row] -> name <= ds. objs [col] -> name)
        if (row != col || diagonal)
          cout         << ds. objs [row] -> name 
               << '\t' << ds. objs [col] -> name 
               << '\t' << matr. get (false, row, col)
               << endl;
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



