// sim2dist.cpp

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


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Convert similarity to distance")
    {
  	  addPositional ("file", dmSuff + "-file without the extension");
  	  addPositional ("attrName", "Attribute name of real-valued object-object table in the " + dmSuff + "-file");
  	  addFlag ("sqr", "Make squared didtances");
  	}



	void body () const
	{
		const string fName    = getArg ("file");
		const string attrName = getArg ("attrName");
		const bool makeSqr    = getFlag ("sqr");
		
		
    Dataset ds (fName);
    
    // sim
    const RealAttr2* sim = nullptr;
    {
      const Attr* attr = ds. name2attr (attrName);
      ASSERT (attr);
      sim = attr->asRealAttr2 ();
    }
    ASSERT (sim);    
    ASSERT (! sim->asPositiveAttr2 ());
    
	
	  // Check sim->matr  
	  {
      Real maxCorrection;
      size_t row_bad, col_bad;
      const_cast <RealAttr2*> (sim) -> matr. symmetrize (maxCorrection, row_bad, col_bad);
      if (maxCorrection > 2 * pow (10, - (Real) sim->decimals))
        ds. comments << "maxCorrection = " + toString (maxCorrection) + " at " + ds. objs [row_bad] -> name + ", " + ds. objs [col_bad] -> name;
    }

	  {
  	  size_t row;
  	  size_t col;
  	  if (sim->matr. existsMissing (false, row, col))
      {
        cout << attrName << '[' << ds. objs [row] -> name << "," << ds. objs [col] -> name << "] is missing" << endl;
        exit (1);
      }  
    }
	

    auto dist = new PositiveAttr2 (sim->name + "_dist", ds, sim->decimals); 
    {
      dist->matr = sim->matr;
      dist->matr. similarity2sqrDistance ();
      if (! makeSqr)
        dist->matr. sqrtAll ();
    }
    
    delete sim;
    
    ds. qc ();

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



