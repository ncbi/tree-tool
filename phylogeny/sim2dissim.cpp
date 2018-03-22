// sim2dist.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
#include "../dm/matrix.hpp"
#include "../dm/dataset.hpp"
using namespace DM_sp;
#include "evolution.hpp"
using namespace DistTree_sp;



namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Convert sequence similarity to evoluiton dissimilarity")
    {
  	  addPositional ("file", dmSuff + "-file without the extension");
  	  addPositional ("attrName", "Attribute name of real-valued object-object table in the " + dmSuff + "-file");
  	}



	void body () const final
	{
		const string fName    = getArg ("file");
		const string attrName = getArg ("attrName");
		
		
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
    dist->matr = sim->matr;
    dist->decimals = 6;  // PAR
    Matrix& matr = dist->matr;
    FOR (size_t, i, matr. rowsSize ())
    {
    	const Real selfScore1 = matr. get (false, i, i);
      FOR (size_t, j, i)
      {
      	// Cf. prots_pair2dissim.cpp
      	const Real score = matr. get (false, i, j);
	    	const Real selfScore2 = matr. get (false, j, j);
	      const Real dissim = score < 0 
	                            ? INF 
	                            : pow (intersection2dissim (selfScore1, selfScore2, score, 0, 0.5, true), 0.35);  // PAR
        matr. putSymmetric (i, j, dissim);
      }
    }    
    FOR (size_t, i, matr. rowsSize ())
     	matr. put (false, i, i, 0.0); 
       
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



