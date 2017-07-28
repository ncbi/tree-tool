// mergeSimilarity.cpp

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;
#include "matrix.hpp"
#include "dataset.hpp"
using namespace DM_sp;



namespace 
{


struct ThisApplication : Application
{
  ThisApplication () 
    : Application ("Print a " + dmSuff + "-file with merged <Similarity> attrributes of a list of " + dmSuff + "-files")
    {
  	  addPositional ("file_list", "List of " + dmSuff + "-files without \"" + dmSuff + "\" with a binary <Similarity> attribute");
  	  addKey ("attr", "Similarity attribute", "Similarity");
  	}



	void body () const
	{
		const string file_list = getArg ("file_list");
		const string attrName  = getArg ("attr");
		
		
		Set<string> objNames;
    {
      // Pass 1
      LineInput li (file_list);
      while (li. nextLine ())
      {
        const Dataset ds (li. line);
        ds. qc ();
        ASSERT (ds. getUnitMult ());
        for (const Obj* obj : ds. objs)
          objNames << obj->name;
      }
    }

		
    Dataset ds_new;
    for (const string& objName : objNames)
      ds_new. appendObj (objName);
    ds_new. setName2objNum ();
    ASSERT (! ds_new. objs. empty ());
    ASSERT (ds_new. objs. size () == objNames. size ());
    

    auto sim_new = new RealAttr2 (attrName, ds_new, 6);  // PAR
    sim_new->matr. putAll (0);
    Matrix count (false, sim_new->matr, false);
    count. putAll (0);
    {
      // Pass 2
      LineInput li (file_list);
      Progress prog;
      while (li. nextLine ())
      {
        prog (li. line);
        const Dataset ds (li. line);
        
        const Attr* attr = ds. name2attr (attrName);
        ASSERT (attr);
        const RealAttr2* sim = attr->asRealAttr2 ();
        ASSERT (sim);

        Matrix& matr = const_cast <RealAttr2*> (sim) -> matr;
        ASSERT (matr. defined ());
        Matrix* rowMean = nullptr;
        Real totalMean;
        matr. centerSimilarity (rowMean, totalMean);
        delete rowMean;
        
        // Normalization
        const Real k = matr. getTrace () / (Real) ds. objs. size ();
        ASSERT (positive (k));
        matr. putProdAll (1 / k);
        
        FOR (size_t, row, ds. objs. size ())
        {
          const size_t row_new = ds_new. getName2objNum (ds. objs [row] -> name);
          ASSERT (row_new != NO_INDEX);
          FOR (size_t, col, ds. objs. size ())
          {
            const size_t col_new = ds_new. getName2objNum (ds. objs [col] -> name);
            ASSERT (col_new != NO_INDEX);
            const Real val = matr. get (false, row, col);
            ASSERT (! isNan (val));
            sim_new->matr. putInc (false, row_new, col_new, val);  
            count. putInc (false, row_new, col_new, 1.0);  
          }
        }
      }
    }


    // Averaging
    FOR (size_t, row, ds_new. objs. size ())
      FOR (size_t, col, ds_new. objs. size ())
      {
        const Real val = sim_new->matr. get (false, row, col);
        if (! isNan (val))
          sim_new->matr. put (false, row, col, val / count. get (false, row, col));  
      }


    ds_new. print (cout);
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



