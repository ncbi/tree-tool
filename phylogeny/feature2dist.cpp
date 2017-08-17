// feature2dist.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/dataset.hpp"
using namespace DM_sp;



namespace 
{


const string distName = "dist";



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Convert feature sets to a distance '" + distName + "' and print a " + dmSuff + "-file")
    {
  	  addPositional ("objects", "File with a list of object files with features");
  	  addPositional ("objects_dir", "Directory with <object> files containing features");
  	}



	void body () const
	{
		const string objectsFName = getArg ("objects");
		const string objects_dir  = getArg ("objects_dir");
		
		
		typedef  Set<string/*feature*/>  Features;
		
    Dataset ds;
    Vector<Features> obj2features;  

    size_t features_min = numeric_limits<size_t>::max();
    size_t features_max = 0;
    {
      LineInput objF (objectsFName);
      while (objF. nextLine ())
      {
        trim (objF. line);
        ds. appendObj (objF. line);
        Features features;
        {
          LineInput featF (objects_dir + "/" + objF. line);
          while (featF. nextLine ())
          {
            trim (featF. line);
            features << featF. line;
          }
        }
        obj2features << features;
        minimize (features_min, features. size  ());
        maximize (features_max, features. size  ());
      }
    }
    ASSERT (ds. objs. size () == obj2features. size ());
    cerr << "# Objects: " << ds. objs. size () << endl;
    cerr << "Min. features: " << features_min << endl;
    cerr << "Max. features: " << features_max << endl;
    ASSERT (features_min);


    auto attr = new PositiveAttr2 (distName, ds, 4);  // PAR
    Progress prog ((uint) ds. objs. size ());
    FOR (size_t, i, ds. objs. size ())
    {
      prog ();
      const Features& f1 = obj2features [i];
      FOR_START (size_t, j, i + 1, ds. objs. size ())
      {
        const Features& f2 = obj2features [j];
        const size_t intersection = f1. intersectSize (f2);
      #if 0
        const Prob c1 = (Real) intersection / (Real) f1. size ();
        const Prob c2 = (Real) intersection / (Real) f2. size ();
        const Real dist = -10 * (log (c1) + log (c2));
      #else
        Features f (f2);
        f << f1;
        const Prob jaccard = (Real) intersection / (Real) f. size ();
        ASSERT (isProb (jaccard));
        const Real dist = -10 * log (jaccard);
      #endif
        attr->putSymm (i, j, dist);
      }
    }
    FOR (size_t, row, ds. objs. size ())
      attr->put (row, row, 0);

    
    ds. qc ();
    ds. saveText (cout);    
	}
};



}  // namespace



int main(int argc, 
         const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



