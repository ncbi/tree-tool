// hash2dissim.cpp  

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
#include "../dm/dataset.hpp"
using namespace DM_sp;
#include "evolution.hpp"
using namespace DistTree_sp;



namespace 
{


const string attrName = "cons";



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Convert hashes to a dissimilarity named \"" + attrName + "\" and print a " + dmSuff + "-file")
    {
    	// Input
  	  addPositional ("objects", "File with a list of objects");
  	  addPositional ("hash_dir", "Directory with hashes for each object");
  	  addKey ("intersection_min", "Min. number of common hashes to compute distance", "50");
  	  addKey ("ratio_min", "Min. ratio of hash sizes (0..1)", "0.5");
  	  // Output
  	  addPositional ("out", "Output " + dmSuff + "-file without " + dmSuff);
  	}



	void body () const final
	{
		const string objectsFName     = getArg  ("objects");
		const string hash_dir         = getArg  ("hash_dir");
		const size_t intersection_min = str2<size_t> (getArg ("intersection_min"));
		const Prob   hashes_ratio_min = str2real (getArg ("ratio_min"));
		const string out              = getArg ("out");
		ASSERT (isProb (hashes_ratio_min));
		ASSERT (! out. empty ());
		
		
    Dataset ds;
    {
  		Set<string> objNames;
      {
        LineInput f (objectsFName);
        while (f. nextLine ())
        {
          trim (f. line);
          objNames << f. line;
        }
      }
      cerr << "# Objects: " << objNames. size () << endl;  
      for (const string& name : objNames)
        ds. appendObj (name);
    }        
    ds. qc ();
    
    
    Vector<Hashes> obj2hashes;  obj2hashes. reserve (ds. objs. size ());
    {
      Progress prog (ds. objs. size ());
      FFOR (size_t, objNum, ds. objs. size ())
      {
        prog (ds. objs [objNum] -> name);
        obj2hashes [objNum] << Hashes (hash_dir + "/" + ds. objs [objNum] -> name);
      }
    }
    
    auto attr = new PositiveAttr2 (attrName, ds, 6);  // PAR
    {
      Progress prog (ds. objs. size ());
      FFOR (size_t, row, ds. objs. size ())
      {
        prog (ds. objs [row] -> name);
        const Hashes& hash1 = obj2hashes [row];
        FOR (size_t, col, row)
        {
          const Hashes& hash2 = obj2hashes [col];
          const Real dissim = hash1. getDissim (hash2, intersection_min, hashes_ratio_min);
          attr->putSymm (row, col, dissim);
        }
        attr->put (row, row, 0);
      }
    }
    
    ds. qc ();
    {
      OFStream f (out + dmSuff);
      ds. saveText (f);    
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



