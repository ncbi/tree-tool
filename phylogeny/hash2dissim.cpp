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
  	  addKey ("hashes_min", "Min. number of common hashes to compute distance", "50");
  	  // Output
  	  addPositional ("out", "Output " + dmSuff + "-file without " + dmSuff);
  	}



	void body () const
	{
		const string objectsFName = getArg  ("objects");
		const string hash_dir     = getArg  ("hash_dir");
		const string out          = getArg  ("out");
		const size_t hashes_min   = str2<size_t> (getArg ("hashes_min"));
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
      Progress prog ((uint) ds. objs. size ());
      FFOR (size_t, objNum, ds. objs. size ())
      {
        prog (ds. objs [objNum] -> name);
        obj2hashes [objNum] << Hashes (hash_dir + "/" + ds. objs [objNum] -> name);
      #if 0
        LineInput hf (hash_dir + "/" + ds. objs [objNum] -> name);
        auto& hashes = obj2hashes [objNum];
        ASSERT (hashes. empty ());
        hashes. reserve (10000);  // PAR
        size_t prev = 0;
        while (hf. nextLine ())
        {
          const size_t hash = str2<size_t> (hf. line);
          ASSERT (hash);
          if (hash <= prev)
            throw runtime_error ("Hash " + hf. line + " is not greater than the previous hash " + toString (prev));
          hashes << hash;
          prev = hash;
        } 
      //if (hashes. empty ())
        //throw runtime_error ("No hashes for " + ds. objs [objNum] -> name);
      #endif
      }
    }
    
    auto attr = new PositiveAttr2 (attrName, ds, 6);  // PAR
    {
      Progress prog ((uint) ds. objs. size ());
      FFOR (size_t, row, ds. objs. size ())
      {
        prog (ds. objs [row] -> name);
        const Hashes& hash1 = obj2hashes [row];
        FOR (size_t, col, row)
        {
          const Hashes& hash2 = obj2hashes [col];
          const Real dissim = hash1. getDissim (hash2, hashes_min);
        #if 0
          const size_t common = hash1.  getIntersectSize (hash2);
          ASSERT (common <= hash1. size ());
          ASSERT (common <= hash2. size ());
          const Real dissim = common >= hashes_min
                                ? - 0.5 * (  log ((Real) common / (Real) hash1. size ()) 
                                           + log ((Real) common / (Real) hash2. size ()) 
                                          )
                                : NAN;
        #endif
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



