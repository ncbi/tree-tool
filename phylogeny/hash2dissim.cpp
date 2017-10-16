// hash2dissim.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
#include "../dm/dataset.hpp"
using namespace DM_sp;



namespace 
{


const string attrName = "cons";



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Convert hashes to a dissimilarity named \"" + attrName + "\" and print a " + dmSuff + "-file")
    {
  	  addPositional ("objects", "File with a list of objects");
  	  addPositional ("hash_dir", "Directory with hashes for each object");
  	  addPositional ("out", "Output " + dmSuff + "-file");
  	}



	void body () const
	{
		const string objectsFName = getArg  ("objects");
		const string hash_dir     = getArg  ("hash_dir");
		const string out          = getArg  ("out");
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
    
    
    typedef Vector<size_t> Hashes;
    
    Vector<Hashes> obj2hashes  (ds. objs. size ());
    {
      Progress prog ((uint) ds. objs. size ());
      FOR (size_t, objNum, ds. objs. size ())
      {
        prog (ds. objs [objNum] -> name);
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
        if (hashes. empty ())
          throw runtime_error ("No hashes for " + ds. objs [objNum] -> name);
      }
    }
    
    auto attr = new PositiveAttr2 (attrName, ds, 6);  // PAR
    {
      Progress prog ((uint) ds. objs. size ());
      FOR (size_t, row, ds. objs. size ())
      {
        prog (ds. objs [row] -> name);
        const Hashes& hash1 = obj2hashes [row];
        FOR (size_t, col, row)
        {
          const Hashes& hash2 = obj2hashes [col];
          const size_t prots = hash1.  getIntersectSize (hash2);
          ASSERT (prots <= hash1. size ());
          ASSERT (prots <= hash2. size ());
        #if 0
          const Real dissim = - log ((Real) prots / (Real) min (hash1. size (), hash2. size ()));
        #else
          const Real dissim = - 0.5 * (  log ((Real) prots / (Real) hash1. size ()) 
                                       + log ((Real) prots / (Real) hash2. size ()) 
                                      );
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



