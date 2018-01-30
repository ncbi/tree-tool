// hash_request2dissim.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;



namespace 
{


typedef Vector<size_t> Hashes;



Hashes getHashes (const string &fName)
{
  LineInput hf (fName);
  Hashes hashes;  hashes. reserve (10000);  // PAR
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
    throw runtime_error ("No hashes for " + fName);
    
  return hashes;
}
    



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Compute hash dissimilarities for requested pairs")
    {
  	  addPositional ("pairs", "File with pairs of objects");
  	  addPositional ("hash_dir", "Directory with hashes for each object");
  	  addPositional ("out", "Output file");
  	}



	void body () const final
	{
		const string pairsFName = getArg  ("pairs");
		const string hash_dir   = getArg  ("hash_dir");
		const string out        = getArg  ("out");
		ASSERT (! hash_dir. empty ());
		ASSERT (! out. empty ());
		
		
    // Cf. hash2dissim.cpp
    LineInput input (pairsFName);
    OFStream output (out);
    ONumber on (output, 6, true);  // PAR
    while (input. nextLine ())
    {
      istringstream iss (input. line);
      string name1, name2;
      iss >> name1 >> name2;
      ASSERT (name1 != name2);
      const Hashes h1 (getHashes (hash_dir + "/" + name1));
      const Hashes h2 (getHashes (hash_dir + "/" + name2));
      const size_t intersection = h1.  getIntersectSize (h2);
      ASSERT (intersection <= h1. size ());
      ASSERT (intersection <= h2. size ());
      const double dissim = - 0.5 * (  log ((double) intersection / (double) h1. size ()) 
                                     + log ((double) intersection / (double) h2. size ()) 
                                    );
      output << name1 << '\t' << name2 << '\t' << dissim;
      if (verbose ())
        output << '\t' << intersection;
      output << endl;
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



