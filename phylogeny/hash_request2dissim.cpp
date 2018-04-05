// hash_request2dissim.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "evolution.hpp"
using namespace DistTree_sp;



namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Compute hash dissimilarities for requested pairs")
    {
    	// Input
  	  addPositional ("pairs", "File with pairs of objects");
  	  addPositional ("hash_dir", "Directory with hashes for each object");
  	  addKey ("intersection_min", "Min. number of common hashes to compute distance", "50");
  	  addKey ("ratio_min", "Min. ratio of hash sizes (0..1)", "0.5");
  	  // Output
  	  addPositional ("out", "Output file with lines: <obj1> <obj2> <dissimlarity>; <obj1> < <obj2>");
  	}



	void body () const final
	{
		const string pairsFName       = getArg  ("pairs");
		const string hash_dir         = getArg  ("hash_dir");
		const size_t intersection_min = str2<size_t> (getArg ("intersection_min"));
		const Prob   hashes_ratio_min = str2real (getArg ("ratio_min"));
		const string out              = getArg  ("out");
		ASSERT (! hash_dir. empty ());
		ASSERT (isProb (hashes_ratio_min));
		ASSERT (! out. empty ());
		
		
    // Cf. hash2dissim.cpp
    OFStream output (out);
    ONumber on (output, 6, true);  // PAR
    map<string/*fName*/,Hashes> name2hashes;
    PairFile input (pairsFName);
    while (input. next ())
    {
      const string fName1 (hash_dir + "/" + input. name1);
      const string fName2 (hash_dir + "/" + input. name2);
      if (! contains (name2hashes, fName1))  name2hashes [fName1] = Hashes (fName1);
      if (! contains (name2hashes, fName2))  name2hashes [fName2] = Hashes (fName2);
      const Hashes& h1 = name2hashes [fName1];
      const Hashes& h2 = name2hashes [fName2];
      const double dissim = h1. getDissim (h2, intersection_min, hashes_ratio_min);
      output << input. name1 << '\t' << input. name2 << '\t' << dissim << endl;
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



