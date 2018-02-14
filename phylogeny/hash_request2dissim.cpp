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
  	  addKey ("hashes_min", "Min. number of common hashes to compute distance", "50");
  	  // Output
  	  addPositional ("out", "Output file");
  	}



	void body () const final
	{
		const string pairsFName = getArg  ("pairs");
		const string hash_dir   = getArg  ("hash_dir");
		const size_t hashes_min = str2<size_t> (getArg ("hashes_min"));
		const string out        = getArg  ("out");
		ASSERT (! hash_dir. empty ());
		ASSERT (! out. empty ());
		
		
    // Cf. hash2dissim.cpp
    LineInput input (pairsFName, 100 * 1024, 1000);  // PAR
    OFStream output (out);
    ONumber on (output, 6, true);  // PAR
    map<string/*fName*/,Hashes> name2hashes;
    while (input. nextLine ())
    {
      istringstream iss (input. line);
      string name1, name2;
      iss >> name1 >> name2;
      ASSERT (name1 != name2);
      const string fName1 (hash_dir + "/" + name1);
      const string fName2 (hash_dir + "/" + name2);
      if (! contains (name2hashes, fName1))  name2hashes [fName1] = Hashes (fName1);
      if (! contains (name2hashes, fName2))  name2hashes [fName2] = Hashes (fName2);
      const Hashes& h1 = name2hashes [fName1];
      const Hashes& h2 = name2hashes [fName2];
      const double dissim = h1. getDissim (h2, hashes_min);
      output << name1 << '\t' << name2 << '\t' << dissim;
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



