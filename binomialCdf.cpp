// binomialCdf.cpp

#undef NDEBUG

#include "common.inc"

#include "common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
using namespace DM_sp;



namespace
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Binomial distribution CDF")
  	{
  	  addPositional ("N", "Number of coin tosses");	
  	  addPositional ("P", "Probability of heads");	
  	  addPositional ("M", "Number of heads");
  	  addFlag("reverse", "1 - CDF");
  	  addFlag("add_pmf", "Add p.m.f.");
  	}



	void body () const final
	{
		const uint n       = str2<uint> (getArg ("N"));
		const Prob p       = str2real (getArg ("P"));
		const uint m       = str2<uint> (getArg ("M"));
		const bool reverse = getFlag ("reverse");
		const bool add_pmf = getFlag ("add_pmf");


    Binomial bin;
    bin. setParam ((int) n, p);
    
    Prob res;
    if (reverse)
      res = 1 - bin. cdf (m - 1);
    else
      res = bin. cdf (m);

    cout << res << endl;
    
    if (add_pmf)
      cout << "pmf = " << bin. pdf (m) << endl;
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


