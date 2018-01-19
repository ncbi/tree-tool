// count.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "numeric.hpp"
using namespace DM_sp;



namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Statistics of a series of numbers from cin")
    {}



	void body () const final
	{
		MeanVar mv;
		for (;;)
		{
			Real x = NAN;
			cin >> x;
			if (cin. eof () && isNan (x))
				break;
			mv << x;
		}
		
		cout << "count" << '\t' << mv. n          << endl
		     << "mean"  << '\t' << mv. getMean () << endl
		     << "var"   << '\t' << mv. getVar ()  << endl
		     << "SD"    << '\t' << mv. getSD ()   << endl;
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



