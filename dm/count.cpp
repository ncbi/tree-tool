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
    : Application ("Statistics of a series of numbers from cin", false)
    {}



	void body () const final
	{
		MeanVar mv;
		string s;
		while (cin >> s)
		{
			const Real x = str2real (s);
			if (! isNan (x))
				mv << x;
		}
		
		cout << "count" << '\t' << mv. n          << endl
		     << "mean"  << '\t' << mv. getMean () << endl
		     << "var"   << '\t' << mv. getVar ()  << endl
		     << "SD"    << '\t' << mv. getSD ()   << endl
		     << "sum"   << '\t' << mv. s          << endl;
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



