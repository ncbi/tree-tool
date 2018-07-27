// count.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "numeric.hpp"
using namespace DM_sp;




namespace 
{
	
	
template <typename T>
	void report (const string &attr,
	             T value)
	  {	cout << attr << '\t' << value << endl; }

	


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Statistics of a sequence of numbers from cin", false)
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
		
		report ("count",   mv. n);
		report ("mean",    mv. getMean ());
		report ("var",     mv. getVar ());
		report ("SD",      mv. getSD ());
		report ("sum",     mv. s);
		report ("mean SD", mv. getSD () / sqrt (mv. n));
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



