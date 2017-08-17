// normalCdf.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
using namespace DM_sp;



namespace
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Normal distribution CDF")
  	{
  	  addPositional ("loc",   "Location parameter");	
  	  addPositional ("scale", "Scale parameter");	  
  	  addPositional ("x",     "Value");
  	}



	void body () const final
	{
		const Real loc   = str2real (getArg ("loc"));
		const Real scale = str2real (getArg ("scale"));
		const Real x     = str2real (getArg ("x"));

    Normal norm;
    norm. setParam (loc, scale);
    cout << norm. cdf (x) << endl;
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


