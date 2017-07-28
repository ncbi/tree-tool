// extractPairs.cpp

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;



namespace 
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("A\\tB1<delim>B2<delim>B3... -> A\\tB1\\nA\\tB2\\nA\\tB3\\n...")
  	{
  	  addPositional("in", "Input file");
  	  addPositional("delim", "Delimiter");
  	}



	void body () const final
	{
		const string fName = getArg ("in");
		const string delim = getArg ("delim");
		ASSERT (! delim. empty ());

    LineInput f (fName, 100 * 1024, 1000);
  	while (f. nextLine ())
  	{ 
  	  trim (f. line);
  	  if (f. line. empty ())
  	    continue;
  	  string lhs (findSplit (f. line, '\t'));
  	  trim (lhs);
  	  replaceStr (f. line, delim, "\t");
  	  while (! f. line. empty ())
  	  {
  	    string rhs (findSplit (f. line, '\t'));
  	    trim (rhs);
  	    cout << lhs << '\t' << rhs << endl;
  	  }
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
