// extractPairs.cpp

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;



namespace 
{
  
  
List<string> getList (string &s,
                      const string &delim)
{
  trim (s);
  replaceStr (s, delim, " ");
  replaceStr (s, "  ", " ");
  return str2list (s);  
}

  
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("A1<delim>A2<delim>A3...\\tB1<delim>B2<delim>B3... -> A1\\tB1\\nA1\\tB2\\nA1\\tB3\\nA2\\tB1\\nA2\\tB2\\nA2\\tB3\\nA3\\tB1\\nA3\\tB2\\nA3\\tB3...")
  	{
  	  addPositional ("in", "Input file");
  	  addKey ("delim", "Delimiter", " ");
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
  	  const List<string> lhsVec (getList (lhs, delim));
  	  const List<string> rhsVec (getList (f. line, delim));
  	  for (const string& a : lhsVec)
    	  for (const string& b : rhsVec)
    	    cout << a << '\t' << b << endl;
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
