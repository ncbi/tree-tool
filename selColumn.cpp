// selColumn.cpp

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;



namespace 
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Select <Column> out of <File>\nSpaces and tabs are automatically included in <delimiters>")
  	{
  	  addPositional ("in", "Text file name or '-'");
  	  addPositional ("col", "Column number (>=1)");
  	  addKey ("delimiters", "Delimiters");
  	}



	void body () const final
  {
	  const string fName   = getArg ("in");
	  const uint targetCol = str2<uint> (getArg ("col"));
	  string delimiters    = getArg ("delimiters");	  
	  ASSERT (targetCol > 0);
	  
	  
	  delimiters += " \t";

    istream* is = & cin;
    Common_sp::AutoPtr <ifstream> ifs;
    if (fName != "-")
    {
      ifs = new ifstream (fName. c_str ());  
      is = ifs. get ();
    }
    uint row = 0;
	  for (;;)
	  {
	  	row++;
		  string s; 
		  FOR (uint, col, targetCol)
		  {
		  	s = getToken (*is, delimiters, "\n");
		  	if (s == "\n")
		  	{
		  		cout << "No column " << targetCol + 1 << " in row " << row << endl;
		  		ERROR;
		  	}
		  }
		  skipLine (*is);  
		  if (s. empty ())
		  	break;
		  cout << s << endl;  
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


