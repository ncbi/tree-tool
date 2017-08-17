// pairs2Boolean.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
using namespace DM_sp;



namespace 
{


void parseLine (string line,
                string &objName,
                string &attrName)
{
//cout << line << endl;  
  const string line_orig (line);
  try 
  {
	  replace (line, '\t', ' ');
	  replaceStr (line, "  ", " ");
	  trim (line);
	  objName  = findSplit (line);
	  attrName = findSplit (line);
	}
	catch (const exception &e)
	{
		throw runtime_error ("Line: " + line_orig + "\n" + e. what ());
	}
}




struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Convert Boolean attributes to a " + dmSuff + "-file")
    {
  	  addPositional ("pairs", "File with lines: <obj> <attr>");
  	}



	void body () const
	{
		const string pairsFName = getArg ("pairs");

		
    Set<string> objNames;
    Set<string> attrNames;
    { 
      cerr << "Pass 1 ..." << endl;
      LineInput f (pairsFName, 10000000, 1000000);  // PAR
      while (f. nextLine ())
      {
        string objName, attrName;
        parseLine (f. line, objName, attrName);
        objNames  << objName;
        attrNames << attrName;
      }
    }
    cerr << "# Objects: " << objNames. size () << endl;
    cerr << "# Boolean attributes: " << attrNames. size () << endl;
    ASSERT (! objNames. empty ());
    ASSERT (! attrNames. empty ());
    
    Dataset ds;
    for (const string& name : objNames)
      ds. appendObj (name);
    ds. setName2objNum ();

    map<string,CompactBoolAttr1*> name2attr;
    for (const string& name : attrNames)
    {
      auto attr = new CompactBoolAttr1 (name, ds);
      name2attr [name] = attr;
      attr->setAll (false);
    }
    
    {
      cerr << "Pass 2 ..." << endl;
      LineInput f (pairsFName, 10000000, 1000000);  // PAR
      while (f. nextLine ())
      {
        string objName, attrName;
        parseLine (f. line, objName, attrName);
        const size_t row = ds. getName2objNum (objName);
        CompactBoolAttr1* attr = name2attr [attrName];
        ASSERT (attr);
        attr->setCompactBool (row, true);
      }
    }
        
    ds. qc ();
    ds. saveText (cout);    
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



