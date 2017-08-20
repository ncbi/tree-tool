// trav.cpp

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;



namespace
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Apply <command> to all items")
  	{
  	  // Iterations n times: use %n only ??
  	  addPositional ("items", "File with items (end-of-line separated) or a directory (in this case items are files in this directory)");
  	  addPositional ("command", "Text with special symbols: \"%d\" = <items>, \"%f\" - an item, \"%n\" - sequential number");
  	  addFlag ("noerror", "Ignore errors");
  	  addFlag ("quote", "Quote %f");
  	  addKey ("blank_lines", "# Blank lines to be printed on the screen after each command", "0");
  	  addKey ("step", "# Items processed to output the progress for", "1000");
  	  addKey ("start", "# Item to start with", "1");
  	}
  


	void body () const
	{
		      string items     = getArg ("items");
		const string cmd_      = getArg ("command");
		const bool ignoreError = getFlag ("noerror");
		const bool quote       = getFlag ("quote");
		const int blank_lines  = str2<int> (getArg ("blank_lines"));
		const uint step        = str2<uint> (getArg ("step"));
		const uint start       = str2<uint> (getArg ("start"));
    ASSERT (! items. empty ());


    // ^C to stop trav if noerror ??

	  // LSFName
	  char lsfName [4096];
	  bool temp = false;
	  if (fileExists (items))
      strcpy (lsfName, items. c_str ());
	  else
    {
	    if (items. at (items. size () - 1) == '/')
	    	items. erase (items. size () - 1);

      ASSERT (directoryExists (items));

    //tmpnam (lsfName);
      strcpy (lsfName, P_tmpdir);
      strcat (lsfName, "/XXXXXX");
      EXEC_ASSERT (mkstemp (lsfName) != -1);
      ASSERT (lsfName [0]);

      const int res = system (("ls " + items + " > " + lsfName). c_str ());
    //printf ("res = %d\n", res);
      ASSERT (! res);
      
      temp = true;
    }
	
	
	  string cmd (cmd_);
	  replaceStr (cmd, "%d", items);
	
	
	  ifstream f (lsfName);
	  ASSERT (f. good ());
	  string item;
	  Progress prog (0, step);
	  while (! f. eof ())
    {
  	  readLine (f, item);
      if (item. empty ())
        continue;
      prog (item);
      if (prog. n < start)
        continue;

      if (temp)
      { 
        const size_t pos = item. rfind ('/');
      	if (pos != string::npos)
          item. erase (0, pos + 1);
      }
      ASSERT (! item. empty ());
      trim (item);

      FOR_REV (size_t, i, item. size ())
        if (item. at (i) == '\\')
        	item. replace (i, 0, "\\");
      if (quote)
        item = "\"" + item + "\"";

      string thisCmd (cmd);
      replaceStr (thisCmd, "%f", item);
      replaceStr (thisCmd, "%n", toString (prog. n));  
      if (verbose ())
      	cerr << thisCmd << endl;
      const int exitStatus = system (thisCmd. c_str ());
      if (exitStatus)
      {
        if (! ignoreError)
        {
        	cout << endl 
        	     << "item=" << item 
        	     << "  status=" << exitStatus 
        	     << endl;
          cout << thisCmd << endl;
        }
        if (exitStatus == -1 || ! ignoreError)
          ERROR;
      }

	    FOR (int, i, blank_lines)
	      cout << " " << endl;
    }
	
	  if (temp)
	    remove (lsfName);
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


