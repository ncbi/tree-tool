// trav.cpp

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;



namespace
{
  
  
struct ItemGenerator
{
  virtual uint steps () const 
    { return 0; }
  virtual bool next (string &item) = 0;
    // Return: false <=> end of items
    // Output: item; may be empty()
};



struct FileItemGenerator : ItemGenerator
{
private:
  const string fName;
  ifstream f;
  const bool temp;
public:
  
  FileItemGenerator (const string& fName_arg,
                     bool temp_arg)
    : fName( fName_arg)
    , f (fName_arg)
    , temp (temp_arg)
    { ASSERT (f. good ()); }
 ~FileItemGenerator ()
    { if (temp)
	      remove (fName. c_str ());
	  }
  
  bool next (string &item) final
    { if (f. eof ())
        return false;
    	readLine (f, item);
      if (temp)
      { const size_t pos = item. rfind ('/');
      	if (pos != string::npos)
          item. erase (0, pos + 1);
      }
      trim (item);
    	return true;
    }    
};

  

struct NumberItemGenerator : ItemGenerator
{
private:
  const uint n;
  uint i {0};
public:
  
  explicit NumberItemGenerator (const string& name)
    : n (str2<uint> (name))
    {}
  
  uint steps () const final
    { return n; }
  bool next (string &item) final
    { if (i == n)
        return false;
      i++;
      item = toString (i);
      return true;
    }
};

  

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Apply <command> to all items")
  	{
  	  addPositional ("items", "File with items (end-of-line separated), a directory (in this case items are files in this directory), or a natural number");
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
		const uint blank_lines = str2<uint> (getArg ("blank_lines"));
		const uint step        = str2<uint> (getArg ("step"));
		const uint start       = str2<uint> (getArg ("start"));
    ASSERT (! items. empty ());


    // ^C to stop trav if noerror ??
    

	  Common_sp::AutoPtr<ItemGenerator> gen;
	  if (fileExists (items))
	    gen. reset (new FileItemGenerator (items, false));
	  else
    {
	    if (items. at (items. size () - 1) == '/')
	    	items. erase (items. size () - 1);
      if (directoryExists (items))
      {
    	  char lsfName [4096] = {'\0'};
        strcpy (lsfName, P_tmpdir);
        strcat (lsfName, "/XXXXXX");
        EXEC_ASSERT (mkstemp (lsfName) != -1);
        ASSERT (lsfName [0]);
  
        const int res = system (("ls " + items + " > " + lsfName). c_str ());
      //printf ("res = %d\n", res);
        ASSERT (! res);
        
	      gen. reset (new FileItemGenerator (lsfName, true));
      }
      else
	      gen. reset (new NumberItemGenerator (items));
    }
    ASSERT (gen. get ());
	
	
	  string cmd (cmd_);
	  replaceStr (cmd, "%d", items);
	
	
	  Progress prog (gen->steps (), step);
	  string item;
	  while (gen->next (item))
    {
      if (item. empty ())
        continue;
	    FOR (uint, i, blank_lines)
	      cout << " " << endl;
      prog (item);
      if (prog. n < start)
        continue;

      // Preparing item for using it in a shell command
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


