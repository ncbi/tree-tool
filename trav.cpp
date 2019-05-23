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
    : Application ("Apply <command> to all <items>")
  	{
  	  addPositional ("items", "File with items (end-of-line separated), a directory (in this case items are files in this directory), or a natural number");
  	  addPositional ("command", "Text with special symbols: \
\"%d\" = <items>, \
\"%f\" - an item, \
\"%n\" - sequential number, \
\"%q\" - single quote, \
\"%Q\" - double quote, \
\"%D\" - $, \
\"%G\" - `");
  	  addKey ("errors", "Ignore errors in running items and save error items into this file");
  	  addFlag ("quote", "Quote %f");
  	  addKey ("blank_lines", "# Blank lines to be printed on the screen after each command", "0");
  	  addKey ("step", "# Items processed to output the progress for", "100");
  	  addKey ("start", "# Item to start with", "1");
  	}
  


	void body () const final
	{
		      string items       = getArg ("items");
		const string cmd_        = getArg ("command");
		const string errorsFName = getArg ("errors");
		const bool quote         = getFlag ("quote");
		const uint blank_lines   = str2<uint> (getArg ("blank_lines"));
		const uint step          = str2<uint> (getArg ("step"));
		const uint start         = str2<uint> (getArg ("start"));
    ASSERT (! items. empty ());


    // ^C to stop trav if noerror ??
    

	  unique_ptr<ItemGenerator> gen;
	  {
  	  const bool isFile = fileExists (items);
  	  const bool isDir = directoryExists (items);
  	  if (isFile || isDir)
  	    gen. reset (new FileItemGenerator (step, isDir, items));
      else 
        if (isdigit (items [0]))
          gen. reset (new NumberItemGenerator (step, items));	  
        else
          throw runtime_error ("File " + strQuote (items) + " does not exist");
    }
    ASSERT (gen. get ());
	
	
	  string cmd (cmd_);
	  replaceStr (cmd, "%d", items);
    replaceStr (cmd, "%q", "'");  
    replaceStr (cmd, "%Q", "\"");  
    replaceStr (cmd, "%D", "$");  
    replaceStr (cmd, "%G", "`");
	

    unique_ptr<OFStream> errors;
    if (! errorsFName. empty ())
    	errors. reset (new OFStream (errorsFName));
	
	  string item;
	  while (gen->next (item))
    {
      if (gen->prog. n < start)
        continue;
	      
      // Preparing item for using it in a shell command
      FOR_REV (size_t, i, item. size ())
        if (item [i] == '\\')
        	item. replace (i, 0, "\\");
      if (quote)
        item = "\"" + item + "\"";

      string thisCmd (cmd);
      replaceStr (thisCmd, "%f", item);
      replaceStr (thisCmd, "%n", toString (gen->prog. n));  
      if (contains (thisCmd, "%"))
        throw runtime_error ("Unprocessed '%' in item=" + item + "\n" + thisCmd);
      if (verbose ())
      	cerr << thisCmd << endl;
      const int exitStatus = system (thisCmd. c_str ());
      if (exitStatus)
      {
        if (errors. get ())
        {
        	*errors << item << endl;
	        if (exitStatus == -1)  // ??
	          ERROR;
        }
        else
          throw runtime_error ("item=" + item + "  status=" + toString (WEXITSTATUS (exitStatus)) + "\n" + thisCmd);
      }

	    FOR (uint, i, blank_lines)
	      cerr << " " << endl;
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


