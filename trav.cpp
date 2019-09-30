// trav.cpp

/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE                          
*               National Center for Biotechnology Information
*                                                                          
*  This software/database is a "United States Government Work" under the   
*  terms of the United States Copyright Act.  It was written as part of    
*  the author's official duties as a United States Government employee and 
*  thus cannot be copyrighted.  This software/database is freely available 
*  to the public for use. The National Library of Medicine and the U.S.    
*  Government have not placed any restriction on its use or reproduction.  
*                                                                          
*  Although all reasonable efforts have been taken to ensure the accuracy  
*  and reliability of the software and data, the NLM and the U.S.          
*  Government do not and cannot warrant the performance or results that    
*  may be obtained by using this software or data. The NLM and the U.S.    
*  Government disclaim all warranties, express or implied, including       
*  warranties of performance, merchantability or fitness for any particular
*  purpose.                                                                
*                                                                          
*  Please cite the author in any work or product based on this material.   
*
* ===========================================================================
*
* Author: Vyacheslav Brover
*
* File Description:
*   Traverse a list of objects and apply a command to each of the objects
*
*/

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
\"%h\" - item hash (0.." + to_string (hash_class_max - 1) + "), \
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
  	  addFlag ("print", "Print command, not execute");
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
		const bool printP        = getFlag ("print");
    ASSERT (! items. empty ());


    // ^C to stop trav if noerror ??
    

	  unique_ptr<ItemGenerator> gen;
	  {
  	  const bool isFile = fileExists (items);
  	  const bool isDir = directoryExists (items);
  	  if (isFile || isDir)
  	    gen. reset (new FileItemGenerator (step, isDir, items));
      else 
        if (isDigit (items [0]))
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
      replaceStr (thisCmd, "%h", to_string (str2hash_class (item)));
      replaceStr (thisCmd, "%n", to_string (gen->prog. n));  
      if (contains (thisCmd, "%"))
        throw runtime_error ("Unprocessed '%' in item=" + item + "\n" + thisCmd);
      if (verbose ())
      	cerr << thisCmd << endl;
      if (printP)
        cout << thisCmd << endl;
      else
      {
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
            throw runtime_error ("item=" + item + "  status=" + to_string (WEXITSTATUS (exitStatus)) + "\n" + thisCmd);
        }
  	    FOR (uint, i, blank_lines)
  	      cerr << " " << endl;
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


