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

#include "common.hpp"
using namespace Common_sp;
#include "version.inc"

#include "common.inc"



namespace
{
  
  
Vector<int> goodStatuses;
unique_ptr<OFStream> errors;
std::mutex errorsMtx;



void executeCommand (const string &cmd,
                     const string &item/*,
                     uint blank_lines*/)
{
  ASSERT (! cmd. empty ());
  ASSERT (! item. empty ());
  
  const int c = system (cmd. c_str ());
    // "set -o pipefail && ..." does not work with dash
  const int exitStatus = WEXITSTATUS (c);
  if (exitStatus && ! goodStatuses. contains (exitStatus))
  {
  //cout << endl << exitStatus << endl;
    if (errors. get ())
    {
    //cout << c << endl;  // always 256 
      errorsMtx. lock ();
    	*errors << item << endl;
      QC_ASSERT (c != -1);
      errorsMtx. unlock ();
    }
    else
      throw runtime_error ("item=" + item + "  status=" + to_string (exitStatus) + "\n" + cmd);
  }
#if 0
  if (isMainThread ())
    FOR (uint, i, blank_lines)
      cerr << " " << endl;
#endif
}
  
  

struct Command
{
  string cmd;
  string item;
};
  
  
  
void executeCommands (size_t from,
                      size_t to,
                      Notype /*&res*/,
                      const Vector<Command> &commands,
                      uint step/*,
                      uint blank_lines*/)
{
  Progress prog (to - from, step);  
  FOR_START (size_t, i, from, to)
  {
    prog ();
    const Command& command = commands [i];
    executeCommand (command. cmd, command. item/*, blank_lines*/);
  }
}

  
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Apply <command> to all <items>")
  	{
      version = VERSION;
  	  addPositional ("items", "Items container: file with items (end-of-line separated), a directory (in this case items are files in this directory in raw order), or a natural number");
  	  addPositional ("command", "Text with special symbols: \
\"%d\" = <items>, \
\"%f\" - item, \
\"%<n>\" - subitem #n (1 <= n <= 9), \
\"%{<n>}\" - subitem #n (n >= 10), \
\"%h\" - item hash (0.." + to_string (hash_class_max - 1) + "), \
\"%q\" - single quote, \
\"%Q\" - double quote, \
\"%D\" - $, \
\"%g\" - `, \
\"%b\" - backslash, \
\"%p\" - %");
      addFlag ("large", "Directory <items> is large: it is subdivided into subdirectories \"0\" .. \"" + to_string (hash_class_max - 1) + "\" which are the hashes of file names");
  	  addKey ("errors", "Ignore errors in running items and save error items into this file");
  	    // Bug: ^C does not stop the program ??
  	//addKey ("blank_lines", "# Blank lines to be printed to stderr after each command", "0");
  	  addKey ("step", "# Items processed to output the progress for", "100");
  	  addKey ("start", "# Item to start with", "1");
  	  addFlag ("zero", "Item numbers are 0-based, otherwise 1-based");
  	  addFlag ("print", "Print command, not execute");
  	  addFlag ("tsv", "<items> is a tsv-file");
  	  addKey ("good_exit_codes", "Comma-delimited list of allowed exit codes", "0");
  	}
  	
  	
 
	void body () const final
	{
		const string itemsName     = getArg ("items");
		const string cmd_orig      = getArg ("command");
		const bool   large         = getFlag ("large");
		const string errorsFName   = getArg ("errors");
//const uint blank_lines       = str2<uint> (getArg ("blank_lines"));
		const uint step            = str2<uint> (getArg ("step"));
		const uint start           = str2<uint> (getArg ("start"));
		const bool zero            = getFlag ("zero");
		const bool printP          = getFlag ("print");
		const bool tsv             = getFlag ("tsv");
		      string goodStatusesS = getArg ("good_exit_codes");

    if (itemsName. empty ())
      throw runtime_error ("Empty items list name");
    if (cmd_orig. empty ())
      throw runtime_error ("Empty command");    
    if (! step)
      throw runtime_error ("-step must be >= 1");
  //if (blank_lines && step > 1)
    //throw runtime_error ("-blank_lines requires -step to be 1");


	  Vector<Command> commands;  commands. reserve (100000);  // PAR
    {
  	  unique_ptr<ItemGenerator> gen;
  	  {
  	    const uint stepItemGen = threads_max > 1 ? 1000 : step;  // PAR
    	  const bool isFile = fileExists (itemsName);
    	  const bool isDir = 
    	    #ifdef _MSC_VER
    	      false
    	    #else
    	      getFiletype (itemsName, true) == Filetype::dir
    	    #endif
    	    ;
    	  if (isDir && tsv)
    	    throw runtime_error ("-tsv cannot be used if <items> is a directory");
    	  if (isFile)
    	    gen. reset (new FileItemGenerator (stepItemGen, itemsName, tsv));
    	  else if (isDir)
    	    gen. reset (new DirItemGenerator (stepItemGen, itemsName, large));
        else 
        {
          size_t n = 0; 
          if (str2<size_t> (itemsName, n))
          {
        	  if (tsv)
        	    throw runtime_error ("-tsv cannot be used if <items> is a number");
            gen. reset (new NumberItemGenerator (stepItemGen, n));	  
          }
          else
            throw runtime_error ("File " + strQuote (itemsName) + " does not exist");
        }
      }
      ASSERT (gen. get ());
  	
  	
      ASSERT (! errors. get ());
      if (! errorsFName. empty ())
      	errors. reset (new OFStream (errorsFName));
  	

 	    constexpr char delChar = '\177';

	    if (contains (cmd_orig, delChar))
	      throw runtime_error ("Command has an ASCII 127 (DEL) character");

	    if (contains (itemsName, delChar))
	      throw runtime_error ("Items container has an ASCII 127 (DEL) character");
	    string itemsName_cleaned (itemsName);  // Can contain '%'
	    replace (itemsName_cleaned, '%', delChar);


  	  string cmd (cmd_orig);
  	  // Constant replacements of '%'
  	  replaceStr (cmd, "%d", itemsName_cleaned);  
      replaceStr (cmd, "%q", "'");  
      replaceStr (cmd, "%Q", "\"");  
      replaceStr (cmd, "%D", "$");  
      replaceStr (cmd, "%g", "`");
      replaceStr (cmd, "%b", "\\");
  	
      bool subitemsP = false;
      FFOR (size_t, i, cmd. size () - 1)
        if (   cmd [i] == '%'
            && ((isDigit (cmd [i + 1]) && cmd [i + 1] != '0') || cmd [i + 1] == '{')
           )
        {
          subitemsP = true;
          break;
        }
        
      commaize (goodStatusesS);
      while (! goodStatusesS. empty ())
      {
        const string s (findSplit (goodStatusesS, ','));
        goodStatuses << str2<int> (s);
      }
    #if 0
      for (const int i : goodStatuses)
        cout << i << endl;  
    #endif

  	  string item;
  	  Istringstream iss;
      StringVector subitems;  
  	  while (gen->next (item))
      {
        if (! tsv)
          trim (item);
        if (item. empty ())
          continue;
        
        ASSERT (gen->prog. n);
        const size_t n = gen->prog. n - (zero ? 1 : 0);
        if (n < start)
          continue;
          
        const string item_orig (item);
        const size_t hash_class = str2hash_class (item_orig);
  	      
  	    if (contains (item, delChar))
  	      throw runtime_error ("Item has an ASCII 127 (DEL) character");
  	    replace (item, '%', delChar);
  	      
        string thisCmd (cmd);
        replaceStr (thisCmd, "%f", item);
        replaceStr (thisCmd, "%h", to_string (hash_class));
        replaceStr (thisCmd, "%n", to_string (n));

        if (subitemsP)
        {
          subitems. clear ();
          if (tsv)
          {
          #ifndef NDEBUG
            const size_t colNum = strCountSet (item, "\t") + 1;
          #endif
            string s (item);
            const bool lastEmpty = isRight (item, "\t");
            while (! s. empty ())
              subitems << findSplit (s, '\t');
            if (lastEmpty)
              subitems << string ();
            ASSERT (subitems. size () == colNum);
          }
          else
          {
            iss. reset (item);
            string s;
            while (! iss. eof ())
            {
              ASSERT (s. empty ());
              iss >> s;
              ASSERT (! s. empty ());
              subitems << std::move (s);
            }
          }
          ASSERT (! subitems. empty ());
          FFOR (size_t, i, subitems. size ())
          {
            const string numS (to_string (i + 1));
            const bool longP = numS. size () > 1;
            replaceStr (thisCmd, string ("%") + (longP ? "{" : "") + numS + (longP ? "}" : ""), subitems [i]); 
          }
        }        

        // Last replacement of '%'
        FFOR (size_t, i, thisCmd. size ())
          if (   thisCmd [i] == '%' 
              && (   i == thisCmd. size () - 1 
                  || thisCmd [i + 1] != 'p'
                 )
             )
            throw runtime_error ("Unprocessed " + strQuote (thisCmd. substr (i, 2)) + " in item: " + item_orig + "\n" + thisCmd);
        replaceStr (thisCmd, "%p", "%");
        
        // Restore original '%'
  	    replace (thisCmd, delChar, '%');

        if (verbose ())
       	  cerr << thisCmd << endl;
       	  
        if (printP)
          cout << thisCmd << endl;
        else
          if (threads_max == 1)
            executeCommand (thisCmd, item_orig/*, blank_lines*/);
      	  else
      	    commands << std::move (Command {std::move (thisCmd), std::move (item_orig)});
      }
    }


    if (threads_max > 1)
    {
      commands. randomOrder ();
      vector<Notype> notypes;
  	  arrayThreads (false, executeCommands, commands. size (), notypes, cref (commands), step/*, blank_lines*/);
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



