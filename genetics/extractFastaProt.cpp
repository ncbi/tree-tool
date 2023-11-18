// extractFastaProt.cpp

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
*   Print protein sequences which are [not] in the target list
*
*/


#undef NDEBUG

#include "../common.hpp"
using namespace Common_sp;
#include "seq.hpp"
using namespace Seq_sp;
#include "../version.inc"

#include "../common.inc"



namespace 
{
	
	
struct Replacement : Named
{
	size_t from {0};
	size_t to {0};
	
	
	Replacement (const string &name_arg,
	             size_t from_arg,
	             size_t to_arg)
	  : Named (name_arg)
	  , from (from_arg)
	  , to (to_arg)
	  { 
	  	QC_ASSERT (goodName (name));
	  	QC_IMPLY (to, from < to);
	  }
	Replacement () = default;
	  
	
	size_t size () const
	  { return to - from; }
};
	



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print protein sequences which are [not] in the target list")
	  {
      version = VERSION;
		  addPositional ("in",  "Input FASTA file with proteins");
		  addPositional ("target", "File with line format: <sequence id> [<replacement id> [<from> <to>]] [...]");
		  addFlag ("remove", "Target list must be removed from the input file, otherwise only the target list is saved");
		  addFlag ("replace", "Replace sequence ids. There can be more than 1 replacement");
		  addFlag ("cut", "Cut out the segment indicated by <from> <to>");
		  addKey ("min_len", "Min. sequence length", "20");
	  }

  
  
	void body () const final
  {
	  const string inFName     = getArg ("in");
	  const string targetFName = getArg ("target");
	  const bool removeTarget  = getFlag ("remove");
	  const bool replaceP      = getFlag ("replace");
	  const bool cut           = getFlag ("cut"); 
	  const size_t len_min     = str2<size_t> (getArg ("min_len"));
	  
	  QC_IMPLY (removeTarget, ! replaceP);
	  QC_IMPLY (cut, replaceP);

  
    map<string/*id*/, Vector<Replacement>> name2replacement;
    {
      LineInput in (targetFName);  
    	Istringstream iss;
    	string name, replacement;
	  	while (in. nextLine ())
	  	{ 
	  	  trim (in. line);
	  	  if (! in. line. empty ())
	  	  {
	  	  	if (replaceP)
	  	    {
	  	    	iss. reset (in. line);
	  	    	iss >> name >> replacement;
	  	    	QC_ASSERT (! replacement. empty ());
	  	    	size_t from = 0;
	  	    	size_t to = 0;
	  	    	if (cut)
	  	    	{
	  	    		iss >> from >> to;
	  	    		if (from > to)
	  	    			throw runtime_error (name + ": from > to: " + in. line);
	  	    		QC_ASSERT (from);
	  	    		from--;
	  	    	}
	  	    	name2replacement [name] << std::move (Replacement (replacement, from, to));
	  	    }
	  	  	else
            name2replacement [in. line];  // = Vector<Replacement> ();
        }
	  	}      
    }


		{
		  Multifasta fa (inFName, true);
		  while (fa. next ())
		  {
		    const Peptide pep (fa, 1000/*PAR*/, true);
		    ASSERT (! strBlank (pep. seq));		
		    if ((! removeTarget) == contains (name2replacement, pep. getId ()))   
		    {
		    	if (replaceP)
		    	{
		    		const Vector<Replacement>& repls = name2replacement [pep. getId ()];
		    		ASSERT (! repls. empty ());
		    		for (const Replacement& repl : repls)
		    		{
		    		  Peptide pep1 (pep);
  		    		if (cut)
  		    		{
  		    			if (repl. to > pep. seq. size ())
  		    				throw runtime_error (pep. name + ": to = " + toString (repl. to) + ", but len = " + toString (pep. seq. size ()));
  		    			pep1. seq = pep. seq. substr (repl. from, repl. size ());
  		    		}
  		    		pep1. name = repl. name + " " + pep. name;
  		    		if (cut)
  		    		  pep1. name += ":" + toString (repl. from + 1) + "-" + toString (repl. to);
  		    	  if (pep1. seq. size () >= len_min)
  		          pep1. saveText (cout);
  		      }
		    	}
		    	else
		    	  if (pep. seq. size () >= len_min)
		          pep. saveText (cout);
		    }
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



