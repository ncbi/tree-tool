// filterFasta.cpp

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
*   Filter FASTA file
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
	
	
struct Replacement final : Named
{
	size_t from {no_index};
	size_t to {no_index};
	
	
	Replacement () = default;
	explicit Replacement (const string &name_arg)
	  : Named (name_arg)
	  { 
	  	QC_ASSERT (goodName (name));
	  }
	Replacement (const string &name_arg,
	             size_t from_arg,
	             size_t to_arg)
	  : Named (name_arg)
	  , from (from_arg)
	  , to (to_arg)
	  { 
	  	QC_ASSERT (goodName (name));
	  	QC_ASSERT (from != no_index);
	  	QC_ASSERT (to != no_index);
	  	QC_ASSERT (from < to);
	  }	
	  
	
	size_t size () const
	  { return to - from; }
};
	



struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Print protein sequences which are [not] in the target list")
	  {
      version = VERSION;
		  addPositional ("in",  "Input FASTA file");
  	  addFlag ("aa", "FASTA file contains protein sequences, otherwise DNA sequences");
		  addFlag ("whole", "Sequence identifiers are whole strings which should not be split by '|'");
		  addKey ("target", "File with line format: <sequence id> [<new sequence id> [<from> <to>]]", "");
		  addFlag ("remove", "Sequences in <target> are removed from the input file, otherwise only the sequences in <target> are saved");
		  addFlag ("replace", "<target> file contains a replacement");
		  addFlag ("cut", "Cut out the segment indicated by <from> <to> in replacement in <target>");
		  addKey ("len_min", "Min. sequence length", "1");
		  addKey ("complexity_min", "Min. sequence complexity", "0");
		  addFlag ("title", "Preserve sequence titles");
		  addFlag ("force", "Do not stop on errors in input data");
	  }

  
  
	void body () const final
  {
	  const string inFName        = getArg ("in");
		const bool   aa             = getFlag ("aa");
	  const bool   whole          = getFlag ("whole");
	  const string targetFName    = getArg ("target");
	        bool   removeP        = getFlag ("remove");
	  const bool   replaceP       = getFlag ("replace");
	  const bool   cutP           = getFlag ("cut"); 
	  const size_t len_min        = str2<size_t> (getArg ("len_min"));
	  const double complexity_min = arg2double ("complexity_min");
	  const bool   titleP         = getFlag ("title");
	  const bool   forceP         = getFlag ("force");

	  
    if (targetFName. empty ())
    {
  	  if (removeP)
  	    throw runtime_error ("-remove requires -target");
      removeP = true;
    }

	  QC_IMPLY (removeP, ! replaceP);
	  QC_IMPLY (cutP, replaceP);

  
    map<string/*id*/, Vector<Replacement>> name2replacement;
    if (! targetFName. empty ())
    {
      LineInput in (targetFName);  
      Set<string> replacements;
    	Istringstream iss;
    	string name;
    	string replacement;
	  	while (in. nextLine ())
	  	{ 
	  	  trim (in. line);
	  	  if (in. line. empty ())
	  	    continue;
  	  	if (replaceP)
  	      try
  	      {
    	    	iss. reset (in. line);
    	    	replacement. clear ();
    	    	iss >> name >> replacement;
    	    	QC_ASSERT (! replacement. empty ());
    	    	if (! replacements. insert (replacement). second)
    	    	  throw runtime_error ("Duplicate replacement sequence id: " + strQuote (replacement));
    	    	size_t from = no_index;
    	    	size_t to   = no_index;
    	    	if (cutP)
    	    	{
    	    		iss >> from >> to;
    	    		QC_ASSERT (to != no_index);
    	    		if (from > to)
    	    			throw runtime_error (name + ": from > to: " + in. line);
    	    		QC_ASSERT (from);
    	    		QC_ASSERT (iss. eof ());
    	    		from--;
      	    	name2replacement [name] << std::move (Replacement (replacement, from, to));
    	    	}
    	    	else
      	    	name2replacement [name] << std::move (Replacement (replacement));
    	    }
    	    catch (const exception &e)
    	    {
    	      reportException (in. lineStr () + "\n" + in. line, e, forceP);
    	    }
  	  	else
          name2replacement [in. line];  // second.empty()
	  	}      
    }
    
    
    const auto process = [&] (Seq &seq)
      {
    	  if (seq. seq. size () < len_min)
    	    return;
        if (seq. getComplexity () < complexity_min)  
          return;
        if (! titleP)
          seq. name. erase (seq. getIdSize ());
        seq. saveText (cout);
      };


		{
		  Multifasta fa (inFName, aa);
		  Set<string> ids;
		  while (fa. next ())
		  {
		    unique_ptr<Seq> seq;
		    if (aa)
		    {
		      auto pep = new Peptide (fa, 1000/*PAR*/, true);
		      pep->pseudo = true;
		      seq. reset (pep);
		    }
		    else
		      seq. reset (new Dna (fa, 100000/*PAR*/, true));
		    ASSERT (seq);
		  //QC_ASSERT (! strBlank (seq->seq));	
		    seq->qc ();	
		    string id (seq->getId ());
		    const string id_ (whole ? id : findSplit (id, '|'));
	    	if (! ids. insert (id_). second)
	    	  throw runtime_error ("Duplicate sequence id: " + strQuote (id_));
		    if (removeP == contains (name2replacement, id_))   
		      continue;
	    	if (replaceP)
	    	{
	    		ASSERT (! name2replacement [id_]. empty ());
	    		for (const Replacement& repl : name2replacement [id_])
	    		  try
  	    		{
    	    		ASSERT (! repl. empty ());
  	    		  unique_ptr<Seq> newSeq (seq->copy ());
  	    		  ASSERT (newSeq);
    	    		if (cutP)
    	    		{
    	    			if (repl. to > seq->seq. size ())
    	    				throw runtime_error (seq->name + ": to = " + to_string (repl. to) + ", but len = " + to_string (seq->seq. size ()));
    	    			newSeq->seq = seq->seq. substr (repl. from, repl. size ());
    	    		}
    	    		newSeq->name = repl. name + " " + seq->name;
    	    		if (cutP)
    	    		  newSeq->name += ":" + to_string (repl. from + 1) + "-" + to_string (repl. to);
      	    	newSeq->qc ();
      	    	process (*newSeq);
    	    	}
      	    catch (const exception &e)
      	    {
      	      reportException (id_ + " -> " + repl. name, e, forceP);
      	    }    	    	
	    	}
	    	else
	    	  process (*seq);
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



