// fasta2len.cpp

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
*   Print lengths of DNA or protein sequences
*
*/


#undef NDEBUG

#include "../common.hpp"
#include "../tsv/tsv.hpp"
using namespace Common_sp;
#include "seq.hpp"
using namespace Seq_sp;
#include "../version.inc"

#include "../common.inc"



namespace 
{
	
	
struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Print lengths of DNA or protein sequences")
    {
      version = VERSION;
  	  addPositional ("in", "FASTA file");
  	  addFlag ("aa", "Protein sequence, otherwise DNA");
  	  addKey ("min_len", "Min. length for output sequences in the file <out>", "0");
  	  addKey ("out", "Output FASTA file with sequences longer than <min_len>");
  	  addKey ("header", "Add comma-separated .tsv-header to the output file");
    }


	
	void body () const final
  {
	  const string inFName  = getArg ("in");
	  const bool aa         = getFlag ("aa");
	  const size_t len_min  = str2<size_t> (getArg ("min_len"));
	  const string outFName = getArg ("out");
	  const string headerS  = getArg ("header");

	  
    Vector<TextTable::Header> header (TextTable::str2header (nvl (headerS, "id,len")));
    if (header. size () != 2)
      throw runtime_error ("Output file should have 2 columns");
    header. back (). numeric = true;

    TextTable tab (true, header);
    tab. saveHeader = ! headerS. empty ();
    {
  	  unique_ptr<OFStream> outF;
  	  if (! outFName. empty ())
  	  	outF. reset (new OFStream (outFName));
  	  Multifasta faIn (inFName, aa);
  	  while (faIn. next ())
  	  {
  	  	unique_ptr<const Seq> seq;
  	  	if (aa)
  	  		seq. reset (new Peptide (faIn, 1024 * 1024, false));
  	  	else
  	  		seq. reset (new Dna     (faIn, 1024 * 1024, false));
  	    seq->qc ();
  	    
  	    StringVector row;
  	    row << seq->getId () << to_string (seq->seq. size ());
  	    tab. rows << std::move (row);
  	      
  	    if (outF && seq->seq. size () >= len_min)
  	    	seq->saveText (*outF);
  	  }
  	}
  	tab. saveText (cout);
  }
};


}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



