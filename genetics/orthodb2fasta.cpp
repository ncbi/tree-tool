// orthodb2fasta.cpp

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
*   Extract FASTA file for given taxids from a gzipped OrthoDB file
*
*/


#undef NDEBUG

#include "../common.hpp"
#include "../gzip.hpp"
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
    : Application ("Extract FASTA file for given taxids from a gzipped OrthoDB file")
    {
      version = VERSION;
  	  addPositional ("orthodb", "OrthoDB gzipped FASTA file");
  	  addPositional ("taxids", "List of taxids");
    }


	
	void body () const final
  {
	  const string orthodbFName = getArg ("orthodb");
	  const string taxidsFName  = getArg ("taxids");


    StringVector taxids;  taxids. reserve (10000);  // PAR
    {
      LineInput f (taxidsFName);
      while (f. nextLine ())
      {
        trim (f. line);
        if (! f. line. empty ())
          taxids << std::move (f. line);
      }
    }
    taxids. sort ();
    taxids. uniq ();
    cerr << "# taxids: " << taxids. size () << '\n';

    GZip f (orthodbFName, 1000000);  // PAR
    bool headerP = true;
    string header;
    while (f. nextLine ())
    {
	    if (headerP)
	    {
	      header = std::move (f. line);
	      QC_ASSERT (! header. empty ());
	      QC_ASSERT (header [0] == '>');
	      replace (header, '\t', ' ');
	    }
	    else
	    {
	      QC_ASSERT (! f. line. empty ());
	      QC_ASSERT (f. line [0] != '>');
	      const Peptide p (header. substr (1), f. line, false);
	      p. qc ();
	      const size_t pos = p. name. find ('_');
	      QC_ASSERT (pos != string::npos);
	      QC_ASSERT (pos);
	      const string taxid (p. name. substr (0, pos));
	      if (taxids. containsFast (taxid))
	        p. saveText (cout);
	    }
	    toggle (headerP);
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



