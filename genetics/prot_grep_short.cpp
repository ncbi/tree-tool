// prot_grep_short.cpp

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
*   "Find proteins in <short> which are substrings of proteins in <in>
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "seq.hpp"
using namespace Seq_sp;
#include "../version.inc"



namespace 
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Find proteins in <short> which are substrings of proteins in <in>")
    {
      version = VERSION;
  	  addPositional ("in", "FASTA file with all proteins");
  	  addPositional ("short", "FASTA file with short proteins");
    }


	
	void body () const final
  {
	  const string inFName    = getArg ("in");
	  const string shortFName = getArg ("short");
  

    VectorOwn<Peptide> needles;  needles. reserve (10000);  // PAR
    {
		  Multifasta trIn (shortFName, true);
		  while (trIn. next ())
		    needles << new Peptide (trIn, 1000, false);
		}
  
  
	  Multifasta fIn (inFName, true);
	  while (fIn. next ())
	  {
	    const Peptide p (fIn, 1000, false);
	    const string& seq = p. seq;
	    for (const Peptide* pep : needles)
	      if (pep->name != p. name)
	      {
	      	// Use fast substring search ??
	      	const size_t pos = seq. find (pep->seq);
	      	if (pos != string::npos)
  	      	cout         << pep->name 
  	      	         << '\t' << p. name 
  	      	         << '\t' << pos + 1
  	      	         << endl;
	      }
	  }
  }
};



}  // namespace



int main(int argc, 
         const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



