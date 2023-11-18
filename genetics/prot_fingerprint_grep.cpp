// prot_fingerprint_grep.cpp

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
*   Find proteins in file <triplet> which are substrings of proteins in <triplet_dir>
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

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Find proteins in file <triplet> which are substrings of proteins in <triplet_dir>")
    {
      version = VERSION;
  	  addPositional ("in", "Input FASTA file of proteins");
  	  addPositional ("triplet", "Triplet of AA");
  	  addPositional ("triplet_dir", "Directory with file <triplet> with proteins containing the triplet");
    }


	
	void body () const final
  {
	  const string inFName     = getArg ("in");
	  const string triplet     = getArg ("triplet");
	  const string triplet_dir = getArg ("triplet_dir");
	  ASSERT (triplet. size () == 3);
  

    const size_t window = 20;  // PAR


    VectorOwn<Peptide> needles;  needles. reserve (10000);  // PAR
    {
		  Multifasta trIn (triplet_dir + "/" + triplet, true);
		  while (trIn. next ())
		  {
		    Peptide* p = new Peptide (trIn, 1000, false);
		    ASSERT (p->seq. size () >= window);
		    needles << p;
		  }
		}
  
  
	  Multifasta fIn (inFName, true);
	  while (fIn. next ())
	  {
	    const Peptide p (fIn, 1000, false);
	    const string& seq = p. seq;
	    if (seq. size () < window)
	    	continue;
	    if (! contains (seq, triplet))
	    	continue;
	    for (const Peptide* needle : needles)
	      if (needle->name != p. name)
	      {
	      	const size_t pos = seq. find (needle->seq);
	      	if (pos != string::npos)
  	      	cout         << needle->name 
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



