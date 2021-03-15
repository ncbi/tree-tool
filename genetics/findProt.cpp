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
*   Print a protein sequence closest to a target
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
    : Application ("Print protein sequence closest to a target")
	  {
      version = VERSION;
		  addPositional ("in",  "Multi-FASTA file with proteins");
		  addPositional ("target", "FASTA file with a target protein");
	  }

  
  
	void body () const final
  {
	  const string inFName     = getArg ("in");
	  const string targetFName = getArg ("target");


    constexpr size_t k = 10;  // PAR
    
    unordered_map<string,size_t/*pos*/> kmer2pos;  kmer2pos. rehash (100000);  // PAR
    size_t len = 0;
    {
      Multifasta targetMF (targetFName, true);
      const Peptide target (targetMF, 1000, false);  // PAR
      if (target. seq. size () < k)
        return;
      len = target. seq. size () - k;
      FFOR (size_t, i, len)
        kmer2pos [target. seq. substr (i, k)] = i;
    }
    ASSERT (len);
    
	  Peptide best;
	  size_t n_max = 0;
    {
  	  Multifasta fa (inFName, true, 1024 * 1024);  // PAR ??
      Vector<char/*bool*/> found (len, 0);
  	  while (fa. next ())
  	  {
  	    const Peptide pep (fa, 1000/*PAR*/, false);
  	    if (pep. seq. size () < k)
  	      continue;
  	    for (char& b : found)
  	      b = 0;
        FFOR (size_t, i, pep. seq. size () - k)
          if (const size_t* pos = findPtr (kmer2pos, pep. seq. substr (i, k)))
  	        found [*pos] = 1;
  	    size_t n = 0;
  	    for (const char b : found)
  	      if (b)
  	        n++;
  	    if (maximize (n_max, n))
  	      best = pep;
  	  }
	  }
	  if (! best. seq. empty ())
	  {
  	  cout << "Matches: " << n_max << " / " << len << endl;
	    best. saveText (cout);
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



