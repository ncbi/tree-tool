// prot_check.cpp

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
*   Check proteins
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
	
	
struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Check proteins")
    {
      version = VERSION;
  	  addPositional ("in", "Protein FASTA file");
    }


	
	void body () const final
  {
	  const string inFName  = getArg ("in");
	  

	  Multifasta faIn (inFName, true);
	  while (faIn. next ())
	  {
  	  const Peptide pep (faIn, 1024 * 1024, false);
	    try
  	  {
  	    pep. qc ();
  	    const size_t stop_pos = pep. seq. find ('*');
  	    if (stop_pos != pep. seq. size () - 1)
  	      throw runtime_error ("stop_codon");
  	    if (pep. seq. front () != 'M')
  	      throw runtime_error ("start_codon");
  	  }
  	  catch (const exception &e)
  	  {
  	    cout << pep. name << '\t' << e. what () << endl;
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



