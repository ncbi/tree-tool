// prot2triplets.cpp

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
*   AA tiplet statistics. Create a file with short proteins
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "seq.hpp"
using namespace Seq_sp;



namespace 
{


const string shortFName = "short";



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("AA tiplet statistics. Create file \"" + shortFName + "\" with short proteins")
    {
  	  addPositional ("in", "Input FASTA file of proteins");
    }


	
	void body () const final
  {
	  const string inFName = getArg ("in");
  
  
    Vector<size_t> aa2index (128, no_index);
    const string alphabet (extPeptideAlphabet);
    FOR (size_t, i, alphabet. size ())
      aa2index [(size_t) alphabet [i]] = i;
  
    typedef Vector<size_t> Aas;
    typedef Vector<Aas> Duplets;
    typedef Vector<Duplets> Triplets;
    Triplets triplets (alphabet. size ());
    ITER (Triplets, it2, triplets)
    {
      it2->resize (alphabet. size ());
      ITER (Duplets, it1, *it2)
        it1->resize (alphabet. size (), 0);
    }  
  
    const size_t window = 20;  // PA
    
    
    OFStream shortF (shortFName);
    

    size_t nSeq = 0;
    size_t nShort = 0;
	  Multifasta fIn (inFName, true);
	  Set<string/*size()=3*/> tripletSet;
	  while (fIn. next ())
	  {
	    Peptide p (fIn, 1000, false);
	    ++nSeq;
	    
	    string& seq = p. seq;
	    if (seq. size () < window)
	    {
	    	++nShort;
    	  p. printCase (shortF, true);
	    	continue;
	    }
	    
	    tripletSet. clear ();
	    ASSERT (seq. size () >= 3);
	    seq. resize (seq. size () + 1, '\0');
	    FOR_REV (size_t, i, seq. size () - 3)
	    {
	      seq [i + 3] = '\0';
	      tripletSet << & seq [i];
	    }
	    
	    CONST_ITER (Set<string>, it, tripletSet)
	      ++ triplets [aa2index [(size_t) (*it) [0]]]
	                  [aa2index [(size_t) (*it) [1]]]
	                  [aa2index [(size_t) (*it) [2]]];
	  }


    cout << "# Sequences = " << nSeq << endl;
    cout << "# Short sequences = " << nShort << endl;
	  
	  const size_t nLong = nSeq - nShort;
	  FOR (size_t, i, alphabet. size ())
	  FOR (size_t, j, alphabet. size ())
	  FOR (size_t, k, alphabet. size ())
	    if (triplets [i] [j] [k])
	    	cout << alphabet [i]
	    	         << alphabet [j]
	    	         << alphabet [k]
	    	         << " " << scientific << (double) triplets [i] [j] [k] / (double) nLong
	    	         << endl;
  }
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



