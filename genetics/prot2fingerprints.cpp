// prot2fingerprints.cpp

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
*   Find the triplet fingerprints of proteins
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

struct TripletF : Named
{
	size_t prots {0};
	OFStream* f {nullptr};
	
	TripletF (const string &dir,
	          const string &triplet)
	  : Named (dir + "/" + triplet)
	  { QC_ASSERT (triplet. size () == 3); }
 ~TripletF ()
    { delete f; }
  
  void save (const Peptide &p)
    { prots++;
    	if (! f)
    		f = new OFStream (name);
    	p. printCase (*f, true);
    }
};



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Find the triplet fingerprints of proteins")
    {
      version = VERSION;
  	  addPositional ("in", "Input file of proteins");
  	  addPositional ("triplets", "Input file with triplet statistics");
  	  addPositional ("out_dir", "Output directory for proteins");
    }


	
	void body () const final
  {
	  const string inFName      = getArg ("in");
	  const string tripletFName = getArg ("triplets");
	  const string out_dir      = getArg ("out_dir");
  
  
    // Cf. prot2triplets.cpp

    
    Vector<size_t> aa2index (128, no_index);
    const string alphabet (extPeptideAlphabet);
    FOR (size_t, i, alphabet. size ())
      aa2index [(size_t) alphabet [i]] = i;
  
    // BLAST database size
    typedef Vector<unsigned long> Aa2size; 
    typedef Vector<Aa2size> Duplet2size;
    Vector<Duplet2size> triplet2size (alphabet. size ());
    for (Duplet2size& ds : triplet2size)
    {
      ds. resize (alphabet. size ());
      for (Aa2size& as : ds)
        as. resize (alphabet. size (), 0);
    }

	  string triplet;  triplet. resize (3);  // Auxiliary

    // # Proteins
    typedef Vector<TripletF> Aa2prots; 
    typedef Vector<Aa2prots> Duplet2prots;
    typedef Vector<Duplet2prots> Triplet2prots;
    Triplet2prots triplet2prots (alphabet. size ());
    FOR (size_t, i, alphabet. size ())
    {
     	triplet [0] = alphabet [i];
    	Duplet2prots& duplet2prots = triplet2prots [i];
      duplet2prots. resize (alphabet. size ());
      FOR (size_t, j, alphabet. size ())
      {
     	  triplet [1] = alphabet [j];
      	Aa2prots& aa2prots = duplet2prots [j];
        FOR (size_t, k, alphabet. size ())
        {
     	    triplet [2] = alphabet [k];
          aa2prots << TripletF (out_dir, triplet);
        }
      }
    }

    // Frequency
    typedef Vector<float> Aa2freq;  
    typedef Vector<Aa2freq> Duplet2freq;
    typedef Vector<Duplet2freq> Triplet2freq;
    Triplet2freq triplet2freq (alphabet. size ());
    for (Duplet2freq& d : triplet2freq)
    {
      d. resize (alphabet. size ());
      for (Aa2freq& a : d)
        a. resize (alphabet. size (), 0);
    }
    // Reading frequencies
    {
      ifstream f (tripletFName. c_str ());
      ASSERT (f. good ());
      while (! f. eof ())
      {
      	if (f. peek () == '#')
      	{
      		skipLine (f);
      		continue;
      	}
      	float freq = 0;
      	triplet. clear ();
        f >> triplet >> freq;
        if (f. eof ())
        {
          ASSERT (! freq);
          break;
        }
        ASSERT (f. peek () == '\n');
        skipLine (f);
        ASSERT (triplet. size () == 3);
        ASSERT (freq > 0);
        ASSERT (freq < 1);
        triplet2freq [aa2index [(size_t) triplet [0]]]
			               [aa2index [(size_t) triplet [1]]]
			               [aa2index [(size_t) triplet [2]]] = freq;
      }
    }
  
  
    const size_t window = 20;  // PAR


    size_t nSeq = 0;
    size_t nShort = 0;
	  Multifasta fIn (inFName, true);
	  Vector<double> weights;  weights. reserve (10000);
	  Set<string/*size()=3*/> fingerprints;
	  string fingerprint;  
	  triplet. resize (3);
	  while (fIn. next ())
	  {
	    const Peptide p (fIn, 1000, false);
	    ++nSeq;
	    
	    const string& seq = p. seq;
	    if (seq. size () < window)
	    {
	    	++nShort;
	    	continue;
	    }
	    
	    ASSERT (seq. size () >= 3);
	    weights. clear ();
	    FOR (size_t, i, seq. size () - 3 + 1)
	    	weights << triplet2freq [aa2index [(size_t) seq [i + 0]]]
										            [aa2index [(size_t) seq [i + 1]]]
										            [aa2index [(size_t) seq [i + 2]]];
	    
	    fingerprints. clear ();
	    fingerprint. clear ();
	    ASSERT (window >= 3);
	    size_t fpIndex_prev = no_index;
	    double weight_min_min = 1;
	    // Heap is faster ??
	    FOR (size_t, i, seq. size () - window + 1) 
	    {
	    	double weight_min = 1;
	    	size_t fpIndex = no_index;
	      FOR_START (size_t, j, i, i + (window - 3) + 1)
	        if (minimize (weight_min, weights [j]))
	        	fpIndex = j;
	      ASSERT (fpIndex != no_index);
	      if (fpIndex_prev != fpIndex)
	      {
	      	fpIndex_prev = fpIndex;
	      	FOR (size_t, k, 3)
	      	  triplet [k] = seq [fpIndex + k];
	      	ASSERT (triplet. size () == 3);
	        fingerprints << triplet;
	        if (minimize (weight_min_min, weight_min))
	        	fingerprint = triplet;
	      }
	    }
	    ASSERT (! fingerprints. empty ());  // => cover of fIn
	    ASSERT (fingerprint. size () == 3);

	    for (const string& s : fingerprints)
	      triplet2size [aa2index [(size_t) s [0]]]
				             [aa2index [(size_t) s [1]]]
				             [aa2index [(size_t) s [2]]] += seq. size ();

		  triplet2prots [aa2index [(size_t) fingerprint [0]]]
			              [aa2index [(size_t) fingerprint [1]]]
					          [aa2index [(size_t) fingerprint [2]]].
			save (p);
	  }


    cout << "# Sequences = " << nSeq << endl;
    cout << "# Short sequences = " << nShort << endl;

	  FOR (size_t, i, alphabet. size ())
	  FOR (size_t, j, alphabet. size ())
	  FOR (size_t, k, alphabet. size ())
    	cout << alphabet [i]
    	         << alphabet [j]
    	         << alphabet [k]
    	         << " " << triplet2size  [i] [j] [k] 
    	         << " " << triplet2prots [i] [j] [k]. prots 
    	         << endl;
  }
};



}  // namespace



int main(int argc, 
         const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}




