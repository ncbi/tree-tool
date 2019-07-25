// testZipf.cpp

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
*   Test of Zipf distribution
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
using namespace DM_sp;



namespace
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Test Zipf estimation")
  	{
  	  addPositional ("size", "Data size");
  	//addPositional ("seed", "Seed for random numbers (> 0)");
  	}



	void body () const final
	{
		const size_t dataSize = str2<size_t> (getArg ("size"));
	//const ulong  seed     = str2<ulong> (getArg ("seed"));


  #if 0
    cout << "zeta(1.5) = " << zeta(1.5) << endl;
    cout << "zeta(2) = " << zeta(2) << endl;
    cout << "zeta(3) = " << zeta(3) << endl;
    cout << "zeta(4) = " << zeta(4) << endl;
  #endif


	  map <uint/*category*/, uint/*freq*/> categs;
    {
	    Zipf zipf;
	    zipf. setSeed (seed_global);
	    zipf. setParam (1.5 /*, 10000*/);  // PAR
	  //ASSERT_EQ (zipf. c, 1/zeta(1.5), 1e-6); 
      FOR (size_t, i, dataSize)  
        categs [(uint) zipf. randDiscrete ()] ++;
	  }


    // Original ranks
    {
	    Dataset ds;
	    auto catAttr = new IntAttr1 ("Category", ds);
	    for (const auto& it : categs)
	    {
	    	const uint key = it. first;
	    	const size_t n = ds. appendObj (toString (key));
	      (*catAttr) [n] = (int) key;
	      const_cast <Obj*> (ds. objs [n]) -> mult = it. second;
	    }
	    ds. qc ();
	  //ds. printCsv (cout);  
	  
	    const Sample sm (ds);	    
	    const UniVariate<IntAttr1> an (sm, *catAttr);
	    
	    Zipf zipf_est;
	    zipf_est. setParam (1.5 /*, 10000*/);  // PAR
	    zipf_est. analysis = & an;
	    zipf_est. estimate ();
	    zipf_est. print (cout);
	    cout << endl;
	    cout << "P-value = " << zipf_est. getFitness_entropy () << endl;
	    cout << endl;
	  }

    
    Vector <uint/*freq*/> freqVec;
    for (const auto& it : categs)
      freqVec << it. second;
    sort (freqVec);
    freqVec. reverse ();


    // Data ranks: Zipf::alpha has a positive bias
    {
	    Dataset ds;
	    auto catAttr = new IntAttr1 ("Category", ds);
	    FOR (size_t, i, freqVec. size ())
	    {
	    	const int key = (int) i + 1;
	    	const size_t n = ds. appendObj (toString (key));
	      (*catAttr) [n] = key;
	      const_cast <Obj*> (ds. objs [n]) -> mult = freqVec [i];
	    }
	    ds. qc ();
	  //ds. printCsv (cout);  
	    
	    const Sample sm (ds);	    
	    const UniVariate<IntAttr1> an (sm, *catAttr);
	    
	    Zipf zipf_est;
	    zipf_est. setParam (1.5 /*, (uint) freqVec. size ()*/); 
	    zipf_est. analysis = & an;
	    zipf_est. estimate ();
	    zipf_est. print (cout);
	    cout << endl;
	    cout << "P-value = " << zipf_est. getFitness_entropy () << endl;
	    cout << endl;
	  }
    
    
  #if 0
    // Data ranks, reverse cumulative
    {
    	Vector <uint> freqCum (freqVec. size ());
    	uint s = 0;
    	FOR_REV (size_t, i, freqCum. size ())
    	{
    		s += freqVec. at (i);
    		freqCum. at (i) = s;
    	}
	    Dataset ds;
	    RealAttr1* catAttr = new RealAttr1 ("Category", ds);
	    FOR (size_t, i, freqCum. size ())
	    {
	    	const uint key = (uint) i + 1;
	    	const size_t n = ds. appendObj (toString (key));
	      (*catAttr) [n] = key;
	      const_cast <Obj*> (ds. objs [n]) -> setMult (freqCum. at (i));
	    }
	    ds. qc ();
	  //ds. printCsv (cout);  
	    
	    const Sample sm (ds);	    
	    const AttrAnalysis an (sm, *catAttr);
	    
	    Zipf zipf_est;
	    zipf_est. setParam (1.5, (uint) freqCum. size ()); 
	    zipf_est. analysis = & an;
	    zipf_est. estimate ();
	    zipf_est. print (cout);
	    cout << endl;
	    cout << "P-value = " << zipf_est. getFitness_entropy () << endl;
	    cout << endl;
	  }
	#endif
    
    
    // Frequency distribution
    map <uint/*freq*/, uint/*freq*/> freqs;
    for (const auto& it : categs)
      freqs [it. second] ++;
    {
	    Dataset ds2;
	    auto freqAttr = new IntAttr1 ("Freq", ds2);
	    for (const auto& it : freqs)
	    {
	    	const uint key = it. first;
	    //if (key > 20)  // PAR
	    	//continue;
	    	const size_t n = ds2. appendObj (toString (key));
	      (*freqAttr) [n] = (int) key;
	      const_cast <Obj*> (ds2. objs [n]) -> mult = it. second;
	    }
	    ds2. qc ();
	  //ds2. printCsv (cout); 
	    
	    const Sample sm (ds2);	    
	    const UniVariate<IntAttr1> an (sm, *freqAttr);
	    
	    Zipf zipf_est;
	    zipf_est. setParam (1.5 /*, (uint) dataSize*/);  
	    zipf_est. analysis = & an;
	    zipf_est. estimate ();
	    zipf_est. print (cout);
	    cout << endl;
	  //cout << "A = " << zipf_est. freqAlpha2alpha () << endl;
	    cout << "P-value = " << zipf_est. getFitness_entropy () << endl;
	    cout << endl;
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


