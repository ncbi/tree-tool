// beta1_test.cpp

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
*   Test beta(1) distribution
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
using namespace DM_sp;
#include "version.inc"



namespace
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Test Beta1 estimation")
  	{
  	  version = VERSION;
  	  addPositional ("size", "Data size");
  	//addPositional ("seed", "Seed for random numbres (> 0)");
  	}



	void body () const final
	{
		const size_t dataSize = str2<size_t> (getArg ("size"));
  //const ulong seed      = str2<ulong>  (getArg ("seed"));


    Dataset ds;
    auto& betaAttr = * new RealAttr1 ("Beta1", ds);
    auto& zipfAttr = * new IntAttr1 ("Zipf", ds);
    const Real beta1alpha = 0.07;  // PAR
    {
	    Beta1 beta1;
	    beta1. setSeed (seed_global);
	    beta1. setParam (beta1alpha);  
	  //cout << beta1. cdf (1.0/(Real)dataSize) << endl;  
	    FOR (size_t, i, dataSize)  
	    {
	    	const size_t n = ds. appendObj (toString (i + 1));
	    	const Real p = beta1. rand ();
	    	ASSERT (isProb (p));
	      betaAttr [n] = p;
	      zipfAttr [n] = (int) floor (p * ((Real) dataSize - 1)) + 1;
	    }
	  }
    ds. qc ();
  //ds. printCsv (cout);  
    
    const Sample sm (ds);

    {
	    const UniVariate<NumAttr1> an (sm, betaAttr);
	    
	    Beta1 beta_est;
	    beta_est. analysis = & an;
	  //beta_est. setParam (...);  // PAR
	    beta_est. estimate ();
	    if (verbose ())
	    {
	      beta_est. print (cout);
  	    cout << endl;
  	    cout << "P-value = " << beta_est. getFitness_entropy () << endl;
  	    cout << endl;
  	  }
	    ASSERT_EQ (beta_est. alpha, beta1alpha, 0.01);  // PAR
	  }

    
    {
	    const UniVariate<IntAttr1> an (sm, zipfAttr);
	    
	    Zipf zipf_est;
	    zipf_est. analysis = & an;
	  //zipf_est. setParam (1.5);  // PAR
      zipf_est. hiBound = (Real) dataSize;  // ??
	    zipf_est. estimate ();
	    if (verbose ())
	    {
  	    zipf_est. print (cout);
  	    cout << endl;
  	    cout << "P-value = " << zipf_est. getFitness_entropy () << endl;  // = 0.98 ??
  	  //cout << "Beta-alpha = " << zipf_est. getBetaAlpha () << endl;
  	    cout << endl;
  	  //cout << zipf_est. pmf (1) << endl;
  	}
	  //ASSERT ??
	  }
	  

  #if 0
    typedef map <uint, uint/*freq*/> Freq;
	  Freq freq;
	  for (Iterator it (zipfAttr); it ();)  
	    freq [(uint) DM_sp::round ((*zipfAttr) [*it])] ++;
	  CONST_ITER (Freq, it, freq)
	    cout << (*it). first << "," << (*it). second << endl;
	#endif
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


