// beta1_test.cpp

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
using namespace DM_sp;



namespace
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Test Beta1 estimation")
  	{
  	  addPositional ("size", "Data size");
  	  addPositional ("seed", "Seed for random numbres (> 0)");
  	}



	void body () const final
	{
		const size_t dataSize = str2<size_t> (getArg ("size"));
		const ulong seed      = str2<ulong>  (getArg ("seed"));


    Dataset ds;
    auto& betaAttr = * new RealAttr1 ("Beta1", ds);
    auto& zipfAttr = * new IntAttr1 ("Zipf", ds);
    const Real beta1alpha = 0.07;  // PAR
    {
	    Beta1 beta1;
	    beta1. setSeed (seed);
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


