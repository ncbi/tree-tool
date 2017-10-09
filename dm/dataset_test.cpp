// dataset_test.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
using namespace DM_sp;



namespace
{

const size_t dataSize = 30000;  // PAR
ulong seed = 0;
VectorOwn <Distribution> distributions;
	
	
	
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Test dataset.cpp")
  	{
  	  addPositional ("seed", "Seed for random numbers");
  	}
	
	
	
	void checkFitness (const string &name,
	                   Distribution &distr) const
	{
	  const Prob pVal = distr. getFitness_entropy ();
	  const bool bad = pVal < 0.01;  // H_0 is rejected
	  if (verbose () || bad)
	    cout << fixed << name << " distribution: P-value of Entropy^ = " << pVal << endl;
	  ASSERT (! bad);
	}
	

	
  template <typename T/*:Distribution*/>
	void estimate (T &distr,
	               Real delta) const
	// Update: distr.randGen
	{
	  if (verbose ())
	    cout << distr. name << endl;
	  
		distr. qc ();
		
    distr. setSeed (seed);

    bool checkEntropyOk = true;
    if (const UniDistribution* distr1 = distr. asUniDistribution ())
	    if (! distr1->stdBounds ())
	    	checkEntropyOk = false;
    if (checkEntropyOk)
    {
    	Unverbose unv;
	    const Real info_mean = distr. getInfoMean ();
	    const Real info_var  = distr. getInfoVar ();
	    if (! (   isNan (info_mean)
	    	     && isNan (info_var)
	    	    )
	    	  )
	    {
		    Real mean, var;
		    distr. info_meanVar_MC (100000, mean, var);
		    if (verbose ())
		    {
		    	cout << endl;
		      cout << "Information:" << endl;
		      cout << "Mean: theory: " << info_mean       << "  experiment: " << mean << endl; 
		      cout << "SD:   theory: " << sqrt (info_var) << "  experiment: " << sqrt (var) << endl; 
		    	cout << endl;
		    }
		    IMPLY (! isNan (info_mean), eqReal (info_mean, mean, 2e-2));
		    IMPLY (! isNan (info_var),  eqReal (sqrt (info_var), sqrt (var), 5e-2));
		  }
		}

    Dataset ds;
    distr. simulate (ds, dataSize);
    ds. qc ();

    Common_sp::AutoPtr<const Analysis1> an (distr. getAnalysisCheck ());
    distr. qc ();
       
    T distr_est (distr);  // Set fixed parameters
    ASSERT (distr_est. getAnalysis ());
    distr_est. qc ();
    distr_est. estimate ();

    if (verbose ())
    {
	    distr_est. print (cout);
	    cout << endl;
	  }
    ASSERT (distr_est. similar (distr, delta));

    checkFitness ("MLE", distr_est);

    // QC of Rand
  	distr_est = distr;  
    distr_est. qc ();
    checkFitness ("Real", distr_est);

    if (verbose ())
      cout << endl;
	}
	
	

	void findDistribution (Distribution &distr,
	                       Real delta) const
	{
		distr. qc ();
		if (const UniDistribution* ud = distr. asUniDistribution ())
      if (! ud->stdBounds ())
    	  return;
		
    distr. setSeed (seed);

    Dataset ds;
    distr. simulate (ds, dataSize);
    ds. qc ();

    Common_sp::AutoPtr<const Analysis1> an (distr. getAnalysisCheck ());
    distr. qc ();


    if (verbose ())
      cout << "Search:" << endl;
    const Distribution* distrBest = nullptr;
    Real entropyBest = INF;
    for (const Distribution* distrTry : distributions)
    {
      ASSERT (& distr != distrTry);
      const_cast <Distribution*> (distrTry) -> shareAnalysis (distr);
	    distrTry->qc ();
	    if (! distrTry->getAnalysis ())
	      continue;
	    if (verbose ())
	      cout << "Trying: " << distrTry->name << endl;
	    const_cast <Distribution*> (distrTry) -> estimate ();
	    const Real entropy_est = distrTry->getEntropy_est ();
		  if (verbose ())
		  {
		    distrTry->print (cout);
     	  cout << "  Entropy^: " << entropy_est 
             << "  P-value = " << distrTry->getFitness_entropy () << endl;
		  }
	    distrTry->qc ();
	    if (minimize (entropyBest, entropy_est))
	    	distrBest = distrTry;
	  }
    IMPLY (distrBest, distrBest->similar (distr, delta));

    if (verbose ())
    	cout << "Altenative distributions:" << endl;
    for (const Distribution* distrTry : distributions)
    if (distrTry != distrBest)
    {
      const_cast <Distribution*> (distrTry) -> shareAnalysis (distr);
	    if (! distrTry->getAnalysis ())
	      continue;
    	const Prob pVal = distrTry->getFitness_entropy (entropyBest);
    	ASSERT (pVal <= 0.001);
		  if (verbose ())
		  {
		    distrTry->print (cout);
     	  cout << "  p-value: " << pVal << endl;
		  }
    }

    if (verbose ())
      cout << endl;
	}



	void body () const final
	{
		seed = str2<ulong> (getArg ("seed"));
		
		
		distributions << new Bernoulli ()
		              << new Categorical ()
		              << new Binomial ()
		              << new UniformDiscrete ()
		              << new Geometric ()
		              << new Zipf ()
		              << new Normal ()
		              << new Cauchy ()
		              << new Beta1 ()
		              << new MultiNormal ()
		              ;
		
		
    // Binomial
    {
	    Binomial bin;
	    bin. setParam (20, 0.3);
	    ASSERT_EQ (bin. cdf ( 0), 0.0008, 1e-4);
	    ASSERT_EQ (bin. cdf ( 1), 0.0076, 1e-4);
	    ASSERT_EQ (bin. cdf ( 2), 0.0355, 1e-4);
	    ASSERT_EQ (bin. cdf ( 3), 0.1071, 1e-4);
	    ASSERT_EQ (bin. cdf ( 4), 0.2375, 1e-4);
	    ASSERT_EQ (bin. cdf ( 5), 0.4164, 1e-4);
	    ASSERT_EQ (bin. cdf ( 6), 0.6080, 1e-4);
	    ASSERT_EQ (bin. cdf ( 7), 0.7723, 1e-4);
	    ASSERT_EQ (bin. cdf ( 8), 0.8867, 1e-4);
	    ASSERT_EQ (bin. cdf ( 9), 0.9520, 1e-4);
	    ASSERT_EQ (bin. cdf (10), 0.9829, 1e-4);
	    ASSERT_EQ (bin. cdf (11), 0.9949, 1e-4);
	    ASSERT_EQ (bin. cdf (12), 0.9987, 1e-4);
	    ASSERT_EQ (bin. cdf (13), 0.9997, 1e-4);
	    ASSERT_EQ (bin. cdf (14), 1.0000, 1e-4);
	    ASSERT_EQ (bin. cdf (20), 1.0000, 1e-4);
	    Binomial* other = bin. copy ();
	    delete other;
	  }
    {
	    Binomial bin;
	    const Real delta = 0.005;
	    bin. loBound = 0;  bin. hiBound = 20;  bin. setParam (20, 0.7);   estimate (bin, delta);
	    bin. loBound = 5;  bin. hiBound = 18;  bin. setParam (20, 0.7);   estimate (bin, delta);
	    bin. loBound = 1;  bin. hiBound = 20;  bin. setParam (20, 0.7);   estimate (bin, delta);
	    bin. loBound = 0;  bin. hiBound = 15;  bin. setParam (20, 0.7);   estimate (bin, delta);
	    bin. loBound = 1;  bin. hiBound = 20;  bin. setParam (20, 0.99);  estimate (bin, delta);
		  if (verbose ())
	      cout << endl;
	  }
    
    // Normal    
    {
    	Normal normal;
    	normal. setParam (0, 1);
    	ASSERT_EQ (normal. cdf (0.155), 0.561589, 1e-6);
    	ASSERT_EQ (normal. cdf (0.344), 0.634577, 1e-6);
    	ASSERT_EQ (normal. cdf (0.545), 0.707123, 1e-6);
    	ASSERT_EQ (normal. cdf (0.696), 0.756786, 1e-6);
    	ASSERT_EQ (normal. cdf (0.853), 0.803170, 1e-6);
    	ASSERT_EQ (normal. cdf (1.045), 0.851989, 1e-6);
    	ASSERT_EQ (normal. cdf (1.245), 0.893434, 1e-6);
    	ASSERT_EQ (normal. cdf (1.393), 0.918190, 1e-6);
    	ASSERT_EQ (normal. cdf (-1.555), 1-0.940027, 1e-6);
    	ASSERT_EQ (normal. cdf (2.444), 0.992737, 1e-6);
	    Normal* other = normal. copy ();
	    delete other;
    }
    {
	    Normal normal;
	    const Real delta = 0.03;
	    normal. setParam (3, 2);  estimate (normal, delta);
	    normal. setParam (4, 1);  estimate (normal, delta);
	    findDistribution (normal, delta);
		  if (verbose ())
	      cout << endl;
    }

    // Zipf
    {
	    Zipf zipf;
	    const Real delta = 0.03;
      zipf. loBound = 1;  zipf. hiBound = 20;  zipf. setParam (1.7);  estimate (zipf, delta);
	    zipf. loBound = 2;  zipf. hiBound = 20;  zipf. setParam (1.7);  estimate (zipf, delta);
	    zipf. loBound = 5;  zipf. hiBound = 20;  zipf. setParam (1.7);  estimate (zipf, delta);
  	  if (verbose ())
	      cout << endl;
	  }
	  
	  // Cauchy
	  {
	  	Cauchy c;
	    const Real delta = 0.05;
	    c. setParam (3, 2);  estimate (c, delta);
	    c. setParam (4, 1);  estimate (c, delta);
	    findDistribution (c, delta);
		  if (verbose ())
	      cout << endl;
	  }
	  
	  // Logistic 
	#if 0
	  // ??
	  {
	  	Logistic lg;
	    const Real delta = 0.05;
	    lg. setParam (3, 2);  estimate (lg, delta);
	    lg. setParam (4, 1);  estimate (lg, delta);
	    findDistribution (lg, delta);
		  if (verbose ())
	      cout << endl;
	  }
	#endif
	
	  // Beta1 vs. Uniform ??
	
	  // MultiNormal
	  {
		  if (verbose ())
	      cout << "MultiNormal" << endl;
	  	const size_t dim = 2;
	  	MultiNormal mn;
	  	mn. setDim (dim);
	  	mn. mu. putAll (1);
	  	mn. sigmaExact. putDiag (0, 1);
	  	mn. sigmaExact. putDiag (1, 1);
	  	mn. sigmaExact. putSymmetric (1, 0, 0.5);
	  	mn. setParam ();
	  	estimate (mn, 0.1);
	    findDistribution (mn, 0.05);  // PAR
		  if (verbose ())
	      cout << endl;
	  	
	    MultiNormal* other = mn. copy ();
	    delete other;

	  	// MultiNormal(1) vs. Normal: same parameters, getInfoMean() etc. ??
	  }
	  
	  // Clustering
	  {
		  if (verbose ())
	      cout << "Clustering" << endl;
	    Dataset ds;
	    const size_t nAttr = 2;
	    Vector<RealAttr1*> attrs;  attrs. reserve (nAttr);
	    FOR (size_t, i, nAttr)
	      attrs << new RealAttr1 ("Attr_" + toString (i + 1), ds, 2);
	    // Component 0
	    {
		  	MultiNormal mn;
		  	mn. setDim (nAttr);
		  	mn. mu [0] = 0;
		  	mn. mu [1] = 0;
		  	mn. sigmaExact. putDiag (0, 1);
		  	mn. sigmaExact. putDiag (1, 1);
		  	mn. sigmaExact. putSymmetric (1, 0, 0);
  	    ASSERT (mn. sigmaExact. psd);
		  	mn. setParam ();
    		mn. setSeed (seed);
		    FOR (size_t, i, 100)  
		    {
		    	const size_t objNum = ds. appendObj ("1_" + toString (i + 1));
		    	mn. randVariable ();
			    FOR (size_t, j, mn. getDim ())
			      (* attrs [j]) [objNum] = mn. variable [j];
		    }
		  }
	    // Component 1
	    {
		  	MultiNormal mn;
		  	mn. setDim (nAttr);
		  	mn. mu [0] = 3;
		  	mn. mu [1] = 3;
		  	mn. sigmaExact. putDiag (0, 1);
		  	mn. sigmaExact. putDiag (1, 1);
		  	mn. sigmaExact. putSymmetric (1, 0, 0);
	  	  ASSERT (mn. sigmaExact. psd);
		  	mn. setParam ();
    		mn. setSeed (seed + 1);
		    FOR (size_t, i, 100)  
		    {
		    	const size_t objNum = ds. appendObj ("2_" + toString (i + 1));
		    	mn. randVariable ();
			    FOR (size_t, j, mn. getDim ())
			      (* attrs [j]) [objNum] = mn. variable [j];
		    }
		  }
		  ds. qc ();
		  const Space1<NumAttr1> sp (ds, true);
		  const Sample sm (ds);
		  if (verbose ())
		  {
		    cout << endl;
		    cout << "Clustering:" << endl << endl;
		  }
		  const Clustering cl (sm, sp, 100, 0.2, false, 0.07);  // PAR
		  cl. qc ();
		  if (verbose ())
		  {
		    cout << endl;
		    cl. print (cout);
		    const Space1<ProbAttr1> spOut (cl. createSpace (ds));
		    cout << endl;
		    ds. print (cout);
		  }
		  ASSERT (cl. mixt. components. size () == 2);
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


