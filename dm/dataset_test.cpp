// dataset_test.cpp

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
*   Test dataset.{hpp,cpp}
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
using namespace DM_sp;
#include "../version.inc"



namespace
{

const size_t dataSize = 30000;  // PAR
VectorOwn <Distribution> distributions;
	
	
	
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Test dataset.cpp")
  	{
  	  version = VERSION;
  	//addPositional ("seed", "Seed for random numbers");
  	  addPositional ("go", "Go");
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
	    cout << "Estimating: " << distr. name << endl;
	  
		distr. qc ();
		
    distr. setSeed (seed_global);

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

    unique_ptr<const Analysis1> an (distr. getAnalysisCheck ());
    distr. qc ();
       
    T distr_est (distr);  // Set fixed parameters
    ASSERT (distr_est. getAnalysis ());
    distr_est. qc ();
    distr_est. estimate ();

    if (verbose ())
    {
	    distr_est. saveText (cout);
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
	  if (verbose ())
	    cout << "Searching: " << distr. name << endl;
		distr. qc ();
		if (const UniDistribution* ud = distr. asUniDistribution ())
      if (! ud->stdBounds ())
    	  return;
		
    distr. setSeed (seed_global);

    Dataset ds;
    distr. simulate (ds, dataSize);
    ds. qc ();

    unique_ptr<const Analysis1> an (distr. getAnalysisCheck ());
    an->qc ();


    if (verbose ())
      cout << "Search:" << endl;
    const Distribution* distrBest = nullptr;
    Real entropyBest = inf;
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
		    distrTry->saveText (cout);
     	  cout << "  Entropy^: " << entropy_est 
             << "  P-value = " << distrTry->getFitness_entropy () << endl;
		  }
	    distrTry->qc ();
	    if (minimize (entropyBest, entropy_est))
	    	distrBest = distrTry;
	  }
    if (distrBest && ! distrBest->similar (distr, delta))
    {
      distr. saveText (cout);
      distrBest->saveText (cout);
      ERROR;
    }

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
  		    distrTry->saveText (cout);
       	  cout << "  p-value: " << pVal << endl;
  		  }
      }

    if (verbose ())
      cout << endl;
	}



	void body () const final
	{
		distributions << new Bernoulli ()
		              << new Categorical ()
		              << new Binomial ()
		              << new UniformDiscrete ()
		              << new Geometric ()
		              << new Zipf ()
		              << new Normal ()
		              << new Exponential ()
		              << new Cauchy ()
		              << new Beta1 ()
		              << new MultiNormal ()
		              ;
		
		
    section ("Binomial", false);
    {
	    Binomial bin;
	    bin. setParam (20, 0.3);
	    bin. qc ();
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
    
    section ("Normal", false);
    {
    	Normal normal;
    	normal. setParam (0.0, 1.0);
    	normal. qc ();
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
    	// getQuantile()
    	for (Real val = -5.0; val <= 5.0; val += 0.1)
    	{
      	const Prob p = normal. cdf (val);  
      	const Real quantile = normal. getQuantile (p);
      	ASSERT_EQ (val, quantile, 1e-6);
      }
	    Normal* other = normal. copy ();
	    delete other;
    }
    {
	    Normal normal;
	    const Real delta = 0.03;
	    normal. setParam (3.0, 2.0);  estimate (normal, delta);
	    normal. setParam (4.0, 1.0);  estimate (normal, delta);
	    findDistribution (normal, delta);
		  if (verbose ())
	      cout << endl;
    }

    section ("Chi2", false);
    {
    	Chi2 chi2;
    	chi2. setParam (1.0);
    	chi2. qc ();
    	ASSERT_EQ (chi2. pdf (1.0), 0.241971, 1e-6);
    	ASSERT_EQ (chi2. pdf (2.0), 0.103777, 1e-6);
    	ASSERT_EQ (chi2. pdf (3.0), 0.0513934, 1e-6);
    	ASSERT_EQ (chi2. cdf (0.455), 0.5, 1e-4);
    	ASSERT_EQ (chi2. cdf (1.642), 0.8, 1e-4);
    	ASSERT_EQ (chi2. cdf (2.706), 0.900, 1e-4);
    	ASSERT_EQ (chi2. cdf (3.841), 0.950, 1e-4);
    	ASSERT_EQ (chi2. cdf (5.024), 0.975, 1e-4);
    	ASSERT_EQ (chi2. cdf (5.412), 0.980, 1e-4);
    	ASSERT_EQ (chi2. cdf (6.635), 0.990, 1e-5);
    	ASSERT_EQ (chi2. cdf (7.879), 0.995, 1e-5);
    	ASSERT_EQ (chi2. cdf (10.828), 0.999, 1e-6);
    	// getQuantileComp()
    	for (Real val = 0.5; val <= 5.0; val += 0.1)
    	{
      	const Prob p = chi2. cdf (val);  
      	const Real quantile = chi2. getQuantileComp (p, 0.0, 1e3);  // PAR
      	ASSERT_EQ (val, quantile, 1e-6);
      }
	    Chi2* other = chi2. copy ();
	    delete other;
    }
    if (false)
    {
	    Chi2 chi2;
	    const Real delta = 0.03;
	    chi2. setParam (3.0);  estimate (chi2, delta);
	    chi2. setParam (4.0);  estimate (chi2, delta);
	    findDistribution (chi2, delta);
		  if (verbose ())
	      cout << endl;
    }

    section ("Zipf", false);
    {
	    Zipf zipf;
	    const Real delta = 0.03;
      zipf. loBound = 1;  zipf. hiBound = 20;  zipf. setParam (1.7);  estimate (zipf, delta);
	    zipf. loBound = 2;  zipf. hiBound = 20;  zipf. setParam (1.7);  estimate (zipf, delta);
	    zipf. loBound = 5;  zipf. hiBound = 20;  zipf. setParam (1.7);  estimate (zipf, delta);
  	  if (verbose ())
	      cout << endl;
	  }
	  
    section ("Cauchy", false);
	  {
	  	Cauchy c;
	    const Real delta = 0.05;
	    c. setParam (3, 2);  estimate (c, delta);
	    c. setParam (4, 1);  estimate (c, delta);
	    findDistribution (c, delta);
		  if (verbose ())
	      cout << endl;
	  }
	  
	#if 0
    section ("Logistic", false);
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
	
	  // MinDistribution
	  section ("MaxDistribution", false);
	  {
	    Chi2 chi2;
    	chi2. setParam (1.0);
    	chi2. qc ();
    	MaxDistribution maxD;
    	maxD. setParam (& chi2, 300);
    	maxD. qc ();
      if (verbose ())
      	for (Real x = 0.0; x < 300.0; x += 10.0)  
      	  cout << x << '\t' << 1.0 - maxD. cdf (x) << endl;
	    MaxDistribution* other = maxD. copy ();
	    delete other;
	  }
	
	  // Beta1 vs. Uniform ??
	
    section ("MultiNormal", false);
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
	  	mn. qc ();
	  	estimate (mn, 0.1);
	    findDistribution (mn, 0.05);  // PAR
		  if (verbose ())
	      cout << endl;
	  	
	    MultiNormal* other = mn. copy ();
	    delete other;

	  	// MultiNormal(1) vs. Normal: same parameters, getInfoMean() etc. ??
	  }
	  
    section ("Clustering", false);
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
    		mn. setSeed (seed_global);
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
    		mn. setSeed (seed_global + 1);
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
		    cl. saveText (cout);
		    const Space1<ProbAttr1> spOut (cl. createSpace (ds));
		    cout << endl;
		    ds. saveText (cout);
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


