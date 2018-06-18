// attrs.cpp

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
    : Application ("Print the description of attributes")
    {
      addPositional ("file", dmSuff + "-file");
      addKey ("corr_min", "Min. correlation between two attributes to report", "nan");
      addKey ("outlier_evalue_max", "Max. outlier e-value", "0.1");
    }
	
	
	
	void body () const final
	{
		const string inFName = getArg ("file");
		const Real corr_min  = str2real (getArg ("corr_min"));
		const Prob outlier_eValue_max = str2real (getArg ("outlier_evalue_max"));
		ASSERT (isProb (outlier_eValue_max));


    const Dataset ds (inFName);
    const Sample sample (ds);
    cout << "mult_sum: " << sample. mult_sum << endl << endl;
    
    cout << "AttrName\tDefintion\tMissings\tStatistics" << endl;
    for (const auto attrRaw : ds. attrs)
    {
    	cout << attrRaw->name << '\t' << attrRaw->getTypeStr () << '\t' << attrRaw->countMissings ();
      if (const NumAttr1* num = attrRaw->asNumAttr1 ())
	    {
	    	const ONumber on (cout, num->decimals + 1, false);
	      Normal normal;
	      const UniVariate<NumAttr1> an (sample, *num);
	      normal. analysis = & an;
	      normal. estimate ();
	      normal. qc ();	      
	      cout << '\t' << normal. loc 
	           << '\t' << normal. scale;
		    if (! isNan (outlier_eValue_max))
		    {
			    Sample sample_pure (ds);
			    for (const bool rightTail : {false, true})
			    {
			    	const Real threshold = num->distr2outlier (sample, normal, rightTail, outlier_eValue_max);
			      cout << '\t' << threshold;
			      FFOR (size_t, i, sample_pure. mult. size ())
			        if (   (  rightTail && (*num) [i] > threshold)
			        	  || (! rightTail && (*num) [i] < threshold)
			        	 )
			        	sample_pure. mult [i] = 0;
			    }
		      const UniVariate<NumAttr1> an_pure (sample_pure, *num);
		      normal. analysis = & an_pure;
		      normal. estimate ();
          const Prob pVal = normal. getFitness_entropy ();
		      cout << '\t' << normal. loc 
		           << '\t' << normal. scale
		           << '\t' << pVal;
			  }
	    }
      else if (const BoolAttr1* boolean = attrRaw->asBoolAttr1 ())
	      cout << '\t' << boolean->getProb (sample);
      else if (const NominAttr1* nomin = attrRaw->asNominAttr1 ())
      {
        Common_sp::AutoPtr<Categorical> cat (nomin->getCategorical (sample));
	      cout << "  " << endl;
        cat->print (cout);
	    }
      cout << endl;
	  //other types ??
	  }

    // Correlations
    // Only for NumAttr1 ??
    // Requires: no missings ??
    if (! isNan (corr_min))
    {
      cout << endl;
      const Space1<NumAttr1> space (ds, true);
      const MultiVariate<NumAttr1> an (sample, space);
      MultiNormal mn; 
      mn. qc ();
      mn. analysis = & an;
      mn. estimate ();
      mn. qc ();
      Progress prog (space. size (), 0);
      FOR (size_t, i, space. size ())
      {
        prog (space [i] -> name);
        FOR (size_t, j, i)
        {
        	const Real corr = mn. sigmaExact. get (false, i, j) / sqrt (mn. sigmaExact. getDiag (i) * mn. sigmaExact. getDiag (j));
          if (geReal (corr, corr_min))
          	cout << "corr (" << space [i] -> name << ", " << space [j] -> name << ") = " << corr << endl;
        }
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


