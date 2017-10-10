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
      addPositional ("corr_min", "Minimum correlation between two attributes to report");
    }
	
	
	
	void body () const final
	{
		const string inFName = getArg ("file");
		const Real corr_min  = str2real (getArg ("corr_min"));


    const Dataset ds (inFName);
    const Sample sample (ds);
    cout << "mult_sum: " << sample. mult_sum << endl << endl;
    
    cout << "AttrName\tDefintion\tMissings\tStatistics" << endl;
    for (const auto attrRaw : ds. attrs)
    {
    	cout << attrRaw->name << '\t' << attrRaw->getTypeStr () << '\t' << attrRaw->countMissings ();
      if (const NumAttr1* num = attrRaw->asNumAttr1 ())
	    {
	      Normal normal;
	      const UniVariate<NumAttr1> an (sample, *num);
	      normal. analysis = & an;
	      normal. estimate ();
	      normal. qc ();
	      cout << '\t' << normal. loc << '\t' << normal. scale;
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
    {
      cout << endl;
      const Space1<NumAttr1> space (ds, true);
      const MultiVariate<NumAttr1> an (sample, space);
      MultiNormal mn; 
      mn. qc ();
      mn. analysis = & an;
      mn. estimate ();
      mn. qc ();
      Progress prog ((uint) space. size (), 0);
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


