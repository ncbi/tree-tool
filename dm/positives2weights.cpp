// positives2weights.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "numeric.hpp"
#include "matrix.hpp"
#include "dataset.hpp"
//#include "prediction.hpp"
using namespace DM_sp;



namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Analyze proportional positive attributes: logarithmize, centralize to 1, find weights")
    {
      // Input
  	  addPositional ("file", dmSuff + "-file without \"" + dmSuff + "\"");
  	  addPositional ("attrs", "List of 2-way <attr> positive attributes in <file>");
  	//addFlag ("print_beta", "Print the coefficients of predictors");
  	}



	void body () const final
	{
		const string inFName    = getArg ("file");
		const string attrsFName = getArg ("attrs");
	//const bool print_beta   = getFlag ("print_beta");


   	const StringVector attrNames (attrsFName, 256);  // PAR
    
    Dataset ds (inFName);
    ds. qc ();
    
    const string suffix ("_log");

	  Space1<RealAttr1> spLog (ds, false);  spLog. reserve (ds. attrs. size ());
	  Vector<Real> means;  means. reserve (ds. attrs. size ());
	  for (const string& attrName : attrNames)
	  {
	  	const Attr* attr_ = ds. name2attr (attrName);
	  	if (! attr_)
	  		throw runtime_error ("Attribute " + attrName + " is not found");
	  	const NumAttr1* attr = attr_->asNumAttr1 ();  // IntAttr1, >= 0
	  	ASSERT (attr);
	  	
	  	auto attr1 = new PositiveAttr1 (attr->name + suffix, ds, attr->decimals / 2 + 1);  // PAR
	  	spLog << attr1;
	  	FFOR (size_t, i, ds. objs. size ())
	  	{
	  		const Real x = attr->getReal (i);
	  	  if (! isNan (x))
	  	  {
	  	  	ASSERT (x >= 0);
	  	  	(*attr1) [i] = log (x + 1);  // >= 0
	  	  }
	  	}
	  	
	    const Sample sm (ds);
	  	Real average = NAN;
	  	Real scatter = NAN;
	  	attr1->getAverageScatter (sm, average, scatter);
	  	ASSERT (average > 0);
	  	means << average;
	  	FFOR (size_t, i, ds. objs. size ())
  	  	(*attr1) [i] /= average; 
	  }
	  ASSERT (spLog. size () == means. size ());
	  ds. qc ();
	  
	  
    auto averageAttr = new PositiveAttr1 ("average", ds, 6);  // PAR
    MVector weights (spLog. size (), 1.0);  // >= 0
    Progress prog;
    for (;;)
    {
	    // *averageAttr
	  	FFOR (size_t, i, ds. objs. size ())
	  	{
	  		(*averageAttr) [i] = 0;
	  		Real weight_sum = 0;
	  		FFOR (size_t, j, spLog. size ())
	  		{
	  			const Real x = (* spLog [j]) [i];
	  			if (isNan (x))
	  				continue;
	  			(*averageAttr) [i] += x * weights [j];
	  			weight_sum += weights [j];
	  		}
	 			(*averageAttr) [i] /= weight_sum;
	 		//sm. mult = ??
	  	}
	
	    // weights[]
	    const MVector weights_old (weights);
	 		FFOR (size_t, j, spLog. size ())
	 		{
	 			Real var = 0;
	 			Real mult_sum = 0;
		  	FFOR (size_t, i, ds. objs. size ())
		  	{
		 			const Real x = (* spLog [j])  [i];
		 			const Real y = (*averageAttr) [i];
		 			if (isNan (x))
		 				continue;
		 			if (isNan (y))
		 				continue;
		 			var += sqr (x - y) * 1;  // PAR
		 			mult_sum += 1;  // PAR
		 		}
		 		if (mult_sum == 0)
		 			throw runtime_error ("Attribute " + spLog [j] -> name + " has no data");
		 		ASSERT (mult_sum > 0);
		 		var /= mult_sum;
		 		ASSERT (var >= 0);
		 		weights [j] = 1 / var;
	 		}
	 		
	 		const Real diff = weights. maxAbsDifferenceVec (weights_old);
      prog (toString (diff));
	 		if (diff < 1e-6)  // PAR
	 			break;
	 	}
	 	
	 	
 		FFOR (size_t, j, spLog. size ())
 		{
	    string name (spLog [j] -> name);
	    EXEC_ASSERT (trimSuffix (name, suffix));
 		  cout         << name
 		       << '\t' << means [j]
 		       << '\t' << weights [j]
 		       << endl;
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



