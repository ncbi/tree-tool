// analyzePositives.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
#include "prediction.hpp"
using namespace DM_sp;



namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Analyze proportional positive attributes: logarithmize, print mean and variance")
    {
      // Input
  	  addPositional ("file", dmSuff + "-file without \"" + dmSuff + "\"");
  	  addPositional ("attrs", "List of 2-way <attr> positive attributes in <file>");
  	  addFlag ("predict", "Predict each attribute by linear regression, print the coefficients of predictors");
  	}



	void body () const final
	{
		const string inFName    = getArg ("file");
		const string attrsFName = getArg ("attrs");
		const bool predict      = getFlag ("predict");


   	const StringVector attrNames (attrsFName, 256);  // PAR
    
    Dataset ds (inFName);
    ds. qc ();

    const string suffix ("_log");

	  Space1<RealAttr1> spLog (ds, false);
	  spLog << ds. addRealAttr1Unit ();  
	  for (const string& attrName : attrNames)
	  {
	  	const Attr* attr_ = ds. name2attr (attrName);
	  	if (! attr_)
	  		throw runtime_error ("Attribute " + attrName + " is not found");
	  	const NumAttr1* attr = attr_->asNumAttr1 ();
	  	ASSERT (attr);
	  	auto attr1 = new RealAttr1 (attr->name + suffix, ds, attr->decimals / 2 + 1);  // PAR
	  	spLog << attr1;
	  	FFOR (size_t, i, ds. objs. size ())
	  	{
	  		const Real x = attr->getReal (i);
	  	  if (! isNan (x))
	  	  {
	  	  	ASSERT (x >= 0);
	  	  	(*attr1) [i] = log (x + 1);
	  	  }
	  	}
	  }
	  

    Progress prog ((uint) spLog. size () - 1);
    for (const RealAttr1* target : spLog)
    {    
      if (target == spLog. front ())  // unit  
        continue; 
      prog (target->name);
    		
	    Space1<RealAttr1> sp (spLog);
	    sp. removeAttr (*target);
	
	    Sample sm (ds);
	    FFOR (size_t, i, ds. objs. size ())
	      sm. mult [i] = target->isMissing (i) ? 0 : 1 /*(1 / (*target) [i])*/;
	    sm. finish ();
	    
	    L2LinearNumPrediction lr (sm, sp, *target);
	    lr. qc ();

	    Real targetVar = NAN;
	    const Real targetMean = lr. getConstTarget (targetVar);
	    ASSERT (targetVar >= 0);
	    string name (target->name);
	    EXEC_ASSERT (trimSuffix (name, suffix));
	    cout         << name 
	         << '\t' << targetMean
	         << '\t' << targetVar;

      if (predict)
      {
		    lr. solveUnconstrained ();
		    ASSERT (! isNan (lr. absCriterion));
		    lr. qc ();
	      cout << '\t' << lr. absCriterion2Error ()
	           << '\t' << lr. getRelTargetCriterion (targetMean); 
		    size_t betaNum = 1;
		    for (const RealAttr1* predictor : spLog)
		    	if (predictor != spLog. front ())  // unit
		    	{
		    		cout << '\t';
		    		if (predictor != target)
		    		{
			    		cout << lr. beta [betaNum];
		    		  betaNum++;
		    		}
		    	}
		    ASSERT (betaNum == lr. beta. size ());	    		
		  }

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



