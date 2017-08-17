// linreg.cpp

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
    : Application ("Linear regression")
  	{
  	  addFlag ("intercept", "Use an intercep");
  	  addPositional ("file", dmSuff + "-file with Number attributes");
  	  addPositional ("target", "Target - Number attribute name");
  	}
	
	
	
	void body () const final
	{
		const bool intercept = getFlag ("intercept");
		const string inFName = getArg ("file");
		const string target  = getArg ("target");

    
    Dataset ds (inFName);
    if (intercept)
      ds. addRealAttr1Unit ();
    ds. qc ();
    
    const auto targetAttr = ds. name2attr (target) -> asRealAttr1 ();
    ASSERT (targetAttr);
    
    Space1<NumAttr1> sp (ds, true);
    sp. removeAttr (*targetAttr);

    const Sample sm (ds);
    
    L2LinearNumPrediction lr (sm, sp, *targetAttr);
    lr. qc ();
    lr. solveUnconstrained ();
    lr. qc ();
    cout << "Abs. Error = " << lr. absCriterion2Error () << endl;
    const Real constTarget = lr. getConstTarget ();
    cout << "Const target = " << constTarget << endl;
    cout << "Rel. Error = " << lr. getRelTargetCriterion (constTarget) << endl;

    cout << endl;
    cout << "beta: value" << endl;
    FOR (size_t, i, lr. space. size ())
      cout << lr. space [i] -> name << ": " << lr. beta [i] << endl;

  #if 0   // From LogisticRegression
    cout << endl;
    cerr << "Importance ..." << endl;
    cout << "beta: importance" << endl;
    lr. setAttrImportance ();
    FOR (size_t, i, lr. space. size ())
      cout << lr. space [i] -> name << ": " << lr. attrImportance [i] << endl;
      

    cout << endl;
    Space1 sp0 (ds, false);
    sp0 << unitAttr << targetAttr;
    LogisticRegression lr0 (sp0, *targetAttr);
    cerr << "Solving control ..." << endl;
    lr0. solve ();
    lr0. qc ();    
    cout << "-logL = " << lr0. negLogLikelihood << endl;
    cout << "P(Correct) = " << lr0. getCorrectPredictionFrac () << endl;
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


