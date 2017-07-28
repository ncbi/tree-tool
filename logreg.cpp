// logreg.cpp

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
#include "prediction.hpp"
using namespace DM_sp;



namespace
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Logistic regression")
  	{
  	  addKey ("importance_min", "Min. relative importance of attributes for attribute selection; 0 - no selection", "0");
  	  addPositional ("file", dmSuff + "-file with Number attributes");
  	  addPositional ("target", "Target - Boolean attribute name");
  	}
	
	
	
	void body () const final
	{
		const Real importance_min = str2real (getArg ("importance_min"));
		const string inFName      = getArg ("file");
		const string target       = getArg ("target");

    
    Dataset ds (inFName);
    RealAttr1* unitAttr = ds. addRealAttr1Unit ();
    ds. qc ();
    
    const BoolAttr1* targetAttr = ds. name2attr (target) -> asBoolAttr1 ();
    ASSERT (targetAttr);
    
    Space1<Attr1> spRaw (ds, true);
    spRaw. removeAttr (*targetAttr);
    const Common_sp::AutoPtr<const Space1<NumAttr1>> sp (spRaw. toNumAttr1 (ds));
      
    const Sample sm (ds);
    
   
    LogisticRegression lr (sm, *sp, *targetAttr);
    lr. qc ();
    cerr << "Solving..." << endl;
    lr. solve ();
    lr. qc ();
    lr. print (cout);
          
    LogisticRegression* lrFinal = & lr;
    Common_sp::AutoPtr<LogisticRegression> lrSel;
    Space1<NumAttr1> spSel (*sp);  // Is used by lrSel
    if (importance_min)
    {
      cerr << "Removing attributes..." << endl;
      lrSel. reset (selectAttrs<LogisticRegression> (sm, spSel, *targetAttr, importance_min));
      if (lrSel. get ())
      { 
        lrFinal = lrSel. get ();
        lrFinal->qc ();
        cout << endl;
        cout << "Selected attributes:" << endl;
        lrFinal->print (cout);
      }
    }
    ASSERT (lrFinal);

    if (verbose ())
    {
      lr. getPredictionAttr ("hat");
      lr. getScoreAttr ("score");
      cout << endl;
      ds. print (cout);  
    }

    cout << endl;
    cout << "beta: importance" << endl;
    lrFinal->setAttrImportance ();
    FOR (size_t, i, lrFinal->space. size ())
      cout << lrFinal->space [i] -> name << ": " << lrFinal->attrImportance [i] << endl;

    cout << endl;
    cout << "Control:" << endl;
    Space1<NumAttr1> sp0 (ds, false);
    sp0 << unitAttr;
    LogisticRegression lr0 (sm, sp0, *targetAttr);
    lr0. solve ();
    lr0. qc ();    
    lr0. print (cout);
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


