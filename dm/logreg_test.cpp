// logreg_test.cpp

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
    : Application ("Test logistic regression")
  	{
  	  addFlag ("lin_dep", "There is a linear dependence between two input variables");
  	//addFlag ("perfect", "There is a separating hyperplane between the true and false objects");
  	//addPositional ("seed", "Seed for random numbers");
  	}
	
	
	
	void body () const final
	{
		const bool lin_dep = getFlag ("lin_dep");
  //const bool perfect = getFlag ("perfect");
	//const ulong seed   = str2<ulong> (getArg ("seed"));
		
		
		// PAR
		const size_t maxObjNum = 1000;  
		const size_t maxAttrNum = 3;  
    
    
    IMPLY (lin_dep, maxAttrNum >= 2);
    
 
    Normal norm;
    norm. setParam (0, 1);
    norm. setSeed (seed_global);

    Dataset ds;
    FOR (size_t, objNum, maxObjNum)
      ds. appendObj ("Obj" + toString (objNum + 1));    
    RealAttr1* a = ds. addRealAttr1Unit ();
    FOR (size_t, attrNum, maxAttrNum - 1/*unit*/)
    {
      const RealAttr1* aPrev = a;
      a = new RealAttr1 ("X" + toString (attrNum + 1), ds, 3);  // PAR
      FOR (size_t, objNum, maxObjNum)
        (*a) [objNum] = lin_dep && attrNum == maxAttrNum - 2 ? (*aPrev) [objNum] : norm. rand ();
    } 
    
    const Space1<NumAttr1> sp (ds, true);

    ExtBoolAttr1* target = new ExtBoolAttr1 ("Target", ds);  
    
    const Sample sm (ds);
    
    LogisticRegression lr (sm, sp, *target);
    lr. beta. putAll (1);
  #if 0
    if (perfect)  // --> test()
      for (Iterator it (sp); it ();)  
        (*target) [*it] = (ebool) (1 <= lr. predict (*it));  // PAR
    else
      lr. simulateTarget (seed_global);
  #endif

    lr. test (seed_global);
  
    if (verbose ())
    {
      cout << "P(target) = " << target->getProb (sm) << endl; 
      cout << "P(correct) = " << lr. getCorrectPredictionFrac () << endl;
      cout << "beta:" << endl;
      lr. beta. saveText (cout);
  
      Unverbose unv;
      if (verbose ())
      {
        lr. getPredictionAttr ("hat");
        ds. print (cout);  
      }
    }
  
    if (!  lin_dep)
      FOR (size_t, i, maxAttrNum)
        ASSERT (fabs (lr. beta [i] - 1) <= 0.5);  // PAR
  
    lr. setAttrImportance ();

  /*
    // JackKnifeAttr
    printf ("\nJack-Knife:\n");
    uint Converged;
    NUM_ATTR1* JackKnifeAttr = LogReg. JackKnife (Converged);
    printf ("# Converged = %u\n", Converged);
    delete JackKnifeAttr;
  */
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



