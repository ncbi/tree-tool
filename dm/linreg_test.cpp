// linreg_test.cpp

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
    : Application ("Test linear regression")
  	{
  	  addFlag("lin_dep", "There is a linear dependence bwteen two input variables");
  	//addFlag("perfect", "There is a separating hyperplane between the true and false objects");
  	  addPositional ("seed", "Seed for random numbers");
  	}
	
	
	
	void body () const final
	{
		const bool lin_dep = getFlag ("lin_dep");
  //const bool perfect = getFlag ("perfect");
		const ulong seed   = str2<ulong> (getArg ("seed"));
		
		
		// PAR
		const size_t maxObjNum = 1000;  
		const size_t maxAttrNum = 3;  
    
    
    IMPLY (lin_dep, maxAttrNum >= 2);
    
 
    Normal norm;
    norm. setParam (0, 1);
    norm. setSeed (seed);

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

    auto target = new RealAttr1 ("Target", ds);  
    ds. qc ();
    
    const Sample sm (ds);
    
    {
      L2LinearNumPrediction lr (sm, sp, *target);
      lr. beta. putAll (1);
      FOR (size_t, i, ds. objs. size ())
        (*target) [i] = lr. predict (i) + norm. rand () / 10;  // PAR
    //lr. setMult ();
      
      lr. solveUnconstrained ();
      lr. qc ();
      if (verbose ())
      {
        cout << "Error = " << lr. absCriterion2Error () << endl;
        cout << "beta:" << endl;
        lr. beta. saveText (cout);
      }
      if (! lin_dep)
        FOR (size_t, i, lr. beta. size ())
          ASSERT (abs (lr. beta [i] - 1) <= 0.05);  // PAR
    }
  
  
    {
      L2LinearNumPrediction lr1 (sm, sp, *target);
      lr1. qc ();
      lr1. beta. putAll (0);  // PAR
      lr1. solveUnconstrainedAlternate (nullptr, false, 10, 0.01); 
      lr1. qc ();
      if (verbose ())
      {
        cout << endl;
        cout << "Error = " << lr1. absCriterion2Error () << endl;
        cout << "beta:" << endl;
        lr1. beta. saveText (cout);
      }
      if (! lin_dep)
        FOR (size_t, i, lr1. beta. size ())
          ASSERT (abs (lr1. beta [i] - 1) <= 0.05);  // PAR    
    }
    
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



