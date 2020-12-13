// logreg_test.cpp

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
*   Test of logistic regression
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
#include "prediction.hpp"
using namespace DM_sp;
#include "../version.inc"



namespace
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Test logistic regression")
  	{
  	  version = VERSION;
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
        ds. saveText (cout);  
      }
    }
  
    if (!  lin_dep)
      FOR (size_t, i, maxAttrNum)
        { QC_ASSERT (fabs (lr. beta [i] - 1) <= 0.5); }  // PAR
  
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



