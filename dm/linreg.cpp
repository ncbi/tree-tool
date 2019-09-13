// linreg.cpp

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
*   Linear regression
*
*/


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
    if (! targetAttr)
      throw runtime_error (strQuote (target) + " is not a Real attribute");
    
    Space1<NumAttr1> sp (ds, true);
    sp. removeAttr (*targetAttr);

    const Sample sm (ds);
    
    L2LinearNumPrediction lr (sm, sp, *targetAttr);
    lr. qc ();
  #if 1
    lr. solveUnconstrained ();
  #else  // for testing
  #if 0  
    FFOR (size_t, attrNum, lr. beta. size ())
      lr. beta [attrNum] = 1.0;  // PAR
  #else
    lr. beta [0] = 0.041179;
    lr. beta [1] = 2.21651;
    lr. beta [2] = 8.04923e-05;
    lr. setAbsCriterion ();
  #endif
    const bool solved = lr. solveUnconstrainedFast (nullptr, true/*??*/, 10, 0.01);  // PAR
    if (! solved)
      cout << "Not solved" << endl;
  #endif    
    lr. qc ();
    cout << "Abs. criterion = " << lr. absCriterion << endl;
    cout << "Abs. error = " << lr. absCriterion2Error () << endl;
    Real scatter = NaN;
    const Real constTarget = lr. getConstTarget (scatter);
    cout << "Const target = " << constTarget << endl;
    cout << "Rel. error = " << lr. getRelTargetCriterion (constTarget) << endl;

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


