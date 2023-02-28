// logreg.cpp

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
*   Logistic regression
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
    : Application ("Logistic regression")
  	{
  	  version = VERSION;
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
    const Space1<NumAttr1> sp (spRaw. toNumAttr1 (ds));
      
    const Sample sm (ds);
    
   
    LogisticRegression lr (sm, sp, *targetAttr);
    lr. qc ();
    section ("Solving", false);
    lr. solve ();
    lr. qc ();
    lr. saveText (cout);
          
    LogisticRegression* lrFinal = & lr;
    unique_ptr<LogisticRegression> lrSel;
    Space1<NumAttr1> spSel (sp);  // Is used by lrSel
    if (importance_min)
    {
      section ("Removing attributes", false);
      lrSel. reset (selectAttrs<LogisticRegression> (sm, spSel, *targetAttr, importance_min));
      if (lrSel. get ())
      { 
        lrFinal = lrSel. get ();
        lrFinal->qc ();
        cout << endl;
        cout << "Selected attributes:" << endl;
        lrFinal->saveText (cout);
      }
    }
    ASSERT (lrFinal);

    if (verbose ())
    {
      lr. getPredictionAttr ("hat");
      lr. getScoreAttr ("score");
      cout << endl;
      ds. saveText (cout);  
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
    lr0. saveText (cout);
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


