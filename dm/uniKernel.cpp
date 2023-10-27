// uniKernel.cpp

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
*   Estimate a Kernel distribution
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "numeric.hpp"
#include "dataset.hpp"
using namespace DM_sp;
#include "../version.inc"



namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
  : Application ("Kernel p.d.f. estimation. Parameters: (prob_uniform, half_window). Print distribution with the step = window")
  { 
    version = VERSION;
	  addPositional ("file", dmSuff + "-file without the extension");
	  addPositional ("attrName", "name of a real attribute");
	  addKey ("window", "half-Window of UniKernel; 0 <=> default", "0");  
  }


	void body () const final
	{
		const string fName    = getArg ("file");  
		const string attrName = getArg ("attrName");
		const Real   window   = str2real (getArg ("window"));
		
		
    Dataset ds (fName);
    ds. qc ();
    
    const Sample sm (ds);
    
    const NumAttr1* attr = nullptr;
    {
      const Attr* attr_ = ds. name2attr (attrName);
      ASSERT (attr_);
      attr = attr_->asNumAttr1 ();
    }
    ASSERT (attr);
    
    const UniVariate<NumAttr1> analysis (sm, *attr);
    
    UniKernel uniKernel (analysis);
    if (window)
      uniKernel. setParam (window);
    else
      uniKernel. estimate ();
    uniKernel. qc ();
    cout << uniKernel. nameParam () << endl;    
    cout << "Mode: " << (*attr) [uniKernel. getModeObjNum ()] << endl;

    if (uniKernel. halfWindow)
    {
      const Real step = 2 * uniKernel. halfWindow;  // PAR
      Real x = uniKernel. attr_min - step;
      while (x <= uniKernel. attr_max + step)
      {
        cout << x << "\t" << uniKernel. pdf (x) << endl;
        x += step;
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



