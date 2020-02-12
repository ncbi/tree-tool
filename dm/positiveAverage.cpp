// calibrateDissims.cpp  // --> dm/positiveAverage.cpp ??

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
*   Positive average model
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
#include "../dm/dataset.hpp"
using namespace DM_sp;
#include "version.inc"



namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Analyze proportional homoscedastic positive attributes. Print attributes statistics")
    {
      version = VERSION;
      // Input
  	  addPositional ("file", dmSuff + "-file without " + strQuote (dmSuff));
  	  addPositional ("power", "Power for raw dissimilarities to be raised to");
  	  addPositional ("outlierSEs", "Number of standard errors for an attribute value to be an outlier");  // explain more ??
  	  addFlag ("optimize_coeff", "Fit the coefficients of the attributes to the average by intercept-free linear regression");
  	  addKey ("iter_max", "Max. number of optimization iterations; 0 - infinity", "0");
  	  // Output
  	  addKey ("output_dissim", "Output merged dissimilarity");
  	}



	void body () const final
	{
		const string inFName            = getArg ("file");
		const Real power                = str2real (getArg ("power"));
		const Real outlierSEs           = str2real (getArg ("outlierSEs"));
		const bool optimize_coeff       = getFlag ("optimize_coeff");
		const size_t iter_max           = str2<size_t> (getArg ("iter_max"));
		const string output_dissimFName = getArg ("output_dissim");
		
		if (outlierSEs <= 0.0)
		  throw runtime_error ("outlierSEs must be positive");
		

    Dataset ds (inFName);
    ds. qc ();
    
    const Sample sm (ds);
	    

    Space<PositiveAttr1> space (ds, false);  space. reserve (ds. attrs. size ());
	  for (const Attr* attr_ : ds. attrs)
	  {
	  	const PositiveAttr1* attr = attr_->asPositiveAttr1 (); 
	  	if (! attr)
	  		throw runtime_error (attr_->name + " should be Positive");
	  		
	    var_cast (attr) -> inf2missing ();
	    
	    if (attr->missingsAll ())
	    {
	      cout << attr->name << ": all values are missing" << endl;
	      continue;
	    }
	  		  	  
	    Real attr_min, attr_max;
	    attr->getMinMax (sm, attr_min, attr_max);
	    ASSERT (attr_min >= 0.0);
	    if (attr_max == 0.0)
	    {
	      cout << attr->name << " is 0.0" << endl;
	      continue;
	    }

  	  FFOR (size_t, objNum, ds. objs. size ())
      {
      	const Real x = (*attr) [objNum];
        var_cast (*attr) [objNum] = pow (x, power);
      }

      space << attr;
	  }
  	if (space. empty ())
  	  throw runtime_error ("No attributes");
	  
	  
	  PositiveAverage pa (sm, space, outlierSEs, optimize_coeff);  
    pa. calibrate (iter_max);
	  pa. qc ();  
   	pa. saveText (cout);   	
    cerr << "Effective number of attributes: " << pa. model. getEffectiveAttrs () << endl;
	  

	 	if (! output_dissimFName. empty ())
	 	{
	 		OFStream f (output_dissimFName + dmSuff);
	 	  sm. save (VectorPtr<Attr> {pa. averageAttr}, f);
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



