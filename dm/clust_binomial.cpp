// clust.cpp

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
*   Clustering of binomial attributes
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
using namespace DM_sp;
#include "../version.inc"



namespace
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Clustering as the decomposition of a mixture of a binomial (with p ~ 1) and uniform distributions")
  	{
  	  version = VERSION;
  	  addPositional ("file", dmSuff + "-file");	 
  	  addPositional ("attr", "name of a non-negative integer attribute"); 
  	  addPositional ("bin_n", "Binomial N parameter");
  	  addKey ("outlier_pValue", "Binomial outlier p-value", "0");
  	}
	
	
	
	void body () const final
	{
		const string inFName      = getArg ("file");
		const string attrName     = getArg ("attr");
		const int bin_n           = str2<int> (getArg ("bin_n"));
    const Prob outlier_pValue = str2real (getArg ("outlier_pValue"));
		ASSERT (bin_n >= 1);
		ASSERT (isProb (outlier_pValue));
		

    Dataset ds (inFName);
    const IntAttr1* attr = nullptr;
    {
    	const Attr* attr_ = ds. name2attr (attrName);
	    if (! attr_)
	    	throw runtime_error ("Attribute " + attrName + " is not found");
	    attr = attr_->asIntAttr1 ();
	    if (! attr)
	    	throw runtime_error ("Attribute " + attrName + " is not integer");
	  }


    const Sample sm (ds);    	    	
    const UniVariate<IntAttr1> an (sm, *attr);
    
    Mixture mixt;

    constexpr Prob bin_prob = 0.9;  // PAR
    ASSERT (isProb (bin_prob));

    auto bin_good = new Binomial;
    bin_good->analysis = & an;
    bin_good->setParam (bin_n, 0.9);  // PAR
    mixt. addComponent (bin_good, bin_prob);  
    
  #if 0
    auto bin_bad = new Binomial;
    bin_bad->analysis = & an;
    bin_bad->setParam (bin_n, 0.5);  // PAR
    mixt. addComponent (bin_bad, (1 - bin_prob) / 2);  
  #endif

    auto ud = new UniformDiscrete;
    ud->analysis = & an;
    ud->setParam (0, bin_n);
    mixt. addComponent (ud, (1 - bin_prob) / 2);  
    
  //mixt. qc ();
    mixt. estimate ();
    mixt. qc ();
    mixt. print (cout);
    
    if (outlier_pValue)
	    FOR_REV (int, val, bin_n)
	      if (bin_good->cdf (val) <= outlier_pValue)
	      {
	        cout << "Max. outlier threshold: " << val << endl;
	      	break;
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


