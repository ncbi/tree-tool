// attrs.cpp

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
*   Print statistics of attributes
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
    : Application ("Print the description of attributes")
    {
      version = VERSION;
      addPositional ("file", dmSuff + "-file");
      addKey ("corr_min", "Min. correlation between two attributes to report", "nan");
      addKey ("outlier_evalue_max", "Max. outlier e-value", "0.1");
      addFlag ("stat", "Compute statistics. For numeric attributes: mean, SD, left outlier threshold, right outlier threshold, mean w/o outliers, SD w/o outliers, normal distribution p-value");
    }
	
	
	
	void body () const final
	{
		const string inFName = getArg ("file");
		const Real corr_min  = str2real (getArg ("corr_min"));
		const Prob outlier_eValue_max = str2real (getArg ("outlier_evalue_max"));
		const bool stat      = getFlag ("stat");
		QC_ASSERT (outlier_eValue_max >= 0.0);


    const Dataset ds (inFName);
    const Sample sample (ds);
  //cout << "mult_sum: " << sample. mult_sum << endl << endl;
    
    cout << "AttrName\tDefinition\tMissings";
    if (stat)
      cout << "\tStatistics\tMain_interval";
    cout << endl;
    {
      Progress prog (ds. attrs. size (), 100);  // PAR
      for (const auto attrRaw : ds. attrs)
      {
        prog ();
      	cout << attrRaw->name << '\t' << attrRaw->getTypeStr () << '\t' << attrRaw->countMissings ();
      	if (stat)
      	{
          if (const NumAttr1* num = attrRaw->asNumAttr1 ())
    	    {
    	    	const ONumber on (cout, num->decimals + 1 + (attrRaw->asBoolAttr1 () ? 3 : 0), false);  // PAR
    	      Normal normal;
    	      const UniVariate<NumAttr1> an (sample, *num);
    	      normal. analysis = & an;
    	      normal. estimate ();
    	      normal. qc ();	      
    	      cout << '\t' << normal. loc 
    	           << '\t' << normal. scale;
    		    if (! isNan (outlier_eValue_max))
    		    {
    			    Real threshold [2] = {NaN, NaN};
    			    for (const bool rightTail : {false, true})
    			    {
    			    	threshold [rightTail] = num->locScaleDistr2outlier (sample, normal, rightTail, outlier_eValue_max);
    			      cout << '\t' << threshold [rightTail];
    			    }
    			    Sample sample_pure (ds);
    		      FFOR (size_t, i, sample_pure. mult. size ())
    		        if (   (*num) [i] > threshold [true]
    		        	  || (*num) [i] < threshold [false]
    		        	  || ! DM_sp::finite ((*num) [i])
    		        	 )
    		        	sample_pure. mult [i] = 0.0;
    		      const UniVariate<NumAttr1> an_pure (sample_pure, *num);
    		      normal. analysis = & an_pure;
    		      normal. estimate ();
              const Prob pVal = max (0.0, 1.0 - normal. getFitness_entropy ());
              if (isNan (pVal))
      		      cout << '\t' << NaN 
      		           << '\t' << NaN 
      		           << '\t' << NaN;
              else
      		      cout << '\t' << normal. loc 
      		           << '\t' << normal. scale
      		           << '\t' << pVal;
    			  }
    	    }
          else if (const BoolAttr1* boolean = attrRaw->asBoolAttr1 ())
    	      cout << '\t' << boolean->getProb (sample);
          else if (const NominAttr1* nomin = attrRaw->asNominAttr1 ())
          {
            unique_ptr<const Categorical> cat (nomin->getCategorical (sample));
    	      cout << "  ";
            cat->saveText (cout);
    	    }
    	  }
        cout << endl;
  	  //other types ??
  	  }
  	}   	 


    // Correlations
    // Only for NumAttr1 ??
    // Requires: no missings ??
    if (! isNan (corr_min))
    {
      cout << endl;
      const Space1<NumAttr1> space (ds, true);
      const MultiVariate<NumAttr1> an (sample, space);
      MultiNormal mn; 
      mn. qc ();
      mn. analysis = & an;
      mn. estimate ();
      mn. qc ();
      Progress prog (space. size (), 0);
      FOR (size_t, i, space. size ())
      {
        prog (space [i] -> name);
        FOR (size_t, j, i)
        {
        	const Real corr = mn. sigmaExact. get (false, i, j) / sqrt (mn. sigmaExact. getDiag (i) * mn. sigmaExact. getDiag (j));
          if (corr >= corr_min)
          	cout << "corr (" << space [i] -> name << ", " << space [j] -> name << ") = " << corr << endl;
        }
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


