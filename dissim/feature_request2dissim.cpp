// feature_request2dissim.cpp

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
*   Print dissimilarities computed from features for requested pairs of objects
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "evolution.hpp"
using namespace DM_sp;
#include "../version.inc"



namespace 
{
  
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print dissimilarities (Hamming distance, SNP dissimilarity, Jaccard index) computed from features for requested pairs of objects")
    {
      version = VERSION;
  	  addPositional ("req", "File with pairs of objects");
  	  addPositional ("dir", "Directory with <object> files containing feature-optional(0/1) pairs");
  	  addKey ("optional_weight", "For Hamming distance: Weight of optional-nonoptional match (0..1)", "NaN");
  	  addKey ("mutation_rate", "For SNP dissimilarity: File where each line has format: <feature_name> <mutation rate>");
  	}



	void body () const final
	{
		const string reqFName        =           getArg ("req");
		const string dir             =           getArg ("dir");
		const Real   optional_weight = str2real (getArg ("optional_weight"));
    const string mutationsFName  =           getArg ("mutation_rate");
		
		if (! isNan (optional_weight))
		{
  		QC_ASSERT (optional_weight >= 0.0);
  		QC_ASSERT (optional_weight <= 1.0);
  		QC_ASSERT (mutationsFName. empty ());
    }
		
		
		Vector<pair<string,Real>> feature2rate;  feature2rate. reserve (10000);  // PAR
		if (! mutationsFName. empty ())
    {
      LineInput f (mutationsFName);
      string name;
      Real rate;
      Istringstream iss;
      while (f. nextLine ())
      {
        iss. reset (f. line);
        rate = NaN;
        iss >> name >> rate;
        QC_ASSERT (! isNan (rate));
        feature2rate << move (pair<string,Real> (name, rate));
      }
      QC_ASSERT (! feature2rate. empty ());
    }


    LineInput reqF (reqFName);
    string obj1, obj2, rest;
    Istringstream iss;
    const ONumber on (cout, 6, true);  // PAR
    while (reqF. nextLine ())
    {
      trim (reqF. line);
      iss. reset (reqF. line);
      rest. clear ();
      iss >> obj1 >> obj2 >> rest;
      QC_ASSERT (! obj1. empty ());
      QC_ASSERT (! obj2. empty ());
      QC_ASSERT (rest. empty ());
      const FeatureVector vec1 (dir + "/" + obj1);
      const FeatureVector vec2 (dir + "/" + obj2);
      cout << obj1 << '\t' << obj2 << '\t' << features2dissim (vec1, vec2, optional_weight, feature2rate) << endl;
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



