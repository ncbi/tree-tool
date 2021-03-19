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
  

void vec2dna_prot (const ObjFeatureVector &vec, 
                   ObjFeatureVector &dnaVec, 
                   ObjFeatureVector &protVec)
// Output: dnaVec, protVec
{
  ASSERT (dnaVec. empty ());
  ASSERT (protVec. empty ());
  dnaVec.  reserve (vec. size());
  protVec. reserve (vec. size());
  for (const ObjFeature& of : vec)
    if (contains (of. name, '-'))  // Protein Mutation's have geneName
      protVec << of;
    else
      dnaVec << of;
  dnaVec. sort ();
  ASSERT (dnaVec. isUniq ());
  protVec. sort ();
  ASSERT (protVec. isUniq ());  
}


  
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print dissimilarities (Hamming distance, SNP dissimilarity, Jaccard index) computed from features for requested pairs of objects")
    {
      version = VERSION;
  	  addPositional ("req", "File with pairs of objects");
  	  addPositional ("dir", "Directory with <object> files containing feature-optional(0/1) pairs");
  	  addKey ("optional_weight", "For Hamming distance: Weight of optional-nonoptional match (0..1)", "NaN");
  	  addKey ("mutation_rate", "For SNP dissimilarity: file where each line has format: <feature_name> <mutation rate 0->1> <mutation rate 1->0>");
  	//addFlag ("freq", "<mutation_rate> file has allele frequencies");
  	  addKey ("dna_weight", "Weight of DNA mutations (>=0). Protein mutations are weighted as 1. Features are DNA or protein mutations, where DNA is not named and proteins are named", "NaN");
  	  addFlag ("virus", "Features are DNA and protein mutations, protein mutations are ignored");
  	}



	void body () const final
	{
		const string reqFName        =           getArg ("req");
		const string dir             =           getArg ("dir");
		const Prob   optional_weight = str2real (getArg ("optional_weight"));
    const string mutationsFName  =           getArg ("mutation_rate");
  //const bool   freq            =           getFlag ("freq");
    const Prob   dnaWeight       = str2real (getArg ("dna_weight"));
    const bool   virus           =           getFlag ("virus");
		
		if (! isNan (optional_weight))
		{
  		QC_ASSERT (isProb (optional_weight));
  		QC_ASSERT (mutationsFName. empty ());
    }    
        
    if (! isNan (dnaWeight))
    { 
      QC_ASSERT (dnaWeight >= 0.0);
      QC_ASSERT (mutationsFName. empty ());
      QC_ASSERT (! virus);
    }
		
		
		unique_ptr<const FeatureVector> feature2rate;
		if (! mutationsFName. empty ())
		{
 		  feature2rate. reset (new FeatureVector (mutationsFName));
      QC_ASSERT (! feature2rate->empty ());
    }
    else
 		  feature2rate. reset (new FeatureVector ());


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
    #if 0
      // ??
      const ObjFeatureVector& vec1 = ObjFeatureVector::getCache (dir + "/" + obj1);
      const ObjFeatureVector& vec2 = ObjFeatureVector::getCache (dir + "/" + obj2);
    #else
      const ObjFeatureVector vec1 (dir + "/" + obj1);
      const ObjFeatureVector vec2 (dir + "/" + obj2);
    #endif
      Real dissim = NaN;
      if (isNan (dnaWeight) && ! virus)
        dissim = features2dissim (vec1, vec2, optional_weight, *feature2rate);
      else
      {
        ObjFeatureVector dnaVec1, protVec1;
        ObjFeatureVector dnaVec2, protVec2;
        vec2dna_prot (vec1, dnaVec1, protVec1);
        vec2dna_prot (vec2, dnaVec2, protVec2);
        const Real dnaDissim  = features2hamming (dnaVec1,  dnaVec2,  optional_weight);
        if (virus)
          dissim = dnaDissim;
        else
        {
          const Real protDissim = features2hamming (protVec1, protVec2, optional_weight);
          dissim = dnaWeight * dnaDissim + protDissim;  
        }
      }
      ASSERT (dissim >= 0.0);
      cout << obj1 << '\t' << obj2 << '\t' << dissim << endl;
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



