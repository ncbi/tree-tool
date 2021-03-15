// evolution.hpp

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
*   Evolution utilities
*
*/


#ifndef EVOLUTION_HPP
#define EVOLUTION_HPP

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
using namespace DM_sp;



namespace DM_sp
{


Real intersection2dissim (Real size1,
                          Real size2,
                          Real intersection,
                          Real intersection_min,
	                        Prob sizes_ratio_min,
	                        bool ave_arithP);
  // Input: ave_arithP; false <=> ave_harm()
  // Return: >= 0; may be NaN
  // Symmetric



struct Hashes : Vector<size_t>
// searchSorted
{
	explicit Hashes (const string &fName);
	Hashes () = default;
	
	Real getDissim (const Hashes &other,
	                size_t intersection_min,
	                Prob hashes_ratio_min) const
		{ return intersection2dissim ( (Real) size ()
			                           , (Real) other. size ()
			                           , (Real) getIntersectSize (other)
			                           , (Real) intersection_min
			                           , hashes_ratio_min
			                           , false  // PAR
			                           ); 
		}
};




// ObjFeature

struct ObjFeature : Named
// Feature of an object
{
  bool optional {false};
  

  explicit ObjFeature (string &&line)
    : Named (move (line))
    { replace (name, '\t', ' ');
      trim (name);
      if (trimSuffix (name, " 0"))
        ;
      else if (trimSuffix (name, " 1"))
        optional = true;
      trim (name);
    }
  ObjFeature (const string &name_arg,
              bool optional_arg)
    : Named (name_arg)
    , optional (optional_arg)
    {}
  ObjFeature () = default;
  
    
  bool operator< (const ObjFeature &other) const
    { return name < other. name; }
  bool operator== (const ObjFeature &other) const
    { return name == other. name; }
};



struct ObjFeatureVector : Vector<ObjFeature>
// sort()'ed, uniq
{
  explicit ObjFeatureVector (const string &fName);
  ObjFeatureVector () = default;

private:
  static map<string,ObjFeatureVector> cache;
public:
  static const ObjFeatureVector& getCache (const string &fName);
};



Real features2hamming (const ObjFeatureVector &vec1,
                       const ObjFeatureVector &vec2,
                       Real optional_weight);

Real features2jaccard (const ObjFeatureVector &vec1,
                       const ObjFeatureVector &vec2);
  // ObjFeature::optional is trated as false


struct Feature : Named
// Boolean attribute
{
  array<Real,2/*bool*/> lambda;
    // >= 0
  

  explicit Feature (const string &line);
    // Format: <name> <rate> <p>
  Feature () = default;
  
  
  bool redundant () const
    { return    ! lambda [0] 
             || ! lambda [1];
    }
  Real lambdaSum () const
    { return lambda [0] + lambda [1]; }
    
  bool operator< (const Feature &other) const
    { return name < other. name; }
  bool operator== (const Feature &other) const
    { return name == other. name; }
};



struct FeatureVector : Vector<Feature>
// sort()'ed, uniq
// !Feature::redundant()
{
  explicit FeatureVector (const string &fName);
  FeatureVector () = default;
};



Real snps2time (const ObjFeatureVector &vec1,
                const ObjFeatureVector &vec2,
                const FeatureVector &feature2lambda);

inline Real features2dissim (const ObjFeatureVector &vec1,
                             const ObjFeatureVector &vec2,
                             Real optional_weight,
                             const FeatureVector &feature2lambda)
  { if (! isNan (optional_weight))
      return features2hamming (vec1, vec2, optional_weight);
    if (feature2lambda. empty ())
      return features2jaccard (vec1, vec2);
    return snps2time (vec1, vec2, feature2lambda);
  }
         

}  // namespace



#endif


