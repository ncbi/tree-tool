// evolution.cpp

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


#undef NDEBUG
#include "../common.inc"

#include "evolution.hpp"

#include "../dm/optim.hpp"




namespace DM_sp
{


Real maps2dissim (Real size1,
                  Real size2,
                  Real maps1,
                  Real maps2,
                  Real maps_min,
                  Prob sizes_ratio_min,
                  bool ave_arithP)
{
	ASSERT (size1 >= 0.0);
	ASSERT (size2 >= 0.0);
  ASSERT (maps1 <= size1);
  ASSERT (maps2 <= size2);
  ASSERT (isProb (sizes_ratio_min));
  
  
  if (verbose ())
    cout << '\t' << maps1 
         << '\t' << maps2
         << '\t' << size1 
         << '\t' << size2 
         << endl;
         
  if (! max (size1, size2))
  	return NaN;

  if (  min (size1, size2) 
  	  / max (size1, size2)
  	    < sizes_ratio_min
  	 )
  	return NaN;
         
  if (maps1 < maps_min)
  	return NaN;
  if (maps2 < maps_min)
  	return NaN;
  	
  const Real dissim1 = log (size1 / maps1);
  const Real dissim2 = log (size2 / maps2);
  ASSERT (dissim1 >= 0.0);
  ASSERT (dissim2 >= 0.0);
  
  if (ave_arithP)
    return ave_arith (dissim1, dissim2);
  //return ave_geom  (dissim1, dissim2);
  else
    return ave_harm  (dissim1, dissim2);
}




// Hashes

Hashes::Hashes (const string &fName)
{
  reserve (10000);  // PAR

  LineInput hf (fName);
  size_t prev = 0;
  while (hf. nextLine ())
  {
    const size_t hash = str2<size_t> (hf. line);
    ASSERT (hash);
    if (hash <= prev)
      throw runtime_error ("Hash " + hf. line + " is not greater than the previous hash " + toString (prev));
    *this << hash;
    prev = hash;
  } 
  searchSorted = true;
/*
  if (empty ())
    throw runtime_error ("No hashes for " + fName);
*/
}



// Read Hashes from a binary file ??




// ObjFeatureVector

ObjFeatureVector::ObjFeatureVector (const string &fName)
{
  reserve (16);  // PAR

  LineInput f (fName);
  while (f. nextLine ())
    (*this) << ObjFeature (std::move (f. line));
    
  sort ();
  QC_ASSERT (isUniq ());
}



map<string,ObjFeatureVector> ObjFeatureVector::cache;



const ObjFeatureVector& ObjFeatureVector::getCache  (const string &fName)
// bad_alloc ??
{
  if (const ObjFeatureVector* fv = findPtr (cache, fName))
    return *fv;

  ObjFeatureVector fv (fName);
  return (cache [fName] = std::move (fv));
}



namespace
{

Real features2hamming_half (const ObjFeatureVector &vec1,
                            const ObjFeatureVector &vec2,
                            Real optional_weight)
{
	ASSERT (optional_weight >= 0.0);
	ASSERT (optional_weight <= 1.0);

  Real n = 0.0;
  for (const ObjFeature& f : vec2)
  {
    const size_t index = vec1. binSearch (f);
    if (f. optional)
    {
      if (index != no_index && ! vec1 [index]. optional)
        n += optional_weight;
    }
    else
      if (index == no_index)
        n += 1.0;
      else
        if (vec1 [index]. optional)
          n += optional_weight;
  }
  return n;
}

}




Real features2hamming (const ObjFeatureVector &vec1,
                       const ObjFeatureVector &vec2,
                       Real optional_weight)
{
  const Real w = optional_weight / 2.0;
  return   features2hamming_half (vec1, vec2, w)
         + features2hamming_half (vec2, vec1, w);
}



Real features2jaccard (const ObjFeatureVector &vec1,
                       const ObjFeatureVector &vec2)
{
  const size_t intersection = vec1. getIntersectionSize (vec2);
  ObjFeatureVector f (vec2);
  f << vec1;
  f. sort ();
  f. uniq ();
  const Prob jaccard = (Real) intersection / (Real) f. size ();
  ASSERT (isProb (jaccard));
  return - log (jaccard);
}




// Feature

Feature::Feature (const string &line)
{ 
  istringstream iss (line);
  iss >> name >> lambda [0] >> lambda [1];
  for (const bool b : {false, true})
    QC_ASSERT (lambda [b] >= 0.0);
}




// FeatureVector

FeatureVector::FeatureVector (const string &fName)
{
  reserve (1000);  // PAR

  LineInput f (fName);
  while (f. nextLine ())
  {
    Feature feature (f. line);
    if (! feature. redundant ())
    {
      (*this) << std::move (feature);
      ASSERT (feature. name. empty ());
    }
  }
    
  sort ();
  QC_ASSERT (isUniq ());
}



//

struct Func_snps2time : Func1
{
  struct Snp 
  {
    array<Real,2/*bool*/> lambda;
      // > 0
    ebool alleles;
      // enull <=> heterozygous
    Real sum () const
      { return   lambda [0] 
               + lambda [1]; 
      }
  };
  Vector<Snp> snps;
  Real lambdaSum_ave {NaN};
  array<size_t,2/*bool*/> num;
  bool same {true};
  

  Func_snps2time (const ObjFeatureVector &vec1,
                  const ObjFeatureVector &vec2,
                  const FeatureVector &feature2lambda)
    {
      snps. reserve (max (vec1. size (), vec2. size ()));
      MeanVar mv;
      for (const bool b : {false, true})
        num [b] = 0;
      Real lambda = NaN;
      for (const Feature& f : feature2lambda)
      {
        // Make time O(1) ??
        const ObjFeature of {f. name, false};
        const size_t i1 = vec1. binSearch (of);
        if (i1 != no_index && vec1 [i1]. optional)
          continue;
        const size_t i2 = vec2. binSearch (of);
        if (i2 != no_index && vec2 [i2]. optional)
          continue;
        const ebool alleles = (i1 != no_index && i2 != no_index
                                 ? etrue
                                 : i1 == no_index && i2 == no_index
                                   ? efalse
                                   : enull
                              );
        snps << Snp {f. lambda, alleles};
        mv << f. lambdaSum ();
        num [(i1 == no_index) == (i2 == no_index)] ++;
        if (same)
        {
          for (const bool b : {false, true})
            if (isNan (lambda))
              lambda = f. lambda [b];
            else
              if (lambda != f. lambda [b])
                same = false;
        }
      }
      lambdaSum_ave = mv. getMean ();
      QC_ASSERT (lambdaSum_ave > 0.0);
    }
  

  Real f (Real t) final
    // Non-monotone
    {
      ASSERT (t >= 0.0);
      if (t == 0.0)
        return -inf;
      Real s = 0.0;  
      for (const Snp& snp : snps)
      {
        const Real lambdaSum = snp. sum ();
        const Real a = exp (t * lambdaSum);
        ASSERT (a > 1.0);
        const Real b = snp. alleles == enull
                         ? - lambdaSum / (a - 1.0)
                         : lambdaSum / ((lambdaSum / snp. lambda [(bool) snp. alleles] - 1.0) * a + 1.0);  
        s += b;
      }
      return s;
    }  
    
  Real f_approx () const
    // Equivalent to Jukes-Cantor model
    { 
    #if 0  // ??
      const size_t diff = num [true] <= num [false] ? 1 : (num [true] - num [false]);
      return - 1.0 / lambdaSum_ave * log ((Real) diff / (Real) (num [true] + num [false]));
    #else
      return num [true] <= num [false]
               ? inf
               : - 1.0 / lambdaSum_ave * log ((Real) (num [true] - num [false]) / (Real) (num [true] + num [false]));
    #endif
    }
};



Real snps2time (const ObjFeatureVector &vec1,
                const ObjFeatureVector &vec2,
                const FeatureVector &feature2lambda)
{
  ASSERT (! feature2lambda. empty ());
  
  Func_snps2time f (vec1, vec2, feature2lambda);

  const Real dissim_init = f. f_approx ();
  if (f. same)
    return dissim_init;
  if (dissim_init == inf)
    return inf;
    
  const Real delta = dissim_init * 1e-7;  // PAR
  const Real dissim = f. findZeroPositive (dissim_init, delta, 1.1);   // PAR
  IMPLY (! isNan (dissim), dissim >= 0.0);
  
  return dissim;
}



}

