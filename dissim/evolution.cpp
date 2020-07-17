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


Real intersection2dissim (Real size1,
                          Real size2,
                          Real intersection,
                          Real intersection_min,
	                        Prob sizes_ratio_min,
	                        bool ave_arithP)
{
	ASSERT (size1 >= 0);
	ASSERT (size2 >= 0);
  ASSERT (intersection <= size1);
  ASSERT (intersection <= size2);
  ASSERT (isProb (sizes_ratio_min));
  
  
  if (verbose ())
    cout << '\t' << intersection 
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
         
  if (intersection < intersection_min)
  	return NaN;
  	
  const Real dissim1 = log (size1 / intersection);
  const Real dissim2 = log (size2 / intersection);
  ASSERT (dissim1 >= 0);
  ASSERT (dissim2 >= 0);
  
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




// Feature

FeatureVector::FeatureVector (const string &fName)
{
  reserve (16);  // PAR

  LineInput f (fName);
  while (f. nextLine ())
    (*this) << Feature (move (f. line));
    
  sort ();
  QC_ASSERT (isUniq ());
}




namespace
{

Real features2hamming_half (const FeatureVector &vec1,
                            const FeatureVector &vec2,
                            Real optional_weight)
{
	ASSERT (optional_weight >= 0.0);
	ASSERT (optional_weight <= 1.0);

  Real n = 0.0;
  for (const Feature& f : vec2)
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




Real features2hamming (const FeatureVector &vec1,
                       const FeatureVector &vec2,
                       Real optional_weight)
{
  return   features2hamming_half (vec1, vec2, optional_weight)
         + features2hamming_half (vec2, vec1, optional_weight);
}



Real features2jaccard (const FeatureVector &vec1,
                       const FeatureVector &vec2)
{
  const size_t intersection = vec1. getIntersectSize (vec2);
  FeatureVector f (vec2);
  f << vec1;
  f. sort ();
  f. uniq ();
  const Prob jaccard = (Real) intersection / (Real) f. size ();
  ASSERT (isProb (jaccard));
  return - log (jaccard);
}



namespace
{

struct Func : Func1
{
  struct Snp 
  {
    const Real rate;
      // > 0
    const bool same; 
  };
  Vector<Snp> snps;
  

  Func (const FeatureVector &vec1,
        const FeatureVector &vec2,
        const Vector<pair<string,Real>> &feature2rate)
    {
      snps. reserve (max (vec1. size (), vec2. size ()));
      for (const pair<string,Real>& f2r : feature2rate)
      {
        QC_ASSERT (f2r. second >= 0.0);
        if (! f2r. second)
          continue;
        const Feature f {f2r. first, false};
        const size_t i1 = vec1. binSearch (f);
        if (i1 != no_index && vec1 [i1]. optional)
          continue;
        const size_t i2 = vec2. binSearch (f);
        if (i2 != no_index && vec2 [i2]. optional)
          continue;
        snps << Snp {f2r. second, (i1 == no_index) == (i2 == no_index)};
      }
    }
  

  Real f (Real t) final
    {
      ASSERT (t >= 0.0);
      if (t == 0.0)
        return -INF;
      Real s = 0.0;
      for (const Snp& snp : snps)
      {
        const Real a = exp (- 2.0 * t * snp. rate);
        ASSERT (a >= 0.0);
        ASSERT (a < 1.0);
        const Prob p = snp. same 
                         ? a * snp. rate / (1.0 + a)
                         : - snp. rate / (1.0 / a - 1.0);  // = - a * snp. rate (1.0 - a)
        s += p;
      }
      return s;
    }  
};

}



Real snps2time (const FeatureVector &vec1,
                const FeatureVector &vec2,
                const Vector<pair<string,Real>> &feature2rate)
{
  ASSERT (! feature2rate. empty ());
  
  Func f (vec1, vec2, feature2rate);

  size_t same = 0;
  size_t diff = 0;
  for (const Func::Snp& snp : f. snps)
    if (snp. same)
      same++;
    else
      diff++;
  if (same < diff)
    return INF;

  const Real delta = 1e-7;  // PAR
  const Real dissim = f. findZero (0.0, 1.0, delta);  // 1.0 = tree length
  ASSERT (dissim >= 0.0);
  ASSERT (dissim <= 1.0);

  if (dissim >= 1.0 - delta)
    return INF;
  return dissim;
}



}


