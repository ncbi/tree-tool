// numeric.cpp

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
*   Numeric utilities
*
*/


#undef NDEBUG
#include "../common.inc"

#include "numeric.hpp"

#include <iostream>



namespace DM_sp
{
 

const Real epsilon = sqrt (numeric_limits<Real>::epsilon ());  // = 1.49012e-08  
const Real pi = acos (-1.0);
Real log_2  = NaN;
Real log_10 = NaN;
Real sqrt_2 = NaN;



bool initNumeric ()
{
  MODULE_INIT

  initCommon ();

#ifdef _MSC_VER
  #pragma warning (disable : 4127)
#endif
  // Real
  ASSERT (numeric_limits<Real>::is_specialized);
  ASSERT (! numeric_limits<Real>::is_integer);
  ASSERT (numeric_limits<Real>::is_signed);
  ASSERT (numeric_limits<Real>::has_infinity);
  ASSERT (numeric_limits<Real>::has_quiet_NaN);

  ASSERT (INF / 2 == INF);
  ASSERT (INF > 0);

  ASSERT (isNan (NaN));
  ASSERT (isNan (INF - INF));
  ASSERT (isNan (INF / INF));
#ifdef _MSC_VER
  #pragma warning (default : 4127)
#endif

#ifndef _MSC_VER
  const Real zero = 0;
  ASSERT (isNan (zero / zero));
#endif

  log_2  = log (2);
  log_10 = log (10);
  sqrt_2 = sqrt (2);
  
  // Initialize local static variables
  lnFactorial (200 - 1);
  
  return true;
}


namespace
{
  const bool init_ = initNumeric ();
} 
 


string real2str (Real x,
                 streamsize decimals,
                 bool scientific_arg)
{
	ostringstream os;
	os << (scientific_arg ? scientific : fixed); os. precision (decimals); os << x;
	return os. str ();
}



Real str2real (const string& s)
{ 
  string s1 (s);
  trim (s1);
  if (s1. empty ())
  	throw runtime_error ("Cannot convert an empty string to number");
  strLower (s1);

  if (s1 == "inf")
	  return INF;
  if (s1 == "-inf")
	  return -INF;
	if (   s1 == "NaN"
	    || s1 == "-NaN"
	    || s1 == "?"
	   )
	  return NaN;

#if 0
  return str2<double> (s1);  // slow
#else
	// 0
	trimSuffix (s1, "e+00");
  while (! s1. empty () && s1 [0] == '0')
    s1. erase (0, 1);
  if (! contains (s1, 'e') && contains (s1, '.'))
	  while (! s1. empty () && s1. back () == '0')
	    s1. erase (s1. size () - 1);
  if (   s1. empty () 
  	  || s1 == "."
  	  || s1 == "-0."
  	 )
  	return 0.0;  	  

  const double res = atof (s1. c_str ());
  if (! res)
  	throw runtime_error ("Cannot convert " + strQuote (s) + " to number");
  	
  return res;
#endif
}



long round (Real a)
{
  if (isNan (a))
    return numeric_limits<long>::min ();
  return (long) floor (a + 0.5);
}




/////////////////////////////////// Sum ///////////////////////////////////

Sum::Sum (Real initMaxAbs)
// Platform dependent ??
: n              (324),
  step           (30),  

  stepInc        (0.9),
  logStep        (log (step)),
  bucket         (NULL),
  maxAbs         (NaN),
  logStep_MaxAbs (NaN)
{     
  ASSERT (n > 0);
  ASSERT (step > 2.0 * 10.0);
  ASSERT (stepInc <= 0.9);
  ASSERT (step * stepInc > 1.0);

  bucket = new Real [n]; 
  reset (initMaxAbs); 
}



Sum::~Sum ()
{
  delete [] bucket;
}



int Sum::absAddendum2index (Real a) const
{
  ASSERT (maxAbs >= 0.0);
  ASSERT (! isNan (logStep_MaxAbs));
  ASSERT (a > 0.0);
  
  
  const int i = (int) floor (logStep_MaxAbs - log (a) / logStep);
  ASSERT (i <= (int) n);
      

  return i;
}



void Sum::stabilize (size_t i) 
{
  ASSERT (i < n);
  
 
  // If i is flagged then bucket [i] == 0.0


  if (bucket [i] == 0.0 || 
      absAddendum2index (fabs (bucket [i])) == (int) i)
    return;
  // Should happen rarely
  
  const Real X = bucket [i];
  bucket [i] = 0.0;
  add (X);
// Flag i
}



void Sum::setMaxAbs (Real initMaxAbs)
{ 
  maxAbs = initMaxAbs;
  if (isNan (maxAbs))
    logStep_MaxAbs = NaN;
  else
    {
      ASSERT (maxAbs >= 0.0);
      logStep_MaxAbs = log (maxAbs) / logStep;
    }
}



void Sum::reset (Real initMaxAbs)
{ 
  FOR (size_t, i, n)
    bucket [i] = 0.0;
    
  hasPlusInf  = false;
  hasMinusInf = false;

  setMaxAbs (initMaxAbs);
}



void Sum::add (Real x)
{
  ASSERT (! isNan (x));
  if (x == INF)
  {
    hasPlusInf = true;
    return;
  }
  if (x == -INF)
  {
    hasMinusInf = true;
    return;
  }
 
 
  if (x == 0.0)
    return;
      
  
  const Real a = fabs (x);
  ASSERT (a != INF);
  
  
  if (isNan (maxAbs))
    {
      setMaxAbs (a * step * stepInc);
      ASSERT (absAddendum2index (a) == 0);
      ASSERT (bucket [0] == 0.0);
    }
  else if (a > maxAbs)
    {  
      setMaxAbs (a * step * stepInc);
      ASSERT (absAddendum2index (a) == 0);
      // Going up may violate the bucket compatibility conditions
      FOR_REV (size_t, i, n)
        stabilize (i);
    }

    
  const int i = absAddendum2index (a);
  ASSERT (i >= 0);     
  if (bucket [i] == 0.0)
    bucket [i] = x;
  else
    {
      // bucket compatibility conditions
      bucket [i] += x;
      stabilize ((size_t) i);
    }
}



Real Sum::get () const
{
  if (hasPlusInf)
    if (hasMinusInf)
      return NaN;
    else
      return INF;
  else
    if (hasMinusInf)
      return -INF;
  
  Real s = 0.0;
  FOR_REV (size_t, i, n)
    s += bucket [i];
  return s;
}




////////////////////////////////// SumLn /////////////////////////////////

void SumLn::reset ()
{
  maxLnX = -INF;
  sum. reset ();
}



void SumLn::addExp (Real lnX)
{
  if (lnX == -INF)
    return;  
 
  if (lnX - maxLnX > 3.0)  // PAR
  {
    const Real s = sum. get ();
    sum. reset ();
    sum. add (s / exp (lnX - maxLnX));
    maxLnX = lnX;
  }
    
  sum. add (exp (lnX - maxLnX));
}




//////////////////////////////////////////////////////////////////////////

Real Series::get (Real maxError) const
{
	ASSERT (maxError > 0);
	
	// integral(x+1) < \sum_{x'=x+1}^\infty f(x') < integral(x)
	// 0 < integral(x) - integral(x+1) < f(x)
	// Sum ??
	Real s = 0;
  FOR_START (uint, x, start, numeric_limits<uint>::max ())
	{
		const Real y = f (x);
		s += y;
		// error = integral (x) - integral (x + 1)
		if (y <= maxError)
			return s + (integral (x) + integral (x + 1)) / 2;
	}
  NEVER_CALL;
  return NaN;
}




/////////////////////////// ContinuedFraction /////////////////////////////

namespace 
{
  inline void setNonZero (Real &target,
		              		    Real source)
    { target = nullReal (source) ? epsilon : source; }
}

  
  
bool ContinuedFraction::lentz (size_t maxIter,
                               Real &result) const
{
  ASSERT (maxIter > 1);
  
  setNonZero (result, getB (0));
  Real c = result;
  Real d = 0.0;
  FOR (size_t, i, maxIter)
  {
    const Real a = getA (i + 1);
    const Real b = getB (i + 1);
    setNonZero (d, b + a * d);
    setNonZero (c, b + a / c);
    d = 1.0 / d;
    const Real delta = c * d;
    if (eqReal (delta, 1.0))
      return true;
    result *=  delta;
  }
        
  return false;
}




///////////////////////////// Special Functions ////////////////////////////

Real lnGamma (Real x)
{
  ASSERT (x > 0.0);

  static constexpr Real c [6] = { 76.18009172947146,
	                               -86.50532032941677,
	                                24.01409824083091,
	                                -1.231739572450155,
	                                 0.001208650973866179,
	                                -0.5395239384953e-5};
                         
  const Real a = x + 4.5;
  const Real b = (x - 0.5) * log (a) - a;
  
  Real s = 1.000000000190015;
  FOR (size_t, i, 6)
    s += c [i] / (x + (Real) i);    
    
  return  b + log (2.5066282746310005 * s);
}



namespace
{
  struct ContinuedFraction_LnComIncGamma : ContinuedFraction
  {
    Real x;
    Real a;
  	
    ContinuedFraction_LnComIncGamma (Real x_arg,
                                     Real a_arg)
      : ContinuedFraction ()
      , x (x_arg)
      , a (a_arg)
      { ASSERT (x > 0.0);
        ASSERT (a > 0.0);
      }

  private:
    Real getA (size_t index) const final
      { ASSERT (index > 0);
        if (index == 1)
          return 1.0;
        const Real j = (Real) (index - 1);
        return  - j * (j - a);
      }
    Real getB (size_t index) const final
      { return index ? x - a + 1.0 + 2.0 * (Real) (index - 1) : 0.0; }
  };
}



Real lnComIncGamma (Real x,
                    Real a)
{
  ASSERT (x >= 0.0);
  ASSERT (a > 0.0);
    
  if (! x)
    return lnGamma (a);

  const ContinuedFraction_LnComIncGamma cf (x, a);  
  Real g;
  if (! cf. lentz (1000, g))  // PAR
    return NaN;
  return - x + a * log (x) + log (g);  
}




/////////////////////////

Prob toProb (Real x, Real delta)
{ 
  ASSERT (delta >= 0.0);
  if (isNan (x))
    return x;
  if (isProb (x))
	  return x;
	if (nullReal (x, delta))
		return 0.0;
	if (eqReal (x, 1, delta))
		return 1.0;
	throw runtime_error ("Not a probability: " + real2str (x));
}
  


string prob2str (Prob x)
{
  ASSERT (isProb (x));  
  const streamsize decimals = 4;  // PAR
  if (x < 0.9)
  	return real2str (x, decimals);
  return "1-" + real2str (1 - x, decimals);
}



Real lnFactorial (uint n)
// William H. Press et al, Numerical Recipes
{
  constexpr size_t resMax = 200;  // Cf. initNumeric()
  static Real res [resMax] = {log(1.), log(1.), log(2.), log(6.), log(24.), log(120.), log(720.)};
                            //     0        1        2        3         4          5          6
  static size_t resN = 6;  // res[0..resN] are defined

  if (n >= resMax)
    return lnGamma (n + 1);
    
  while (resN < n)
  {
    resN++;
    res [resN] = res [resN - 1] + log ((Real) resN);
  }
  return res [n];
}



Real erf (Real x) 
// Numeric Recipes on C, p. 221
{
  const Real z = fabs (x);
  const Prob t = 1 / (1 + z/2);
  const Real ans = t * exp (- sqr (z) - 1.26551223 
                             + t * (1.00002368 
                                    + t * (0.37409196
                                           + t * (0.09678418
                                                  + t * (-0.18628806
                                                         + t * (0.27886807
                                                                + t * (-1.13520398
                                                                       + t * (1.48851587
                                                                              + t * (-0.82215223
                                                                                     + t * 0.17087277
                                                                                    )
                                                                             )
                                                                      )
                                                               )
                                                        )
                                                 )
                                          )
                                    )
                              );
  return 1 - (x >= 0 ? ans : (2 - ans)); 
}



namespace
{
	
struct ZetaSeries : Series
{
	Real s;
	ZetaSeries (Real s_arg,
	            uint start_arg)
	  : s (s_arg)
	  { ASSERT (s > 1);
	  	start = start_arg; 
	  	ASSERT (start >= 1);
	  }
  Real f (uint x) const 
    { return 1 / pow (x, s); }
  Real integral (uint x) const 
	  { return 1 / (pow (x, s - 1) * (s - 1)); }
};
	
}



Real zeta (Real alpha,
           uint from)
{
  ASSERT (greaterReal (alpha, 1));
  
  constexpr size_t size = 1000;  // PAR
	static Real res [size];
	{
	  static bool init = false;
		if (! init)
		{
			FOR (size_t, i, size)
			  res [i] = 0;
			init = true;
		}
	}
	size_t index = NO_INDEX;
  if (from == 1 && alpha > 1 && alpha < 2)
  {
	  index = (size_t) round ((alpha - 1) * size);
	  if (index < size)
	  	if (const Real res_ = res [index])
	  	  return res_;
  }
  
	ZetaSeries zs (alpha, from);
	const Real res_ = zs. get (1e-4);  // PAR
	
  if (from == 1 && index < size)
  	res [index] = res_;
	return res_;
}



Real zeta (Real alpha,
           uint from,
           uint to)
{
  Real s = 0;
  FOR_REV_END (uint, i, from, to + 1)
    s += pow (i, - alpha);
  return s;
}



namespace
{
	
struct ZetaLnSeries : Series
{
	Real s;
	ZetaLnSeries (Real s_arg,
	              uint start_arg)
	  : s (s_arg)
	  { ASSERT (s > 1);
	  	start = start_arg; 
	  	ASSERT (start >= 1);
	  }
  Real f (uint x) const 
    { return 1 / pow (x, s) * log (x); }
  Real integral (uint x) const 
	  { const Real b = 1 - s;
	    return pow (x, b) * (1 / sqr (b) - log (x) / b);
	  } 
};
	
}



Real zetaLn (Real alpha,
             uint from)
{
  ASSERT (greaterReal (alpha, 1));

  constexpr size_t size = 1000;  // PAR
	static Real res [size];
	{
	  static bool init = false;
		if (! init)
		{
			FOR (size_t, i, size)
			  res [i] = 0;
			init = true;
		}
	}
	size_t index = NO_INDEX;
  if (from == 1 && alpha > 1 && alpha < 2)
  {
	  index = (size_t) round ((alpha - 1) * size);
	  if (index < size)
	  	if (const Real res_ = res [index])
	  	  return res_;
  }
  
  
  ZetaLnSeries zs (alpha, from);
	const Real res_ = zs. get (1e-4);  // PAR
	
  if (from == 1 && index < size)
  	res [index] = res_;
	return res_;
}



Real zetaLn (Real alpha,
             uint from,
             uint to)
{
  Real s = 0;
  FOR_REV_END (uint, i, from, to + 1)
    s += pow (i, - alpha) * log (i);
  return s;
}




// WeightedMeanVar

void WeightedMeanVar::add (Real x,
                           Real weight)
{ 
  ASSERT (! isNan (weight));  
  if (weight == 0)
    return;

  ASSERT (! isNan (x));
  const Real x2 = sqr (x);  
  weightedSum  += x  * weight;
  weightedSum2 += x2 * weight;
  weights      +=      weight;
}



void WeightedMeanVar::add (const WeightedMeanVar& other)
{ 
  weightedSum  += other. weightedSum;
  weightedSum2 += other. weightedSum2;
  weights      += other. weights;
}



void WeightedMeanVar::subtract (const WeightedMeanVar& other)
{ 
  ASSERT (geReal (weights, other. weights));
  
  weightedSum  -= other. weightedSum;
  weightedSum2 -= other. weightedSum2;
  weights      -= other. weights;
}




// MeanVar

void MeanVar::qc () const
{
  if (! qc_on)
    return;
  QC_ASSERT (n >= 0);
  QC_ASSERT (s2 >= 0);
  QC_IMPLY (s2 > 0, n > 0);
}




// Histogram

Histogram::Histogram (Real start_arg,
                      Real stop_arg,
                      Real binRange_arg)
: start (start_arg)
, stop (stop_arg)
, binRange (binRange_arg)
, bins ((size_t) round (ceil ((stop - start) / binRange)), 0)
{ 
  ASSERT (lessReal (start, stop)); 
  ASSERT (binRange > 0.0);
  ASSERT (bins. size ());      
}



void Histogram::saveText (ostream &os) const
{ 
  os << "bin\tcount" << endl;
  FFOR (size_t, i, bins. size ())
    os << start + binRange * (Real) i << '\t' << bins [i] << endl;
}



size_t Histogram::getBin (Real x) const
{ 
  if (! betweenEqualReal (x, start, stop))
    return NO_INDEX;
  size_t bin = (size_t) round (floor ((x - start) / binRange));
  if (bin == bins. size ())
    bin--;
  ASSERT (bin < bins. size ());
  return bin;
}



Histogram& Histogram::operator<< (Real x)
{ 
  const size_t bin = getBin (x);
  if (bin != NO_INDEX)
    bins [bin] ++;
  return *this;
}




}




