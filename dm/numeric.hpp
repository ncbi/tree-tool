// numeric.hpp

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


#ifndef NUMERIC_HPP_69862
#define NUMERIC_HPP_69862


#include <limits>
#include <cmath>
#include "../common.hpp"
using namespace Common_sp;



namespace DM_sp
{


bool initNumeric ();
  // Module initialization
  // Invoked automatically
  // Invokes: initCommon()


template <typename T /*:integer*/>  
  inline T divide (T n, 
                   T divisor,
                   T &remainder)
  { static_assert (numeric_limits<T>::is_integer(), "T must be integer");
    if (n >= 0)
    { remainder = n % divisor;
      return n / divisor;
    }
    const T absN = abs (n);
    const T result = - ((absN - 1) / divisor + 1);
    remainder = n - result * divisor;
    return result;
  }
  
template <typename T /*: number*/>  
  inline int num2sign (T t)
    { return t > 0 ? 1 : (t < 0 ? -1 : 0); }

template <typename T /*: number*/>  
  inline T sqr (T x)
    { return x * x; }



typedef  double  Real;
  // The floating-point datatype (short and fast)
  // Platform-specific
static_assert (numeric_limits<Real>::is_specialized);
static_assert (! numeric_limits<Real>::is_integer);
static_assert (numeric_limits<Real>::is_signed);
static_assert (numeric_limits<Real>::has_infinity);
static_assert (numeric_limits<Real>::has_quiet_NaN);

string real2str (Real x,
                 streamsize decimals = 2,
                 bool scientific_arg = true);

Real str2real (const string& s);

constexpr Real inf = numeric_limits<Real>::infinity ();  
static_assert (inf / 2 == inf);
static_assert (inf > 0.0);

inline bool finite (Real x)
  { return -inf < x && x < inf; }

inline Real nvlReal (Real x, Real val)
  { return isNan (x) ? val : x; }

constexpr Real epsilon = sqrt (numeric_limits<Real>::epsilon ());  // = 1.49012e-08  
  // Precision

constexpr Real pi = acos (-1.0);
extern Real log_2;
extern Real log_10;
extern Real sqrt_2;


inline bool eqReal (Real x, Real y, Real delta = epsilon)
  { if (x == y)  // inf
      return true;
    const Real m = fabs (max (x, y));
    if (m == inf)
      return false;
    return fabs (x - y) <= max (m, 1.0) * delta; 
	}

inline bool sameReal (Real x, Real y, Real delta = epsilon)
  { return    (isNan (x) && isNan (y))
	         || eqReal (x, y, delta);
	}

inline bool nullReal (Real x, Real delta = epsilon)
  { return eqReal (x, 0.0, delta); } 

inline bool lessReal (Real x, Real y, Real delta = epsilon)
  { return x < y && ! eqReal (x, y, delta); }

inline bool greaterReal (Real x, Real y, Real delta = epsilon)
  { return x > y && ! eqReal (x, y, delta); }

inline bool leReal (Real x, Real y, Real delta = epsilon)
  { return x <= y || eqReal (x, y, delta); }

inline bool geReal (Real x, Real y, Real delta = epsilon)
  { return x >= y || eqReal (x, y, delta); }

inline bool positive (Real x, Real delta = epsilon)
  { return greaterReal (x, 0.0, delta); }

inline bool negative (Real x, Real delta = epsilon)
  { return lessReal (x, 0.0, delta); }
  
inline bool leRealRel (Real x, Real y, Real delta = epsilon)
  { return    leReal (x, y) 
           || leReal (x / y, 1.0, delta);
  }
  
inline bool geRealRel (Real x, Real y, Real delta = epsilon)
  { return    geReal (x, y) 
           || geReal (x / y, 1.0, delta);
  }
  
inline bool betweenReal (Real x, Real lo, Real hi, Real delta = epsilon) 
  { return geReal (x, lo, delta) && lessReal (x, hi, delta); }

inline bool betweenEqualReal (Real x, Real lo, Real hi, Real delta = epsilon) 
  { return geReal (x, lo, delta) && leReal (x, hi, delta); }

long round (Real x);
  // NaN -> numeric_limits::min()

inline bool isInteger (Real x)
  { return fabs (x - (Real) round (x)) < epsilon; }

inline bool maximizeReal (Real &a,
				                  Real b,
				                  Real delta)
  { if (a >= b - delta)
  	  return false;
 	  a = b;
 	  return true;
  }

inline bool minimizeReal (Real &a,
				                  Real b,
				                  Real delta)
  { if (a <= b + delta)
  	  return false;
 	  a = b;
 	  return true;
  }

inline void sqrEquation (Real a,
                         Real b,
                         Real c,
                         Real &x1,
                         Real &x2)
  { if (a)
    { const Real det = sqrt (sqr (b) - 4.0 * a * c);
      x1 = (- b - det) / (2.0 * a);
      x2 = (- b + det) / (2.0 * a);
    }
    else
    { x1 = x2 = - c / b; }
  }


inline Real ave_arith (Real x, Real y)
  { return 0.5 * (x + y); }
inline Real ave_geom (Real x, Real y)
  { return sqrt (x * y); }
inline Real ave_harm (Real x, Real y)   // harmonic
  { return 2.0 / (1.0 / x + 1.0 / y); }



struct LogReal : Root
{
  bool sign {true};
    // false <=> negative
  Real n {0.0};
  
  LogReal () = default;
  explicit LogReal (Real r)
    : sign (r >= 0)
    , n (log (fabs (r)))
    {}
    
  LogReal& operator*= (Real x)
    { if (x < 0.0)
      { toggle (sign);
        x = - x;
      }
      n += log (x);
      return *this;
    }
  Real get () const
    { return exp (n) * (sign ? 1.0 : -1.0); }
  bool isZero () const
    { return n == -inf; }
  bool nullReal () const
    { return DM_sp::nullReal (exp (n)); }
};



struct Sum final : Root, Nocopy
// Summation with better precision
{
private:
  size_t n {0};
    // # Buckets
  Real step {NaN};
    // >> 2.0
  Real stepInc {NaN};
  Real logStep {NaN};
    // = log(step)
  Real* bucket {nullptr};
    // Size = n
    // if bucket[i] != 0 then Threshold(i+1) < |bucket[i]| <= Threshold(i)
    //   where Threshold(i) = maxAbs / step^i
  Real maxAbs {NaN};
    // >= 0.0
  Real logStep_MaxAbs {NaN};
    // = log (maxAbs) / logStep
  bool hasPlusInf {false};
  bool hasMinusInf {false};
public:

  explicit Sum (Real initMaxAbs = NaN);
    // Invokes: reset()
 ~Sum ();

private:
  int absAddendum2index (Real a) const;
    // Return: < n
  void stabilize (size_t i);
    // Time: O(n)
  void setMaxAbs (Real intiMaxAbs);
public:
  void reset (Real initMaxAbs = NaN);
    // Faster if initMaxAbs >= max |X|
  void add (Real x);
    // More precise if |x|'s increase
    // Invokes: stabilize()
  Real get () const;
    // Return: sum bucket[]
};



struct SumLn final : Root, Nocopy
// Summation of natural logarithms with better precision
{
private:
  Real maxLnX {NaN};
  Sum sum;
public:

  SumLn ()
	   { reset (); }

  void reset ();
  void addExp (Real lnX);
    // More precise if |lnX|'s increase
  Real getLn () const
    { return log (sum. get ()) + maxLnX; }
};



struct Series : Root
{
	uint start {0};
	
	// Input: x >= start
  virtual Real f (uint x) const = 0;
    // Return: > 0
    // Requires: f(x+1) < f(x)
  virtual Real integral (uint x) const = 0;
    // Return: \int_x^\infty f(x') dx'
  Real get (Real maxError) const;
    // REturn: \sum_{x=start}^\infty f(x) with error <= maxError
};



struct ContinuedFraction : Root
{
  virtual Real getA (size_t index) const = 0;
  virtual Real getB (size_t index) const = 0;  
  bool lentz (size_t maxIter,
              Real &result) const;
    // Output: result
    // Return: true if converges
    // Recommended: maxIter = 100 
};



// Special Functions 

Real lnGamma (Real x);
  // Logarithm of Gamma function
  // Requires: x > 0.0
  // Author: Lanczos, C. 1964, SIAM Journal on Numerical Analysis, ser. B, vol. 1, pp. 86-96

Real lnComIncGamma (Real x,
                    Real a);
  // Return: Logarithm of complementary incomplete Gamma function, may be NaN if not converged
  // Requires: x >= 0.0, a > 0.0
  // Invokes: lnGamma (), ContinuedFraction::lentz ()
  // lnComIncGamma (x, 0.0) = lnGamma (X)



// Probability
typedef Real Prob;

inline bool isProb (Real p)
  { return p >= 0.0 && p <= 1.0; }
  
inline bool probIsLogic (Prob p)
  { return p == 0.0 || p == 1.0; }
  
Prob toProb (Real x, 
             Real delta = 1e-5);  // PAR
  
inline void makeProb (Real &x, 
                      Real delta = 1e-5)  // PAR
  { x = toProb (x, delta); }

inline Prob negateProb (Prob p,
                        bool act)
  { return act ? 1.0 - p : p; }
  // negateProb(negateProb(p,a),b) = negateProb(p,(a+b)%2)

inline Real prob2info (Prob p)
  { return - log2 (p); }

string prob2str (Prob x);

inline Prob ebool2prob (ebool b)
  { switch (b)
		{ 
			case etrue:  return 1.0;
			case efalse: return 0.0;
			case enull:  return 0.5;
		}
		throw runtime_error ("Never call");
  }

inline ebool prob2ebool (Prob p)
  { if (eqReal (p, 0.5))
  	  return enull;
  	if (p < 0.5)
  		return efalse;
  	return etrue;
  }

inline Real convexCombination (Prob p,
                               Real a,
                               Real b)
  { return p == 0.0 ? b : 
  	       p == 1.0 ? a :
  	       p * a + (1.0 - p) * b; 
  }

inline Prob probit (Real a)
  { const Real b = exp (a);
    return b / (b + 1.0);
  }

inline Real logit (Prob p)
  { return log (p / (1.0 - p)); }

Real lnFactorial (uint n);

inline Real multiplyLog (Real a, Real logB)
  { if (a == 0.0 /*&& logB == -inf*/)
  	  return 0.0;
  	return a * logB;
  }



Real erf (Real x);
  // Fractional error < 1.2e-7 for N(0,1)


// Riemann zeta function
// Input: from >= 1
Real zeta (Real alpha,
           uint from = 1);
  // Return: \sum_{i=from}^\infty i^{-alpha}
  // Input: alpha > 1
Real zeta (Real alpha,
           uint from,
           uint to);
  // Return: \sum_{i=from}^to i^{-alpha}
Real zetaLn (Real alpha,
             uint from = 1);
  // Return: \sum_{i=from}^\infty i^{-alpha} * ln i
  // Input: alpha > 1
Real zetaLn (Real alpha,
             uint from,
             uint to);
  // Return: \sum_{i=from}^to i^{-alpha} * ln i
 


struct WeightedMeanVar
// Weighted mean
{
  Real weightedSum {NaN};
  Real weightedSum2 {NaN};
  Real weights {NaN};

  
  WeightedMeanVar ()
    { clear (); }

    
  void clear ()
    { weightedSum = 0.0;
      weightedSum2 = 0.0;
      weights = 0.0;
    }
  void saveText (ostream &os) const
    { os         << weightedSum 
         << '\t' << weightedSum2 
         << '\t' << weights
         << endl; 
    }
  void add (Real x,
            Real weight = 1.0);
    // Input: weight: may be < 0
    // Invariant: { add (x, weight); add (x, - weight); }
  void add (const WeightedMeanVar& other);
  void subtract (const WeightedMeanVar& other);
  void addValue (Real x)
    { weightedSum2 += 2.0 * x * weightedSum + sqr (x) * weights;
      weightedSum  += x * weights;
    }
    // {for (i) add (x_i, w_i); addValue (x); } = {for (i) add (x_i + x, w_i); }
  Real getMean () const
    { return weightedSum ? weightedSum / weights : 0.0; }
  Real getVar () const
    { return (weightedSum2 ? weightedSum2 / weights : 0.0) - sqr (getMean ()); }
  Real getSD () const
    { return sqrt (getVar ()); }
  Real mean2var (Real mean) const
    { const Real bias = getMean () - mean;
    	return getVar () + sqr (bias); 
    }
  Real getOutlier_min (Real zScore = 3.0) const
    { return getMean () + zScore * getSD (); }
};



struct MeanVar final : Root
{
  streamsize decimals {2};
  int n {0};
  Real s {NaN};
  Real s2 {NaN};
  Real v_min {NaN};
  Real v_max {NaN};
  

  explicit MeanVar (streamsize decimals_arg = 2)
    : decimals (decimals_arg)
    { clear (); }
  void qc () const override;
  void saveText (ostream &os) const override
    { os << fixed; os. precision (decimals);
      os << "[" << v_min << " .. " << getMean () << " +- " << getSD () << " .. " << v_max << "] N = " << n;
    }
  void clear () override
    { n = 0;
    	s = 0.0;
    	s2 = 0.0;
    	v_min =  inf;
    	v_max = -inf;
    }


  MeanVar& operator<< (Real x)
    { n++;
    	s += x;
    	s2 += sqr (x);
    	minimize (v_min, x);
    	maximize (v_max, x);
    	return *this;
    }
  void add (const MeanVar& other,
            Real delta)
    { add_ (other, delta, false); }
  void subtract (const MeanVar& other,
                 Real delta)
    { add_ (other, delta, true); }
private:
  void add_ (const MeanVar& other,
             Real delta,
             bool reverse)
    { const int mode = reverse ? -1 : 1;
      n  += mode * other. n;
      s  += mode * (other. s + other. n * delta);
      s2 += mode * (other. s2 + 2 * other. s * delta + other. n * sqr (delta)); 
    }
public:

  // Estimates    
  Real getMean () const
    { return n ? s / n : NaN; }
  Real getVar (bool biased = true) const
    { return max (0.0, s2 - n * sqr (getMean ())) / (biased ? n : (n - 1)); }
    // biased = MLE
  Real getSD (bool biased = true) const
    { return sqrt (getVar (biased)); }
  Real mean2var (Real mean) const
    { const Real bias = getMean () - mean;
    	return getVar () + sqr (bias); 
    }
    // Unbiased
};



struct Correlation
{
  MeanVar a_mv;
  MeanVar b_mv;
  Real ab {NaN};

  
  Correlation ()
    { clear (); }


  void clear ()
    { a_mv. clear ();
      b_mv. clear ();
      ab = 0.0;
    }
  void add (Real a,
            Real b)
    { a_mv << a;
      b_mv << b;
      ab += a * b;
    }
  Real getCovariance () const
    { return a_mv. n ? ab / a_mv. n - a_mv. getMean () * b_mv. getMean () : 0.0; }
  Real getCorrelation () const
    { const Real a_sd = a_mv. getSD ();
      const Real b_sd = b_mv. getSD ();
      return (a_sd && b_sd) ? getCovariance () / (a_sd * b_sd) : 0.0; 
    }
};



struct Histogram final : Root
{
  Real start {NaN};
  Real stop {NaN};
  Real binRange {NaN};
  Vector<size_t> bins;
  
  Histogram (Real start_arg,
             Real stop_arg,
             Real binRange_arg);
  void saveText (ostream &os) const override;
  
  size_t getBin (Real x) const;
  Histogram& operator<< (Real x);
};

 

}  



#endif
