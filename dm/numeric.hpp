// numeric.hpp

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


template <class T /*:integer*/>  
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
  
template <class T /*: number*/>  
  inline int num2sign (T t)
    { return t > 0 ? 1 : (t < 0 ? -1 : 0); }

template <class T>  
  inline T sqr (T x)
    { return x * x; }

template <class T>
  inline T absDiffUnsigned (T x, T y)
    { if (x > y) return x - y; return y - x; }



///////////// Platform-specific /////////////////////
typedef double Real;
  // The floating-point datatype (short and fast)
  // See fstFloat.cpp  
/////////////////////////////////////////////////////

string real2str (Real x,
                 streamsize decimals = 2,
                 bool scientific_arg = true);

Real str2real (const string& s);

const Real INF = numeric_limits<Real>::infinity ();  
  // No "isinf()"

inline bool finite (Real x)
  { return -INF < x && x < INF; }

const Real NaN = numeric_limits<Real>::quiet_NaN ();  

inline bool isNan (Real x)
  { return x != x; }

inline Real nvlReal (Real x, Real val)
  { return isNan (x) ? val : x; }

extern const Real epsilon;
  // Precision
extern const Real pi;
extern Real log_2;
extern Real log_10;
extern Real sqrt_2;


inline bool eqReal (Real x, Real y, Real delta = epsilon)
  { if (x == y)  // INF
      return true;
    const Real m = fabs (max (x, y));
    if (m == INF)
      return false;
    return fabs (x - y) <= max (m, 1.0) * delta; 
	}

inline bool sameReal (Real x, Real y, Real delta = epsilon)
  { return    (isNan (x) && isNan (y))
	         || eqReal (x, y, delta);
	}

inline bool nullReal (Real x, Real delta = epsilon)
  { return eqReal (x, 0, delta); } 

inline bool lessReal (Real x, Real y, Real delta = epsilon)
  { return x < y && ! eqReal (x, y, delta); }

inline bool greaterReal (Real x, Real y, Real delta = epsilon)
  { return x > y && ! eqReal (x, y, delta); }

inline bool leReal (Real x, Real y, Real delta = epsilon)
  { return x <= y || eqReal (x, y, delta); }

inline bool geReal (Real x, Real y, Real delta = epsilon)
  { return x >= y || eqReal (x, y, delta); }

inline bool positive (Real x, Real delta = epsilon)
  { return greaterReal (x, 0, delta); }

inline bool negative (Real x, Real delta = epsilon)
  { return lessReal (x, 0, delta); }
  
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

inline Real nonNegative (Real x)
  { if (negative (x))
      throw logic_error ("negative x = " + toString (x));
    else
      return max (0.0, x);
  }

inline void sqrEquation (Real a,
                         Real b,
                         Real c,
                         Real &x1,
                         Real &x2)
  { if (nullReal (a))
    { x1 = x2 = - c / b; }
    else
    { const Real det = sqrt (sqr (b) - 4 * a * c);
      x1 = (- b - det) / (2 * a);
      x2 = (- b + det) / (2 * a);
    }
  }


inline Real ave_arith (Real x, Real y)
  { return 0.5 * (x + y); }
inline Real ave_geom (Real x, Real y)
  { return sqrt (x * y); }
inline Real ave_harm (Real x, Real y)   // harmonic
  { return 2 / (1 / x + 1 / y); }



struct LogReal : Root
{
  bool sign {true};
    // false <=> negative
  Real n {0};
  
  LogReal () = default;
  explicit LogReal (Real r)
    : sign (r >= 0)
    , n (log (fabs (r)))
    {}
    
  LogReal& operator*= (Real x)
    { if (x < 0)
      { toggle (sign);
        x = - x;
      }
      n += log (x);
      return *this;
    }
  Real get () const
    { return exp (n) * (sign ? 1 : -1); }
  bool isZero () const
    { return n == -INF; }
  bool nullReal () const
    { return DM_sp::nullReal (exp (n)); }
};



struct Sum : Root, Nocopy
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



struct SumLn : Root, Nocopy
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

bool lnComIncGamma (Real x,
                    Real start,
                    size_t  maxIter,
                    Real &result);
  // Logarithm of complementary incomplete Gamma function
  // Output: result
  // Return: true if converges
  // Requires: X > 0.0, Start >= 0.0
  // Invokes: lnGamma (), CONTINUED_FRACTION::Lentz ()
  // Result of LnComIncGamma (X, 0.0) = LnGamma (X)



// Probability
typedef Real Prob;

inline bool isProb (Real p)
  { return p >= 0 && p <= 1; }
  
inline bool probIsLogic (Prob p)
  { return p == 0 || p == 1; }
  
Prob toProb (Real x);
  
inline void makeProb (Real &x)
  { x = toProb (x); }

inline Prob negateProb (Prob p,
                        bool act)
  { return act ? 1 - p : p; }
  // negateProb(negateProb(p,a),b) = negateProb(p,(a+b)%2)

inline Real prob2info (Prob p)
  { return - log2 (p); }

string prob2str (Prob x);

inline Prob ebool2prob (ebool b)
  { switch (b)
		{ 
			case ETRUE:  return 1;
			case EFALSE: return 0;
			case UBOOL:  return 0.5;
		}
		throw runtime_error ("Never call");
  }

inline ebool prob2ebool (Prob p)
  { if (eqReal (p, 0.5))
  	  return UBOOL;
  	if (p < 0.5)
  		return EFALSE;
  	return ETRUE;
  }

inline Real convexCombination (Prob p,
                               Real a,
                               Real b)
  { return p == 0 ? b : 
  	       p == 1 ? a :
  	       p * a + (1 - p) * b; 
  }

inline Prob probit (Real a)
  { const Real b = exp (a);
    return b / (b + 1);
  }

inline Real logit (Prob p)
  { return log (p / (1 - p)); }

Real lnFactorial (uint n);

inline Real multiplyLog (Real a, Real logB)
  { if (a == 0 /*&& logB == -INF*/)
  	  return 0;
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
    { weightedSum = 0;
      weightedSum2 = 0;
      weights = 0;
    }
  void add (Real x,
            Real weight = 1);
    // Input: weight: may be < 0
    // Invariant: { add (x, weight); add (x, - weight); }
  void add (const WeightedMeanVar& other);
  void subtract (const WeightedMeanVar& other);
  void addValue (Real x)
    { weightedSum2 += 2 * x * weightedSum + sqr (x) * weights;
      weightedSum  += x * weights;
    }
    // {for (i) add (x_i, w_i); addValue (x); } = {for (i) add (x_i + x, w_i); }
  Real getMean () const
    { return weightedSum ? weightedSum / weights : 0; }
  Real getVar () const
    { return (weightedSum2 ? weightedSum2 / weights : 0) - sqr (getMean ()); }
  Real getSD () const
    { return sqrt (getVar ()); }
  Real mean2var (Real mean) const
    { const Real bias = getMean () - mean;
    	return getVar () + sqr (bias); 
    }
  Real getOutlier_min (Real zScore = 3) const
    { return getMean () + zScore * getSD (); }
};



struct MeanVar : Root
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
    	s = 0;
    	s2 = 0;
    	v_min =  INF;
    	v_max = -INF;
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
    { return a_mv. n ? ab / a_mv. n - a_mv. getMean () * b_mv. getMean () : 0; }
  Real getCorrelation () const
    { const Real a_sd = a_mv. getSD ();
      const Real b_sd = b_mv. getSD ();
      return nullReal (a_sd) || nullReal (a_sd) ? 0 : getCovariance () / (a_sd * b_sd); 
    }
};



struct Histogram : Root
{
  Real start {NaN};
  Real stop {NaN};
  Real binRange {NaN};
  Vector<size_t> bins;
  
  Histogram (Real start_arg,
             Real stop_arg,
             Real binRange_arg);
  void saveText (ostream &os) const;
  
  size_t getBin (Real x) const;
  Histogram& operator<< (Real x);
};

 

}  



#endif
