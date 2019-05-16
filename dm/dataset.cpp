// dataset.cpp

#undef NDEBUG
#include "../common.inc"

#include "dataset.hpp"

#include "optim.hpp"



namespace DM_sp
{ 


const string dmExt ("dm");
const string dmSuff ("." + dmExt);




// Obj

void Obj::qc () const
{ 
  if (! qc_on)
    return;
  Root::qc ();
		
  QC_IMPLY (! name. empty (), goodName (name));
  QC_ASSERT (mult >= 0.0); 
  QC_ASSERT (mult < INF);
}




// RealScale

string RealScale::getAverageStrValue (StringVector &&valuesStr)
{
  if (valuesStr. empty ())
    return "nan";
    
  if (valuesStr. size () == 1)
    return move (valuesStr [0]);
   
  MeanVar mv;
  for (string s : valuesStr)
  {
    replaceStr (s, ",", "");
    mv << str2real (s);
  }
   
  return toString (mv. getMean ());  
}






// Attr

const char* missingStr = "?";



Attr::Attr (const string &name_arg,
		        Dataset &ds_arg,
		        bool rightAlign_arg)
: ds (ds_arg)
, dsIt (ds_arg. attrs. end ())
, name (name_arg)
, rightAlign (rightAlign_arg)
{ 
  ds_arg. addAttr (this);  
}



Attr::~Attr ()
{
  ASSERT (dsIt != ds. attrs. end ());
  var_cast (ds). attrs. erase (dsIt);
  var_cast (ds). name2attr_. erase (name);
}



void Attr::qc () const
{
  if (! qc_on)
    return;
  Root::qc ();
  
  QC_ASSERT (*dsIt == this);
}



void Attr::moveAfter (const Attr* pred)
{
  IMPLY (pred, & pred->ds == & ds);
  ASSERT (pred != this);
  
  if (pred && std::next (pred->dsIt) == dsIt)
    return;
  if (! pred && ds. attrs. begin () == dsIt)
    return;
  
  Dataset& ds_ = var_cast (ds);

  ds_. attrs. erase (dsIt);
  
  ASSERT (! ds. attrs. empty ());
  auto it = pred ? pred->dsIt : ds_. attrs. begin ();
  it++;  // Visual C++: "expression: list iterator not incrementable" <= sentinel of an empty list
  dsIt = ds_. attrs. insert (it, this);
}



size_t Attr::countMissings () const
{
  size_t n = 0;
  FFOR (size_t, objNum, ds. objs. size ())
    if (isMissing (objNum))
      n++;
  return n;
}



bool Attr::existsMissing (size_t &objNum) const
{
  for (objNum = 0; objNum < ds. objs. size (); objNum++)
    if (ds. objs [objNum] -> mult > 0 && isMissing (objNum))
      return true;
      
  objNum = NO_INDEX;
  return false;
}



bool Attr::missingsAll () const
{
  FFOR (size_t, objNum, ds. objs. size ())
    if (! isMissing (objNum))
      return false;
  return true;
}



#if 0
void Attr::summary (ostream &f) const
{
  f << name << endl;
  f << getTypeStr () << endl;
  if (isConstant ())
    f << "Constant" << endl;
}
#endif




// Attr1 

void Attr1::setMissingAll ()
{
  FFOR (size_t, objNum, ds. objs. size ())
    setMissing (objNum);
}



size_t Attr1::getWidth_max () const
{
  size_t width = string (missingStr). size ();  
  FFOR (size_t, i, ds. objs. size ())
    if (! isMissing (i))
      maximize (width, value2str (i). size ());
   
  return width;
}




// NumAttr1

namespace
{

struct ObjNum_Real
{
  size_t objNum;
  Real value;
  
  bool operator< (const ObjNum_Real &other) const
    { return value < other. value; }
};

}



NumAttr1::NumAttr1 (const string &name_arg,
                    Dataset &ds_arg,
                    const NumAttr1 &from)
: Attr1 (name_arg, ds_arg, true)
, RealScale (from. decimals)
{ 
  ASSERT (ds_arg. objs. size () == from. ds. objs. size ());
}



NumAttr1& NumAttr1::operator= (const NumAttr1& other)
{ 
  ASSERT (& other. ds == & ds);
  ASSERT (other. decimals == decimals);
  if (! RealScale::operator== (other))
    throw runtime_error (FUNC "Different RealScle in " + name + " and in " + other. name);
  return *this;
}



NumAttr1::Value& NumAttr1::operator[] (size_t objNum) 
{ 
  if (const RealAttr1* a = asRealAttr1 ())
    return var_cast (a) -> operator[] (objNum);
  throw runtime_error (FUNC "Assignment to NumAttr1"); 
}



void NumAttr1::getAverageScatter (const Sample &sample,
                                  Real &average,
	                                Real &scatter) const
{
	Real n = 0.0;
	Real s = 0.0;
	Real s2 = 0.0;
  for (Iterator it (sample); it ();)  
    if (! isMissing (*it))
    {
    	const Real x = getReal (*it);
    	n += it. mult;
    	s += it. mult * x;
    	s2 += it. mult * sqr (x);
    }
  
  average = s / n;
  scatter = s2 / n - sqr (average);
  IMPLY (! isNan (scatter), scatter >= 0.0);
}



Real NumAttr1::getSqr_ave (const Sample &sample) const
{
	Real n = 0.0;
	Real s = 0.0;
  for (Iterator it (sample); it ();)  
    if (! isMissing (*it))
    {
    	const Real x = getReal (*it);
    	n += it. mult;
    	s += it. mult * sqr (x);
    }
  
  return s / n;
}



Vector<NumAttr1*> NumAttr1::toNumAttr1 (Dataset &ds_arg) const 
{
  Vector<NumAttr1*> vec;
  auto* a = new RealAttr1 (name + "_num", ds_arg);  
  vec << a;
  FOR (size_t, i, ds. objs. size ())
    (*a) [i] = getReal (i);
  return vec;
}



Vector<RealAttr1*> NumAttr1::standardize (Dataset &ds_arg,
	                                        const Sample &sample) const
{
	ASSERT (& ds_arg == & ds);
	
  Vector<RealAttr1*> vec;
  Real average, scatter;
  getAverageScatter (sample, average, scatter);
  if (! isNan (scatter) && ! nullReal (scatter))
  {
    auto* a = new RealAttr1 (name + "_std", ds_arg);  
    vec << a;
    ASSERT (positive (scatter));
    const Real sd = sqrt (scatter);
    for (Iterator it (sample); it ();)  
      if (! isMissing (*it))
        (*a) [*it] = (getReal (*it) - average) / sd;
  }
  return vec;
}



bool NumAttr1::existsLessThan (const Sample &sample,
                               Value minValue,
                               size_t &objNum) const
{
  for (Iterator it (sample); it ();)
    if (getReal (*it) < minValue)
    {
      objNum = *it;
      return true;
    }
  return false;
}



void NumAttr1::getMinMax (const Sample &sample,
                          Value &min, 
                          Value &max) const
{
  min =  INF;
  max = -INF;
  for (Iterator it (sample); it ();)
  { 
    minimize (min, getReal (*it));
    maximize (max, getReal (*it));
  }
}



NumAttr1::Value NumAttr1::getMedian (const Sample &sample) const
{
  Vector<Value> vec;  vec. reserve (ds. objs. size ());
  for (Iterator it (sample); it ();)
  {
  	const Real x = getReal (*it);
  	if (! isNan (x))
    	vec << x;
  }
  if (vec. empty ())
  	return NaN;
  // Time: --> linear ??
  vec. sort ();
  return vec [vec. size () / 2];
}



Real NumAttr1::distr2outlier (const Sample &sample,
                              LocScaleDistribution &distr,
                              bool rightTail,
                              Real outlier_EValue_max) const
{
  Vector<ObjNum_Real> vec;  vec. reserve (ds. objs. size ());
  for (Iterator it (sample); it ();)  
    if (! isMissing (*it))
      vec << ObjNum_Real {*it, getSign (rightTail) * getReal (*it)};
  if (vec. size () <= 2)
    return NaN;
    
  vec. sort ();
    
  Real mult_sum = 0.0;
  Real s = 0.0;
  Real s2 = 0.0;
  Real mean = NaN;
  Real var = 0.0;
  for (const ObjNum_Real& it : vec)
  {
    const Real x = it. value;  // Next object

    if (positive (var) && mult_sum >= 0.5 * sample. mult_sum)  // PAR
    {
      distr. setMeanVar (mean, var);
      const Prob p = 1.0 - distr. cdf (x);
      if (leReal (p * mult_sum, outlier_EValue_max))
        return getSign (rightTail) * x;
    }

    const Real mult = ds. objs [it. objNum] -> mult;
    mult_sum += mult;
    ASSERT (geReal (sample. mult_sum, mult_sum));
    s  += mult * x;
    s2 += mult * sqr (x);
    mean = s / mult_sum;
    var = s2 / mult_sum - sqr (mean);
    ASSERT (! negative (var));
  }
  if (verbose ())
    cout << "mean = " << mean << ' ' << "SD = " << sqrt (var) << endl;
  
  return getSign (rightTail) * INF;  
}




// RealScale

const RealScale::Value RealScale::missing = NaN;

	


// RealAttr1 

RealAttr1::RealAttr1 (const string &name_arg,
						          Dataset &ds_arg,
						          streamsize decimals_arg)
: NumAttr1 (name_arg, ds_arg, decimals_arg)
, values (ds. objs. size (), missing)
{ 
  values. reserve (ds. objs. capacity ());
}



RealAttr1::RealAttr1 (const string &name_arg,
                      Dataset &ds_arg,
                      const RealAttr1 &from)
: NumAttr1 (name_arg, ds_arg, from)
, values (from. values)
{ 
  ASSERT (ds_arg. objs. size () == from. ds. objs. size ());
}



RealAttr1& RealAttr1::operator= (const RealAttr1& other)
{ 
  NumAttr1::operator= (other);
  values = other. values;
  return *this;
}



void RealAttr1::qc () const
{
  if (! qc_on)
    return;
	NumAttr1::qc ();
	
  // values[]
  QC_ASSERT (values. size () == ds. objs. size ());
}



bool RealAttr1::isConstant () const
{
  Real x_min = INF;
  Real x_max = -INF;
  FFOR (size_t, i, ds. objs. size ())
    if (! isMissing (i))
    {
    	minimize (x_min, values [i]);
    	maximize (x_max, values [i]);
    }
  return x_max - x_min < pow (10.0, - (int) decimals) / 2.0;
}



void RealAttr1::appendObj ()
{
  values << missing;
}



size_t RealAttr1::getInfCount () const
{
  size_t n = 0;
  FFOR (size_t, objNum, ds. objs. size ())
    if (! finite ((*this) [objNum]))
      n++;
  return n;
}



size_t RealAttr1::inf2missing ()
{
	size_t n = 0;
  FFOR (size_t, row, ds. objs. size ())
    if (! finite ((*this) [row]))
    {
      if (verbose ())
        cout << "Infinity:" 
             << ' ' << ds. objs [row] -> name 
             << ' ' << (*this) [row]
             << endl;
      setMissing (row);
      n++;
    }
  return n;
}



void RealAttr1::multiply (Real coeff)
{
  FFOR (size_t, objNum, ds. objs. size ())
    if (! isMissing (objNum))
      (*this) [objNum] = (*this) [objNum] * coeff;
  decimals += (streamsize) max<long> (0, DM_sp::round (- log10 (coeff)));
}



Real RealAttr1::normal_likelihood2max (const Sample &sample) const
{
  Vector<ObjNum_Real> vec;  vec. reserve (ds. objs. size ());
  for (Iterator it (sample); it ();)  
    if (! isMissing (*it))
      vec << ObjNum_Real {*it, getReal (*it)};
  if (vec. size () <= 2)
    return NaN;
    
  vec. sort ();
    
  size_t objNum_threshold = NO_INDEX;
  Real negLogLikelihood_min = INF;
  const Real x_max = vec. back (). value;
  Real mult_sum = 0;
  Real s = 0;
  Real s2 = 0;
  Normal normal;
  const Real coeff = log (2 * pi) + 1;
  for (const auto& it : vec)
  {
    const Real mult = ds. objs [it. objNum] -> mult;
    mult_sum += mult;
    const Real x = it. value;
    s  += mult * x;
    s2 += mult * sqr (x);
    const Real mean = s / mult_sum;
    const Real var = s2 / mult_sum - sqr (mean);
    ASSERT (! negative (var));
    if (! positive (var))
      continue;
    normal. setParam (mean, sqrt (var));
    ASSERT (geReal (sample. mult_sum, mult_sum));
    const Real rest = sample. mult_sum - mult_sum;
    const Real negLogLikelihood = 
        mult_sum / 2 * (coeff + log (var)) 
      + (nullReal (rest)
           ? 0
           : eqReal (x, x_max) 
             ? INF
             : rest * (log ((x_max - x) / (1 - normal. cdf (x))))
        );
    if (minimize (negLogLikelihood_min, negLogLikelihood))
      objNum_threshold = it. objNum;
  }
  
  return objNum_threshold == NO_INDEX ? NaN : getReal (objNum_threshold);
}



#if 0
Real RealAttr1::getSkew (Real mult,
                         Real mean,
                         Real variance) const
// Biased -??
{
  Real sumDiff3 = 0.0;
  for (Iterator it (this); it ();)
    sumDiff3 += it. mult * pow (values [*it] - mean, 3);

  const Real sd3 = variance * sqrt (variance);
  return sumDiff3 / (sqrt (15.0 * mult) * sd3);
}



Real RealAttr1::getKurt (Real mult,
                         Real mean,
                         Real variance) const
// Biased -??
{
  Real sumDiff4 = 0;
  for (Iterator it (this); it ();)
    sumDiff4 += it. mult * pow (values [*it] - mean, 4);

  const Real a = sumDiff4 / (mult * sqr (variance)) - 3.0;
  return a / (sqrt (96.0 / mult));
}



namespace {
Real histogramLikelihood (const Vector <Real> &vec,
                          Real delta)
// Return: Average negative log-likelihood
// Requires: vec is sorted ascending
{
  ASSERT (vec. size () >= 2); 
  ASSERT (positive (delta));
  

  Real negLogLikelihood = 0;
  size_t left  = 0;  // = min left  s.t. vec.at(left)  >= vec.at(center) - delta
  size_t right = 0;  // = min right s.t. vec.at(right) >= vec.at(center) + delta
  for (size_t center = 0; center < vec. size (); center++)
  {
    const Real centerValue = vec [center];

    while (left < center)
      if (vec [left] >= centerValue - delta)
        break;
      else
        left++;

    while (right < vec. size ())
      if (vec [right] >= centerValue + delta)
        break;
      else
        right++;
    
    ASSERT (left <= center);
    ASSERT (center < right);    
    // center is excluded
    const Real density =   (Real) (right - left - 1) 
                          / (delta * 2 * (Real) (vec. size () - 1));
    negLogLikelihood += - log (density);
  }

  
  return negLogLikelihood / (Real) vec. size ();
}
}



void RealAttr1::histogram (ostream &f) const
{
  Vector<Real> vec (ds->objs. size ());
  size_t i = 0;
  for (Iterator it (this); it (); )
  {
    vec [i] = values [*it];
    i++;
  }
  vec. resize (i);
  if (i <= 2)
    return;

  vec. sort ();

  const Real range = vec [vec. size () - 1] - vec [0];
  Real delta = range * 1.1;
  // ??
  while (delta >= range / (Real) vec. size ())
  {
    f << delta << " " << histogramLikelihood (vec, delta) << endl;
    delta /= 2;
  }
}
#endif




// Type transformation

/*
ORD_ATTR1* RealAttr1::MakeOrdAttr (const char* Suffix,
                                   Real       Min,
                                   Real       Max,
                                   size_t        CategNum,
                                   size_t        &IntervalObjNum) const
{
  ASSERT (! StrIsEmpty (Suffix));
  ASSERT (LessEqualFloat (Min, Max));
  ASSERT (CategNum > 1);
  
 
  // Bounds
  Matrix Bounds (false, CategNum + 1, 1);
  const Real Step = (Max - Min) / CategNum;
  ASSERT (! Negative (Step));
  For (i, CategNum + 1)
    Bounds. put (false, i, 0, Min + i * Step);
  ASSERT (EqualFloat (Bounds. Get (false, CategNum, 0), Max));
  

  // Attr
  char* S = StrNew (Name);
  StrAppend (S, Suffix);
  auto* Attr = new ORD_ATTR1 (S, ds);
  ASSERT (GoodObject (Attr));
  free (S);
  
  For (i, CategNum)
    {
      char CategName [256];
      sprintf (CategName, "%0.*f-%0.*f", (int) Decimals, Bounds. Get (false, i,     0),
                                         (int) Decimals, Bounds. Get (false, i + 1, 0));
      if (! Attr->AddCateg (CategName))
        {
          delete Attr;
          return NULL;
        }
    }
  
  IntervalObjNum = 0;
  FOR (objNum, ds->objs)
    if (! isMissing (objNum))
      For (i, CategNum)
        if (BetweenFloat (values [objNum], Bounds. Get (false, i,     0),
                                         Bounds. Get (false, i + 1, 0)))
          {
            Attr->values [objNum] = i;
            IntervalObjNum++;
            break;
          }
  
  
  return Attr;
}
*/



// Time series

#if 0
RealAttr1* RealAttr1::smoothUniform (const string &suffix,
                                   int start,
                                   int end) const
{
  ASSERT (! suffix. empty ());
  ASSERT (start <= end);
  

  auto* smoothed = new RealAttr1 (name + suffix, * var_cast (ds), decimals);
  
  
  for (Iterator it (ds); it ();)
  {
    Real s = 0;
    Real w = 0;
    for (int i = (int) *it + start; 
             i < (int) *it + end; 
             i++)
      if (between (i, 0, (int) ds->objs. size ()))
      {
        const Obj* obj = ds->objs [(size_t) i];
        if (   obj->active () 
            && ! isMissing ((size_t) i)
           )
        {
          s += values [(size_t) i];
          w += obj->mult;
        }
      }
         
    if (w)
      smoothed->values [*it] = s / w;
  }
         
  
  return smoothed;
}
#endif




// PositiveAttr1

void PositiveAttr1::qc () const
{
  if (! qc_on)
    return;
  RealAttr1::qc ();

  // values[]
  FFOR (size_t, i, ds. objs. size ())
    if (! isMissing (i))
    {
      QC_ASSERT (values [i] >= 0);
    }
}



Vector<RealAttr1*> PositiveAttr1::standardize (Dataset &ds_arg,
	                                             const Sample &sample) const 
{
	ASSERT (& ds_arg == & ds);
	
	RealAttr1* attr = logarithmize (ds_arg, "_log_std");
	
  Vector<RealAttr1*> vec;
  Real average, scatter;
  attr->getAverageScatter (sample, average, scatter);
  if (isNan (scatter) || nullReal (scatter))
  	delete attr;
  else
  {
    vec << attr;
    ASSERT (positive (scatter));
    const Real sd = sqrt (scatter);
    for (Iterator it (sample); it ();)  
      if (! isMissing (*it))
        (*attr) [*it] = ((*attr) [*it] - average) / sd;
  }
  return vec;
}



PositiveAttr1* PositiveAttr1::standardizePositive (Dataset &ds_arg,
	                                                 const Sample &sample) const 
{
	ASSERT (& ds_arg == & ds);
	
  Real average, scatter;
  getAverageScatter (sample, average, scatter);
  if (isNan (average) || nullReal (average))
  	return nullptr;
  else
  {
		auto* attr = new PositiveAttr1 (name + "_pos_std", ds_arg, decimals);
    ASSERT (positive (average));
    for (Iterator it (sample); it ();)  
      if (! isMissing (*it))
        (*attr) [*it] = (*this) [*it] / average;
	  return attr;
  }
}



RealAttr1* PositiveAttr1::logarithmize (Dataset &ds_arg,
                                        const string &suffix) const
{
	ASSERT (& ds_arg == & ds);
		
	auto* attr = new RealAttr1 (name + suffix, ds_arg, decimals + 1);  // PAR
  FFOR (size_t, objNum, ds. objs. size ())
  {
  	const Real x = (*this) [objNum];
    if (x > 0)
      (*attr) [objNum] = log (x);
  }
  return attr;
}




// ProbAttr1

void ProbAttr1::qc () const
{
  if (! qc_on)
    return;
  RealAttr1::qc ();

  // values[]
  FFOR (size_t, i, ds. objs. size ())
    QC_IMPLY (! isMissing (i), isProb (values [i]));
}



Prob ProbAttr1::getProb (size_t objNum) const
{ 
  const Prob r = (*this) [objNum];
  IMPLY (! isNan (r), isProb (r));
  return r;
}




// IntAttr1

const IntAttr1::Value IntAttr1::missing = numeric_limits<IntAttr1::Value>::min ();



IntAttr1::IntAttr1 (const string &name_arg,
                    Dataset &ds_arg)
: NumAttr1 (name_arg, ds_arg, 0)
, values (ds. objs. size (), missing)
{ 
  values. reserve (ds. objs. capacity ());
}


                    
IntAttr1::IntAttr1 (const string &name_arg,
                    Dataset &ds_arg,
                    const IntAttr1 &from)
: NumAttr1 (name_arg, ds_arg, from)
, values (from. values)
{ 
  ASSERT (ds_arg. objs. size () == from. ds. objs. size ());
}


                    
IntAttr1& IntAttr1::operator= (const IntAttr1& other)
{ 
  NumAttr1::operator= (other);
  values = other. values;
  return *this;
}


  
void IntAttr1::qc () const
{
  if (! qc_on)
    return;
	NumAttr1::qc ();
	
  // values[]
  QC_ASSERT (values. size () == ds. objs. size ());
  FFOR (size_t, i, ds. objs. size ())
    if (! isMissing (i))
    {
      QC_ASSERT (finite (values [i]));
    }
}



bool IntAttr1::isConstant () const
{
  Value x_min = numeric_limits<Value>::max ();
  Value x_max = numeric_limits<Value>::min ();
  FFOR (size_t, i, ds. objs. size ())
    if (! isMissing (i))
    {
    	minimize (x_min, values [i]);
    	maximize (x_max, values [i]);
    }
  return x_max == x_min;
}





// BoolAttr1

const BoolAttr1::Value BoolAttr1::missing = UBOOL;



bool BoolAttr1::str2bool (const string &s) const
{
  if (s. size () == 1)
  {
    const char c = toUpper (s [0]);
    if (   c == '1'
        || c == 'T'
        || c == 'Y'
       ) 
      return true;
    if (   c == '0'
        || c == 'F'
        || c == 'N'
       ) 
      return false;
  }
  
  throw runtime_error (FUNC "Unknown Boolean value " + strQuote (s));
}



BoolAttr1::Value& BoolAttr1::operator[] (size_t objNum) 
{ 
  if (const ExtBoolAttr1* a = asExtBoolAttr1 ())
    return var_cast (a) -> operator[] (objNum);
  throw runtime_error (FUNC "Assignment to BoolAttr1"); 
}



void BoolAttr1::getStat (const Sample &sample,
                         array<size_t,3/*ebool*/> &stat) const
{
	stat. fill (0);
  for (Iterator it (sample); it ();)  
  	switch (getBool (*it))
  	{
  	  case EFALSE: stat [0] ++; break;
  	  case ETRUE:  stat [1] ++; break;
  	  case UBOOL:  stat [2] ++; break;
  	  default: ERROR;
  	}
}



Prob BoolAttr1::getProb (const Sample &sample) const
{
	Real n = 0;
	Real s = 0;
  for (Iterator it (sample); it ();)  
  {
  	const ebool x = getBool (*it);
  	n += it. mult;
  	if (x == ETRUE)
  	  s += it. mult;
  }
  return s / n;
}




// ExtBoolAttr1

ExtBoolAttr1::ExtBoolAttr1 (const string &name_arg,
                            Dataset &ds_arg) 
: BoolAttr1 (name_arg, ds_arg)
, values (ds. objs. size (), UBOOL)
{ 
  values. reserve (ds. objs. capacity ());
}



ExtBoolAttr1::ExtBoolAttr1 (const string &name_arg,
                            Dataset &ds_arg,
                            const ExtBoolAttr1 &from)
: BoolAttr1 (name_arg, ds_arg)
, values (from. values)
{
  ASSERT (ds_arg. objs. size () == from. ds. objs. size ());
}



void ExtBoolAttr1::qc () const
{
  if (! qc_on)
    return;
	BoolAttr1::qc ();
	  
  QC_ASSERT (values. size () == ds. objs. size ());
}



bool ExtBoolAttr1::isConstant () const
{ 
  ebool x_min = ETRUE;
  ebool x_max = EFALSE;
  FFOR (size_t, i, ds. objs. size ())
    if (! isMissing (i))
    {
    	minimize (x_min, values [i]);
    	maximize (x_max, values [i]);
    }
  return x_min == x_max;
}



ExtBoolAttr1& ExtBoolAttr1::operator= (const ExtBoolAttr1& other)
{ 
  ASSERT (& other. ds == & ds);
  values = other. values;
  return *this;
}




// CompactBoolAttr1

CompactBoolAttr1::CompactBoolAttr1 (const string &name_arg,
                                    Dataset &ds_arg,
                                    const CompactBoolAttr1 &from)
: BoolAttr1 (name_arg, ds_arg)
, values (from. values)
{
  ASSERT (ds_arg. objs. size () == from. ds. objs. size ());
}



void CompactBoolAttr1::qc () const
{
  if (! qc_on)
    return;
	BoolAttr1::qc ();
	  
  QC_ASSERT (values. size () == ds. objs. size ());
}



bool CompactBoolAttr1::isConstant () const
{ 
  bool x_min = true;
  bool x_max = false;
  for (const bool b : values)
  {
  	minimize (x_min, b);
  	maximize (x_max, b);
  }
  return x_min == x_max;
}



CompactBoolAttr1& CompactBoolAttr1::operator= (const CompactBoolAttr1& other)
{ 
  ASSERT (& other. ds == & ds);
  ASSERT (values. size () == other. values. size ());
  values = other. values;
  return *this;
}



void CompactBoolAttr1::setAll (bool value)
{
//values. resize (0, value);  
  values. resize (ds. objs. size (), value); 
}




// NominAttr1::Dependene

void NominAttr1::Dependence::qc () const
{
  if (! qc_on)
    return;
  QC_ASSERT (isProb (pValue));
}
      
      


// NominAttr1

const size_t NominAttr1::missing = NO_INDEX;
	


NominAttr1::NominAttr1 (const string &name_arg,
                        Dataset &ds_arg) 
: Attr1 (name_arg, ds_arg, false)
, values (ds. objs. size ())
{ 
  values. reserve (ds. objs. capacity ());
	setMissingAll (); 
}



NominAttr1::NominAttr1 (const string &name_arg,
                        Dataset &ds_arg,
                        const NominAttr1 &from)
: Attr1 (name_arg, ds_arg, false)
, categories (from. categories)
, values (from. values)
{
  ASSERT (ds_arg. objs. size () == from. ds. objs. size ());
}



NominAttr1& NominAttr1::operator= (const NominAttr1& other)
{ 
  ASSERT (& other. ds == & ds);
  ASSERT (categories == other. categories);
  values = other. values;
  return *this;
}



void NominAttr1::qc () const
{
  if (! qc_on)
    return;
	Attr1::qc ();
	
//ASSERT (! categories. empty ());
	for (const string &cat : categories)
	  QC_ASSERT (! cat. empty ());
	  
  FFOR (size_t, i, ds. objs. size ())
    if (! isMissing (i) && ! ((*this) [i] < categories. size ()))
    {
      cout << i << ": " << (*this) [i] << ' ' << categories. size () << endl;
      ERROR;
    }
    
  QC_ASSERT (categories. size () == categMap. size ());
  FFOR (size_t, i, categories. size ())
    QC_ASSERT ((* categMap. find (categories [i])). second == i);
}



string NominAttr1::getTypeStr () const
{
	string s;
	s = "Nominal";
	for (const string& cat : categories)
	  s += " " + cat;
	  
	return s;
}



#if 0
void NominAttr1::summary (ostream & /*f*/) const
{
	NOT_IMPLEMENTED;
}
#endif



void NominAttr1::str2value (size_t objNum,
                            const string &s)
{
	size_t index = missing;
  FFOR (size_t, i, categories. size ())
    if (categories [i] == s)
    {
    	index = i;
    	break;
    }
  if (index == missing)
  {
    cout << name << "[" << ds. objs [objNum] -> name << "] = !" << s << '!' << endl;
    for (const string& s1 : categories)
      cout << s1 << endl;
    ERROR;
  }
  
	values [objNum] = index;
}



Vector<NumAttr1*> NominAttr1::toNumAttr1 (Dataset &ds_arg) const
{
	ASSERT (& ds_arg == & ds);
	
  Vector<NumAttr1*> vec;
  if (categories. size () >= 2)
  {
    vec. reserve (categories. size () - 1);
    FFOR (size_t, i, categories. size () - 1)
      vec << new ExtBoolAttr1 (name + "_" + categories [i], ds_arg);
    FFOR (size_t, objNum, ds. objs. size ())
      if (! isMissing (objNum))
      {
        const size_t value = (*this) [objNum];
        ASSERT (value <= vec. size ());
        FFOR (size_t, i, vec. size ())
        {
          ExtBoolAttr1* a = var_cast (vec [i] -> asExtBoolAttr1 ());
          ASSERT (a);
          (*a) [objNum] = (value == i ? ETRUE : EFALSE);
        }
      }
  }
  return vec;
}



Vector<RealAttr1*> NominAttr1::standardize (Dataset &/*ds_arg*/,
	                                          const Sample &/*sample*/) const
{
  NOT_IMPLEMENTED;
  return Vector<RealAttr1*> ();
}



size_t NominAttr1::category2index (const string &category)
{
	ASSERT (categMap. size () == categories. size ());
	
  const CategMap::const_iterator it = categMap. find (category);
  if (it != categMap. end ())
  	return (*it). second;
  
  const size_t index = categories. size ();
  categMap [category] = index;
  categories. push_back (category);
  
  return index;
}



void NominAttr1::rebuildCategMap ()
{
  categMap. clear ();
  FFOR (size_t, i, categories. size ())
    categMap [categories [i]] = i;
}



Categorical* NominAttr1::getCategorical (const Sample &sample) const
{
  auto* cat = new Categorical ();

  MVector categ2num (categories. size (), 0);
//categ2num. putAll (0);
  for (Iterator it (sample); it ();)
    if (! isMissing (*it))
      categ2num [(*this) [*it]] += it. mult;
  categ2num. balanceRow (true, 0, 1);

  cat->probs. resize (categories. size ());
  FFOR (size_t, i, categories. size ())
    cat->probs [i] = categ2num [i];
    
  cat->setParam ();
  
  return cat;
}



bool NominAttr1::isDegenerate (const Sample &sample) const
{
  unique_ptr<const Categorical> cat (getCategorical (sample));
  return cat->getUniqueCategory () != missing;
}



void NominAttr1::deleteEmptyCategories (const Sample &sample)
{
  size_t catNum = 0;
  Vector<size_t> newCat (categories. size (), missing);
  {
    unique_ptr<const Categorical> cat (getCategorical (sample));
    ASSERT (categories. size () == cat->probs. size ());
    FFOR (size_t, i, cat->probs. size ())
      if (! nullReal (cat->probs [i]))
      {
        newCat [i] = catNum;
        catNum++;
      }
  }
    
  if (catNum == categories. size ())
    return;

  for (size_t &val : values)
    if (val != missing)
      val = newCat [val];
  
  FOR_REV (size_t, i, categories. size ())
    if (newCat [i] == missing)
      categories. eraseAt (i);
  ASSERT (catNum == categories. size ());
  
  rebuildCategMap ();
}



size_t NominAttr1::missing2category (const string &missingName)
{
  if (! existsMissing ())
    return NO_INDEX;

  const size_t missingCategory = categories. size ();
  ASSERT (missingCategory != NO_INDEX);
  EXEC_ASSERT (category2index (missingName) == missingCategory);
  for (size_t &val : values)
    if (val == missing)
      val = missingCategory;

  ASSERT (! existsMissing ());

  rebuildCategMap ();
  
  return missingCategory;
}



void NominAttr1::renameCategories ()
{
  FOR_REV (size_t, i, categories. size ())
    categories [i] = "C" + toString (i + 1);
  rebuildCategMap ();
}



NominAttr1::Dependence NominAttr1::getDependence (const Sample &sample,
                                                  size_t value,
                                                  const RealAttr1* attr) const
{
  ASSERT (attr);
  ASSERT (value < categories. size ());
    
  WeightedMeanVar cluster;
  WeightedMeanVar other;
  for (Iterator it (sample); it ();)
    if (   ! attr->isMissing (*it)
        && !       isMissing (*it)
       )
    {
      const Real r = (*attr) [*it];
      if ((*this) [*it] == value)
        cluster. add (r, it. mult); 
      else
        other. add (r, it. mult); 
    }
  
  const Dependence::NormalParam paramCluster (cluster. getMean (), cluster. getSD (), attr->decimals + 1);
  const Dependence::NormalParam paramOther   (other.   getMean (), other.   getSD (), attr->decimals + 1);
      
  Normal normCluster;
  normCluster. setParam ( paramCluster. mean
                        , paramCluster. sd
                        );
  
  Normal normOther;
  normOther. setParam ( paramOther. mean
                      , paramOther. sd
                      );

  const Prob pValue =    isNan (normCluster. loc)
                      || isNan (normOther.   loc)
                        ? 1
                        : min ( normCluster. pValue_2tail (paramOther.   mean)
                              , normOther.   pValue_2tail (paramCluster. mean)
                              );

  const Dependence dep (pValue, paramCluster, paramOther);
  dep. qc ();
  
  return dep;
}




// Attr2 

#if 0
void Attr2::saveText (ostream &os) const
{ 
  const AttrAnalysis an (this);  
  save (os, an. sample); 
}
#endif



bool Attr2::isMissing (size_t objNum) const
{
  FFOR (size_t, col, ds. objs. size ())
    if (isMissing2 (objNum, col))
      return true;
  return false;
}



void Attr2::setMissingAll ()
{
  FFOR (size_t, row, ds. objs. size ())
  FFOR (size_t, col, ds. objs. size ())
    setMissing (row, col);
}



bool Attr2::existsMissing2 (size_t &row,
                            size_t &col) const
{
  for (row = 0; row < ds. objs. size (); row++)
    if (ds. objs [row] -> mult > 0)
      for (col = 0; col < ds. objs. size (); col++)
        if (ds. objs [col] -> mult > 0)
          if (isMissing2 (row, col))
            return true;
      
  row = NO_INDEX;
  col = NO_INDEX;
  return false;
}



size_t Attr2::getWidth_max () const
{
  size_t width = string (missingStr). size ();  
  FFOR (size_t, row, ds. objs. size ())
	  FFOR (size_t, col, ds. objs. size ())
	    if (! isMissing2 (row, col))
		    maximize (width, value2str (row, col). size ());
  return width;
}



void Attr2::save (ostream &os,
                  const Sample &sample) const
{
  for (Iterator row (sample); row ();)  
  {
    for (Iterator col (sample); col ();)  
	  {
	  	if (*col)
	  		os << ' ';
	    os << value2str (*row, *col);
	  }
	  os << endl;
	}
}




// RealAttr2 

RealAttr2::RealAttr2 (const string &name_arg,
						          Dataset &ds_arg,
						          streamsize decimals_arg)
: Attr2 (name_arg, ds_arg, true)
, RealScale (decimals_arg)
, matr (ds. objs. size ())
{ 
	setMissingAll (); 
}



RealAttr2::RealAttr2 (const string &name_arg,
					            Dataset &ds_arg,
					            const RealAttr2 &from)
: Attr2 (name_arg, ds_arg, true)
, RealScale (from. decimals)
, matr (from. matr)
{
  ASSERT (ds_arg. objs. size () == from. ds. objs. size ());
}



void RealAttr2::qc () const
{
  if (! qc_on)
    return;
	Attr2::qc ();
	
  // matr[][]
  QC_ASSERT (matr. isSquare ());
  QC_ASSERT (matr. rowsSize (false) == ds. objs. size ());
#if 0
  FFOR (size_t, row, ds. objs. size ())
  FFOR (size_t, col, ds. objs. size ())
    if (   ! isMissing2 (row, col)
        && ! finite (get (row, col))
       )
    {
      cout << name << ' ' << ds. objs [row] -> name << ' ' << ds. objs [col] -> name << ' ' << get (row, col) << endl;
      ERROR; 
    }
#endif
}



bool RealAttr2::isConstant () const
{
  Real x_min = INF;
  Real x_max = -INF;
  FFOR (size_t, row, ds. objs. size ())
  FFOR (size_t, col, ds. objs. size ())
    if (! isMissing2 (row, col))
    {
    	minimize (x_min, get (row, col));
    	maximize (x_max, get (row, col));
    }
  return x_max - x_min < pow (10.0, - (int) decimals) / 2.0;
}



void RealAttr2::appendObj ()
{
	const size_t objNum_max_old = matr. rowsSize (false);
	FOR (char, b, 2)
    matr. insertRows (b, objNum_max_old, 1);
  FFOR (size_t, objNum, objNum_max_old + 1)
  {
    setMissing (objNum, objNum_max_old);
    setMissing (objNum_max_old, objNum);
  }
}



RealAttr2& RealAttr2::operator= (const RealAttr2& other)
{ 
  ASSERT (& other. ds == & ds);
  RealScale::operator= (other);
  matr = other. matr;
  return *this;
}



bool RealAttr2::existsLessThan (Real minValue,
                                size_t &row,
                                size_t &col) const
{
  for (row = 0; row < ds. objs. size (); row++)
  for (col = 0; col < ds. objs. size (); col++)
	  if (get (row, col) < minValue)
	    return true;
	row = NO_INDEX;
	col = NO_INDEX;
  return false;
}



void RealAttr2::setDiag (Real value)
{
  FFOR (size_t, row, ds. objs. size ())
    put (row, row, value);
}


size_t RealAttr2::getInfCount () const
{
  size_t n = 0;
  FFOR (size_t, row, ds. objs. size ())
    FFOR (size_t, col, ds. objs. size ())
      if (! finite (get (row, col)))
        n++;
  return n;
}



size_t RealAttr2::inf2missing ()
{
	size_t n = 0;
  FFOR (size_t, row, ds. objs. size ())
    FFOR (size_t, col, ds. objs. size ())
      if (! finite (get (row, col)))
      {
        if (verbose ())
          cout << "Infinity:" 
               << ' ' << ds. objs [row] -> name 
               << ' ' << ds. objs [col] -> name 
               << ' ' << get (row, col)
               << endl;
        setMissing (row, col);
        n++;
      }
  return n;
}




// PositiveAttr2

void PositiveAttr2::qc () const
{
  if (! qc_on)
    return;
  RealAttr2::qc ();

  // values[]
  FFOR (size_t, row, ds. objs. size ())
    FFOR (size_t, col, ds. objs. size ())
      if (! isMissing2 (row, col))
      {
        QC_ASSERT (get (row, col) >= 0);
      }
}




/////////////////////////////////////// Dataset ////////////////////////////////////

namespace 
{
  
void loadObjPair (const string &objName1, 
                  const string &objName2, 
                  StringVector &&values, 
                  Attr2 &attr)
{
  ASSERT (objName1. empty () == objName1. empty ());
  
  if (objName1. empty ())
    return;

  const size_t objNum1 = attr. ds. getName2objNum (objName1);
  const size_t objNum2 = attr. ds. getName2objNum (objName2);
  if (objNum1 == NO_INDEX)
    throw runtime_error (FUNC "Unknown object " + strQuote (objName1) + " while reading two-way attribute " + strQuote (attr. name));
  if (objNum2 == NO_INDEX)
    throw runtime_error (FUNC "Unknown object " + strQuote (objName2) + " while reading two-way attribute " + strQuote (attr. name));

  const string averageS (attr. getAverageStrValue (move (values)));
  
  attr. str2value (objNum1, objNum2, averageS); 
}

  
}



void Dataset::load (istream &is)
{
  if (! is. good ())
    throw runtime_error (FUNC "Cannot load dataset");
  
  
  string s;
  
  
  // Comments
  while (! is. eof () && s. empty ())
  {
    is >> s; 
    if (s. empty ())
    	;  // eof()
    else if (s [0] == '#')
    {
      const size_t len = 1024 * 10;  // PAR
      char line [len];  
      is. getline (line, len, '\n');
      comments << s. substr (1) + line;
      s. clear ();
    }
  }
  

  // Objects
  strUpper (s);
  ASSERT (s == "OBJNUM");
  size_t maxObjNum;
  is >> maxObjNum;
  FOR (size_t, i, maxObjNum)
  {
    auto* obj = new Obj (toString (i + 1));
    objs. push_back (obj);
  }

  bool named = false;
  is >> s;
  strUpper (s);
  if (s == "NAME")
  	named = true;
  else 
  	ASSERT (s == "NONAME");

  bool multP = false;
  is >> s;
  strUpper (s);
  if (s == "MULT")
  	multP = true;
  else 
  	ASSERT (s == "NOMULT");

  
  // attrs
  is >> s;
  strUpper (s);
  ASSERT (s == "ATTRIBUTES");
  {
    Progress prog (0, 10000);  // PAR
    string attrName;
    string dataS;
    string type;
    for (;;)
    {
      // attrName
      is >> attrName;
      if (! isAlpha (attrName [0]))
      	throw runtime_error (FUNC "Bad protein name: " + strQuote (attrName));
      dataS = attrName;
      strUpper (dataS);
      if (dataS == "DATA")
        break;
        
      prog ();
  
      // Type name
      is >> type;
      strUpper (type);
  
  
      // Attr
      Attr* attr = nullptr;
      streamsize decimals;
      if (type == "REAL")
      {
        is >> decimals;
        attr = new RealAttr1 (attrName, *this, decimals);
      }
      else if (type == "REAL2")
      {
        is >> decimals;
        attr = new RealAttr2 (attrName, *this, decimals);
      }
      else if (type == "POSITIVE")
      {
        is >> decimals;
        attr = new PositiveAttr1 (attrName, *this, decimals);
      }
      else if (type == "POSITIVE2")
      {
        is >> decimals;
        attr = new PositiveAttr2 (attrName, *this, decimals);
      }
      else if (type == "PROBABILITY")
      {
        is >> decimals;
        attr = new ProbAttr1 (attrName, *this, decimals);
      }
      else if (type == "INTEGER")
        attr = new IntAttr1 (attrName, *this);
      else if (type == "BOOLEAN")
        attr = new ExtBoolAttr1 (attrName, *this);
      else if (type == "COMPACTBOOLEAN")
        attr = new CompactBoolAttr1 (attrName, *this);
      else if (   type == "NOMINAL"
             //|| type == "ORDINAL"
              )
      {
        NominAttr1* nominAttr = /*type == "NOMINAL" 
                                  ?*/ new NominAttr1 (attrName, *this)
                                  /*: new OrdAttr1   (attrName, *this)*/;
        attr = nominAttr;
        string categoriesS;
        bool more = true;
        while (more)
        {
          getline (is, s);
          trim (s);
          more = trimSuffix (s, "\\");
          categoriesS += s;
        }
        replace (categoriesS, '\t', ' ');
        replaceStr (categoriesS, "  ", " ");
        const List<string> categories (str2list (categoriesS));
        for (const string& cat : categories)
          nominAttr->category2index (cat);
      //ASSERT (! nominAttr->categories. empty ());
      }
      if (! attr)
        throw runtime_error (FUNC "Attribute " + strQuote (attrName) + " of unknown type " + strQuote (type));
    }
  }
  if (attrs. empty ())
  	throw runtime_error (FUNC "No attributes");
    
    
  // objs, Obj::data
  {
    Progress prog (objs. size (), max<size_t> (1, 1000000 / attrs. size ()));  // PAR
    string val;
    FFOR (size_t, i, objs. size ())
    {
      prog ();
      
      // objs
      const Obj* obj = objs [i];
      if (named)
        is >> var_cast (obj) -> name;
      if (multP)
        is >> var_cast (obj) -> mult;
      
      // Obj::data
      for (const Attr* a : attrs)
        if (Attr1* attr = var_cast (a->asAttr1 ()))
        {
          is >> val;
          if (val. empty ())
            throw runtime_error (FUNC " end-of-file");
          if (val == missingStr)
            attr->setMissing (i);
          else
            attr->str2value (i, val);
        }
    }
  }
  
  
  setName2objNum ();
  
  
  // Attr2, Obj::comment
  string s1;
  while (! is. eof ())
  {
  	s. clear ();
	  is >> s;
	  if (s. empty ())
	  	break;
	  s1 = s;
	  strUpper (s1);
	  if (s1 == "COMMENT")
	  {
	    ASSERT (! objCommented ());
      is. ignore (numeric_limits<streamsize>::max(), '\n');	    
      for (const Obj* obj_ : objs)
      {
        Obj* obj = var_cast (obj_);
        getline (is, obj->comment);
        replace (obj->comment, '\t', ' ');
        trim (obj->comment);
        if (obj->comment == missingStr)
          obj->comment. clear ();
      }
	  //ASSERT (objCommented ());
	  }
	  else if (s1 == "PAIR_DATA")
	  {
	    size_t pairs = 0;
	    size_t attrsNum = 0; 
	    is >> pairs >> attrsNum;
	    VectorPtr<Attr2> attr2s;  attr2s. reserve (attrsNum);
	    string attrName;
	    FOR (size_t, col, attrsNum)
	    {
	      attrName. clear ();
	      is >> attrName;
	      if (attrName. empty ())
	        throw runtime_error (FUNC "PAIR_DATA: No attribute name");
	      const Attr* attr = name2attr (attrName);
	      if (! attr)
	        throw runtime_error (FUNC "PAIR_DATA: Unknown attribute " + strQuote (attrName));
	      const Attr2* attr2 = attr->asAttr2 ();
	      if (! attr2)
	        throw runtime_error (FUNC "PAIR_DATA: Not a two-way attribute " + strQuote (attrName));
	      attr2s << attr2;
	    }
	    string objName1, objName2, value;
	    Verbose verb (1);
	    Progress prog (pairs, 1000);  // PAR
	    FOR (size_t, row, pairs)
	    {
	      prog ();
	      objName1. clear ();
	      objName2. clear ();
	      is >> objName1 >> objName2;
        const size_t objNum1 = getName2objNum (objName1);
        if (objNum1 == NO_INDEX)
          throw runtime_error (FUNC "PAIR_DATA: Unknown object " + strQuote (objName1));
        const size_t objNum2 = getName2objNum (objName2);
        if (objNum2 == NO_INDEX)
          throw runtime_error (FUNC "PAIR_DATA: Unknown object " + strQuote (objName2));
        FOR (size_t, col, attrsNum)
        {
          is >> value;
          if (value. empty ())
            throw runtime_error (FUNC "end-of-file");
          Attr2* attr2 = var_cast (attr2s [col]);
          if (value == missingStr)
            attr2->setMissing (objNum1, objNum2);
          else
            attr2->str2value (objNum1, objNum2, value);
        }
	    }
	  }
	  else
	  {
  	  const Attr* attr = name2attr (s);
  	  if (! attr)
  	    throw runtime_error (FUNC "Two-way attribute " + strQuote (s) + " is not found");
  	  Attr2* attr2 = var_cast (attr->asAttr2 ());
  	  if (! attr2)
  	    throw runtime_error (FUNC + strQuote (s) + " is not a two-way attribute");
  	  const string twoWayAttrName ("two-way attribute " + strQuote (attr2->name));
  	  is >> s;
  	  strUpper (s);
  	  if (s == "FULL")
  	  {
        string value;
    	  FFOR (size_t, row, objs. size ())
      	  FFOR (size_t, col, objs. size ())
          {
            is >> value;
            if (value. empty ())
              throw runtime_error (FUNC "end-of-file");
            if (value == missingStr)
              attr2->setMissing (row, col);
            else
              try { attr2->str2value (row, col, value); }
                catch (const exception &e) 
                {
                  throw runtime_error (e. what () + string (": ") + twoWayAttrName 
                                       + " [" + objs [row] -> name 
                                       + ", " + objs [col] -> name 
                                       + "]");
                }
          }
      }
      else if (s == "PARTIAL")
      {
        size_t matrixObjNum;
        is >> matrixObjNum;
        Vector<StringVector> matr;  matr. resize (matrixObjNum);
        Vector<size_t> objNums;  objNums. reserve (matrixObjNum);
        FOR (size_t, i, matrixObjNum)
        {
          is >> s;
          size_t objNum = NO_INDEX;
          if (find (name2objNum, s, objNum))
            objNums << objNum;
          else
            throw runtime_error (FUNC "Unknown object " + strQuote (s) + " in " + twoWayAttrName);
          matr [i]. resize (matrixObjNum);
          string value;
          FOR (size_t, j, matrixObjNum)
          {
            is >> value;
            if (value. empty ())
              throw runtime_error (FUNC "end-of-file in " + twoWayAttrName);
            matr [i] [j] = move (value);
          }
        }
        FOR (size_t, i, matrixObjNum)
          FOR (size_t, j, matrixObjNum)
            if (matr [i] [j] == missingStr)
              attr2->setMissing (objNums [i], objNums [j]);
            else
              attr2->str2value (objNums [i], objNums [j], matr [i] [j]);
      }
      else if (s == "PAIRS")
      {
        size_t pairsNum = NO_INDEX;
        is >> pairsNum;
        if (pairsNum == NO_INDEX)
          throw runtime_error (FUNC "End-of-file reading " + twoWayAttrName);
      //is. ignore (numeric_limits<streamsize>::max (), '\n');
        skipLine (is);
        string line;
        Istringstream iss;
        string objName1;
        string objName2;
        string value;
        StringVector values;
        string objName1_old;
        string objName2_old;
        string averageS;
        FFOR (size_t, i, pairsNum)
        {
          ASSERT (objName1_old. empty () == objName2_old. empty ());
          readLine (is, line);
          iss. reset (line);
          value. clear ();
          iss >> objName1 >> objName2 >> value;
        //cout << objName1 << ' ' << objName2 << ' ' << value << endl;  
          const string atLine (" at line " + toString (i + 1));
          if (value. empty ())
            throw runtime_error (FUNC "End-of-file reading " + twoWayAttrName + atLine);
          if (objName1 < objName1_old)
            throw runtime_error (FUNC + objName1 + " is not ordered in " + twoWayAttrName + atLine);
          if (objName1 == objName1_old && objName2 < objName2_old)
            throw runtime_error (FUNC + objName2 + " is not ordered in " + twoWayAttrName + atLine);
          if (   objName1 != objName1_old
              || objName2 != objName2_old
             )
          {
          #if 0
            if (values. size () >= 2)  
            {
              cerr << attr2->name << ' '<< objName1_old << ' ' << objName2_old << ": ";
              for (const string &value : values)
                 cerr << ' ' << value;
              cerr << endl;  
            }
          #endif
            loadObjPair (objName1_old, objName2_old, move (values), *attr2);
            objName1_old = objName1;
            objName2_old = objName2;
            values. clear ();
          }
          values << move (value);
        }
        loadObjPair (objName1_old, objName2_old, move (values), *attr2);
      }
      else throw runtime_error (FUNC "Unknown representation " + strQuote (s) + " for " + twoWayAttrName);
    }
	}
    
  
  qc ();
}



Dataset::Dataset (const Eigens &eigens)
{
  FFOR (size_t, i, eigens. getDim ())
    appendObj ();
  addRealAttr1Unit ();
  auto* orderAttr = new RealAttr1 ("Order", *this);  
  auto* eigenValueFracAttr = new RealAttr1 ("EigenValueFrac", *this, 4);  // PAR
  FFOR (size_t, i, eigens. getDim ())
  {
    (*orderAttr)          [i] = log (i + 1);
    (*eigenValueFracAttr) [i] = log (eigens. explainedFrac (i));
  }
}



void Dataset::qc () const
{
  if (! qc_on)
    return;
	Root::qc ();

  for (const Obj* obj: objs)
    obj->qc ();

  Set<string> attrNames;
  for (const Attr* attr : attrs)
  {
    QC_ASSERT (& attr->ds == this);
    attr->qc ();
    if (! attr->name. empty () && attrNames. contains (attr->name))
    {
      cout << "Duplicate: " << attr->name << '!' << endl;
      ERROR;
    }
    attrNames << attr->name;
  }
  
  QC_IMPLY (! name2objNum. empty (), name2objNum. size () == objs. size ());
  QC_ASSERT (attrs. size () == name2attr_. size ());
}



void Dataset::saveText (ostream &os) const
{ 
  VectorPtr<Attr> vec;  vec. reserve (attrs. size ());
  insertAll (vec, attrs);
  const Sample sample (*this);
  sample. save (vec, os);
}



JsonArray* Dataset::comments2Json (JsonContainer* parent,
                                   const string& name) const
{
  if (comments. empty ())
    return nullptr;

  auto* jComm = new JsonArray (parent, name);
  for (const auto& s : comments)
    new JsonString (s, jComm); 
    
  return jComm;
}




// objs

void Dataset::setName2objNum ()
{
  name2objNum. clear ();
  name2objNum. reserve (objs. size ());  
  FFOR (size_t, i, objs. size ())
  {
    const string& objName = objs [i] -> name;
    if (contains (name2objNum, objName))
      throw runtime_error (FUNC "Duplicate name " + objName);
    name2objNum [objName] = i;
  }
}



size_t Dataset::getName2objNum (const string &objName) const
{
  const auto& it = name2objNum. find (objName);
	if (it == name2objNum. end ())
		return NO_INDEX;
	return it->second; 
}



size_t Dataset::appendObj (const string &objName)
{
  string realObjName = objName;
  if (realObjName. empty ())
    realObjName = toString (objs. size () + 1);

  objs << new Obj (realObjName);
  for (const Attr* attr : attrs)
    var_cast (attr) -> appendObj (); 
    
  return objs. size () - 1;
}

 

void Dataset::list2ObjNames (const VectorPtr<Named> &names)
{
  ASSERT (names. size () == objs. size ());
  FFOR (size_t, i, objs. size ())
    var_cast (objs [i]) -> name = names [i] -> name;
}



void Dataset::setMultAll (Real mult)
{
  for (const Obj* obj : objs)
    var_cast (obj) -> mult = mult;
}



void Dataset::attr2mult (const RealAttr1* attr)
{
  ASSERT (attr);
  ASSERT (& attr->ds == this);
  FFOR (size_t, i, objs. size ())
    var_cast (objs [i]) -> mult = (*attr) [i];
}



void Dataset::objInterval2Active (size_t startObjNum,
                                  size_t endObjNum)
{
  FFOR (size_t, i, objs. size ())
    var_cast (objs [i]) -> mult = between (i, startObjNum, endObjNum);
}



bool Dataset::getUnitMult () const
{
  FFOR (size_t, i, objs. size ())
    if (! eqReal (objs [i] -> mult, 1))
    	return false;
  return true;
}



bool Dataset::objCommented () const
{
  for (const Obj* obj : objs)
    if (! obj->comment. empty ())
      return true;
  return false;
}



void Dataset::deleteAttrs ()
{
  while (! attrs. empty ())
  {
    const Attr* attr = attrs. front ();
    delete attr;
  }
}



string Dataset::findNewAttrName (const string &namePrefix) const
{
  ASSERT (! namePrefix. empty ());
  
  string attrName;
  FOR (size_t, i, numeric_limits<size_t>::max ())
  {
    attrName = namePrefix + "_" + to_string (i);
    if (! name2attr (attrName))
      return attrName;
  }
  
  throw runtime_error (FUNC "No new attribute names");
}



void Dataset::addAttr (Attr* attr)
{
  ASSERT (attr);
  ASSERT (! name2attr (attr->name));
  attrs << attr;  
  name2attr_ [attr->name] = attr;

  attr->dsIt = attrs. end (); 
  attr->dsIt--;
  ASSERT (attr->dsIt != attrs. end ()); 
}





//////////////////////////// Sample //////////////////////////////

Sample::Sample (const Dataset &ds_arg)
: ds   (& ds_arg)
, mult (ds_arg. objs. size (), NaN)
{
  const VectorOwn<Obj>& objs = ds_arg. objs;
  FFOR (size_t, i, objs. size ())
    mult [i] = objs [i] -> mult;

  finish ();
}



void Sample::finish ()
{
  mult_sum = 0.0;
  nEffective = 0.0;
  for (const Real m : mult)
  {
    mult_sum += m;
    if (m)
      nEffective++;
  }
  ASSERT (mult_sum >= 0.0);
//if (! mult. empty () && ! mult_sum)
  //throw runtime_error (FUNC "mult_sum is 0");
}



void Sample::qc () const
{
  if (! qc_on)
    return;
  QC_ASSERT (ds);
  QC_ASSERT (mult. size () == ds->objs. size ());
  for (const Real m : mult)
    QC_ASSERT (m >= 0.0);
  QC_ASSERT (positive (mult_sum));  
}



Real Sample::getMaxMult () const
{
  Real n = 0.0;
  for (Iterator it (*this); it ();)  
    maximize (n, it. mult);
  return n;
}



size_t Sample::getObjNameLen_max () const
{
	size_t len = 0;
  for (Iterator it (*this); it ();)  
    maximize (len, ds->objs [*it] -> name. size ());
  return len;
}



void Sample::missing2mult (const Attr1* attr1)
{
  ASSERT (attr1);
  ASSERT (& attr1->ds == ds);
  ASSERT (attr1->ds. objs. size () == mult. size ());
  
  FFOR (size_t, i, mult. size ())
    if (attr1->isMissing (i))
    	mult [i] = 0;
}



void Sample::save (const VectorPtr<Attr> &attrs,
                   ostream &os) const
{
  ASSERT (os. good ());
  
  const bool unitMult = ds->getUnitMult ();

  size_t activeObjs = 0;    
  for (Iterator it (*this); it ();)  
    activeObjs++;

  os << "OBJNUM " << activeObjs << " name ";
  if (unitMult)
    os << " no";
  os << "mult" << endl;
  
  // Attributes
  os << "ATTRIBUTES" << endl;
  for (const Attr* attr : attrs)
  	os << "  " << attr->name << " " << attr->getTypeStr () << endl;


  // Data
  os << "DATA" << endl;
  
  Vector<size_t> widths;
  for (const Attr* attr : attrs)
    if (const Attr1* attr1 = attr->asAttr1 ())
   	  widths << attr1->getWidth_max ();
  
  const size_t objWidth = getObjNameLen_max ();
  const size_t multWidth = toString (getMaxMult ()). size ();  // integer ??
  for (Iterator it (*this); it ();)  
    {
    	os. width ((streamsize) objWidth);
      os << std::left << ds->objs [*it] -> name << " ";
    	os. width ((streamsize) multWidth);
    	if (! unitMult)
    	{
  	    os << std::right << fixed; os. precision (3); os << it. mult << " ";   // PAR
  	  }
      size_t col = 0;
      for (const Attr* attr : attrs)
  	    if (const Attr1* attr1 = attr->asAttr1 ())
  	    {
  	      os << " ";
  	      os. width ((streamsize) widths [col]);
  	      os << (attr1->rightAlign ? std::right : std::left) 
  	      	 << (attr1->isMissing (*it) ? missingStr : attr1->value2str (*it));
  	      col++;
  	    }
  	  ASSERT (col == widths. size ());
      os << endl;
    }

  for (const Attr* attr : attrs)
    if (const Attr2* attr2 = attr->asAttr2 ())
    {
    	os << attr2->name << " FULL" << endl;  // ??
    	attr2->save (os, *this);
    }

  // COMMENT
  if (ds->objCommented ())
  {
    os << "COMMENT" << endl;
    for (Iterator it (*this); it ();)  
      os << nvl (ds->objs [*it] -> comment, missingStr) << endl;
  }
}




//////////////////////////////////////////////////////////////////////////////////////

PositiveAttr2* getDist2 (const Space1<RealAttr1> &space,
                         const string &attrName,
                         Dataset &ds)
{
  ASSERT (! space. empty ());

  auto* dist = new PositiveAttr2 (attrName, ds, 4);  // PAR
  FFOR (size_t, row, ds. objs. size ())
  {
  	dist->matr. putDiag (row, 0);
    FOR (size_t, col, row)
    {
    	Real diff = 0;  
    	for (const RealAttr1* attr : space)
    	{
    	  ASSERT (& attr->ds == & ds);
    	  if (   ! attr->isMissing (row)
    	      && ! attr->isMissing (col)
    	     )
    	    diff += sqr ((*attr) [row] - (*attr) [col]); 
    	}
      dist->matr. putSymmetric (row, col, diff);
    }
  }
  
  return dist;
}



PositiveAttr2* getHammingDist (const Space1<ProbAttr1> &space,
                              const string &attrName,
                              Dataset &ds)
{
  ASSERT (! space. empty ());

#if 0
  // Attributes with most (*attr)[row] \approx 0 have too big weight
	MVector c (getDim ());  
	if (normalize)
	{
  	FFOR (size_t, attrNum, getDim ())
  	{
  	  const RealAttr1* attr = getRealAttr1 (attrNum);
  	  c [attrNum] = 0;
		  FFOR (size_t, row, ds. objs. size ())
  	    c [attrNum] += (*attr) [row];
  	  c [attrNum] /= (Real) ds. objs. size ();
  	}
  	FFOR (size_t, attrNum, getDim ())
  	{
  	  c [attrNum] -= sqr (c [attrNum]);
  	  ASSERT (isProb (c [attrNum]));
  	}
	}
#endif
		
	
  auto* dist = new PositiveAttr2 (attrName, ds, 4);  // PAR
  FFOR (size_t, row, ds. objs. size ())
  {
  	dist->matr. putDiag (row, 0);
    FOR (size_t, col, row)
    {
    	Real diff = 0;  
    	for (const ProbAttr1* attr : space)
    	{
    	  ASSERT (& attr->ds == & ds);
    	  if (   attr->isMissing (row)
    	      || attr->isMissing (col)
    	     )
    	    continue;
    	  // p_i = P(attr for object i = 1) 
    	  const Prob p1 = (*attr) [row];  
    	  const Prob p2 = (*attr) [col];  
    	  Real diff1 = p1 + p2 - 2 * p1 * p2;
    	  diff += diff1;
    	}
      dist->matr. putSymmetric (row, col, diff);
    }
  }
  
  return dist;
}




RealAttr2* getSimilarity (const Space1<RealAttr1> &space,
                          const string &attrName,
                          Dataset &ds)
{
  ASSERT (! space. empty ());

  auto* sim = new RealAttr2 (attrName, ds, 4);  // PAR
  Progress prog (ds. objs. size (), max<size_t> (1, 10000 / space. size ()));  // PAR
  FFOR (size_t, row, ds. objs. size ())
  {
    prog ();
    FFOR_START (size_t, col, row, ds. objs. size ())
    {
    	Real s = 0;  
      for (const RealAttr1* attr : space)
    	{
    	  ASSERT (& attr->ds == & ds);
        // Missings are treated as 0, which should be the attribute mean
        const Real rowValue = (*attr) [row]; 
        if (attr->isMissingValue (rowValue))
          continue;
    	  const Real colValue = (*attr) [col]; 
    	  if (attr->isMissingValue (colValue))
    	    continue;  
    	  s += rowValue * colValue; 
    	}
      sim->matr. putSymmetric (row, col, s);
    }
  }
  
  sim->matr. psd = true;
  
  return sim;
}





// Distribution

void Distribution::qc () const 
{
  if (! qc_on)
    return;
  Named::qc ();

  if (getParamSet ())
  {  	
	  QC_ASSERT (getDim ());
	//QC_IMPLY (analysis, getDim () == analysis->space. size ());  
	}
}



size_t Distribution::getModeObjNum () const
{
  size_t objNum = NO_INDEX;
  Real logPdf_max = -INF; 
  for (Iterator it (getAnalysisCheck () -> sample); it ();)  
  {
  	data2variable (*it);
  	if (maximize (logPdf_max, logPdfVariable ()))
  	  objNum = *it;
  }
  ASSERT (objNum != NO_INDEX);
  return objNum;
}



void Distribution::simulate (Dataset &ds,
                             size_t objsSize)
{
  ASSERT (ds. empty ());

  auto* an_ = createAnalysis (ds);
  
  FOR (size_t, i, objsSize)  
  {
  	EXEC_ASSERT (ds. appendObj (toString (i + 1)) == i);
  	randVariable ();
  	variable2data (i);
  }
  
  an_->sample = Sample (ds);
  an_->qc ();
}



void Distribution::info_meanVar_MC (uint sampleSize,
                                    Real &mean,
                                    Real &var) const
{
	// sampleSize --> errorMax ??
  MeanVar mv;
  FOR (uint, i, sampleSize)
  {
  	randVariable ();
    mv << logPdfVariable ();
  }
  mean = - mv. getMean ();
  var  =   mv. getVar (true);  
}



void Distribution::info_meanVar (uint sampleSize,
                                 Real &mean,
                                 Real &var) const
{
	mean = getInfoMean ();
	var  = getInfoVar ();
	if (   ! isNan (mean) 
		  && ! isNan (var)
		 )
		return;

//if (! isNan (var))  ??

  Real mean_est;
  info_meanVar_MC (sampleSize, mean_est, var);		
  
  if (isNan (mean))
	  mean = mean_est;
	else
	{
	  const Real bias = mean - mean_est;
		var += sqr (bias);
  }
}



Real Distribution::getLogLikelihood () const
{
	Sum s;
  for (Iterator it (getAnalysisCheck () -> sample); it ();)  
  {
  	data2variable (*it);
  	s. add (it. mult * logPdfVariable ());
 /*if (! finite (s))
  	{
  	  cout << an. space. ds. objs [*it] -> name << endl;
  	  ERROR;
  	}*/
  }
/*if (verbose ())
    cout << "logLikelihood = " << s << endl;  */
  return s. get ();
}



Normal Distribution::getEntropyDistribution () const
{
  Unverbose unv;

	// mean, sigma
	Real mean, var;
  info_meanVar (100000/*PAR*/, mean, var); 
  const Real sigma = sqrt (var / getAnalysisCheck () -> sample. mult_sum);

  Normal normal;
  normal. setParam (mean, sigma);
  if (verbose ())
	  cout << "Entropy distribution: mean = " << mean << ", SD = " << sigma << endl;
  
  return normal;
}



Prob Distribution::getFitness_entropy (Real entropy_est_best) const
{	
	return getEntropyDistribution (). cdf (entropy_est_best); 
}



Prob Distribution::getFitness_entropy () const
{ 
  Unverbose unv;

	const Real entropy_est = getEntropy_est ();
  if (verbose ())
	  cout << "Entropy^: " << entropy_est << endl;

  return getEntropyDistribution (). pValue_2tail (entropy_est);
}



size_t Distribution::getWeakestObjNum () const
{
  size_t objNum = NO_INDEX; 
  Real pdf_min = INF;
  for (Iterator it (getAnalysisCheck () -> sample); it ();)  
  {
  	data2variable (*it);
    if (minimize (pdf_min, pdfVariable ()))
    	objNum = *it;
  }
  ASSERT (objNum != NO_INDEX);

  return objNum;
}




// Bernoulli

void Bernoulli::qc () const
{
  if (! qc_on)
    return;
	Distribution::qc ();
		
  QC_IMPLY (getParamSet (), isProb (p));
}



Analysis1* Bernoulli::createAnalysis (Dataset &ds)
{
  auto* attr = new ExtBoolAttr1 ("X", ds);  
  const Sample sm (ds);
  auto* analysis_ = new An (sm, *attr);
  analysis = analysis_;
  
  return analysis_;
}



bool Bernoulli::similar (const Distribution &distr,
                         Real delta) const 
{ 
	if (const Bernoulli* b = distr. asBernoulli ())
	  return eqReal (p, b->p, delta);
	// Convert to Bernoulli !??
	if (const Categorical* cat = distr. asCategorical ())
	  return    cat->probs. size () == 2 
	         && eqReal (p, cat->probs [1], delta);
	if (const Binomial* bin = distr. asBinomial ())
	  return    bin->n == 2 
	         && eqReal (p, bin->p, delta);

  return false;
}




// Categorical

void Categorical::qc () const
{
  if (! qc_on)
    return;
	Distribution::qc ();
		
  if (getParamSet ())
  {		
	  Prob s = 0;
	  for (const Prob p : probs)
	    s += p;
	  QC_ASSERT (eqReal (s, 1));
	
		QC_ASSERT (eqReal (probSum [probSum. size () - 1], 1));	
  }
}



void Categorical::saveText (ostream& os) const
{ 
	os << name << "(" << probs. size () << ") ";
	bool first = true;
  for (const Prob p : probs)
  {
  	if (! first)
  		os << ", ";
    os << p;
    first = false;
  }
}



void Categorical::setParamFunc ()
{
  balanceProb ();
		
  // probSum
	probSum. resize (probs. size ());
	Prob s = 0;
	FFOR (size_t, i, probs. size ())
  {
	  s += probs [i];
	  probSum [i] = s;
	}
	probSum. searchSorted = true;
}



void Categorical::setParam (const DiscreteDistribution &distr,
		                        int x_from,
		                        int x_to)
{
	ASSERT (x_from <= x_to);
	
	const int n = x_to - x_from + 1;
	
	probs. resize ((size_t) n);
	FOR (int, i, n)
	  probs [(size_t) i] = distr. pmf (x_from + i);

	setParam ();
}



Analysis1* Categorical::createAnalysis (Dataset &ds)
{
  auto* attr = new NominAttr1 ("X", ds);  
  FFOR (size_t, i, probs. size ())
    attr->category2index ("C" + toString (i + 1));

  const Sample sm (ds);
  auto* analysis_ = new An (sm, *attr);    
  analysis = analysis_;
  
  return analysis_;
}



bool Categorical::similar (const Distribution &distr,
	                         Real delta) const
{
	if (const Categorical* cat = distr. asCategorical ())
	{
	  if (cat->probs. size () != probs. size ())
	  	return false;
	  FFOR (size_t, i, probs. size ())
	    if (! eqReal (probs [i], cat->probs [i], delta))
	    	return false;
	  return true;
	}

  return false;
  
}



void Categorical::estimate ()
{ 
  ASSERT (analysis);
  
	const NominAttr1& attr = analysis->attr;

  // probs
  for (Prob& p : probs)
    p = 0;
  for (Iterator it (analysis->sample); it ();)  
   	probs [attr [*it]] += it. mult;
    
  setParam ();
}



void Categorical::balanceProb ()
{
  Real s = 0;
  for (const Prob p : probs)
    s += p;
  for (Prob& p : probs)
    p /= s;
}



size_t Categorical::getUniqueCategory () const
{
  size_t unique = NO_INDEX;
  FFOR (size_t, i, probs. size ())
    if (! nullReal (probs [i]))
    {
      if (unique == NO_INDEX)
        unique = i;
      else
        return NO_INDEX;
    }
  return unique;
}




// UniDistribution

void UniDistribution::qc () const
{
  if (! qc_on)
    return;
  Distribution::qc ();

  if (getParamSet ())
  {
    QC_ASSERT (loBound <= hiBound);
    QC_ASSERT (isProb (p_supp));
    QC_ASSERT (isProb (p_ltSupp));
    QC_ASSERT (isProb (p_supp + p_ltSupp));
  }
  
#if 0
  // estimate() takes care of this
  if (analysis)
  {
  	const RealAttr1* attr = analysis->space. getRealAttr1 ();
  	const Real lo = getLoBoundEffective ();
  	const Real hi = getHiBoundEffective ();
  	if (   lo != -INF
  		  && hi !=  INF
  		 )
		  for (Iterator it (*analysis); it ();)  
		  {
		   	QC_ASSERT (geReal ((*attr) [*it], lo));
		   	QC_ASSERT (leReal ((*attr) [*it], hi));  
		   	  // Binomial::n in a Mixture ??
		  }
  }
#endif
}




// DiscreteDistribution

void DiscreteDistribution::qc () const
{
  if (! qc_on)
    return;
	UniDistribution::qc ();
		
  if (getParamSet ())
  {
  	QC_IMPLY (loBound != -INF, isInteger (loBound));
  	QC_IMPLY (hiBound !=  INF, isInteger (hiBound));
  }
}



Analysis1* DiscreteDistribution::createAnalysis (Dataset &ds)
{
  auto* attr = new IntAttr1 ("X", ds);  
  const Sample sm (ds);
  auto* analysis_ = new An (sm, *attr);
  analysis = analysis_;
  
  return analysis_;
}


int DiscreteDistribution::real2int (Real x)
{
  ASSERT (isInteger (x));
	return (int) round (x);
}



#if 0
Prob DiscreteDistribution::cdfDiscrete_ (int x) const 
{
	Prob s = 0;
	FOR_REV (int, i, x + 1)  // Usually pmf(i) > pmf(i+1)
	  s += pmf (i);
	  	 
  return s;
}
#endif




// Binomial

void Binomial::qc () const
{
  if (! qc_on)
    return;
	DiscreteDistribution::qc ();		
		
	QC_ASSERT (n >= 0);
  if (getParamSet ())
  {
    QC_ASSERT (n > 0);
	  QC_ASSERT (isProb (p));
	}

	bernoulli. qc ();
}



void Binomial::setParamFunc ()
{ 
  n_re = (Real) n;
	bernoulli. setParam (p);
	lnFacN = lnFactorial ((uint) n);
	lnP = log (p);
	lnPCompl = log (1 - p);
	
  const int lo = real2int (getLoBoundEffective ());
  const int hi = real2int (getHiBoundEffective ());
	p_ltSupp = lo ? cdfDiscrete_ (lo - 1) : 0;
	const Real p_gtSupp = hi == n ? 0 : (1 - cdfDiscrete_ (hi));
	p_supp = 1 - (p_ltSupp + p_gtSupp);
}



void Binomial::estimate ()
{
  ASSERT (analysis);
	const IntAttr1& attr = analysis->attr;


#if 0
  if (! nFixed && ! n)
	  for (Iterator it (*analysis); it ();)  
    	maximize (n, real2int ((*attr) [*it]));
  if (! n )
  	return;
#endif


  // p
	if (stdBounds ())
  {
		Real s = 0;
		Real sCompl = 0;
	  for (Iterator it (analysis->sample); it ();)  
	  {
	  	const An::Value x = attr [*it];
	  	ASSERT (x >= 0);
	  	ASSERT (x <= n);
	  	s      += it. mult * x;
	  	sCompl += it. mult * (n - x);
	  }  
	  const Real a = sCompl / s;	  
    setParam (n, 1 / (a + 1));
	}
	else
	{
		Prob p_ = 0.5;  // Init
		// Convergence ??
		do
		{
			setParam (n, p_);
			if (   eqReal (p, 0)
				  || eqReal (p, 1)
				 )
				break;
			Real c;
			{
				const int loBound_ = real2int (getLoBoundEffective ());
				const int hiBound_ = real2int (getHiBoundEffective ());
				Real sum = 0;
				Real sumLog = 0;
				FFOR_START (int, x, loBound_, hiBound_ + 1)
				{
					const Real y = pmf_ (x);
				//if (verbose ())
					//cout << "pmf(" << x << ")=" << y << endl;
		      sumLog += y * ((Real) x / p - (Real) (n - x) / (1 - p));
				  sum    += y;
			    ASSERT (! isNan (sumLog));
				}
				ASSERT (geReal (sum, 0));
				c = sumLog / sum;
				ASSERT (! isNan (c));
			}
			Real x_ave;
			{
				Real mx = 0;
				Real m  = 0;
			  for (Iterator it (analysis->sample); it ();)  
			  {
			  	const An::Value x = attr [*it];
			  	mx += it. mult * x;
			  	m  += it. mult;
			  }  
			  if (nullReal (m))
			  {
			  	setParam (n, 0);
			  	break;
			  }
			  x_ave = mx / m;
			  ASSERT (x_ave >= 0);
				ASSERT (leReal (x_ave, n_re));
				minimize (x_ave, n_re);
			}
		  if (fabs (c) <= 0.1)  // PAR
		  	p_ = x_ave / n_re;
		  else
		  {
				// c*p_^2 - (n+c)*p_ + x_ave = 0
			  const Real d = sqr (n_re + c) - 4 * c * x_ave;  // Determinant
			    // If c >= 0 then d = (n+c)^2 - 4*c*x_ave >= (n+c)^2 - 4*c*n = (n-c)^2 >= 0
			    // If c <  0 then d = (n+c)^2 - 4*c*x_ave <= (n+c)^2 - 4*c*n = (n-c)^2 
			  ASSERT (d >= 0);
				const Real d_sqrt = sqrt (d);
				IMPLY (c >= 0, geReal (d_sqrt, fabs (n_re - c)));
				IMPLY (c >= 0, leReal (d_sqrt, fabs (n_re + c)));
				IMPLY (c <  0, leReal (d_sqrt, fabs (n_re - c)));
				IMPLY (c <  0, geReal (d_sqrt, fabs (n_re + c)));
				  // If c >= 0 & n-c >= 0 then n+c - d_sqrt <= n+c - (n-c) = 2*c
				  // If c >= 0 & n-c <  0 then n+c + d_sqrt <= n+c + (c-n) = 2*c
				  // If c <  0 & n-c >= 0 then n+c - d_sqrt >= n+c - (n-c) = 2*c
				  // If c <  0 & n-c <  0 then n+c + d_sqrt >= n+c + (c-n) = 2*c
				p_ = (n_re + c + ((n_re - c >= 0) ? -1 : 1) * d_sqrt) / (2 * c);  // The other solution: isProb() => second local MLE
				ASSERT (isProb (p_));
			}
		//cout << "c=" << c << " p_=" << p_ << endl;  
		}
		while (fabs (p - p_) > 1e-5);  // PAR
	}


#if 0
  // --> clusterSizeDist.cpp ??
 	// n, p
  if (! nFixed)
  {
	  int n_best = n;
	  Prob p_best = p;
	  Real entropyEst_best = INF;
	  while (minimizeEntropy (entropyEst_best))
	  {
		 	n_best = n;
		 	p_best = p;
		 	const int n_new = n + 1;
		 	const Prob p_new = pFixed ? p : (p * n / n_new);
	    setParam (n_new, p_new);
	  }
	  setParam (n_best, p_best);
  }
#endif
}



Real Binomial::logPmf_ (int x) const
{ 
	if (x > n)
		return -INF;
	return   lnFacN - lnFactorial ((uint) x) - lnFactorial ((uint) (n - x)) 
	       + multiplyLog ((Real) x, lnP)
	       + multiplyLog ((Real) (n - x), lnPCompl); 
}



Prob Binomial::cdfDiscrete_ (int x) const
{ 
	ASSERT (x <= n);
	
	int x_start = 0;
	int x_end = x;
	if (x_end > n / 2)
	{
		x_start = x_end + 1;
		x_end   = n;
		ASSERT (x_start);
	}

	Sum sum; 
	FFOR_START (int, i, x_start, x_end + 1)
	  sum. add (pmf_ (i));
	Prob s = min (sum. get (), 1.0);
	 
	if (x_start)
		s = 1 - s;

  return s;
}



int Binomial::randDiscrete_ () const
{
	int r = 0;
  FOR (int, i, n)
  {
    bernoulli. randVariable ();
    r += bernoulli. variable;
  }
  return r;
}




// UniformDiscrete

void UniformDiscrete::qc () const
{
  if (! qc_on)
    return;
	DiscreteDistribution::qc ();
		
  if (getParamSet ())
  {
  	QC_ASSERT (stdBounds ());  // ??
  	QC_ASSERT (min <= max);
  }
}




// Geometric

void Geometric::qc () const
{
  if (! qc_on)
    return;
	DiscreteDistribution::qc ();

  if (getParamSet ())
  {
    QC_ASSERT (stdBounds ());  // ??
    QC_ASSERT (isProb (p));
  }
}



void Geometric::estimate ()
// TOTEST
{
  ASSERT (analysis);
	const IntAttr1& attr = analysis->attr;
	
	Real s = 0;
	Real n = 0;
  for (Iterator it (analysis->sample); it ();)  
  {
  	const An::Value x = attr [*it];
  	s += it. mult * x;
  	n += it. mult;
  }
  
  setParam (n / (s + n));
}



Real Geometric::logPmf_ (int x) const
{
	ASSERT (x >= 1);
	return lnP + ((Real) x - 1) * lnPCompl;
}




// Zipf

void Zipf::qc () const
{
  if (! qc_on)
    return;
	DiscreteDistribution::qc ();

  if (getParamSet ())
  {
    QC_ASSERT (alpha > 1);
  //IMPLY (hiBound == INF, alpha > 1);  
  	QC_ASSERT (c > 0);
  }
  
	cat. qc ();
}



void Zipf::setParamFunc ()
{
 	const uint lo = (uint) real2int (getLoBoundEffective ());
	const Real zeta_gtSupp = hiBound == INF ? 0 : zeta (alpha, (uint) real2int (hiBound));
	const Real zeta_geSupp = zeta (alpha, lo);
  const Real zeta_ltSupp = lo ? zeta (alpha, 1, lo - 1) : 0;
	
	c = 1 / (zeta_ltSupp + zeta_geSupp);
	lnC = log (c);

 	p_supp = (zeta_geSupp - zeta_gtSupp) * c;
 	p_ltSupp = zeta_ltSupp * c;
 	
//cout << lo << " " << alpha << " " << zeta (alpha, lo) << " " << zeta_gtSupp << " " << zeta_geSupp << " " << zeta_ltSupp << endl;  
 	
 	if (hiBound < INF)
 		cat. setParam (*this, (int) lo, real2int (hiBound));
}



namespace {

struct Zipf_estimate_Func : Func1
{
  uint from;
  uint to;
    // 0 <=> unlimited
  Real goal;

  Real f (Real x) 
  { 
  	const Real alpha = - x;
  	if (to)
      return zetaLn (alpha, from, to) / zeta (alpha, from, to) - goal;
    return zetaLn (alpha, from) / zeta (alpha, from) - goal;
  }  
};

}



void Zipf::estimate ()
{
  ASSERT (analysis);  
  
  Unverbose unv;

	const IntAttr1& attr = analysis->attr;
	
	//   \sum_x m_x log x \over \sum_x m_x 
	// = \sum_i=1^\infty i^alpha log i \over \sum_i=1^\infty i^alpha
	
  Zipf_estimate_Func zf;
  {
	  zf. from = (uint) real2int (getLoBoundEffective ());
  	zf. to   = (uint) (hiBound == INF ? 0 : real2int (hiBound));
		Real s = 0;
		Real w = 0;
	  for (Iterator it (analysis->sample); it ();)  
	  {
	  	const An::Value x = attr [*it];
	  	s += it. mult * log (x);
	  	w += it. mult;
	  }
	  zf. goal = s / w;  // = \sum_x m_x log x \over \sum_x m_x 
  }
#if 0
  CONST_ITER (Data, it, data)
    cout << (*it). first << " " << (*it). second << endl;
#endif
  const Real alpha_max = -1 /*(hiBound == INF ? -1 : 0)*/ - 1e-5;  // PAR
  const Real alpha_min = -10;  // PAR
  if (verbose ())
  {
	  cerr << "goal = " << zf. goal << endl;  
	  cerr << "f(alpha_min) = " << zf. f (alpha_min) << endl;  
	  cerr << "f(alpha_max) = " << zf. f (alpha_max) << endl;  
	  cerr << "f(-0.5)  = " << zf. f (-0.5) << endl;  
	}
  Real alpha_;  
  if      (geReal (zf. f (alpha_min), 0))  // ??
  	alpha_ = - alpha_min;
  else if (leReal (zf. f (alpha_max), 0))  // ??
  	alpha_ = - alpha_max;  
  else
    alpha_ = - zf. findZero (alpha_min, alpha_max, 0.001);  // PAR
  
#if 0
  FOR (uint, i, 30)
  {
  	setParam (alpha_ + 0.001 * i);
  	cout << alpha << " " << getEntropy_est () << endl;
  }
#endif   
  setParam (alpha_);
}



Real Zipf::logPmf_ (int x) const
{ 
  ASSERT (x > 0);
	return lnC - alpha * log (x); 
}



int Zipf::randDiscrete_ () const
{ 
 	if (hiBound == INF)
 		NOT_IMPLEMENTED;
  cat. randVariable ();
	return (int) cat. variable + real2int (getLoBoundEffective ());
}




// ContinuousDistribution

Analysis1* ContinuousDistribution::createAnalysis (Dataset &ds)
{
  auto* attr = new RealAttr1 ("X", ds);  
  const Sample sm (ds);
  auto* analysis_ = new An (sm, *attr);
  analysis = analysis_;
  
  return analysis_;
}




// LocScaleDistribution

void LocScaleDistribution::qc () const
{
  if (! qc_on)
    return;
  ContinuousDistribution::qc ();

  QC_IMPLY (! isNan (scale), scale >= 0.0);
}




// Normal

const Real Normal::coeff = log (2.0 * pi);
	
	

void Normal::qc () const
{
  if (! qc_on)
    return;
	LocScaleDistribution::qc ();

	QC_IMPLY (getParamSet (), stdBounds ());  // ??
}
	
	

void Normal::estimate ()
{
  ASSERT (analysis);
	const NumAttr1& attr = analysis->attr;
	
	Real n = 0;
	Real s = 0;
	Real s2 = 0;
  for (Iterator it (analysis->sample); it ();)  
    if (! attr. isMissing (*it))
    {
    	const An::Value x = attr. getReal (*it);
    	n += it. mult;
    	s += it. mult * x;
    	s2 += it. mult * sqr (x);
    }
  
  const Real mu_ = s / n;
  
  setParam (mu_, sqrt (s2 / n - sqr (mu_)));
}



Real Normal::rand_ () const
// Numeric Recipes in C, p. 289
{
  Real a, b, d2;
  do
  {
    a = 2 * randProb () - 1;
    b = 2 * randProb () - 1;
    d2 = sqr (a) + sqr (b);
  }
  while (d2 >= 1 || d2 == 0);
  const Real f = sqrt (-2 * log (d2) / d2);
  return unstnd (a * f);
}



Real Normal::getQuantile (Prob p) const 
// https://gist.github.com/kmpm/1211922/:

/// Original C++ implementation found at http://www.wilmott.com/messageview.cfm?catid=10&threadid=38771
/// C# implementation found at http://weblogs.asp.net/esanchez/archive/2010/07/29/a-quick-and-dirty-implementation-of-excel-norminv-function-in-c.aspx
/*
 *     Compute the quantile function for the normal distribution.
 *
 *     For small to moderate probabilities, algorithm referenced
 *     below is used to obtain an initial approximation which is
 *     polished with a final Newton step.
 *
 *     For very large arguments, an algorithm of Wichura is used.
 *
 *  REFERENCE
 *
 *     Beasley, J. D. and S. G. Springer (1977).
 *     Algorithm AS 111: The percentage points of the normal distribution,
 *     Applied Statistics, 26, 118-121.
 *
 *      Wichura, M.J. (1988).
 *      Algorithm AS 241: The Percentage Points of the Normal Distribution.
 *      Applied Statistics, 37, 477-484.
 */
{
  ASSERT (isProb (p));
 
  if (p == 0)
    return -INF;
  if (p == 1)
    return INF;
  if (scale == 0)
    return loc;

  const Real q = p - 0.5;
  /*-- use AS 241 --- */
  /* double ppnd16_(double *p, long *ifault)*/
  /*      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
          Produces the normal deviate Z corresponding to a given lower
          tail area of P; Z is accurate to about 1 part in 10**16.
  */
  Real val = NaN;
  if (fabs(q) <= .425)
  {
    // 0.075 <= p <= 0.925 
    const Real r = .180625 - q * q;
    val = q * (((((((r * 2509.0809287301226727 +
                      33430.575583588128105) * r + 67265.770927008700853) * r +
                    45921.953931549871457) * r + 13731.693765509461125) * r +
                  1971.5909503065514427) * r + 133.14166789178437745) * r +
                3.387132872796366608)
           / (((((((r * 5226.495278852854561 +
                    28729.085735721942674) * r + 39307.89580009271061) * r +
                  21213.794301586595867) * r + 5394.1960214247511077) * r +
                687.1870074920579083) * r + 42.313330701600911252) * r + 1);
  }
  else
  { 
    // closer than 0.075 from {0,1} boundary 

    /* r = min(p, 1-p) < 0.075 */
    Real r = (q > 0) ? 1 - p : p;
    r = sqrt (-log (r));
    /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */

    if (r <= 5)
    { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
      r += -1.6;
      val = (((((((r * 7.7454501427834140764e-4 +
                 .0227238449892691845833) * r + .24178072517745061177) *
               r + 1.27045825245236838258) * r +
              3.64784832476320460504) * r + 5.7694972214606914055) *
            r + 4.6303378461565452959) * r +
           1.42343711074968357734)
          / (((((((r *
                   1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                  r + .0151986665636164571966) * r +
                 .14810397642748007459) * r + .68976733498510000455) *
               r + 1.6763848301838038494) * r +
              2.05319162663775882187) * r + 1);
    }
    else
    { /* very close to  0 or 1 */
      r += -5;
      val = (((((((r * 2.01033439929228813265e-7 +
                 2.71155556874348757815e-5) * r +
                .0012426609473880784386) * r + .026532189526576123093) *
              r + .29656057182850489123) * r +
             1.7848265399172913358) * r + 5.4637849111641143699) *
           r + 6.6579046435011037772)
          / (((((((r *
                   2.04426310338993978564e-15 + 1.4215117583164458887e-7) *
                  r + 1.8463183175100546818e-5) * r +
                 7.868691311456132591e-4) * r + .0148753612908506148525)
               * r + .13692988092273580531) * r +
              .59983220655588793769) * r + 1);
    }
    if (q < 0.0)
      val = -val;
  }

  return loc + scale * val;
}
  
  


// Exponential

void Exponential::qc () const
{
  if (! qc_on)
    return;
	LocScaleDistribution::qc ();

	if (getParamSet ())
	{
	  QC_ASSERT (stdBounds ());  // ??
  	QC_ASSERT (loc >= 0);
  	QC_ASSERT (scale == loc);
  }
}
	
	

void Exponential::estimate ()
{
  ASSERT (analysis);
	const NumAttr1& attr = analysis->attr;
	
	Real n = 0;
	Real s = 0;
  for (Iterator it (analysis->sample); it ();)  
    if (! attr. isMissing (*it))
    {
    	const An::Value x = attr. getReal (*it);
    	n += it. mult;
    	s += it. mult * x;
    }
  
  setParam (s / n);
}




// Cauchy

void Cauchy::qc () const
{
  if (! qc_on)
    return;
	LocScaleDistribution::qc ();

	QC_IMPLY (getParamSet (), stdBounds ());  // ??
}



namespace {

struct CauchyScaleFunc : Func1
{
	const Cauchy::An& an;
  Real loc;
  
  CauchyScaleFunc (const Cauchy::An& an_arg,
                   Real loc_arg)
    : an (an_arg)
    , loc (loc_arg)
    {}

  Real f (Real x) 
  { 
  	const NumAttr1& attr = an. attr; 
  	const Real x2 = sqr (x);
		Real n = 0;
		Real s = 0;
	  for (Iterator it (an. sample); it ();)  
	  {
	  	const Cauchy::An::Value y = attr. getReal (*it);
	  	n += it. mult;
	  	s += 1 / (x2 + sqr (y - loc)) * it. mult;
	  }
	  
	  return (nullReal (x2) ? 0 : (x2 * s)) - n / 2;
  }  
};

}



void Cauchy::estimate ()
{
  ASSERT (analysis);
	const NumAttr1& attr = analysis->attr;
	
  setParam (NaN, NaN);
	Real loc_ = attr. getMedian (analysis->sample);
	Real scale_ = INF;
	
	// Convergence ??
	do
	{
	  setParam (loc_, scale_);

		// scale_
		{
			CauchyScaleFunc f (*analysis, loc_);
			Real minX, maxX;
			attr. getMinMax (analysis->sample, minX, maxX);
			const Real maxScale = max ( fabs (loc_ - minX)
						                    , fabs (loc_ - maxX)
						                    );
		  scale_ = f. findZero (0, maxScale, 1e-4);  // PAR
		  ASSERT (positive (scale_));
		}
		
	  // loc_	
	  {
			Real n = 0;
			Real s = 0;
		  for (Iterator it (analysis->sample); it ();)  
		  {
		  	const An::Value x = attr. getReal (*it);
		  	const Real weight = it. mult / (sqr (scale_) + sqr (x - loc_));
		  	n += weight;
		  	s += weight * x;
		  }
		  loc_ = s / n;
		}
	}
	while (   fabs (loc   - loc_)   > 1e-5   // PAR
	       || fabs (scale - scale_) > 1e-5   // PAR
	      );
}	
	



// Beta1

void Beta1::qc () const
{
  if (! qc_on)
    return;
	ContinuousDistribution::qc ();

  if (getParamSet ())
  {
  	QC_ASSERT (stdBounds ());  // ??
  	QC_ASSERT (alpha > 0);
  }
}



void Beta1::estimate ()
{
  ASSERT (analysis);
	const NumAttr1& attr = analysis->attr;
	
	Real s = 0;
	Real n = 0;
  for (Iterator it (analysis->sample); it ();)  
  {
  	const An::Value x = attr. getReal (*it);
  	ASSERT (isProb (x));
  	s += it. mult * log (x);
  	n += it. mult;
  }
  
  setParam (- n / s);
}




// UniKernel

void UniKernel::Point::qc () const
{
  if (! qc_on)
    return;
  QC_ASSERT (! isNan (value));
  QC_ASSERT (positive (mult));
}



void UniKernel::qc () const
{
  if (! qc_on)
    return;
  ContinuousDistribution::qc ();
    
    
  QC_ASSERT (analysis);
  QC_ASSERT (positive (analysis->sample. mult_sum));
  QC_ASSERT (points.  size () <= analysis->sample. size ());
  QC_ASSERT (multSum. size () <= analysis->sample. size ());
  QC_ASSERT (points. size () == multSum. size ());
  QC_ASSERT (! points. empty ());


  if (getParamSet ())
  {
    QC_ASSERT (! isNan (attr_min));
    QC_ASSERT (! isNan (attr_max));
    QC_ASSERT (getRange () >= 0);
    QC_ASSERT (isProb (uniform_prob));
    
    Real prevValue = NaN;
    for (const Point& p : points)
    {
      p. qc ();
      IMPLY (! isNan (prevValue), prevValue <= p. value);      
      prevValue = p. value;
    }

    Real prevMult = NaN;
    for (const Real& mult : multSum)
    {
      QC_ASSERT (positive (mult));
      QC_IMPLY (! isNan (prevMult), prevMult < mult);      
      prevMult = mult;
    }

    QC_ASSERT (halfWindow >= 0);
    QC_ASSERT (positive (halfWindow) == positive (getRange ()));
  }
}



void UniKernel::estimate ()
{
  uniform_prob = 0;
  halfWindow = 0;
  
  const Real sd = setPoints ();

  if (! nullReal (getRange ()))
  {
    // uniform_prob, halfWindow  
    // Strange optimum: large uniform_prob, small halfWindow ??
    Real halfWindow_lo = sd * 0.01;  // PAR
    Real halfWindow_hi = getRange () / 2;  // otherwise uniform distribution will have a higher likelihood
    Real step = pow (halfWindow_hi / halfWindow_lo, 0.01);  // PAR
    if (verbose ())
      cout << halfWindow_lo << " <= halfWindow <= " << halfWindow_hi << "  step = " << step << endl;
    Progress prog;
    while (greaterReal (halfWindow_hi, halfWindow_lo, /*getRange () * */ 1e-4))  // PAR
    {
      estimate_ (halfWindow_lo, halfWindow_hi, step);
      halfWindow_lo = halfWindow / pow (step, 10);
      halfWindow_hi = halfWindow * pow (step, 10);
      step = sqrt (step);
      prog (toString (halfWindow_hi - halfWindow_lo));
    }
    ASSERT (positive (halfWindow));
  }  

  ASSERT (getParamSet ());
}



Real UniKernel::setPoints ()
{
	const NumAttr1& attr = getAttr ();

  points. clear ();
  for (Iterator it (analysis->sample); it ();)  
    if (   ! attr. isMissing (*it)
        && positive (it. mult)
        && ! attr. isMissing (*it)
       )
    	points << Point (attr. getReal (*it), it. mult);
  ASSERT (! points. empty ());
  points. sort ();

  multSum. clear ();
  WeightedMeanVar mv;
  {
    Sum s;
    for (const Point& p : points)
    {
      s. add (p. mult);
      multSum << s. get ();
      mv. add (p. value, p. mult);
    }
  }
  if (verbose ())
    cout << "SD = " << mv. getSD () << endl;  
    
  var_cast (analysis) -> sample. mult_sum = multSum. back ();

  // attr_min, attr_max
  attr. getMinMax (analysis->sample, attr_min, attr_max);
  if (verbose ())
    cout << attr_min << " " << attr_max << endl;
    
  return mv. getSD ();
}



void UniKernel::estimate_ (Real halfWindow_lo,
                           Real halfWindow_hi,
                           Real step)
{
  ASSERT (positive (halfWindow_lo));
  ASSERT (halfWindow_lo < halfWindow_hi);
  ASSERT (step > 1);
  
  
  Real halfWindow_best = NaN;   
  Prob uniform_prob_best = NaN;   
  Real logLikelihood_max = -INF;  
  halfWindow = halfWindow_lo;
  while (leReal (halfWindow, halfWindow_hi))
  {  
    halfWindow *= step;    
    const Real height = getHeight ();
    set_uniform_prob ();

    // Jack-knife
    Sum logLikelihood;
    for (const Point& point : points)
    {
      const Real pdfValue = pdf (point. value) - (1 - uniform_prob) * point. mult * height;
        // Correct for: checkPtr (analysis) -> sample. mult_sum -= p. mult ??
      ASSERT (pdfValue >= 0);
      logLikelihood. add (log (pdfValue) * point. mult);
    }
        
    Unverbose unv;
    if (verbose ())
      cout << "uniform_prob = " << uniform_prob << "  halfWindow = " << halfWindow << "  logLikelihood = " << logLikelihood. get () << endl;
    
    if (maximize (logLikelihood_max, logLikelihood. get ()))
    {
      halfWindow_best   = halfWindow;
      uniform_prob_best = uniform_prob;
    }
  }
  ASSERT (halfWindow_best);
  
  
  halfWindow   = halfWindow_best;
  uniform_prob = uniform_prob_best;
  if (verbose ())
    cout << "step = " << step << "  best: uniform_prob = " << uniform_prob << "  halfWindow = " << halfWindow << "  logLikelihood = " << logLikelihood_max << endl;
}



void UniKernel::set_uniform_prob ()
{
  ASSERT (positive (halfWindow));
  
#if 1  
  Real outlierMult = 0;
  FFOR (size_t, i, points. size ())
    if (   (i == 0                   || geReal (points [i].     value - points [i - 1]. value, halfWindow))
        && (i == points. size () - 1 || geReal (points [i + 1]. value - points [i].     value, halfWindow))
       )
      outlierMult += points [i]. mult;
  uniform_prob = outlierMult / checkPtr (analysis) -> sample. mult_sum;
#else
  for (;;)
  {
    Real uniform_mult = 0;
    for (const Point& point : points)
    {
    //const Real pdfValue = pdf (point. value);
      const Real pdfValue = pdf (point. value) - (1 - uniform_prob) * point. mult * height;
      ASSERT (pdfValue > 0);
      const Prob p = (uniform_prob * getUniformHeight ()) / pdfValue;
      uniform_mult += p * point. mult;
    }
    const Prob uniform_prob_new = uniform_mult / checkPtr (analysis) -> sample. mult_sum;
    ASSERT (isProb (uniform_prob_new));
    if (eqReal (uniform_prob, uniform_prob_new, 1e-5))  // PAR
      break;
    uniform_prob = uniform_prob_new;
  }    
  ASSERT (isProb (uniform_prob));
  {
    ONumber on (cout, 6, true);
    cout << "uniform_prob = " << uniform_prob << endl;
  }
#endif

  ASSERT (uniform_prob >= 0);
}



size_t UniKernel::findIndex (Point x) const
{
  ASSERT (! points. empty ());
  const size_t index = points. binSearch (x, false);  
    // = min {index : x <= points[index].value}
  if (index == NO_INDEX)
    if (points. empty ())
      return NO_INDEX;
    else
      return points. size () - 1;    
  else
    if (eqReal (x. value, points [index]. value))
      return index;
    else
      if (index == 0)
        return NO_INDEX;
      else
        return index - 1;
}



Real UniKernel::pdf_ (Real x) const
{
  if (halfWindow == 0)
    return eqReal (x, attr_min) ? INF : 0;
  
  const size_t left  = findIndex (Point (x - halfWindow));
  const size_t right = findIndex (Point (x + halfWindow));
  
  Real multSumLeft = 0;
  if (left == NO_INDEX)
  {
    if (right == NO_INDEX)
      return 0;
  }
  else
    multSumLeft = multSum [left];
      
  ASSERT (right != NO_INDEX);
  ASSERT (right < multSum. size ());
  
  const Real windowMultsum = multSum [right] - multSumLeft;
  ASSERT (windowMultsum >= 0);
  
  const Real uniformPdf = betweenEqualReal (x, attr_min, attr_max) ? getUniformHeight () : 0;
  return        uniform_prob  * uniformPdf
         + (1 - uniform_prob) * windowMultsum * getHeight ();  // = sum_i Kernel_i(x)
}



Prob UniKernel::cdf_ (Real x) const
{
  const size_t index = findIndex (Point (x));
  if (index == NO_INDEX)
    return 0;
  ASSERT (index < multSum. size ());
  return multSum [index] / multSum. back ();
}



Real UniKernel::getMean () const
{
	Real average = NaN;
	Real scatter = NaN;
  getAttr (). getAverageScatter (analysis->sample, average, scatter);
  return average;
}



Real UniKernel::getVar () const
{
	Real average = NaN;
	Real scatter = NaN;
  getAttr (). getAverageScatter (analysis->sample, average, scatter);
  return scatter;
}




// MultiDistribution

void MultiDistribution::qc () const
{
  if (! qc_on)
    return;
  Distribution::qc ();
    
  if (getParamSet ())
  {
    QC_ASSERT (variable. size () == getDim ());
    QC_ASSERT (x_field. size () == getDim ());
  }
}



Analysis1* MultiDistribution::createAnalysis (Dataset &ds)
{
  Space1<NumAttr1> space (ds, false);
  FFOR (size_t, i, getDim ())
    space << new RealAttr1 ("X" + toString (i + 1), ds);  
  const Sample sm (ds);
  auto* analysis_ = new An (sm, space);
  analysis = analysis_;
  
  return analysis_;
}




// MultiNormal

void MultiNormal::qc () const
{
  if (! qc_on)
    return;
  MultiDistribution::qc ();

  if (getParamSet ())
  {  	
	  QC_ASSERT (mu. size () == getDim ());
	  QC_ASSERT (mu. defined ());
	  
	  QC_ASSERT (sigmaExact. defined ());
	  QC_ASSERT (sigmaExact. isSimilarity ());
	  QC_ASSERT (sigmaExact. rowsSize (false) == getDim ());
	  QC_ASSERT (sigmaExact. psd);

	  QC_ASSERT (sigmaInflated. defined ());
	  QC_ASSERT (sigmaInflated. isSimilarity ());
	  QC_ASSERT (sigmaInflated. rowsSize (false) == getDim ());
	  QC_ASSERT (sigmaInflated. psd);

	  QC_ASSERT (variance_min. size () == getDim ());
    QC_ASSERT (variance_min. defined ());
	  QC_ASSERT (variance_min. min () >= 0);
	  
	  QC_IMPLY (variance_min. min () == 0, sigmaExact. maxAbsDiff (false, sigmaInflated, false) == 0);

    if (coeff != INF)
    { 
  	  QC_ASSERT (sigmaInv. defined ());
  	  QC_ASSERT (sigmaInv. isSimilarity ());
  	  QC_ASSERT (sigmaInv. rowsSize (false) == getDim ());
  	  QC_ASSERT (sigmaInv. psd);
  	}

	  QC_ASSERT (zs. size () == getDim ());
	}

	for (const Normal& n : zs)
	  n. qc ();
}



void MultiNormal::saveText (ostream& os) const
{
  MultiDistribution::saveText (os);
  os << endl;    	

  os << "Mean:" << endl;    	
  mu. saveText (os);

  os << "VC(exact):" << endl;    	
  sigmaExact. saveText (os);

  if (variance_min. min ())
  {  
    os << "VC(inflated):" << endl;    	
    sigmaInflated. saveText (os);

    os << "variance_min:" << endl;
    variance_min. saveText (os);  
  }
  
  if (verbose ())
  {  
    os << "VC-1:" << endl;    	
    sigmaInv. saveText (os);
    
    os << "coeff: " << coeff << endl;
  }
}



void MultiNormal::setDim (size_t dim)
{
	ASSERT (dim);
	
	MultiDistribution::setDim (dim);

  mu. resize (dim);
	mu. putAll (0);

  sigmaExact. resize (false, dim, dim);
	sigmaExact. putAll (0);
	
  sigmaInflated. resize (false, dim, dim);
	sigmaInflated. putAll (0);
	
  if (! variance_min. size ())
  {
    variance_min. resize (dim);
    variance_min. putAll (0);
  }

  sigmaInv. resize (false, dim, dim);
	sigmaInv. putAll (NaN);
	
  zs. resize (dim);
	for (Normal &n : zs)
	  n. setParam (0, 1);
}



void MultiNormal::setParamFunc ()
{ 
  Unverbose unv;
  
  // sigma
	if (verbose ())
	{
		cout << endl;
		cout << "variance_min:" << endl;
		variance_min. saveText (cout);
		cout << "Mu:" << endl;
		mu. saveText (cout);
		cout << "VC:" << endl;
		sigmaExact. saveText (cout);
	}
  const bool inflated = inflateSigma ();
	if (verbose ())
	{
		cout << endl;
		cout << "Inflated = " << inflated << endl;
		cout << "VC:" << endl;
		sigmaInflated. saveText (cout);
	}
  	
  coeff = INF;

	cholesky = sigmaInflated. getCholesky (false);
	if (! cholesky. defined ())
	  return;

  sigmaInv = sigmaInflated; 
	const Determinant det (sigmaInv. inverse ());
	if (   det. get () <= 0
	    || nullReal (det. getAverage () /*, 1e-5*/)  // PAR  
	   )  
	{
	#if 0
	  cout << "det = " << det. get () << "  ave = " << det. getAverage () << "  abs = " << sigmaInv. maxAbs () << endl;
	  ERROR;
	#endif
	  return;
	}
	
	coeff = 0.5 * ((Real) getDim () * Normal::coeff + det. n);		

	if (verbose ())
	{
		cout << endl;
		cout << "VC-1:" << endl;
		sigmaInv. saveText (cout);
		cout << endl;
		cout << "Cholesky:" << endl;
	  cholesky. saveText (cout);
	}
}



void MultiNormal::setSeed (ulong seed) const
{
	MultiDistribution::setSeed (seed);
	FFOR (ulong, i, zs. size ())
	  zs [i]. setSeed (seed + i + 1);
}



void MultiNormal::randVariable () const
{
  MVector zVec (zs. size ());
	FFOR (size_t, i, zs. size ())
	  zVec [i] = zs [i]. rand ();
	x_field. multiply (false, cholesky, false, zVec, false);
	x_field. add (false, mu, false, 1);
	
  FFOR (size_t, i, variable. size ())	
	  variable [i] = x_field [i];
}



void MultiNormal::estimate ()
{
  ASSERT (analysis);
    
  setDim (analysis->space. size ());
  
  mu. putAll (0);
  sigmaExact. putAll (0);

  MVector muCount (mu);
  Matrix sigmaCount (sigmaExact);
  {
    Progress prog (analysis->sample. nEffective, analysis->size (1, 2) > 1e7);   // PAR
    for (Iterator it (analysis->sample); it ();)  
    {
      prog ();
    	data2variable (*it);
    	const Real m = it. mult;
    	FFOR (size_t, row, mu. size ())
    	{
    	  const Real a = variable [row];
    	  if (isNan (a))
    	    continue;
    	  mu      [row] += a * m;
    	  muCount [row] += m;
      	FFOR (size_t, col, row + 1)
      	{
  	  	  const Real b = variable [col];
  	  	  if (isNan (b))
  	  	    continue;
      	  sigmaExact. putInc (false, row, col, a * b * m);
      	  sigmaCount. putInc (false, row, col, m);
      	}
      }
    }
  }
  
	FFOR (size_t, row, sigmaExact. rowsSize (false))
	{
	  {
  	  const Real n = muCount [row];
  	  if (nullReal (n))
  	    mu [row] = 0;
  	  else
  	    mu [row] /= n;
  	}
  	FFOR (size_t, col, row + 1)
  	{
  	  const Real n = sigmaCount. get (false, row, col);
  	  if (nullReal (n))
  	    sigmaExact. put (false, row, col, 0);
  	  else
  	    sigmaExact. putProd (false, row, col, 1 / n);
  	}
  }
  
  sigmaExact. lower2upper (false);
  
  // Centering
  Matrix mu2 (getDim ());
  mu2. multiply (false, mu, false, mu, true);
  sigmaExact. add (false, mu2, false, -1);
  
  sigmaExact. psd = true;
  
  setParam ();
  if (verbose ())
    if (isNan (coeff))
      cout << "det = 0" << endl;
}



void MultiNormal::getSigmaRaw (Matrix &sigmaRaw) const
{
  ASSERT (sigmaRaw. equalLen (false, sigmaExact, false));
  
  sigmaRaw. multiply (     false
                     , mu, false
                     , mu, true
                     );       
  sigmaRaw. add (false, sigmaExact, false, 1);
}



bool MultiNormal::inflateSigma ()
{
  sigmaInflated = sigmaExact;
  
  if (nullReal (variance_min. maxAbs ()))
	  return false;


  if (verbose ())
    cout << "inflateSigma" << endl;
  ASSERT (sigmaExact. psd);
	Eigens eigens (sigmaExact, sigmaExact. rowsSize (false), 1, 0, 1e-4 * sqrt (variance_min. min ()), 2000);  // PAR
	{
	  Unverbose unv;
	  if (verbose ())
	    eigens. saveText (cout);
	}
	if (isNan (eigens. explainedVarianceFrac))
	{
  	FFOR (size_t, i, sigmaExact. rowsSize (false))
  	  sigmaInflated. putDiag (i, variance_min [i]);
  	return true;
	}
  eigens. makeSquare ();		
	eigens. qc ();

	bool changed = false;
	ASSERT (sigmaInflated. rowsSize (false) == eigens. getDim ());
  FFOR (size_t, i, sigmaInflated. rowsSize (false))
    if (maximize (eigens. values [i], variance_min [i]))
    {
    	changed = true;
    	if (verbose ())
    	  cout << "eigens.values[" << i << "] = " << variance_min [i] << endl;
    }
  if (! changed)
  	return false;
	
	eigens. restore (sigmaInflated);

	
	return true;
}




// Mixture

void Mixture::Component::qc () const
{
  if (! qc_on)
    return;
	QC_ASSERT (distr. get ());
	distr->qc ();
	QC_ASSERT (isProb (prob));
}



Real Mixture::Component::getMult () const
{
	Real s = 0;
  for (Iterator it (distr->getAnalysisCheck () -> sample); it ();)  
    s += objProb [*it] * it. mult;
  return s;
}



Real Mixture::Component::getDeltaness () const
{
	Real s = 0;
  for (Iterator it (distr->getAnalysisCheck () -> sample); it ();)  
    maximize (s, objProb [*it] * it. mult);
    
  return getMult () - s;
}



void Mixture::Component::estimate ()
{
  Sample& sample = var_cast (distr->getAnalysisCheck ()) -> sample;
  const Component2Sample c2s (sample, *this);  
  if (verbose ())
    cout << "multSum = " << sample. mult_sum << endl;
//if (sample. mult_sum < 0.001)  // PAR
  //return false;
	distr->estimate (); 
	ASSERT (distr->getParamSet ());
//return true;
}



Real Mixture::Component::getFitness_entropy () const
{ 
  const Component2Sample c2s (var_cast (distr->getAnalysisCheck ()) -> sample, * var_cast (this)); 
	return distr->getFitness_entropy (); 
}



void Mixture::Component::merge (const Component* comp)
{
  ASSERT (comp);
  ASSERT (comp != this);  
  
  prob += comp->prob;
  makeProb (prob);
  
  ASSERT (objProb. size () == comp->objProb. size ());
  FFOR (size_t, i, objProb. size ())
  {
    objProb [i] += comp->objProb [i];
    makeProb (objProb [i]);
  }
}



void Mixture::qc () const
{
  if (! qc_on)
    return;
	Distribution::qc ();

	if (! getParamSet ())
  	return;

	Prob s = 0.0;
	const Distribution* distr0 = components [0] -> distr. get ();
	const bool discrete = distr0->asDiscreteDistribution ();
	const size_t objs = getAnalysis () ? getAnalysis () -> sample. mult. size () : 0;
	for (const Component* comp : components)
	{
		comp->qc ();
	  QC_ASSERT (comp->distr->getAnalysis () == getAnalysis ());
		QC_ASSERT ((bool) comp->distr->asDiscreteDistribution () == discrete);
		QC_ASSERT (comp->distr->getDim () == getDim ());
	  const UniDistribution* u = comp->distr->asUniDistribution ();
    if (const UniDistribution* u0 = distr0->asUniDistribution ())
		{
			QC_ASSERT (u);
			QC_ASSERT (eqReal (u0->loBound, u->loBound));
			QC_ASSERT (eqReal (u0->hiBound, u->hiBound));
		}
		else
		  { QC_ASSERT (! u); }
		QC_ASSERT (comp->objProb. size () == objs);
		QC_IMPLY (comp->distr->getAnalysis (), objs == comp->distr->getAnalysis () -> sample. mult. size ());
	  s += comp->prob;      
	}
	QC_ASSERT (eqReal (s, 1.0));
	
  FOR (size_t, objNum, objs)  
  {
    Prob objS = 0.0;
  	for (const Component* comp : components)
  	  objS += comp->objProb [objNum];
  	QC_IMPLY (! isNan (objS), eqReal (objS, 1.0));
  }
	
	cat. qc ();
}



void Mixture::saveText (ostream& os) const
{
	FFOR (size_t, i, components. size ())
	{
		os << i + 1 << ": ";
		const Component* comp = components [i];
		comp->saveText (os);
		if (! comp->distr->asUniDistribution ())
		  os << endl;
	}
}



Mixture::Component* Mixture::addComponent (Distribution* distr,
                                           Prob prob)
{ 
  ASSERT (distr)
//ASSERT (distr->analysis);
  
  auto c = new Component (distr, prob);
  if (const Analysis1* an = distr->getAnalysis ())
    c->objProb. resize (an->sample. ds->objs. size (), NaN); 
  components << c; 
  
  return c;
}



void Mixture::setParamFunc ()
{ 
	balanceProb (); 
	
	cat. probs. clear ();  
	cat. probs. reserve (components. size ());
	for (const Component* comp : components)
	  cat. probs << comp->prob;
	cat. setParam ();
}



void Mixture::setSeed (ulong seed) const
{
  Distribution::setSeed (seed);
    
	cat. setSeed (seed + 1);
	FFOR (size_t, i, components. size ())
	  components [i] -> distr. get () -> setSeed (seed + 2 + i);
}



size_t Mixture::paramCount () const
{
	size_t n = components. size () - 1;  // Component::prob
	for (const Component* comp : components)
	  n += comp->distr->paramCount ();
	  
	return n;
}



bool Mixture::getParamSet () const 
{ 
	if (components. empty ())
		return false;
	for (const Component* comp : components)
  {
    ASSERT (comp);
  	if (isNan (comp->prob))
  		return false;
	  if (! comp->distr->getParamSet ())
	  	return false;
	}
	return true; 
}



bool Mixture::similar (const Distribution &distr,
                       Real delta) const
{ 
	if (const Mixture* mixt = distr. asMixture ())
	{
		if (components. size () != mixt->components. size ())
			return false;
		FFOR (size_t, i, components. size ())
		  if (! components [i] -> distr->similar (* mixt->components [i] -> distr, delta))
		  	return false;
		return true;
	}
		
  return false;
}



Real Mixture::pdfVariable () const
{
	Real r = 0;
	for (const Component* comp : components)
	  r += comp->pdfProb ();
	return r;
}



Real Mixture::logPdfVariable () const
{
#if 0
	const Real r = log (pdfVariable (x));
	if (r > 0 || finite (r))
	  return r;
#endif
	  
	SumLn s;
	for (const Component* comp : components)
	  s. addExp (comp->logPdfProb ());
	return s. getLn ();
}



#if 0
void Mixture::getMeanVec (MVector &mean) const 
{
	ASSERT (mean. size () == getDim ());
	
	mean. putAll (0);
	MVector x (getDim ());
	for (const Component* comp : components)
	{
		comp->distr->getMeanVec (x);
		mean. add (false, x, false, comp->prob);
	}
}



void Mixture::getVC (Matrix &vc) const
{
	ASSERT (vc. rowsSize (false) == getDim ());
	ASSERT (vc. rowsSize (true)  == getDim ());

#if 1
  NOT_IMPLEMENTED;
#else	
	vc. putAll (0);
	Matrix vcComp (getDim ());
	for (const Component* comp : components)
	{
		comp->distr->getVC (vcComp);
		vc. add (false, vcComp, false, sqr (comp->prob));
	}
#endif
}
#endif



namespace 
{
  bool componentComp (const Mixture::Component* a,
                      const Mixture::Component* b)
  {
    ASSERT (a);
    ASSERT (b);
    return a->distr->getSortingValue () < b->distr->getSortingValue ();
  }
}



void Mixture::estimate ()
{
  ASSERT (! components. empty ());
  
  const Analysis1* analysis = getAnalysisCheck ();
  
  bool objProbStart = getParamSet ();
//cout << "objProbStart = " << objProbStart << endl;  

#ifndef NDEBUG	
	if (! objProbStart)
    for (const Component* c : components)
    {
    	ASSERT (c->objProb. size () == analysis->sample. mult. size ());
      for (Iterator it (analysis->sample); it ();)      	
    	  ASSERT (isProb (c->objProb [*it]));
    }
#endif
  	

  Real entropyEst_best = INF;
  do
  {
    if (verbose ())
      cout << endl << "entropyEst_best = " << entropyEst_best << endl;
      
	  // Component::objProb
	  if (objProbStart)
	  {
		  MVector vec (components. size ());
		  for (Iterator it (analysis->sample); it ();)  
		    if (components. size () == 1)
		      var_cast (components [0]) -> objProb [*it] = 1;
		    else
  		  {
  		  	data2variable (*it);
  		  //ASSERT (! variable. contains (NaN));
  		  	// Bayes' theorem
  				FFOR (size_t, i, components. size ())
  				  vec [i] = components [i] -> logPdfProb ();				
  				vec. expBalanceLogRow (true, 0, 1);
  				vec. expRow (true, 0);
  				ASSERT (eqReal (vec. sumRow (true, 0), 1));
  				FFOR (size_t, i, components. size ())
  				{
  					const Prob newObjProb = vec [i];
  					ASSERT (isProb (newObjProb));
  					var_cast (components [i]) -> objProb [*it] = newObjProb;
  			  }
  		  }  
	  }
	  else
	  	objProbStart = true;
	
	  // Component::distr->Parameters
	  FOR_REV (size_t, i, components. size ())  // --> for ??
	  {
	    if (verbose ())
	      cout << "Estimating component " << (i + 1) << endl;
	  #if 1
	    var_cast (components [i]) -> estimate ();
	    if (verbose ())
	      components [i] -> print (cout);
	  #else
	    if (   ! var_cast (components [i]) -> estimate ()
	        && components. size () >= 2
	       )
	    {
	    	components. erasePtr (i);
  	    if (verbose ())
  	      cout << "Component is erased" << endl;
	    }
	  #endif
	  }
	
	  // Component::prob
	  {
		  MVector vec (components. size ());
			FFOR (size_t, i, components. size ())
		    vec [i] = components [i] -> getMult ();
	  	vec. balanceRow (true, 0, 1);
	  	FFOR (size_t, i, components. size ())
	  	{
	  	  var_cast (components [i]) -> prob = toProb (vec [i]);
	  	  if (verbose ())
	  	  {
	  	    ONumber on (cout, 6, true);
	  	    cout << "P(component " << (i + 1) << ") = " << vec [i] << endl;
	  	  }
	  	}
		}
  }
  while (minimizeEntropy (entropyEst_best, 1e-4));  // PAR  
  

  components. sort (componentComp);
  
  
  setParam ();
}



#if 0
Prob Mixture::getFitness_entropy_min () const
{
	Prob p = 1;
	for (const Component* comp : components)
	  minimize (p, comp->getFitness_entropy ());
  return p;
}
#endif



Prob Mixture::getConfusion () const
{
  const Analysis1* analysis = getAnalysisCheck ();

	Real s = 0;
  for (Iterator it (analysis->sample); it ();)  
  {
  	Prob p_max = 0;
  	for (const Component* comp : components)
  	  maximize (p_max, comp->objProb [*it]);
  	s += p_max * it. mult;
  }
  
  const Prob p = 1 - s / analysis->sample. mult_sum;
  ASSERT (isProb (p));
  
  return p;
}



Prob Mixture::getOverlap (size_t compNum1,
                          size_t compNum2) const
{
  const Analysis1* analysis = getAnalysisCheck ();

  const Component* comp1 = components [compNum1];
  const Component* comp2 = components [compNum2];
  Real confused = 0;
  Real total = 0;
  for (Iterator it (analysis->sample); it ();)  
  {
  	data2variable (*it);
    const Real r1 = comp1->pdfProb ();
    const Real r2 = comp2->pdfProb ();
    total += (r1 + r2) * it. mult;
    confused += min (r1, r2) * it. mult;
  }
  
  const Prob p = confused / total;
  ASSERT (isProb (p));
  ASSERT (leReal (p, 0.5));
  
  return p;
}



void Mixture::balanceProb ()
{
	Real s = 0;
	for (const Component* comp : components)
	  s += comp->prob;
	ASSERT (s > 0);
	for (const Component* comp : components)
	  var_cast (comp) -> prob /= s;
}



void Mixture::mergeComponents ()
{
  FFOR (size_t, i, components. size ())
  {
	  const Component* comp = components [i];
  	FOR_REV_END (size_t, j, i + 1, components. size ())
  	{
  	  const Component* other = components [j];
	    if (comp->distr->similar (* other->distr, 1e-3))  // PAR  
	    {
        var_cast (comp) -> merge (other);
  			components. erasePtr (j);
  	  }
  	}
  }
  
  setParam ();
}



bool Mixture::deleteComponent (size_t num)
{
  const Analysis1* analysis = getAnalysisCheck ();

	const Component* comp_del = components [num];
	
	const Prob p_del = comp_del->prob;

	if (eqReal (p_del, 1))
	  return false;
  for (Iterator it (analysis->sample); it ();)  
  	if (eqReal (comp_del->objProb [*it], 1))
  	  return false;

  for (const Component* c : components)
    if (c != comp_del)
    {
      Prob& p = var_cast (c) -> prob;
      p /= 1 - p_del;
      makeProb (p);
    }
  for (Iterator it (analysis->sample); it ();)  
  {
    const Prob p_obj_del = comp_del->objProb [*it];
  	ASSERT (lessReal (p_obj_del, 1));
	  for (const Component* c : components)
      if (c != comp_del)
  	  {
  	    Prob& p = var_cast (c) -> objProb [*it];
  	    p /= 1 - p_obj_del;
  	    makeProb (p);
  	  }
  }
	
	components. erasePtr (num);	  
  setParam ();

	return true;
}




// PrinComp

void PrinComp::qc () const
{ 
  if (! qc_on)
    return;
  P::qc ();
    
  mn. qc ();
  QC_ASSERT (space. size () == mn. getDim ());
  QC_IMPLY (mn. getAnalysis (), mn. getAnalysis () == this);
  
	eigens. qc ();
	QC_ASSERT (! isNan (eigens. explainedVarianceFrac));
	QC_ASSERT (eigens. getInitSize () == mn. getDim ());
	QC_ASSERT (eigens. getDim () <= mn. getDim ());
}



void PrinComp::project (const PrinComp::P::Value &x,
                        MVector &projection) const
{	
  ASSERT (x. size () == space. size ());

  MVector x_centered (x);
  x_centered. addRow (         true, 0
                     , mn. mu, true, 0
                     , -1.0
                     );
  projection. multiply ( false
                       , eigens. basis, true
                       , x_centered, false
                       ); 
}



Space1<RealAttr1> PrinComp::createSpace (const string &attrPrefix,
                                         Dataset &ds) const
{
	ASSERT (! attrPrefix. empty ());
	ASSERT (& ds == sample. ds);
	
	Space1<RealAttr1> sp (ds, false);
	FFOR (size_t, i, getOutDim ())
	  sp << new RealAttr1 (attrPrefix + toString (i + 1), ds);
	
  MVector projection (getOutDim ());  
  FFOR (size_t, i, ds. objs. size ())
  {
    project (i, projection);
  	FFOR (size_t, col, projection. size ())
  	  (* var_cast (sp [col])) [i] = projection [col];
  }
	
	return sp;
}



Dataset PrinComp::createAttrMds (const string &attrPrefix,
                                 Vector<Real> &quality) const
{
  ASSERT (quality. empty ());
  
	Dataset ds;
	FFOR (size_t, attrNum, mn. mu. size ())
	  ds. appendObj (space [attrNum] -> name);
	FFOR (size_t, eigenNum, getOutDim ())
  {
	  auto* attr = new RealAttr1 (attrPrefix + toString (eigenNum + 1), ds);
	  FFOR (size_t, attrNum, mn. mu. size ())
	    (*attr) [attrNum] = getAttrMds (attrNum, eigenNum);
	  quality << eigens. explainedFrac (eigenNum);
	}
	
	return ds;
}



Real PrinComp::getChi2 (size_t objNum) const
{
  MVector projection (getOutDim ());  
  project (objNum, projection);
  Real chi2 = 0;
  FFOR (size_t, i, getOutDim ())
  {
    const Real coeff = eigens. values [i];
    ASSERT (positive (coeff));
    chi2 += sqr (projection [i]) / coeff;
  }
  return chi2;
}




// Clustering

Clustering::Clustering (const Sample &sample_arg,
                        const Space1<NumAttr1> &space_arg,
	                      size_t clusters_max,
	                      Real sd_min,
	                      bool sd_min_is_relative,
	                      Real entropyDimensionPrecision)
: P (sample_arg, space_arg)
, variance_min (space. size ())
{
  ASSERT (clusters_max >= 2);
  ASSERT (sd_min >= 0);
  ASSERT (entropyDimensionPrecision >= 0);
  
  
  // variance_min
  {
    const Real var_rel_min = sqr (sd_min);
    FFOR (size_t, i, space. size ())
    {
      Real var = 1;
      if (sd_min_is_relative)
      {
        const UniVariate<NumAttr1> an (sample, * space [i]);
        Normal normal;
        normal. analysis = & an;
        normal. estimate ();
        normal. qc ();
        var = normal. getVar ();
      }
      variance_min [i] = var * var_rel_min;
    }
  }


  ASSERT (mixt. components. empty ());
  mixt. qc ();
  
  
  Unverbose unv;
  
  
  // Global cluster
  {
  	if (verbose ())
  	  cout << "Global Cluster" << endl; 
  	Mixture::Component* comp = addCluster (1);
	  comp->distr->estimate ();
	  mixt. setParam ();
	}
  mixt. estimate (); 
	if (verbose ())
    mixt. qc ();
  if (! mixt. getParamSet ())  
    return;


  {
    Progress prog;
	  Real entropy_best = mixt. getEntropy_est ();
  	FOR (uint, iter, 10000)  // PAR
	  {
	  	if (verbose ())
	  	{
	  	  cout << endl << "# Clusters = " << mixt. components. size () << "  Entropy_best = " << entropy_best << endl; 
	  	  Unverbose unv1;
		  	if (verbose ())
		  	{
  	      mixt. print (cout);  
		  	  cout << endl;  
  	    }
	  	}

	  	Mixture workBest;  
	    ASSERT (! workBest. getParamSet ());

	  	// Splitting a cluster
	  	if (mixt. components. size () < clusters_max)
	  	{
  	  	FFOR (size_t, i, mixt. components. size ())
  	  	{
  	  		Clustering work (*this);
  	  		work. setAnalysis ();  // --> Clustering copy constructor ??
  	  	//work. print (cout);  
  	  	  work. splitCluster (i);
  	  	  work. mixt. estimate ();
  	  	  if (verbose ())
  	  	  {
  	  	    cout << "Split cluster " << i + 1 << ":"
  	  	             << " entropy = " << work. mixt. getEntropy_est ()  
  	  	             << "  confusion = " << work. mixt. getConfusion ()
  	  	             << "  paramCout = " << work. mixt. paramCount ()
  	  	             << "  mult_sum = " << sample. mult_sum
  	  	             << "  entropy_best = " << entropy_best
  	  	             << endl;
  	  	    work. mixt. print (cout);
  	  	  }
  	  	  if (   work. mixt. getParamSet ()
  	  	      && work. mixt. minimizeEntropy (entropy_best, entropyDimensionPrecision) 
  	  	     )
  	  	  	workBest = move (work. mixt);
  	  	}
  	  	if (verbose ())
  	      workBest. qc ();
  	    prog (toString (mixt. components. size ()));
  	  }

	    // Deleting a cluster
	    if (mixt. components. size () >= 2)
	    {
  	  	FFOR (size_t, i, mixt. components. size ())
  	  	{
  	  		Clustering work (*this);
  	  		work. setAnalysis ();
  	  	//work. print (cout);  
  	  	  if (! work. deleteCluster (i))
  	  	    continue;
  	  	  work. mixt. estimate ();
  	  	  if (verbose ())
  	  	  {
  	  	    cout << "Deleted cluster " << i + 1 << ":"
  	  	             << " entropy = " << work. mixt. getEntropy_est ()  
  	  	             << "  confusion = " << work. mixt. getConfusion ()
  	  	             << "  paramCout = " << work. mixt. paramCount ()
  	  	             << "  mult_sum = " << sample. mult_sum
  	  	             << "  entropy_best = " << entropy_best
  	  	             << endl;
  	  	    work. mixt. print (cout);
  	  	  }
  	  	  if (   work. mixt. getParamSet ()
  	  	      && work. mixt. minimizeEntropy (entropy_best, entropyDimensionPrecision) 
  	  	     )
  	  	  	workBest = /*move*/ (work. mixt);
  	  	}
  	  	if (verbose ())
	        workBest. qc ();
  	    prog (toString (mixt. components. size ()));
	    }

	    if (! workBest. getParamSet ()) 
	      break;
	    mixt = /*move*/ (workBest);
	    setAnalysis ();
	  }
	}


  ASSERT (mixt. getParamSet ());
	mixt. mergeComponents ();
	mixt. qc ();
}



Clustering::Clustering (const Clustering &clustering,
                        const vector<bool> &toMerge)
: P (clustering)
, variance_min (clustering. variance_min)
{
  ASSERT (clustering. getOutDim () == toMerge. size ());
  mixt = clustering. mixt;
  
  size_t compNum1 = NO_INDEX;
  FFOR (size_t, i, mixt. components. size ())
    if (toMerge [i])
    {
      compNum1 = i;
      break;
    }
    
  if (compNum1 != NO_INDEX)    
    FOR_REV (size_t, compNum2, mixt. components. size ())
    {
      if (compNum2 == compNum1)
        break;
      else
        if (toMerge [compNum2])
          merge (compNum1, compNum2);        
    }
        
  mixt. setParam ();  
  setAnalysis ();
}
  


void Clustering::qc () const
{
  if (! qc_on)
    return;
  P::qc ();
   
  mixt. qc ();
  for (const Mixture::Component* c : mixt. components)
  	QC_ASSERT (c->distr->asMultiNormal ());

  QC_ASSERT (getOutDim ()); 
  QC_ASSERT (variance_min. min () >= 0);
}



void Clustering::saveText (ostream &os) const
{ 
  if (mixt. getDim () == 1)
  {
    {
      TabDel td;
      td << "Cluster" << "Mean" << "Var" << "P";
      os << td. str () << endl;
    }
    FFOR (size_t, i, mixt. components. size ())
    {
      TabDel td;
      const MultiNormal* mn = getMultiNormal (i);
      ASSERT (mn);
      td << (i + 1) << mn->mu [0] << mn->sigmaInflated. get (false, 0, 0) << mixt. components [i] -> prob;      
      os << td. str () << endl;
    }
  }
  else
    mixt. saveText (os); 
}



void Clustering::setAnalysis ()
{ 
  FFOR (size_t, i, getOutDim ())
    var_cast (getMultiNormal (i)) -> analysis = this;
}



void Clustering::splitCluster (size_t num)
{
  if (verbose ())
    cout << "Splitting cluster " << num + 1 << endl;

	Mixture::Component* oldComp = var_cast (mixt. components [num]);
	oldComp->prob /= 2;  // ??
	const MultiNormal* oldMn = oldComp->distr->asMultiNormal ();
	ASSERT (oldMn);
	ASSERT (oldMn->getParamSet ());
	
	Mixture::Component* newComp = addCluster (oldComp->prob);

  // oldComp->objProb[], newComp->objProb[]
	const PrinComp pc (sample, space, *oldMn, 1, 1.0, 0.0, 0.0); 
	ASSERT (pc. getOutDim () <= 1);
  MVector projection (1);
  for (Iterator it (sample); it ();)  
  {
  	if (pc. getOutDim () == 1)
  	  pc. project (*it, projection);
  	else
  	  projection [0] = 0;
  	if (projection [0] < 0)  
  		newComp->objProb [*it] = 0;
  	else
		{
  		newComp->objProb [*it] = oldComp->objProb [*it];
  		oldComp->objProb [*it] = 0;
		}
  }
  
  ASSERT (! mixt. getParamSet ());

  if (verbose ())
    cout << "Cluster is split" << endl;
}



Space1<ProbAttr1> Clustering::createSpace (Dataset &ds) const
{
  ASSERT (& ds == & space. ds);
  
	Space1<ProbAttr1> sp (ds, false);
	FFOR (size_t, i, getOutDim ())
	  sp << new ProbAttr1 ( "Cluster" + toString (i + 1), ds, 2/*PAR*/);
	
  for (Iterator it (sample); it ();)  
  	FFOR (size_t, col, getOutDim ())
  	  (* var_cast (sp [col])) [*it] = mixt. components [col] -> objProb [*it];
	
	return sp;
}



NominAttr1* Clustering::createNominAttr (const string &attrName,
	                                       Prob prob_min,
	                                       Dataset &ds) const
{
	ASSERT (isProb (prob_min));
//ASSERT (prob_min >= 0.5);
	
	auto* attr = new NominAttr1 (attrName, ds);
 	FFOR (size_t, col, getOutDim ())
 	  EXEC_ASSERT (attr->category2index ("C" + toString (col + 1)) == col);
  ASSERT (! attr->categories. empty ());

  for (Iterator it (sample); it ();)  
  {
  	size_t col_best = NO_INDEX;
  	Prob p_max = 0;
  	FFOR (size_t, col, getOutDim ())
  	  if (maximize (p_max, mixt. components [col] -> objProb [*it]))
  	  	col_best = col;
  	ASSERT (col_best != NO_INDEX);
  	if (geReal (p_max, prob_min))
 	  	(*attr) [*it] = col_best;
  }
  
  return attr;
}



ProbAttr1* Clustering::createProbAttr (const string &attrName,
                                       streamsize decimals,
                                       Dataset &ds) const
{
	auto* attr = new ProbAttr1 (attrName, ds, decimals);  

  for (Iterator it (sample); it ();)  
  {
  	Prob p_max = 0;
  	FFOR (size_t, col, getOutDim ())
  	  maximize (p_max, mixt. components [col] -> objProb [*it]);
  	(*attr) [*it] = p_max;
  }
  	  
  return attr;
}



namespace 
{

struct CategoryCluster : DisjointCluster
{
  size_t category;
  explicit CategoryCluster (size_t category_arg)
    : category (category_arg)
    {}
};  
  
}



bool Clustering::mergeClose (NominAttr1 &nominAttr,
                             ProbAttr1 &probAttr,
                             Prob confused_max) const
{
  ASSERT (& nominAttr. ds == sample. ds);
  ASSERT (& probAttr.  ds == sample. ds);

  VectorOwn<CategoryCluster> clusters;  clusters. reserve (mixt. components. size ());
  FFOR (size_t, i, mixt. components. size ())
    clusters << new CategoryCluster (i);
  ASSERT (clusters. size () == mixt. components. size ());
  FFOR (size_t, i, mixt. components. size ())
  	FOR (size_t, j, i)
  	{
  	  const Prob confused = mixt. getOverlap (i, j);
  	  if (verbose ())
  	    cout << "Cluster " << i + 1 << " vs. cluster " << j + 1 << ": confused = " << confused << endl;
 	    if (geReal (confused, confused_max))
 	      var_cast (clusters [i]) -> merge (* var_cast (clusters [j]));
    }

  unordered_map<const DisjointCluster*, Vector<size_t>> cluster2categories;  cluster2categories. rehash (clusters. size ());
  FFOR (size_t, i, clusters. size ())
    cluster2categories [var_cast (clusters [i]) -> getDisjointCluster ()] << i;
  ASSERT (! cluster2categories. empty ());

  if (cluster2categories. size () == 1)
    return false;

  bool merged = false;    
  for (const auto& clusterIt : cluster2categories)
  {
    const Vector<size_t>& categories = clusterIt. second;
    ASSERT (! categories. empty ());
    if (categories. size () == 1)
      continue;
    merged = true;
    const size_t mainCategory = categories [0];
    for (Iterator it (sample); it ();)  
    {
      if (! categories. contains (nominAttr [*it]))
        continue;
      ASSERT (nominAttr [*it] >= mainCategory);
      nominAttr [*it] = mainCategory;
      Prob p = 0;
      FFOR (size_t, i, categories. size ())
        p += mixt. components [categories [i]] -> objProb [*it];
      makeProb (p);
      probAttr [*it] = p;
    }
  }

  nominAttr. deleteEmptyCategories (sample);
  nominAttr. renameCategories ();
      
  nominAttr. qc ();
  probAttr. qc ();
  
  return merged;
}



void Clustering::merge (size_t compNum1,
                        size_t compNum2)
{
  ASSERT (compNum1 < getOutDim ());
  ASSERT (compNum2 < getOutDim ());  
  ASSERT (compNum1 != compNum2);
    
  Mixture::Component* comp1 = var_cast (mixt. components [compNum1]);
  Mixture::Component* comp2 = var_cast (mixt. components [compNum2]);
  ASSERT (comp1 != comp2);

  MultiNormal* mn1 = var_cast (comp1->distr->asMultiNormal ());
  MultiNormal* mn2 = var_cast (comp2->distr->asMultiNormal ());
  ASSERT (mn1);
  ASSERT (mn2);
  ASSERT (mn1 != mn2);
  
  const Prob p1 = comp1->prob;
  const Prob p2 = comp2->prob;
  
  MVector mu (mn1->mu. size (), 0);
//mu. putAll (0);
  mu. add (false, mn1->mu, false, p1);
  mu. add (false, mn2->mu, false, p2);
  mu. putProdAll (1 / (p1 + p2));
  
  Matrix sigma (false, mn1->sigmaExact, false, 0);
  {
  //sigma. putAll (0);
    Matrix raw (false, mn1->sigmaExact, false);
    mn1->getSigmaRaw (raw);
    sigma. add (false, raw, false, p1);
    mn2->getSigmaRaw (raw);
    sigma. add (false, raw, false, p2);
    sigma. putProdAll (1 / (p1 + p2));
    raw. multiply (     false
                  , mu, false
                  , mu, true
                  );       
    ASSERT (sigma. psd);
    sigma. add (false, raw, false, -1);  
    sigma. psd = true;
  }
  
  mn1->setDim (mn1->getDim ());
  mn1->mu = mu;
  mn1->sigmaExact = sigma;
  mn1->setParam ();
  
  comp1->merge (comp2);
  
	mixt. components. erasePtr (compNum2);	  
  mixt. setParam ();
}



void Clustering::processSubclusters (const string &clusterAttrName,
                                     bool merge_close,
                                     const VectorPtr<Attr> &attrs_orig,
                                     JsonArray* jClusts,
                                     Prob attr_pvalue,
                                     const string &outDir,
                                     const string &outGenericDm,
                                     Dataset &ds,
                                     Space1<Attr1> &sp1) const
{
  ASSERT (! outGenericDm. empty ());
  ASSERT (& ds == sample. ds);
  ASSERT (& ds == & sp1. ds);
  

  NominAttr1* clustNominAttr = createNominAttr (clusterAttrName, 0, ds);  
  ProbAttr1* probAttr = createProbAttr (clusterAttrName + "_prob", 2, ds);  // PAR

  if (merge_close)
    mergeClose (*clustNominAttr, *probAttr, 1e-2);  // PAR 
    
  const size_t missingCategory = clustNominAttr->missing2category ("Outlier");
  
  clustNominAttr->qc ();
  probAttr->qc ();
  
  sp1 << clustNominAttr;
  sp1 << probAttr;           

  if (getOutDim () < 2)
    return;


  if (! merge_close)
  {
    cerr << "Canonical ..." << endl;  
    Unverbose unv;
    const Canonical can (*this);
    can. qc ();
    cout << "Canonical projections = " << can. getOutDim () << endl;
    if (can. getOutDim ())
      sp1 << can. createSpace ("CC_", ds);
  }
  

  FFOR (size_t, i, clustNominAttr->categories. size ())
  {
    ostringstream oss;
    {
      Sample subset (sample);
      FFOR (size_t, row, space. ds. objs. size ())
        if ((*clustNominAttr) [row] != i)
          subset. mult [row] = 0.0;
      subset. finish ();    
      subset. save (attrs_orig, oss);
    }

    JsonMap* jClust = nullptr;
    if (jClusts)
    {
      jClust = new JsonMap (jClusts);
      new JsonString (clustNominAttr->categories [i], jClust, "name");
      auto* jDeps = new JsonArray (jClust, "dependencies");
      for (const Attr* attr : attrs_orig)
        if (const RealAttr1* attr1 = attr->asRealAttr1 ())
        {
          const NominAttr1::Dependence dep (clustNominAttr->getDependence (sample, i, attr1));
          if (leReal (dep. pValue, attr_pvalue))
          {
            JsonMap* jDep = dep. toJson (jDeps);
            new JsonString (attr1->name, jDep, "attr");
            new JsonInt    ((int) attr1->decimals, jDep, "decimals");
          }
        }
      if (i != missingCategory && oss. str (). size () < 10000000)  // PAR
        new JsonString (oss. str (), jClust, "dm");
    }

    if (! outDir. empty ()) 
    {
      // Can be used as fName for mds.cpp or pca.cpp
      OFStream f (outDir, outGenericDm + "." + toString (i + 1), dmExt);
      f << oss. str ();
    }

    if (   ! merge_close
        && i != missingCategory
       )
    {
      if (verbose ())
        cerr << "Canonical for cluster " << i + 1 << " ..." << endl;  

      Vector<bool> toMerge; 
      FFOR (size_t, j, getOutDim ())
        toMerge. push_back (i != j);
                
      Unverbose unv;

      const Clustering cl (*this, toMerge); 
      ASSERT (cl. getOutDim () == 2);
      
      const Canonical can (cl);
      can. qc ();


      // Between-center projection
      {
        const string attrPrefix ("BC" + toString (i + 1) + "_");
        
      	auto* attr = new RealAttr1 (attrPrefix + "1", ds);
      	sp1 << attr;
        MVector center (getMultiNormal (i) -> mu);
        center. normalizeRow (true, 0);
        for (Iterator it (sample); it ();)  
        {
          data2variable (*it);
          MVector vec (variable);
          vec. addRow (                   true, 0
                      , can. between. mu, true, 0  // global center
                      , -1.0
                      );
      	  (*attr) [*it] = vec. multiplyVec (center);
        }
        
        MVector center_diff (cl. getMultiNormal (0) -> mu);
        center_diff. addRow (                               true, 0
                            , cl. getMultiNormal (1) -> mu, true, 0
                            , -1.0
                            );

        if (jClust)
          new JsonDouble (center_diff. sumSqr (), numeric_limits<uint>::max(), jClust, "between_center_d2"); 

        const UniVariate<NumAttr1> analysis (sample, *attr);
        UniKernel uniKernel (analysis);
        uniKernel. estimate ();
        uniKernel. qc ();
      //cout << uniKernel. nameParam () << endl;
        auto* pdf = new PositiveAttr1 (attrPrefix + "2", ds, 3);  // PAR
        sp1 << pdf;
        for (Iterator it (sample); it ();)  
          (*pdf) [*it] = uniKernel. pdf ((*attr) [*it]);
      }
  
  
      if (can. getOutDim ())
      {
        ASSERT (can. getOutDim () == 1);
        const string attrPrefix ("CC" + toString (i + 1) + "_");
        
        const Space1<RealAttr1> spCan (can. createSpace (attrPrefix, ds));
        ASSERT (spCan. size () == 1);
        sp1 << spCan;
        if (verbose ())
        {
          cout << "Separating cluster " << i + 1 << ": " << can. eigenValues [0] << endl;
          can. print (cout);
        }
      #if 0
        Real average = NaN;
        Real scatter = NaN;
        spCan [0] -> getAverageScatter (sample, average, scatter);
        ASSERT (eqReal (average, 0, 1e-2));
      //ASSERT (eqReal (scatter, eigenValues [0] + 1, 1e-2));
      #endif
        if (jClust)
          new JsonDouble (/*scatter*/ can. eigenValues [0], 1, jClust, "canonical");  // PAR

        const RealAttr1& attr = * spCan [0]; 
        const UniVariate<NumAttr1> analysis (sample, attr);
        UniKernel uniKernel (analysis);
        uniKernel. estimate ();
        uniKernel. qc ();
      //cout << uniKernel. nameParam () << endl;
        auto* pdf = new PositiveAttr1 (attrPrefix + "2", ds);
        sp1 << pdf;
        for (Iterator it (sample); it ();)  
          (*pdf) [*it] = uniKernel. pdf (attr [*it]);
      }
    }
  }


#if 0
  if (classAttr)
  {
    Matrix contTab (false, clustNominAttr->categories. size (), classAttr->categories. size (), 0);
  //contTab. putAll (0);
    AttrAnalysis an (clustNominAttr, classAttr);
    for (Iterator it (an); it ();)  
      contTab. putInc (false, (*clustNominAttr) [*it], (*classAttr) [*it], 1);
    // Quality ??
    // Rename clustNominAttr->categories ??
    cout << endl;
    cout << clusterAttrName << " Class N" << endl;
    ONumber on (cout, 0, false);
    FFOR (size_t, row, contTab. rowsSize (false))
      FFOR (size_t, col, contTab. rowsSize (true))
        if (const Real v = contTab. get (false, row, col))
          cout << clustNominAttr->categories [row] << " " << classAttr->categories [col] << " " << v << endl;
  }
#endif
}




// Canonical

Canonical::Canonical (const Clustering &clustering)
: P (clustering. sample, clustering. space)
, requires (getRequirement (clustering))
, between (getBetween (clustering))
, within (getWithin (clustering))
, choleskyInv (getCholeskyInv (within. sigmaExact))
{ 
  if (! choleskyInv. defined ())
    return;

  MultiNormal mn (between);
  mn. sigmaExact. multiplyBilinear ( false
                                   , between. sigmaExact, false
                                   , choleskyInv, false
                                   );

  const size_t dim_max = clustering. getOutDim () - 1;
  ASSERT (dim_max);

  const Eigens eigens (mn. sigmaExact, dim_max, 1, 1e-3, 1e-4, 5000);  // PAR
  eigens. qc ();
  ASSERT (eigens. getInitSize () == choleskyInv. rowsSize (false));  
  ASSERT (eigens. values. size () <= dim_max);  
  
  basis. resize (false, between. getDim (), eigens. basis. rowsSize (true));
  basis. multiply ( false
                  , choleskyInv, false
                  , eigens. basis, false
                  );  

  eigenValues. resize (eigens. values. size ());
  eigenValues = eigens. values;

  basis_norm. resize (false, between. getDim (), eigens. basis. rowsSize (true));
  basis_norm = basis;
  {
    Real sumSqrNorma = NaN;
    EXEC_ASSERT (basis_norm. normalize (true, sumSqrNorma));
  }
}



MultiNormal Canonical::getBetween (const Clustering &clustering)
{
  MultiNormal between;

  Dataset ds;
  Vector<RealAttr1*> attrs;
  
  for (const Mixture::Component* comp : clustering. mixt. components)
    var_cast (ds. objs [ds. appendObj ()]) -> mult = comp->prob;
  FFOR (size_t, i, clustering. space. size ())
    attrs << new RealAttr1 (clustering. space [i] -> name, ds); 

  size_t objNum = 0;
  for (const Mixture::Component* comp : clustering. mixt. components)
  {
    const MultiNormal* compMn = comp->distr->asMultiNormal ();
    ASSERT (ds. attrs. size () == compMn->mu. size ());
    FFOR (size_t, i, clustering. space. size ())
      (*attrs [i]) [objNum] = compMn->mu [i];
    objNum++;
  }
  ASSERT (objNum == clustering. getOutDim ());

  {    
    const Sample dsSample (ds);
    const Space1<NumAttr1> dsSp (ds, true);
    const MultiVariate<NumAttr1> analysis (dsSample, dsSp);
    between. analysis = & analysis;
    between. estimate ();
  }
  
  return between;
}



MultiNormal Canonical::getWithin (const Clustering &clustering)
{
  MultiNormal within;
  within. setDim (clustering. space. size ());
  for (const Mixture::Component* comp : clustering. mixt. components)
  {
    const MultiNormal* compMn = comp->distr->asMultiNormal ();
    ASSERT (compMn);
    within. sigmaExact. add (false, compMn->sigmaExact, false, comp->prob);
  }
  within. setParam ();
  
  return within;
}



Matrix Canonical::getCholeskyInv (const Matrix &withinSigma) 
{
  Matrix m (withinSigma. getCholesky (true));
  if (! m. defined ())
    return m;
  m. qc ();
  ASSERT (m. isSquare ());
  
  if (! m. inverse (). get ())
    m. clear ();
    
  return m;
}



void Canonical::qc () const
{ 
  if (! qc_on)
    return;
  P::qc ();
    
  between. qc ();
  within. qc ();
  QC_ASSERT (between. getDim () == within. getDim ());
  QC_IMPLY (getOutDim (), choleskyInv. defined ());

  if (getOutDim ())
  {
    QC_ASSERT (basis. rowsSize (false) == between. getDim ());
    QC_ASSERT (basis. rowsSize (true) == eigenValues. size ());
    QC_ASSERT (basis_norm. equalLen (false, basis, false));
    QC_ASSERT (eigenValues. size () >= 1);

    choleskyInv. qc ();
    QC_ASSERT (choleskyInv. isSquare ());
    QC_ASSERT (between. getDim () == choleskyInv. rowsSize (false));
    basis. qc ();
    QC_ASSERT (basis. defined ());
    basis_norm. qc ();
    QC_ASSERT (basis_norm. defined ());
    eigenValues. qc ();
    eigenValues. defined ();
    QC_ASSERT (positive (eigenValues. min ()));
    
    if (false)  // --> getError () ??
    {
      Matrix m (getOutDim ());
      
      m. multiplyBilinear (false, within. sigmaExact, false, basis, false);
      FFOR (size_t, row, getOutDim ())
        FFOR (size_t, col, getOutDim ())
        {
          const Real a = m. get (false, row, col);
          if (row == col) 
            { if (! eqReal (a, 1, 1e-1))  // PAR
              {
                cout << a << endl;
                ERROR;
              }
            }
          else
            { if (! eqReal (a, 0, 1e-1))  // PAR
              {
                cout << a << endl;  
                ERROR;
              }
            }
        }   
        
      m. multiplyBilinear (false, between. sigmaExact, false, basis, false);
      FFOR (size_t, row, getOutDim ())
        FFOR (size_t, col, getOutDim ())
        {
          const Real a = m. get (false, row, col);
          if (row == col) 
            { QC_ASSERT (eqReal (a, eigenValues [row], 1e-2)); }
          else
            { QC_ASSERT (eqReal (a, 0, 1e-1)); }
        }   
    }
  }
}  



void Canonical::saveText (ostream &os) const
{ 
  os << endl << "Between:" << endl;
  between. saveText (os);

  os << endl << "Within:" << endl;
  within. saveText (os);
  
  if (choleskyInv. defined ())
  { 
    os << endl << "Cholesky:" << endl;
    choleskyInv. saveText (os);

    os << endl << "basis:" << endl;
    basis. saveText (os);

    os << endl << "eigenValues:" << endl;
    eigenValues. saveText (os);
  }
}


void Canonical::project (const Canonical::P::Value &x,
                         MVector &projection) const
{	
  ASSERT (getOutDim ());
  ASSERT (x. size () == space. size ());  

  MVector x_centered (x);
  x_centered. addRow (              true, 0
                     , between. mu, true, 0
                     , -1.0
                     );
  projection. multiply ( false
                       , basis_norm, true
                       , x_centered, false
                       ); 
}



Space1<RealAttr1> Canonical::createSpace (const string &attrPrefix,
                                          Dataset &ds) const
{
	ASSERT (! attrPrefix. empty ());
	ASSERT (& ds == sample. ds);
	
	Space1<RealAttr1> sp (ds, false);
	FFOR (size_t, i, getOutDim ())
	  sp << new RealAttr1 (attrPrefix + toString (i + 1), ds);
	
  MVector projection (getOutDim ());  
  FFOR (size_t, i, ds. objs. size ())
  {
    project (i, projection);  
  	FFOR (size_t, col, projection. size ())
  	  (* var_cast (sp [col])) [i] = projection [col];
  }
	
#ifndef NDEBUG
  if (sp. size () == 1)
  {
    Real average = NaN;
    Real scatter = NaN;
    sp [0] -> getAverageScatter (sample, average, scatter);
    ASSERT (eqReal (average, 0, 1e-2));
  //ASSERT (eqReal (scatter, eigenValues [0] + 1, 1e-2));
      // Requires: basis is used in project()
  }
#endif

	return sp;
}




// Mds

Mds::Mds (const Sample &sample_arg,
          const RealAttr2 &attr2,
		      size_t outDim_max,
	        Prob totalExplainedFrac_max,
	        Prob explainedFrac_min)
: Analysis (sample_arg)
, eigens ( attr2. matr
 	       , outDim_max
	       , totalExplainedFrac_max
	       , explainedFrac_min
	       , 1e-5  // PAR
	       , 10000  // PAR
	       )
{
  ASSERT (sample. ds == & attr2. ds);
  ASSERT (sample. ds->getUnitMult ());
//ASSERT (! attr2. matr. psd);  // PCA via MDS
}



void Mds::qc () const
{ 
  if (! qc_on)
    return;
  Analysis::qc ();
  
	eigens. qc (); 
	QC_IMPLY (getOutDim (), ! isNan (eigens. explainedVarianceFrac));
}



Space1<RealAttr1> Mds::createSpace (const string &realAttrPrefix,
                                    const string &imaginaryAttrPrefix,
                                    Dataset &ds) const
{
  ASSERT (! realAttrPrefix. empty ());
  ASSERT (! imaginaryAttrPrefix. empty ());
  ASSERT (realAttrPrefix != imaginaryAttrPrefix);
	ASSERT (& ds == sample. ds);
	
	Space1<RealAttr1> sp (ds, false);
	Vector<RealAttr1*> realAttrs;  realAttrs. reserve (getOutDim ());
	Vector<RealAttr1*> imaginaryAttrs (getOutDim (), nullptr);
	FFOR (size_t, i, getOutDim ())
	{
	  auto* realAttr = new RealAttr1 (realAttrPrefix + toString (i + 1), ds);;
	  sp << realAttr;
	  realAttrs << realAttr;
	  if (eigens. imaginaryMdsAttr (i))
	  {
  	  auto* imaginaryAttr = new RealAttr1 (imaginaryAttrPrefix + toString (i + 1), ds);;
  	  sp << imaginaryAttr;
	    imaginaryAttrs [i] = imaginaryAttr;
	  }
	}
	
  for (Iterator it (sample); it ();)  
  	FFOR (size_t, attrNum, eigens. values. size ())
  	{
  	  const Real x = get (*it, attrNum, false);
  	  if (isNan (x))
  	  {
  	    cout << "Object=" << *it << " Attr=" << attrNum << endl;
  	    cout << eigens. basis. get (false, *it, attrNum) << " " << eigens. values [attrNum] << endl;
  	    ERROR;
  	  }
  	  (*realAttrs [attrNum]) [*it] = x;
  	  if (imaginaryAttrs [attrNum])
  	    (*imaginaryAttrs [attrNum]) [*it] = get (*it, attrNum, true);
  	}
	
	return sp;
}




///////////////////////////// PositiveAverage //////////////////////////////////

// PositiveAverageModel::Component

PositiveAverageModel::Component::Component (const PositiveAverageModel& pam_arg,
                                            const string &line,
		                                        bool loadStat)
: pam (pam_arg)
{
	// Cf. saveText()
	istringstream iss (line);
	iss >> name; 
	if (loadStat)
	{
		iss >> coeff >> var;
		if (! (coeff >= 0.0))
			throw runtime_error (FUNC "Coefficient should be >= 0");
    if (! (coeff < INF))
			throw runtime_error (FUNC "Coefficient should be finite");
		if (! (var >= 0.0))
			throw runtime_error (FUNC "Variance should be >= 0");
  	setVar (var);
	}
}



void PositiveAverageModel::Component::qc () const
{
	if (! qc_on)
		return;

	Named::qc ();
	
	QC_ASSERT (goodName (name));
	QC_ASSERT (coeff >= 0.0);
	QC_ASSERT (var >= 0.0);
	QC_IMPLY (coeff == 0.0, var == INF);
}



void PositiveAverageModel::Component::setVar (Real var_arg)
{ 
  ASSERT (var_arg >= 0.0);

  var = var_arg;
  if (nullReal (var))
    throw runtime_error (FUNC "name " + strQuote (name) + ": too small variance");
    
  sd = sqrt (var);
  weight = 1.0 / var;
}



void PositiveAverageModel::Component::setOutlier (Real value_target) const
{ 
	outlier = (value == INF || /*fabs*/ (value - value_target) / sd > pam. outlierSEs);
}




// PositiveAverageModel

PositiveAverageModel::PositiveAverageModel (const string &fName,
                                            bool loadStat)
{ 
	LineInput f (fName);
	while (f. nextLine ())
	{
	  if (f. line. empty ())
	    continue;
	    
    if (loadStat && isNan (outlierSEs))
    {
      outlierSEs = stod (f. line);
      ASSERT (! isNan (outlierSEs));
    }
    else
		{ 
			Component comp (*this, f. line, loadStat);
		  components << move (comp);
		}
  }
}    



void PositiveAverageModel::qc () const
{
	if (! qc_on)
		return;

  QC_ASSERT (outlierSEs >= 0.0);
	for (const Component& comp : components)
		comp. qc ();
}



Real PositiveAverageModel::get () const
{
	Real ave = NaN;
  for (const Component& comp : components)
   	comp. outlier = false;  // To converge
	for (;;)
	{
		const Real ave_prev = ave;
		
    Real sum = 0.0;
    Real weights = 0.0;
    for (const Component& comp : components)
    	if (   Component::validValue (comp. value)
    		  && ! comp. outlier
    		 )
    	{
    		ASSERT (comp. value >= 0.0);
    		if (comp. value == INF)
    		{
    		  ASSERT (comp. var == INF);
    		  continue;
    		}
    		const Real weight = comp. weight;
    		ASSERT (weight >= 0.0);
    		ASSERT (weight < INF);
    		sum     += weight * comp. value;
    		weights += weight;
    	}  		
    ave = sum / weights;
    
    if (isNan (ave) || fabs (ave - ave_prev) < 1e-6)  // PAR
    	break;
    
    for (const Component& comp : components)
    	comp. setOutlier (ave);
  }
  maximize (ave, 0.0);
    
    
  return ave;
}



Matrix PositiveAverageModel::getParam () const
{ 
	Matrix param (components. size (), 2);
	FFOR (size_t, i, components. size ())
  {
	  param. put (false, i, 0, components [i]. coeff);
	  param. put (false, i, 1, components [i]. var);
	}
	return param;
}




// PositiveAverage

PositiveAverage::PositiveAverage (const Sample &sample_arg,
                                  const Space1<PositiveAttr1> &space_arg,
                                  Real outlierSEs_arg,
                                  bool optimizeCoeff_arg)
: P (sample_arg, space_arg)
, model (outlierSEs_arg)
, optimizeCoeff (optimizeCoeff_arg)
{
  // model
  model. components. reserve (space. size ());
  streamsize decimals_max = 0;
  for (const PositiveAttr1* attr : space)
  {
   	model. components << move (PositiveAverageModel::Component (attr->name, model));
    maximize (decimals_max, attr->decimals);
  }
  if (model. components. empty ())
    throw runtime_error (FUNC "PositiveAverage: No attributes");
    	  
  averageAttr = new PositiveAttr1 ("average", var_cast (* sample. ds), decimals_max + 1);

#if 0  
  outliers. resize (space. size ());
  for (Vector<bool>& vec : outliers)
    vec. resize (sample. ds->objs. size (), false);
#endif
  
  calibrate ();
}



void PositiveAverage::calibrate ()
{
  ASSERT (averageAttr);
  
  
  // Initialization 
  FFOR (size_t, attrNum, space. size ())
  {
    PositiveAverageModel::Component& comp = model. components [attrNum];
  	{
      const Real center = space [attrNum] -> getMedian (sample);
      ASSERT (center >= 0.0);
      if (center == 0.0)
        throw runtime_error (FUNC + space [attrNum] -> name + ": center is 0");
      comp. coeff = 1.0 / center;
    }
  	comp. setVar (0.1);  // PAR
  }

  
  for (;;)
  {
    // averageAttr
    // PositiveAverageModel::Component::outliers
  	FFOR (size_t, i, sample. ds->objs. size ())  
  	{
      FFOR (size_t, attrNum, space. size ())
		    model. components [attrNum]. setValue ((* space [attrNum]) [i]);
  		var_cast (*averageAttr) [i] = model. get ();
  	#if 0
      FFOR (size_t, attrNum, space. size ())
		    outliers [attrNum] [i] = model. components [attrNum]. outlier;
		#endif
  	}  	

    if (optimizeCoeff)
    {
      // Make average of averageAttr = 1
      Real average = NaN;
      Real scatter = NaN;
      averageAttr->getAverageScatter (sample, average, scatter);
      ASSERT (average > 0.0);
    	FFOR (size_t, i, sample. ds->objs. size ())
    		var_cast (*averageAttr) [i] /= average;  /* use errWeight ?? */
    }

    const Matrix param_old (model. getParam ());

    FFOR (size_t, attrNum, space. size ())
      setComponent (model. components [attrNum], * space [attrNum] /*, outliers [attrNum]*/);
 		
 		const Real diff = model. getParam (). maxAbsDiff (false, param_old, false);
    cerr << diff << endl;  /*+ " " + toString (getVar ())*/
 		if (diff < 1e-6)  // PAR
 			break;
 			
 	//saveText (cout);   
 	}
}



void PositiveAverage::setComponent (PositiveAverageModel::Component &comp,
                                    const PositiveAttr1 &attr/*,
                                    const Vector<bool> &attrOutliers*/)
{
	ASSERT (averageAttr);


	if (optimizeCoeff)
	{
  	// coeff
  	Real x2_sum = 0.0;
  	Real xy_sum = 0.0;
  	FFOR (size_t, i, sample. ds->objs. size ())
  	{
  	  const Real x = attr [i];
  		const Real y = (*averageAttr) [i];  
  		if (! PositiveAverageModel::Component::validValue (x))
  			continue;
  		if (! PositiveAverageModel::Component::validValue (y))
  			continue;
  		ASSERT (y >= 0.0);
  	#if 0
  		if (attrOutliers [i])
  		  continue;
  	#endif
  	#if 0
  		if (y == 0.0)
  		  continue;
  	  const Real errWeight = pow (y, - model. varPower);
  	#endif
  		const Real mult = sample. ds->objs [i] -> mult /* * errWeight*/;
  		ASSERT (mult >= 0.0);
  		x2_sum += mult * sqr (x);
  		xy_sum += mult * x * y;
  	}
  	
   	ASSERT (x2_sum >= 0.0);
  	if (x2_sum == 0.0)
  	  comp. coeff = INF;
  	else
  	{
    	if (xy_sum <= 0.0)
    	  comp. coeff = INF;
    	else
    	  comp. coeff = xy_sum / x2_sum;
    }
  }
	ASSERT (comp. coeff > 0.0);


  // var
  {
  	Real s = 0.0;
  	Real mult_sum = 0.0;
  	FFOR (size_t, i, sample. ds->objs. size ())
  	{
  	  const Real x = comp. setValue (attr [i]);  // y_hat
  		const Real y = (*averageAttr) [i];
  		if (! PositiveAverageModel::Component::validValue (x))
  			continue;
  		if (! PositiveAverageModel::Component::validValue (y))
  			continue;
  		ASSERT (y >= 0.0);
  	#if 0
  		if (attrOutliers [i])
  		  continue;
  	#endif
  	#if 0
  		if (y == 0.0)
  		  continue;
  		const Real errWeight = pow (y, - model. varPower);
  	#endif
  		const Real mult = sample. ds->objs [i] -> mult /* * errWeight*/;
  		ASSERT (mult >= 0.0);
  		s        += mult * sqr (x - y);  
  		mult_sum += mult;  
  	}

  	ASSERT (mult_sum >= 0.0);
  	if (mult_sum == 0.0)
  		comp. setVar (INF);
  	else
    	comp. setVar (s / mult_sum);
  }
	

	comp. qc ();
}



void PositiveAverage::qc () const
{
	if (! qc_on)
		return;

  P::qc ();
  model. qc ();
  QC_ASSERT (model. components. size () == space. size ());
//QC_ASSERT (model. components. size () == outliers. size ());  
  FFOR (size_t, i, space. size ())
  {
    QC_ASSERT (model. components [i]. name == space [i] -> name);
  //QC_ASSERT (outliers [i]. size () == sample. ds->objs. size ());
  }
  QC_ASSERT (averageAttr);
  QC_ASSERT (& averageAttr->ds == sample. ds);
}




}
