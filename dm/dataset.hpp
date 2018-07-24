// dataset.hpp
// "Data Master" library

#ifndef DATASET_HPP_63472  
#define DATASET_HPP_63472

#include "../common.hpp"
using namespace Common_sp;
#include "numeric.hpp"
#include "matrix.hpp"
using namespace DM_sp;



namespace DM_sp
{



extern const string dmExt;
extern const string dmSuff;  // = "." + dmExt
extern const char* missingStr;



// For time: n = Dataset::objs.size()



struct Obj : Root, DisjointCluster
// Object
{
  string name;
    // Me be empty()
  Real mult {1};
    // >= 0
    // Multiplicity (absolute frequency)
    // --> Sample Dataset::sample, remove extra Sample's ??
  string comment;


  explicit Obj (const string &name_arg) 
    : name (name_arg)
    {}
  void qc () const override;
};



struct RealScale
// Real numbers, finite()
{
  typedef  Real  Value;
  static const Value missing;  
    // NAN
  streamsize decimals;
    // Purposes: display, data precision


protected:    
  explicit RealScale (streamsize decimals_arg)
    : decimals (decimals_arg)
    {}
public:
  bool operator== (const RealScale &other) const
    { return decimals == other. decimals; }
  
  Real getPrecision () const
    { return pow (10.0, - (Real) decimals); }
};



constexpr streamsize decimals_def = 4;  // PAR

    
  

// Attr

struct Dataset;  

struct Sample;

struct Attr1;
  struct NumAttr1;
    struct RealAttr1;
      struct PositiveAttr1;
      struct ProbAttr1;
    struct IntAttr1;
    struct BoolAttr1;
      struct ExtBoolAttr1;
      struct CompactBoolAttr1;
  struct NominAttr1;
//struct OrdAttr1;
struct Attr2;
  struct RealAttr2;
  struct PositiveAttr2;


// Distribution's
struct Categorical;
struct LocScaleDistribution;



struct Attr : Root, Nocopy 
// Attribute
{
//typedef ... Value;
//static const Value missing;
  const Dataset &ds;  
private:
  List<const Attr*>/*=Dataset::Attrs*/::iterator dsIt;
    // In ds.attrs
  friend struct Dataset;
public:
  string name;
    // May be empty()
  bool rightAlign;
    // true iff the value strings are to be right-aligned 
/*
  static const DataType missing;
    // Missing value
  values of DataType
    // Init: missing
*/


protected:
  // Insert into ds.attrs
  Attr (const string &name_arg,
        Dataset &ds_arg,
        bool rightAlign_arg);
    // Invokes: Dataset::addAttr()
    // Descendants invoke: setMissingAll()
/*
  Attr (const string &name_arg,
        Dataset &ds_arg,
        const Attr &from);
    // Copy
*/
 ~Attr ();
public:
  void qc () const;
  void saveText (ostream &os) const
    { os << getTypeStr (); }
    

  virtual const Attr1* asAttr1 () const
    { return nullptr; }
  virtual const NumAttr1* asNumAttr1 () const
    { return nullptr; }
  virtual const RealAttr1* asRealAttr1 () const
    { return nullptr; }
  virtual const PositiveAttr1* asPositiveAttr1 () const
    { return nullptr; }
  virtual const ProbAttr1* asProbAttr1 () const
    { return nullptr; }
  virtual const IntAttr1* asIntAttr1 () const
    { return nullptr; }
  virtual const BoolAttr1* asBoolAttr1 () const
    { return nullptr; }
  virtual const ExtBoolAttr1* asExtBoolAttr1 () const
    { return nullptr; }
  virtual const CompactBoolAttr1* asCompactBoolAttr1 () const
    { return nullptr; }
  virtual const NominAttr1* asNominAttr1 () const
    { return nullptr; }
  virtual const Attr2* asAttr2 () const
    { return nullptr; }
  virtual const RealAttr2* asRealAttr2 () const
    { return nullptr; }
  virtual const PositiveAttr2* asPositiveAttr2 () const
    { return nullptr; }
  static const Attr* as (const Attr* attr)
    { return attr; }

  // Type
  virtual string getTypeStr () const = 0;
    // Requires: Corresponds to <dmSuff>-format
    
  virtual Attr* copyAttr (const string &name_arg) const = 0;
  virtual Attr* copyToDataset (Dataset &other) const = 0;
    
  // Position
  void moveAfter (const Attr* pred);
    // Update:: ds.attrs
    
  // Missings
  virtual bool isMissing (size_t objNum) const = 0;
  size_t countMissings () const;
  bool existsMissing (size_t &objNum) const;
    // Output: objNum; !Return <=> NO_INDEX
  bool existsMissing () const
    { size_t objNum;
      return existsMissing (objNum); 
    }
  bool missingsAll () const;
  virtual void setMissingAll () = 0;
  
  // Input: values
  virtual bool isConstant () const = 0;
    // Disregarding missing values
    // Requires: !missingAll()
  virtual size_t getWidth_max () const = 0; 
    // Maximum width of a value string
  // Update: values
  virtual void appendObj () = 0;
};



struct Attr1 : Attr
// Attribute of 1 object (unary attribute)
{
  // Data
//Vector<DataType> values;
  // size() == ds->objs.size() 


protected:
  Attr1 (const string &name_arg,
         Dataset &ds_arg,
         bool rightAlign_arg) 
    : Attr (name_arg, ds_arg, rightAlign_arg)
    {}
public:


  const Attr1* asAttr1 () const final
    { return this; }
  static const Attr1* as (const Attr* attr)
    { return attr->asAttr1 (); }

  void setMissingAll () final;
    // Invokes: setMissing()
  virtual void setMissing (size_t objNum) = 0;

  // values: I/O
  size_t getWidth_max () const; 
  virtual string value2str (size_t objNum) const = 0;
    // Requires: !isMissing(objNum)
  virtual void str2value (size_t objNum,
                          const string &s) = 0;
    // Output: values[objNum]
    // Requires: s does not encode missing
  virtual Json* value2json (JsonContainer* parent, 
                            size_t objNum) const = 0;
    // Return: new
    // Requires: !isMissing(objNum)

  virtual Vector<NumAttr1*> toNumAttr1 (Dataset &ds_arg) const = 0;
    // Invokes: new NumAttr1(,ds,)
  virtual Vector<RealAttr1*> standardize (Dataset &ds_arg,
                                          const Sample &sample) const = 0;
    // Make attributes centered and normalized
    // Invokes: new RealAttr1(,ds,)
};



struct NumAttr1 : Attr1, RealScale
{
protected:
  NumAttr1 (const string &name_arg,
            Dataset &ds_arg,
            streamsize decimals_arg)
    : Attr1 (name_arg, ds_arg, true)
    , RealScale (decimals_arg)
    {}
  NumAttr1 (const string &name_arg,
            Dataset &ds_arg,
            const NumAttr1 &from);
public:
  NumAttr1& operator= (const NumAttr1& other);


  const NumAttr1* asNumAttr1 () const final
    { return this; }
  static const NumAttr1* as (const Attr* attr)
    { return attr->asNumAttr1 (); }
    
  Vector<NumAttr1*> toNumAttr1 (Dataset &ds_arg) const final;
  Vector<RealAttr1*> standardize (Dataset &ds_arg,
                                  const Sample &sample) const override;
  virtual Value getReal (size_t objNum) const = 0;

  Value operator[] (size_t objNum) const
    { return getReal (objNum); }
  Value& operator[] (size_t objNum);

  // Properties
  bool existsLessThan (const Sample &sample,
                       Value minValue,
                       size_t &objNum) const;
    // Output: objNum: valid if return is true
  void getMinMax (const Sample &sample,
                  Value &min, 
                  Value &max) const;
  void getAverageScatter (const Sample &sample,
                          Real &average,
                          Real &scatter) const;
  Real getSqr_ave (const Sample &sample) const;
    // Return: >= 0
  Value getMedian (const Sample &sample) const;
  Real distr2outlier (const Sample &sample,
                      LocScaleDistribution &distr,
                      bool rightTail,
                      Real outlier_EValue_max) const;
    // Return: (if rightTail) min(x) s.t. (1-\Phi(x)) * mult_sum <= outlier_EValue_max; may be NAN or INF
    // Idempotent after removing outliers
    // Time: O(n log(n))
};



struct RealAttr1 : NumAttr1
// Missing = NAN
{
protected:
  Vector<Value> values;
public:

  
  RealAttr1 (const string &name_arg,
             Dataset &ds_arg,
             streamsize decimals_arg = decimals_def); 
  RealAttr1 (const string &name_arg,
             Dataset &ds_arg,
             const RealAttr1 &from);
  RealAttr1& operator= (const RealAttr1& other);
  void qc () const override;


  const RealAttr1* asRealAttr1 () const final
    { return this; }
  static const RealAttr1* as (const Attr* attr)
    { return attr->asRealAttr1 (); }
    
  RealAttr1* copyAttr (const string &name_arg) const override
    { return new RealAttr1 (name_arg, const_cast <Dataset&> (ds), *this); }
  RealAttr1* copyToDataset (Dataset &other) const override
    { return new RealAttr1 (name, other, *this); }

  string getTypeStr () const override
    { return "Real " + toString (decimals); }
  bool isConstant () const final;  
  void appendObj () final;
  bool isMissing (size_t objNum) const final
    { return isNan (values [objNum]); }
  void setMissing (size_t objNum) final
    { values [objNum] = missing; }
  string value2str (size_t objNum) const final
    { return real2str (values [objNum], decimals); }
  void str2value (size_t objNum,
                  const string &s) final
    { string s1 (s); 
      replaceStr (s1, ",", "");
      values [objNum] = str2real (s1); 
    }
  JsonDouble* value2json (JsonContainer* parent,
                          size_t objNum) const final
    { return new JsonDouble ((*this) [objNum], decimals, parent, name); }
  Value getReal (size_t objNum) const final
    { return (*this) [objNum]; }

  Value operator[] (size_t objNum) const
    { return values [objNum]; }
  Value& operator[] (size_t objNum) 
    { return values [objNum]; }
  
  static bool isMissingValue (Value value)
    { return isNan (value); }
  void setAll (Real value)
    { values. setAll (value); }
  size_t getInfCount () const;
  size_t inf2missing ();
    // Return: number of replacements
    
  void multiply (Real coeff);
  Real normal_likelihood2max (const Sample &sample) const;
    // Model: distribution of *this = if value <= t then truncated Normal else Uniform
    // Return: t, may be NAN
    // Time: O(n log(n))
};



struct PositiveAttr1 : RealAttr1
{
  PositiveAttr1 (const string &name_arg,
                 Dataset &ds_arg,
                 streamsize decimals_arg = decimals_def)  
    : RealAttr1 (name_arg, ds_arg, decimals_arg)
    {}
  PositiveAttr1 (const string &name_arg,
                 Dataset &ds_arg,
                 const PositiveAttr1 &from)
    : RealAttr1 (name_arg, ds_arg, from)
    {}
  PositiveAttr1& operator= (const PositiveAttr1& other)
    { RealAttr1::operator= (other);
      return *this;
    }
  void qc () const override;


  const PositiveAttr1* asPositiveAttr1 () const final
    { return this; } 
  static const PositiveAttr1* as (const Attr* attr)
    { return attr->asPositiveAttr1 (); }
    
  PositiveAttr1* copyAttr (const string &name_arg) const final
    { return new PositiveAttr1 (name_arg, const_cast <Dataset&> (ds), *this); }
  PositiveAttr1* copyToDataset (Dataset &other) const final
    { return new PositiveAttr1 (name, other, *this); }
    
  string getTypeStr () const final
    { return "Positive " + toString (decimals); }
  Vector<RealAttr1*> standardize (Dataset &ds_arg,
                                  const Sample &sample) const final;
    // Invokes: logarithmize()

  PositiveAttr1* standardizePositive (Dataset &ds_arg,
	                                    const Sample &sample) const;
	  // Return: nullptr <=> no data
	  //         average = 1
  RealAttr1* logarithmize (Dataset &ds_arg,
                           const string &suffix = "_log") const;
    // 0 --> missing
};



struct ProbAttr1 : RealAttr1
{
  ProbAttr1 (const string &name_arg,
             Dataset &ds_arg,
             streamsize decimals_arg = decimals_def)
    : RealAttr1 (name_arg, ds_arg, decimals_arg)
    {}
  ProbAttr1 (const string &name_arg,
             Dataset &ds_arg,
             const ProbAttr1 &from)
    : RealAttr1 (name_arg, ds_arg, from)
    {}
  ProbAttr1& operator= (const ProbAttr1& other)
    { RealAttr1::operator= (other);
      return *this;
    }
  void qc () const final;


  const ProbAttr1* asProbAttr1 () const final
    { return this; }
  static const ProbAttr1* as (const Attr* attr)
    { return attr->asProbAttr1 (); }
    
  ProbAttr1* copyAttr (const string &name_arg) const final
    { return new ProbAttr1 (name_arg, const_cast <Dataset&> (ds), *this); }
  ProbAttr1* copyToDataset (Dataset &other) const final
    { return new ProbAttr1 (name, other, *this); }

  string getTypeStr () const final
    { return "Probability " + toString (decimals); }

  Prob getProb (size_t objNum) const;
};



struct IntAttr1 : NumAttr1
{
  typedef int Value;
  static const Value missing;  
private:
  Vector<Value> values;
public:


  IntAttr1 (const string &name_arg,
            Dataset &ds_arg);
  IntAttr1 (const string &name_arg,
            Dataset &ds_arg,
            const IntAttr1 &from);
  IntAttr1& operator= (const IntAttr1& other);
  void qc () const final;


  const IntAttr1* asIntAttr1 () const final
    { return this; }
  static const IntAttr1* as (const Attr* attr)
    { return attr->asIntAttr1 (); }
    
  IntAttr1* copyAttr (const string &name_arg) const final
    { return new IntAttr1 (name_arg, const_cast <Dataset&> (ds), *this); }
  IntAttr1* copyToDataset (Dataset &other) const final
    { return new IntAttr1 (name, other, *this); }

  string getTypeStr () const final
    { return "Integer"; }
  bool isConstant () const final;
  void appendObj () final
    { values. push_back (missing); }
  bool isMissing (size_t objNum) const final
    { return values [objNum] == missing; }
  void setMissing (size_t objNum) final
    { values [objNum] = missing; }
  string value2str (size_t objNum) const final
    { return toString (values [objNum]); }
  JsonInt* value2json (JsonContainer* parent,
                       size_t objNum) const final
    { return new JsonInt ((*this) [objNum], parent, name); }

  Value operator[] (size_t objNum) const
    { return values [objNum]; }
  Value& operator[] (size_t objNum) 
    { return values [objNum]; }

  void setAll (Value value)
    { values. setAll (value); }

  void str2value (size_t objNum,
                  const string &s)
    { (*this) [objNum] = str2<Value> (s); }
  Real getReal (size_t objNum) const 
    { const Value v = (*this) [objNum];
      return v == missing ? NumAttr1::missing : (Real) v; 
    }
};



struct BoolAttr1 : NumAttr1
{
  typedef ebool Value;
  static const Value missing;  


protected:
  BoolAttr1 (const string &name_arg,
             Dataset &ds_arg)
    : NumAttr1 (name_arg, ds_arg, 0)
    {}
public:
  

  const BoolAttr1* asBoolAttr1 () const final
    { return this; }
  static const BoolAttr1* as (const Attr* attr)
    { return attr->asBoolAttr1 (); }
    
  void str2value (size_t objNum,
                  const string &s) final
    { setBool (objNum, (ebool) str2bool (s)); }
  JsonBoolean* value2json (JsonContainer* parent,
                           size_t objNum) const final
    { return new JsonBoolean (getBool (objNum), parent, name); }
  Real getReal (size_t objNum) const final
    { const ebool v = getBool (objNum);
      return v == UBOOL ? NumAttr1::missing : (Real) v; 
    }
    
  virtual Value getBool (size_t objNum) const = 0;
  virtual void setBool (size_t objNum,
                        Value value) = 0;
  bool str2bool (const string &s) const;

  Value operator[] (size_t objNum) const
    { return getBool (objNum); }
  Value& operator[] (size_t objNum);

  // Properties
  void getStat (const Sample &sample,
                array<size_t,3/*ebool*/> &stat) const;
    // Output: stat
  Prob getProb (const Sample &sample) const;
};



struct ExtBoolAttr1 : BoolAttr1
{
private:
  Vector<Value> values;
public:
    

  ExtBoolAttr1 (const string &name_arg,
                Dataset &ds_arg);
  ExtBoolAttr1 (const string &name_arg,
                Dataset &ds_arg,
                const ExtBoolAttr1 &from);
  ExtBoolAttr1& operator= (const ExtBoolAttr1& other);
  void qc () const final;


  const ExtBoolAttr1* asExtBoolAttr1 () const final
    { return this; }
  static const ExtBoolAttr1* as (const Attr* attr)
    { return attr->asExtBoolAttr1 (); }
    
  ExtBoolAttr1* copyAttr (const string &name_arg) const final
    { return new ExtBoolAttr1 (name_arg, const_cast <Dataset&> (ds), *this); }
  ExtBoolAttr1* copyToDataset (Dataset &other) const final
    { return new ExtBoolAttr1 (name, other, *this); }

  string getTypeStr () const final
    { return "Boolean"; }
  bool isConstant () const final;
  void appendObj () final
    { values. push_back (UBOOL); }
  bool isMissing (size_t objNum) const final
    { return values [objNum] == UBOOL; }
  void setMissing (size_t objNum) final
    { values [objNum] = UBOOL; }
  string value2str (size_t objNum) const final
    { return values [objNum] ? "1" : "0"; }
  Value getBool (size_t objNum) const final 
    { return (*this) [objNum]; }
  void setBool (size_t objNum,
                Value value) final 
    { (*this) [objNum] = value; }

  Value operator[] (size_t objNum) const
    { return values [objNum]; }
  Value& operator[] (size_t objNum) 
    { return values [objNum]; }

  void setAll (Value value)
    { values. setAll (value); }
};



struct CompactBoolAttr1 : BoolAttr1
{
private:
  vector<bool> values;
    // Init: false
public:
    

  CompactBoolAttr1 (const string &name_arg,
                    Dataset &ds_arg)
    : BoolAttr1 (name_arg, ds_arg)
    { setAll (false);}
  CompactBoolAttr1 (const string &name_arg,
                    Dataset &ds_arg,
                    const CompactBoolAttr1 &from);
  CompactBoolAttr1& operator= (const CompactBoolAttr1& other);
  void qc () const final;


  const CompactBoolAttr1* asCompactBoolAttr1 () const final
    { return this; }
  static const CompactBoolAttr1* as (const Attr* attr)
    { return attr->asCompactBoolAttr1 (); }
        
  CompactBoolAttr1* copyAttr (const string &name_arg) const final
    { return new CompactBoolAttr1 (name_arg, const_cast <Dataset&> (ds), *this); }
  CompactBoolAttr1* copyToDataset (Dataset &other) const final
    { return new CompactBoolAttr1 (name, other, *this); }

  string getTypeStr () const final
    { return "CompactBoolean"; }
  bool isConstant () const final;
  void appendObj () final
    { values. push_back (false); }
  bool isMissing (size_t /*objNum*/) const final
    { return false; }
  void setMissing (size_t /*objNum*/) final
    { throw runtime_error ("setMissing for CompactBoolAttr1"); }
  string value2str (size_t objNum) const final
    { return values [objNum] ? "1" : "0"; }
  Value getBool (size_t objNum) const final 
    { return (Value) (*this) [objNum]; }
  void setBool (size_t objNum,
                Value value) final
    { if (value == missing)
        throw runtime_error ("setBool(missing) for CompactBoolAttr1");
      setCompactBool (objNum, value);
    }

  bool operator[] (size_t objNum) const
    { return values [objNum]; }
  void setCompactBool (size_t objNum,
                       bool value)
    { values [objNum] = value; }
  void setAll (bool value);
};



struct NominAttr1 : Attr1
// Nominal (aka categorical, qualitative)
// TOTEST
{
  typedef size_t Value;
    // Category index
  static const Value missing;  
    // NO_INDEX
  StringVector categories;
    // Index: Value
private:
  typedef map<string,Value> CategMap;
  CategMap categMap;
  Vector<Value> values;
    // values[] < categories.size()
public:
    

  NominAttr1 (const string &name_arg,
              Dataset &ds_arg);
  NominAttr1 (const string &name_arg,
              Dataset &ds_arg,
              const NominAttr1 &from);
  NominAttr1& operator= (const NominAttr1& other);
  void qc () const final;


  const NominAttr1* asNominAttr1 () const final
    { return this; }
  static const NominAttr1* as (const Attr* attr)
    { return attr->asNominAttr1 (); }
    
  NominAttr1* copyAttr (const string &name_arg) const final
    { return new NominAttr1 (name_arg, const_cast <Dataset&> (ds), *this); }
  NominAttr1* copyToDataset (Dataset &other) const final
    { return new NominAttr1 (name, other, *this); }

  string getTypeStr () const final;
  bool isConstant () const final
    { return categories. size () == 1; }
  void appendObj () final
    { values. push_back (missing); }
  bool isMissing (size_t objNum) const final
    { return values [objNum] == missing; }
  void setMissing (size_t objNum) final
    { values [objNum] = missing; }
  string value2str (size_t objNum) const final
    { return categories [values [objNum]]; }
  void str2value (size_t objNum,
                  const string &s) final;
  JsonString* value2json (JsonContainer* parent,
                          size_t objNum) const final
    { return new JsonString (value2str (objNum), parent, name); }
  Vector<NumAttr1*> toNumAttr1 (Dataset &ds_arg) const final;
  Vector<RealAttr1*> standardize (Dataset &ds_arg,
                                  const Sample &sample) const final;

  Value operator[] (size_t objNum) const
    { return values [objNum]; }
  Value& operator[] (size_t objNum) 
    { return values [objNum]; }

  Value category2index (const string &category);
    // Update: categMap
private:
  void rebuildCategMap ();
    // Input: categories
    // Output: categMap
public:
  void setAll (size_t value)
    { values. setAll (value); }
  Categorical* getCategorical (const Sample &sample) const;
  bool isDegenerate (const Sample &sample) const;
  void deleteEmptyCategories (const Sample &sample);
  Value missing2category (const string &missingName);
    // Return: new category or NO_INDEX
    // Update: categories
  void renameCategories ();
    // Update: categories: C<i>
  
  
  struct Dependence : Root  // --> Analysis ??
  {
    struct NormalParam : Root
    { 
      Real mean; 
      Real sd; 
      streamsize decimals;
      NormalParam (Real mean_arg, 
                   Real sd_arg,
                   streamsize decimals_arg)
        : mean (mean_arg)
        , sd (sd_arg)
        , decimals (decimals_arg)
        {}
      void saveText (ostream &os) const final
        { os << mean << " +- " << sd; }
      JsonMap* toJson (JsonContainer* parent_arg,
                       const string& name_arg = noString) const final
        { auto jMap = new JsonMap (parent_arg, name_arg);
          new JsonDouble (mean, decimals, jMap, "mean");
          new JsonDouble (sd,   decimals, jMap, "sd");
          return jMap;
        }
    };
    Prob pValue;
    NormalParam cluster;
    NormalParam other;
    Dependence (Prob pValue_arg,
                NormalParam cluster_arg,
                NormalParam other_arg)
      : pValue (pValue_arg)
      , cluster (cluster_arg)
      , other (other_arg)
      {}
    void qc () const final;
    void saveText (ostream &os) const final
      { os << "pValue=" << pValue << "  cluster ";
        cluster. saveText (os);
        os << "  other ";
        other. saveText (os);
      }
    JsonMap* toJson (JsonContainer* parent_arg,
                     const string& name_arg = noString) const final
      { auto jMap = new JsonMap (parent_arg, name_arg);
        new JsonDouble (pValue, 3, jMap, "pValue");  // PAR
        cluster. toJson (jMap, "cluster");
        other.   toJson (jMap, "other");
        return jMap;
      }
  };
  Dependence getDependence (const Sample &sample,
                            size_t value,
                            const RealAttr1* attr) const;
};



#if 0
struct OrdAttr1 : NominAttr1
{
  OrdAttr1 (const string &name_arg,
            Dataset &ds_arg)
    : NominAttr1 (name_arg, ds_arg)
    {}
  OrdAttr1 (const string &name_arg,
            Dataset &ds_arg,
            const OrdAttr1 &from)
    : NominAttr1 (name_arg, ds_arg, from)
    {}
};
#endif



struct Attr2 : Attr
// Attribute of 2 objects (two-way attribute), object-object table
{
//Vector<Vector<DataType> > values;
  // size() == ds->objs.size() * ds->objs.size()


protected:
  Attr2 (const string &name_arg,
         Dataset &ds_arg,
         bool rightAlign_arg) 
    : Attr (name_arg, ds_arg, rightAlign_arg)
    {}
public:


  const Attr2* asAttr2 () const final
    { return this; }
  static const Attr2* as (const Attr* attr)
    { return attr->asAttr2 (); }

  bool isMissing (size_t objNum) const final;
  void setMissingAll () final;
    // Invokes: setMissing()
  size_t getWidth_max () const final; 

  // Missings
  virtual bool isMissing2 (size_t row,
                           size_t col) const = 0;
  bool existsMissing2 (size_t &row,
                       size_t &col) const;
    // Output: row,col; !Return <=> NO_INDEX
  virtual void setMissing (size_t row,
                           size_t col) = 0;
  // values: I/O
  virtual string value2str (size_t row,
                            size_t col) const = 0;
    // Requires: !isMissing2(row,col)
  virtual void str2value (size_t row,
                          size_t col,
                          const string &s) = 0;
    // Output: values[row][col]
    // Requires: s does not encode missing

  void save (ostream &os,
             const Sample &sample) const;
};



struct RealAttr2 : Attr2, RealScale
{
  Matrix matr;

  
  RealAttr2 (const string &name_arg,
             Dataset &ds_arg,
             streamsize decimals_arg = decimals_def);
  RealAttr2 (const string &name_arg,
             Dataset &ds_arg,
             const RealAttr2 &from);
  RealAttr2& operator= (const RealAttr2& other);
  void qc () const override;
  void saveText (ostream &os) const final
    { Attr2::saveText (os);
      os << endl;
      matr. print (os);
    }


  const RealAttr2* asRealAttr2 () const final
    { return this; }
  static const RealAttr2* as (const Attr* attr)
    { return attr->asRealAttr2 (); }
    
  RealAttr2* copyAttr (const string &name_arg) const override
    { return new RealAttr2 (name_arg, const_cast <Dataset&> (ds), *this); }
  RealAttr2* copyToDataset (Dataset &other) const override
    { return new RealAttr2 (name, other, *this); }

  string getTypeStr () const override
    { return "Real2 " + toString (decimals); }
  bool isConstant () const final;  
  void appendObj () final;
  bool isMissing2 (size_t row,
                   size_t col) const final
    { return isNan (get (row, col)); }
  void setMissing (size_t row,
                   size_t col) final
    { put (row, col, missing); }
  string value2str (size_t row,
                    size_t col) const final
    { return real2str (get (row, col), decimals); }
  void str2value (size_t row,
                  size_t col,
                  const string &s) final
    { string s1 (s);
      replaceStr (s1, ",", "");
      put (row, col, str2real (s1)); 
    }
    
  Value get (size_t row,
             size_t col) const
    { return matr. get (false, row, col); }
  void put (size_t row,
            size_t col, 
            Value value) 
    { matr. put (false, row, col, value); }
  void putSymm (size_t row,
                size_t col, 
                Value value) 
    { put (row, col, value); 
      if (row != col)
        put (col, row, value); 
    }
  
  bool existsLessThan (Value minValue,
                       size_t &row,
                       size_t &col) const;
    // Output: row,col: valid if return is true
  void setAll (Value value)
    { matr. putAll (value); }
  size_t getInfCount () const;
  size_t inf2missing ();
    // Return: number of replacements
};



struct PositiveAttr2 : RealAttr2
{
  PositiveAttr2 (const string &name_arg,
                 Dataset &ds_arg,
                 streamsize decimals_arg = decimals_def)
    : RealAttr2 (name_arg, ds_arg, decimals_arg)
    {}
  PositiveAttr2 (const string &name_arg,
                 Dataset &ds_arg,
                 const PositiveAttr2 &from)
    : RealAttr2 (name_arg, ds_arg, from)
    {}
  PositiveAttr2& operator= (const PositiveAttr2& other)
    { RealAttr2::operator= (other);
      return *this;
    }
  void qc () const override;


  const PositiveAttr2* asPositiveAttr2 () const final
    { return this; }
  static const PositiveAttr2* as (const Attr* attr)
    { return attr->asPositiveAttr2 (); }
    
  PositiveAttr2* copyAttr (const string &name_arg) const final
    { return new PositiveAttr2 (name_arg, const_cast <Dataset&> (ds), *this); }
  PositiveAttr2* copyToDataset (Dataset &other) const final
    { return new PositiveAttr2 (name, other, *this); }

  string getTypeStr () const final
    { return "Positive2 " + toString (decimals); }

  bool areClose (size_t row,
                 size_t col,
                 Real distance_max) const
    // Input: distance_max: >= 0
    { const Real d = get (row, col);
      return DM_sp::finite (d) && (! distance_max || d <= distance_max);
    }
};



struct Dataset : Root
// Object-attribute table
// Time series: the greater objNum the later time point
{
  List<string> comments;
  VectorOwn<Obj> objs;  // --> Vector<Obj> ??
    // Requires: Obj::name's are different
  List<const Attr*> attrs;
    // Attr::name's are different or empty()
private:
  map<string/*Attr::name*/,const Attr*> name2attr_;
  map<string/*Obj::name*/,size_t/*Obj index*/> name2objNum;
public:


  struct Error : runtime_error
    { explicit Error (const string &str) : runtime_error ("Dataset: " + str) {} };


  Dataset () 
    {}
  explicit Dataset (const string &fName)
    { { ifstream is (fName + dmSuff);
        load (is);
      }
      qc ();
    }
  explicit Dataset (istream &is)
    { load (is);
      qc ();
    }
private:
  void load (istream &is);
    // Loading from a text file in the <dmSuff>-format:
    //
    //    {# <Comment>}
    //    ObjNum <# objects> [no]name [no]mult
    //    Attributes
    //      {<Attribute name> <Type>}
    //    Data
    //      {[<Object name>] [<mult>] {<Value>}}
    //   [COMMENT
    //      {<Object comment>}
    //   ]
    //
    // where <Type> ::= {Real | Positive | Probability} <# decimals> |
    //                  Nominal {<Category name>}+   // line may wrap if ended with "\"
    //                  {Integer | Boolean | CompactBoolean}
    //
    // All keywords are case-insensitive
    // Missings are coded by <missing>
    // Empty lines are allowed
    // Invokes: qc()
public:
  explicit Dataset (const Eigens &eigens);
    // RealAttr1: unit, "Order" (log), "EigenValueFrac" (log)
 ~Dataset ()
    { deleteAttrs (); }
  void qc () const;
  void saveText (ostream &os) const;


  bool empty () const
    { return objs. empty () && attrs. empty (); }
  JsonArray* comments2Json (JsonContainer* parent,
                            const string& name) const;   
    // Return: may be nullptr                        

  // objs
  bool name2objNumSet () const
    { return ! name2objNum. empty (); }
  void setName2objNum ();
    // Output: name2objNum
    // If duplicate names then throw
  size_t getName2objNum (const string &objName) const;
    // Return: NO_INDEX <=> not found
    // Input: name2objNum
  StringVector getObjNames () const
    { StringVector names;  names. reserve (name2objNum. size ());
      for (const auto it : name2objNum)
        names << it. first;
      names. sort ();
      return names;
    }
  size_t appendObj (const string &objName = string ());
    // Return: objNum
  void list2ObjNames (const VectorPtr<Named> &names);
  bool objCommented () const;
  
  string pair2name (size_t objNum1,
                    size_t objNum2) const
    { string s1 (objs [objNum1] -> name);
      string s2 (objs [objNum2] -> name);
      if (s1 > s2)
        swap (s1, s2);
      return s1 + "-" + s2;
    }

  // --> Sample ??
  // Obj::mult
  void setMultAll (Real mult);
  void attr2mult (const RealAttr1* attr);
  void objInterval2Active (size_t startObjNum,
                           size_t endObjNum);
  bool getUnitMult () const;
    // Return: true <=> for all objNum obj[objNum]->mult = 1

  // attrs
  void deleteAttrs ();
  const Attr* name2attr (const string &attrName) const
    { return attrName. empty () ? nullptr : findPtr (name2attr_, attrName); }
    // Return: nullptr <=> not found
private:
  void addAttr (Attr* attr);  
  friend struct Attr;
public:
  RealAttr1* addRealAttr1Unit (const string &attrName = "unit")
    { auto attr = new RealAttr1 (attrName, *this, 0);
      attr->setAll (1);
      attr->moveAfter (nullptr);
      return attr;
    }
};


  
struct Sample : Root
// Usage: Sample s(); [s.mult *= ...; s.finish();]  // *= --> setMult() ??
{
  const Dataset* ds;
    // !nullptr
  Vector<Real> mult;
    // size() = ds.objs.size()
  Real mult_sum {0};
    // >= 0
  size_t nEffective {0};


  explicit Sample (const Dataset &ds_arg);
    // Invokes: finish()
  void finish ();
    // Input: mult[]
    // Output: mult_sum, nEffective
  void qc () const override;
  bool empty () const override
    { return mult. empty (); }


  size_t size () const
    { return mult. size (); } 
  Real getMult_avg () const
    { return mult_sum / (Real) nEffective; }
  Real getMaxMult () const;
  size_t getObjNameLen_max () const;  
  void missing2mult (const Attr1* attr1);
  void save (const VectorPtr<Attr> &attrs,
             ostream &os) const;
    // <dmSuff>-format
};



struct SaveSample : Sample
/* Usage:
   struct X : SaveSample
     { X (Sample &sample_arg, ...) { <modify sample_arg> } };
   { X x (sample);
     <use sample>
   }
*/
{
protected:
  Sample& sample; 

  explicit SaveSample (Sample &sample_arg)
    : Sample (sample_arg)
    , sample (sample_arg)
    {}
public:
 ~SaveSample ()
    { sample = *this; }
};



struct Iterator : Nocopy  // --> FOR ??
{
private:
  const Vector<Real>& mult_;
  // Output
  size_t objNum {NO_INDEX};
public:
  Real mult {NAN};

  explicit Iterator (const Sample &sample)
    : mult_ (sample. mult)
    {}

  bool operator() ()
    // Update: objNum, mult
    // Return: true => mult
		{ if (objNum == mult_. size ())
		    return false;
		  for (objNum = objNum == NO_INDEX ? 0 : (objNum + 1); 
		       objNum < mult_. size (); 
		       objNum++
		      )
		  { mult = mult_ [objNum];
		    if (mult)
		      return true;
		  }
		  return false;
		}
  size_t operator* () const 
    { return objNum; }
    // Valid if operator()
};



template <typename T/*:Attr*/>
  struct Space : VectorPtr<T>
  // For running time: n = ds.objs.size(), p = size()
  // Elements: in ds.attrs, !nullptr, all are different
  {
    typedef  VectorPtr<T>  P;
    const Dataset& ds;
  
  
    Space (const Dataset &ds_arg,
           bool copyAttrs)
      : ds (ds_arg)
      { if (copyAttrs)
        { P::reserve (ds. attrs. size ());
          for (const Attr* attr : ds. attrs)
            if (const T* a = T::as (attr))
              *this << a;
        }
      }
    template <typename U/*:T*/>
      Space (const Space<U> &other)
        : P (other)
        , ds (other. ds)
        {}
    void qc () const
      { if (! qc_on)
          return;
        ASSERT (! P::empty ());
        Set<const T*> attrSet;
        for (const T* attr : *this)
        { ASSERT (attr);
          if (& attr->ds != & ds)
          { cout << & attr->ds << " " << & ds << " " << attr->name << endl;
            ERROR;
          }
          attrSet << attr;
        }
        ASSERT (P::size () == attrSet. size ());
      }
  
  
    const T* name2attr (const string &attrName) const
      { if (! attrName. empty ())
          for (const T* attr : *this)
            if (attr->name == attrName)
              return attr;
        return nullptr;
      }
    size_t getMaxNameLen () const
      { size_t len = 0;
        for (const Attr* attr : *this)
          maximize (len, attr->name. size ());
        return len;
      }
  
    template <typename U/*:T*/>
      Space<T>& operator<< (const U* attr)
        { ASSERT (attr);
          ASSERT (& attr->ds == & ds);
        //if (name2attr (attr->name))
          //ERROR_MSG ("Attribute " + strQuote (attr->name) + " already exists in the space");
          P::operator<< (attr);  
          return *this;
        }
    template <typename U/*:T*/>
      Space<T>& operator<< (const Space<U> &other)
        { for (const U* attr : other)
            *this << attr;
          return *this;
        }
    
    const T* existsMissing (size_t &badObjNum) const
      { for (const T* attr : *this)
          if (attr->existsMissing (badObjNum))
            return attr;
        badObjNum = NO_INDEX;
        return nullptr;
      }
    bool existsMissing () const
      { size_t badObjNum;
      	return existsMissing (badObjNum); 
      }
    VectorPtr<T> removeConstants ()
      { VectorPtr<T> erased;
        for (Iter<VectorPtr<T> > iter (*this); iter. next (); )
          if ((*iter)->isConstant ())
          { erased << *iter;
            iter. erase ();
          }
        return erased;
      }
    Space<T>& removeAttr (const T &attr)
      { for (Iter<VectorPtr<T> > iter (*this); iter. next (); )
          if (*iter == & attr)
            iter. erase ();
        return *this;
      }
      
  protected:
    const Dataset* dsPtr () const
      { return & ds; }
  };
  
    

template <typename T/*:Attr1*/>
  struct Space1 : Space<T>
  {
    typedef  Space<T>  P;
    
        
    Space1 (const Dataset &ds_arg,
            bool copyAttrs)
      : P (ds_arg, copyAttrs)
      {}
    template <typename U/*:T*/>
      Space1 (const Space<U> &other)
        : P (other)
        {}


    void printCsv (const Sample &sample,
                   ostream &os) const
      { ASSERT (sample. ds == P::dsPtr ());
        const bool unitMult = P::ds. getUnitMult ();
        // Header
        os << "ObjName";
        if (! unitMult)
          os << ",Mult";
        for (const T* attr : *this) 
          os << "," << attr->name;
        os << endl;
        // Data
        for (Iterator it (sample); it ();)  
        { os << P::ds. objs [*it] -> name;
          if (! unitMult)
            os << "," << fixed << sample. mult [*it];
          for (const T* attr : *this)
            os << "," << (attr->isMissing (*it) ? missingStr : attr->value2str (*it)); 
          os << endl;
          ASSERT (os. good ());
        }
      }
    JsonArray* toJson (const Sample &sample,
                       JsonContainer* parent,
                       const string& name = noString) const 
      // Element format: {objName:S,mult:R,comment:S,attr:(<attrName>:S)*}
      { ASSERT (sample. ds == P::dsPtr ());
        const bool unitMult = P::ds. getUnitMult ();
        const bool commented = P::ds. objCommented ();
        auto jObjs = new JsonArray (parent, name);
        for (Iterator it (sample); it ();)  
        { const Obj* obj = P::ds. objs [*it];
          auto jObj = new JsonMap (jObjs);
          new JsonString (obj->name, jObj, "objName");
          if (! unitMult)
            new JsonDouble (sample. mult [*it], 3, jObj, "mult");  // PAR
          if (commented)
            new JsonString (obj->comment, jObj, "comment");
          auto jAttr = new JsonMap (jObj, "attr");
          for (const T* attr : *this)
          { Json* j = nullptr;
            if (attr->isMissing (*it))
              j = new JsonNull (jAttr, attr->name);
            else
              j = attr->value2json (jAttr, *it); 
            ASSERT (j);
          }
        }
        return jObjs;
      }

    Space1<NumAttr1> toNumAttr1 (Dataset &ds_arg) const
      { ASSERT (P::dsPtr () == & ds_arg);
        Space1<NumAttr1> sp (ds_arg, false);
        for (const T* attr : *this)
        { const Vector<NumAttr1*> vec (attr->toNumAttr1 (ds_arg));
          for (const NumAttr1* attr_new : vec)
            sp << attr_new;
        }
        return sp;
      }
    Space1<RealAttr1> standardize (const Sample &sample,
                                   Dataset &ds_arg) const
      { ASSERT (P::dsPtr () == sample. ds);
        ASSERT (P::dsPtr () == & ds_arg);
        Space1<RealAttr1> sp (ds_arg, false);
        for (const T* attr : *this)
          { const Vector<RealAttr1*> vec (attr->standardize (ds_arg, sample));
            for (const RealAttr1* attr_new : vec)
              sp << attr_new;
          }
        return sp;
      }
  };




// Return: !nullptr
// Requires: space.defined()

PositiveAttr2* getDist2 (const Space1<RealAttr1> &space,
                         const string &attrName,
                         Dataset &ds);
  // Return: distance in L_2
  
PositiveAttr2* getHammingDist (const Space1<ProbAttr1> &space,
                               const string &attrName,
                               Dataset &ds);
  // Return: distance in L_1
  // Input: if values \in {0,1} then Return = getDist2()
  
RealAttr2* getSimilarity (const Space1<RealAttr1> &space,
                          const string &attrName,
                          Dataset &ds);



  
///////////////////////////////////////////////////////////////////////////////

struct Analysis : Root
{
  Sample sample;
//const Attr* <attributes>
//typedef ... Value;
    // Matches <attributes>
//mutable Value variable;


protected:
  explicit Analysis (const Sample &sample_arg)
    : sample (sample_arg)
    {}
public:
  void qc () const override
    { if (! qc_on)
        return;
      sample. qc (); 
    }
};



struct Analysis1 : Analysis  
//const Attr1* <attributes>;
{
protected:
  explicit Analysis1 (const Sample &sample_arg)
    : Analysis (sample_arg)
    {}
public:


  virtual void data2variable (size_t objNum) const = 0;
    // Input: Dataset
    // Output: variable
  virtual void variable2data (size_t objNum) = 0;
    // Input: variable
    // Output: Dataset
};



template <typename T/*:Attr1*/>
  struct UniVariate : Analysis1
  {
    const T& attr;

    typedef  typename T::Value  Value;
    mutable Value variable;
    
  
    UniVariate (const Sample &sample_arg,
                const T &attr_arg)
      : Analysis1 (sample_arg)
      , attr (attr_arg)
      , variable ()
      {}
    void qc () const override
      { if (! qc_on)
          return;
        Analysis1::qc();
        attr. qc ();
        ASSERT (& attr. ds == sample. ds);
      }


    void data2variable (size_t objNum) const final
      { variable = attr [objNum]; }
    void variable2data (size_t objNum) final
      { const_cast <T&> (attr) [objNum] = variable; }
  };
  
  

template <typename T/*:Attr1*/>
  struct MultiVariate : Analysis1
  {
    Space1<T> space;
    
    typedef  Vector<typename T::Value>  Value;
    mutable Value variable;

  
    template <typename U/*:T*/>
      MultiVariate (const Sample &sample_arg,
                    const Space1<U> &space_arg)
        : Analysis1 (sample_arg)
        , space (space_arg)
        , variable (space. size ())
        {}
    void qc () const override
      { if (! qc_on)
          return;
        Analysis1::qc();
        space. qc ();
        ASSERT (& space. ds == sample. ds);
        ASSERT (space. size () == variable. size ()); 
        ASSERT (! space. empty ());
      }


    Real size (uint samplePower,
               uint spacePower) const
      { return   pow ((Real) sample. nEffective, samplePower) 
               * pow ((Real) space. size (),     spacePower); 
      }

    void data2variable (size_t objNum) const override
      { FOR (size_t, i, variable. size ())
          variable [i] = (* space [i]) [objNum]; 
      }
    void variable2data (size_t objNum) override
      { FOR (size_t, i, variable. size ())
          const_cast <T&> (* space [i]) [objNum] = variable [i]; 
      }

    void saveData (ostream &os) const
      { const VectorPtr<Attr> attrs (space);
        sample. save (attrs, os);
      }


    // Requires: T:NumAttr1
  private:
	  static void setAttrSim_ (size_t from,
					                   size_t to,
					                   const Space1<T> &space,  
	                           const Sample &sample,
	                           Matrix &attrSim,
				                     bool vc) 
	    // Input: vc: variance-covariance: allows missing values
		  {	ASSERT (attrSim. rowsSize () == space. size ());
		    FFOR_START (size_t, attrNum1, from, to)
		    { const NumAttr1& a1 = * space [attrNum1];
		  	  FFOR_START (size_t, attrNum2, attrNum1, space. size ())
		      { const NumAttr1& a2 = * space [attrNum2];
		        Real s = 0;
		        Real mult_sum = 0;
		        FFOR (size_t, i, sample. mult. size ())
		          if (const Real mult = sample. mult [i])
		          { const Real x1 = a1. getReal (i);
		          	const Real x2 = a2. getReal (i);
		          	if (vc)
		            {
			          	if (isNan (x1))
			          		continue;
			          	if (isNan (x2))
			          		continue;
			            mult_sum += mult;
			          }
		            s += x1 * x2 * mult;
		          }
		        ASSERT (! isNan (s));
		        attrSim. putSymmetric (attrNum1, attrNum2, vc ? (mult_sum ? s / mult_sum : 0) : s);
		      }
		   	}
		  }		
		public:
	  void setAttrSim (Matrix &attrSim,
	                   bool vc) const
	    // Input: vc: variance-covariance: allows missing values
		  {	ASSERT (attrSim. rowsSize () == space. size ());
		    Unverbose unv;
		  #if 0
		    // thread_num depends on sqrt(sample.mult.size ())
		    arrayThreads (setAttrSim_, space. size (), std::cref (space), std::cref (sample), std::ref (attrSim), vc);
		  #else
		    Progress prog (space. size ());  
		    FFOR (size_t, attrNum1, space. size ())
		    { prog ();
		      const NumAttr1& a1 = * space [attrNum1];
		  	  FFOR_START (size_t, attrNum2, attrNum1, space. size ())
		      { const NumAttr1& a2 = * space [attrNum2];
		        Real s = 0;
		        Real mult_sum = 0;
		        FFOR (size_t, i, sample. mult. size ())
		          if (const Real mult = sample. mult [i])
		          { const Real x1 = a1. getReal (i);
		          	const Real x2 = a2. getReal (i);
		          	if (vc)
		            {
			          	if (isNan (x1))
			          		continue;
			          	if (isNan (x2))
			          		continue;
			            mult_sum += mult;
			          }
		            s += x1 * x2 * mult;
		          }
		        ASSERT (! isNan (s));
		        attrSim. putSymmetric (attrNum1, attrNum2, vc ? (mult_sum ? s / mult_sum : 0) : s);
		      }
		   	}
			#endif
		  }
  };    




//////////////////////// Distribution ////////////////////////

//struct Distribution
  struct Bernoulli;
  struct Categorical;
  struct UniDistribution;
    struct DiscreteDistribution;  
      struct Binomial;
      struct UniformDiscrete;
      struct Geometric;
      struct Zipf;
    struct ContinuousDistribution;
      struct LocScaleDistribution;
        struct Normal;
        struct Exponential;
        struct Cauchy;
      struct Beta1;
      struct UniKernel;
  struct MultiDistribution;
    struct MultiNormal;
  struct Mixture;



struct Distribution : Named  
// Random variable
{ 
//typedef  <Analysis1>  An;
//const An* analysis;
//mutable An::Value variable;
    // Instantiation of the random variable
    // missings ??
//Parameters
  //Confidence interval of Parameters ??
//Functions of parameters
  //Sufficient statistics ??

protected:
  mutable Rand randGen;
    // Random numbers generator
public:
    

protected:
  explicit Distribution (const string &name_arg)
    : Named (name_arg)
    {}
public:
  Distribution* copy () const override = 0;
  void qc () const override;


  virtual const UniDistribution* asUniDistribution () const
    { return nullptr; }
  virtual const DiscreteDistribution* asDiscreteDistribution () const
    { return nullptr; }
  virtual const Bernoulli* asBernoulli () const
    { return nullptr; }
  virtual const Categorical* asCategorical () const
    { return nullptr; }
  virtual const Binomial* asBinomial () const
    { return nullptr; }
  virtual const UniformDiscrete* asUniformDiscrete () const
    { return nullptr; }
  virtual const Geometric* asGeometric () const
    { return nullptr; }
  virtual const Zipf* asZipf () const
    { return nullptr; }
  virtual const ContinuousDistribution* asContinuousDistribution () const
    { return nullptr; }
  virtual const LocScaleDistribution* asLocScaleDistribution () const
    { return nullptr; }
  virtual const Normal* asNormal () const
    { return nullptr; }
  virtual const Exponential* asExponential () const
    { return nullptr; }
  virtual const Cauchy* asCauchy () const
    { return nullptr; }
  virtual const Beta1* asBeta1 () const
    { return nullptr; }
  virtual const UniKernel* asUniKernel () const
    { return nullptr; }
  virtual const MultiDistribution* asMultiDistribution () const
    { return nullptr; }
  virtual const MultiNormal* asMultiNormal () const
    { return nullptr; }
  virtual const Mixture* asMixture () const
    { return nullptr; }
    

  virtual const Analysis1* getAnalysis () const = 0;
  const Analysis1* getAnalysisCheck () const
    { if (const Analysis1* an = getAnalysis ())
        return an;
      throw runtime_error ("No data to analyze a distribution");
    }
  // Output: <analysis>
  virtual Analysis1* createAnalysis (Dataset &ds) = 0;
    // Return: new <Analysis1>
    // Output: analysis = Return
    // Update: ds: new <Attr>
  virtual void removeAnalysis () = 0;
    // Output: analysis = nullptr    
  virtual void shareAnalysis (const Distribution &distr) = 0;
    // Output: analysis = distr.analysis if possible: may be nullptr
  // Input: *getAnalysisCheck()
  virtual void variable2analysis () const = 0;
    // Input: variable
    // Output: analysis->variable
  virtual void analysis2variable () const = 0;
    // Input: analysis->variable
    // Output: variable
  void variable2data (size_t objNum)
    { variable2analysis ();
      const_cast <Analysis1*> (getAnalysisCheck ()) -> variable2data (objNum);
    }
  void data2variable (size_t objNum) const
    { getAnalysisCheck () -> data2variable (objNum);
      analysis2variable ();
    }

  // Update: randGen
  virtual void setSeed (ulong seed) const
    { randGen. setSeed (seed); }
  Prob randProb () const
    { return randGen. getProb (); }
    // Return: Uniform on Prob
    // Requires: after setSeed()

//void setParam (<Parameters>_arg);    
    // Do:
    //   <Parameters> = <Parameters>_arg;
    //   finish();
protected:
  virtual void finish ()
    { setParamFunc ();
      qc ();
    }
  virtual void setParamFunc () 
    // Input: <Parameters>
    // Output: <Functions of parameters>
    {}
public:
  virtual void estimate () = 0;
    // Input: *getAnalysisCheck()
    // Update: <Parameters>
    // MLE
    // Requires: qc()
    // Postcondition: getParamSet()
    // Invokes: setParam()
    // Skip missing values ??!

  virtual bool getParamSet () const = 0;
    // Input: <Parameters>
  
  // Requires: getParamSet()

  virtual size_t getDim () const 
    { return 1; }
    // Return: # dimensions (= # scalar random variables); >= 1
  virtual size_t paramCount () const
    { return getDim (); }
  virtual bool similar (const Distribution &distr,
                        Real delta) const = 0;
    // a.similar(b) <=> a ~ b
  virtual Real getSortingValue () const = 0;
    // To sort() Distribution's

  // Input: variable
  virtual Real pdfVariable () const
    { return exp (logPdfVariable ()); }
    // Return: >= 0
  virtual Real logPdfVariable () const
    { return log (pdfVariable ()); }
  size_t getModeObjNum () const;
    
  virtual void randVariable () const = 0;
    // Output: variable
    // Invokes: randProb()    
  void simulate (Dataset &ds,
                 size_t objsSize);
    // Output: ds, analysis
    // Requires: ds.empty()
    // Invokes: createAnalysis(), randVariable(), variable2data()

  // Information = -log(pdfVariable(X))
  virtual Real getInfoMean () const
    { return NAN; }
  Real getEntropy () const  // synonym
    { return getInfoMean (); }  
    // If X ~ Categorical with n equally probable categories then getEntropy() = log(n)
  virtual Real getInfoVar () const
    { return NAN; }    

  // Input: *getAnalysisCheck()
  
  void info_meanVar_MC (uint sampleSize,
                        Real &mean,
                        Real &var) const;
    // Monte-Carlo
    // Output: mean, var
    // Invokes: randVariable()    
  void info_meanVar (uint sampleSize,
                     Real &mean,
                     Real &var) const;
    // Output: mean, var
    // Invokes: getInfoMean(), getInfoVar(), info_meanVar_MC()
    
  Real getLogLikelihood () const;
  Real getEntropy_est () const
    { return - getLogLikelihood () / getAnalysisCheck () -> sample. mult_sum; }
  bool minimizeEntropy (Real &entropy_best,
                        Real delta) const
    { const Real entropy = getEntropy_est ();
      const Real deltaAdj = delta * (Real) getDim ();
      if (verbose ())
        cout << "new entropy = " << entropy << "  deltaAdj = " << deltaAdj << endl;
      return minimizeReal (entropy_best, entropy, deltaAdj); 
    }
    // Entropy depends on variance => a relative measurement 
  Normal getEntropyDistribution () const;
    // Invokes: info_meanVar()
  // Entropy fitness test of *this
  // Return: p-value 
  // Invoke: getEntropyDistribution()
  Prob getFitness_entropy (Real entropy_est_best) const;
    // Return: P(Entropy <= entropy_est_best)
  Prob getFitness_entropy () const; 
    // 2-tail Normal test
    // Return: P(getEntropy_est() is unlikely for Entropy)
  //
//Prob getFitness_KolmogorovSmirnov ();  // ??

  size_t getWeakestObjNum () const;
};



struct Bernoulli : Distribution 
{
  typedef  UniVariate<BoolAttr1>  An;
  const An* analysis {nullptr};
  mutable An::Value variable {BoolAttr1::missing};
  // Parameters
  Prob p {NAN};


  Bernoulli () 
    : Distribution ("Bernoulli") 
    {}
  Bernoulli* copy () const final
    { return new Bernoulli (*this); }
  void qc () const override;
  void saveText (ostream &os) const override
    { os << name << "(" << prob2str (p) << ")"; }
  
  
  const Bernoulli* asBernoulli () const final
    { return this; }

  const Analysis1* getAnalysis () const final
    { return analysis; }
  Analysis1* createAnalysis (Dataset &ds) final;
  void removeAnalysis () final
    { analysis = nullptr; }
  void shareAnalysis (const Distribution &distr) final
    { analysis = nullptr;
      if (const Bernoulli* b = distr. asBernoulli ())
        analysis = b->analysis;
    }
  void variable2analysis () const final
    { checkPtr (analysis) -> variable = variable; }
  void analysis2variable () const final
    { variable = checkPtr (analysis) -> variable; }
  // Parameters
  void setParam (Prob p_arg)    
    { p = p_arg; 
      finish ();
    }
  void estimate () final
    { NOT_IMPLEMENTED; }
  bool getParamSet () const 
    { return ! isNan (p); }
  bool similar (const Distribution &distr,
                Real delta) const;
  Real getSortingValue () const 
    { return p; }
  Real pdfVariable () const
    { return pmf (variable); }
  void randVariable () const
    { variable = randProb () <= p ? ETRUE : EFALSE; }

  Prob pmf (An::Value x) const
    { switch (x)
      { case ETRUE:  return p;
        case EFALSE: return 1 - p;
        case UBOOL:  return NAN;
      }
      throw runtime_error ("Unknown ebool value");
    }
  Real getMean () const  
    { return p; }
  Real getVar () const  
    { return p * (1 - p); }
};



struct Categorical : Distribution 
// Zipf property: pmf(x) = c * pow(rank(pmf(x)),-alpha)
//                <=> rank(pmf(X)) ~ Zipf
{
  typedef  UniVariate<NominAttr1>  An; 
  const An* analysis {nullptr};
  mutable An::Value variable {NominAttr1::missing};
  // Parameters
  Vector<Prob> probs;
    // sum = 1
  // Functions of parameters
private:
  Vector<Prob> probSum;
    // size() = probs.size()
    // at(0) = 0; at(i) < at(i+1)
public:


  Categorical () 
    : Distribution ("Categorical") 
    {}
  Categorical* copy () const final
    { return new Categorical (*this); }
  void qc () const;
  void saveText (ostream& os) const;

  
  const Categorical* asCategorical () const final
    { return this; }

  const Analysis1* getAnalysis () const 
    { return analysis; }
  Analysis1* createAnalysis (Dataset &ds);
  void removeAnalysis () 
    { analysis = nullptr; }
  void shareAnalysis (const Distribution &distr) 
    { analysis = nullptr;
      if (const Categorical* d = distr. asCategorical ())
        if (d->probs. size () == probs. size ())
          analysis = d->analysis;
    }
  void variable2analysis () const
    { checkPtr (analysis) -> variable = variable; }
  void analysis2variable () const
    { variable = checkPtr (analysis) -> variable; }
  // Parameters
  void setParam ()
    // Input: probs: proportional to the probabilities
    { finish (); }
private:
  void setParamFunc ();
public:
  void setParam (const DiscreteDistribution &distr,
                 int x_from,
                 int x_to);
    // Input: x_from, x_to: inclusive
  void estimate () final;
  bool getParamSet () const 
    { return ! probs. empty (); }
  size_t paramCount () const
    { return probs. size () - 1; }
  bool similar (const Distribution &distr,
                Real delta) const;
  Real getSortingValue () const 
    { return getNominalVar (); }
  //
  Real pdfVariable () const
    { return pmf (variable); }
  void randVariable () const
    { variable = probSum. binSearch (randProb (), false); }
    
  Prob pmf (size_t x) const
    { return probs [x]; }
  Prob getNominalVar () const 
    { Prob s = 0;
      for (const Prob p : probs)
        s += sqr (p);
      return s;
    }
  size_t getUniqueCategory () const;
    // Return: category i if it is the only one with a non-zero probs[i], otherwise NO_INDEX
private:
  void balanceProb ();
    // Update: probs
};



struct UniDistribution : Distribution
// Univariate distribution
// Standard support (unconditional distribution): [stdLoBound(), stdHiBound()]
// Bounded support (truncated distribution): [getLoBoundEffective(), getHiBoundEffective()]
{
  // Parameters
  Real loBound {-INF};
  Real hiBound {INF};

  // Functions of parameters
  Prob p_supp {1};
    // = P(getLoBoundEffective() <= X <= getHiBoundEffective())
  Prob p_ltSupp {0};
    // = P(X < getLoBoundEffective())
  // p_supp + p_ltSupp <= 1
protected:
  Real log_p_supp {0};
    // <= 0
public:
    

protected:
  explicit UniDistribution (const string &name_arg)
    : Distribution (name_arg) 
    {}
public:
  void qc () const;
  void saveText (ostream& os) const
    { os << nameParam (); 
      if (! stdBounds ())
        os << " on [" << getLoBoundEffective () << ',' << getHiBoundEffective () << ']';
    }
  virtual string nameParam () const = 0;
  
  
  const UniDistribution* asUniDistribution () const
    { return this; }

  // Parameters
protected:
  void finish ()
    { Distribution::finish ();
      log_p_supp = log (p_supp);
    }
  void setParamFunc () 
    { p_supp = 1;
      p_ltSupp = 0;
    }
public:
  bool similar (const Distribution &distr,
                Real delta) const 
    { if (const UniDistribution* u = distr. asUniDistribution ())
        return    eqReal (getLoBoundEffective (), u->getLoBoundEffective (), delta)
               && eqReal (getHiBoundEffective (), u->getHiBoundEffective (), delta);
      return false;
    }
  Real getSortingValue () const 
    { return getMean (); }

  // Input: bounds
  virtual Real stdLoBound () const
    { return -INF; }
  virtual Real stdHiBound () const
    { return INF; }
  bool stdBounds () const
    { return    leReal (loBound, stdLoBound ())
             && geReal (hiBound, stdHiBound ());
    }
  Real getLoBoundEffective () const
    { return max (loBound, stdLoBound ()); }
  Real getHiBoundEffective () const
    { return min (hiBound, stdHiBound ()); }
  bool supported (Real x) const
    { return    geReal (x, getLoBoundEffective ())
             && leReal (x, getHiBoundEffective ());
    }
    
  // Requires: getParamSet()

  // Standard support
private:
  virtual Real pdf_ (Real x) const
    { return exp (logPdf_ (x)); }
    // >= 0
  virtual Real logPdf_ (Real x) const
    { return log (pdf_ (x)); }
  virtual Prob cdf_ (Real x) const = 0;
    // Return: P(X<=x)
    // cdf_(stdHiBound()) = 1
  virtual Real rand_ () const = 0;  
    // Return: standard or bounded support
    // Update: randGen
public:

  // Bounded support
  Real pdf (Real x) const
    { if (p_supp == 0)
        return 0;
      if (supported (x))
        return pdf_ (x) / p_supp; 
      return 0;
    }
    // Return: >= 0
  Real logPdf (Real x) const
    { if (p_supp == 0)
        return -INF;
      if (supported (x))
        return logPdf_ (x) - log_p_supp; 
      return -INF;
    }
  Prob cdf (Real x) const
    { if (lessReal (x, getLoBoundEffective ()))
        return 0;
      if (geReal (x, getHiBoundEffective ()))
        return 1;
      return (cdf_ (x) - p_ltSupp) / p_supp;
    }
  Real rand () const
    { Real x;
      do x = rand_ ();
        while (! supported (x));
      return x;
    }

  // Standard support
  // Moments
  virtual Real getMean () const = 0;
  virtual Real getVar () const = 0;
  Real getSd () const
    { return sqrt (getVar ()); }    
};



struct DiscreteDistribution : UniDistribution
{
  typedef  UniVariate<IntAttr1>  An; 
  const An* analysis {nullptr};
  mutable An::Value variable {IntAttr1::missing};


protected:
  explicit DiscreteDistribution (const string &name_arg) 
    : UniDistribution (name_arg) 
    {}
public:
  void qc () const override;

    
  const DiscreteDistribution* asDiscreteDistribution () const final
    { return this; }

  const Analysis1* getAnalysis () const final
    { return analysis; }
  Analysis1* createAnalysis (Dataset &ds) final;
  void removeAnalysis () final
    { analysis = nullptr; }
  void shareAnalysis (const Distribution &distr) final
    { analysis = nullptr;
      if (const DiscreteDistribution* d = distr. asDiscreteDistribution ())
        if (   eqReal (d->getLoBoundEffective (), getLoBoundEffective ())
            && eqReal (d->getHiBoundEffective (), getHiBoundEffective ())
           )
          analysis = d->analysis;
    }
  void variable2analysis () const final
    { checkPtr (analysis) -> variable = variable; }
  void analysis2variable () const final
    { variable = checkPtr (analysis) -> variable; }
  // Parameters
  Real pdfVariable () const final
    { return pdf (variable); }
  Real logPdfVariable () const final
    { return logPdf (variable); }
  void randVariable () const final
    { variable = randDiscrete (); }

  Real stdLoBound () const override
    { return 0; }
private:
  Real pdf_ (Real x) const final
    { return pmf_ (real2int (x)); }
    // Mixture of delta-functions
  Real logPdf_ (Real x) const
    { return logPmf_ (real2int (x)); }
  Prob cdf_ (Real x) const
    { return cdfDiscrete_ (real2int (isInteger (x) ? x : floor (x))); }
  Real rand_ () const
    { return randDiscrete_ (); }
public:

protected:
  static int real2int (Real x);
    // Input: x is integer  
  // Standard support
  virtual Prob pmf_ (int x) const
    { return exp (logPmf_ (x)); }
  virtual Real logPmf_ (int x) const
    { return log (pmf_ (x)); }
    // Return: <= 0
  virtual Prob cdfDiscrete_ (int x) const = 0;
    // Return: P(X <= x)
  virtual int randDiscrete_ () const = 0;
    // Return: standard or bounded support
    // Update: randGen
public:
  
  // Bounded support
  Prob pmf (int x) const
    { return pdf (x); }
  Real logPmf (int x) const
    { return logPdf (x); }
    // Return: <= 0
  Prob cdfDiscrete (int x) const
    { return cdf (x); }
  int randDiscrete () const
    { return real2int (rand ()); }
};



struct Binomial : DiscreteDistribution 
// Sum of n Berboulli(p)
{
  // Parameters
  An::Value n {0};
    // > 0
    // 0 <=> unknown
  Prob p {NAN};
private:
  // Functions of parameters
  Real n_re {NAN};
  Bernoulli bernoulli;
  Real lnFacN {NAN};
  Real lnP {NAN};
  Real lnPCompl {NAN};
public:


  Binomial () 
    : DiscreteDistribution ("Binomial") 
    {}
  Binomial* copy () const final
    { return new Binomial (*this); }
  void qc () const final;
    

  const Binomial* asBinomial () const
    { return this; }

  void setSeed (ulong seed) const override
    { DiscreteDistribution::setSeed (seed);
      bernoulli. setSeed (seed + 1);
    }
  // Parameters
  void setParam (int n_arg,
                 Prob p_arg)    
    { n = n_arg;
      p = p_arg;
      finish ();
    }
private:
  void setParamFunc ();
public:
  void estimate () final;
  bool getParamSet () const 
    { return n && ! isNan (p); }
  size_t paramCount () const
    { return 2; }
  bool similar (const Distribution &distr,
                Real delta) const
    { if (! DiscreteDistribution::similar (distr, delta))
        return false;
      if (const Binomial* bin = distr. asBinomial ())
        return    n == bin->n
               && eqReal (p, bin->p, delta);
      return false;
    }
  // UniDistribution
  string nameParam () const final
    { return name + "(" + toString (n) + "," + prob2str (p) + ")"; }
  Real stdHiBound () const final
    { return n_re; }
  Real getMean () const  
    { return n_re * bernoulli. getMean (); }
  Real getVar () const  
    { return n_re * bernoulli. getVar (); }
private:
  Real logPmf_ (int x) const;
  Prob cdfDiscrete_ (int x) const;
  int randDiscrete_ () const;
};



struct UniformDiscrete : DiscreteDistribution 
// Continuous: Uniform
{
  // Parameters
  An::Value min {1};
  An::Value max {0};
  // min <= max


  UniformDiscrete () 
    : DiscreteDistribution ("Uniform(Discrete)") 
    {}
  void qc () const;
  UniformDiscrete* copy () const final
    { return new UniformDiscrete (*this); }


  const UniformDiscrete* asUniformDiscrete () const
    { return this; }

  // Parameters
  void setParam (int min_arg,
                 int max_arg)    
    { min = min_arg;
      max = max_arg;
      finish ();
    }
  void estimate () final
    { if (! getParamSet ())  // ??
    	  NOT_IMPLEMENTED; 
    } 
  bool getParamSet () const 
    { return min <= max; }
  size_t paramCount () const
    { return 2; }
  bool similar (const Distribution &distr,
                Real delta) const
    { if (! DiscreteDistribution::similar (distr, delta))
        return false;
      if (const UniformDiscrete* ud = distr. asUniformDiscrete ())
        return    min == ud->min
               && max == ud->max;
      return false;
    }
  // UniDistribution
  string nameParam () const final
    { return name + "(" + toString (min) + "," + toString (max) + ")"; }
  Real stdLoBound () const final
    { return (Real) min; }
  Real stdHiBound () const final
    { return (Real) max; }
  Real getMean () const  
    { return (Real) (min + max) / 2.0; }
  Real getVar () const  
    { return (Real) (range () + 1) * ((Real) range () - 1) / 12.0; }
private:
  Real pmf_ (int /*x*/) const
    { return 1.0 / range (); }          
  int randDiscrete_ () const
    { return min + (int) randGen. get ((ulong) range ()); }
  Prob cdfDiscrete_ (int x) const
    { return (Real) (x - min + 1) * 1.0 / range (); }
public:

  int range () const
    { return max - min + 1; }
};



struct Geometric : DiscreteDistribution 
// Values: >= 1
// Continuous: Exponential
{
  // Parameters
  Prob p {NAN};
  // Functions of parameters
private:
  Real lnP {NAN};
  Real lnPCompl {NAN};
public:


  Geometric () 
    : DiscreteDistribution ("Geometric") 
    {}
  Geometric* copy () const final
    { return new Geometric (*this); }
  void qc () const final;


  const Geometric* asGeometric () const final
    { return this; }

  // Parameters
  void setParam (Prob p_arg)    
    { p = p_arg;
      finish ();
    }
private:
  void setParamFunc ()
    { lnP = log (p);
      lnPCompl = log (1 - p);
    }
public:
  void estimate () final;
  bool getParamSet () const 
    { return ! isNan (p); }
  bool similar (const Distribution &distr,
                Real delta) const
    { if (! DiscreteDistribution::similar (distr, delta))
        return false;
      if (const Geometric* geom = distr. asGeometric ())
        return eqReal (p, geom->p, delta);
      return false;
    }
  // UniDistribution
  string nameParam () const final
    { return name + "(" + prob2str (p) + ")"; }
  Real stdLoBound () const final
    { return 1; }
  Real getMean () const  
    { return 1 / p; }
  Real getVar () const  
    { return (1 - p) / sqr (p); }
private:
  Real logPmf_ (int x) const;
  Prob cdfDiscrete_ (int /*x*/) const
    { NOT_IMPLEMENTED; return NAN; }
  int randDiscrete_ () const
    { NOT_IMPLEMENTED; return -1; }
};



struct Zipf : DiscreteDistribution   // not a distribution ??
// = Zeta 
// Zipf -> data -> frequency ranks ~ Zipf ??
// Zipf -> data -> frequencies !~ Zipf
// Beta1 -> discrete data !~ Zipf
// Continuous: Pareto
{
  // Parameters
  Real alpha {NAN};
    // > 1
private:
  // Functions of parameters
  Real c {NAN};
    // > 0
  Real lnC {NAN};
  Categorical cat;
    // Valid if hiBound < INF
public:


  Zipf () 
    : DiscreteDistribution ("Zipf") 
    {}
  Zipf* copy () const final
    { return new Zipf (*this); }
  void qc () const final;


  const Zipf* asZipf () const final
    { return this; }

  void setSeed (ulong seed) const override
    { DiscreteDistribution::setSeed (seed);
      cat. setSeed (seed + 1);
    }
  // Parameters
  void setParam (Real alpha_arg)
    { alpha = alpha_arg;
      finish ();
    }
private:
  void setParamFunc ();
public:   
  void estimate () final;
  bool getParamSet () const 
    { return ! isNan (alpha); }
  bool similar (const Distribution &distr,
                Real delta) const
    { if (! DiscreteDistribution::similar (distr, delta))
        return false;
      if (const Zipf* zipf = distr. asZipf ())
        return eqReal (alpha, zipf->alpha, delta);
      return false;
    }
  // UniDistribution
  string nameParam () const final
    { return name + "(" + real2str (alpha, 3) + ")"; }
  Real stdLoBound () const final
    { return 1; }
  Real getMean () const  
    { NOT_IMPLEMENTED; return NAN; }
  Real getVar () const  
    { NOT_IMPLEMENTED; return NAN; }
private:
  Prob pmf_ (int x) const
    { return c * pow (x, - alpha); }
  Real logPmf_ (int x) const;
  Prob cdfDiscrete_ (int /*x*/) const
    { NOT_IMPLEMENTED; return NAN; }
  int randDiscrete_ () const;
public:

#if 0
  // WRONG    
  Real freqAlpha2alpha () const  
    { return 1 / (alpha - 1); }
  Real getBetaAlpha () const
    // Return: alpha of the Beta(alpha,1), alpha > 0
    // pmf(1) = Beta1::cdf(1/upBound) = (1/upBound)^alpha 
    { return - log (c) / log (upBound); }
#endif
};



struct ContinuousDistribution : UniDistribution  
// Stable distributions ??
{
  typedef  UniVariate<NumAttr1>  An;
  const An* analysis {nullptr};
  mutable An::Value variable {RealScale::missing};
    

protected:
  explicit ContinuousDistribution (const string &name_arg)
    : UniDistribution (name_arg) 
    {}
public:
  
  
  const ContinuousDistribution* asContinuousDistribution () const final
    { return this; }

  const Analysis1* getAnalysis () const final
    { return analysis; }
  Analysis1* createAnalysis (Dataset &ds) final;
  void removeAnalysis () final
    { analysis = nullptr; }
  void shareAnalysis (const Distribution &distr) final
    { analysis = nullptr;
      if (const ContinuousDistribution* d = distr. asContinuousDistribution ())
        if (   eqReal (d->getLoBoundEffective (), getLoBoundEffective ())
            && eqReal (d->getHiBoundEffective (), getHiBoundEffective ())
           )
          analysis = d->analysis;
    }
  void variable2analysis () const final
    { checkPtr (analysis) -> variable = variable; }
  void analysis2variable () const final
    { variable = checkPtr (analysis) -> variable; }
  // Parameters
  Real pdfVariable () const final
    { return pdf (variable); }
  Real logPdfVariable () const final
    { return logPdf (variable); }
  void randVariable () const final
    { variable = rand (); }
  Real getInfoMean () const final
    { return stdBounds () ? getInfoMean_ () : NAN; }
  Real getInfoVar () const final
    { return stdBounds () ? getInfoVar_ () : NAN; }

  // Unbounded support
  virtual Real getQuantile (Prob /*p*/) const 
    { NOT_IMPLEMENTED; return NAN; }
    // Return: inverse cdf()
  virtual Real getInfoMean_ () const
    { return NAN; }
  virtual Real getInfoVar_ () const
    { return NAN; }
};



struct LocScaleDistribution : ContinuousDistribution
{
  // Parameters
  An::Value loc {NAN};
    // Location
    // If stdLoBound() = 0 then 0
  An::Value scale {NAN};
    // > 0
//Lower bound on scale ??
  // <Others>

  // Functions of parameters
protected:
  Real log_scale {NAN};
  // <Others>
public:
    

protected:
  explicit LocScaleDistribution (const string &name_arg)
    : ContinuousDistribution (name_arg) 
    {}
public:
  void qc () const override;


  const LocScaleDistribution* asLocScaleDistribution () const final
    { return this; }

protected:
  void finish () override
    { UniDistribution::finish ();
      log_scale = log (scale);
    }
public:
  bool getParamSet () const override
    { return    ! isNan (loc)
             && ! isNan (scale);
    }
  size_t paramCount () const override
    { return 2; }   
  bool similar (const Distribution &distr,
                Real delta) const override
    { if (! ContinuousDistribution::similar (distr, delta))
        return false;
      if (const LocScaleDistribution* cd = distr. asLocScaleDistribution ())
        return    sameReal (loc,   cd->loc,   delta)
               && sameReal (scale, cd->scale, delta);
      return false;
    }
  // UniDistribution
  string nameParam () const override
    { return  name + "(" + real2str (loc, 3) + "," + real2str (scale, 3) + ")"; }
private:
  Real pdf_ (Real x) const final
    { return pdfStnd (stnd (x)) / scale; }
    // >= 0
  Real logPdf_ (Real x) const final
    { return logPdfStnd (stnd (x)) - log_scale; }
public:

  Real stnd (Real x) const
    { return (x - loc) / scale; }
  Real unstnd (Real x) const
    { return x * scale + loc; }
protected:
  // loc = 0, scale = 1
  virtual Real pdfStnd (Real xStnd) const
    { return exp (logPdfStnd (xStnd)); }
    // >= 0
  virtual Real logPdfStnd (Real xStnd) const
    { return log (pdfStnd (xStnd)); }

  // Unbounded support
  Real getInfoMean_ () const final
    { return getInfoMeanStnd () + log_scale; }
  virtual Real getInfoMeanStnd () const 
    { return NAN; }
public:
    
  virtual void setMeanVar (Real /*mean*/,
                           Real /*var*/) 
    { NOT_IMPLEMENTED; }
};



struct Normal : LocScaleDistribution 
// getEntropy_est() = 0.5
{
  static const Real coeff;
    // = log(2*pi)
  
  
  Normal ()
    : LocScaleDistribution ("Normal")
    {}
  Normal* copy () const final
    { return new Normal (*this); }
  void qc () const final;


  const Normal* asNormal () const final
    { return this; }

  // Parameters
  void setParam (Real loc_arg,
                 Real scale_arg)
    { loc   = loc_arg;
      scale = scale_arg;
      finish ();
    }
  void estimate () final;
  // UniDistribution
private:
  Prob cdf_ (Real x) const final
    { return (erf (stnd (x) / sqrt_2) + 1) / 2; }
  Real rand_ () const final;
public:
  Real getMean () const final
    { return loc; }
  Real getVar () const final
    { return sqr (scale); }
  // ContinuousDistribution
private:
  Real logPdfStnd (Real xStnd) const final
    { return coeff == INF 
               ? nullReal (xStnd)
                 ? 0  // PAR
                 : -INF
               : - 0.5 * (coeff + sqr (xStnd)); 
    }
public:
  Real getQuantile (Prob p) const final;
private:
  Real getInfoVar_ () const final
    { return 0.5; }
  Real getInfoMeanStnd () const final
    { return 0.5 * (1 + coeff); }
public:
  
  void setMeanVar (Real mean,
                   Real var) final
    { setParam (mean, sqrt (var)); }
  Prob pValue_2tail (Real x) const
    { const Real diff = abs (getMean () - x);
      if (nullReal (diff))
        return 1;
      if (nullReal (scale))
        return 0;
      return 2 * cdf (getMean () - diff); 
    }
};



struct Exponential : LocScaleDistribution 
{
  Exponential ()
    : LocScaleDistribution ("Exponential")
    {}
  Exponential* copy () const final
    { return new Exponential (*this); }
  void qc () const final;


  const Exponential* asExponential () const final
    { return this; }

  // Parameters
  void setParam (Real loc_arg)
    { loc = loc_arg;
      scale = loc;
      finish ();
    }
  void estimate () final;
  // UniDistribution
private:
  Prob cdf_ (Real x) const final
    { return 1 - exp (- x / loc); }
  Real rand_ () const final
    { return - loc * log (1 - randProb ()); }
public:
  Real getMean () const final
    { return loc; }
  Real getVar () const final
    { return sqr (scale); }
  // ContinuousDistribution
private:
  Real logPdfStnd (Real xStnd) const final
    { return - xStnd; }
  Real getInfoMeanStnd () const final
    { return 1; }
  Real getInfoVar_ () const final
    { return 1; } 
  string nameParam () const final
    { return  name + "(" + real2str (loc, 3) + ")"; }
  Real stdLoBound () const final
    { return 0; }
public:

  void setMeanVar (Real mean,
                   Real /*var*/) final
    { setParam (mean); }
};



struct Cauchy : LocScaleDistribution 
// Cauchy(0,1) ~ Normal(0,1) / Normal(0,1)
// Cauchy(loc,scale) ~ Student_1(loc,scale)
// Cauchy(0,1) ~ tan(pi*(Uniform(0,1) - 0.5))
{
  Cauchy ()
    : LocScaleDistribution ("Cauchy")
    {}
  Cauchy* copy () const final
    { return new Cauchy (*this); }
  void qc () const final;


  const Cauchy* asCauchy () const final
    { return this; }

  // Parameters
  void setParam (Real loc_arg,
                 Real scale_arg)
    { loc   = loc_arg;
      scale = scale_arg;
      finish ();
    }
  void estimate () final;
  Real getSortingValue () const final
    { return loc; }
  // UniDistribution
private:
  Prob cdf_ (Real x) const final
    { return 0.5 + atan (stnd (x)) / pi; }
  Real rand_ () const final
    { return unstnd (tan (pi * (randProb () - 0.5))); }
public:
  Real getMean () const final
    { return NAN; }
  Real getVar () const final
    { return NAN; }
  // ContinuousDistribution
private:
  Real pdfStnd (Real xStnd) const final
    { return 1 / (pi * (1 + sqr (xStnd))); }
  Real getInfoMeanStnd () const final
    { return log (4 * pi); }
  // Solve: \int \log (x^2 + a) / (x^2 + a) dx ??
//Real getInfoVar_ () const
  //{ return NAN; }    
  // Solve: \int \log^2 (x^2 + a) / (x^2 + a) dx ??
};



struct Beta1 : ContinuousDistribution 
// = Beta(alpha,1)
// loc, scale: use ??
{  
  // Parameters
  Real alpha {NAN};
    // 0..1 (Prob ??)
  // Functions of parameters
private:
  Real logAlpha {NAN};
public:


  Beta1 ()
    : ContinuousDistribution ("Beta")
    {}
  Beta1* copy () const final
    { return new Beta1 (*this); }
  void qc () const final;


  const Beta1* asBeta1 () const final
    { return this; }

  // Parameters
  void setParam (Real alpha_arg)
    { alpha = alpha_arg;
      finish ();
    }
private:
  void setParamFunc ()
    { logAlpha = log (alpha); }
public:
  void estimate () final;
  bool getParamSet () const 
    { return ! isNan (alpha); }
  bool similar (const Distribution &distr,
                Real delta) const
    { if (! ContinuousDistribution::similar (distr, delta))
        return false;
      if (const Beta1* beta1 = distr. asBeta1 ())
        return eqReal (alpha, beta1->alpha, delta);
      return false;
    }
  Real getInfoMean_ () const final
    { return logAlpha + (1 - alpha) / alpha; }
  // UniDistribution
  string nameParam () const
    { return  name + "(" + real2str (alpha) + ",1)"; }
  Real stdLoBound () const final
    { return 0; }
  Real stdHiBound () const final
    { return 1; }
private:
  Real pdf_ (Real x) const
    { return alpha * pow (x, alpha - 1); }
  Prob cdf_ (Real x) const
    { return pow (x, alpha); }
  Real rand_ () const
    { return pow (randProb (), 1 / alpha); }
    // Return: isProb()
public:
  Real getMean () const
    { return alpha / (alpha + 1); }
  Real getVar () const
    { return alpha / (sqr (alpha + 1) * (alpha + 2)); }
  // ContinuousDistribution
private:
  Real pdfStnd (Real xStnd) const
    { return pdf_ (xStnd); }
};



struct UniKernel : ContinuousDistribution 
// loc, scale: use ??
// MultiVariate: use for clustering (single-linkage using halfWindow), cf. dbscan ??
{  
  // Parameters
private:
  struct Point : Root
  { Real value;
    Real mult;
    Point (Real value_arg,
           Real mult_arg)
      : value (value_arg)
      , mult (mult_arg)
      { qc (); }
    explicit Point (Real value_arg)
      : value (value_arg)
      , mult (NAN)
      {}
    Point ()
      : value (NAN)
      , mult (NAN)
      {}
    void qc () const;
    bool operator== (const Point &other) const
      { return value == other. value; }
    bool operator< (const Point &other) const
      { return value < other. value; }
  };
  Vector<Point> points;
    // Sorted ascending
  Vector<Real> multSum;
public:
  Real attr_min {NAN};
  Real attr_max {NAN};
  Prob uniform_prob {NAN};
  Real halfWindow {NAN};
    // >= 0


  explicit UniKernel (const UniVariate<NumAttr1> &analysis_arg)
    : ContinuousDistribution ("UniKernel")
    { analysis = & analysis_arg; }
  UniKernel* copy () const final
    { return new UniKernel (*this); }
  void qc () const final;


  const UniKernel* asUniKernel () const final
    { return this; }

  // Parameters
  void setParam (Real window)
    { halfWindow = window / 2;
      setPoints ();
      set_uniform_prob ();
      finish ();
    }
  void estimate () final;
    // Invokes: pdf()
private:
  Real setPoints ();
    // Return: SD
  void estimate_ (Real halfWindow_lo,
                  Real halfWindow_hi,
                  Real step);
    // Output: uniform_prob, halfWindow
  void set_uniform_prob ();
    // Input: halfWindow
    // Output: uniform_prob
  size_t findIndex (Point x) const;
    // Return: max{index : x >= points[index].value}; may be NO_INDEX
    // Time: O(log(analysis->sample->nEffective))
  Real getHeight () const
    { return 1 / (checkPtr (analysis) -> sample. mult_sum * 2 * halfWindow); }
    // Kernel
    // Integral over (2 * halfWindow) = 1 / analysis->sample. mult_sum
  Real getUniformHeight () const
    { return 1 / getRange (); }
    // Integral over getRange() = 1
public:
  bool getParamSet () const
    { return isProb (uniform_prob) && halfWindow >= 0; }
  bool similar (const Distribution &/*distr*/,
                Real /*delta*/) const
    { return false; }
  // UniDistribution
  string nameParam () const
    { return name + " (" + prob2str (uniform_prob) + "," + real2str (2 * halfWindow, getAttr (). decimals + 1) + ")"; }  
private:
  Real pdf_ (Real x) const;
    // Invokes: findIndex()
    // halfWindow = 0 => INF or 0
  Prob cdf_ (Real x) const;
  Real rand_ () const
    { return points [(size_t) round (randProb () * ((Real) points. size () - 1))]. value; }
    // Bootstrap
public:
  Real getMean () const;
  Real getVar () const;
  // ContinuousDistribution
private:
  Real pdfStnd (Real xStnd) const
    { return pdf_ (xStnd); }

public:
  const NumAttr1& getAttr () const
    { return checkPtr (analysis) -> attr; }
  Real getRange () const
    { return attr_max - attr_min; }
};

  

struct MultiDistribution : Distribution
{
  typedef  MultiVariate<NumAttr1>  An;
  const An* analysis {nullptr};
  mutable An::Value variable;
  mutable MVector x_field;


protected:
  MultiDistribution (const string &name_arg)
    : Distribution (name_arg)
    {}
public:
  void qc () const override;


  const MultiDistribution* asMultiDistribution () const final
    { return this; }    

  const Analysis1* getAnalysis () const 
    { return analysis; }
  Analysis1* createAnalysis (Dataset &ds);
  void removeAnalysis () 
    { analysis = nullptr; }
  void shareAnalysis (const Distribution &distr) 
    { analysis = nullptr;
      if (const MultiDistribution* d = distr. asMultiDistribution ())
        if (getDim () == d->getDim ())
          analysis = d->analysis;
    }
  void variable2analysis () const
    { checkPtr (analysis) -> variable = variable; }
  void analysis2variable () const
    { variable = checkPtr (analysis) -> variable; }
  // Parameters
//setDim(); set <Parameters>; setParam();
  void setDim (size_t dim)
    { variable. resize (dim); 
      x_field.  resize (dim); 
    }
  size_t getDim () const 
    { return variable. size (); }
  bool similar (const Distribution &distr,
                Real /*delta*/) const
    { if (const MultiDistribution* d = distr. asMultiDistribution ())
        return getDim () == d->getDim ();
      return false;
    }
  Real getSortingValue () const 
    { MVector v (getDim ());
      getMeanVec (v);
      return v. sumSqr ();
    }

  // Input: <Parameters>
  virtual void getMeanVec (MVector &mean) const = 0;
  virtual void getVC (Matrix &vc) const = 0;
    
};



struct MultiNormal : MultiDistribution
{
  // Parameters
  MVector mu;
  Matrix sigmaExact;
    // VC matrix
    // size() = (getDim(), getDim())
    // det >= 0 
  Matrix sigmaInflated;
    // VC matrix
    // size() = (getDim(), getDim())
    // det >= 0 
    // Inflated if variance_min > 0, otherwise = sigmaExact
  MVector variance_min;
    // size() = mu.size()
    // Init: 0
    // Is not estimate()'ed
    
  // Functions of parameters
  Matrix sigmaInv;
  Real coeff {NAN};
    // May be NAN
    // = 0.5 * (getDim() * Normal::coeff + log(det(sigmaInflated)))
private: 
  Common_sp::AutoPtr <Matrix> cholesky;
  // Random numbers generator
  Vector<Normal> zs;
    // size() == getDim()
public: 


  MultiNormal ()
    : MultiDistribution ("Normal")
    { /*sigma. psd = true;*/ }
  MultiNormal* copy () const final
    { return new MultiNormal (*this); }
  void qc () const override;
  void saveText (ostream& os) const override;


  const MultiNormal* asMultiNormal () const final
    { return this; }    

  void setSeed (ulong seed) const;
  // Parameters
  void setDim (size_t dim);
    // Output: mu = 0, sigmaInflated = 0, sigmaExact = 0
  void setParam ()
    // Input: mu, sigmaInflated, sigmaExact
    { finish (); }
private:
  void setParamFunc ();
    // Invokes: inflateSigma()
public:
  void estimate () final;
  bool getParamSet () const
    { return    getDim () 
             && mu. defined () 
             && sigmaExact. defined () 
             && sigmaInflated. defined () 
             && ! isNan (coeff)
             ; 
    }
  size_t paramCount () const
    { return getDim () + (getDim () * (getDim () + 1)) / 2; }
  bool similar (const Distribution &distr,
                Real delta) const
    { if (! MultiDistribution::similar (distr, delta))
        return false;
      if (const MultiNormal* mn = distr. asMultiNormal ())
        return    mu.            maxAbsDiff (false, mn->mu,            false) <= delta
               && sigmaInflated. maxAbsDiff (false, mn->sigmaInflated, false) <= delta;
      return false;
    }
  Real logPdfVariable () const
    { if (isNan (coeff))
        throw runtime_error ("MultiNormal::coeff is NAN"); 
      x_field = variable;
      if (coeff == INF)
        return nullReal (x_field. maxAbsDiff (false, mu, false)) ? 0 : -INF;  // PAR
      const Real mah = sigmaInv. getMahalanobis ( false
                                                , x_field, true, 0
                                                , mu,      true, 0
                                                );
      return -0.5 * mah - coeff;
    }
  void randVariable () const;
    // variable = mu + cholesky * zs
  Real getInfoMean () const final
    { if (isNan (coeff))
        throw runtime_error ("MultiNormal::coeff is NAN"); 
      return 0.5 * (Real) getDim () + coeff; 
    }
  Real getInfoVar () const final
    { return (Real) getDim () / 2.0; }    
  // MultiDistribution
  void getMeanVec (MVector &mean) const 
    { mean = mu; }
  void getVC (Matrix &vc) const
    { vc = sigmaExact; }
  
  void getSigmaRaw (Matrix &sigmaRaw) const;
    // Output: sigmaRaw: sigmaExact + mu * mu'
private:
  bool inflateSigma ();
    // Make the variance in the space of PC >= variance_min
    // Input: sigmaExact
    // Output: sigmaInflated
};



struct Mixture : Distribution  
{
  // Parameters
  struct Component : Root
  {
    // Parameters
    Common_sp::AutoPtr<Distribution> distr;
      // !nullptr
      // getParamSet()
    Prob prob;
      // Functions of parameters
    Vector<Prob> objProb;
      // Size = Analysis1::samlpe.size()
      // Requires: (bool)distr->getAnalysisCheck()
    
  private:
    Component (Distribution* distr_arg,
               Prob prob_arg)
      : distr (distr_arg)
      , prob (prob_arg)
      {}
    friend struct Mixture;
  public:
    Component* copy () const final
      { return new Component (*this); }
    void qc () const override;
    void saveText (ostream& os) const override
      { distr->saveText (os); 
        os << " P=" << prob << endl;
      }
      
    // Input: Parameters
    Real pdfProb () const
      { return prob ? prob * distr->pdfVariable () : 0; }
    Real logPdfProb () const
      { return log (prob) + distr->logPdfVariable (); }
    Real getMult () const;
    Real getDeltaness () const;
      // = getMult() - max_obj obj->mult
      // Return: >= 0
      //         0 <=> delta-function
    void estimate ();
    Real getFitness_entropy () const;
    void merge (const Component* comp);     
      // Update: prob, objProb[]
  protected:
    struct Component2Sample : SaveSample
    {
      Component2Sample (Sample &sample_arg,
                        const Component &comp)
        // Update: sample_arg
        : SaveSample (sample_arg)  
        { for (Iterator it (sample_arg); it ();)  
            sample_arg. mult [*it] *= comp. objProb [*it];
          sample_arg. finish ();
        }
    };
  public:
  #if 0
    const MultiNormal* getMultiNormal () const
      { return distr->asMultiNormal (); }
      // Return: may be nullptr
  #endif
  };
  VectorOwn<Component> components;
    // Same: Distribution::{<analysis>,<variable>}
    // sum_j(components[j]->prob) = 1
    // For each i: sum_j(components[j]->objProb[i]) = 1
    // All are continuous or all are discrete 
    // sort'ed by Distribution::getSortingValue()
    
  // Functions of parameters
protected:
  Categorical cat;
public:


  Mixture ()
    : Distribution ("Mixture")
    {}
  Mixture* copy () const final
    { return new Mixture (*this); }
  void qc () const override;
  void saveText (ostream& os) const override;


  const Mixture* asMixture () const final
    { return this; }

  const Analysis1* getAnalysis () const final
    { return components. empty () ? nullptr : components [0] -> distr->getAnalysis (); }
  Analysis1* createAnalysis (Dataset &ds) final
    { if (components. empty ())
        return nullptr;
      return const_cast <Component*> (components [0]) -> distr->createAnalysis (ds);
    }
  void removeAnalysis () final
    { for (const Component* comp : components)
        const_cast <Component*> (comp) -> distr->removeAnalysis ();
    }
  void shareAnalysis (const Distribution &distr) final
    { for (const Component* comp : components)
      { const_cast <Component*> (comp) -> distr->shareAnalysis (distr);
        if (! comp->distr->getAnalysis ())
        { removeAnalysis ();
          break;
        }
      }
    }
  void variable2analysis () const final
    { if (! components. empty ())
        components [0] -> distr->variable2analysis ();
    }
  void analysis2variable () const final
    { for (const Component* comp : components)
        comp->distr->analysis2variable ();
    }
  void setSeed (ulong seed) const override;
  // Parameters
  Component* addComponent (Distribution* distr,
                           Prob prob = NAN);
    // Return: !nullptr
    // Append: components (may be re-ordered later)
//distr->setParam()
  void setParam ()
    { finish (); }
private:
  void setParamFunc ();
public:   
  void estimate () final;
    // Input: Component::{parameters, prob} or Component::objProb
    // Invokes: components.sort()
  bool getParamSet () const;
    // Input: Component::{parameters, prob}
  size_t getDim () const 
    { return components. empty () ? SIZE_MAX : components [0] -> distr->getDim (); }
  size_t paramCount () const;
  bool similar (const Distribution &distr,
                Real delta) const;
  Real getSortingValue () const
    { Real s = 0;
      for (const Component* comp : components)
        s += comp->prob * comp->distr->getSortingValue ();
      return s;
    }
  Real pdfVariable () const;
    // Return: sum(Component::pdfProb(x))
  Real logPdfVariable () const;
  void randVariable () const
    { cat. randVariable ();
      Distribution* distr = components [cat. variable] -> distr. get ();
      distr->randVariable (); 
      variable2analysis ();
      analysis2variable ();
    }

//Prob getFitness_entropy_min () const;
    // Over components
  Prob getConfusion () const;
    // Return: 1 - average_obj max_cluster P(obj in cluster)
  Prob getOverlap (size_t compNum1,
                   size_t compNum2) const;
    // Return: <= P(X is assigned to a wrong cluster): <= 0.5

private:
  void balanceProb ();
    // Update: Component::prob
public:
  void mergeComponents ();
    // Merge Component's with identical Distribution's
  bool deleteComponent (size_t num);
    // Return: success    
};




////////////////////////////////// Construction ///////////////////////////////

struct PrinComp : MultiVariate<NumAttr1>
// Principal components
{
  typedef  MultiVariate<NumAttr1>  P;
  
  MultiNormal mn;  
  Eigens eigens;
    // len = mn.getDim()


  PrinComp (const Sample &sample_arg,
            const Space1<NumAttr1> &space_arg,
            const MultiNormal &mn_arg,
            size_t maxCount,
            Prob totalExplainedFrac_max,
            Prob explainedFrac_min,
            Real error)
    : P (sample_arg, space_arg)
    , mn (mn_arg)
    , eigens ( mn. sigmaExact
             , maxCount
             , totalExplainedFrac_max
             , explainedFrac_min
             , error
             , 5000  // PAR
             )
    { mn. analysis = nullptr; }
  PrinComp* copy () const final
    { return new PrinComp (*this); }
  void qc () const override;
  void saveText (ostream &os) const override
    { if (verbose ())
        mn. saveText (os);  
      eigens. saveText (os);
    }


  size_t getOutDim () const
    { return eigens. getDim (); }
  Real getQuality () const
    { return eigens. totalExplainedFrac (); }
  Real getVariance (size_t attrNum) const
    { return eigens. values [attrNum]; }
  void project (const P::Value &x,
                MVector &projection) const;
    // Input: x: size() = space.size()
    // Output: projection: size() = getOutDim()
  void project (size_t objNum,
                MVector &projection) const
    { data2variable (objNum);
      project (variable, projection);
    }
  Space1<RealAttr1> createSpace (const string &attrPrefix,
                                 Dataset &ds) const;
    // Sum d^2 of createSpace() is equal to sum d^2 of the original space
    // Invokes: project()
  Real getAttrMds (size_t attrNum,
                   size_t eigenNum) const
    { return eigens. getMds (attrNum, eigenNum, false) / sqrt_2; }
  Dataset createAttrMds (const string &attrPrefix,
                         Vector<Real> &quality) const;
    // Return RealAttr1->name = attrPrefix + <i>
    // Output: quality; size() = Return->attrs.size()
    // Invokes: getAttrMds(), new RealAttr1
  Real getChi2 (size_t objNum) const;
    // Return: squared standardized principal components
    //         ~ Chi2(getOutDim()) if data has multinormal distribution
};



struct Clustering : MultiVariate<NumAttr1>
{
  typedef  MultiVariate<NumAttr1>  P;
  
  Mixture mixt;
    // Component::distr->asMultiNormal()
    // components[i] - cluster i
    // May be: !getParamSet()
  MVector variance_min;  
    // Of each component of MultiNormal
    // Function of data and {sd_min,sd_min_is_relative} of constructor
  
  
  Clustering (const Sample &sample_arg,
              const Space1<NumAttr1> &space_arg,
              size_t clusters_max,
              Real sd_min,
              bool sd_min_is_relative,
              Real entropyDimensionPrecision = 0.001);  // PAR
    // Input: sd_min: (relative) min. SD for each attribute in each cluster
  Clustering (const Clustering &clustering,
              const vector<bool> &toMerge);
    // Requires: clustering.getOutDim() = toMerge.size()
private:
  void setAnalysis ();
    // Invokes: = this
  Mixture::Component* addCluster (Prob p)
    { MultiNormal* distr = new MultiNormal ();
      distr->setDim (space. size ());
      distr->variance_min = variance_min;
      distr->analysis = this;
      return mixt. addComponent (distr, p);
    }
  void splitCluster (size_t num);
    // Invokes: Mixture::addComponent(); PrinComp
  bool deleteCluster (size_t num)
    { return mixt. deleteComponent (num); }
public:
  Clustering* copy () const final
    { return new Clustering (*this); }
  void qc () const override;
  void saveText (ostream &os) const override;


  size_t getOutDim () const
    { return mixt. components. size (); }
    // Return: >= 1
  const MultiNormal* getMultiNormal (size_t i) const
    { return mixt. components [i] -> distr->asMultiNormal (); }
    // Return !nullptr
    
  // Return !nullptr
  // Update: mixt.space->ds
  Space1<ProbAttr1> createSpace (Dataset &ds) const;
    // new ProbAttr1() for each cluster
  NominAttr1* createNominAttr (const string &attrName,
                               Prob prob_min,
                               Dataset &ds) const;
    // Input: prob_min: min. class membership probability to produce a non-missing value
  ProbAttr1* createProbAttr (const string &attrName,
                             streamsize decimals,
                             Dataset &ds) const;

  bool mergeClose (NominAttr1 &nominAttr,
                   ProbAttr1 &probAttr,
                   Prob confused_max) const;
    // Does nothing if all clusters overlap
    // Return: merge is done
    // Input: confused_max: max. allowed Mixture::getOverlap()
    // Requires: nominAttr = *createNominAttr()
    //           probAttr  = *createProbAttr()
    // Invokes: mixt.getOverlap(), DisjointCluster::merge()
  void merge (size_t compNum1,
              size_t compNum2);
    // Invokes: mixt.components.erasePtr(compNum2)
  void processSubclusters (const string &clusterAttrName,
                           bool merge_close,
                           const VectorPtr<Attr> &attrs_orig,
                           JsonArray* jClusts,
                           Prob attr_pvalue,
                           const string &outDir,
                           const string &outGenericDm,
                           Dataset &ds,
                           Space1<Attr1> &sp1) const;
    // Create files <outDir>/<outGenericDm>.<N><dmSuff>
    // Input: jClusts: may be nullptr
    //        outDir: "/dev/null" <=> no file is created
    // Update: sp1 (append)
    // Invokes: createNominAttr(clusterAttrName), createProbAttr(clusterAttrName + "_prob"), mergeClose()
};



#if 0
struct OneCluster : MultiVariate<NumAttr1>
{
  typedef  MultiVariate<NumAttr1>  P;
  
  Mixture mixt;
    // Of MultiNormal and MultiUniform (= "outliers")
};
#endif



struct Canonical : MultiVariate<NumAttr1>
{
  typedef  MultiVariate<NumAttr1>  P;
  
  Notype requires;
  MultiNormal between;  
  MultiNormal within;  
  Common_sp::AutoPtr<const Matrix> choleskyInv;
    // May be nullptr ??
  Matrix basis;
  Matrix basis_norm;
  MVector eigenValues;
    // size() <= clustering. getOutDim () - 1
  

  explicit Canonical (const Clustering &clustering);
private:
  static Notype getRequirement (const Clustering &clustering)
    { if (clustering. getOutDim () <= 1)
        throw logic_error ("Canonical: Too few clusters");
      return Notype ();
    }
  static MultiNormal getBetween (const Clustering &clustering);
  static MultiNormal getWithin (const Clustering &clustering);
  static Matrix* getCholeskyInv (const Matrix &withinSigma);
public:
  Canonical* copy () const final
    { return new Canonical (*this); }
  void qc () const override;
  void saveText (ostream &os) const override;


  size_t getOutDim () const
    { return eigenValues. size (); }
    // Return: may be 0
  // Requires: getOutDim()
  void project (const P::Value &x,
                MVector &projection) const;
    // Input: x: size() = space.size()
    // Output: projection
    // Requires: projection.size() = getOutDim()
  void project (size_t objNum,
                MVector &projection) const
    { data2variable (objNum);
      project (variable, projection);
    }
  Space1<RealAttr1> createSpace (const string &attrPrefix,
                                 Dataset &ds) const;
};



struct Mds : Analysis 
// Linear multi-dimensional scaling
{
  Eigens eigens;
  
  
  Mds (const Sample &sample_arg,
       const RealAttr2 &attr2,
       size_t outDim_max,
       Prob totalExplainedFrac_max,
       Prob explainedFrac_min);
    // Requires: attr2.matr: defined(), symmetric, centered
    //           space.ds.getUnitMult() ??
    // Time: O(n^2 outDim_max)
  Mds* copy () const final
    { return new Mds (*this); }
  void qc () const final;
  void saveText (ostream &os) const override
    { eigens. saveText (os); }


  size_t getOutDim () const
    { return eigens. getDim (); }
  Real getQuality () const
    { return eigens. totalExplainedFrac (); }
  Real getVariance (size_t attrNum) const
    { return eigens. values [attrNum] / sample. mult_sum; }
  // attrNum < eigens.values.size()
  Real get (size_t objNum,
            size_t attrNum,
            bool imaginary) const
    { return eigens. getMds (objNum, attrNum, imaginary); }
  Real project (const Matrix &similarity,
                bool t,
                size_t row,
                size_t attrNum) const
    { return multiplyVec (eigens. basis, true, attrNum, similarity, t, row) / sqrt (eigens. values [attrNum]); }
  Space1<RealAttr1> createSpace (const string &realAttrPrefix,
                                 const string &imaginaryAttrPrefix,
                                 Dataset &ds) const;
    // Update: ds
};



}  



#endif
