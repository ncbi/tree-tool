// prediction.hpp 
 
#ifndef PREDICTNUM_HPP_78455
#define PREDICTNUM_HPP_78455

#include "../common.hpp"
using namespace Common_sp;
#include "numeric.hpp"
#include "matrix.hpp"
#include "dataset.hpp"
using namespace DM_sp;

#define LIN_PROG 0
#if LIN_PROG
  #include "../C/CPLEX.hpp"  
#endif



namespace DM_sp
{



template <typename Pred/*Attr1*/, typename Targ/*Attr1*/>
  struct Prediction : MultiVariate<Pred>
  {
    typedef  MultiVariate<Pred>  P;
    typedef  Pred  Predictor;
    typedef  Targ  Target;
    
    const Target& target;
    mutable typename Target::Value targetVariable;
  
  
  protected:
    Prediction (const Sample &sample_arg,
                const Space1<Predictor> &space_arg,
                const Target& target_arg)
      : P (sample_arg, space_arg)
      , target (target_arg)
      , targetVariable ()
      {}
      // To invoke: resize()
  public:
    void qc () const override
      { if (! qc_on)
          return;
        P::qc ();
        target. qc ();
        ASSERT (& target. ds == P::sample. ds);
      }
  
  
    void data2variable (size_t objNum) const override
      { P::data2variable (objNum);
        targetVariable = target [objNum];
      }
    void variable2data (size_t objNum) override
      { P::variable2data (objNum);
        const_cast <Target&> (target) [objNum] = targetVariable;
      }
      
    bool isOverfit () const
      { return P::space. size () > P::space. ds. objs. size (); }
    size_t getMaxPredictorNameLen () const;
    bool getExistsMissing () const
      { size_t badObjNum;
      	return P::space. existsMissing () || target. existsMissing (badObjNum); 
      }
  };




struct LogisticRegression : Prediction<NumAttr1,BoolAttr1>
{
  // Output
  Real negLogLikelihood_ave {NAN};
    // !isNan() <=> solve() is successul
    // Init: NAN
    // >= 0
  MVector beta;
    // rowsSize(true) = space.size()
    // Precision ??
    // getSeparated(), beta *= k, k > 1 => negLogLikelihood_ave is improved
  Real target_score_min {NAN};
  Real non_target_score_max {NAN};
  MVector attrImportance;
    // = negLogLikelihood_ave(without attr) / negLogLikelihood_ave(with attr) - 1
    // May be < 0 if separated()
    // rowsSize(true) = space.size()
    // at(i) = relative change of negLogLikelihood_ave if space[i] is erase()'d

    
  LogisticRegression (const Sample &sample_arg,
                      const Space1<NumAttr1> &space_arg,
                      const BoolAttr1 &target_arg)
    : Prediction<NumAttr1,BoolAttr1> (sample_arg, space_arg, target_arg)
    { resize (); }
  LogisticRegression* copy () const final
    { return new LogisticRegression (*this); }
  void qc () const override;
    // Requires: space[0] = ds.addNumAttrUnit()
  void print (ostream &os) const override
    { os << "-logL = " << negLogLikelihood_ave << endl;
      os << "P(Correct) = " << getCorrectPredictionFrac () << endl;
      os << "Target_score_min = " << target_score_min << endl;
      os << "Other_score_max  = " << non_target_score_max << endl;
      os << "------------------------------" << endl;
      os << "beta: value" << endl;
      FOR (size_t, i, space. size ())
        os << space [i] -> name << ": " << beta [i] << endl;
    }


  void resize () 
    { beta.           resize (space. size ());
      attrImportance. resize (space. size ());
      negLogLikelihood_ave = NAN;
      target_score_min     = NAN;
      non_target_score_max = NAN;
    }

  Real getScore (size_t objNum) const;
    // Input: beta
  Prob predict (size_t objNum) const
    { return probit (getScore (objNum)); }
    // Prediction of target
  Real getNegLogLikelihood_ave () const;
    // Estimator of the conditional entropy
    // Invokes: getScore()
  Prob getCorrectPredictionFrac () const
    { return exp (- negLogLikelihood_ave); }
  RealAttr1* getScoreAttr (const string &targetNameSuffix)
    { auto a = new RealAttr1 (target. name + "_" + targetNameSuffix, const_cast <Dataset &> (space. ds));  
      for (Iterator it (sample); it ();)  
        (*a) [*it] = getScore (*it);
      return a;
    }
  ProbAttr1* getPredictionAttr (const string &targetNameSuffix)
    { auto a = new ProbAttr1 (target. name + "_" + targetNameSuffix, const_cast <Dataset &> (space. ds));  
      for (Iterator it (sample); it ();)  
        (*a) [*it] = predict (*it);
      return a;
    }
  bool getSeparated () const
    { return lessReal (non_target_score_max, target_score_min); }
  void simulateTarget (ulong seed);
    // Output: target[]
  void solve ();
    // Output: negLogLikelihood_ave, beta
  void test (ulong seed);
    // Input: beta
    // Invokes: solve()
  void setAttrImportance ();
    // Input: negLogLikelihood_ave: !isNan()
    // Output: attrImportance
    // Time: p * time(solve)
    // Invokes: solve()
    // --> Prediction ??
};



template <typename T /*:Prediction*/>
T* selectAttrs (const Sample &sample,
                Space1<typename T::Predictor> &space,
                const typename T::Target &target,
                Real attrImportance_min)
// Return: final T
// Update: space
// Time: O(p * time(setAttrImportance))
{ Progress prog ((uint) space. size ());
  for (;;)
  { Common_sp::AutoPtr<T> t (new T (sample, space, target));
    t->solve ();
    t->qc ();
    if (verbose ())
      cerr << t->negLogLikelihood_ave << endl;  
    t->setAttrImportance ();
    const typename T::Predictor* attr_bad = nullptr;
    size_t i_bad = NO_INDEX;
    Real importance_min = INF;
    FOR (size_t, i, t->space. size ())
      if (minimize (importance_min, t->attrImportance [i]))
      { attr_bad = t->space. at (i);
        i_bad = i;
      }
    if (importance_min >= attrImportance_min)
    { t->solve ();
      return t. release ();
    }
    prog (attr_bad->name);
    if (verbose ())
      cerr << attr_bad->name << ": " << importance_min << endl;  
    space. eraseAt (i_bad);
  }
  return nullptr;
}




struct LinearNumPrediction; 



// Constraint

struct Constraint : Root
// Linear constraint on LinearNumPrediction::beta
{
  const LinearNumPrediction& pred;
  MVector lhs;
  Real rhs;
  char sense;
    // 'L': <=
    // 'E': ==
    // 'G': >=

  Constraint (const LinearNumPrediction &pred_arg);
};




//////////////////// LinearNumPrediction /////////////////////

struct LinearNumPrediction : Prediction<NumAttr1,RealAttr1>
// target = sum_i beta_i * x_i
{
  typedef  Prediction<NumAttr1,RealAttr1>  P;
  
  // INPUT
  VectorOwn<Constraint> constraints;  

  // OUTPUT
  Real absCriterion {NAN};
    // >= 0, -> min
    // May be NAN
    // Init: NAN
  MVector beta;
    // size = predictors_.size()
    // Init: NAN
  MVector betaSD;
    // Uncertainty of beta
    // Disregards constraints
    // size = predictors_.size()
	  // Init: NAN
	  
#if LIN_PROG
protected:
  Commpn_sp::AutoPtr<CPLEX> Problem;
    // Init: nullptr
    // Requries: beta [] are the first variables
public:
#endif


protected:
  LinearNumPrediction (const Sample &sample_arg,
                       const Space1<NumAttr1> &space_arg,
                       const RealAttr1 &target_arg)
    : P (sample_arg, space_arg, target_arg)
    , absCriterion (NAN)
  #if LIN_PROG
    , Problem (nullptr)
  #endif
    {}
public:
  void qc () const override;
  void saveText (ostream &os) const override;


  void resize () 
    { beta.   resize (space. size ());
      betaSD. resize (space. size ());
    }

  virtual void solveUnconstrained () = 0;
    // Output: absCriterion: may be isNan()
    //         beta: if isNan(absCriterion) then unchanged

#if LIN_PROG
  // Problem: solve with constraints
  void setProblemVerbose (bool On);
  virtual bool setProblem () = 0;
    // Input: Table() w/o OBJ::GetMult ()
    // Output: Problem
    // Return: true if success
    // Requires: !Problem
  void closeProblem ();
    // Output: !Problem
  void setVarSelected ();
    // Update: *Problem
  virtual bool setWeight () = 0;
   	// Return: true if success
    // Input: OBJ::W ()
    // Update: *Problem
   	//          multSum = GetObjSumW () !
  bool solve ();
    // Update: *Problem
    // Output: absCriterion, beta - if failed then NAN
    // Return: true if success
    // Requires: after SetVarSelected (), SetWeight () = true
    // Invokes: Problem2Result ()
protected:
  virtual void problem2Result ();
    // Output: absCriterion, beta
public:
#endif

  // Requires: after Solve[Unconstrained] ()
  Real predict (size_t objNum) const;
    // Time: O(pn)
  Real getResidual (size_t objNum) const
    { return target [objNum] - predict (objNum); }
  virtual void setAbsCriterion () = 0;
    // Input: beta
    // Output: absCriterion
    // Invokes: getResidual()
  virtual void setAbsCriterion (const RealAttr1& residual) = 0;
    // Input: beta, residual
    // Output: absCriterion
  virtual Real absCriterion2Error () const = 0;
    // Return: >= 0
    // Input: absCriterion
  virtual Real getConstTarget (Real &scatter) const = 0;
    // Return: Prediction of target as a constant
    // Output: scatter
  virtual Real getRelTargetCriterion (Real constTarget) const = 0;
    // Input: absCriterion
    // Return: Criterion relative to the prediction of target by constTarget
    //         = absCriterion / maxAbsCriterion
    //         If maxAbsCriterion = 0 then INF
  Real getRelCriterion () const
    { Real scatter;
    	return getRelTargetCriterion (getConstTarget (scatter)); 
    }
   	// "R2"
  virtual LinearNumPrediction* makeLinearNumPrediction (const Sample &sample_arg,
                                                        const Space1<NumAttr1> &space_arg,
                                                        const RealAttr1 &target_arg) const = 0;
  bool solveUnconstrainedAlternate (const RealAttr1* predictionAttr,
                                    bool betaNonNegative,
                                    uint maxIter,
                                    Real errorRelDiff);
    // Approximation by alternate optimization
    // Return: converged
    // Input: predictionAttr: may be nullptr
    // Update: beta
    // Output: absCriterion: !isNan()
    // Depends on space ordering
    // Requires: &predictionAttr->ds = &space.ds
    // Invokes: solveUnconstrained()
    // Time: O(maxIter * p * n)
  bool solveUnconstrainedFast (const RealAttr1* predictionAttr,
                               bool betaNonNegative,
                               uint maxIter,
                               Real errorRelDiff);
    // Return: converged
    // Input: predictionAttr: may be nullptr
    // Update: beta
    // Output: absCriterion: !isNan()
    // Invokes: solveUnconstrained() and/or solveUnconstrainedAlternate()
    // Time: O(maxIter * p * n)

  RealAttr1& makePredictionAttr (const string &nameSuffix,
                                 Dataset &ds) const;
  RealAttr1& makeResidualAttr (const string &attrName,
                               const RealAttr1* predictionAttr,
                               Dataset &ds) const;
    // Input: predictionAttr: may be nullptr
    // Requires: & predictionAttr->ds = & space. ds

#if 0
  // Graphics
  TIME_SERIES_WINDOW* drawTimeSeriesPrediction (const char* WindowName) const;
    // Invokes: MakePredictionAttr ()
    //          new GRAPH_WINDOW ()
#endif

protected:
//Constraint* getBetaConstraint (size_t BCNum) const
  //{ return (Constraint*) constraints. AtCheck (BCNum); }
#if LIN_PROG
  void setBetaConstraints (SPARSE_MATRIX &lhs,
	                  		   size_t         InitRow) const;
    // Input: constraints
    // Output: lhs, Problem
    // Requries: beta [] are in the first columns of lhs
    // InitRow is the first row of lhs containing a beta constraint
#endif
};



#if LIN_PROG
  LinearNumPrediction* GetLR (METRIC Metric,
              						    size_t MaxRow,
              						    size_t predictorsSize,
              						    bool WithIntercept);
    // Return: !=nullptr
#endif




#if LIN_PROG
class L1_LINEAR_NUM_PREDICT : LinearNumPrediction
// absCriterion = weighted sum of absolute values of residuals
// betaSD -??
{
typedef LinearNumPrediction inherited;

public:
  L1_LINEAR_NUM_PREDICT (DataSet* initTable,
                         size_t initSpaceNum,
                         const RealAttr1* initTarget,
           							 bool WithIntercept);


  bool SetProblem ();
  bool SetWeight ();
  void SetAbsCriterion ();
  Real AbsCriterion2Error () const;
  Real GetConstTarget () const;
  Real GetRelCriterion (Real ConstTarget) const;
};
#endif




struct L2LinearNumPrediction : LinearNumPrediction
// absCriterion = weighted sum of squares of residuals ("MSD")
// targetConst = average of target
{
  Matrix attrSim;
    // Attribute similarity matrix: X^t * X
    // size = (space.size(), space.size())
  Matrix betaCovariance;
    // size = (space.size(), space.size())
#if LIN_PROG
  Real absCriterionUnconstrained {NAN};
    // Init: NAN
  MVector betaUnconstrained;
    // size = space.size()
#endif


  L2LinearNumPrediction (const Sample &sample_arg,
                         const Space1<NumAttr1> &space_arg,
                         const RealAttr1 &target_arg)
    : LinearNumPrediction (sample_arg, space_arg, target_arg)
  #if LIN_PROG
    , betaUnconstrained (space_arg. size ())
  #endif
    { resize (); }


  void resize () 
    { LinearNumPrediction::resize ();
      attrSim.        resize (false, space. size (), space. size ());
      betaCovariance. resize (false, space. size (), space. size ());
    }
  void solveUnconstrained () final;
    // Output: attrSim, betaCovariance, betaSD
    // n < p => isNan(absCriterion)
    // Time: O(p^2 (n + p))
#if LIN_PROG
  bool setProblem ();
  bool setWeight ();
    // Input: attrSim
  void problem2Result ();
#endif
  void setAbsCriterion () final;
  void setAbsCriterion (const RealAttr1& residual) final;
  Real absCriterion2Error () const final
    { return sqrt (absCriterion / sample. mult_sum); }
    // MSE
  Real getConstTarget (Real &scatter) const final;
  Real getRelTargetCriterion (Real constTarget) const final;
  L2LinearNumPrediction* makeLinearNumPrediction (const Sample &sample_arg,
                                                  const Space1<NumAttr1> &space_arg,
                                                  const RealAttr1 &target_arg) const final
    { return new L2LinearNumPrediction (sample_arg, space_arg, target_arg); }
};




struct EigensLinearRegression : Root
{
  Dataset ds;
  const RealAttr1& target;
  L2LinearNumPrediction lr;
    
  
  explicit EigensLinearRegression (const Eigens &eigens)
    : ds (eigens)
    , target (* checkPtr (ds. name2attr ("EigenValueFrac")) -> asRealAttr1 ())
    , lr (Sample (ds), Space1<NumAttr1> (ds, true). removeAttr (target), target)
    { lr. solveUnconstrained (); }
  void qc () const override;
  void saveText (ostream &os) const override
    { lr. saveText (os); }
  JsonMap* toJson (JsonContainer* parent,
                   const string& name) const override;
};




#if LIN_PROG
class L_INF_LINEAR_NUM_PREDICT : LinearNumPrediction
// betaSD -??
{
typedef LinearNumPrediction inherited;

public:
  L_INF_LINEAR_NUM_PREDICT (DataSet* initTable,
                            size_t            initSpaceNum,
                            const RealAttr1* initTarget,
              							bool            WithIntercept);


  bool SetProblem ();
  bool SetWeight ();  
  void SetAbsCriterion ();
  Real AbsCriterion2Error () const;
  Real GetConstTarget () const;
  Real GetRelCriterion (Real ConstTarget) const;
};
#endif





/////////////////////// ExprNumPredict /////////////////////

#if 0
struct ExprNumPredict : Prediction
// Non-linear regression
{
  // INPUT
  const RealAttr1 &target;
    // In space
  size_t maxParamNum;
  Matrix BetaDelta;
    // Precision of beta
    // Size = (MaxParamNum, 1)
    // Init: NAN

  // OUTPUT
  MVector beta;
    // Coefficients of prediction function
    // Size = (MaxParamNum, 1)
    // Init: NAN
  Real absCriterion;
    // Absolute criterion
    // >= 0, -> min
    // Init: NAN
    // Metric L_2


  ExprNumPredict (Space1<NumAttr1> &space,
                  const RealAttr1 &target_arg
                  size_t initMaxParamNum);
  void qc () const;
    // Requires: No missing values
  
  
  virtual Real GetTargetPredicted (size_t         objNum,
                                    const Matrix &TrialBeta) const = 0;
  virtual void GetTargetPredictedDeriv (size_t         objNum,
                                        const Matrix &TrialBeta,
                                        Matrix       &TargetPredictedDeriv) const;
    // Output: TargetPredictedDeriv
    // Requires: TargetPredictedDeriv: Size = (MaxParamNum, 1)
    // Default: RunError
  bool Solve (size_t MaxIter);
    // Update: beta
    // Output: absCriterion
    // Return: true if converged
    // Requires: beta, BetaDelta are not NAN
    // Invokes: FUNCTION_MULT::OptimizeMarquardt (),
    //          GetTargetPredicted (), GetTargetPredictedDeriv ()
};
#endif



}



#endif

