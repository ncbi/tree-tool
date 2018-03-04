// prediction.cpp

#undef NDEBUG
#include "../common.inc"

#include "prediction.hpp"

#include "optim.hpp"




namespace DM_sp
{
using namespace std;
using namespace Common_sp;




// LogisticRegression

void LogisticRegression::qc () const
{
  if (! qc_on)
    return;
  Prediction::qc ();
    
  ASSERT (& target);

#ifndef NDEBUG
  FOR (size_t, i, space. size ())
  {
    const NumAttr1* a = space [i];
    ASSERT (a);
    if (! i)
    {
      ASSERT (a->isConstant ());
      ASSERT (a->getReal (0) == 1);
    }
  } 
  
  ASSERT (beta. size () == space. size ());   
  ASSERT (attrImportance. size () == space. size ());   
  
  if (! isNan (negLogLikelihood_ave))
  {
    ASSERT (! isNan (target_score_min));
    ASSERT (! isNan (non_target_score_max));
  }
#endif
}



Real LogisticRegression::getScore (size_t objNum) const
{
  Real s = 0;
  FOR (size_t, attrNum, space. size ())
    s += space [attrNum] -> getReal (objNum) * beta [attrNum];
  return s;
}



Real LogisticRegression::getNegLogLikelihood_ave () const
{
  Real s = 0;
  for (Iterator it (sample); it ();)  
  {
    const Real z = getScore (*it);
    s += it. mult * (log (exp (z) + 1) - z * target. getBool (*it));
  }
  ASSERT (s >= 0);
  return s / sample. mult_sum;
}



namespace 
{

struct LogRegFuncMult : FuncMult
{
  LogisticRegression& lr;
  

  LogRegFuncMult (LogisticRegression &lr_arg)
    : FuncMult (lr_arg. space. size ())
    , lr (lr_arg)
    {}



  Real f (const MVector &x)
    { 
      lr. beta = x; 
      return lr. getNegLogLikelihood_ave ();
    }

  void getGradient (const MVector &x,
         				    MVector &gradient)
    {
      lr. beta = x; 
      gradient. putAll (0);
      for (Iterator it (lr. sample); it ();)  
      {
        const Real a = it. mult * (lr. predict (*it) - lr. target. getBool (*it));
        FOR (size_t, argNum, maxArgNum)
      		gradient. putInc (false, argNum, 0, a * lr. space [argNum] -> getReal (*it));
      }
      gradient. putProdAll (1 / lr. sample. mult_sum);
    }
    

  void getHessian (const MVector &x,
        			     Matrix &hessian)
    { 
      lr. beta = x; 
      hessian. putAll (0.0);
      for (Iterator it (lr. sample); it ();)  
      {
        const Prob p = lr. predict (*it) ;
        const Real a = it. mult * p * (1 - p);
        FOR (size_t, argNum1, maxArgNum)
        {
          const Real x1 = lr. space [argNum1] -> getReal (*it);
          FOR (size_t, argNum2, maxArgNum)
        		hessian. putInc (false, argNum1, argNum2, a * x1 * lr. space [argNum2] -> getReal (*it));
        }
      }
      hessian. putProdAll (1 / lr. sample. mult_sum);
    }
};

}



void LogisticRegression::simulateTarget (ulong seed)
{ 
  Rand rand (seed);
  FOR (size_t, i, space. ds. objs. size ())
    const_cast <BoolAttr1&> (target). setBool (i, (ebool) (rand. getProb () <= predict (i)));
}



void LogisticRegression::solve ()
{
  // Check for trivially separating attributes ??
  
  ASSERT (! target. existsMissing ());
  ASSERT (beta. size () >= 1);

  target_score_min = INF;
  non_target_score_max = -INF;

  if (beta. size () == 1)
  {
    beta [0] = logit (target. getProb (sample));
    negLogLikelihood_ave = getNegLogLikelihood_ave ();
    return;
  }

  negLogLikelihood_ave = NAN;
  beta. putAll (0);  // PAR

  LogRegFuncMult f (*this);

  MVector x (beta);

  MVector dX (space. size ());
  FOR (size_t, attrNum, space. size ())
  {
    Real average, scatter;
    space [attrNum] -> getAverageScatter (sample, average, scatter);
    ASSERT (! isNan (scatter));
    dX. put (false, attrNum, 0,
      max (1e-6 * sqrt (scatter), 1e-6 /*Intercept precision*/));  // PAR
  }

  if (f. optimizeMarquardt (true, 1e-6, x, dX, 1000))  // PAR ??
  {
    negLogLikelihood_ave = getNegLogLikelihood_ave ();
    ASSERT (negLogLikelihood_ave >= 0);
  }
  beta = x;    
  
  for (Iterator it (sample); it ();)  
  {
    const Real score = getScore (*it);
    if (target. getBool (*it))
      minimize (target_score_min, score);
    else
      maximize (non_target_score_max, score);
  }
}



void LogisticRegression::test (ulong seed)
{ 
  simulateTarget (seed);
//setMult ();
  const Real negLogLikelihood_ave_old = getNegLogLikelihood_ave ();
  if (verbose ())
  {
    cout << "-logL original = " << negLogLikelihood_ave_old << endl;
    cout << endl;
  }
  
  solve ();
  
  if (verbose ())
    cout << "-logL optimized = " << negLogLikelihood_ave << endl;

  ASSERT (geReal (negLogLikelihood_ave_old, negLogLikelihood_ave));
}



void LogisticRegression::setAttrImportance ()
{
  ASSERT (! isNan (negLogLikelihood_ave));
  
//Temorarily set Obj::mult to 0 if !Obj::active() ??

#ifndef NDEBUG
  const VectorPtr<NumAttr1> tmp (space); 
#endif

  attrImportance [0] = NAN;
  const Real negLogLikelihood_ave_old = negLogLikelihood_ave;
  for (Iter <VectorPtr<NumAttr1>> iter (space); iter. next (); )
  {
    const size_t i = iter. getIndex ();
    if (i == 0)
      continue;

    const NumAttr1* attr = iter. erase ();  
    resize ();  
    
    solve ();  
    const Real negLogLikelihood_ave_without_attr = negLogLikelihood_ave;
    const bool separated = getSeparated ();
    
    iter. insert (attr);
    resize ();  

    attrImportance [i] = (negLogLikelihood_ave_without_attr - negLogLikelihood_ave_old) / negLogLikelihood_ave_old;
    IMPLY (! separated, ! negative (attrImportance [i]));
  }
  
#ifndef NDEBUG
  qc ();  
  ASSERT (tmp == space);  
#endif
}




///////////////////////////////////////////////////////////////////

Constraint::Constraint (const LinearNumPrediction &pred_arg)
: pred (pred_arg)
, lhs (pred_arg. space. size ())
, rhs (0)
, sense ('E')
{ 
  lhs. putAll (0); 
}




/////////////////////////// LinearNumPrediction /////////////////////////////

void LinearNumPrediction::qc () const
{
  if (! qc_on)
    return;
  Prediction::qc ();

  ASSERT (& target);
  
  for (const auto attr : space)
    ASSERT (attr);

  for (const Constraint* cons : constraints)
    ASSERT (& cons->pred == this);
}



void LinearNumPrediction::saveText (ostream &os) const
{ 
  os << "absCriterion = " << absCriterion << endl;
  FOR (size_t, i, beta. size ())
    os << "beta[" << i << "] = " << beta [i] << " +- " << betaSD [i] << endl;
}



bool LinearNumPrediction::solveUnconstrainedAlternate (const RealAttr1* predictionAttr,
                                                       bool betaNonNegative,
                                                       uint maxIter,
                                                       Real errorRelDiff)
{
  ASSERT (maxIter);
  ASSERT (errorRelDiff > 0);
  ASSERT (beta. defined ());
#ifndef NDEBUG
  if (betaNonNegative)
    FOR (size_t, attrNum, space. size ())
      ASSERT (beta [attrNum] >= 0);
#endif

  absCriterion = NAN;
  
  Dataset& ds = const_cast <Dataset&> (space. ds);
  
  Common_sp::AutoPtr<RealAttr1> residual (& makeResidualAttr ("_res", predictionAttr, ds));  

  setAbsCriterion (*residual);
  const Real errorPrev_init = absCriterion2Error ();
  ASSERT (! isNan (errorPrev_init));  
  
  Common_sp::AutoPtr<RealAttr1> target1 (new RealAttr1 ("_target", ds)); 

    
  bool ok = true;
  Progress prog;
  Space1<NumAttr1> sp (ds, false);
  Real errorPrev = errorPrev_init;
  FOR (uint, iter, maxIter)
  {
    prog (real2str (absCriterion /*errorPrev*/, 6));  // PAR
    FOR (size_t, attrNum, space. size ())
    {
      // beta[attrNum]
      // Assume beta[j] are correct for all j != attrNum
      const NumAttr1& predictor = * space [attrNum];
      sp. clear ();
      sp << & predictor;
      for (Iterator it (sample); it ();)  
        (*target1) [*it] = (*residual) [*it] + beta [attrNum] * predictor. getReal (*it);
      Unverbose unv;
      Common_sp::AutoPtr<LinearNumPrediction> lp (makeLinearNumPrediction (sample, sp, *target1));
      ASSERT (eqReal (lp->sample. mult_sum, sample. mult_sum));
      lp->solveUnconstrained ();
      lp->qc ();
      if (isNan (lp->absCriterion))
        continue;
      if (betaNonNegative)
        maximize (lp->beta [0], 0.0);
      const Real betaDelta = beta [attrNum] - lp->beta [0];
      for (Iterator it (sample); it ();)  
        (*residual) [*it] = (*residual) [*it] + betaDelta * predictor. getReal (*it);
      beta [attrNum] = lp->beta [0];
    //cout << "attrNum = " << attrNum << "  beta = " << beta [attrNum] << endl;
      // Save lp->absCriterion as a property of beta[attrNum] ??
    }
    setAbsCriterion (*residual);
    const Real error = absCriterion2Error ();
    {
      Unverbose unv;
      if (verbose ())
        cout << error << " " << errorPrev << endl;
    }
    if (errorPrev / error - 1 <= errorRelDiff) 
    {
      if (error >= errorPrev_init)
        ok = false;
      goto quit;
    }
    errorPrev = error;
  }
  ok = false;

    
quit:
  setAbsCriterion ();
  ASSERT (! isNan (absCriterion));
    
  return ok;
}




bool LinearNumPrediction::solveUnconstrainedFast (const RealAttr1* predictionAttr,
                                                  bool betaNonNegative,
                                                  uint maxIter,
                                                  Real errorRelDiff)
{
  bool solved = false;

  const MVector beta_init (beta);
  ASSERT (beta_init. defined ());
  IMPLY (betaNonNegative, ! negative (beta. min ()));
  
  absCriterion = NAN;
  
  if (beta. size () /*p*/ <= 2 * maxIter)  
  {
    if (verbose ())
      cout << "Direct solution ..." << endl;
    solveUnconstrained ();   // Time = O(p^2 n / 2)
    if (! isNan (absCriterion))
    {
      solved = true;
      if (betaNonNegative)
        FOR (size_t, i, beta. size ())
          if (maximize (beta [i], 0.0))
            solved = false;
    }
  }

  if (! solved)
  {
    if (verbose ())
      cout << "Alternate optimization ..." << endl;
    beta = beta_init;
  #if 0
    // Use corrected beta after solveUnconstrained() ??
    setAbsCriterion ();
    const Real absCriterion_old = absCriterion;
  #endif
    solved = solveUnconstrainedAlternate (predictionAttr, betaNonNegative, maxIter, errorRelDiff);  // Time = O(maxIter * p * n)
  #if 0
    if (absCriterion > absCriterion_old)
    {
      beta = beta_init;
      solved = solveUnconstrainedAlternate (predictionAttr, betaNonNegative, maxIter, errorRelDiff);  // Time = O(maxIter * p * n)
    }
  #endif
  }
  
  return solved;
}



RealAttr1& LinearNumPrediction::makePredictionAttr (const string &nameSuffix,
                                                    Dataset &ds) const
{ 
  ASSERT (& ds == sample. ds);

  auto& prediction = * new RealAttr1 (target. name + "_" + nameSuffix, ds, target. decimals); 
  for (Iterator it (sample); it ();)  
    prediction [*it] = predict (*it);

  return prediction;  
}



RealAttr1& LinearNumPrediction::makeResidualAttr (const string &attrName,
                                                  const RealAttr1* predictionAttr,
                                                  Dataset &ds) const
{ 
  IMPLY (predictionAttr, & predictionAttr->ds == & space. ds);
  ASSERT (& ds == sample. ds);

  auto& residual = * new RealAttr1 (attrName, ds, target. decimals); 
  for (Iterator it (sample); it ();)  
    residual [*it] = predictionAttr 
                       ? target [*it] - (*predictionAttr) [*it] 
                       : getResidual (*it);

  return residual;
}



#if 0
void LinearNumPrediction::removeTargetOutlierValues ()
{
	??
  Normal /*Exponential*/ distr;
  outlier_min = exp (residualAttr->distr2outlier (sample, distr, outlier_EValue_max));  
}
#endif



#if LIN_PROG
// Problem

void LinearNumPrediction::setProblemVerbose (bool On)
{
  ASSERT (Problem);
  
  
  if (! Problem->SetVerbose (On))
    RunError;
}



void LinearNumPrediction::closeProblem ()
{
  delete Problem;
  Problem = nullptr;
}



void LinearNumPrediction::setVarSelected ()
{
  ASSERT (Problem);


  For (col, space. size ())
  {
 		Problem->ChangeBound (col, 'L', - CPX_INFBOUND);
 		Problem->ChangeBound (col, 'U',   CPX_INFBOUND);
  }
}



bool LinearNumPrediction::solve ()
{
  ASSERT (Problem);


  beta. clear ();
  absCriterion = NAN;

  if (! Problem->Solve ())
    return false;

  ASSERT (Problem->GetCurMaxRow () == Table () -> GetActiveObjNum ());
  ASSERT (Problem->GetCurMaxCol () == Problem->space. size ());

  Problem2Result ();
  

  return true;
}



void LinearNumPrediction::problem2Result ()
{
  ASSERT (Problem);


  absCriterion = Problem->ObjVal;

  For (Col, space. size ())
    beta. put (false, Col, 0, Problem->X [Col]);
}
#endif



#if 0
Real LinearNumPrediction::getIntercept () const
{
  ASSERT (! isNan (absCriterion));
  

  if (interceptAttr)
    return beta. get (false, space. size () - 1, 0);
  return 0.0;
}
#endif



Real LinearNumPrediction::predict (size_t objNum) const
{
  Real s = 0;
  FOR (size_t, col, space. size ())
    s += beta [col] * space [col] -> getReal (objNum);
  return s;
}




#if 0
TIME_SERIES_WINDOW* LinearNumPrediction::DrawTimeSeriesPrediction (const char* WindowName) const
{
  ASSERT (! StrIsEmpty (WindowName));
  

  // PredAttr
  RealAttr1* PredAttr = MakePredictionAttr ("_pred");
  ASSERT (PredAttr);
  
 
  // TSW
  TIME_SERIES_WINDOW* TSW = new TIME_SERIES_WINDOW (WindowName, Table ());
  ASSERT (TSW);
  
  TSW->AddAttr (target);
  TSW->AddAttr (PredAttr);
  
  TSW->ForceSameYAxis ();
  
  
  delete PredAttr;
  
  
  return TSW;
}
#endif



#if LIN_PROG
void LinearNumPrediction::SetBetaConstraints (SPARSE_MATRIX &lhs,
                               										    size_t          InitRow) const
{
  ASSERT (Problem);
  

  ForCollection (BCNum, constraints)
    {
      Constraint* BC = GetBetaConstraint (BCNum);
      
      const size_t Row = InitRow + BCNum;
      
      For (Col, space. size ())
      		lhs. Append (Row, Col, BC->lhs [Col]);
      
      Problem->rhs         [Row] = BC->rhs;
      Problem->ConstrSense [Row] = BC->sense;
    }
}
#endif



#if LIN_PROG
LinearNumPrediction* GetLR (METRIC  Metric,
                   						  size_t    MaxRow,
                   						  size_t    predictorsSize,
                   						  bool WithIntercept)
{
  LinearNumPrediction* LR = nullptr;
  switch (Metric)
    {
      case 0: LR = new L1_LINEAR_NUM_PREDICT    (MaxRow, space. size (), WithIntercept);
      Case 1: LR = new L2LinearNumPrediction    (MaxRow, space. size (), WithIntercept);
      Case 2: LR = new L_INF_LINEAR_NUM_PREDICT (MaxRow, space. size (), WithIntercept);
      Default: RunError;
    }
  ASSERT (LR);


  return LR;
}
#endif





#if LIN_PROG
////////////////////////// L1_LINEAR_NUM_PREDICT ////////////////////////////

L1_LINEAR_NUM_PREDICT::L1_LINEAR_NUM_PREDICT (DataSet*  initTable,
                                              size_t             initSpaceNum,
                                              const RealAttr1* initTarget,
                                						  bool             WithIntercept):
  inherited ("LINEAR PREDICTION IN L_1",
             initTable, initSpaceNum, initTarget, WithIntercept),
{
}

					    

bool L1_LINEAR_NUM_PREDICT::SetProblem ()
/* LP Problem:

   sum_i (Weight_i * C_i + Weight_i * D_i) -> min
     s.t.
   beta^t * X_i + C_i - D_i = Y_i
   C_i >= 0
   D_i >= 0
   Constraints on beta
   beta are unbounded
   
   In an optimal solution C_i = 0 or D_i = 0

*/
{
  ASSERT (! Problem);


  // Problem setup
  // Columns for LP: beta, C, D
  const int LPMaxCol = space. size () + 2 * MaxRow;

  // Rows for LP: MaxRow, constraints
  const int LPMaxRow = MaxRow + constraints. size ();


  Problem = new CPLEX ("LR_L1", LPMaxRow, LPMaxCol, LPMaxRow * LPMaxCol /*??*/);
  ASSERT (Problem);


  Problem->SetVerbose (false);


  SPARSE_MATRIX lhs;


  // Constraints
  For (Row, MaxRow)
  {
		For (Col, space. size ())
		  lhs. Append (Row, Col, X [Row] [Col]);
		lhs. Append (Row, s[ace size () + 2 * Row + 0,  1.0);
		lhs. Append (Row, space size () + 2 * Row + 1, -1.0);

		Problem->ConstrSense [Row] = 'E';
		Problem->rhs         [Row] = Y [Row];
  }
  SetBetaConstraints (lhs, MaxRow);


  bool Ok = Problem->SetLP (CPX_MIN, LPMaxRow, LPMaxCol, lhs, 0, 0);


  return Ok;
}



bool L1_LINEAR_NUM_PREDICT::SetWeight ()
{
  ASSERT (Problem);


  For (Row, MaxRow)
    For (i, 2)
      if (! Problem->ChangeObjCoeff (space. size () + 2 * Row + i, Weight [Row]))
      		RunError;

  multSum = GetObjSumW ();
  

  return true;
}



void L1_LINEAR_NUM_PREDICT::SetAbsCriterion ()
{
  absCriterion = 0.0;
  For (Row, MaxRow)
    absCriterion += AbsFloat (Residual [Row]) * Weight [Row];
}



Real L1_LINEAR_NUM_PREDICT::AbsCriterion2Error () const
{
  return absCriterion / multSum;
}



class Y_WEIGHT: public ROOT
{
typedef ROOT inherited;

public:
  Real Y;
  Real Weight;

  Y_WEIGHT (Real initY,
		          Real initWeight);
};



Y_WEIGHT::Y_WEIGHT (Real initY,
            				    Real initWeight):
  inherited (),

  Y         (initY),
  Weight    (initWeight)
{
  ASSERT (Weight >= 0.0);
}



static bool YWeightGreaterEqual (const ROOT* A,
                  						            const ROOT* B)
{
  return ((Y_WEIGHT*) A) -> Y >= ((Y_WEIGHT*) B) -> Y;
}



Real L1_LINEAR_NUM_PREDICT::GetConstTarget () const
{
  Real WSum = 0.0;
  {
    For (Row, MaxRow)
      WSum += Weight [Row];
  }
  const Real WThreshold = WSum / 2;

      
  HEAP YWeightHeap (YWeightGreaterEqual);
  {
    For (Row, MaxRow)
      YWeightHeap. Insert (new Y_WEIGHT (Y [Row], Weight [Row]));
  }
  YWeightHeap. Sort ();


  ConstTarget = NAN;
  Real WSortSum = 0.0;
  {
    For (Row, MaxRow)
    {
   	  const Y_WEIGHT* YWeight = (Y_WEIGHT*) YWeightHeap. Array. AtCheck (Row);
   	  WSortSum += YWeight->Weight;
   	  if (GreaterEqualFloat (WSortSum, WThreshold))
     	  {
     	    ConstTarget = YWeight->Y;
     	    break;
     	  }
    }
  }
  ASSERT (! isNan (ConstTarget));
  
  
  return ConstTarget;
}



Real L1_LINEAR_NUM_PREDICT::GetRelCriterion (Real ConstTarget) const
{
  if (isNan (absCriterion) ||
      isNan (ConstTarget))
    return NAN;


  // MaxAbsCriterion
  Real MaxAbsCriterion = 0.0;
  For (Row, MaxRow)
    MaxAbsCriterion += Weight [Row] * AbsFloat (Y [Row] - ConstTarget);


  if (nullFloat (MaxAbsCriterion))
    return  -INF;
  else
    return  1.0 - absCriterion / MaxAbsCriterion;
}
#endif




////////////////////////// L2LinearNumPrediction ////////////////////////////

void L2LinearNumPrediction::solveUnconstrained ()
{
  absCriterion = NAN;
  const bool existsMissing = getExistsMissing ();


  setAttrSim (attrSim, existsMissing);

  // betaCovariance: temporary
  betaCovariance = attrSim;
  if (   betaCovariance. inverse (). isZero ()
    //|| ! attrSim. checkInverse (betaCovariance, false)
     )
    return;
        
  // xy
  MVector xy (space. size ());
  FFOR (size_t, attrNum, space. size ())
  {
    const NumAttr1& a = * space [attrNum];
    Real s = 0;
    Real mult_sum = 0;
    FFOR (size_t, i, sample. mult. size ())
      if (const Real mult = sample. mult [i])
      {
      	const Real x = a. getReal (i);
      	if (existsMissing)
        {
	      	if (isNan (x))
	      		continue;
	      	if (isNan (target [i]))
	      		continue;
	        mult_sum += mult;
	      }
        s += x * target [i] * mult;
      }
    if (existsMissing)
	    xy [attrNum] = mult_sum > 0 ? s / mult_sum : 0;
	  else
	  {
	    ASSERT (! isNan (s));
	    xy [attrNum] = s;
	  }
 	}

  // beta
  beta. multiply (false, 
                  betaCovariance, false, 
                  xy,             false);


  // absCriterion
  const Real yHat2 = beta. multiplyVec (xy);
  ASSERT (yHat2 >= 0);

  Real y2 = 0;
  {
	  Real mult_sum = 0;
	  FFOR (size_t, i, sample. mult. size ())
	    if (const Real mult = sample. mult [i])
	    	if (! isNan (target [i]))
	    	{
	    		mult_sum += mult;
			    y2 += sqr (target [i]) * mult;
			  }
		if (existsMissing)
			y2 /= mult_sum;
	}
	ASSERT (y2 >= 0);

  absCriterion = y2 - yHat2;
  if (existsMissing)
  {
  	absCriterion *= sample. mult_sum;
  	maximize (absCriterion, 0.0);
  }
  else
  {
	  if (! positive (absCriterion))
	    setAbsCriterion ();
	}
  ASSERT (absCriterion >= 0);


  // betaCovariance  
  const Real epsilonVar = absCriterion / sample. mult_sum;  // biased
  betaCovariance. putProdAll (existsMissing ? absCriterion : epsilonVar);  
  
  // betaSD
  betaSD. diag2row (true, 0, betaCovariance);
  betaSD. sqrtAll ();
}



#if LIN_PROG
bool L2LinearNumPrediction::SetProblem ()
/* QP Problem:

   beta^t X^t diag (Weight) X beta -> min
     s.t.
   Constraints on beta (derived from constraints on beta~)
   beta is unbounded,
     where
   beta = beta~ - beta^,
   beta^ = unconstrained weighted L2 LR coefficients of Y on X
   beta~ is the result

*/
{
  ASSERT (! Problem);



  // Problem setup
  // Columns for LP: beta
  const int LPMaxCol = space. size ();

  // Rows for LP: constraints
  const int LPMaxRow = constraints. size ();


  Problem = new CPLEX ("LR_L2", LPMaxRow, LPMaxCol, LPMaxRow * LPMaxCol);
  ASSERT (Problem);

  Problem->SetVerbose (false);


  SPARSE_MATRIX lhs;


  // Constraints
  SetBetaConstraints (lhs, 0);


  Real** QM = (Real**) NewMatrix (LPMaxCol, LPMaxCol, sizeof (Real));


  bool Ok = Problem->SetLP (CPX_MIN, LPMaxRow, LPMaxCol, lhs, QM, 0);


  FreeMatrix ((void**) QM, LPMaxCol);


  return Ok;
}



bool L2LinearNumPrediction::SetWeight ()
{
  ASSERT (Problem);


  // betaUnconstrained (= beta^)
  SolveUnconstrained ();
  betaUnconstrained. Copy (false, Col, false);
  absCriterionUnconstrained = absCriterion;
  if (isNan (absCriterionUnconstrained))
   	return false;
  ASSERT (absCriterionUnconstrained >= 0.0);


  // Constraints on beta~ -> constraints on beta
  ForCollection (BCNum, constraints)
  {
    const Constraint* BC = GetBetaConstraint (BCNum);
    
    Real C = 0.0;
    For (col, space. size ())
    		if (VarSelected [col])
    		  C += BC->lhs [col] * betaUnconstrained [col];

    if (! Problem->ChangeRHSCoeff (BCNum, BC->rhs - C))
    		RunError;
  }
  

  // QM
  For (i, space. size ())
    FOR_START (size_t, j, i, space. size ())
      if (! Problem->ChangeQMCoeff (i, j, attrSim [i] [j]))
      		RunError;


  multSum = GetObjSumW ();


  return true;
}



void L2LinearNumPrediction::Problem2Result ()
{
  ASSERT (Problem);


  absCriterion = Problem->ObjVal + absCriterionUnconstrained;

  beta. Copy (false, betaUnconstrained, false);
  For (attrNum, space. size ())
    beta. putInc (false, attrNum, 0, Problem->X [attrNum]);
}
#endif



void L2LinearNumPrediction::setAbsCriterion ()
{
  absCriterion = 0;
  Real mult_sum = 0;
  FFOR (size_t, i, sample. mult. size ())
    if (const Real mult = sample. mult [i])
    	{
    		const Real x = getResidual (i);
    		if (! isNan (x))
    		{
	        absCriterion += sqr (x) * mult; 
	        mult_sum += mult;
	      }
	    }
  absCriterion *= sample. mult_sum / mult_sum;
}


void L2LinearNumPrediction::setAbsCriterion (const RealAttr1& residual)
{
  ASSERT (& residual. ds == & space. ds);
  
  absCriterion = 0;
  Real mult_sum = 0;
  FFOR (size_t, i, sample. mult. size ())
    if (const Real mult = sample. mult [i])
    	{
    		const Real x = residual [i];
    		if (! isNan (x))
    		{
	        absCriterion += sqr (x) * mult; 
	        mult_sum += mult;
	        ASSERT (! isNan (absCriterion));
	        ASSERT (! isNan (mult_sum));	        
	      }
	    }
  absCriterion *= sample. mult_sum / mult_sum;
}
  


Real L2LinearNumPrediction::getConstTarget (Real &scatter) const
{ 
  Real average; 
  target. getAverageScatter (sample, average, scatter);
  ASSERT (! isNan (scatter));
  return average;
}



Real L2LinearNumPrediction::getRelTargetCriterion (Real constTarget) const
{
  if (   isNan (absCriterion) 
      || isNan (constTarget)
     )
    return NAN;

  // maxAbsCriterion
  Real maxAbsCriterion = 0;
  for (Iterator it (sample); it ();)  
    maxAbsCriterion += it. mult * sqr (target [*it] - constTarget);

  if (nullReal (maxAbsCriterion))
    return INF;
  else
    return absCriterion / maxAbsCriterion;
}




///////////////////////// EigensLinearRegression ///////////////////////////


void EigensLinearRegression::qc () const
{ 
  if (! qc_on)
    return;
  ds. qc ();
  lr. qc ();
  ASSERT (lr. beta. size () == 2);  
}


JsonMap* EigensLinearRegression::toJson (JsonContainer* parent,
                                         const string& name) const
{
  if (ds. objs. size () <= 2)
    return nullptr;

  const size_t decimals = 4;  // PAR
  JsonMap* jLr = new JsonMap (parent, name);
  new JsonDouble (1 - lr. getRelCriterion (), decimals, jLr, "log_scale_R2");  
  new JsonDouble (exp (lr. beta [0]), decimals, jLr, "intercept");  
  new JsonDouble (lr. beta [1], decimals, jLr, "coeff");  
  
  return jLr;
}




#if LIN_PROG
///////////////////////// L_INF_LINEAR_NUM_PREDICT ///////////////////////////

L_INF_LINEAR_NUM_PREDICT::L_INF_LINEAR_NUM_PREDICT (DataSet*  initTable,
                                                    size_t             initSpaceNum,
                                                    const RealAttr1* initTarget,
                                      						  bool             WithIntercept):
  inherited ("LINEAR PREDICTION IN L_INF",
             initTable, initSpaceNum, initTarget, WithIntercept),
{
}

					    

bool L_INF_LINEAR_NUM_PREDICT::SetProblem ()
/* LP Problem:

   Alpha -> min
     s.t.
   beta^t * X_i + C_i - D_i = Y_i
   Weight_i * C_i - Alpha <= 0.0
   Weight_i * D_i - Alpha <= 0.0
   C_i >= 0.0
   D_i >= 0.0
   Alpha >= 0.0  (redundant)
   Constraints on beta
   beta are unbounded

*/
{
  ASSERT (! Problem);
  

  // Problem setup
  // Columns for LP: beta, C, D, Alpha
  const int LPMaxCol = space. size () + 2 * MaxRow + 1;

  // Rows for LP: X, C, D, constraints
  const int LPMaxRow = 3 * MaxRow + constraints. size ();


  Problem = new CPLEX ("LR_LInf", LPMaxRow, LPMaxCol, LPMaxRow * LPMaxCol /*??*/);
  ASSERT (Problem);

  Problem->SetVerbose (false);


  SPARSE_MATRIX lhs;


  Problem->ObjCoeff [LPMaxCol - 1] = 1.0;
  
  
  // Constraints
  {
    For (Row, MaxRow)
      {
      		For (Col, space. size ())
      	      lhs. Append (Row, Col, X [Row] [Col]);
      		lhs. Append (Row, space. size () +          Row,  1.0);  // C
      		lhs. Append (Row, space. size () + MaxRow + Row, -1.0);  // D
      		Problem->ConstrSense [Row] = 'E';
      		Problem->rhs         [Row] = Y [Row];
      
      		// C
      		lhs. Append (MaxRow + Row, LPMaxCol - 1, -1.0);
      		Problem->ConstrSense [MaxRow + Row] = 'L';
      
      		// D
      		lhs. Append (2 * MaxRow + Row, LPMaxCol - 1, -1.0);
      		Problem->ConstrSense [2 * MaxRow + Row] = 'L';
      }
  }
  SetBetaConstraints (lhs, 3 * MaxRow);


  bool Ok = Problem->SetLP (CPX_MIN, LPMaxRow, LPMaxCol, lhs, 0, 0);



  return Ok;
}



bool L_INF_LINEAR_NUM_PREDICT::SetWeight ()
{
  ASSERT (Problem);
  

  For (Row, MaxRow)
    For (i, 2)
      if (! Problem->ChangeLHSCoeff ((i + 1) * MaxRow + Row, 
                        								     space. size () + i * MaxRow + Row, 
                        								     Weight [Row]))
      		RunError;

  multSum = GetObjSumW ();


  return true;
}



void L_INF_LINEAR_NUM_PREDICT::SetAbsCriterion ()
{
  absCriterion = 0.0;
  for (Iterator it (sample); it ();)  
    maximize (absCriterion, AbsFloat (GetResidual (*it) * it->weight);
}



Real L_INF_LINEAR_NUM_PREDICT::AbsCriterion2Error () const
{
  return absCriterion;
}



Real L_INF_LINEAR_NUM_PREDICT::GetConstTarget () const
{ 
  ConstTarget = NAN;
  Real MaxAbsCriterion = -INF;
  For (Row1, MaxRow)
    For (Row2, MaxRow)
      if (Y [Row2] >= Y [Row1])
      		if (maximize (MaxAbsCriterion, 
		                    (Y [Row2] - Y [Row1]) * Weight [Row1] * Weight [Row2] 
              		      / (Weight [Row1] + Weight [Row2])))
		        ConstTarget = (Y [Row1] * Weight [Row1] + Y [Row2] * Weight [Row2]) /
	                       (Weight [Row1] + Weight [Row2]);
  ASSERT (! isNan (ConstTarget));
  
  
  return ConstTarget;
}



Real L_INF_LINEAR_NUM_PREDICT::GetRelCriterion (Real ConstTarget) const
{
  if (isNan (absCriterion) ||
      isNan (ConstTarget))
    return NAN;


  // MaxAbsCriterion
  Real MaxAbsCriterion = 0.0;
  for (Iterator it (sample); it ();)  
    maximize (MaxAbsCriterion, it->weight * AbsFloat (target [*it] - ConstTarget));
  

  if (nullFloat (MaxAbsCriterion))
    return  -INF;
  else
    return  1.0 - absCriterion / MaxAbsCriterion;
}
#endif





//////////////////////////// ExprNumPredict //////////////////////////

#if 0
ExprNumPredict::ExprNumPredict (const Space1<NumAttr1>* initAttrSpace,
                                const RealAttr1*  initTarget,
                                size_t              initMaxParamNum):
  space    (initAttrSpace),
  target       (initTarget),
  MaxParamNum  (initMaxParamNum),
  BetaDelta    (false, MaxParamNum, 1),

  beta         (false, MaxParamNum, 1),
  absCriterion (NAN)
{
  ASSERT (space);
  ASSERT (target);
  ASSERT (contains <const Attr*> (space->list, target));
  

  beta. clear ();
  BetaDelta. clear ();
}



void ExprNumPredict::GetTargetPredictedDeriv (size_t         /*objNum*/,
                                                const Matrix &/*TrialBeta*/,
                                                Matrix       &/*TargetPredictedDeriv*/) const
{
  throw NeverCall ();
}



namespace
{

class FUNCTION_MULT_PREDICT_NUM : public FuncMult
{
  const ExprNumPredict* pred;
    // !=nullptr
  Matrix TargetPredictedDeriv;
public:
    
    
  FUNCTION_MULT_PREDICT_NUM (const ExprNumPredict* initPred);
  
  
  Real f (const Matrix &X);
    // Return: pred->absCriterion
  void GetGradient (const Matrix &X,
             			    Matrix       &Gradient);
  void GetHessian (const Matrix &X,
          			      Matrix       &Hessian);
    // Approximation. Marquardt Hessian, p.s.d.
};




FUNCTION_MULT_PREDICT_NUM::FUNCTION_MULT_PREDICT_NUM (const ExprNumPredict* initPred):
  FuncMult             (initPred->MaxParamNum),
  
  pred                  (initPred),
  TargetPredictedDeriv  (false, MaxArgNum, 1)
{
  ASSERT (pred);
}



Real FUNCTION_MULT_PREDICT_NUM::f (const Matrix &X)
{
  ASSERT (X. rowSize (false) == MaxArgNum);
  
 
  Real S = 0.0;
  for (Iterator it (pred->sample); it ();)  
    {
      const Real TargetPredicted = pred->GetTargetPredicted (*it, X);
      if (isNan (TargetPredicted))
        {
          printf ("\n%s: \n", it->name);
          X. print (true, false, 0);
          ERROR;
        }
      S += it->weight * sqr (pred->target [*it] - TargetPredicted);
    }
  ASSERT (! isNan (S));
    
    
  return S / pred->multSum;
}



void FUNCTION_MULT_PREDICT_NUM::GetGradient (const Matrix &X,
             			                             Matrix       &Gradient)
{
  ASSERT (X.        rowSize (false) == MaxArgNum);
  ASSERT (Gradient. rowSize (false) == MaxArgNum);
  
 
  Gradient. putAll (0.0);

  for (Iterator it (pred->sample); it ();)  
    {
      const Real TargetPredicted = pred->GetTargetPredicted (*it, X);
      const Real Residual = it->weight * (pred->target [*it] - TargetPredicted);
      pred->GetTargetPredictedDeriv (*it, X, TargetPredictedDeriv);
      ASSERT (TargetPredictedDeriv. Defined ());
      Gradient. Add (false, Residual, TargetPredictedDeriv, false);
    }
  ASSERT (Gradient. Defined ());

  Gradient. putProdAll (- 2.0 / pred->multSum);
}



void FUNCTION_MULT_PREDICT_NUM::GetHessian (const Matrix &X,
                                   			      Matrix       &Hessian)
{
  ASSERT (X.        rowSize (false) == MaxArgNum);
  ASSERT (Hessian.  rowSize (false) == MaxArgNum);
  ASSERT (Hessian. isSquare ());
  
 
  Hessian. putAll (0.0);

  for (Iterator it (pred->sample); it ();)  
    {
      pred->GetTargetPredictedDeriv (*it, X, TargetPredictedDeriv);
      ASSERT (TargetPredictedDeriv. Defined ());
      Hessian. addVecVecT (false, TargetPredictedDeriv, true, 0, it->weight);
    }
  ASSERT (Hessian. Defined ());

  Hessian. putProdAll (2.0 / pred->multSum);
}

}



bool ExprNumPredict::Solve (size_t MaxIter)
{
  ASSERT (beta. Defined ());
  ASSERT (BetaDelta. Defined ());
  ASSERT (! isNan (multSum));


  FUNCTION_MULT_PREDICT_NUM Search (this);

  const bool Ok = Search. optimizeMarquardt (true, 0.0, beta, BetaDelta, (unsigned short) MaxIter);
  absCriterion = Search. f (beta);
  
  
  return Ok;
}
#endif



}
