// optim.cpp

#undef NDEBUG
#include "../common.inc"

#include "optim.hpp"



namespace DM_sp
{



bool practicallyOptimal (bool  Min,
                         Real PracticalOptimum,
                         Real Y)
{
  ASSERT (! isNan (PracticalOptimum));
  
  if (isNan (Y))
    return false;

  if (Min)
    return Y <= PracticalOptimum;
  else
    return Y >= PracticalOptimum;
}




////////////////////////////// Func1 ////////////////////////////////////

Real Func1::findZero (Real x_min,
                      Real x_max,
                      Real precision) 
{
	ASSERT (positive (precision));	
	ASSERT (! positive (f (x_min)));
	ASSERT (! negative (f (x_max)));
	
	Unverbose unv;
	
	for (;;)
	{
	  ASSERT (leReal (x_min, x_max));
  	const Real x = (x_min + x_max) / 2;
		if (x_max - x_min <= precision)
			return x;
  	const Real y = f (x);
  	ASSERT (! isNan (y));
    if (verbose ())
      cout << "x =" << x << "  y = " << y << endl;  
  	if (fabs (y) <= precision)  
  		return x;
  	if (y < 0)
  		x_min = x;
  	else
  		x_max = x;
  }
}



Real Func1::findZeroPositive (Real x_init,
                              Real precision) 
{
	ASSERT (x_init > 0);
	
  Real x_min = x_init;
  while (positive (f (x_min)))
    x_min /= 2;
  if (x_min != x_init)
  	return findZero (x_min, 2 * x_min, precision);
  if (verbose ())
    cout << "x_init=" << x_init << "  x_min=" << x_min << endl;  
    
  Real x_max = x_init;
  while (negative (f (x_max)))
    x_max *= 2;
  if (x_max != x_init)
  	return findZero (x_max / 2, x_max, precision);
  if (verbose ())
    cout << "x_init=" << x_init << "  x_max=" << x_max << endl;   
    
  return findZero (x_min, x_max, precision);
}



Real Func1::optimizeFibonacci (bool Min, 
                               Real x_min, 
                               Real YXMin, 
                               Real x_max, 
                               Real YXMax, 
                               Real dX)
{
  ASSERT (! isNan (YXMin));
  ASSERT (! isNan (YXMax));
  ASSERT (x_min <= x_max);
  

	Unverbose unv;


  if (x_min == x_max)
    return x_min;
    

  ASSERT (dX > 0.0);


  // Golden Section (= 0.618...)
  const Real Psi = (sqrt (5.0) - 1.0) / 2.0;

  
  // x_min <= A <= B <= x_max
  Real Range = x_max - x_min;
  Real A = x_max - Psi * Range;
  Real B = x_min + Psi * Range;
  Real YA = f (A);
  Real YB = f (B);
  ASSERT (! isNan (YA));
  ASSERT (! isNan (YB));
  Real Res = ((YXMin < YXMax) == Min) ? x_min : x_max;
  Real YRes = Min ? min (YXMin, YXMax) : max (YXMin, YXMax);
  const Real InitRes  = Res;
  const Real InitYRes = YRes;
  while (Range > dX)
  {
    ASSERT (leReal (x_min, A));
    ASSERT (leReal (A, B));
    ASSERT (leReal (B, x_max));
    
  /*
    printf ("x_min = %f,  YXMin = %f,  x_max = %f,  YXMax = %f\n",
      x_min, YXMin, x_max, YXMax);
  */

    if (Min)
      {
        // To allow non-convex|concave functions
      //ASSERT (leReal (max (YA, YB), max (YXMin, YXMax)));
        Res = (YA < YB) ? A : B;
        YRes = min (YA, YB);
      }
    else
      {
        // To allow non-convex|concave functions
      //ASSERT (geReal (min (YA, YB), min (YXMin, YXMax)));
        Res = (YA > YB) ? A : B;
        YRes = max (YA, YB);
      }
    
    if ((YA > YB) == Min)
      {
        x_min  = A;
        YXMin = YA;
        
        A  = B;
        YA = YB;
        
        Range = x_max - x_min;
        
        B = x_min + Psi * Range;
        YB = f (B);
        ASSERT (! isNan (YB));
      }
    else
      {
        x_max  = B;
        YXMax = YB;
        
        B  = A;
        YB = YA;
        
        Range = x_max - x_min;
       
        A = x_max - Psi * Range;
        YA = f (A);
        ASSERT (! isNan (YA));
      }
  }
    

  if ((InitYRes < YRes) == Min)
    Res = InitRes;
 

  return Res;
}



bool Func1::minimizeConvex (Real x_min,
                            Real x_max,
                			      Real x_init,
                            Real dX,
                            unsigned short iter_max,
                            Real &x)
{
  ASSERT (dX > 0.0);

  
	Unverbose unv;


  bool Ok = false;
  Real x1 = x_init;  
  Real f1 = f (x1);
  ASSERT (! isNan (f1));
  FOR (unsigned short, iter, iter_max)
  {
    ASSERT (x_min <= x_max);
    ASSERT (betweenReal (x1, x_min, x_max));
    if (verbose ())
    	cout << "f1=" << f1 << endl;

    if (x_max - x_min <= dX)
    {
      Ok = true;
      break;
    }

    // x2
    Real x2;
    if (x1 - x_min > x_max - x1)
      x2 = x_min + (x1   - x_min) / 2.0;  
    else
      x2 = x1   + (x_max - x1)   / 2.0;
    ASSERT (betweenReal (x2, x_min, x_max));
    
    // f2
    const Real f2 = f (x2);
    ASSERT (! isNan (f2));

    if (f1 > f2)
    {
      swap (x1, x2);
      f1 = f2;
    }
    // f(x1) <= f(x2)
    
    if (x1 < x2)
      x_max = x2;
    else
      x_min = x2;
  }
  x = x1;  
  
  
  return Ok;
}




//////////////////////////// FuncMult //////////////////////////////////

// MoveFunc

MoveFunc::MoveFunc (FuncMult* initMF,
                              const MVector*  initMFX0,
                              const MVector*  initMFGradient)
: MF         (initMF)
, MFX0       (initMFX0)
, MFGradient (initMFGradient)
, MFX        (MF->maxArgNum)
{
  ASSERT (MF);

  ASSERT (MFX0);
  ASSERT (MFX0->rowsSize (false) /*Len [false]*/ == MF->maxArgNum);
  ASSERT (MFX0->rowsSize (true)  /*Len [true]*/  == 1);

  ASSERT (MFGradient);
  ASSERT (MFGradient->rowsSize (false) /*Len [false]*/ == MF->maxArgNum);
  ASSERT (MFGradient->rowsSize (true)  /*Len [true]*/  == 1);
}



Real MoveFunc::f (Real X)
{ 
  ASSERT (! isNan (X));
  ASSERT (X >= 0.0); 
  

  if (SetMFX (X))
    return MF->f (MFX); 
  else
    return NaN;
}




// FuncMult

FuncMult::FuncMult (size_t maxArgNum_arg)
: maxArgNum (maxArgNum_arg)
{
  ASSERT (maxArgNum);
}



void FuncMult::getGradient (const MVector & /*X*/,
                            MVector & /*Gradient*/)
{
  NOT_IMPLEMENTED;
}



void FuncMult::getHessian (const MVector & /*X*/,
                           Matrix       & /*Hessian*/)
{
  NOT_IMPLEMENTED;
}



Real FuncMult::optimizeLinear (bool           Min,
                                Real          PracticalOptimum,
                                MoveFunc* MF,
                                const Real    dK,
                                bool           Rough)
{
  ASSERT (MF);
  ASSERT (dK > 0.0);
  ASSERT (finite (dK));


  Unverbose unv;

  // ??
  const Real dKCoeff = 10.0;
  const unsigned short  MaxIter = 1000;


  Real KPrevPrev = 0.0;
  Real YPrevPrev = MF->f (KPrevPrev);
  ASSERT (! isNan (YPrevPrev));
  ASSERT (finite (YPrevPrev));
  if (practicallyOptimal (Min, PracticalOptimum, YPrevPrev))
    return KPrevPrev;
  
  Real KPrev = KPrevPrev;
  Real YPrev = YPrevPrev;
  
  Real K = dK * dKCoeff;
  Real Y = NaN;
  while (isNan (Y))
  {
    K /= dKCoeff;
    Y = MF->f (K);
    if (practicallyOptimal (Min, PracticalOptimum, Y))
      return K;
  }
  
  if (verbose ())      
    cerr << endl << "  Y-YPrev=" << Y - YPrev << "  K=" << K << endl;

  unsigned short iter;
  for (iter = 0; iter < MaxIter; iter++)
  {
    if ((YPrev < Y) == Min)
      break;
   
    KPrevPrev = KPrev;
    YPrevPrev = YPrev;

    KPrev = K;
    YPrev = Y;

    K *= dKCoeff;
    Y = MF->f (K);
    if (isNan (Y))
      return KPrev;
    if (practicallyOptimal (Min, PracticalOptimum, Y))
      return K;

    if (verbose ())      
      cerr << "  Y=" << Y << "  K=" << K << endl;
  }
  ASSERT (finite (YPrev));
  
  
  if (Rough)
    if (iter > 2)  // ??
      return KPrev;
  

  if (! finite (Y))
  {
    Real KL = KPrev;
    Real KU = K;
    while (KU - KL > dK)
    {
      const Real M = (KU + KL) / 2.0;
      const Real YM = MF->f (M);
      if (verbose ())      
        cerr << "  YM=" << YM << "  M=" << M << "  KL=" << KL << "  KU=" << KU << endl;
      ASSERT (! isNan (YM));
      if (finite (YM))
        KL = M;
      else
        KU = M;
    }
    K = KL;
    Y = MF->f (K);
    if (verbose ())      
      cerr << "  Y=" << Y << "  K=" << K << endl;
  }
  ASSERT (! isNan (Y));
  ASSERT (finite (Y));
  if (practicallyOptimal (Min, PracticalOptimum, Y))
    return K;


  const Real KDiff = K - KPrevPrev;
  ASSERT (KDiff >= 0.0);
  return MF->optimizeFibonacci (Min, KPrevPrev, YPrevPrev, K, Y, 
                                max (max (min (dK, KDiff / 10.0), KDiff / 1e6), 1e-20));
}




// optimizeGradient

bool FuncMult::optimizeGradient (bool         Min,
                                 Real        PracticalOptimum,
                                 MVector       &X,
                                 const MVector &dX,
                                 unsigned short         MaxIter)
{
  ASSERT (X.  rowsSize (false) /*Len [false]*/ == maxArgNum);
  ASSERT (dX. rowsSize (false) /*Len [false]*/ == maxArgNum);
  ASSERT (dX. min () >= 0.0);


  Unverbose unv;

  MVector Gradient (maxArgNum);


  bool Stop = false;
  unsigned short iter;
  for (iter = 0; iter < MaxIter; iter++)
  {
    if (verbose ())
      cerr << "\rGradient: " << iter + 1 << " / " << MaxIter << ": " << f (X) << " ";

    if (Stop)
      break;

    // X -> Gradient
    getGradient (X, Gradient);
    if (nullReal (Gradient. maxAbs ()))
      break;
    if (Min)
      Gradient. negate ();

    // dX, Gradient -> dK
    Real dK = INF;
    FOR (unsigned short, ArgNum, maxArgNum)
      minimize (dK, dX [ArgNum] / fabs (Gradient [ArgNum]));
    ASSERT (finite (dK));
    ASSERT (dK > 0.0);

    // X, Gradient, dK -> K
    const Real K = GetGradientMove (Min, PracticalOptimum, X, Gradient, dK);
    ASSERT (K >= 0.0);

    // Gradient, K -> Gradient, X
    X. add (false, Gradient, false, K);

    // Stop      
    Stop = true;
    FOR (unsigned short, ArgNum, maxArgNum)
      if (fabs (Gradient [ArgNum]) > dX [ArgNum])
        Stop = false; 

    if (verbose ())
    {
      cerr << endl << iter + 1 << ": Y = " << f (X) << "  K = " << K << ";  ";

      cerr << "X [] =";
      X. printRow (true, 0, cerr);

      cerr << ";  Gradient [] =";
      Gradient. printRow (true, 0, cerr);

      cerr << endl;
    }
  }
  if (verbose ())
    cerr << endl;
  

  return (iter < MaxIter);
}



struct GradientMoveFunc : public MoveFunc
{
  GradientMoveFunc (FuncMult* initMF,
                    const MVector*  initMFX0,
                    const MVector*  initMFGradient)
    : MoveFunc (initMF, initMFX0, initMFGradient) 
    {}


  bool SetMFX (Real X)
  // MFX = MFX0 + MFGradient * X
  {
    ASSERT (! isNan (X));
    ASSERT (X >= 0.0);
    
   
    MFX = *MFX0;
    MFX. add (false, *MFGradient, false, X);
    
    
    return true;
  }
};




Real FuncMult::GetGradientMove (bool         Min,
                                 Real        PracticalOptimum,
                                 const MVector &X,
                                 const MVector &Gradient,
                                 const Real  dK)
{
  ASSERT (X.        rowsSize (false) /*Len [false]*/ == maxArgNum);
  ASSERT (Gradient. rowsSize (false) /*Len [false]*/ == maxArgNum);
  ASSERT (finite (dK));
  ASSERT (dK > 0.0);
  
  GradientMoveFunc MF (this, & X, & Gradient);
  return optimizeLinear (Min, PracticalOptimum, & MF, dK, false);
}




// optimizeMarquardt

bool FuncMult::optimizeMarquardt (bool         Min,
                                  Real        PracticalOptimum,
                                  MVector       &X,
                                  const MVector &dX,
                                  unsigned short         MaxIter)
{
  ASSERT (X.  rowsSize (false) == maxArgNum);
  ASSERT (dX. rowsSize (false) == maxArgNum);
  ASSERT (positive (dX. min ()));


	Unverbose unv;


#if 0
  cout << "dX:" << endl;
  dX. saveText (cout);
#endif


  MVector Gradient (maxArgNum);
  Matrix Hessian (maxArgNum);
  MVector XNew (maxArgNum);
  bool Stop = false;
  unsigned short iter;
  Real y_prev = NaN;
  for (iter = 0; iter < MaxIter; iter++)
  {
    if (verbose ())
      cerr << "\rMarquardt: " << iter + 1 << " / " << MaxIter << ": " << f (X) << " ";

    if (Stop)
      break;
    const Real y = f (X);
    if (practicallyOptimal (Min, PracticalOptimum, y))
      break;
    if (eqReal (y, y_prev))
      break;
    y_prev = y;

    // XNew
    getGradient (X, Gradient);
    if (nullReal (Gradient. maxAbs ()))
      break;
    getHessian (X, Hessian);

  #if 0
    cout << "gradient:" << endl;
    Gradient. saveText (cout);
    cout << "Hessian:" << endl;
    Hessian. saveText (cout);
    cout << "dX:" << endl;
    dX. saveText (cout);
  #endif

    if (! GetMarquardtX (Min, PracticalOptimum, X, Gradient, Hessian, dX, XNew))
    {
      if (verbose ())
        cerr << endl << "Marquardt is inapplicable" << endl;
      return false;
    }

    // Stop      
    Stop = true;
    FOR (unsigned short, ArgNum, maxArgNum)
      if (fabs (XNew [ArgNum] - X [ArgNum]) > dX [ArgNum])
        Stop = false; 

    X = XNew;

    if (verbose ())
    {
      cerr << endl << iter + 1 << ": Y = " << f (X) << ";  ";

      cerr << ";  Gradient [] =";
      Gradient. printRow (true, 0, cerr);

      cerr << ";  X [] =";
      X. printRow (true, 0, cerr);

      cerr << endl;
    }
  }
  if (verbose ())
    cerr << endl;
  

  return (iter < MaxIter);
}



struct MarquardtMoveFunc : public MoveFunc
{
  const Matrix* const MFHessian;
    // Requires: Size = (MF->maxArgNum, MF->maxArgNum)
protected:
  Matrix Q;
    // Size = (MF->maxArgNum, MF->maxArgNum + 1)
public:


  MarquardtMoveFunc (FuncMult* initMF,
                     const MVector*  initMFX0,
                     const MVector*  initMFGradient,
                     const Matrix*  initMFHessian)
    : MoveFunc (initMF, initMFX0, initMFGradient)
    , MFHessian (initMFHessian)
    , Q         (false, MF->maxArgNum, MF->maxArgNum + 1)
    {
      ASSERT (MFHessian);
      ASSERT (MFHessian->isSquare ());
      ASSERT (MFHessian->rowsSize (false) /*Len [false]*/ == MF->maxArgNum);
    }


  bool SetMFX (Real X)
  // X -> 0.0: Linear (gradient) optimization
  // X = 0.0: MFX = MFX0
  // X -> INF: Quadratic optimization
  {
    ASSERT (! isNan (X));
    ASSERT (X >= 0.0);
    
    // delta MFX
    if (X == 0.0)
      MFX. putAll (0.0); 
    else
    {
      // Q
      Q. copyShift (false, 0, MF->maxArgNum, *MFGradient, false);
      Q. negateRow (true, MF->maxArgNum);
      
      Q. copyShift (false, 0, 0, *MFHessian, false);
      const Real Lambda = 1.0 / X;
      const Real x1 = 1.0 + Lambda;
      ASSERT (finite (x1));
      FOR (size_t, Row, Q. rowsSize (false) /*Len [false]*/)
        Q. putProd (false, Row, Row, x1);

      
      if (! MFX. solveSystem (false, 0, Q, false))
        return false;
    }

    MFX. add (false, *MFX0, false, 1);
      
    return true;
  }
};



bool FuncMult::GetMarquardtX (bool         Min,
                              Real        PracticalOptimum,
                              const MVector &X0, 
                              const MVector &Gradient, 
                              const Matrix &Hessian, 
                              const MVector &dX,
                              MVector       &XNew)
{
  ASSERT (X0.       rowsSize (false) /*Len [false]*/ == maxArgNum);
  ASSERT (Gradient. rowsSize (false) /*Len [false]*/ == maxArgNum);
  ASSERT (Hessian.  rowsSize (false) /*Len [false]*/ == maxArgNum);
  ASSERT (Hessian. isSquare ());
  ASSERT (dX.   rowsSize (false) /*Len [false]*/ == maxArgNum);
  ASSERT (XNew. rowsSize (false) /*Len [false]*/ == maxArgNum);
  
  
  Unverbose unv;


  MarquardtMoveFunc MF (this, & X0, & Gradient, & Hessian);


  const Real Y0 = f (X0);
  ASSERT (! isNan (Y0));
  const Real YNew = MF. f (INF);
  if (isNan (YNew) ||
      (YNew > Y0) == Min)
  {
    // dK
    Real A = 0.0;
    FOR (unsigned short, ArgNum, maxArgNum)
      maximize (A, fabs (Gradient [ArgNum] / (Hessian. getDiag (ArgNum) * dX [ArgNum])));
    const Real dK = 1.0 / A;

  //cout << "dK = " << dK << endl;

    if (! finite (dK))  // Marquardt method is inapplicable
      return false;  
      
    const Real K = optimizeLinear (Min, PracticalOptimum, & MF, dK, true);
    const Real Y = MF. f (K);
    ASSERT (! isNan (Y));
  }


  XNew = MF. MFX;
    
    
  return true;
}




// optimizeHypercube

bool FuncMult::optimizeHypercube (bool         Min,
                                  MVector       &X,
                                  const MVector &XStep,
                                  unsigned short         MaxIter,
                                  Real        &YBest)
{
  ASSERT (X.  rowsSize (false) /*Len [false]*/ == maxArgNum);
  ASSERT (XStep. rowsSize (false) /*Len [false]*/ == maxArgNum);
  ASSERT (XStep. min () > 0);
  ASSERT (MaxIter > 0);


	Unverbose unv;


  // XInc [ArgNum] = X [ArgNum] + XStep [ArgNum] * Inc [ArgNum]
  MVector XInc (maxArgNum);

  // = -1, 0, 1
  int* Inc = new int [maxArgNum];
  int* IncBest = new int [maxArgNum];
  // Increment to get to the previous point
  int* IncPrev = new int [maxArgNum];


  // Initializing YBest, IncBest []
  YBest = Min ? INF : -INF;
  FOR (unsigned short, ArgNum, maxArgNum)
    IncBest [ArgNum] = 0;


  if (verbose ())
    cerr << endl;
  bool Changed = false;
  FOR (unsigned short, iter, MaxIter)
  {
    if (verbose () &&
        iter > 0)
    {
      cerr << "\rHypercube: " << iter << ": Y = " << YBest << "   X [] =";
      X. printRow (true, 0, cerr);
      cerr << "   XStep [] =";
      XStep. printRow (true, 0, cerr);
      cerr << endl;
    }


    // Searching the hypercube: 3^maxArgNum
    // YBest, IncBest []
    FOR (unsigned short, ArgNum, maxArgNum)
      Inc [ArgNum] = -1;
    for (;;)
    {
      if (verbose ())
      {
        cerr << "\r";
        FOR (unsigned short, ArgNum, maxArgNum)
          cerr << (Inc [ArgNum] == -1 ? 'v' : 
                   Inc [ArgNum] ==  1 ? '^' : '-');
      }


      bool Visited;
      if (iter == 0)
        Visited = false;
      else
      {
        Visited = true;
        FOR (unsigned short, ArgNum, maxArgNum)
          if (abs (Inc [ArgNum] - IncPrev [ArgNum]) >= 2)
            {
              Visited = false;
              break;
            }
      }


      if (! Visited)
      {
        // Y = f (XInc)
        XInc = X;
        FOR (unsigned short, ArgNum, maxArgNum)
          XInc. putInc (false, ArgNum, 0, XStep [ArgNum] * Inc [ArgNum]);
        const Real Y = f (XInc);
        
        if (! isNan (Y) &&
            ! eqReal (Y, YBest) &&
            (Y < YBest) == Min)
        {
          YBest = Y;
          FOR (unsigned short, ArgNum, maxArgNum)
            IncBest [ArgNum] = Inc [ArgNum];
        }
      }


      // Incrementing Inc []
      {
        unsigned short ArgNum;
        for (ArgNum = 0; ArgNum < maxArgNum; ArgNum++)
          if (Inc [ArgNum] == 1)
            Inc [ArgNum] = -1;
          else
          {
            Inc [ArgNum] ++;
            break;
          }
        if (ArgNum == maxArgNum)
          break;
      }
    }
    
    
    Changed = false;
    FOR (unsigned short, ArgNum, maxArgNum)
    {
      IncPrev [ArgNum] = - IncBest [ArgNum];
      if (IncBest [ArgNum] != 0)
      {
        X. putInc (false, ArgNum, 0, XStep [ArgNum] * IncBest [ArgNum]);
        IncBest [ArgNum] = 0;
        Changed = true;
      }
    }
    if (! Changed)
      break;
  }
  if (verbose ())
    cerr << "\r";


  delete [] Inc;  
  delete [] IncBest;
  delete [] IncPrev;


  return ! Changed;
}



bool FuncMult::optimizeHypercubeDepth (bool Min,
                                       MVector &X,
                                       const MVector &XStep,
                                       Real XStepCoeff,
                                       unsigned short MaxIter,
                                       unsigned short MaxDepth,
                                       Real &YBest)
// Theory -??
{
  ASSERT (X.  rowsSize (false) == maxArgNum);
  ASSERT (XStep. rowsSize (false) == maxArgNum);
  ASSERT (XStep. min () > 0);
  ASSERT (XStepCoeff > 0);
  ASSERT (XStepCoeff < 1);
  ASSERT (MaxDepth > 0);


  // _XStep
  MVector _XStep (XStep);

/* ??
  Matrix XMod (false, X, false);
  Matrix XDiff (false, X, false);
*/
  
  FOR (unsigned short, Depth, MaxDepth)
  {
    if (! optimizeHypercube (Min, X, _XStep, MaxIter, YBest))
      return false;

  /* ??
    // XDiff
    XMod. CopyData (false, X, false);
    FOR (unsigned short, ArgNum, maxArgNum)
    {
      const Real Delta = _XStep [ArgNum];
      XMod. putInc (false, ArgNum, 0, Delta);
      const Real Y1 = f (XMod);
      XMod. putInc (false, ArgNum, 0, - 2 * Delta);
      const Real Y2 = f (XMod);
      XDiff. put (false, ArgNum, 0, max (fabs (Y1 - YBest), fabs (Y2 - YBest)));
    }
    XDiff. StandardizeRow (true, 0);
    FOR (unsigned short, ArgNum, maxArgNum)
      _XStep. putProd (false, ArgNum, 0, XStepCoeff * (0.5 + 0.5 * XDiff [ArgNum]));
  */
    _XStep. putProdAll (XStepCoeff);
      
    if (verbose ())
    {
      cerr << Depth + 1 << ": YBest = " << YBest << "  X [] =";
      X. printRow (true, 0, cerr);
      cerr << "  _XStep [] =";
      _XStep. printRow (true, 0, cerr);
      cerr << endl;
    }
  }
      

  return true;
}



}


