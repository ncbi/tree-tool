// optim.hpp

#ifndef OPTIM_HPP_78785
#define OPTIM_HPP_78785


#include "../common.hpp"
using namespace Common_sp;
#include "numeric.hpp"
#include "matrix.hpp"



namespace DM_sp
{


bool practicallyOptimal (bool min,
                         Real PracticalOptimum,
                         Real Y);
  // Return: true if Y is practically optimal w.r.t. PracticalOptimum




//////////////////////// Func1 /////////////////////////

struct Func1 : Root
// Function of 1 argument
{
  virtual Real f (Real x) = 0;
    // Return: May be NaN

  Real findZero (Real x_min,
	               Real x_max,
	               Real precision);
	  // Return: x s.t. f(x)=0, x in [x_min, x_max]
	  // Requires: f is increasing
	  // Binary search
  Real findZeroPositive (Real x_init,
                         Real precision);
    // Requires: the domain of f() >= 0
    // Invokes: findZero()
  Real optimizeFibonacci (bool Min, 
         		      	      Real XMin,
        		       	      Real YXMin, 
        		       	      Real XMax, 
        		       	      Real YXMax, 
        		       	      Real dX);
    // Find a local extremum by Fibonacci search
    // If f is a unimodal function, extremum is between XMin and XMax,
    //   and f' is 0 only at extremum then the local extremum 
    //   is global
    // Return: x = arg LocalExtr_X f(x) within dX, x \in [XMin, XMax]
    //         If cannot optimize then XMin or XMax
    // Requires: XMin <= XMax
    //           YXMin|Max = f (XMin|Max)
    //           dX > 0
    //           f(x) != NaN for x \in [XMin, XMax]
    // # iterations = log ((XMax - XMin) / dX) / log (Psi), 
    //   where Psi = Golden Section ((sqrt (5) - 1) / 2)
  bool minimizeConvex (Real x_min,
                			 Real x_max,
                			 Real x_init,
                			 Real dX,
                			 unsigned short iter_max,
                			 Real &x);
    // Difference from optimizeFibonacci ()  ??
    // Return: success
    // Output: x = arg min_x f(x), x \in [x_min, x_max]
    // Return: true if x is minimum within dX
    // Requires: x_min <= x_max
    //           dX > 0
    //           f(x) != NaN for x \in [x_min, x_max]
    // If f is convex then local minimum is global
};




///////////////////////// FuncMult ///////////////////////

struct FuncMult;  



struct MoveFunc : Func1, Nocopy
// Optimization of FuncMult in 1 dimension in the neighborhood of MFX0
{
  FuncMult* const MF;
    // != NULL
  const MVector* const MFX0;
    // != NULL
    // Requires: Size = (MF->maxArgNum, 1)
  const MVector* const MFGradient;
    // != NULL
    // Requires: Size = (MF->maxArgNum, 1)
  MVector MFX;
    // x for MF->f 
    // Size = (MF->maxArgNum, 1)


protected:
  MoveFunc (FuncMult* initMF,
            const MVector* initMFX0,
            const MVector* initMFGradient);
public:


  Real f (Real x);
    // Output: MFX
    // Return: MF->f (MFX); may be NaN
    // Requires: x >= 0.0; f(0.0) != NaN
    // Invokes: SetMFX (x)
  virtual bool SetMFX (Real x) = 0;
    // Output: MFX
    // Return: true if success
    // Requires: x >= 0.0; x = 0.0 => MFX = MFX0
};



struct FuncMult : Root
// Function of many arguments
// Newton method -??
// Constraints -??
// Find a global extremum ("Simulated annealing", "Genetic algorithm") -??
{
  size_t maxArgNum;
    // # Arguments
    // > 0


protected:
  explicit FuncMult (size_t maxArgNum_arg);
public:

  // Requires: Size of x, dX, gradient = (maxArgNum, 1)

  virtual Real f (const MVector &x) = 0;
    // Return: May be NaN
  virtual void getGradient (const MVector &x,
                     			  MVector &gradient);
    // Output: gradient
  virtual void getHessian (const MVector &x,
                  			   Matrix &hessian);
    // Output: hessian; Size = (maxArgNum, maxArgNum); Symmetric
    // If hessian is p.s.d. for all x, then f is convex 
    // Quality of optimization ??


protected:
  Real optimizeLinear (bool       Min,
                       Real       PracticalOptimum,
          				     MoveFunc*  MF,
          				     const Real dK,
          				     bool       Rough);
    // Optimization in the neighborhood of MF->MFX0
    // Input: dK - corresponds to dX[]
    // Return: arg LocalExtr_K MF->f(K), s.t. K >= 0.0;
    //         very large K if MF->f(K) is improving
    // Invokes: MF->optimizeFibonacci()
    // Requires: MF->f(K) != NaN in the neighborhood of MF->MFX0
public:


  // Find a local extremum
  // If f is convex|concave then local extremum is global
  // Update: x - initial point and the resulting extremum

  bool optimizeGradient (bool         Min,
                         Real         PracticalOptimum,
                     		 MVector       &x,
                      	 const MVector &dX,
                     		 unsigned short MaxIter);
    // Gradient method to find a local extremum (method of steepest ascent|descent)
    // Input: dX[] - precision of gradient; > 0.0
    // Return: true if converged
    // Invokes: getGradient(), GetGradientMove()
private:
  Real GetGradientMove (bool         Min,
                        Real         PracticalOptimum,
                     		const MVector &x,
                     		const MVector &Gradient,
                     		const Real   dK);
    // x + K * Gradient -> Min, s.t. K >= 0
    // Return: K
    // Invokes: optimizeLinear(not Rough)
public:

  bool optimizeMarquardt (bool           Min,
                          Real           PracticalOptimum,
                					MVector         &x,
                					const MVector   &dX,
                					unsigned short MaxIter);
    // Join with optimizeGradient() -??
    // Marquardt method to find a local extremum
    // Input: dX[] - precision of x; > 0.0
    // Return: true if converged
    // Update: x
    // Invokes: getGradient(), getHessian(), GetMarquardtX()
private:
  bool GetMarquardtX (bool         Min,
                      Real         PracticalOptimum,
                      const MVector &X0, 
                      const MVector &Gradient, 
                      const Matrix &Hessian, 
            					const MVector &dX,
                      MVector       &XNew);
    // Output: XNew
    // Return: true if success
    // Invokes: optimizeLinear(Rough)
    // Requires: dX [] > 0.0; f(X0) != NaN
public:

  bool optimizeHypercube (bool           Min,
                   				MVector         &x,
                   				const MVector   &XStep,
                   				unsigned short MaxIter,
                   				Real           &YBest);	
    // Hypercube search
    // Hypercube = x[i] + {0, XStep[i], -XStep[i]}
    // Input: XStep; > 0
    // Output: YBest = f(x)
    // Return: true if converged
  bool optimizeHypercubeDepth (bool           Min,
                      				 MVector         &x,
                      				 const MVector   &XStep,
                      				 Real           XStepCoeff,
                      				 unsigned short MaxIter,
                      				 unsigned short MaxDepth,
                      				 Real           &YBest);
    // Invokes: optimizeHypercube () MaxDepth times 
    //   multiplying XStep by XStepCoeff each time
    // Update: x
    // Output: YBest
    // Requires: 0 < XStepCoeff < 1
};


}



#endif

