// expr.cpp

#undef NDEBUG
#include "../common.inc"

#include "expr.hpp"



namespace Algebra_sp
{
using namespace std;
using namespace Common_sp;



namespace
{

Expr* cutOut (VectorOwn <Expr> &args,
	            size_t argNum) 
{ 
	const Expr* e = args [argNum];  
	args [argNum] = nullptr;
	return const_cast <Expr*> (e);
}

}




// Operation

#ifndef NDEBUG
namespace
{
  Vector <const Operation*> operations;
}
#endif



Operation::Operation (const string &name_arg,
						          size_t nArgs_arg)
: Named (name_arg)
, nArgs (nArgs_arg)
{ 
  ASSERT (! name. empty ());
	ASSERT (nArgs >= 1);  // ??
	operations << this;
}



void Operation::qc () const
{
  if (! qc_on)
    return;
	Named::qc ();
	IMPLY (rank (), prefix () || infix ());
	IMPLY (rank (), name. size () == 1);
}



bool Operation::leftAssociate (VectorOwn <Expr> &args) const
{
	ASSERT (args. size () == 2);
	
	OperExpr* e1 = const_cast <OperExpr*> (args [1] -> isOper (*this));
	if (! e1)
		return false;

	swap (args [0], args [1]);
	ASSERT (! e1->args [1] -> isOper (*this));
	swap (args [1], e1->args [1]);
	swap (e1->args [0], e1->args [1]);
	args [0] = e1->stndizeArgs ();

	return true;
}



bool Operation::orderAssociated (VectorOwn <Expr> &args) const
{
	ASSERT (args. size () == 2);

	// Ordering by simpler() 
  // Bubble sort
	if (OperExpr* e0 = const_cast <OperExpr*> (args [0] -> isOper (*this)))
  {		
		// e0->args are already ordered
		if (args [1] -> simpler (* e0->args [1]))
		{ 
			swap (args [1], e0->args [1]);
		  args [0] = e0->stndizeArgs ();
		  return true;
		}
	}
	else if (args [1] -> simpler (* args [0]))
	{
		swap (args [0], args [1]);  
	//return true;   // ??
  }
		
	return false;
}




Plus plus;
Minus minus;
Times times;
Over over;
Power power;




// Plus

namespace
{

struct TimesSeries : Root
{
private:
  const Expr* *e;
    // !0
    // !isNumber()
    // isOper(times) => !args[1]->isNumber()
public:
  

  explicit TimesSeries (const Expr* &e_arg)
    : e (& e_arg)
    { qc (); }
  void qc () const
    { if (! qc_on)
        return;
      ASSERT (e);
	  	ASSERT (*e);
    	ASSERT (! (*e)->isNumber ()); 
    }


	bool tailIsNumber () const
	  { if (const OperExpr* oe = (*e)->isOper (times))
	  	  return oe->args [0] -> isNumber ();
	  	return true;  // <=> Rational(1)
	  }
  void next ()
    { ASSERT (! tailIsNumber ());
    	const OperExpr* oe = (*e)->isOper (times);
    	ASSERT (oe);
      e = & const_cast <OperExpr*> (oe) -> args [0];
      qc ();
    }
  // Return: !0
	const Expr* getMultiplicand () const
	  { if (const OperExpr* oe = (*e)->isOper (times))
	  	  return oe->args [1];
	  	return *e;
	  }
	Expr* getNumber () 
	  { ASSERT (tailIsNumber ());
	  	if (const OperExpr* oe = (*e)->isOper (times))
	  	  return cutOut (const_cast <OperExpr*> (oe) -> args, 0);
	  	return new NatExpr (1);
	  }
	const Expr** getNumberArg () 
	  { ASSERT (tailIsNumber ());
	  	if (const OperExpr* oe = (*e)->isOper (times))
	  	  return & const_cast <OperExpr*> (oe) -> args [0];
	  	return e;
	  }
};

}



Expr* Plus::stndize (VectorOwn <Expr> &args) const
{ 
  ASSERT (args. size () == nArgs);

  if (leftAssociate (args))
		return stndize (args);

  if (orderAssociated (args)) 
		return stndize (args);

  Rational rat0, rat1;		

  if (   args [0] -> isRational (rat0)
  	  && ! rat0. numerator
  	 )
  	return cutOut (args, 1);

  // Compute numbers
  if (   args [0] -> isRational (rat0)
  	  && args [1] -> isRational (rat1)
  	 )
  {
  	rat0. add (rat1);
  	return rat0. getExpr ();
  }
  
  // n * x + m * x --> (n + m) * x  
  if (! args [0] -> isNumber ())
  {
  	const Expr* *e0 = & args [0];
  	if (const OperExpr* oe = args [0] -> isOper (plus))
  	  e0 = & const_cast <OperExpr*> (oe) -> args [1];
	  TimesSeries ts0 (*e0);
	  TimesSeries ts1 (args [1]);
    while (ts0. getMultiplicand () -> equal (* ts1. getMultiplicand ()))
      if (   ! ts0. tailIsNumber ()
          && ! ts1. tailIsNumber ()
         )
	    {
	    	ts0. next ();
	    	ts1. next ();
	    }
	    else
		    if (   ts0. tailIsNumber ()
		    	  && ts1. tailIsNumber ()
		    	 )
		    {
		    	const Expr* *e = ts0. getNumberArg ();
		    	ASSERT (e);
		    	ASSERT (*e);
		    	const bool eNumber = (*e)->isNumber ();
		    	Expr* pe = new OperExpr ( plus
		                              , ts0. getNumber ()
		                              , ts1. getNumber ()
		                              );
		    	*e = eNumber ? pe : new OperExpr (times, *e, pe);
		    	return cutOut (args, 0) -> stndize ();
		    }
		    else
		    	break;
	}

	return nullptr;
}




// Minus

Expr* Minus::stndize (VectorOwn <Expr> &args) const
{
  ASSERT (args. size () == nArgs);

	if (OperExpr* e0 = const_cast <OperExpr*> (args [0] -> isOper (*this)))
		return cutOut (e0->args, 0);
		
	if (OperExpr* e0 = const_cast <OperExpr*> (args [0] -> isOper (plus)))
		return (new OperExpr ( plus
							           , (new OperExpr (*this, cutOut (e0->args, 0))) -> stndizeArgs ()
							           , (new OperExpr (*this, cutOut (e0->args, 1))) -> stndizeArgs ()
							           )
					 ) -> stndizeArgs ();

  Rational rat;
  if (args [0] -> isRational (rat))
  {
  	if (! rat. numerator)
  	  return new NatExpr (0);
  }
  else
		return (new OperExpr ( times
							           , new OperExpr (minus, new NatExpr (1))
							           , cutOut (args, 0)
							           )
					 ) -> stndizeArgs ();
  		
	return nullptr;
}




// Times

Expr* Times::stndize (VectorOwn <Expr> &args) const
{
  ASSERT (args. size () == nArgs);

  // Move plus out of args[] up 
  if (args [0] -> isOper (plus))
  	swap (args [0], args [1]);
  if (OperExpr* sum = const_cast <OperExpr*> (args [1] -> isOper (plus)))
  {
  	const Expr* a0 = sum->args [0];
  	const Expr* a1 = sum->args [1];
  	sum->args [0] = (new OperExpr (times, args [0] -> copy (), a0)) -> stndizeArgs ();
  	sum->args [1] = (new OperExpr (times, args [0],            a1)) -> stndizeArgs ();  // args[0] may be delete'd
  	cutOut (args, 0);
  	cutOut (args, 1);
  	return sum->stndizeArgs ();
  }
  
  if (leftAssociate (args))
		return stndize (args);
  if (orderAssociated (args))
		return stndize (args);

  Rational rat0, rat1;

  if (args [0] -> isRational (rat0))
  {
  	if (! rat0. numerator)
  	  return new NatExpr (0);
    if (rat0. isOne ())
  	  return cutOut (args, 1);
  }

  // Compute
  if (   args [0] -> isRational (rat0)
  	  && args [1] -> isRational (rat1)
  	 )
  {
  	rat0. multiply (rat1);
  	return rat0. getExpr ();
  }
  
  // args[]->isOper(power): cf. Plus ??
  	
	return nullptr;
}




// Over

Expr* Over::stndize (VectorOwn <Expr> &args) const
{
  ASSERT (args. size () == nArgs);

	if (OperExpr* e0 = const_cast <OperExpr*> (args [0] -> isOper (*this)))
		return cutOut (e0->args, 0);
		
	if (OperExpr* e0 = const_cast <OperExpr*> (args [0] -> isOper (minus)))
		return (new OperExpr ( minus
		                     , (new OperExpr (*this, cutOut (e0->args, 0))) -> stndizeArgs ()
		                     )
		       ) -> stndizeArgs ();
		
	if (OperExpr* e0 = const_cast <OperExpr*> (args [0] -> isOper (times)))
		return (new OperExpr ( times
							           , (new OperExpr (*this, cutOut (e0->args, 0))) -> stndizeArgs ()
							           , (new OperExpr (*this, cutOut (e0->args, 1))) -> stndizeArgs ()
							           )
					 ) -> stndizeArgs ();

  Rational rat;
  if (args [0] -> isRational (rat))
  {
  	ASSERT (rat. denominator == 1);
  	ASSERT (! rat. negative);
  	if (rat. numerator == 1)
  	  return cutOut (args, 0);
  }
  		
	return nullptr;
}




// Power

Expr* Power::stndize (VectorOwn <Expr> &args) const
{
  ASSERT (args. size () == nArgs);

  // ??

  Rational rat0, rat1;

  if (args [0] -> isRational (rat0))
  {
	  if (args [1] -> isRational (rat1))
	  {
		  // Compute  	
	  	rat0. power (rat1);
	  	return rat0. getExpr ();
	  }
	  else
	  {
	  	if (! rat0. numerator)
	 	    return cutOut (args, 0);  
	    if (rat0. isOne ())
	  	  return cutOut (args, 0);
	  }
	}
  if (args [1] -> isRational (rat1))
  	if (! rat1. numerator)  
  		return new NatExpr (1);
  		
	return nullptr;
}




// Rational

Rational::Rational (uint numerator_arg,
									  uint denominator_arg,
									  bool negative_arg)
: numerator (numerator_arg)
, denominator (denominator_arg)
, negative (negative_arg)
{
  ASSERT (denominator);
  IMPLY (! numerator, ! negative);  
  simplify ();
}



Expr* Rational::getExpr () const
{
  Expr* e = new NatExpr (numerator);
	if (numerator && negative)
		e = new OperExpr (minus, e);
	if (denominator != 1)
		e = new OperExpr (times, e, new OperExpr (over, new NatExpr (denominator)));
	return e;
}



void Rational::add (const Rational &rat)
{
	uint a =      denominator;
	uint b = rat. denominator;
	for (;;)
	{
	  const uint c = gcd (a, b);
	  ASSERT (c);
	  if (c == 1)
	  	break;
	  a /= c;
	  b /= c;
	}
	
	Rational rat1 (rat);
			  numerator   *= b;
			  denominator *= b;
	rat1. numerator   *= a;
	rat1. denominator *= a;
	ASSERT (denominator == rat1. denominator);	

	if (negative == rat1. negative)
		numerator += rat1. numerator;  // Overflow ??
	else
		if (numerator >= rat1. numerator)
			numerator -= rat1. numerator;
		else
		{
			numerator = rat1. numerator - numerator;
			toggle (negative);
		}
		
  simplify ();
}



void Rational::multiply (const Rational &rat)
{
  // Overflow ??
  numerator   *= rat. numerator;
  denominator *= rat. denominator;
  if (rat. negative)
  	toggle (negative);

  simplify ();
}



void Rational::power (const Rational &rat)
{
	ASSERT (rat. denominator == 1);  // ??

	Rational rat1 (rat);
	rat1. simplify ();

	if (rat1. negative)
		swap (numerator, denominator);

	simplify ();

  // Overflow ??
	numerator   = powInt (numerator,   rat1. numerator);
	denominator = powInt (denominator, rat1. numerator);
	if (divisible (rat1. numerator, 2))
		negative = false;
}



void Rational::simplify ()
{
	for (;;)
	{
		const uint c = gcd (numerator, denominator);
		ASSERT (c);
		if (c == 1)
			break;
		numerator   /= c;
		denominator /= c;
	}
}




// Expr

bool Expr::simpler (const Expr &expr) const 
{ 
	if (isNumber ())
		if (expr. isNumber ())
			;
		else
			return true;
	else
		if (expr. isNumber ())
			return false;

  return simpler_ (expr);
}



const OperExpr* Expr::isOper (const Operation &oper) const
{ 
	if (const OperExpr* oe = asOperExpr ())
	  if (& oe->oper == & oper)
	  	return oe;
	return nullptr;
}



bool Expr::isRational (Rational &rat) const
{
	uint numerator, denominator;
	bool negative;
	
	if (getNumerator (numerator, negative))
	{
		rat = Rational (numerator, 1, negative);
		return true;
	}

  if (getDenominator (denominator))
	{
		rat = Rational (1, denominator, false);
		return true;
	}

	if (const OperExpr* e = isOper (times))
  	if (e->args [0] -> getNumerator (numerator, negative))
		  if (e->args [1] -> getDenominator (denominator))
		  {
		    rat = Rational (numerator, denominator, negative);
		  	return true;
		  }

  rat. denominator = 0;
  return false;
}



bool Expr::getNumerator (uint &numerator,
	                       bool &negative) const
{
	if (const NatExpr* n = asNatExpr ())
	{
		numerator = n->num;
		negative = false;
		return true;
	}
	
	if (const OperExpr* e = isOper (minus))
		if (const NatExpr* n = e->args [0] -> asNatExpr ())
		{
			numerator = n->num;
			negative = true;
			return numerator > 0;
		}

  return false;
}
	


bool Expr::getDenominator (uint &denominator) const
{
	if (const OperExpr* e = isOper (over))
		if (const NatExpr* n = e->args [0] -> asNatExpr ())
		{
			denominator = n->num;
	  	if (! denominator)
	 		  throw logic_error ("Zero division");  
			return true;
		}
	
  return false;
}
	



// OperExpr

void OperExpr::qc () const
{ 
  if (! qc_on)
    return;
	ASSERT (oper. nArgs == args. size ());
	CONST_ITER (VectorOwn <Expr>, it, args)
  { 
  	ASSERT (*it);
	  (*it)->qc ();
	}
}



Expr* OperExpr::copy () const 
{ 
	OperExpr* e = new OperExpr (*this); 
	for (auto &arg : e->args)
	  arg = arg->copy ();
	  	  
	return e;
}



string OperExpr::str_ (size_t parentRank,
	                     bool arg0,
	                     const Operation* leftOper) const 
{ 
	const bool braces =    oper. rank () > parentRank
	                    || (   oper. rank () == parentRank
	                        && oper. infix ()
	                        && ! oper. associative () 
	                        && arg0
	                       );

  string s;
	if (oper. infix ())
		s =         args [0] -> str_ (oper. rank (), true, braces ? 0 : leftOper) 
	      + " " + args [1] -> str_ (oper. rank (), false, & oper);
	else
		if (oper. prefix ())  
		{
			s += leftOper 
			       ? leftOper->rank () == oper. rank ()
							 // leftOper is suppressed: "+ -" --> "-", "* /" --> "/"
			         ? oper. name
							 // "/" --> "1 /"
			         : leftOper->name + " " + oper. prefixName ()
			       : oper. prefixName ();
			s += " " + args [0] -> str_ (oper. rank () - 1, true, 0);
		}
		else
		{
			ASSERT (! braces);
			s += prepend (leftOper) + oper. name + " (";
			FOR (size_t, i, args. size ())
		  { 
		  	if (i)
			  	s += ", ";
			  s += args [i] -> str_ (SIZE_MAX, (bool) i, 0);
			} 	
			s += ")";
	}
	
	return (braces ? (prepend (leftOper) + "(") : "") + s + (braces ? ")" : "");
}



bool OperExpr::simpler_ (const Expr &expr) const 
{ 
	if (const OperExpr* oe = expr. asOperExpr ())
	{ 
		if (& oper == & oe->oper)
	  { 
	  	FOR_REV (size_t, i, args. size ())
	  	  if (args [i] -> simpler (* oe->args [i]))
	  	  	return true;
	  	  else if (oe->args [i] -> simpler (* args [i]))
	  	  	return false;
	  	return false;
	  }
	  else
	    return oper. complexity () < oe->oper. complexity ();
	}
	return false;

} 



Expr* OperExpr::stndize () 
{ 
	Offset ofs;
  if (verbose ())
  {
		ofs. newLn (cout);
  	cout << "Before: " << str ();
  }
  	
	ITER (VectorOwn <Expr>, it, args)
	  *it = const_cast <Expr*> (*it) -> stndize ();
	Expr* res = stndizeArgs ();
	
  if (verbose ())
  {
  	ofs. newLn (cout);
  	cout << "After:  " << res->str ();
  }
  	
	return res;
}



bool OperExpr::isNumber () const
{
	CONST_ITER (VectorOwn <Expr>, it, args)
	  if (! (*it)->isNumber ())
	  	return false;
	 return true;
}



Expr* OperExpr::stndizeArgs ()
{
	if (Expr* expr = oper. stndize (args))
	{
		ASSERT (expr != this);
		delete this;
		expr->qc ();
	  return expr;
	}
	return this;
}




// VarExpr

VarExpr::VarExpr (const string &name_arg)
: name (name_arg)
{ 
	ASSERT (! name. empty ()); 
	ASSERT (! strCountSet (name, " "));
}




//

#ifndef NDEBUG
namespace
{

bool qc ()
{
	map <size_t/*complexity*/, const Operation*> complexity2oper;
	CONST_ITER (Vector <const Operation*>, it, operations)
	{
		ASSERT (! complexity2oper [(*it)->complexity ()]);
		complexity2oper [(*it)->complexity ()] = *it;
		(*it)->qc ();  // !qc_on => always true ??
	}

  // For Rational	
	ASSERT (minus. complexity () < over. complexity ());
	
  return true;
}


  const bool qc_ = qc ();
}
#endif



}
