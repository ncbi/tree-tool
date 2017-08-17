// expr.hpp

#ifndef EXPR_HPP
#define EXPR_HPP

#include "../common.hpp"
using namespace Common_sp;



namespace Algebra_sp
{



struct Expr;



struct Operation : Named
// Function: Number, Number, ... -> Number
{
  size_t nArgs;
	  // >= 1
	

protected:
	Operation (const string &name_arg,
	           size_t nArgs_arg);
    // Descendants: singletons
public:
	void qc () const;


  virtual size_t complexity () const = 0;	
	  // Different for all descendants
	virtual Expr* stndize (VectorOwn <Expr> &/*args*/) const 
	  { return nullptr; }
	  // Standardize and simplify
	  // Return: !0 <=> OperExpr of args will be delete'd; stndize()'d 
	  // Update: args
	  // Requires: args.size() == nArgs
	  //           args: stndize()'d
protected:
	bool leftAssociate (VectorOwn <Expr> &args) const;
	bool orderAssociated (VectorOwn <Expr> &args) const;
	  // Requires: args: leftAssociate()'d
public:
	virtual size_t rank () const
	  { return 0; }
	bool prefix () const
	  { return rank () && nArgs == 1; }
	bool infix () const
	  { return rank () && nArgs == 2; }
	virtual string prefixName () const
	  { return name; }
	  // Valid if prefix()
	virtual bool associative () const
	  { return infix (); }
	  // Presupposition: infix()
};



struct Plus : Operation
{
	Plus () : Operation ("+", 2)  {}

  size_t complexity () const 
    { return 3; }
	Expr* stndize (VectorOwn <Expr> &args) const;
	size_t rank () const
	  { return 3; }
};
extern Plus plus;
		


struct Minus : Operation
{
	Minus () : Operation ("-", 1)  {}

  size_t complexity () const 
    { return 1; }
  Expr* stndize (VectorOwn <Expr> &args) const;
	size_t rank () const
	  { return 3; }
};
extern Minus minus;
		


struct Times : Operation
{
	Times () : Operation ("*", 2)  {}

  size_t complexity () const 
    { return 4; }
	Expr* stndize (VectorOwn <Expr> &args) const;
	size_t rank () const
	  { return 2; }
};
extern Times times;



struct Over : Operation 
{
//args[0] != NatExpr(0) ??

	Over () : Operation ("/", 1)  {}

  size_t complexity () const 
    { return 2; }
  Expr* stndize (VectorOwn <Expr> &args) const;
	size_t rank () const
	  { return 2; }
	string prefixName () const
	  { return "1 " + name; }
};
extern Over over;
		


struct Power : Operation
{
//args[0] != NatExpr(0) && args[1] != NatExpr(0) ??

	Power () : Operation ("^", 2)  {}

  size_t complexity () const 
    { return 5; }
	Expr* stndize (VectorOwn <Expr> &args) const;
	size_t rank () const
	  { return 1; }
	bool associative () const
	  { return false; }
};
extern Power power;



struct Rational  // Redundant ??
{
  uint numerator;
  uint denominator;
  bool negative;


  Rational (uint numerator_arg,
					  uint denominator_arg,
					  bool negative_arg);
	  // Invokes: simplify()
	Rational ()
	  : numerator (0)
	  , denominator (0)
	  , negative (false)
	  {}
	void print (ostream &os) const
	  { os << (negative ? '-' : ' ') << numerator << '/' << denominator; }
					  

  bool valid () const
    { return (bool) denominator; }
  Expr* getExpr () const;
    // Opposite to Expr::isRational()
    // Complete Expr: times (minus NatExpr, over NatExpr)
  void add (const Rational &rat);
  void multiply (const Rational &rat);
  void power (const Rational &rat);
  void simplify ();
  bool isOne () const
    { return numerator == 1 && denominator == 1 && ! negative; }
};

	

// Ordered by simpler() descending
//struct Expr
	struct OperExpr;
	struct VarExpr;
	struct NatExpr;

	
	
struct Expr : Root
{
protected:
	Expr ()
	  {}
public:
	virtual Expr* copy () const = 0;

	
	virtual const OperExpr* asOperExpr () const
	  { return 0; }
	virtual const VarExpr* asVarExpr () const
	  { return 0; }
	virtual const NatExpr* asNatExpr () const
	  { return 0; }

  string str () const
    { return str_ (SIZE_MAX, true, nullptr); }	
private:
	virtual string str_ (size_t parentRank,
	                     bool arg0,
	                     const Operation* leftOper) const = 0;
    // Input: arg0 <=> this == parent->arg[0]; valid if parent->infix()
    //        leftOper: may be 0
	friend struct OperExpr;
protected:
	static string prepend (const Operation* leftOper)
	  { return leftOper ? (leftOper->name + " ") : ""; }
public:
	bool simpler (const Expr &expr) const;
	  // Antisymmetric (=> irreflexive)
	  // isNumber() is simplest
	  // Invokes: simpler_()
private:
	virtual bool simpler_ (const Expr &expr) const = 0;
	  // Antisymmetric (=> irreflexive)
	  // Requires: this->isNumber() == expr.isNumber()
public:
	bool equal (const Expr &expr) const
	  { return    !       simpler (expr) 
	  	       && ! expr. simpler (*this);
	  }
	virtual Expr* stndize () 
	  { return this; }
	  // Standardize and simplify
	  // Return: !0
	  // Idempotent
	  // May invoke: delete this
	  // Time: Exponential (cf. CNF to DNF transformation)
	  /* +,*,-,/:
	     <sum> ::= <summand>*
	     <summand> ::= <rational> <prod>
	     <prod> ::= <multiplicand>*
	     <multiplicand> ::= <variable> | <operation>(<expr>*)  -- !isNumber()
	     
	     <integer> ::= [-] <natural>
	     <rational> ::= <integer> [/ <natural>]
	  */
	const OperExpr* isOper (const Operation &oper) const;
	  // Return: May be 0
	virtual bool isNumber () const = 0;
  bool isRational () const
    { Rational rat;
    	return isRational (rat);
    }
  bool isRational (Rational &rat) const;
    // Output: rat: valid if Return
    // Invokes: getNumerator(), getDenominator()
private:
	// Output: all: valid if Return
  bool getNumerator (uint &numerator,
	                   bool &negative) const;
  bool getDenominator (uint &denominator) const;
    // denominator > 0
};
	
	
	
struct OperExpr : Expr
{
	const Operation &oper;
	VectorOwn <Expr> args;
	  // size() = oper.nArgs
	  // !0
	  // Tree

  OperExpr (const Operation &oper_arg,
            const Expr* arg1)
    : oper (oper_arg)
    { args << arg1; }
  OperExpr (const Operation &oper_arg,
            const Expr* arg1,
            const Expr* arg2)
    : oper (oper_arg)
    { args << arg1 << arg2; }
  OperExpr (const Operation &oper_arg,
            const Expr* arg1,
            const Expr* arg2,
            const Expr* arg3)
    : oper (oper_arg)
    { args << arg1 << arg2 << arg3; }
  void qc () const;
	Expr* copy () const;


	const OperExpr* asOperExpr () const
	  { return this; }

private:
	string str_ (size_t parentRank,
               bool arg0,
	             const Operation* leftOper) const;
	bool simpler_ (const Expr &expr) const;
public:
	Expr* stndize ();
	bool isNumber () const;

  Expr* stndizeArgs ();
    // Return: !0
    // Invokes: oper.stndize(args)
    //          May: delete this
};
	
	

struct VarExpr : Expr
{
	string name;
	  // !empty()

	  
	explicit VarExpr (const string &name_arg);
	Expr* copy () const 
	  { return new VarExpr (*this); }


	const VarExpr* asVarExpr () const
	  { return this; }
	  
private:
	string str_ (size_t /*parentRank*/,
               bool /*arg0*/,
	             const Operation* leftOper) const
	  { return prepend (leftOper) + name; }
	bool simpler_ (const Expr &expr) const 
	  { if (expr. asOperExpr ())
	  	  return true;
	  	if (const VarExpr* v = expr. asVarExpr ())
	  		return name < v->name;
	  	return false;
	  } 
public:
	bool isNumber () const
	  { return false; }
};
	
	

struct NatExpr : Expr
{
	uint num;

	
	explicit NatExpr (uint num_arg)
	  : num (num_arg)
	  {}
	Expr* copy () const 
	  { return new NatExpr (*this); }


	const NatExpr* asNatExpr () const
	  { return this; }

private:
	string str_ (size_t /*parentRank*/,
               bool /*arg0*/,
	             const Operation* leftOper) const
	  { return prepend (leftOper) + toString (num); }
	bool simpler_ (const Expr &expr) const 
	  { if (expr. asOperExpr ())
	  	  return true;
	  	if (expr. asVarExpr ())
	  	  return true;
	  	if (const NatExpr* i = expr. asNatExpr ())
	  		return num < i->num;
	  	return false;
	  } 
public:
	bool isNumber () const
	  { return true; }
};
	


#if 0
??
struct Function 
{
	Vector <string> args;
	const Expr* expr;
	  // !0
};
#endif



}  
#endif



// ??
// To do:
// parsing
// A * A -> A^2

