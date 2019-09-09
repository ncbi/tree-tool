// language.hpp
// Natural language

#ifndef LANGUAGE_HPP
#define LANGUAGE_HPP

#include "../common.hpp"
using namespace Common_sp;



namespace Lang_sp
{


// SExpr
struct SExprNamed;
  struct SExprGeneral;   // text
  struct SExprContent;   // text
  struct SExprVariable;
struct SExprEndOfList;
struct SExprTruncated;   // text

typedef  map <string, const SExprGeneral*>  Var2SExpr;
  // !nullptr
typedef  map <string, string/*SExprGeneral::name*/>  Var2name;

struct Transformation;



struct SExpr : Root
{
protected:
  SExpr () = default;
public:
  virtual SExpr* copy () const = 0;


  virtual const SExprNamed* asSExprNamed () const
    { return nullptr; }
  virtual const SExprGeneral* asSExprGeneral () const
    { return nullptr; }
  virtual const SExprContent* asSExprContent () const
    { return nullptr; }
  virtual const SExprVariable* asSExprVariable () const
    { return nullptr; }
  virtual const SExprEndOfList* asSExprEndOfList () const
    { return nullptr; }
  virtual const SExprTruncated* asSExprTruncated () const
    { return nullptr; }

protected:
  static SExpr* parse (TokenInput &ti);
    // Return: nullptr <=> EOF
    //         descendants: SExprNamed
public:

  virtual bool unify (const SExprNamed &/*pattern*/,
                      Var2SExpr &/*var2SExpr*/) const 
    { throw logic_error ("Unifying a bad SExpr"); }
    // Return: success
    // Output: var2SExpr::second
    // Requires: *this has no SExprVariable's
  bool equal (const SExprNamed &pattern) const
    { Var2SExpr var2SExpr;
      return unify (pattern, var2SExpr);
    }
};



struct SExprNamed : SExpr
{
  string name;


protected:
  explicit SExprNamed (const string &name_arg)
    : name (name_arg)
    {}
public:
  void qc () const override;
  void saveText (ostream &os) const override
    { os << name; }


  const SExprNamed* asSExprNamed () const final
    { return this; }

  virtual Var2name getVarNames (const string &/*parentName*/) const 
    { return Var2name (); }
  virtual const SExprNamed* copySubstitute (const Var2SExpr &var2SExpr) const = 0;
    // Return: !nullptr
};



struct SExprGeneral : SExprNamed
{
  VectorPtr<SExpr> children;
    // !nullptr
    // Tree
protected:
  mutable bool referred {false};
public:


  SExprGeneral (const string &name_arg,
                size_t childrenNum) 
    : SExprNamed (name_arg)
    { children. resize (childrenNum, nullptr); }
 ~SExprGeneral ()
    { deleteChildren (); }
private:
  void deleteChildren ();
public:
  void qc () const override;
  void saveText (ostream &os) const override;
  SExprGeneral* copy () const final;


  const SExprGeneral* asSExprGeneral () const final
    { return this; }

  static SExprGeneral* parse (TokenInput &ti)
    { if (SExpr* s = SExpr::parse (ti))
        if (const SExprGeneral* sg = s->asSExprGeneral ())
          return var_cast (sg);
        else
          throw logic_error ("Parsing non-SExprGeneral");
      return nullptr;
    }

  Var2name getVarNames (const string &/*parentName*/) const final;
  const SExprGeneral* copySubstitute (const Var2SExpr &var2SExpr) const final;
  bool unify (const SExprNamed &pattern,
              Var2SExpr &var2SExpr) const final;
  bool apply (const Transformation &tr);
    // Return: success
  void applyDown (const Transformation &tr);
};



struct SFalse : SExprGeneral
  { SFalse () : SExprGeneral ("false", 0) {} };

struct STrue : SExprGeneral
  { STrue () : SExprGeneral ("true", 0) {} };



struct SExprContent : SExprNamed
// name: arbitrary text
{
  explicit SExprContent (const string &name_arg) 
    : SExprNamed (name_arg)
    {}
  void saveText (ostream &os) const override
    { os << strQuote (name); }
  SExprContent* copy () const final
    { return new SExprContent (name); }


  const SExprContent* asSExprContent () const final
    { return this; }

  const SExprContent* copySubstitute (const Var2SExpr &/*var2SExpr*/) const final
    { return new SExprContent (name); }
  bool unify (const SExprNamed &pattern,
              Var2SExpr &var2SExpr) const final
    { return    pattern. asSExprContent () 
             && name == pattern. name; 
    }
};



struct SExprVariable : SExprNamed
{
  explicit SExprVariable (const string &name_arg) 
    : SExprNamed (name_arg)
    {}
  SExprVariable* copy () const final
    { return new SExprVariable (name); }


  const SExprVariable* asSExprVariable () const final
    { return this; }

  Var2name getVarNames (const string &parentName) const final
    { Var2name var2name;
      var2name [name] = parentName;
      return var2name; 
    }
  const SExprVariable* copySubstitute (const Var2SExpr &/*var2SExpr*/) const final
    { throw logic_error ("SExprVariable::copySubstitute()"); }
};



struct SExprEndOfList : SExpr
// Auxiliary
{ 
  SExprEndOfList () = default;
  SExprEndOfList* copy () const final
    { throw logic_error ("SExprEndOfList::copy()"); }


  const SExprEndOfList* asSExprEndOfList () const final
    { return this; }
};



struct SExprTruncated : SExpr
// To be replaced by a non-SExprTruncated SExpr
{ 
  SExprTruncated () = default;
  void saveText (ostream &os) const override
    { os << "__truncated__"; }
  SExprTruncated* copy () const final
    { throw logic_error ("SExprTruncated::copy()"); }


  const SExprTruncated* asSExprTruncated () const final
    { return this; }

  bool unify (const SExprNamed &/*pattern*/,
              Var2SExpr &/*var2SExpr*/) const final
    { return false; }
};




// Transformation

struct Transformation : Root
{
  const SExprGeneral* lhs {nullptr};
  const SExprGeneral* rhs {nullptr};
  // !nullptr


  explicit Transformation (TokenInput &ti);
  void qc () const override;
  void saveText (ostream &os) const override
    { lhs->saveText (os);
      os << " -> ";
      rhs->saveText (os);
      os << ';' << endl;
    }
  bool empty () const final
    { return ! lhs && ! rhs; }
};




// Unicode

typedef  uint  Codepoint;
  // 0 - non-existing code point



struct Utf8 
{
private:
  ifstream is;
public:

  explicit Utf8 (const string &fName);

  bool get (Codepoint &v);
    // Return: false <=> EOF
    // Output: v
private:
  static size_t nextBytes (char &c);
    // Update: c - leading byte
};



struct Language : Root
/* 
  CF-grammar:
    char_star -> char char_star | nil
    char -> delimiter | cap_letter
    delimiter -> content
    cap_letter -> letter capital  # bold underscore italics ...
    letter -> content
    capital -> false | true
  text = char_star where first and last char's are delimiter('#')
*/
{
  // delimiters
  Vector<Codepoint> delimiters;
  StringVector delimiterNames;

  // digits
  Vector<Codepoint> digits;
  StringVector digitNames;

  // Letters
  Codepoint startCapital;
  Codepoint startSmall;
  size_t size;
  Vector<Codepoint> extraCapital;
  Vector<Codepoint> extraSmall;
  //
  StringVector letterNames;
    // Of small letters

  Vector<Transformation> trs;


protected:
  Language (Codepoint startCapital_arg,
            Codepoint startSmall_arg,
            size_t size_arg);
public:
  void qc () const override;


private:
  size_t capitalShift () const
    { return startSmall - startCapital; }
  size_t extraCapitalShift () const
    { return extraCapital. empty () ? 0 : (extraSmall [0] - extraCapital [0]); }
  Codepoint toSmall (Codepoint c) const;
    // Return: 0 <=> not in Alphabet
    //         c is capital => Return != c
  const string& small2name (Codepoint c) const;
public:
  SExprGeneral* codepoint2SExpr (Codepoint c) const;
    // Return: !nullptr
    // Input: 0: text start/end
  SExprGeneral* utf8_2SExpr (Utf8 &text) const;
    // Return: !nullptr

  void readTransformations (const string &tfmFName);
    // Output: trs
  void transform (SExprGeneral& root) const
    { Progress prog;
      for (const Transformation& tr : trs)
      { prog ();
        root. applyDown (tr);
      }
    }
};



struct RussianAlphabet : Language  
{
  RussianAlphabet ()
    : Language (0x410, 0x430, 32)
    { extraCapital << 0x401;  // Yo
      extraSmall   << 0x451;  // yo
      letterNames << "a" << "b" << "v" << "g" << "d" << "je" << "zh" << "z" << "ji" << "iy" << "k" << "l" << "m" << "n" << "o" << "p"
            << "r" << "s" << "t" << "u" << "f" << "x" << "c" << "ch" << "sh" << "sch" << "'" << "i" << "j" << "e" << "ju" << "ja"
            << "yo";
    }
};



#if 0
struct Category : Named
{
  typedef VectorPtr<Category> Expansion;
    // And
    // !empty()
  Vector<Expansion> expansions;
    // Or
    // empty() <=> terminal

  explicit Category (const string &name_arg)
    : Named (name_arg)
    {}
};



struct Construction;



struct SubConstruction
{
  const Construction* cons {nullptr};
    // !nullptr
  Vector<SExpr<Category>> args;
    // Match cons->args
};



struct Choice
{
  SubConstruction lhs;

  typedef Vector<SubConstruction> Expansion;
  Vector<Expansion> expansions;
};



struct Construction : Named
{
  VectorPtr<Category> args;
//Vector<Address> ??
  Vector<Choice> choices;
    // empty() <=> terminal (<=> "phoneme")

  explicit Construction (const string &name_arg)
    : Named (name_arg)
    {}
};



struct Language : Root
{
  static constexpr char* arrowS {"->"};
  map <string/*Category::name*/, const Category*> categories;
  map <string/*Construction::name*/, const Construction*> constructions;
  const Construction* rootConstruction {nullptr};
    // !nullptr


  explicit Language (const string &fName);
private:
  void parseGrammar (LineInput &f,
                     uchar pass);
  void parseRule (const string &line,
                  uchar pass);
};
#endif


}  // namespace



#endif

