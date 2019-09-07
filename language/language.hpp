// language.hpp
// Natural language

#ifndef LANGUAGE_HPP
#define LANGUAGE_HPP

#include "../common.hpp"
using namespace Common_sp;



namespace Lang_sp
{


// SExpr
struct SExprGeneral;
struct SContent;
struct SVariable;
struct SEndOfList;
struct STruncated;



struct SExpr : Named
{
protected:
  explicit SExpr (const string &name_arg) 
    : Named (name_arg)
    {}
public:
  void qc () const override;
  void saveText (ostream &os) const override
    { os << name; }


  virtual const SExprGeneral* asSExprGeneral () const
    { return nullptr; }
  virtual const SContent* asSContent () const
    { return nullptr; }
  virtual const SVariable* asSVariable () const
    { return nullptr; }
  virtual const SEndOfList* asSEndOfList () const
    { return nullptr; }
  virtual const STruncated* asSTruncated () const
    { return nullptr; }

  static SExpr* parse (TokenInput &ti);
    // Return: nullptr <=> EOF

  virtual StringVector getVarNames () const 
    { return StringVector (); }

  typedef  map <string, const SExpr*>  Var2SExpr;
  virtual bool unify (const SExpr &pattern,
                        Var2SExpr &var2SExpr) const 
    { throw logic_error ("Matching a bad SExpr"); }
    // Return: success
    // Output: var2SExpr::second
    // Requires: *this has no SVariable's
  bool equal (const SExpr &pattern) const
    { Var2SExpr var2SExpr;
      return unify (pattern, var2SExpr);
    }
};



struct SExprGeneral : SExpr
{
  VectorOwn<SExpr> children;
    // !nullptr
    // Tree


  SExprGeneral (const string &name_arg,
                size_t childrenNum) 
    : SExpr (name_arg)
    { children. resize (childrenNum, nullptr); }
  void qc () const override;
  void saveText (ostream &os) const override
    { os << name;
      if (children. empty ())
        return;
      os << " (";
      bool first = true;
      for (const SExpr* child : children)
      { if (! first)
          os << ' ';
        if (child)
          child->saveText (os);
        else
          os << "nil";
        first = false;
      }
      os << ')';
    }


  const SExprGeneral* asSExprGeneral () const final
    { return this; }

  StringVector getVarNames () const final
    { StringVector s;
      for (const SExpr* child : children)
        s << child->getVarNames ();
      s. sort ();
      s. uniq ();
      return s;
    }
  bool unify (const SExpr &pattern,
                Var2SExpr &var2SExpr) const final;
};



struct SNil : SExprGeneral
  { SNil () : SExprGeneral ("nil", 0) {} };
  // Also repressnts "false"

struct STrue : SExprGeneral
  { STrue () : SExprGeneral ("true", 0) {} };



struct SContent : SExpr
// name: arbitrary text
{
  explicit SContent (const string &name_arg) 
    : SExpr (name_arg)
    {}

  const SContent* asSContent () const final
    { return this; }
  bool unify (const SExpr &pattern,
                Var2SExpr &var2SExpr) const final;
};



struct SVariable : SExpr
{
  explicit SVariable (const string &name_arg) 
    : SExpr (name_arg)
    {}

  const SVariable* asSVariable () const final
    { return this; }

  StringVector getVarNames () const final
    { return StringVector ({name}); }

  bool assign (const SExpr &target,
               Var2SExpr &var2SExpr) const;
};



struct SEndOfList : SExpr
// Auxiliary
{ 
  SEndOfList () 
    : SExpr ("end_of_list_") 
    {} 

  const SEndOfList* asSEndOfList () const final
    { return this; }
};



struct STruncated : SExpr
// To be replacde by a non-STruncated SExpr
{ 
  STruncated () 
    : SExpr ("truncated_") 
    {} 

  const STruncated* asSTruncated () const final
    { return this; }

  bool unify (const SExpr &/*pattern*/,
                Var2SExpr &/*var2SExpr*/) const final
    { return false; }
};




// Transformation

struct Transformation : Root
{
  const SExpr* lhs {nullptr};
  const SExpr* rhs {nullptr};
  // !nullptr


  explicit Transformation (TokenInput &ti);
  void qc () const override;


  bool empty () const final
    { return ! lhs && ! rhs; }
  void saveText (ostream &os) const override
    { lhs->saveText (os);
      os << " -> ";
      rhs->saveText (os);
      os << ';' << endl;
    }
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
    char -> delimiter | dia_letter
    delimiter -> content
    dia_letter -> letter capital  # bold underscore italics ...
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
  // Optional: map into the main letters ??
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
  SExpr* codepoint2SExpr (Codepoint c) const;
    // Return: !nullptr
    // Input: 0: text start/end
  SExpr* utf8_2SExpr (Utf8 &text) const;
    // Return: !nullptr

  void readTransformations (const string &tfmFName);
    // Output: trs
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
struct Codestring : Root
{
  const Alphabet& alphabet;
  Vector<Codepoint> arr;


  explicit Codestring (const Alphabet &alphabet_arg)
    : alphabet (alphabet_arg)
    {}
  void saveText (ostream &os) const final
    { for (const Codepoint c : arr)
        os << alphabet. small2name (c) << ' ';
    }
  bool empty () const final
    { return arr. empty (); }
  void clear () final
    { arr. clear (); }


  bool operator< (const Codestring &other) const
    { return arr < other. arr; }
  Codestring& operator<< (Codepoint c)
    { if (! c)
        throw logic_error ("0 code point");
      arr << c;
      return *this;
    }
  size_t size () const
    { return arr. size (); }
  Codepoint front () const
    { return arr. front (); }
  Codepoint back () const
    { return arr. back (); }
  bool contains (Codepoint c) const
    { return arr. contains (c); }
};



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



// -> common.hpp ??
template <typename T /*:Named*/>
  struct SExpr
  {
    const T* t {nullptr};
      // !nullptr
    Vector<SExpr<T>> children;
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

