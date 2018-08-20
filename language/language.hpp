// language.hpp
// Natural language

#ifndef LANGUAGE_HPP
#define LANGUAGE_HPP

#include "../common.hpp"
using namespace Common_sp;



namespace Lang
{


// --> utf8.{hpp,cpp} ??
typedef  uint  Codepoint;
  // Unicode
  // 0 - non-existing code point



struct Utf8 
{
private:
  ifstream is;
public:

  explicit Utf8 (const string &fName);

  bool get (Codepoint &v);
    // Return: false <=> EOF
private:
  static size_t nextBytes (char &c);
    // Update: c - leading byte
};



struct Alphabet
{
  static const Vector<Codepoint> delimiters;
    // Of words
  // Main
  Codepoint startCapital;
  Codepoint startSmall;
  size_t size;
  // Dispensable letters: map into the main code points
  Vector<Codepoint> extraCapital;
  Vector<Codepoint> extraSmall;
  Vector<string> names;
  //
  Vector<Codepoint> vowels;
    // Small


  Alphabet (Codepoint startCapital_arg,
            Codepoint startSmall_arg,
            size_t size_arg)
    : startCapital (startCapital_arg)
    , startSmall (startSmall_arg)
    , size (size_arg)
    {}
  void qc () const;


private:
  size_t capitalShift () const
    { return startSmall - startCapital; }
  size_t extraCapitalShift () const
    { return extraCapital. empty () ? 0 : (extraSmall [0] - extraCapital [0]); }
  Codepoint toSmall_ (Codepoint c) const;
public:
  Codepoint toSmall (Codepoint c) const;
    // Return: 0 <=> not in Alphabet
    //         c is capital => Return != c
  string small2name (Codepoint c) const;
};



struct Russian : Alphabet  // ??
{
  Russian ()
    : Alphabet (0x410, 0x430, 32)
    { extraCapital << 0x401;  // Yo
      extraSmall   << 0x451;  // yo
      vowels << 0x430 << 0x435 << 0x438 << 0x43E << 0x443 << 0x44B << 0x44D << 0x44E << 0x44F << 0x451;
      names << "a" << "b" << "v" << "g" << "d" << "je" << "zh" << "z" << "ji" << "iy" << "k" << "l" << "m" << "n" << "o" << "p"
            << "r" << "s" << "t" << "u" << "f" << "x" << "c" << "ch" << "sh" << "sch" << "'" << "i" << "j" << "e" << "ju" << "ja"
            << "yo";
      qc ();
    }
};



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


}  // namespace



#endif

