// Russian_test.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "language.hpp"
using namespace Lang_sp;



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
#endif



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Text processing")
    { 
      addPositional ("tfm",  "Language transformations"); 
      addPositional ("text", "Input txet file"); 
    }



  void body () const final
  {
    const string tfmFName  = getArg ("tfm");
    const string textFName = getArg ("text");
    ASSERT (! tfmFName. empty ());
    ASSERT (! textFName. empty ());


    RussianAlphabet russian;
    russian. readTransformations (tfmFName);
    russian. qc ();
    
    Utf8 textF (textFName);
    unique_ptr<SExprGeneral> root (russian. utf8_2SExpr (textF, 100));  // PAR
    ASSERT (root. get ());
    root->qc ();
  //root->saveText (cout);
  //cout << endl;

    russian. transform (*root);
    root->saveText (cout);
    cout << endl;


#if 0
    map <Codepoint, size_t> char2count;
    map <Codestring, size_t> consSeqs;  // small
    {
      Utf8 text (textFName);
      Codepoint v;
      Codestring word (russian);
      bool badWord = false;
      while (text. get (v))
      {
        char2count [v] ++;
        const bool delimiter = RussianAlphabet::delimiters. contains (v);
        const Codepoint sc = russian. toSmall (v);
        if (delimiter || russian. vowels. contains (sc))
        {
          if (delimiter)
            word << '#';
          if (! word. empty () && ! badWord)
            consSeqs [word] ++;
          word. clear ();
          badWord = false;
          if (delimiter)
            word << '#';
        }
        else
          if (sc)
            word << sc;
          else
            badWord = true;
      }
    }


    // Report
    Set<Codepoint> consonants;
    for (const auto& it : consSeqs)
      for (const Codepoint c : it. first. arr)
        consonants << c;
    for (const Codepoint c : consonants)
    {
      cout << endl << russian. small2name (c) << ":" << endl;
      //
      for (const auto& it : consSeqs)
      {
        const Codestring& s = it. first;
        if (   s. contains (c)
            && s. front () == '#'
            && s. size () > 2  // non-trivial
           )
          cout << s << ' ' << it. second << endl;
      }
      cout << endl;
      //
      for (const auto& it : consSeqs)
      {
        const Codestring& s = it. first;
        if (   s. contains (c)
            && s. back () == '#'
            && s. size () > 2  // non-trivial
           )
          cout << s << ' ' << it. second << endl;
      }
      cout << endl;
      //
      for (const auto& it : consSeqs)
      {
        const Codestring& s = it. first;
        if (   s. contains (c)
            && s. front () != '#'
            && s. back ()  != '#'
            && s. size () > 2  // non-trivial
           )
          cout << s << ' ' << it. second << endl;
      }
      cout << endl;
    }

  /*
    for (const auto it : char2count)
    {
      const Codepoint v = it. first;
      cout << hex << v << " " << (v >= ' ' && v < 255 ? (char) v : ' ') << " " << dec << it. second << endl;
    }
    cout << endl;
  */

  /*
    cout << endl;
    for (const auto it : consSeqs)
    {
      const Codestring& s = it. first;
      if (s. front () == '#' && s. back () == '#' && s. size () > 3)
        cout << s << ' ' << it. second << endl;
    }

    cout << endl;
    for (const auto it : consSeqs)
    {
      const Codestring& s = it. first;
      if (s. front () == '#' && s. back () != '#' && s. size () > 2)
        cout << s << ' ' << it. second << endl;
    }

    cout << endl;
    for (const auto it : consSeqs)
    {
      const Codestring& s = it. first;
      if (s. front () != '#' && s. back () == '#' && s. size () > 2)
        cout << s << ' ' << it. second << endl;
    }

    cout << endl;
    for (const auto it : consSeqs)
    {
      const Codestring& s = it. first;
      if (s. front () != '#' && s. back () != '#' && s. size () > 2)
        cout << s << ' ' << it. second << endl;
    }
  */
  #endif
  }
};



int main (int argc,
          const char* argv [])
{
  ThisApplication app;
  return app. run (argc, argv);
}
