// language.cpp

#undef NDEBUG
#include "../common.inc"
#include "cfgrammar.hpp"

#include "language.hpp"



namespace Lang
{


// Utf8

Utf8::Utf8 (const string &fName)
: is (fName, ios_base::binary | ios_base::in)
{ 
  static_assert (sizeof (Codepoint) >= 4, "Too small Codepoint");
  ASSERT (is. good ());
  char c;
  EXEC_ASSERT (getChar (is, c));
  ASSERT (c == '\xEF');  // ??
  EXEC_ASSERT (getChar (is, c));
  ASSERT (c == '\xBB');  // ??
  EXEC_ASSERT (getChar (is, c));
  ASSERT (c == '\xBF');  // ??
}



bool Utf8::get (Codepoint &v)
{
  char c;
  if (! getChar (is, c))
    return false;
  const size_t next = nextBytes (c);
  ASSERT (next < sizeof (Codepoint));
  v = static_cast <unsigned char> (c);
  // Continuation bytes
  FOR (size_t, i, next)
  {
    EXEC_ASSERT (getChar (is, c));
    ASSERT ((c & '\xC0') == '\x80');
    v *= 64;
    v += static_cast <unsigned char> (c & '\x3F');
  }
  return true;
}



size_t Utf8::nextBytes (char &c)
{
  if (! (c & '\x80'))
    return 0;
  size_t n = 0;
  while (c & '\x80')
  {
    c = c << 1;
    n++;
  }
  ASSERT (n > 1);
  c = c >> n;
  return n - 1;
}




// Alphabet

const Vector<Codepoint> Alphabet::delimiters 
  ({ 0x0D , 0x0A, ' ', 0x09, '.', ',', ';', ':', '!', '?', '&', '(', ')' 
  , '-', '+', '*', '/', '=', '[', ']', '{', '}', '|', '\\', '<', '>'
  , '"', '\'', 0xAB/*<<*/, 0xBB/*>>*/, 0xA9/*Copyright*/, 0x2014/*dash*/
  , '#', '$', '%', '@'
  });
                 


void Alphabet::qc () const
{ 
  if (! qc_on)
    return;

  QC_ASSERT (size);
  QC_ASSERT (startCapital > 32);
  QC_ASSERT (startSmall >= startCapital + (int) size);
  QC_ASSERT (extraCapital. size () == extraSmall.size ());
  if (! extraCapital. empty ())
  { 
    QC_ASSERT (extraSmall [0] > extraCapital [0]);
    FOR (size_t, i, extraCapital. size ())
    {
      QC_ASSERT (extraCapital [i]);
      QC_ASSERT (! toSmall_ (extraCapital [i]));
      QC_ASSERT (! toSmall_ (extraSmall   [i]));
      QC_ASSERT (extraSmall [i] == extraCapital [i] + extraCapitalShift ());
    }
  }
  for (const Codepoint v : vowels)
  {
    QC_ASSERT (v);
    QC_ASSERT (toSmall (v) == v);
  }
  // delimiters
  Set<Codepoint> delimiterSet;
  for (const Codepoint d : delimiters)
  {
    QC_ASSERT (d);
    QC_ASSERT (! toSmall (d));
    delimiterSet << d;
  }
  QC_ASSERT (delimiters. size () == delimiterSet. size ());
  // names
  QC_ASSERT (names. size () == size + extraCapital. size ());
  Set<string> nameSet;
  for (const string& name : names)
  {
    QC_ASSERT (! name. empty ());
    nameSet << name;
  }
  QC_ASSERT (nameSet. size () == names. size());
}



Codepoint Alphabet::toSmall_ (Codepoint c) const
{ 
  ASSERT (c >= 0);
  if (   c >= startSmall
      && c < startSmall + (int) size
      )
    return c;
  if (   c >= startCapital
      && c < startCapital + (int) size
      )
    return c + capitalShift ();
  return 0;
}



Codepoint Alphabet::toSmall (Codepoint c) const
{ 
  if (const Codepoint sc = toSmall_ (c))
    return sc;
  if (extraSmall. contains (c))
    return c;
  if (extraCapital. contains (c))
    return c + extraCapitalShift ();
  return 0;
}



string Alphabet::small2name (Codepoint c) const
{ 
  if (c < 128)  // ASCII
    return string (1, c);
  if (c >= startSmall && c < startSmall + size)
    return names [c - startSmall];
  size_t index;
  EXEC_ASSERT (extraSmall. find (c, index));
  return names [size + index];
}




// Language

namespace
{

Cfgr::Grammar* cfgr = nullptr;

void addCharRule (string &s, 
                  char c)
  { s += string ("char -> \"") + c + "\"\n"; }

}



Language::Language (const string &fName)
{
  if (! cfgr)
  {
    string s (R"cfgr(
sigma -> categ_rule
sigma -> cons_rule

categ_rule -> w sp "->" categ_choices
categ_choices -> categ_seq
categ_choices -> categ_seq sp "|" categ_choices
categ_seq -> w
categ_seq -> w space categ_seq

cons_rule -> cons sp "->" cons_choices
cons_choices -> cons_seq
cons_choices -> cons_seq sp "|" cons_choices
cons_seq -> cons
cons_seq -> cons space cons_seq

cons -> w sp "[" args sp "]" [pointer]
args -> arg
args -> arg sp "," args

arg -> SExpr [space w]
arg -> pointer

SExpr -> w
SExpr -> sp "(" SExprList sp ")"
SExprList -> SExpr
SExprList -> SExpr space SExprList

pointer -> sp "^" w

w -> sp char+
space -> " "
sp -> space*
)cfgr");
    for (char c = 'a'; c <= 'z'; c++)
      addCharRule (s, c);
    for (char c = 'A'; c <= 'Z'; c++)
      addCharRule (s, c);
    for (char c = '0'; c <= '9'; c++)
      addCharRule (s, c);
    addCharRule (s, '_');
    if (verbose ())
      cout << s << endl;  

    istringstream iss (s);
    CharInput in (iss);
    cfgr = new Cfgr::Grammar (in);
    cfgr->qc();
    if (verbose ())
      cfgr->print (cout);  
    cfgr->prepare ();  
  }


  LineInput f (fName, 100 * 1024, 1);  // PAR
  f. commentStart = "#";
  FOR (uchar, pass, 2)
    parseGrammar (f, pass);
}



void Language::parseGrammar (LineInput &f,
                             uchar pass)
{
  f. reset ();
  string line;
  while (f. nextLine ())
  {
    if (f. line. empty ())
      continue;
    replace (f. line, '\t', ' ');
    trim (f. line);
    if (f. line == "END")
      break;
    if (contains (f. line, arrowS))
    {
      parseRule (line, pass);
      line = f. line;
    }
    else
      line += " " + f. line;
  }
  parseRule (line, pass);
}



void Language::parseRule (const string &line,
                          uchar pass)
{
  ASSERT (cfgr);

  if (line. empty ())
    return;
  ASSERT (! isspace (line [0]));

  const size_t pos = line. find (arrowS);
  ASSERT (pos != string::npos);
  {
    const size_t rpos = line. rfind (arrowS);
    ASSERT (pos <= rpos);
    if (pos < rpos)
      throw runtime_error ("There can be only one \"->\" on a line");
  }

  string name (line. substr (0, pos));
  trim (name);

  if (pass ==  0)
  {
    if (contains (line, "["))
    {
      if (! findPtr (categories, name))
        categories [name] = new Category (name);
    }
    else
    {
      if (! findPtr (constructions, name))
        constructions [name] = new Construction (name);
    }
    return;
  }


  ASSERT (pass == 1);

  Cfgr::Sentence text (*cfgr, line);
  const Cfgr::Syntagms& syntagms = cfgr->parseSentence (text);

  const Cfgr::Syntagm* syntagm = nullptr;
  for (const Cfgr::Syntagm* syntagm_ : syntagms)
    if (syntagm_->right)
    {
      if (syntagm)
        throw runtime_error ("Ambiguous: " + line);
      syntagm = syntagm_;
    }
  if (! syntagm)
    throw runtime_error ("Not parsed: " + line);
  if (verbose ())
  {
    syntagm->print (cout);
    cout << endl;
  }
}



}  // namespace

