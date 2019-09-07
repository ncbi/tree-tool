// language.cpp

#undef NDEBUG
#include "../common.inc"
#include "cfgrammar.hpp"

#include "language.hpp"



namespace Lang_sp
{


// SExpr

void SExpr::qc () const 
{
 if (! qc_on)
   return;

  QC_ASSERT (! name. empty ());
  QC_ASSERT (! contains (name, ' '));
}



SExpr* SExpr::parse (TokenInput &ti)
{
  const Token t1 (ti. get ());
  if (t1. empty ())
    return nullptr;

  SExpr* s = nullptr;
  switch (t1. type)
  {
    case Token::eName: 
      if (isUpper (t1. name [0]))
        s = new SVariable (t1. name); 
      else
        s = new SExprGeneral (t1. name, 0); 
      break;
    case Token::eText: s = new SContent (t1. name); break;
    case Token::eDelimiter: 
      if (t1. isDelimiter (')'))
        s = new SEndOfList ();
      else if (   t1. isDelimiter ('-')
               || t1. isDelimiter (';')
              )
        throw runtime_error (FUNC "')' is expected");
      else
        throw runtime_error (FUNC "S-expression cannot start with a delimiter");
      break;
    default: throw runtime_error (FUNC "bad Token::type");
  }
  ASSERT (s);
  if (! s->asSExprGeneral ())
    return s;

  const Token t2 (ti. get ());
  if (! t2. isDelimiter ('('))
  {
    ti. last = t2;
    return s;
  }

  SExprGeneral* sg = var_cast (s->asSExprGeneral ());
  ASSERT (sg);
  ASSERT (sg->children. empty ());
  for (;;)
  {
    SExpr* child = SExpr::parse (ti);
    if (! child)
      throw runtime_error (FUNC "EOF");
    if (child->asSEndOfList ())
    {
      delete child;
      break;
    }
    sg->children << child;
  }

  return sg;
}




// SExprGeneral

void SExprGeneral::qc () const 
{
 if (! qc_on)
   return;

  SExpr::qc ();

  for (const SExpr* child : children)
  {
    QC_ASSERT (child); 
    QC_ASSERT (! child->asSEndOfList ());
    child->qc ();
    QC_IMPLY (child->asSVariable (), children. size () == 1);
  }
}



bool SExprGeneral::unify (const SExpr &pattern,
                          Var2SExpr &var2SExpr) const 
{
  if (const SExprGeneral* p = pattern. asSExprGeneral ())
  {
    if (children. size () != p->children. size ())
      return false;
    FFOR (size_t, i, children. size ())
      if (! children [i] -> unify (* p->children [i], var2SExpr))
        return false;
    return true;
  }

  if (const SVariable* v = pattern. asSVariable ())
    return v->assign (*this, var2SExpr);

  return false;
}




// SContent

bool SContent::unify (const SExpr &pattern,
                      Var2SExpr &var2SExpr) const
{ 
  if (const SVariable* v = pattern. asSVariable ())
    return v->assign (*this, var2SExpr);
  return    pattern. asSContent () 
         && name == pattern. name; 
}




// SVariable

bool SVariable::assign (const SExpr &target,
                        Var2SExpr &var2SExpr) const
{
  auto& it = var2SExpr. find (name);
  if (it == var2SExpr. end ())
    throw logic_error ("Unknown variable " + strQuote (name) + " in pattern");
  if (it->second)
    return target. equal (* it->second);
  it->second = & target;

  return true;
}




// Transformation

Transformation::Transformation (TokenInput &ti)
{
  try 
  {
    lhs = SExpr::parse (ti);
    if (! lhs)
      return;
    Token (ti. ci, '-');
    Token (ti. ci, '>');
    rhs = SExpr::parse (ti);
    QC_ASSERT (rhs);
    Token (ti. ci, ';');
  }
  catch (const exception &e)
  {  
    throw runtime_error ("Line " + to_string (ti. ci. lineNum + 1) + ", position " + to_string (ti. ci. charNum + 1) + ": " + e. what ());
  }
}



void Transformation::qc () const
{
  if (! qc_on)
    return;

  if (empty ())
    return;

  QC_ASSERT (lhs);
  QC_ASSERT (rhs);
  lhs->qc ();
  rhs->qc ();

  const StringVector lhsVars (lhs->getVarNames ());
  const StringVector rhsVars (rhs->getVarNames ());
  QC_ASSERT (lhsVars == rhsVars);
}



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




// Language

Language::Language (Codepoint startCapital_arg,
                    Codepoint startSmall_arg,
                    size_t size_arg)
: startCapital (startCapital_arg)
, startSmall (startSmall_arg)
, size (size_arg)
{
  // delimiters, delimiterNames
  // spaces
  delimiters     << 0x0D << 0x0A << ' '     << 0x09  << 0x0C; 
  delimiterNames << "cr" << "lf" << "space" << "tab" << "ff";
  // punctuation
  delimiters     << '.'      << ','     << ';'         << ':'     << '!'           << '?'        << '-'      << 0x2014;
  delimiterNames << "period" << "comma" << "semicolon" << "colon" << "exclamation" << "question" << "hyphen" << "dash";
  // quotes
  delimiters     << '"'       << '\''      << 0xAB/*<<*/   << 0xBB/*>>*/;
  delimiterNames << "d_quote" << "s_quote" << "open_quote" << "close_quote";
  // parentheses
  delimiters     << '('           << ')'            << '['            << ']'             << '{'          << '}';
  delimiterNames << "open_parens" << "close_parens" << "open_bracket" << "close_bracket" << "open_brace" << "close_brace";
  // slashes
  delimiters     << '/'     << '\\';
  delimiterNames << "slash" << "back_slash";
  // signs
  delimiters     << '&'        << '+'    << '*'    << '='     << '<'    << '>'       << '#'     << '$'      << '%'       << '@'  << 0xA9;
  delimiterNames << "ampesand" << "plus" << "star" << "equal" << "less" << "greater" << "pound" << "dollar" << "percent" << "at" << "copyright";

  // digits, digitNames
  digits     << '0'    << '1'   << '2'   << '3'     << '4'    << '5'    << '6'   << '7'     << '8'     << '9';
  digitNames << "zero" << "one" << "two" << "three" << "four" << "five" << "six" << "seven" << "eight" << "nine";
}
                 


void Language::qc () const
{ 
  if (! qc_on)
    return;

  // delimiters
  Set<Codepoint> delimiterSet;
  for (const Codepoint d : delimiters)
  {
    QC_ASSERT (d);
    delimiterSet << d;
  }
  QC_ASSERT (delimiters. size () == delimiterSet. size ());
  QC_ASSERT (delimiters. size () == delimiterNames. size ());

  // digits
  Set<Codepoint> digitSet;
  for (const Codepoint d : digits)
  {
    QC_ASSERT (d);
    digitSet << d;
  }
  QC_ASSERT (digits. size () == digitSet. size ());
  QC_ASSERT (digits. size () == digitNames. size ());

  QC_ASSERT (size);
  QC_ASSERT (startCapital > 32);
  QC_ASSERT (startSmall >= startCapital + (int) size);

  QC_ASSERT (extraCapital. size () == extraSmall.size ());
  if (! extraCapital. empty ())
  { 
    QC_ASSERT (extraSmall [0] > extraCapital [0]);
    FFOR (size_t, i, extraCapital. size ())
    {
      QC_ASSERT (extraCapital [i]);
      QC_ASSERT (toSmall (extraSmall [i]) == extraSmall [i]);
      QC_ASSERT (extraSmall [i] == extraCapital [i] + extraCapitalShift ());
    }
  }

#if 0
  for (const Codepoint v : vowels)
  {
    QC_ASSERT (v);
    QC_ASSERT (toSmall (v) == v);
  }
#endif

  // letterNames
  QC_ASSERT (letterNames. size () == size + extraCapital. size ());
  Set<string> nameSet;
  for (const string& name : letterNames)
  {
    QC_ASSERT (! name. empty ());
    nameSet << name;
  }
  QC_ASSERT (nameSet. size () == letterNames. size ());


  for (const Transformation& tr : trs)
    tr. qc ();
}



Codepoint Language::toSmall (Codepoint c) const
{ 
  ASSERT (c > 0);

  if (   c >= startSmall
      && c < startSmall + (int) size
     )
    return c;
  if (   c >= startCapital
      && c < startCapital + (int) size
     )
    return c + capitalShift ();

  if (extraSmall. contains (c))
    return c;
  if (extraCapital. contains (c))
    return c + extraCapitalShift ();

  return 0;
}



const string& Language::small2name (Codepoint c) const
{ 
  if (c >= startSmall && c < startSmall + size)
    return letterNames [c - startSmall];
  size_t index;
  EXEC_ASSERT (extraSmall. find (c, index));
  return letterNames [size + index];
}



SExpr* Language::codepoint2SExpr (Codepoint c) const
{ 
  SExpr* character = nullptr;
  if (c)
  {
    const size_t i = delimiters. indexOf (c);
    if (i != NO_INDEX)
    { 
      auto delimContent = new SContent (delimiterNames [i]);
      auto delim = new SExprGeneral ("delimiter", 1);
      delim->children [0] = delimContent;
      character = delim;
    }
    else
      if (const Codepoint cs = toSmall (c))
      { 
        auto letterContent = new SContent (small2name (cs));
        auto letter = new SExprGeneral ("letter", 1);
        letter->children [0] = letterContent;
        auto capital = new SExprGeneral ("capital", 1);
        if (cs == c)
          capital->children [0] = new SNil ();
        else
          capital->children [0] = new STrue ();
        auto diaLetter = new SExprGeneral ("dia_letter", 2);
        diaLetter->children [0] = letter;
        diaLetter->children [1] = capital;  // diacritic
        character = diaLetter;
      }
      else
      { 
        auto letter = new SContent (to_string (c));
        auto foreign = new SExprGeneral ("foreign", 1);
        foreign->children [0] = letter;
        character = foreign;
      }
  }
  else
  {
    auto delimContent = new SContent ("#");
    auto delim = new SExprGeneral ("delimiter", 1);
    delim->children [0] = delimContent;
    character = delim;
  }
  ASSERT (character);

  auto s = new SExprGeneral ("char", 1);
  s->children [0] = character;

  return s;
}



SExpr* Language::utf8_2SExpr (Utf8 &text) const
{
  auto root = new SExprGeneral ("char_star", 2);
  root->children [0] = codepoint2SExpr (0);
  Codepoint c = 0;
  size_t i = 0; // ??
  SExprGeneral* parent = root;
  bool truncated = false;
  while (text. get (c))
  {
    ASSERT (c);
    auto char_star = new SExprGeneral ("char_star", 2);
    char_star->children [0] = codepoint2SExpr (c);
    parent->children [1] = char_star;
    parent = char_star;
    i++;
    if (i >= 100)  // PAR
    {
      truncated = true;
      break;
    }
  }
  if (truncated)
    parent->children [1] = new STruncated ();
  else
  {
    auto char_star = new SExprGeneral ("char_star", 2);
    char_star->children [0] = codepoint2SExpr (0);
    char_star->children [1] = new SNil ();
    parent->children [1] = char_star;
  }

  return root;
}



void Language::readTransformations (const string &tfmFName)
{
  ASSERT (trs. empty ());

  TokenInput ti (tfmFName, '#');
  for (;;)
  {
    const Transformation tr (ti);
    tr. qc ();
    if (tr. empty ())
      break;
    trs << tr;
  }
}



#if 0
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
#endif



}  // namespace

