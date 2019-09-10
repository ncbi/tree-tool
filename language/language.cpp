// language.cpp

#undef NDEBUG
#include "../common.inc"
//#include "cfgrammar.hpp"

#include "language.hpp"




namespace Lang_sp
{


// SExpr

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
        s = new SExprVariable (t1. name); 
      else
        s = new SExprGeneral (t1. name, 0); 
      break;
    case Token::eText: s = new SExprContent (t1. name); break;
    case Token::eDelimiter: 
      if (t1. isDelimiter (')'))
        s = new SExprEndOfList ();
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
    if (child->asSExprEndOfList ())
    {
      delete child;
      break;
    }
    ASSERT (child->asSExprNamed ());
    sg->children << child;
  }

  return sg;
}




// SExprNamed

void SExprNamed::qc () const 
{
 if (! qc_on)
   return;

  QC_ASSERT (! name. empty ());
  QC_ASSERT (! contains (name, ' '));
}




// SExprGeneral

void SExprGeneral::deleteChildren ()
{ 
  ASSERT (! referred);

  for (const SExpr* child : children)
    if (child)
    { 
      if (const SExprGeneral* sg = child->asSExprGeneral ())
        if (sg->referred)
        { 
          sg->referred = false;
          continue;
        }
      delete child;
    }
}



void SExprGeneral::qc () const 
{
 if (! qc_on)
   return;

  SExpr::qc ();

  for (const SExpr* child : children)
  {
    QC_ASSERT (child); 
    QC_ASSERT (! child->asSExprEndOfList ());
    child->qc ();
    QC_IMPLY (child->asSExprVariable (), children. size () == 1);
  }
}



void SExprGeneral::saveText (ostream &os) const
{ 
  os << name;
  if (children. empty ())
    return;
  if (isSmall ())
  {
    os << " (";
    bool first = true;
    for (const SExpr* child : children)
    { 
      if (! first)
        os << ' ';
      child->saveText (os);
      first = false;
    }
    os << ')';
  }
  else
  {
    const ebool left = isLeftList ();
    switch (left)
    {
      case EFALSE: os << " [*LIST]"; break;
      case ETRUE:  os << " [LIST*]"; break;
      default: break;
    }
    Offset::newLn (os);
    os << '(';
    {
      const Offset ofs;
      if (left == UBOOL)
        for (const SExpr* child : children)
        { 
          Offset::newLn (os);
          child->saveText (os);
        }
      else
        saveList (os, (bool) left);
    }
    Offset::newLn (os);
    os << ')';
  }
}



void SExprGeneral::saveList (ostream &os,
                             bool left) const
{ 
  const SExprGeneral* sg = this;
  for (;;)
  {
    if (sg->children. empty ())
      break;
    ASSERT (sg->children. size () == 2);
    Offset::newLn (os);
    sg->children [left] -> saveText (os);
    if (const SExprTruncated* st = sg->children [! left] -> asSExprTruncated ())
    {
      Offset::newLn (os);
      st->saveText (os);
      break;
    }
    sg = sg->children [! left] -> asSExprGeneral ();
    ASSERT (sg);
  }
}



SExprGeneral* SExprGeneral::copy () const
{
  auto sg = new SExprGeneral (name, children. size ());
  FFOR (size_t, i, children. size ())
    sg->children [i] = children [i] -> copy ();
  return sg;
}



Var2name SExprGeneral::getVarNames (const string &/*parentName*/) const
{ 
  Var2name var2name;
  for (const SExpr* child : children)
    if (const SExprNamed* sn = child->asSExprNamed ())
    {
      const Var2name other (sn->getVarNames (name));
      for (const auto& otherIt : other)
      {
        const auto& it = var2name. find (otherIt. first);
        if (it == var2name. end ())
          var2name [otherIt. first] = otherIt. second;
        else if (it->second != otherIt. second)
          throw logic_error ("Different SExpr names of variable " + strQuote (it->first));
      }
    }
    else
      throw logic_error ("getVarNames of a non-SExprNamed");

  return var2name;
}



const SExprGeneral* SExprGeneral::copySubstitute (const Var2SExpr &var2SExpr) const
{
  if (children. size () == 1)
    if (const SExprVariable* v = children [0] -> asSExprVariable ())
    {
      const auto& it = var2SExpr. find (v->name);
      if (it == var2SExpr. end ())
        throw logic_error (FUNC "Variable " + strQuote (v->name) + " is missing");  // Cannot happen ??
      const SExprGeneral* sg = it->second;
      ASSERT (sg);
      if (sg->referred)
        sg = sg->copy ();
      else
        sg->referred = true;
      ASSERT (name == sg->name);
      return sg;
    }

  auto sg = new SExprGeneral (name, children. size ());
  FFOR (size_t, i, children. size ())
    if (const SExprNamed* sn = children [i] -> asSExprNamed ())
      sg->children [i] = sn->copySubstitute (var2SExpr);
    else
      throw logic_error (FUNC "not SExprNamed");
  return sg;
}



bool SExprGeneral::unify (const SExprNamed &pattern,
                          Var2SExpr &var2SExpr) const 
{
  const SExprGeneral* p = pattern. asSExprGeneral ();
  if (! p)
    return false;

  if (name != pattern. name)
    return false;

  if (p->children. size () == 1)
    if (const SExprVariable* v = p->children [0] -> asSExprVariable ())
    {
      const auto& it = var2SExpr. find (v->name);
      if (it == var2SExpr. end ())
        var2SExpr [v->name] = this;
      else
      {
        ASSERT (it->second);
        return equal (* it->second);
      }
      return true;
    }
      
  if (children. size () != p->children. size ())
    return false;
  FFOR (size_t, i, children. size ())
  {
    const SExprNamed* p_child = p->children [i] -> asSExprNamed ();
    if (! p_child)
      throw logic_error (FUNC "Non-SExprNamed in pattern");
    if (! children [i] -> unify (*p_child, var2SExpr))
      return false;
  }

  return true;
}



bool SExprGeneral::apply (const Transformation &tr)
{
  Var2SExpr var2SExpr;
  if (! unify (* tr. lhs, var2SExpr))
    return false;

  SExprGeneral* sg = var_cast (tr. rhs->copySubstitute (var2SExpr));
  ASSERT (sg->name == name);
  deleteChildren ();
  children = move (sg->children);
  ASSERT (sg->children. empty ());
  delete sg;

  return true;
}



void SExprGeneral::applyDown (const Transformation &tr)
{
  while (apply (tr));

  for (const SExpr* child : children)
    if (const SExprGeneral* sg = child->asSExprGeneral ())
      var_cast (sg) -> applyDown (tr);
}




// Transformation

Transformation* Transformation::parse (TokenInput &ti)
{
  SExprGeneral* lhs = SExprGeneral::parse (ti);
  if (! lhs)
    return nullptr;
  Token (ti. ci, '-');
  Token (ti. ci, '>');
  SExprGeneral* rhs = SExprGeneral::parse (ti);
  QC_ASSERT (rhs);
  Token (ti. ci, ';');
  return new Transformation (lhs, rhs);
}



void Transformation::qc () const
{
  if (! qc_on)
    return;

  QC_ASSERT (lhs);
  QC_ASSERT (rhs);
  QC_ASSERT (lhs != rhs);
  lhs->qc ();
  rhs->qc ();
  QC_ASSERT (lhs->name == rhs->name);

  const Var2name lhsVars (lhs->getVarNames (string ()));
  const Var2name rhsVars (rhs->getVarNames (string ()));
#undef TRANSFORMATION_REVERSIBLE  // ??
#ifdef TRANSFORMATION_REVERSIBLE
  if (lhsVars. size () != rhsVars. size ())
    throw logic_error (FUNC "Different number of variables in LHS and RHS");
#endif
  for (const auto& it : lhsVars)
  {
    const auto& other = rhsVars. find (it. first);
    if (other == rhsVars. end ())
    #ifdef TRANSFORMATION_REVERSIBLE
      throw logic_error (FUNC "LHS variable " + strQuote (it. first) + " is missing in RHS");
    #else
      ;
    #endif
    else
      if (it. second != other->second)
        throw logic_error (FUNC "Variable " + strQuote (it. first) + " is under " + strQuote (it. second) + " and " + strQuote (other->second));
  }
  for (const auto& it : rhsVars)
  {
    const auto& other = lhsVars. find (it. first);
    if (other == lhsVars. end ())
    #ifdef TRANSFORMATION_REVERSIBLE
      throw logic_error (FUNC "RHS variable " + strQuote (it. first) + " is missing in LHS");
    #else
      ;
    #endif
    else 
      if (it. second != other->second)
        throw logic_error (FUNC "Variable " + strQuote (it. first) + " is under " + strQuote (other->second) + " and " + strQuote (it. second));
  }
}




// Utf8

Utf8::Utf8 (const string &fName)
: is (fName, ios_base::binary | ios_base::in)
{ 
  static_assert (sizeof (Codepoint) >= 4, "Too small Codepoint");
  ASSERT (is. good ());
  char c;
  EXEC_ASSERT (getChar (is, c));
  ASSERT (c == '\xEF');  
  EXEC_ASSERT (getChar (is, c));
  ASSERT (c == '\xBB');  
  EXEC_ASSERT (getChar (is, c));
  ASSERT (c == '\xBF');  
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
    c = (char) (2 * c);
    n++;
  }
  ASSERT (n > 1);
  FOR (size_t, i, n)
    c = (char) (c / 2);
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
  QC_ASSERT (startSmall >= startCapital + size);

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


  for (const Transformation* tr : trs)
  {
    QC_ASSERT (tr);
    tr->qc ();
  }
}



Codepoint Language::toSmall (Codepoint c) const
{ 
  ASSERT (c > 0);

  if (   c >= startSmall
      && c < startSmall + size
     )
    return c;
  if (   c >= startCapital
      && c < startCapital + size
     )
    return c + (Codepoint) capitalShift ();

  if (extraSmall. contains (c))
    return c;
  if (extraCapital. contains (c))
    return c + (Codepoint) extraCapitalShift ();

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



SExprGeneral* Language::codepoint2SExpr (Codepoint c) const
{ 
  SExpr* character = nullptr;
  if (c)
  {
    const size_t i = delimiters. indexOf (c);
    if (i != NO_INDEX)
    { 
      auto delimContent = new SExprContent (delimiterNames [i]);
      auto delim = new SExprGeneral ("delimiter", 1);
      delim->children [0] = delimContent;
      character = delim;
    }
    else
      if (const Codepoint cs = toSmall (c))
      { 
        auto letterContent = new SExprContent (small2name (cs));
        auto letter = new SExprGeneral ("letter", 1);
        letter->children [0] = letterContent;
        auto capital = new SExprGeneral ("capital", 1);
        if (cs == c)
          capital->children [0] = new SFalse ();
        else
          capital->children [0] = new STrue ();
        auto diaLetter = new SExprGeneral ("cap_letter", 2);
        diaLetter->children [0] = letter;
        diaLetter->children [1] = capital;  // diacritic
        character = diaLetter;
      }
      else
      { 
        auto letter = new SExprContent (to_string (c));
        auto foreign = new SExprGeneral ("foreign", 1);
        foreign->children [0] = letter;
        character = foreign;
      }
  }
  else
  {
    auto delimContent = new SExprContent ("#");
    auto delim = new SExprGeneral ("delimiter", 1);
    delim->children [0] = delimContent;
    character = delim;
  }
  ASSERT (character);

  auto s = new SExprGeneral ("char", 1);
  s->children [0] = character;

  return s;
}



SExprGeneral* Language::utf8_2SExpr (Utf8 &text,
                                     size_t char_max) const
{
  ASSERT (char_max);

  auto root = new SExprGeneral ("char_list", 2);
  root->children [0] = codepoint2SExpr (0);
  Codepoint c = 0;
  size_t i = 0; 
  SExprGeneral* parent = root;
  bool truncated = false;
  while (text. get (c))
  {
    ASSERT (c);
    auto char_list = new SExprGeneral ("char_list", 2);
    char_list->children [0] = codepoint2SExpr (c);
    parent->children [1] = char_list;
    parent = char_list;
    i++;
    if (i >= char_max)  
    {
      truncated = true;
      break;
    }
  }
  if (truncated)
    parent->children [1] = new SExprTruncated ();
  else
  {
    auto char_list = new SExprGeneral ("char_list", 2);
    char_list->children [0] = codepoint2SExpr (0);
    char_list->children [1] = new SExprGeneral ("char_list", 0);
    parent->children [1] = char_list;
  }

  return root;
}



void Language::readTransformations (const string &tfmFName)
{
  ASSERT (trs. empty ());

  TokenInput ti (tfmFName, '#');
  for (;;)
  {
    Transformation* tr = Transformation::parse (ti);
    if (! tr)
      break;
    try { tr->qc (); }
      catch (const exception &e)
      {
        throw runtime_error ("Line " + to_string (ti. ci. lineNum + 1) + ", position " + to_string (ti. ci. charNum + 1) + ": " + e. what ());
      }
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

