// grammar.cpp

#undef NDEBUG
#include "common.inc"

#include "grammar.hpp"



namespace Grammar_sp
{


bool initGrammar ()
{
  MODULE_INIT

  static_assert (numeric_limits<Char>::is_integer, "Char must be integer");

  return true;
}


namespace
{
  const bool init_ = initGrammar ();
}




string chars2str (const Vector<Char> &vec)
{
  string s;  s. reserve (vec. size ());
  for (const Char c : vec)
  {
    ASSERT (c < 128);
    s += (char) c;
  }
  
  return move (s);
}




// Sentence

Sentence::Sentence (const Grammar &grammar,
                    const string &s) throw (BadPos)
{
  seq. resize (s. size () + 2);
  seq [0]. init (grammar, eot). right = true;
  FOR (size_t, i, s. size ())
    try 
    { 
      const char c = s [i];
      ASSERT (c != eot);
      seq [i + 1]. init (grammar, (Char) c); 
    }
    catch (const exception &e) 
      { throw BadPos (i, e. what ()); }
  seq [s. size () + 1]. init (grammar, eot). right = true;
}



void Sentence::qc () const
{
  ASSERT (! seq. empty ());
#ifndef NDEBUG
  FOR (size_t, i, seq. size ())
  {
    const Position& pos = seq [i];
    if (pos. c == eot)
      { ASSERT (i == 0 || i == seq. size () - 1); }
    else
      pos. qc ();
  }
#endif  
}



size_t Sentence::countWrongTerminalSyntagms () const
{
  size_t n = 0;
  for (const Position& pos : seq)
    for (const auto mapIt : pos. terminal2syntagms)
      for (const Syntagm* syntagm : * mapIt. second)
        if (! syntagm->right)  
          n++;
  return n;
}



size_t Sentence::countWrongNonTerminalSyntagms () const
{
  size_t n = 0;
  for (const Position& pos : seq)
    for (const auto mapIt : pos. nonTerminal2syntagms)
      for (const Syntagm* syntagm : * mapIt. second)
        if (! syntagm->right)  
          n++;
  return n;
}



size_t Sentence::countRollbacks () const
{
  size_t n = 0;
  for (const Position& pos : seq)
    n += pos. rollbacks;
  return n;
}



const NonTerminalSyntagm* Sentence::getLastSyntagm () const
{
  const Syntagm* lastSyntagm = nullptr;
  ptrdiff_t index = 0;
  for (const Position& pos : seq) 
    for (const auto it : pos. nonTerminal2syntagms)
      for (const Syntagm* syntagm : *(it. second))
      {
        const ptrdiff_t diff = & syntagm->end - & seq [0];
        ASSERT (diff > 0);
        if (maximize (index, diff))
          lastSyntagm = syntagm;
      }
        
  if (! lastSyntagm)
    return nullptr;  
  const NonTerminalSyntagm* res = lastSyntagm->asNonTerminalSyntagm ();
  return res;
}



void Sentence::reportError (ostream &os) const
{
  const NonTerminalSyntagm* lastSyntagm = getLastSyntagm ();
  if (! lastSyntagm)
    return;
  size_t index = 0;
  bool reached = false;
  for (const Position& pos : seq)
  {
    if (& pos == & lastSyntagm->end)
    {
      reached = true;
      os << '|';
    }
    os << (char) pos. c;
    if (! reached)
      index++;
  }
  os << endl;
  os << "  At: " << index << endl;
}



void Sentence::printWrongSyntagms (ostream &os,
                                   size_t size_min) const
{
  for (const Position& pos : seq) 
    for (const auto it : pos. nonTerminal2syntagms)
      for (const Syntagm* syntagm : *(it. second))
        if (! syntagm->right 
            && syntagm->size () >= size_min
           )
        {
          syntagm->Syntagm::saveText (os);
          os << endl;
        }
}



VectorPtr<NonTerminalSyntagm> Sentence::findSyntagms (const string& symbolName) const
{
  VectorPtr<NonTerminalSyntagm> vec;
  for (const Position& pos : seq)
    for (const auto it : pos. nonTerminal2syntagms)
      for (const Syntagm* syntagm : *(it. second))
        if (   syntagm->right 
            && syntagm->symbol. name == symbolName
           )
          vec << syntagm->asNonTerminalSyntagm ();
  return move (vec);
}




// Syntagm

void Syntagm::qc () const
{
  ASSERT (begin. c != eot);
}



void Syntagm::saveText (ostream &os) const
{ 
  os << symbol. name << ": " << chars2str (str ()); 
}



size_t Syntagm::size () const
{ 
  const long int diff = & end - & begin;
  ASSERT (diff >= 0);
  return (size_t) diff; 
}




// TerminalSyntagm

#if 0
TerminalSyntagm::TerminalSyntagm (const Position& begin_arg,
                                  const Position& end_arg,
                                  const Terminal& terminal)
: Syntagm (begin_arg, end_arg, terminal)
{
  init (begin_arg, terminal);
}
#endif



TerminalSyntagm::TerminalSyntagm (const Position& begin_arg,
                                  const Terminal& terminal)
: Syntagm (begin_arg, * (& begin_arg + 1), terminal)
{ 
  init (begin_arg, terminal);
}



void TerminalSyntagm::init (const Position& begin_arg, 
                            const Terminal& terminal)
{
  const_cast <Syntagms&> (findMake (const_cast <Position&> (begin_arg). terminal2syntagms, & terminal)) << this;
}



void TerminalSyntagm::qc () const
{
  Syntagm::qc ();

  ASSERT (symbol. asTerminal ());
  ASSERT (findPtr (begin. terminal2syntagms, symbol. asTerminal ()));
  ASSERT (& begin != & end);
}



Vector<Char> TerminalSyntagm::str () const 
{ 
  Vector<Char> vec (1, begin. c); 
  return move (vec);
}




// NonTerminalSyntagm

NonTerminalSyntagm::NonTerminalSyntagm (const Position& begin_arg,
                                        const Position& end_arg,
                                        const Rule &rule_arg,
                                        const VectorPtr<Syntagm> &children_arg)
: Syntagm (begin_arg, end_arg, rule_arg. lhs)
, rule (rule_arg)
, children (children_arg)  // move() ??
{
  const Syntagms* syntagms = findPtr (begin. nonTerminal2syntagms, & rule.lhs);
  ASSERT (syntagms);
  * const_cast <Syntagms*> (syntagms) << this;
}



void NonTerminalSyntagm::qc () const
{
  Syntagm::qc ();

  ASSERT (symbol. asNonTerminal ());
  ASSERT (findPtr (begin. nonTerminal2syntagms, symbol. asNonTerminal ()));
  ASSERT (& symbol == & rule. lhs);
  ASSERT (children. size () == rule. rhs. size ());
#ifndef NDEBUG
  const Syntagm* prev = nullptr;
  FOR (size_t, i, children. size ())
  {
    const Syntagm* syntagm = children [i];
    ASSERT (syntagm);
    ASSERT (& syntagm->symbol == rule. rhs [i]);
    IMPLY (i == 0, & syntagm->begin == & begin);
    IMPLY (i == children. size () - 1, & syntagm->end == & end);
    IMPLY (prev, & prev->end == & syntagm->begin);
  //syntagm->qc ();   // Invoked in Position::qc()
    prev = syntagm;
  }
#endif
}



void NonTerminalSyntagm::saveText (ostream &os) const
{
  rule. saveText (os);
  Offset ofs;
  for (const Syntagm* syntagm : children)
  {
    ofs. newLn (os);
    syntagm->saveText (os);
  }
}



Vector<Char> NonTerminalSyntagm::str () const
{
  Vector<Char> vec;  vec. reserve (size ());
  for (const Syntagm* syntagm : children)
    vec << syntagm->str ();
/*
  for (const Position* pos = & begin; pos != & end; pos++)
    vec << pos->c;    
*/
  return move (vec);
}




void NonTerminalSyntagm::setRight () 
{ 
  Syntagm::setRight ();
  for (const Syntagm* syntagm : children)
    const_cast <Syntagm*> (syntagm) -> setRight ();
}




// Position

void Position::qc () const
{
  ASSERT (! terminal2syntagms. empty ()); 
#ifndef NDEBUG
  const Syntagm* prev = nullptr;
  for (const auto mapIt : terminal2syntagms)
  {
    ASSERT (mapIt. first);
    ASSERT (mapIt. second);
    Set<const Syntagm*> dSet;
    for (const Syntagm* syntagm : * mapIt. second)
    {
      ASSERT (syntagm);
      syntagm->qc ();
      ASSERT (syntagm->asTerminalSyntagm ());
      ASSERT (& syntagm->begin == this);
      ASSERT (c == syntagm->symbol. asTerminal () -> c);
      ASSERT (mapIt. first == & syntagm->symbol);
      dSet << syntagm;
      IMPLY (prev, & prev->begin == & syntagm->begin);
      prev = syntagm;
    }
    ASSERT (dSet. size () == mapIt. second->size ());
  }
  for (const auto mapIt : nonTerminal2syntagms)
  {
    ASSERT (mapIt. first);
    ASSERT (mapIt. second);
    Set<const Syntagm*> dSet;
    for (const Syntagm* syntagm : * mapIt. second)
    {
      ASSERT (syntagm);
      syntagm->qc ();
      ASSERT (syntagm->asNonTerminalSyntagm ());
      ASSERT (& syntagm->begin == this);
      ASSERT (mapIt. first == & syntagm->symbol);
      dSet << syntagm;
      IMPLY (prev, & prev->begin == & syntagm->begin);
      prev = syntagm;
    }
    ASSERT (dSet. size () == mapIt. second->size ());
  }
#endif
}



TerminalSyntagm& Position::init (const Grammar &grammar,
                                 Char c_arg)
{
  c = c_arg;
/*if (isspace (c))
    c = ' '; */
  if (const Terminal* t = findPtr (grammar. char2terminal, c))
    return * new TerminalSyntagm (*this, *t); 
  else
    throw runtime_error ("Bad character");
}




// Rule::Occurrences

Set<const Symbol*> Rule::Occurrences::getSymbols () const
{
  Set<const Symbol*> res;
  for (const Rule::Occurrence ro : *this)
    res << ro. getSymbol ();
  return move (res);
}




// Rule::Occurrence

void Rule::Occurrence::qc () const
{
  ASSERT (*rhsIt);
  ASSERT (rhsIt - rule. rhs. begin () >= 0);
}



void Rule::Occurrence::saveText (ostream &os) const
{
  os << rule. num << '.' << getIndex ();
}



size_t Rule::Occurrence::getIndex () const
{ 
  const long int diff = rhsIt - rule. rhs. begin ();
  ASSERT (diff >= 0);
  return (size_t) diff; 
}



string Rule::Occurrence::getName () const
{
  ostringstream oss;
  oss << getSymbol () -> name << '/';
  saveText (oss);
  return oss. str ();
}



bool Rule::Occurrence::operator< (const Occurrence &other) const
{ 
  LESS_PART (rule, other. rule, num);
  return getIndex () < other. getIndex (); 
}



bool Rule::Occurrence::ruleCanStartWith () const
{
  CONST_ITER (Rhs, it, rule. rhs)
    if (it == rhsIt)
      return true;
    else if (! (*it)->erasable)
      return false;
  NEVER_CALL;
}



bool Rule::Occurrence::ruleCanEndWith () const
{
  RhsIt it = rule. rhs. end ();
  do
  {
    it--;
    if (it == rhsIt)
      return true;
    if (! (*it)->erasable)
      return false;
  }
  while (it != rule. rhs. begin ());
  NEVER_CALL;
}



Rule::Occurrences Rule::Occurrence::getFirstROs () const
{
  Occurrences res (getSymbol () -> firstROs);
  res << *this;
  return move (res);
}



Rule::Occurrences Rule::Occurrence::getLastROs () const
{
  Occurrences res (getSymbol () -> lastROs);
  res << *this;
  return move (res);
}




// Rule

void Rule::finish (Grammar &grammar)
{
  ASSERT (rhs. empty ());
  rhs. resize (rhsS. size ());
  FOR (size_t, i, rhsS. size ())
  {
    const string& name = rhsS [i];
    Symbol* s = grammar. string2Symbol [name];
    if (! s)
      s = grammar. getSymbol<Terminal /*Letter*/> (name);
    rhs [i] = s;
    s->ruleOccurrences << Rule::Occurrence (*this, rhs. begin () + (int) i);
  }
  ASSERT (rhs. size () == rhsS. size ());
  rhsS. clear ();
}



void Rule::qc () const
{
  Root::qc ();

  ASSERT (! rhs. contains (nullptr));
  IMPLY (empty (), isErasable ());
  IMPLY (isErasable (), isTerminable ());
  IMPLY (isErasable (), lhs. erasable);
  IMPLY (isErasable (), erasable);
  ASSERT (firstTerminals. empty () == erasable);
  ASSERT (! firstTerminals. contains (nullptr));
  ASSERT (singleRhs. size () == rhs. size ());
#if 0
#ifndef NDEBUG
  for (const Symbol* s : rhs)
    ASSERT (lhs. grammar == s->grammar);
#endif
#endif

  // Non-redundant
  IMPLY (rhs. size () == 1, & lhs != rhs [0]);
}



void Rule::saveText (ostream &os) const
{ 
  os << "Rule# " << num << ":   " << lhs. name << ' ' << Grammar::arrowS;
  for (const Symbol* s : rhs)  
    os << ' ' << s->name; 

#if 0
  os << "  1st:";
  for (const Terminal* t : firstTerminals)
    os << ' ' << t->name;
#endif
}



bool Rule::isLeftRecursive () const
{ 
  return    ! empty () 
         && & lhs == rhs [0]; 
}



bool Rule::isRightRecursive () const
{ 
  return    ! empty () 
         && & lhs == rhs [rhs. size () - 1]; 
}



bool Rule::isTerminable () const
{
  for (const Symbol* s : rhs)
    if (! s->terminable)
      return false;
  return true;
}



bool Rule::isErasable () const
{
  for (const Symbol* s : rhs)
    if (! s->erasable)
      return false;
  return true;
}



#if 0
bool Rule::isLeftRecursiveErasable () const
{
  ASSERT (isLeftRecursive ());

  bool first = true;
  for (const Symbol* s : rhs)
  {
    if (   ! first 
        && ! s->erasable
       )
      return false;
    first = false;
  }
  return true;
}
#endif



#if 0
Rule::Occurrences Rule::canStartWith (const Symbol* s) const
{
//ASSERT (rhs. contains (s));

  Occurrences ros;
  CONST_ITER (Rhs, it, rhs)
  {
    if (*it == s)
      ros << Occurrence (*this, it);
    if (! (*it)->erasable)
      break;
  }

  return move (ros);
}



Rule::Occurrences Rule::canEndWith (const Symbol* s) const
{
//ASSERT (rhs. contains (s));

  Occurrences ros;
  RhsIt it = rhs. end (); 
  do
  {
    it--;
    if (*it == s)
      ros << Occurrence (*this, it);
    if (! (*it)->erasable)
      break;
  }
  while (it != rhs. begin ());

  return move (ros);
}
#endif



void Rule::setFirstTerminalsErasable ()
{
  ASSERT (firstTerminals. empty ());
  ASSERT (! erasable);

  const bool leftRecursive = isLeftRecursive ();
  bool first = true;
  for (const Symbol* s : rhs)
  {
    if (leftRecursive && first)
      ;
    else
    {
      for (const Symbol* f : s->getFirstSymbols ())
        if (const Terminal* t = f->asTerminal ())
          firstTerminals << t;
      if (! s->erasable)
        return;
    }
    first = false;
  }

  firstTerminals. clear ();
  erasable = true;
}




void Rule::setSingleRhs ()
{
  ASSERT (! singleRhs. contains (true));
  FOR (size_t, i, rhs. size ())
  {
    bool single = true;
    FOR (size_t, j, rhs. size ())
    if (   j != i
        && ! rhs. at (j) -> erasable
       )
      single = false;
    singleRhs [i] = single;
  }
}
  


void Rule::parseIt (RhsIt rhsIt,
                    const Position &pos,
                    const VectorPtr<Syntagm> &children) const
{
  IMPLY (! children. empty (), & children. back () -> end == & pos);
  
  const Position& begin = children. empty () ? pos : children [0] -> begin;

#if 0
  if (verbose ())
  {
    cout << endl;
    cout << & pos << "  " << rhsIt - rhs. begin () << "  " << children. size () << "  " << & begin << "  ";
    print (cout);
    cout << endl;
    for (const Syntagm* syntagm : children)
      cout << syntagm->symbol. name << ": " << chars2str (syntagm->str ()) << '|' << endl;
  }
#endif
      
  
  if (rhsIt == rhs. end ())
  {
    if (   lhs. name != NonTerminal::sigmaS 
        && ! lhs. canEnd (pos)
       )
      return;
    
    NonTerminalSyntagm* nts = new NonTerminalSyntagm (begin, pos, *this, children);

    // Rule::isLeftRecursive()
    if (lhs. name != NonTerminal::sigmaS)  // Otherwise pos has no content
    {
      VectorPtr<Syntagm> children_new;  children_new. reserve (1); 
      children_new << nts;
      for (const auto it : pos. terminal2syntagms)
        if (const VectorPtr<Rule>* rules = findPtr (lhs. terminal2rules [true], it. first))
          for (const Rule* rule : *rules)
          {
            ASSERT (rule->isLeftRecursive ());
            rule->parseIt (std::next (rule->rhs. begin ()), pos, children_new); 
          }
    }
  }
  else if (const Syntagms* syntagms = (*rhsIt)->parse (pos))
  {
    rhsIt++;
    for (const Syntagm* syntagm : *syntagms)
    {
      VectorPtr<Syntagm> children_new;  children_new. reserve (children. size () + 1);
      children_new = children;
      children_new << syntagm;
      parseIt (rhsIt, syntagm->end, children_new);
    }
  }
}



#if 0
// Deterministic parser
void Rule::parse (Pos &pos,
                  VectorPtr<Syntagm> &children) const
{
  Pos::CheckPoint cp (pos);
  size_t pushed = 0;
  const bool leftRecursive = isLeftRecursive ();
  bool first = true;
  for (const Symbol* s : rhs)
  {
    if (leftRecursive && first)
      ;
    else if (const Syntagm* syntagm = s->parse (pos))
    {
      children << syntagm;
      pushed++;
    }
    else
    {
      cp. rollback ();  
      children. pop (pushed);
      return false;
    }
    first = false;
  }

  return true;
}
#endif




// Symbol

void Symbol::qc () const
{
  Named::qc ();
    
  ASSERT (! isLeft  (name, "["));
  ASSERT (! isRight (name, "]"));

  IMPLY (erasable, terminable);

//IMPLY (! grammar, asGrammar());
#ifndef NDEBUG
  for (const Rule::Occurrence ro : ruleOccurrences)
  {
    ro. qc ();
    ASSERT (ro. getSymbol () == this);
  }
#endif

  ASSERT (! terminals. contains (nullptr));
#ifndef NDEBUG
  for (const Rule::Occurrence ro : firstROs)
    ro. qc ();
  for (const Rule::Occurrence ro : lastROs)
    ro. qc ();
  bool found = false;
  for (const Symbol* s : getFirstSymbols ())
  {
    ASSERT (s);
  //ASSERT (s->grammar == grammar);
    if (const Terminal* t = s->asTerminal ())
    {
      found = true;
      ASSERT (terminals. contains (t));
    }
  }
  // Non-redundant
  ASSERT (/*(bool) grammar ==*/ found)

  found = false;
  for (const Symbol* s : getLastSymbols ())
  {
    ASSERT (s);
  //ASSERT (s->grammar == grammar);
    if (const Terminal* t = s->asTerminal ())
    {
      found = true;
      ASSERT (terminals. contains (t));
    }
  }
  // Non-redundant
  ASSERT (/*(bool) grammar ==*/ found)

  ASSERT (! terminals. contains (nullptr));
#endif

  // Non-redundant
#ifndef NDEBUG
  for (const bool out : Bool)
  {
    IMPLY (name != NonTerminal::sigmaS, ! neighbors [out]. empty ());
    ASSERT (! neighbors [out]. contains (nullptr));
  }
#endif


  ASSERT (! replacements. contains (nullptr));
  ASSERT (replacements. contains (this));
}



void Symbol::saveText (ostream &os) const
{
  os << name << ":" << endl;
  
  for (const bool out : Bool)
  {
    os << (out ? "Next" : "Prev") << ":";
    for (const Symbol* s : neighbors [out])
      os << ' ' << s->name;
    os << endl;
  }
  
  os << "Replacements:";
  for (const Symbol* s : replacements)
    os << ' ' << s->name;
  os << endl;  
}



bool Symbol::isRoot () const
{ 
  return    ruleOccurrences. size () == 1 
         && ruleOccurrences. front (). rule. lhs. name == NonTerminal::sigmaS;
}



void Symbol::setErasable ()
{
  if (erasable)
    return;
  erasable = true;
  for (const Rule::Occurrence ro : ruleOccurrences)
    if (ro. rule. isErasable ())
      const_cast <NonTerminal&> (ro. rule. lhs). setErasable ();
}



void Symbol::setTerminable ()
{
  if (terminable)
    return;
  terminable = true;
  for (const Rule::Occurrence ro : ruleOccurrences)
    if (ro. rule. isTerminable ())
      const_cast <NonTerminal&> (ro. rule. lhs). setTerminable ();
}



void Symbol::addFirstRO (Rule::Occurrence ro)
{
  if (firstROs. contains (ro))
    return;
  firstROs << ro;
  for (const Rule::Occurrence symbolRo : ruleOccurrences)
    if (symbolRo. ruleCanStartWith ())
      const_cast <NonTerminal&> (symbolRo. rule. lhs). addFirstRO (ro);
}



void Symbol::addLastRO (Rule::Occurrence ro)
{
  if (lastROs. contains (ro))
    return;
  lastROs << ro;
  for (const Rule::Occurrence symbolRo : ruleOccurrences)
    if (symbolRo. ruleCanEndWith ())
      const_cast <NonTerminal&> (symbolRo. rule. lhs). addLastRO (ro);
}



void Symbol::addTerminal (const Terminal* t)
{
  ASSERT (t);

  if (terminals. contains (t))
    return;
  terminals << t;
  for (const Rule::Occurrence ro : ruleOccurrences)
    const_cast <NonTerminal&> (ro. rule. lhs). addTerminal (t);
}



Set<const Symbol*> Symbol::getRulePrevs () const
{
  Set<const Symbol*> s;
  for (const Rule::Occurrence ro : ruleOccurrences)
    if (ro. isFirst ())
      s << nullptr;
    else
      s << * std::prev (ro. rhsIt);
  return s;
}



Set<const Symbol*> Symbol::getRuleNexts () const
{
  Set<const Symbol*> s;
  for (const Rule::Occurrence ro : ruleOccurrences)
    if (ro. isLast ())
      s << nullptr;
    else
      s << * std::next (ro. rhsIt);
  return s;
}




// Terminal

const string Terminal::eotName ("EOT");



Terminal::Terminal (const string &name_arg)
: Symbol (name_arg)
{
  ASSERT (! name. empty ());
  try 
  {
    if (isDigit (name [0]))
    {
      c = str2<Char> (name);
      ASSERT (c < ' ' || c >= 127);
    }
    else
    {    
      ASSERT (name. size () == 3);
      ASSERT (name [0] == Token::quote);
      ASSERT (name [2] == Token::quote);
      c = (Char) name [1];
      ASSERT (c >= ' ')
      ASSERT (c < 127);
    }
    ASSERT (c != eot);
  }
  catch (const exception &e)
  {
    throw runtime_error ("Bad Terminal: " + name + "\n" + e. what ());
  }
}



void Terminal::qc () const
{
  Symbol::qc ();

#ifndef NDEBUG
//if (grammar)
  {
    ASSERT (! erasable);
    ASSERT (terminable);
    ASSERT (firstROs. empty ());
    ASSERT (lastROs. empty ());
    ASSERT (terminals. size () == 1);
    ASSERT (* terminals. begin () == this);
    for (const bool out : Bool)
      ASSERT (Set<Rule::Occurrence> (terminalNeighbors [out]) == ruleOccurrences);
  }
#endif
}



bool Terminal::differentiated (Rule::Occurrence ro,
                               bool out) const
{
  const Occurrence2Terminals& o2t = terminalNeighbors [out];
  ASSERT (! o2t. empty ());
  if (o2t. size () == 1)
    return false;
  const Set<const Terminal*>* s = findPtr (o2t, ro);
  ASSERT (s);
  for (const auto it : o2t)
    if (   ! (it. first == ro)
        && s->intersects (it. second)
       )
      return false;
  return true;      
}




// NonTerminal

const string NonTerminal::sigmaS = "Sigma";



void NonTerminal::qc () const
{
  Symbol::qc ();

  ASSERT (isAlpha (name [0]));

//ASSERT (grammar);
#ifndef NDEBUG
  size_t leftRecursive = 0;
  size_t rightRecursive = 0;
  size_t nonRecursive = 0;
  Set<const Rule*> ruleSet;
  for (const Rule* r : lhsRules)
  {
    ASSERT (r);
    r->qc ();
    ASSERT (& r->lhs == this);
    ruleSet << r;
    if (r->isLeftRecursive ())
      leftRecursive++;
    else if (r->isRightRecursive ())
      rightRecursive++;
    else
      nonRecursive++;
  }
  ASSERT (ruleSet. size () == lhsRules. size ());
  // Non-redundant
  IMPLY (leftRecursive, nonRecursive);
  IMPLY (rightRecursive, nonRecursive);
#endif

//ASSERT (terminal2rules [false]. size () == terminal2rules [true] . size ());
#ifndef NDEBUG
  Set<const Rule*> all;
  for (const bool b : Bool)
    for (const auto it : terminal2rules [b])
    {
      const Terminal* t = it. first;
      ASSERT (t);
    //ASSERT (t->grammar == grammar);
      const VectorPtr<Rule>& rules = it. second;
      Set<const Rule*> s;
      s. insertAll (rules);
      ASSERT (s. size () == rules. size ());
      for (const Rule* r : rules)
      {
        ASSERT (r);
      //ASSERT (r->lhs. grammar == grammar);
        ASSERT (! r->erasable);
        ASSERT (r->isLeftRecursive () == (bool) b);
      }
      all << s;
    }
  Set<const Rule*> s;
  s. insertAll (erasableRules);
  ASSERT (s. size () == erasableRules. size ());
  for (const Rule* r : erasableRules)
  {
    ASSERT (r);
  //ASSERT (r->lhs. grammar == grammar);
    ASSERT (r->erasable);
    ASSERT (! r->isLeftRecursive ());  // Otherwise dummy left-recursivity
    all << r;
  }
  ASSERT (all == ruleSet);
#endif

  IMPLY (erasable, ! erasableRules. empty ());

  // Non-redundant
  ASSERT (! lhsRules. empty ());  
  IMPLY (lhsRules. size () == 1, ! lhsRules [0] -> empty ());
#if 0
  if (isRoot () != isTransient ())
  {
    cout << name << " " << (int) isRoot () << " " << (int) isTransient () << endl;
    ERROR;
  }
#endif
  IMPLY (getEmptyRule (), lhsRules. size () > 1);
  ASSERT (! firstROs. empty ());
  ASSERT (! lastROs. empty ());
}



void NonTerminal::saveText (ostream &os) const
{
  Symbol::saveText (os);

  os << "Complexity: " << getComplexity () << endl;

  for (const Rule* r : lhsRules)
  {
    r->saveText (os);
    os << endl;
  }

  if (erasable)
    os << "Erasable" << endl;

  os << "1st:";
  for (const Rule::Occurrence ro : firstROs)
    os << ' ' << ro. getName ();
  os << endl;

  os << "Last:";
  for (const Rule::Occurrence ro : lastROs)
    os << ' ' << ro. getName ();
  os << endl;

  os << "Terminals:";
  for (const Terminal* t : terminals)
    os << ' ' << t->name;
  os << endl;
  
  os << "RegularStar:";
  for (const Symbol* s : regularStar)
    os << ' ' << s->name;
  os << endl;
}



const Rule* NonTerminal::getEmptyRule () const
{
  for (const Rule* r : lhsRules)
    if (r->empty ())
      return r;
  return nullptr;
}



void NonTerminal::setTerminal2rules ()
{
#ifndef NDEBUG
  for (const bool b : Bool)
    { ASSERT (terminal2rules [b]. empty ()); }
#endif
  ASSERT (erasableRules. empty ());
  for (const Rule* r : lhsRules)
    if (r->erasable)
      erasableRules << r;
    else
      for (const Terminal* t : r->firstTerminals)
        terminal2rules [r->isLeftRecursive ()] [t] << r;
}



void NonTerminal::setRegularStar ()
{
  ASSERT (replacements. contains (this));
  ASSERT (regularStar. empty ());
  for (const Symbol* s : replacements)
    if (const NonTerminal* nt = s->asNonTerminal ())
      for (const Rule* r : nt->lhsRules)
        FOR (size_t, i, r->rhs. size ())
          if (r->singleRhs [i])
            FOR (size_t, j, r->rhs. size ())
              if (   j != i
                  && r->rhs [j] -> replacements. contains (this)
                 )
                regularStar << r->rhs [i] -> replacements;
}



double NonTerminal::getComplexity () const
{
  // Relation to running time ??
  size_t n = erasableRules. size ();
  for (const bool b : Bool)
    for (const auto it : terminal2rules [b])
    {
      const VectorPtr<Rule>& rules = it. second;
      maximize (n, rules. size () /*+ (b ? 0 : erasableRules. size ())*/);
    }
  ASSERT (n);
#if 0
  if (getEmptyRule ())
    n--;
  ASSERT (n);  
#endif
  return log (n);
}




bool NonTerminal::canEnd (const Position &pos) const
{
  for (const auto it : pos. terminal2syntagms)
    if (neighbors [true]. contains (it. first))
      return true;
  return false;
}



const Syntagms* NonTerminal::parse (const Position& pos) const
// Use prefix tree of Rule::Rhs's ??
{
  if (const Syntagms* res = findPtr (pos. nonTerminal2syntagms, this))
    return res;

#if 0
  if (verbose ())
    cout << & pos << ' ' << (char) pos. c << ' ' << name << endl;
#endif
#if 0
  Common_sp::AutoPtr<Offset> ofs;
  if (verbose ())
  {
    ofs. reset (new Offset ());
    ofs->newLn (cout);
    cout << name;
  }
#endif

  const Syntagms& res = findMake (const_cast <Position&> (pos). nonTerminal2syntagms, this);
  ASSERT (res. empty ());

  const VectorPtr<Syntagm> children;
  for (const auto it : pos. terminal2syntagms)
    if (const VectorPtr<Rule>* rules = findPtr (terminal2rules [false], it. first))
      for (const Rule* rule : *rules)
      {
        ASSERT (! rule->isLeftRecursive ());
        rule->parseIt (rule->rhs. begin (), pos, children); 
      }
  for (const Rule* rule : erasableRules)
  {
    ASSERT (! rule->isLeftRecursive ());
    rule->parseIt (rule->rhs. begin (), pos, children); 
  }

  return & res;

#if 0
  // Deterministic parser
  NonTerminalSyntagm* syntagm = nullptr;
  VectorPtr<Syntagm> children;  children. reserve (6);  // PAR

  const Pos pos_init (pos);
  const VectorPtr<Terminal> terminals (pos. getTerminals ());
  for (const Terminal* t : terminals)
  {
    for (const Rule* rule : terminal2rules [false]. at (t))
      if (rule->parse (pos, children))
      {
        syntagm = new NonTerminalSyntagm (pos_init, pos, * rule, children);  // commitment
        break;
      }
    if (syntagm)
      break;
  }
  if (! syntagm)
    for (const Rule* rule : erasableRules)
      if (rule->parse (pos, children))
      {
        syntagm = new NonTerminalSyntagm (pos_init, pos, * rule, children);  // commitment
        break;
      }
  if (! syntagm)
    return nullptr;
   
  bool found = true;
  while (found)
  {
    found = false;
    const VectorPtr<Terminal> terminals (pos. getTerminals ());
    for (const Terminal* t : terminals)
    {
      children. clear ();   
      children << syntagm;
      for (const Rule* rule : terminal2rules [true]. at (t))
        if (rule->parse (pos, children))
        {
          syntagm = new NonTerminalSyntagm (pos_init, pos, * rule, children);  // commitment
          found = true;
          break;
        }
      if (found)
        break;
    }
  }

  return syntagm;    
#endif
}



#if 0
// Letter

const string Letter::space ("space");



void Letter::qc () const
{
  Terminal::qc ();

  ASSERT (grammar);
  ASSERT (! endOfSentence (c));
  ASSERT (! erasable);
}
#endif




#if 0
// SymbolTree

void SymbolTree::qc () const
{
  IMPLY (children. empty (), ! rules. empty ());
#ifndef NDEBUG
  for (const auto it : children)
  {
    ASSERT (it. first);
    const SymbolTree* st = it. second;
    ASSERT (st);
    st->qc ();
  }
#endif

  ASSERT (! rules. contains (nullptr));
}



void SymbolTree::saveText (ostream &os) const
{
  for (const Rule* r : rules)
    os << ' ' << r->num;
  Offset ofs;
  for (const auto it : children)
  {
    ofs. newLn (os);
    os << it. first->name << ':';
    it. second->saveText (os);
  }
}
#endif




// Grammar

const string Grammar::arrowS ("->");



Grammar::Grammar (const string &fName)
: /*Terminal (commentC, nullptr)
,*/ ruleNum (0)
, startSymbol (nullptr)
, eotSymbol (nullptr)
//, symbolTree (* new SymbolTree ())
{
  string startS;

  CharInput in (fName);
  try 
  {
    List<Token> tokens;
    for (;;)
    {
      const Token token (in);
      token. qc ();
      if (token. isDelimiter (commentC))
      {
        in. getLine ();
        continue;
      }
      if (   token. charNum == 0
          || token. empty ()
         )
      {
        if (! tokens. empty ())
        {
          ASSERT (tokens. size () >= 3);
          ASSERT (tokens. front (). type == Token::eName);
          const string lhs (tokens. popFront (). name);
          ASSERT (arrowS. size () == 2);
          ASSERT (tokens. popFront (). isDelimiter (arrowS [0]));
          ASSERT (tokens. popFront (). isDelimiter (arrowS [1]));
        #ifndef NDEBUG
          for (const Token& t : tokens)
            ASSERT (! contains (t. name, commentC));
        #endif
          List<Token>::const_iterator it = tokens. cbegin ();
          addRule (lhs, it, tokens. cend (), false);
          if (startS. empty ())
            startS = lhs;
        }
        tokens. clear ();
      }
      if (token. empty ())
        break;
      tokens << token;
    }
  }
  catch (...)
  {
    cout << "At line " << in. lineNum + 1 << endl;
    throw;
  }

  finish (startS);
}



Grammar::Grammar (const Grammar &other,
                  const VectorPtr<NonTerminal> &newTerminals)
: /*Terminal (other. name, nullptr)
,*/ ruleNum (0)
, startSymbol (nullptr)
, eotSymbol (nullptr)
//, symbolTree (* new SymbolTree ())
{
  ASSERT (! newTerminals. empty ());
#ifndef NDEBUG
  other. qc ();
#endif

  DerivedSymbolGraph dsg (other);      
  for (const NonTerminal* nt : newTerminals)
  {
    SymbolGraph::Node* node = dsg. symbol2node [nt];
    ASSERT (node);
    node->deleteNeighborhood (false);
  }
  const_cast <SymbolGraph::Node*> (dsg. getStartSymbol ()) -> setReachable ();

  for (const Symbol* s : other. symbols)
    if (const NonTerminal* nt = s->asNonTerminal ())
    {
      SymbolGraph::Node* node = dsg. symbol2node [nt];
      ASSERT (node);
      if (! node->reachable)
        continue;
      if (newTerminals. contains (nt))
        continue;
      for (const Rule* r : nt->lhsRules)
      {
        Vector<string> rhs;
        for (const Symbol* child : r->rhs)
          rhs << child->name;
        addRule (r->lhs. name, rhs);
      }
    }

  finish (other. startSymbol->name);
}



Grammar::Grammar (const Grammar &other,
                  const Terminal* universalDelimiter)
: /*Terminal (other. name, nullptr)
,*/ ruleNum (0)
, startSymbol (nullptr)
, eotSymbol (nullptr)
//, symbolTree (* new SymbolTree ())
{
  ASSERT (other. symbols. contains (universalDelimiter));
#ifndef NDEBUG
  other. qc ();
#endif

  for (const Symbol* s : other. symbols)
    if (const NonTerminal* nt = s->asNonTerminal ())
      for (const Rule* r : nt->lhsRules)
      {
        Vector<string> rhs;
        for (const Symbol* child : r->rhs)
          if (child != universalDelimiter)
            rhs << child->name;
        addRule (r->lhs. name, rhs);
      }

  finish (other. startSymbol->name);
}



void Grammar::addRule (const string &lhs,
                       const Vector<string> &rhs)
{
//ASSERT (! startSymbol);

  NonTerminal* nt = getSymbol<NonTerminal> (lhs);
  ASSERT (nt);
  Rule* r = new Rule (ruleNum, *nt, rhs);
  nt->lhsRules << r;
  ruleNum++;
}



void Grammar::addRule (const string &lhs,
                       List<Token>::const_iterator &rhsIt,
                       const List<Token>::const_iterator &rhsEnd,
                       bool optional)
{
  ASSERT (! lhs. empty ());
  ASSERT (lhs != NonTerminal::sigmaS);
  
  string ascii ("   ");
  ASSERT (ascii. size () == 3);
  ascii [0] = Token::quote;
  ascii [2] = ascii [0];

  Vector<string> rhs;
  size_t rhsPart = 1;
  bool stop = false;
  bool prevName = false;
  while (! stop && rhsIt != rhsEnd)
  {
    const Token token (*rhsIt);
    token. qc ();
    rhsIt++;
    ASSERT (! token. empty ());
    bool isName = false;
    switch (token. type)
    {
      case Token::eName: 
        isName = true;
      case Token::eNumber: 
        rhs << token. name; 
        rhsPart++;
        break;
      case Token::eText: 
        for (const char c : token. name)
        {
          ascii [1] = c;
          rhs << ascii;
        }
        rhsPart++;
        break;
      case Token::eDelimiter: 
        {
          ASSERT (token. name. size () == 1);
          const char last = token. name [0];
          switch (last)
          {
            case '[': {
                        const string optionalName (lhs + commentC + toString (rhsPart));
                        addRule (optionalName, rhsIt, rhsEnd, true); 
                        addRule (optionalName, Vector<string> ());
                        rhs << optionalName;
                        rhsPart++;
                      }
                      break;
            case ']': ASSERT (optional);   
                      stop = true;                   
                      break;
            case '*':
            case '+': {
                        ASSERT (prevName);
                        ASSERT (! rhs. empty ());
                        const string name (rhs. back ());
                        const string nameMod (name + last);
                        rhs. back () = nameMod;
                        if (! string2Symbol [nameMod])
                        {
                          getSymbol<NonTerminal> (nameMod);
                          List<Token> tokens;
                          List<Token>::const_iterator it = tokens. cbegin ();
                          switch (last)
                          {
                            case '+': tokens << name << name << '*';
                                      break;
                            case '*': addRule (nameMod, it, tokens. cend (), false);
                                      tokens << name << '+';
                                      break;
                            default: ERROR;
                          }
                          it = tokens. cbegin ();
                          addRule (nameMod, it, tokens. cend (), false);
                        }
                      }
                      break;
            default:  ERROR_MSG ("Unknown special grammar symbol: \'" + token. name + "\'");
          }
        }
        break;
      default: ERROR;
    }    
    prevName = isName;
  }
  ASSERT (stop == optional);
  
  addRule (lhs, rhs);
}



void Grammar::finish (const string &startS)
{
  ASSERT (! string2Symbol. empty ());

  ASSERT (! eotSymbol);
  eotSymbol = new Terminal ();
  ASSERT (! string2Symbol [eotSymbol->name]);
  string2Symbol [eotSymbol->name] = const_cast <Terminal*> (eotSymbol);
  symbols << eotSymbol;

  ASSERT (! startSymbol);
  startSymbol = new NonTerminal ();
  ASSERT (! string2Symbol [startSymbol->name]);
  string2Symbol [startSymbol->name] = const_cast <NonTerminal*> (startSymbol);
  symbols << startSymbol;
  {
    Symbol* s = string2Symbol [startS /*name*/];
    ASSERT (s);
    Vector<string> rhs;  rhs. reserve (3);
    rhs << eotSymbol->name << s->name << eotSymbol->name;
    addRule (startSymbol->name, rhs);
  }

  // Rule::finish()
  {
    VectorPtr<NonTerminal> nonTerminals;  
    for (const Symbol* s : symbols)
      if (const NonTerminal* nt = s->asNonTerminal ())
        nonTerminals << nt;
    for (const NonTerminal* nt : nonTerminals)
      for (const Rule* r : nt->lhsRules)
        const_cast <Rule*> (r) -> finish (*this);  // updates symbols
  }

  ASSERT (symbols. size () == string2Symbol. size ());
  string2Symbol. clear ();

  setTerminable ();
  setErasable ();
  setRuleSingleRhs ();
  setFirstROs ();
  setLastROs ();
  setParsingTable ();
  setTerminals ();
  for (const bool out : Bool)
    setNeighbors (out);
//setSymbolTree ();
  setReplacements ();
  setRegularStar ();
}



void Grammar::qc () const
{
//Terminal::qc ();

//ASSERT (grammar != this);
  ASSERT (! symbols. empty ());
  ASSERT (startSymbol);
//ASSERT (symbols [0] == startSymbol);

  // eotSymbol
  ASSERT (eotSymbol);
//ASSERT (eotSymbol->ruleOccurrences. empty ());

#ifndef NDEBUG
  Set<const Symbol*> symbolSet;
  Set<Char> charSet;
  size_t terminals = 0;
  for (const Symbol* s : symbols)
  {
    ASSERT (s);
  //ASSERT (s->grammar == this);
    try { s->qc (); }
      catch (...) { cout << s->name << endl; throw; }
    symbolSet << s;
    if (const Terminal /*Letter*/ * let = s->asTerminal /*Letter*/ ())
    {
      terminals++;
      charSet << let->c;
    }
    // Non-redundant
    IMPLY (/*s != eotSymbol &&*/ s->ruleOccurrences. empty (), s == startSymbol);  
  }
  ASSERT (symbolSet. size () == symbols. size ());
  ASSERT (charSet. size () == terminals);

  for (const auto it : char2terminal)
    ASSERT (symbols. contains (it. second));      

  // Valid CF-grammar
  if (const Symbol* s = getNonReachable ())
    ERROR_MSG (s->name + " is not reachable");
  if (const Symbol* s = getNonTerminable ())
    ERROR_MSG (s->name + " is not terminable");
#endif

#if 0
#ifndef NDEBUG
  for (const Rule* r : symbolTree. rules)
    ASSERT (r->erasable);
#endif
  symbolTree. qc ();
#endif
}



void Grammar::saveText (ostream &os) const
{
  ASSERT (startSymbol);

  os << "Start: " << startSymbol->name << endl;
  os << endl;

  os << "Non-terminals:" << endl;
  for (const Symbol* s : symbols)
    if (const NonTerminal* nt = s->asNonTerminal ())
    {
      nt->saveText (os);
      os << endl;
    }
  os << endl;

  os << "Terminals:" << endl;
  for (const Symbol* s : symbols)
    if (const Terminal* t = s->asTerminal ())
    {
      t->saveText (os);
      os << endl;
    }
  os << endl;

  os << "Complexity: " << getComplexity () << endl;

  os << endl << "Unique prevs: " << endl;
  map<const Symbol* /*prev*/, Set<const Symbol*> /*nexts*/> mainNexts;
  for (const Symbol* s : symbols)
  {
    const Set<const Symbol*> prevs (s->getRulePrevs ());
    if (prevs. size () == 1)
      if (const Symbol* prev = * prevs. begin ())
      {
        mainNexts [prev] << s;
        os << s->name << " after " << prev->name << endl;
      }
  }
  for (const auto it : mainNexts)
  {
    Set<const Symbol*> nexts (it. first->getRuleNexts ());
    nexts. erase (nullptr);
  //os << it. first->name << " " << it. second. size () << " " << nexts. size () << endl;  
    if (it. second == nexts)
      os << "Implied: " << it. first->name << endl;    
  }

  os << endl << "Unique nexts: " << endl;
  map<const Symbol* /*next*/, Set<const Symbol*> /*prevs*/> mainPrevs;
  for (const Symbol* s : symbols)
  {
    const Set<const Symbol*> nexts (s->getRuleNexts ());
    if (nexts. size () == 1)
      if (const Symbol* next = * nexts. begin ())
      {
        mainPrevs [next] << s;
        os << s->name << " before " << Symbol::getName (next) << endl;
      }
  }
  for (const auto it : mainPrevs)
  {
    Set<const Symbol*> prevs (it. first->getRulePrevs ());
    prevs. erase (nullptr);
    if (it. second == prevs)
      os << "Implied: " << it. first->name << endl;    
  }

#if 0
  os << endl;
  os << "Symbol tree:" << endl;
  symbolTree. saveText (os);
  os << endl;
#endif
}



void Grammar::setTerminable ()
{
  for (const Symbol* s : symbols)
    if (const Terminal* t = s->asTerminal ())
      const_cast <Terminal*> (t) -> setTerminable ();
    else 
    {
      const NonTerminal* nt = s->asNonTerminal ();
      ASSERT (nt);
      if (nt->getEmptyRule ())
        const_cast <NonTerminal*> (nt) -> setTerminable ();
    }
}



void Grammar::setErasable ()
{
  for (const Symbol* s : symbols)
    if (const NonTerminal* nt = s->asNonTerminal ())
      if (nt->getEmptyRule ())
        const_cast <NonTerminal*> (nt) -> setErasable ();
}



void Grammar::setRuleSingleRhs ()
{
  for (const Symbol* s : symbols)
    if (const NonTerminal* nt = s->asNonTerminal ())
      for (const Rule* r : nt->lhsRules)
        const_cast <Rule*> (r) -> setSingleRhs ();
}



void Grammar::setFirstROs ()
{ 
  for (const Symbol* s : symbols)
    for (const Rule::Occurrence ro : s->ruleOccurrences)
      if (ro. ruleCanStartWith ())
        const_cast <NonTerminal&> (ro. rule. lhs). addFirstRO (ro);
}



void Grammar::setLastROs ()
{ 
  for (const Symbol* s : symbols)
    for (const Rule::Occurrence ro : s->ruleOccurrences)
      if (ro. ruleCanEndWith ())
        const_cast <NonTerminal&> (ro. rule. lhs). addLastRO (ro);
}



void Grammar::setParsingTable ()
{
  for (const Symbol* s : symbols)
    if (const NonTerminal* nt = s->asNonTerminal ())
    {
      for (const Rule* r : nt->lhsRules)
        const_cast <Rule*> (r) -> setFirstTerminalsErasable ();
      const_cast <NonTerminal*> (nt) -> setTerminal2rules ();
    }
}



void Grammar::setTerminals ()
{
  for (const Symbol* s : symbols)
    if (const Terminal* t = s->asTerminal ())
      const_cast <Terminal*> (t) -> addTerminal (t);
}



void Grammar::setNeighbors (bool out)
{
  const FollowROGraph frg (*this);
  frg. qc ();
//frg. print (cout);
  for (const DiGraph::Node* node : frg. nodes)
  {
    const Rule::Occurrence ro = static_cast <const ROGraph::Node*> (node) -> ro;
    const Symbol* symbol = ro. getSymbol ();
    if (const Terminal* t = symbol->asTerminal ())
      const_cast <Terminal*> (t) -> terminalNeighbors [out] [ro] = Set<const Terminal*> ();
    for (const DiGraph::Arc* arc : node->arcs [out])
    {
      const Symbol* neighbor = static_cast <const ROGraph::Node*> (arc->node [out]) -> ro. getSymbol ();
      const_cast <Symbol*> (symbol) -> neighbors [out] << neighbor;
      if (const Terminal* t = symbol->asTerminal ())
        if (const Terminal* neighborT = neighbor->asTerminal ())
          const_cast <Terminal*> (t) -> terminalNeighbors [out] [ro] << neighborT;
    }
  }
}



void Grammar::setReplacements ()
{
  SingleDerivedSymbolGraph sdsg (*this);
  sdsg. qc ();
  for (const DiGraph::Node* node : sdsg. nodes)
  {
    const SymbolGraph::Node* parent = static_cast <const SymbolGraph::Node*> (node);
    ASSERT (parent->symbol. replacements. empty ());
    sdsg. clearReachable ();
		const_cast <SymbolGraph::Node*> (parent) -> setReachable ();
    for (const DiGraph::Node* childNode : sdsg. nodes)
    {
      const SymbolGraph::Node* child = static_cast <const SymbolGraph::Node*> (childNode);
      if (child->reachable)
        const_cast <Symbol&> (parent->symbol). replacements << & child->symbol;
    }
  }
}



void Grammar::setRegularStar ()
{
  for (const Symbol* s : symbols)
    if (const NonTerminal* nt = s->asNonTerminal ())
      const_cast <NonTerminal*> (nt) -> setRegularStar ();
}



#if 0
void Grammar::setSymbolTree ()
{
  ASSERT (symbolTree. children. empty ());
  for (const Symbol* s : symbols)
    if (const NonTerminal* nt = s->asNonTerminal ())
      for (const Rule* r : nt->lhsRules)
      {
        SymbolTree* st = const_cast <SymbolTree*> (& symbolTree);
        for (const Symbol* rs : r->rhs)
        {
          ASSERT (st);
          SymbolTree* st_new = const_cast <SymbolTree*> (st->children [rs]);
          if (! st_new)
          {
            st_new = new SymbolTree ();
            st->children [rs] = st_new;
          }
          st = st_new;
        }
        st->rules << r;
      }
}
#endif



const Symbol* Grammar::getNonReachable () const
{
  const DerivedSymbolGraph dsg (*this);    
  const_cast <SymbolGraph::Node*> (dsg. getStartSymbol ()) -> setReachable ();
  if (const SymbolGraph::Node* n = dsg. getNonReachable ())
    return & n->symbol;
  return nullptr;
}



const Symbol* Grammar::getNonTerminable () const
{
  for (const Symbol* s : symbols)
    if (! s->terminable)
      return s;
  return nullptr;
}



void Grammar::setChar2terminal ()
{
  for (const Symbol* s : symbols)
    if (const Terminal* t = s->asTerminal ())
      char2terminal [t->c] = t;
}



const Symbol* Grammar::findParseCycle () const
{
  FirstDerivedSymbolGraph fdsg (*this);
  fdsg. scc ();
  for (DiGraph::Node* n : fdsg. nodes)
    if (n->getConnectedComponent () != n)  // There exists a cycle which is not a self-loop
      return  & static_cast <FirstDerivedSymbolGraph::Node*> (n) -> symbol;
  return nullptr;
}



const Rule* Grammar::getLeftRecursiveErasable () const
{
  for (const Symbol* s : symbols)
    if (const NonTerminal* nt = s->asNonTerminal ())
      for (const Rule* r : nt->lhsRules)
        if (r->isLeftRecursive () && r->erasable)
          return r;
  return nullptr;
}



void Grammar::findAmbiguity () const throw (Ambiguity)
{
  for (const Symbol* s : symbols)
    if (const NonTerminal* nt = s->asNonTerminal ())
      if (! nt->regularStar. empty ())
        for (const Symbol* nextSymbol : nt->neighbors [true])
          if (const NonTerminal* next = nextSymbol->asNonTerminal ())
          {
            Set<const Symbol*> intersection (next->regularStar);
            intersection. intersect (nt->regularStar);
            if (! intersection. empty ())
              throw Ambiguity (intersection. front () -> name + "* tandem");
          }
  // tandem of the same erasable symbol ??
  // Other cases ??
}



void Grammar::prepare () throw (StdParserError, Ambiguity)
{
  if (const Symbol* s = findParseCycle ())
    throw StdParserError ("First symbol cycle: " + s->name);

  if (const Rule* r = getLeftRecursiveErasable ())
  {
    r->print (cout);
    throw StdParserError ("Left-recursive erasable rule");
  }

  setChar2terminal ();
  
  findAmbiguity ();
}



double Grammar::getComplexity () const
{
  double comp = 0;
  for (const Symbol* s : symbols)
    if (const NonTerminal* nt = s->asNonTerminal ())
      comp += nt->getComplexity ();
  return comp;
}



const Syntagms& Grammar::parseSentence (Sentence &sentence) const
{ 
  ASSERT (! sentence. seq. empty ());
  const Syntagms* syntagms = startSymbol->parse (sentence. seq [0 /*1*/]);
  ASSERT (syntagms);
  sentence. qc ();
  for (const Syntagm* syntagm : *syntagms)
    const_cast <Syntagm*> (syntagm) -> setRight ();
  return *syntagms;
}



VectorPtr<NonTerminal> Grammar::getCutSymbols () const
{
  DerivedSymbolGraph dsg (*this);    

  VectorPtr<NonTerminal> res;
  for (DiGraph::Node* node : dsg. nodes)
  {
    const VectorPtr<DiGraph::Node> children (node->getChildren ());
    if (children. empty ())
      continue;

    const NonTerminal* nt = static_cast <SymbolGraph::Node*> (node) -> symbol. asNonTerminal ();
    ASSERT (nt);
    if (nt == startSymbol)
      continue;

    node->deleteNeighborhood (false);

    dsg. connectedComponents ();
    const DiGraph::Node* root = const_cast <SymbolGraph::Node*> (dsg. getStartSymbol ()) -> getConnectedComponent ();
    size_t terminalNum = 0;
    for (DiGraph::Node* node1 : dsg. nodes)
      if (const Terminal* t = static_cast <SymbolGraph::Node*> (node1) -> symbol. asTerminal ())
        if (/*t != eotSymbol &&*/ node1->getConnectedComponent () != root)
        {
          ASSERT (nt->terminals. contains (t));
          terminalNum++;
        }
    if (nt->terminals. size () == terminalNum)
      res << nt;

    // Restore neighborhood
    for (const DiGraph::Node* child : children)
      new DiGraph::Arc (const_cast <DiGraph::Node*> (child), node);
  }

  return move (res);
}



const Terminal* Grammar::getUniversalDelimiter () const
{
  const FollowSymbolGraph fsg (*this);
  fsg. qc ();
//fsg. print (cout);
  const Terminal* ud = nullptr;
  for (const DiGraph::Node* node : fsg. nodes)
  {
    const SymbolGraph::Node* symbolNode = static_cast <const SymbolGraph::Node*> (node);
    if (const Terminal* t = symbolNode->symbol. asTerminal ())
    {
      const Terminal* prev = symbolNode->getNextTerminal (false);
      const Terminal* next = symbolNode->getNextTerminal (true);
      if (prev && next && prev == next && prev != t)
        if (ud && ud != prev)
          return nullptr;
        else
          ud = prev;
      else
        if (ud && ud != t)
          return nullptr;
        else
          ud = t;
    }
  }
  return ud;
}



void Grammar::splitTerminals () const
{
  for (const Symbol* symbol : symbols)
    if (const Terminal* t = symbol->asTerminal ())
      for (const bool out : Bool)
        for (const Rule::Occurrence ro : t->ruleOccurrences)
          if (t->differentiated (ro, out))
            cout << " " << ro. getName () << '/' << (out ? "next" : "prev");  // ??
}



void Grammar::terminals4scc () const
{
  DerivedSymbolGraph dsg (*this);
  dsg. scc ();
//dsg. contractScc ();
//dsg. print (cout);

  map<const DiGraph::Node* /*scc*/, Set<const Symbol*> /*terminals for scc*/> scc2terminals;
  for (const DiGraph::Node* n : dsg. nodes)
  {
    const DiGraph::Node* scc = n->scc;
    if (scc != n)  // Non-singleton
      scc2terminals [scc] = Set<const Symbol*> ();
  }
  for (const DiGraph::Node* n : dsg. nodes)
  {
    const DiGraph::Node* scc = n->scc;
    if (const Set<const Symbol*>* terminals = findPtr (scc2terminals, scc))
      for (const DiGraph::Arc* arc : n->arcs [false])
      {
        const DiGraph::Node* other = arc->node [false];
        if (other->scc != scc)
          const_cast <Set<const Symbol*>&> (*terminals) << & static_cast <const SymbolGraph::Node*> (other) -> symbol;
      }
  }

  // Only print ??
  for (const auto it : scc2terminals)
  {
    cout << static_cast <const SymbolGraph::Node*> (it. first) -> symbol. name << ":";
    for (const Symbol* s : it. second)
      cout << ' ' << s->name;
    cout << endl;
  }
}




// SymbolGraph::Node

void SymbolGraph::Node::setReachable ()
{
  if (reachable)
    return;
  reachable = true;
	for (Arc* arc : arcs [false])
		static_cast <Node*> (arc->node [false]) -> setReachable ();
}



const Terminal* SymbolGraph::Node::getNextTerminal (bool out) const
{
  const Terminal* res = nullptr;
  for (const DiGraph::Arc* arc : arcs [out])
    if (const Terminal* t = static_cast <const SymbolGraph::Node*> (arc->node [out]) -> symbol. asTerminal ())
    {
      if (res)
        return nullptr;
      else
        res = t;
    }
  return res;
}




// SymbolGraph

const SymbolGraph::Node* SymbolGraph::getNonReachable () const
{
  for (DiGraph::Node* n : nodes)
  {
    const SymbolGraph::Node* n_ = static_cast<SymbolGraph::Node*> (n);
    if (! n_->reachable)
      return n_;
  }
  return nullptr;
}




// DerivedSymbolGraph

void DerivedSymbolGraph::initArcs ()
{
  for (const Symbol* s : grammar. symbols)
    if (const NonTerminal* nt = s->asNonTerminal ())
    {
      Node* parent = symbol2node [nt];
      ASSERT (parent);
      for (const Rule* r : nt->lhsRules)
        for (const Symbol* sChild : r->rhs)
        {
          Node* child = symbol2node [sChild];
          ASSERT (child);
          if (! parent->isIncident (child, false))
            new DiGraph::Arc (child, parent);
        }
    }
}




// FirstDerivedSymbolGraph

void FirstDerivedSymbolGraph::initArcs ()
{
  for (const Symbol* s : grammar. symbols)
    if (const NonTerminal* nt = s->asNonTerminal ())
    {
      Node* parent = symbol2node [nt];
      ASSERT (parent);
      for (const Rule* r : nt->lhsRules)
        for (const Symbol* sChild : r->rhs)
        {
          Node* child = symbol2node [sChild];
          ASSERT (child);
          if (! parent->isIncident (child, false))
            new DiGraph::Arc (child, parent);
          if (! sChild->erasable)
            break;
        }
    }
}




// SingleDerivedSymbolGraph

void SingleDerivedSymbolGraph::initArcs ()
{
  for (const Symbol* s : grammar. symbols)
    if (const NonTerminal* nt = s->asNonTerminal ())
    {
      Node* parent = symbol2node [nt];
      ASSERT (parent);
      for (const Rule* r : nt->lhsRules)
        FOR (size_t, i, r->rhs. size ())
          if (r->singleRhs [i])
          {
            Node* child = symbol2node [r->rhs [i]];
            ASSERT (child);
            if (! parent->isIncident (child, false))
              new DiGraph::Arc (child, parent);
          }
    }
}




// FollowSymbolGraph

void FollowSymbolGraph::initArcs ()
{
  // use FollowROGraph ??
  for (const Symbol* s : grammar. symbols)
    if (const NonTerminal* nt = s->asNonTerminal ())
      for (const Rule* r : nt->lhsRules)
      {
        const VectorPtr<Symbol>& rhs = r->rhs;
        FOR (size_t, j, rhs. size ())
          FOR (size_t, i, j)
          {
            bool erasable = true;
            FOR_START (size_t, k, i + 1, j)
              if (! rhs [k] -> erasable)
              {
                erasable = false;
                break;
              }
            if (! erasable)
              continue;
            for (const Symbol* prevS : rhs [i] -> getLastSymbols ())
            {
              Node* prev = symbol2node [prevS];
              ASSERT (prev);
              for (const Symbol* nextS : rhs [j] -> getFirstSymbols ())
              {
                Node* next = symbol2node [nextS];
                ASSERT (next);
                if (! prev->isIncident (next, true))
                  new DiGraph::Arc (prev, next);
              }
            }
          }
      }
}




// ROGraph::Node

Set<const Terminal*> ROGraph::Node::getNextTerminals (bool out) const
{
  Set<const Terminal*> res;
  for (const DiGraph::Arc* arc : arcs [out])
    if (const Terminal* t = static_cast <const ROGraph::Node*> (arc->node [out]) -> ro. getSymbol () -> asTerminal ())
      res << t;
  return res;
}




// FollowROGraph

void FollowROGraph::initArcs ()
{
  for (const Symbol* s : grammar. symbols)
    if (const NonTerminal* nt = s->asNonTerminal ())
      for (const Rule* r : nt->lhsRules)
      {
        const VectorPtr<Symbol>& rhs = r->rhs;
        FOR (size_t, j, rhs. size ())
          FOR (size_t, i, j)  // i < j
          {
            bool erasable = true;
            FOR_START (size_t, k, i + 1, j)  // i < k < j
              if (! rhs [k] -> erasable)
              {
                erasable = false;
                break;
              }
            if (! erasable)
              continue;
            for (const Rule::Occurrence prevRo : Rule::Occurrence (*r, i). getLastROs () /*rhs [i] -> lastROs*/)
            {
              Node* prev = ro2node [prevRo];
              ASSERT (prev);
              for (const Rule::Occurrence nextRo : Rule::Occurrence (*r, j). getFirstROs () /*rhs [j] -> firstROs*/)
              {
                Node* next = ro2node [nextRo];
                ASSERT (next);
                if (! prev->isIncident (next, true))
                  new DiGraph::Arc (prev, next);
              }
            }
          }
      }
}



}  // namespace



/* TO DO: ??
  canEnd() for last rule occurrence
*/
