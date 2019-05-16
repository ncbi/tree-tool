// cfgrammar.cpp

#undef NDEBUG
#include "../common.inc"

#include "cfgrammar.hpp"



namespace Cfgr
{


// Sentence

Sentence::Sentence (const Grammar &grammar,
                    const string &s) /*throw (BadPos)*/
{
  seq. resize (s. size () + 2);
  seq [0]. init (grammar, eot). right = true;
  FOR (size_t, i, s. size ())
    try 
    { 
      const char c = s [i];
      ASSERT (c != eot);
      seq [i + 1]. init (grammar, c); 
    }
    catch (const exception &e) 
      { throw BadPos (i, e. what ()); }
  seq [s. size () + 1]. init (grammar, eot). right = true;
}



void Sentence::qc () const
{
  if (! qc_on)
    return;

  QC_ASSERT (! seq. empty ());
  FOR (size_t, i, seq. size ())
  {
    const Position& pos = seq [i];
    if (pos. c == eot)
      { QC_ASSERT (i == 0 || i == seq. size () - 1); }
    else
      pos. qc ();
  }
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



const NonTerminalSyntagm* Sentence::getLastLongestSyntagm () const
{
  const Syntagm* lastlongestSyntagm = nullptr;
  ptrdiff_t diff_max = 0;
  size_t size_max = 0;
  for (const Position& pos : seq) 
    for (const auto& it : pos. nonTerminal2syntagms)
      for (const Syntagm* syntagm : *it. second)
      {
        const ptrdiff_t diff = & syntagm->end - & seq [0];
        ASSERT (diff > 0);
        if (maximize (diff_max, diff))
          size_max = 0;
        if (   diff_max == diff 
            && maximize (size_max, syntagm->size ())
           )
          lastlongestSyntagm = syntagm;
      }
        
  if (! lastlongestSyntagm)
    return nullptr;  
  const NonTerminalSyntagm* res = lastlongestSyntagm->asNonTerminalSyntagm ();
  ASSERT (res);
  return res;
}



void Sentence::reportError (ostream &os) const
{
  const NonTerminalSyntagm* lastlongestSyntagm = getLastLongestSyntagm ();
  if (! lastlongestSyntagm)
    return;
  size_t index = 0;
  bool reached = false;
  for (const Position& pos : seq)
  {
    if (& pos == & lastlongestSyntagm->end)
    {
      reached = true;
      os << '|';
    }
    os << (char) pos. c;
    if (! reached)
      index++;
  }
  ASSERT (reached);
  os << endl;
  os << "  At: " << index << endl;
  
  const Position* pos = & seq. at (1);
  while (pos <= & seq. back ())
  {
    const NonTerminalSyntagm* longestSyntagm = pos->getLongestSyntagm ();
    if (! longestSyntagm)
      break;
    longestSyntagm->print (cout);
    pos = & longestSyntagm->end;
  }
  cout << endl;
}



void Sentence::printWrongSyntagms (ostream &os,
                                   size_t size_min) const
{
  for (const Position& pos : seq) 
    for (const auto& it : pos. nonTerminal2syntagms)
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
    for (const auto& it : pos. nonTerminal2syntagms)
      for (const Syntagm* syntagm : *(it. second))
        if (   syntagm->right 
            && syntagm->symbol. name == symbolName
           )
          vec << syntagm->asNonTerminalSyntagm ();
  return vec;
}




// Syntagm

void Syntagm::qc () const
{
  if (! qc_on)
    return;
  QC_ASSERT (begin. c != eot);
}



void Syntagm::saveText (ostream &os) const
{ 
  os << symbol. name << ": " << str (); 
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
  if (! qc_on)
    return;
  Syntagm::qc ();

  QC_ASSERT (symbol. asTerminal ());
  QC_ASSERT (findPtr (begin. terminal2syntagms, symbol. asTerminal ()));
  QC_ASSERT (& begin != & end);
}



string TerminalSyntagm::str () const 
{ 
  return string (1, begin. c); 
}




// NonTerminalSyntagm

NonTerminalSyntagm::NonTerminalSyntagm (const Position& begin_arg,
                                        const Position& end_arg,
                                        const Rule &rule_arg,
                                        const VectorPtr<Syntagm> &children_arg)
: Syntagm (begin_arg, end_arg, rule_arg. lhs)
, rule (rule_arg)
, children (children_arg)  // move ??
{
  const Syntagms* syntagms = findPtr (begin. nonTerminal2syntagms, & rule.lhs);
  ASSERT (syntagms);
  * const_cast <Syntagms*> (syntagms) << this;
}



void NonTerminalSyntagm::qc () const
{
  if (! qc_on)
    return;
  Syntagm::qc ();

  QC_ASSERT (symbol. asNonTerminal ());
  QC_ASSERT (findPtr (begin. nonTerminal2syntagms, symbol. asNonTerminal ()));
  QC_ASSERT (& symbol == & rule. lhs);
  QC_ASSERT (children. size () == rule. rhs. size ());
  const Syntagm* prev = nullptr;
  FOR (size_t, i, children. size ())
  {
    const Syntagm* syntagm = children [i];
    QC_ASSERT (syntagm);
    QC_ASSERT (& syntagm->symbol == rule. rhs [i]);
    QC_IMPLY (i == 0, & syntagm->begin == & begin);
    QC_IMPLY (i == children. size () - 1, & syntagm->end == & end);
    QC_IMPLY (prev, & prev->end == & syntagm->begin);
  //syntagm->qc ();   // Invoked in Position::qc()
    prev = syntagm;
  }
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



string NonTerminalSyntagm::str () const
{
  string vec;  vec. reserve (size ());
  for (const Syntagm* syntagm : children)
    vec += syntagm->str ();
/*
  for (const Position* pos = & begin; pos != & end; pos++)
    vec << pos->c;    
*/
  return vec;
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
  if (! qc_on)
    return;
    
  QC_ASSERT (! terminal2syntagms. empty ()); 
  const Syntagm* prev = nullptr;
  for (const auto mapIt : terminal2syntagms)
  {
    QC_ASSERT (mapIt. first);
    QC_ASSERT (mapIt. second);
    Set<const Syntagm*> dSet;
    for (const Syntagm* syntagm : * mapIt. second)
    {
      QC_ASSERT (syntagm);
      syntagm->qc ();
      QC_ASSERT (syntagm->asTerminalSyntagm ());
      QC_ASSERT (& syntagm->begin == this);
      QC_ASSERT (c == syntagm->symbol. asTerminal () -> c);
      QC_ASSERT (mapIt. first == & syntagm->symbol);
      dSet << syntagm;
      QC_IMPLY (prev, & prev->begin == & syntagm->begin);
      prev = syntagm;
    }
    QC_ASSERT (dSet. size () == mapIt. second->size ());
  }
  for (const auto mapIt : nonTerminal2syntagms)
  {
    QC_ASSERT (mapIt. first);
    QC_ASSERT (mapIt. second);
    Set<const Syntagm*> dSet;
    for (const Syntagm* syntagm : * mapIt. second)
    {
      QC_ASSERT (syntagm);
      syntagm->qc ();
      QC_ASSERT (syntagm->asNonTerminalSyntagm ());
      QC_ASSERT (& syntagm->begin == this);
      QC_ASSERT (mapIt. first == & syntagm->symbol);
      dSet << syntagm;
      QC_IMPLY (prev, & prev->begin == & syntagm->begin);
      prev = syntagm;
    }
    QC_ASSERT (dSet. size () == mapIt. second->size ());
  }
}



TerminalSyntagm& Position::init (const Grammar &grammar,
                                 char c_arg)
{
  c = c_arg;
/*if (isspace (c))
    c = ' '; */
  if (const Terminal* t = findPtr (grammar. char2terminal, c))
    return * new TerminalSyntagm (*this, *t); 
  else
    throw runtime_error ("Bad character");
}



const NonTerminalSyntagm* Position::getLongestSyntagm () const
{
  const Syntagm* longestSyntagm = nullptr;
  size_t size_max = 0;
  for (const auto& it : nonTerminal2syntagms)
    for (const Syntagm* syntagm : *it. second)
      if (maximize (size_max, syntagm->size ()))
        longestSyntagm = syntagm;
        
  if (! longestSyntagm)
    return nullptr;  
  const NonTerminalSyntagm* res = longestSyntagm->asNonTerminalSyntagm ();
  ASSERT (res);
  return res;
}




// Rule::Occurrences

Set<const Symbol*> Rule::Occurrences::getSymbols () const
{
  Set<const Symbol*> res;
  for (const Rule::Occurrence ro : *this)
    res << ro. getSymbol ();
  return res;
}




// Rule::Occurrence

void Rule::Occurrence::qc () const
{
  if (! qc_on)
    return;
  QC_ASSERT (*rhsIt);
  QC_ASSERT (rhsIt - rule. rhs. begin () >= 0);
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
  return false;
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
  return false;
}



Rule::Occurrences Rule::Occurrence::getFirstROs () const
{
  Occurrences res (getSymbol () -> firstROs);
  res << *this;
  return res;
}



Rule::Occurrences Rule::Occurrence::getLastROs () const
{
  Occurrences res (getSymbol () -> lastROs);
  res << *this;
  return res;
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
  if (! qc_on)
    return;
  Root::qc ();

  QC_ASSERT (! rhs. contains (nullptr));
  QC_IMPLY (empty (), isErasable ());
  QC_IMPLY (isErasable (), isTerminable ());
  QC_IMPLY (isErasable (), lhs. erasable);
  QC_IMPLY (isErasable (), erasable);
  QC_ASSERT (firstTerminals. empty () == erasable);
  QC_ASSERT (! firstTerminals. contains (nullptr));
  QC_ASSERT (singleRhs. size () == rhs. size ());
#if 0
  for (const Symbol* s : rhs)
    QC_ASSERT (lhs. grammar == s->grammar);
#endif

  // Non-redundant
  QC_IMPLY (rhs. size () == 1, & lhs != rhs [0]);
}



void Rule::saveText (ostream &os) const
{ 
  os << "Rule# " << num << ":   " << lhs. name << ' ' << Grammar::arrowS;
  for (const Symbol* s : rhs)  
    os << ' ' << s->name; 
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

  return ros;
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

  return ros;
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
      for (const auto& it : pos. terminal2syntagms)
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
  if (! qc_on)
    return;
  Named::qc ();
    
  QC_ASSERT (! isLeft  (name, "["));
  QC_ASSERT (! isRight (name, "]"));

  QC_IMPLY (erasable, terminable);

//IMPLY (! grammar, asGrammar());
  for (const Rule::Occurrence ro : ruleOccurrences)
  {
    ro. qc ();
    QC_ASSERT (ro. getSymbol () == this);
  }

  QC_ASSERT (! terminals. contains (nullptr));
  for (const Rule::Occurrence ro : firstROs)
    ro. qc ();
  for (const Rule::Occurrence ro : lastROs)
    ro. qc ();
  bool found = false;
  for (const Symbol* s : getFirstSymbols ())
  {
    QC_ASSERT (s);
  //QC_ASSERT (s->grammar == grammar);
    if (const Terminal* t = s->asTerminal ())
    {
      found = true;
      QC_ASSERT (terminals. contains (t));
    }
  }
  // Non-redundant
  QC_ASSERT (/*(bool) grammar ==*/ found)

  found = false;
  for (const Symbol* s : getLastSymbols ())
  {
    QC_ASSERT (s);
  //QC_ASSERT (s->grammar == grammar);
    if (const Terminal* t = s->asTerminal ())
    {
      found = true;
      QC_ASSERT (terminals. contains (t));
    }
  }
  // Non-redundant
  QC_ASSERT (/*(bool) grammar ==*/ found)

  QC_ASSERT (! terminals. contains (nullptr));

  // Non-redundant
  for (const bool out : {false, true})
  {
    QC_IMPLY (name != NonTerminal::sigmaS, ! neighbors [out]. empty ());
    QC_ASSERT (! neighbors [out]. contains (nullptr));
  }


  QC_ASSERT (! replacements. contains (nullptr));
  QC_ASSERT (replacements. contains (this));
}



void Symbol::saveText (ostream &os) const
{
  os << name << ":" << endl;
  
  for (const bool out : {false, true})
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

Terminal::Terminal (const string &name_arg)
: Symbol (name_arg)
{
  ASSERT (! name. empty ());
  try 
  {
    if (isDigit (name [0]))
    {
      const int i = str2<int> (name);
      ASSERT (i > 0);
      ASSERT (i < ' ' || i >= 127);
      ASSERT (i < 256);
      c = (char) i;
    }
    else
    {    
      ASSERT (name. size () == 3);
      ASSERT (name [0] == '\"');
      ASSERT (name [2] == '\"');
      c = name [1];
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
  if (! qc_on)
    return;
  Symbol::qc ();

//if (grammar)
  {
    QC_ASSERT (! erasable);
    QC_ASSERT (terminable);
    QC_ASSERT (firstROs. empty ());
    QC_ASSERT (lastROs. empty ());
    QC_ASSERT (terminals. size () == 1);
    QC_ASSERT (* terminals. begin () == this);
    for (const bool out : {false, true})
      QC_ASSERT (Set<Rule::Occurrence> (terminalNeighbors [out]) == ruleOccurrences);
  }
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
  for (const auto& it : o2t)
    if (   ! (it. first == ro)
        && s->intersects (it. second)
       )
      return false;
  return true;      
}




// NonTerminal

void NonTerminal::qc () const
{
  if (! qc_on)
    return;
  Symbol::qc ();

  QC_ASSERT (isAlpha (name [0]));

//QC_ASSERT (grammar);
  size_t leftRecursive = 0;
  size_t rightRecursive = 0;
  size_t nonRecursive = 0;
  Set<const Rule*> ruleSet;
  for (const Rule* r : lhsRules)
  {
    QC_ASSERT (r);
    r->qc ();
    QC_ASSERT (& r->lhs == this);
    ruleSet << r;
    if (r->isLeftRecursive ())
      leftRecursive++;
    else if (r->isRightRecursive ())
      rightRecursive++;
    else
      nonRecursive++;
  }
  QC_ASSERT (ruleSet. size () == lhsRules. size ());
  // Non-redundant
  QC_IMPLY (leftRecursive, nonRecursive);
  QC_IMPLY (rightRecursive, nonRecursive);

//QC_ASSERT (terminal2rules [false]. size () == terminal2rules [true] . size ());
  Set<const Rule*> all;
  for (const bool b : {false, true})
    for (const auto& it : terminal2rules [b])
    {
      const Terminal* t = it. first;
      QC_ASSERT (t);
    //QC_ASSERT (t->grammar == grammar);
      const VectorPtr<Rule>& rules = it. second;
      Set<const Rule*> s;
      s. insertAll (rules);
      QC_ASSERT (s. size () == rules. size ());
      for (const Rule* r : rules)
      {
        QC_ASSERT (r);
      //QC_ASSERT (r->lhs. grammar == grammar);
        QC_ASSERT (! r->erasable);
        QC_ASSERT (r->isLeftRecursive () == (bool) b);
      }
      all << s;
    }
  Set<const Rule*> s;
  s. insertAll (erasableRules);
  QC_ASSERT (s. size () == erasableRules. size ());
  for (const Rule* r : erasableRules)
  {
    QC_ASSERT (r);
  //QC_ASSERT (r->lhs. grammar == grammar);
    QC_ASSERT (r->erasable);
    QC_ASSERT (! r->isLeftRecursive ());  // Otherwise dummy left-recursivity
    all << r;
  }
  QC_ASSERT (all == ruleSet);

  QC_IMPLY (erasable, ! erasableRules. empty ());

  // Non-redundant
  QC_ASSERT (! lhsRules. empty ());  
  QC_IMPLY (lhsRules. size () == 1, ! lhsRules [0] -> empty ());
#if 0
  if (isRoot () != isTransient ())
  {
    cout << name << " " << (int) isRoot () << " " << (int) isTransient () << endl;
    ERROR;
  }
#endif
  QC_IMPLY (getEmptyRule (), lhsRules. size () > 1);
  QC_ASSERT (! firstROs. empty ());
  QC_ASSERT (! lastROs. empty ());
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
  for (const bool b : {false, true})
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
  for (const bool b : {false, true})
    for (const auto& it : terminal2rules [b])
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
  for (const auto& it : pos. terminal2syntagms)
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
  unique_ptr<Offset> ofs;
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
  for (const auto& it : pos. terminal2syntagms)
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
  if (! qc_on)
    return;
  Terminal::qc ();

  QC_ASSERT (grammar);
  QC_ASSERT (! endOfSentence (c));
  QC_ASSERT (! erasable);
}
#endif




#if 0
// SymbolTree

void SymbolTree::qc () const
{
  if (! qc_on)
    return;
  QC_IMPLY (children. empty (), ! rules. empty ());
  for (const auto& it : children)
  {
    QC_ASSERT (it. first);
    const SymbolTree* st = it. second;
    QC_ASSERT (st);
    st->qc ();
  }

  QC_ASSERT (! rules. contains (nullptr));
}



void SymbolTree::saveText (ostream &os) const
{
  for (const Rule* r : rules)
    os << ' ' << r->num;
  Offset ofs;
  for (const auto& it : children)
  {
    ofs. newLn (os);
    os << it. first->name << ':';
    it. second->saveText (os);
  }
}
#endif




// Grammar

Grammar::Grammar (const Grammar &other,
                  const VectorPtr<NonTerminal> &newTerminals)
// : Terminal (other. name, nullptr)
//, symbolTree (* new SymbolTree ())
{
  ASSERT (! newTerminals. empty ());
  other. qc ();

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
        StringVector rhs;
        for (const Symbol* child : r->rhs)
          rhs << child->name;
        addRule (r->lhs. name, rhs);
      }
    }

  finish (other. startSymbol->name);
}



Grammar::Grammar (const Grammar &other,
                  const Terminal* universalDelimiter)
// : Terminal (other. name, nullptr)
//, symbolTree (* new SymbolTree ())
{
  ASSERT (other. symbols. contains (universalDelimiter));
  other. qc ();

  for (const Symbol* s : other. symbols)
    if (const NonTerminal* nt = s->asNonTerminal ())
      for (const Rule* r : nt->lhsRules)
      {
        StringVector rhs;
        for (const Symbol* child : r->rhs)
          if (child != universalDelimiter)
            rhs << child->name;
        addRule (r->lhs. name, rhs);
      }

  finish (other. startSymbol->name);
}



void Grammar::init (CharInput &in)
{
  string startS;

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
          ASSERT (string (arrowS). size () == 2);
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



void Grammar::addRule (const string &lhs,
                       const StringVector &rhs)
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
  if (lhs == NonTerminal::sigmaS)
    throw runtime_error ("Symbol " + strQuote (NonTerminal::sigmaS) + " is reserved");
  
  string ascii ("   ");
  ASSERT (ascii. size () == 3);
  ascii [0] = '\"';
  ascii [2] = ascii [0];

  StringVector rhs;
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
        rhs << token. name; 
        rhsPart++;
        break;
      case Token::eInteger: 
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
                        addRule (optionalName, StringVector ());
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
                            case '+': tokens << move (Token (name)) 
                                             << move (Token (name))
                                             << move (Token ('*'));
                                      break;
                            case '*': addRule (nameMod, it, tokens. cend (), false);
                                      tokens << move (Token (name)) 
                                             << move (Token ('+'));
                                      break;
                            default: ERROR;
                          }
                          it = tokens. cbegin ();
                          addRule (nameMod, it, tokens. cend (), false);
                        }
                      }
                      break;
            default:  ERROR_MSG ("Unknown special grammar symbol: " + strQuote (token. name));
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
    StringVector rhs;  rhs. reserve (3);
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
  for (const bool out : {false, true})
    setNeighbors (out);
//setSymbolTree ();
  setReplacements ();
  setRegularStar ();
}



void Grammar::qc () const
{
  if (! qc_on)
    return;
//Terminal::qc ();

//QC_ASSERT (grammar != this);
  QC_ASSERT (! symbols. empty ());
  QC_ASSERT (startSymbol);
//QC_ASSERT (symbols [0] == startSymbol);

  // eotSymbol
  QC_ASSERT (eotSymbol);
//QC_ASSERT (eotSymbol->ruleOccurrences. empty ());

  Set<const Symbol*> symbolSet;
  Set<char> charSet;
  size_t terminals = 0;
  for (const Symbol* s : symbols)
  {
    QC_ASSERT (s);
  //QC_ASSERT (s->grammar == this);
    try { s->qc (); }
      catch (...) { cout << s->name << endl; throw; }
    symbolSet << s;
    if (const Terminal /*Letter*/ * let = s->asTerminal /*Letter*/ ())
    {
      terminals++;
      charSet << let->c;
    }
    // Non-redundant
    QC_IMPLY (/*s != eotSymbol &&*/ s->ruleOccurrences. empty (), s == startSymbol);  
  }
  QC_ASSERT (symbolSet. size () == symbols. size ());
  QC_ASSERT (charSet. size () == terminals);

  for (const auto& it : char2terminal)
    QC_ASSERT (symbols. contains (it. second));      

  // Valid CF-grammar
  if (const Symbol* s = getNonReachable ())
    ERROR_MSG (s->name + " is not reachable");
  if (const Symbol* s = getNonTerminable ())
    ERROR_MSG (s->name + " is not terminable");

#if 0
  for (const Rule* r : symbolTree. rules)
    QC_ASSERT (r->erasable);
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
  for (const auto& it : mainNexts)
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
  for (const auto& it : mainPrevs)
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
    const Rule::Occurrence prevRo = static_cast <const ROGraph::Node*> (node) -> ro;
    const Symbol* symbol = prevRo. getSymbol ();
    if (const Terminal* t = symbol->asTerminal ())
      const_cast <Terminal*> (t) -> terminalNeighbors [out] [prevRo] = Set<const Terminal*> ();
    for (const DiGraph::Arc* arc : node->arcs [out])
    {
      const Rule::Occurrence nextRo = static_cast <const ROGraph::Node*> (arc->node [out]) -> ro;
      const Symbol* neighbor = nextRo. getSymbol ();
      const_cast <Symbol*> (symbol) -> neighbors [out] << neighbor;
      if (const Terminal* t = symbol->asTerminal ())
        if (const Terminal* neighborT = neighbor->asTerminal ())
          const_cast <Terminal*> (t) -> terminalNeighbors [out] [prevRo] << neighborT;
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
    if (n->getDisjointCluster () != n)  // There exists a cycle which is not a self-loop
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



void Grammar::findAmbiguity () const /*throw (Ambiguity)*/
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
            {
              Set<string> tandemSet;
              const FollowROGraph frg (*this);
              frg. qc ();
              for (const DiGraph::Node* node : frg. nodes)
                if (static_cast <const ROGraph::Node*> (node) -> ro. getSymbol () == s)
                  for (const DiGraph::Arc* arc : node->arcs [true])
                    if (static_cast <const ROGraph::Node*> (arc->node [true]) -> ro. getSymbol () == nextSymbol)
                      tandemSet << static_cast <const ROGraph::Arc*> (arc) -> toString ();
              ASSERT (! tandemSet. empty ());
              List<string> tandems;
              insertAll (tandems, tandemSet);
              throw Ambiguity (intersection. front () -> name + "* tandem in rule(s):" + "\n" + list2str (tandems, "\n"));
            }
          }
  // tandem of the same erasable symbol ??
  // Other cases ??
}



void Grammar::prepare () /*throw (StdParserError, Ambiguity)*/
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
    const DiGraph::Node* root = const_cast <SymbolGraph::Node*> (dsg. getStartSymbol ()) -> getDisjointCluster ();
    size_t terminalNum = 0;
    for (DiGraph::Node* node1 : dsg. nodes)
      if (const Terminal* t = static_cast <SymbolGraph::Node*> (node1) -> symbol. asTerminal ())
        if (/*t != eotSymbol &&*/ node1->getDisjointCluster () != root)
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

  return res;
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
      for (const bool out : {false, true})
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
  for (const auto& it : scc2terminals)
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




// ROGraph::Arc

void ROGraph::Arc::qc () const
{
  if (! qc_on)
    return;
  DiGraph::Arc::qc ();
    
  QC_ASSERT (rhs_pos1 < rhs_pos2)
  QC_ASSERT (rhs_pos2 < rule. rhs. size ());
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
            for (const Rule::Occurrence prevRo : Rule::Occurrence (*r, i). getLastROs ())
            {
              Node* prev = ro2node [prevRo];
              ASSERT (prev);
              for (const Rule::Occurrence nextRo : Rule::Occurrence (*r, j). getFirstROs ())
              {
                Node* next = ro2node [nextRo];
                ASSERT (next);
                if (! prev->isIncident (next, true))
                  new Arc (prev, next, *r, i, j);
              }
            }
          }
      }
}



}  // namespace



/* TO DO: ??
  canEnd() for last rule occurrence
*/
