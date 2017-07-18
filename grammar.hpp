// grammar.h

#ifndef GRAMMAR_HPP
#define GRAMMAR_HPP

#include "common.hpp"



namespace Grammar_sp
{
using namespace std;
using namespace Common_sp;


// Forward
struct Grammar;
struct Symbol;
  struct Terminal;
  struct NonTerminal;
struct Rule;
struct Syntagm;
  struct TerminalSyntagm;
  struct NonTerminalSyntagm;
struct Position;



typedef  unsigned int  Char;
  // Unicode
const Char eot = 0;
  // "End of text"
typedef  Set<Char>  Alphabet;
  // !contains(eot)
string chars2str (const Vector<Char> &vec);
  // Requires: Char in ASCII




// Syntagm

struct Sentence : Root 
// DAG of Syntagm's
// Ambiguous text
{
  Vector<Position> seq;
//const Alphabet* alphabet;  ??

  struct BadPos : runtime_error  
    { size_t pos;
      BadPos (size_t pos_arg, 
              const string &what_arg) 
        : runtime_error (what_arg + " at position " + toString (pos_arg + 1)) 
        , pos (pos_arg)
        {}
    };
  Sentence (const Grammar &grammar,
            const string &s) throw (BadPos);
    // After: grammar.setChar2terminal() 
    // Requires: Char in ASCII
  void qc () const;

  size_t countWrongTerminalSyntagms () const;
    // Input: Syntagm::right
  size_t countWrongNonTerminalSyntagms () const;
    // Input: Syntagm::right
  size_t countRollbacks () const;
    // Input: Position::rollbacks
  const NonTerminalSyntagm* getLastLongestSyntagm () const;
    // Return: may be nullptr
  void reportError (ostream &os) const;
    // Requires: Char in ASCII
  void printWrongSyntagms (ostream &os,
                           size_t size_min) const;
  VectorPtr<NonTerminalSyntagm> findSyntagms (const string& symbolName) const;
};



struct Syntagm : Root
// Symbol application
{
  // In the same Sentence
  const Position& begin;
  const Position& end;

  const Symbol& symbol;
    // In Grammar::symbols

  // Analysis
  bool right;

protected:
  Syntagm (const Position& begin_arg,
           const Position& end_arg,
           const Symbol &symbol_arg)
    : begin (begin_arg)
    , end (end_arg)
    , symbol (symbol_arg)
    , right (false)
    {}
public:
  void qc () const;
  void saveText (ostream &os) const;

  virtual const TerminalSyntagm* asTerminalSyntagm () const 
    { return nullptr; }
  virtual const NonTerminalSyntagm* asNonTerminalSyntagm () const 
    { return nullptr; }

  size_t size () const;
  virtual Vector<Char> str () const = 0;
  virtual void setRight () 
    { right = true; }
};



struct TerminalSyntagm : Syntagm
{
  TerminalSyntagm (const Position& begin_arg,
                   const Terminal& terminal);
private:
  void init (const Position& begin_arg,
             const Terminal& terminal);
public:
  void qc () const;

  const TerminalSyntagm* asTerminalSyntagm () const 
    { return this; }

  Vector<Char> str () const;
};



struct NonTerminalSyntagm : Syntagm
{
  const Rule& rule;
    // In Grammar
  VectorPtr<Syntagm> children;  
    // Tree
    // In Sentence
    // empty() <=> rule.lhs is erased

  NonTerminalSyntagm (const Position& begin_arg,
                      const Position& end_arg,
                      const Rule &rule_arg,
                      const VectorPtr<Syntagm> &children_arg);
  void qc () const;
  void saveText (ostream &os) const;

  const NonTerminalSyntagm* asNonTerminalSyntagm () const 
    { return this; }

  Vector<Char> str () const;
  void setRight ();
};



typedef  VectorOwn<Syntagm>  Syntagms;
  // !nullptr
  // All Syntagm's of a Symbol* starting at Position 
  // Different
  // (\exists s1, s2 \in nonTerminal2syntagms : s1 != s2 & s1->symbol = s2->symbol & s1.end = s2.end) <=> ambiguity



struct Position : Root
// In Sentence
{
  Char c;

  // Output
  // !nullptr
  map<const Terminal*, const Syntagms*> terminal2syntagms;  // of TerminalSyntagm*
    // !contains(t) <=> Syntagm's for t do not exist
  map<const NonTerminal*, const Syntagms*> nonTerminal2syntagms;  // of NonTerminalSyntagm*
    // !contains(nt) <=> nt has not been parse()'d
    // May be empty()
    // Memoization of NonTerminal::parse()

  size_t rollbacks;  // not maintaiend ??
    // --> incomplete Syntagm's in Syntagms ??


  Position ()
    : c (eot)
    , rollbacks (0)
    {}
 ~Position ()
    { for (auto it : terminal2syntagms)  
        delete it. second;
      for (auto it : nonTerminal2syntagms)  
        delete it. second;
    }
  TerminalSyntagm& init (const Grammar &grammar,
                         Char c_arg);
  void qc () const;

  const NonTerminalSyntagm* getLongestSyntagm () const;
    // Return: may be nullptr
};




// Grammar

struct Rule : Root
{
  typedef  VectorPtr<Symbol>  Rhs;
  typedef  Rhs::const_iterator  RhsIt;


  struct Occurrence;
  struct Occurrences : Set<Occurrence>  
  { Set<const Symbol*> getSymbols () const;
      // Return: !contains(nullptr)
  };


  struct Occurrence : Root
  // "Morphoneme"
  {
    const Rule& rule;
    RhsIt rhsIt;
      // In rule.rhs

    Occurrence (const Rule& rule_arg,
                RhsIt rhsIt_arg)
      : rule (rule_arg)
      , rhsIt (rhsIt_arg)
      {}
    Occurrence (const Rule& rule_arg,
                size_t index)
      : rule (rule_arg)
      , rhsIt (rule_arg. rhs. begin () + (int) index)
      {}
    bool operator== (const Occurrence &other) const
      { return    & rule == & other. rule
               && rhsIt == other. rhsIt;
      }
    void qc () const;
    void saveText (ostream &os) const;

    size_t getIndex () const;
    string getName () const;
    const Symbol* getSymbol () const
      { return *rhsIt; }
    bool operator< (const Occurrence &other) const;
    bool isFirst () const
      { return rhsIt == rule. rhs. begin (); }
    bool isLast () const
      { return rhsIt + 1 == rule. rhs. end (); }
    bool ruleCanStartWith () const;
    bool ruleCanEndWith () const;
    // Return: contains(*this)
    Occurrences getFirstROs () const;
    Occurrences getLastROs () const;
  };


  size_t num;
    // Unique within Rule's of the Grammar

private:
  friend struct Grammar;
  Vector<string> rhsS;  
    // Temporary
public:
  const NonTerminal& lhs;
    // In Grammar::symbols
  Rhs rhs;  
    // In lhs.grammar->symbols

  // Analysis
  // isLeftRecursive() => rhs[0] is skipped
  Set<const Terminal*> firstTerminals;  
    // In lhs.grammar->symbols
  bool erasable;
  Vector<bool> singleRhs;
    // size() = rhs.size()
    // true <=> all other symbols in rhs are erasable


  Rule (size_t num_arg,
        const NonTerminal& lhs_arg,
        const Vector<string> &rhsS_arg)
    : num (num_arg)
    , rhsS (rhsS_arg)
    , lhs (lhs_arg)
    , erasable (false)
    , singleRhs (rhsS_arg. size (), false)
    {}
private:
  void finish (Grammar &grammar);
    // Update: *grammar
public:
  void qc () const;
  void saveText (ostream &os) const;


  bool empty () const
    { return rhs. empty (); }
  bool isLeftRecursive () const;
  bool isRightRecursive () const;
  bool isTerminable () const;
  bool isErasable () const;
//bool isLeftRecursiveErasable () const;
    // Requires: isLeftRecursive()

private:
  void setFirstTerminalsErasable ();  
    // Input: Symbol::firstROs
  void setSingleRhs ();
    // Input: Symbol::erasable
public:
  void parseIt (RhsIt rhsIt,
                const Position &pos, 
                const VectorPtr<Syntagm> &children) const;
    // Invokes: new NonTerminalSyntagm(begin,pos,...)
};



struct Symbol : Named
{
  // Analysis
  Rule::Occurrences ruleOccurrences;
    // Occurrence::getSymbol() = this
  bool erasable;
  bool terminable;
  // In grammar->symbols
  Rule::Occurrences firstROs;  
  Rule::Occurrences lastROs;  
  // !nullptr
  Set<const Terminal*> terminals;
  Set<const Symbol*> neighbors [2/*bool next*/];
  Set<const Symbol*> replacements;
    // contains(x) <=> <x> \in D(<this>)


protected:
  explicit Symbol (const string &name_arg)
    : Named (name_arg)
    , erasable (false)
    , terminable (false)
    {}
public:    
  void qc () const;
  void saveText (ostream &os) const;


  virtual const Terminal* asTerminal () const
    { return nullptr; }
  virtual const NonTerminal* asNonTerminal () const
    { return nullptr; }
#if 0
  virtual const Letter* asLetter () const
    { return nullptr; }
  virtual const Grammar* asGrammar () const
    { return nullptr; }
#endif

  static string getName (const Symbol* s)   
    { return s ? s->name : "#"; }

  bool isRoot () const;
  Set<const Symbol*> getFirstSymbols () const
    { return move (firstROs. getSymbols () << this); }
  Set<const Symbol*> getLastSymbols () const
    { return move (lastROs. getSymbols () << this); }
private:
  friend struct Grammar;
  void setErasable ();
  void setTerminable ();
  void addFirstRO (Rule::Occurrence ro);
  void addLastRO (Rule::Occurrence ro);
  void addTerminal (const Terminal* t);
  Set<const Symbol*> getRulePrevs () const;
    // nullptr <=> *this is the first Symbol in rule.rhs
  Set<const Symbol*> getRuleNexts () const;
    // nullptr <=> *this is the last Symbol in rule.rhs
public:
  virtual const Syntagms* parse (const Position &pos) const = 0;
    // Non-deterministic (ambiguous) prefix parser
    // Return: pos::Syntagms
    // Update: pos.nonTerminal2syntagms(this) - all Syntagm's of *this starting at pos
};



typedef  map <Rule::Occurrence, Set<const Terminal*> >  Occurrence2Terminals;



struct Terminal : Symbol
// name: '<ASCII>' or <Unicode>
{
  static const string eotName;
  Char c;

  // Analysis
  Occurrence2Terminals terminalNeighbors [2/*bool next*/];  
    // keys() = ruleOccurrences 


  Terminal ()
    : Symbol (eotName)  
    , c (eot)
    {}
    // Initial or final symbol of a sentence
  explicit Terminal (const string &name_arg);
  void qc () const;


  const Terminal* asTerminal () const
    { return this; }
  static const Terminal* as (const Symbol* s) 
    { return s->asTerminal (); }

  bool differentiated (Rule::Occurrence ro,
                       bool out) const;
    // Input: ro in terminalNeighbors[out].keys()

  const Syntagms* parse (const Position &pos) const
    { return findPtr (pos. terminal2syntagms, this); }
    // Return: TerminalSyntagm*'s
};



struct NonTerminal : Symbol
// name: no quotes
{
  static const string sigmaS;
  VectorOwn<Rule> lhsRules;
    // *this == Rule::lhs 

  // Analysis
  // Cover of lhsRules 
  map <const Terminal*, VectorPtr<Rule> > terminal2rules [2/*Rule::isLeftRecursive()*/];
    // Parsing table
    // Index by Terminal::num ??
    // !Rule::erasable
  VectorPtr<Rule> erasableRules;
    // Rule::erasable

  // !nullptr
  Set<const Symbol*> regularStar;
    // contains(x) <=> *this = x* or x+ (regular expression) depending on erasable


  NonTerminal ()
    : Symbol (sigmaS)
    {}
  explicit NonTerminal (const string &name_arg/*,
                        const Grammar &grammar_arg*/)
    : Symbol (name_arg /*, & grammar_arg*/)
    {}
  void qc () const;
  void saveText (ostream &os) const;


  const NonTerminal* asNonTerminal () const
    { return this; }
  static const NonTerminal* as (const Symbol* s) 
    { return s->asNonTerminal (); }

  const Rule* getEmptyRule () const;
    // Return: some empty rule
  bool isTransient () const
    { return    lhsRules.        size () == 1 
             && ruleOccurrences. size () == 1; 
    }


private:
  friend struct Grammar;
  void setTerminal2rules ();
  void setRegularStar ();
    // Input: replacements
  double getComplexity () const;
    // Return: >= 0
public:
  bool canEnd (const Position &pos) const;
    // Requires: pos != Sentence::seq.end()
  const Syntagms* parse (const Position &pos) const;
    // Return: !nullptr, VectorOwn* of NonTerminalSyntagm*
};


  
#if 0
struct Letter : Terminal
{
  static const string space;
  Char c;


  Letter (const string &name_arg,
          const Grammar &grammar_arg)
    : Terminal (name_arg, & grammar_arg)
    , c (getChar (name_arg))
    {}
private:
  static Char getChar (const string &name) 
    { return name == space
               ? ' ' 
               : name. size () == 1
                 ? name [0]
                 : throw runtime_error ("Bad Terminal " + name);
    }
public:
  void qc () const;

  const Letter* asLetter () const
    { return this; }
  static const Letter* as (const Symbol* s) 
    { return s->asLetter (); }
};
#endif



#if 0
struct SymbolTree : Root
// Rule::Rhs: left to right
{
  map <const Symbol*, const SymbolTree* /*!nullptr*/> children;
  VectorPtr<Rule> rules;
    // children.empty() => !rules.empty()

 ~SymbolTree ()
    { for (auto it : children)
        delete it. second;
    }
  void qc () const;
  void saveText (ostream &os) const;
};
#endif



struct Grammar : Root
// grammar  ??
  // Tree
  // nullptr <=> *this is a root grammar
  // !nullptr <=> *this is a symbol grammar
{
  static const string arrowS;
  static const char commentC = '#';
private:
  friend struct Rule;
  size_t ruleNum;
  map<string,Symbol*> string2Symbol;
public:

  VectorOwn<Symbol> symbols;  
  // In symbols
  const NonTerminal* startSymbol;
  const Terminal* eotSymbol;

#if 0
  const SymbolTree& symbolTree;
    // rules = erasable rules
#endif

  map <Char, const Terminal*> char2terminal; 
    // size() == charSize
    // symbols.contains(Terminal*)


  explicit Grammar (const string &fName);
    // Normal CF-grammar
    // Input: fName: sequence of rules 
    //   <rule> ::= <NonTerminal> <arrowS> <Symbol>*
    //     space-delimited
    //     starts from the first position of a line
    //   '#' starts a comment which ends at the line end
    //   Special characters: [ ] * +
    //   Automatically added rules for symbols A* or A+ in r.h.s.:
    //     A+ -> A A*
    //     A* ->
    //     A* -> A+
    //   <Terminal symbol> ::= <Unicode> | "<sequence of printable ASCII symbols>"
    // Output: startSymbol = the <NonTerminal> of the first rule
  // Subgrammars ??
//explicit Grammar (const NonTerminal &nt);
  Grammar (const Grammar &other,
           const VectorPtr<NonTerminal> &newTerminals);
    // Input: newTerminals: in other, !empty()
  Grammar (const Grammar &other,
           const Terminal* universalDelimiter);
    // Input: other.symbols.contains(universalDelimiter)
private:
  template <typename T>
    T* getSymbol (const string &lhs)
      // Return: !nullptr
      { Symbol* s = string2Symbol [lhs];
        if (! s)
        { s = new T (lhs /*, *this*/);
          symbols << s;
          string2Symbol [lhs] = s;
        };
        const T* res = T::as (s);
        ASSERT (res);
        return const_cast <T*> (res);
      }
  void addRule (const string &lhs,
                const Vector<string> &rhs);
    // Update: ruleNum
  void addRule (const string &lhs,
                List<Token>::const_iterator &rhsIt,
                const List<Token>::const_iterator &rhsEnd,
                bool optional);
    // Update: rhsIt
    // Invokes: addRule(lhs,rhs)
  void finish (const string &startS);
    // Invokes: Rule::finish(), see below
public:
  void qc () const;
    // Invokes: see below
  void saveText (ostream &os) const;
    // Invokes: getComplexity()


//const Grammar* asGrammar () const
  //{ return this; }

  // Invoked in finish()
  void setTerminable ();
  void setErasable ();
  void setRuleSingleRhs ();
  void setFirstROs ();
    // Input: Symbol::erasable
  void setLastROs ();
    // Input: Symbol::erasable
  void setParsingTable ();
    // Input: Symbol::firstROs
  void setTerminals ();
  void setNeighbors (bool out);
    // Output: Symbol::neighbors[], Terminal::terminalNeighbors[]
//void setSymbolTree ();
  void setReplacements ();
    // Output: Symbol::replacements
  void setRegularStar ();
    // Input:  Symbol::replacements
    // Output: NonTerminal::regularStar

  // Invoked in qc():
  const Symbol* getNonReachable () const;
    // Return: nullptr <=> all symbols are reachable
  const Symbol* getNonTerminable () const;
    // Return: !terminable; nullptr <=> all symbols are terminable

  void setChar2terminal ();
    // Update: char2terminal
  const Symbol* findParseCycle () const;
    // Return: Symbol in a parse() cycle; may be nullptr 
  const Rule* getLeftRecursiveErasable () const;
    // Return: may be nullptr
  struct Ambiguity : runtime_error  { Ambiguity (const string &what_arg) : runtime_error (what_arg) {} };
  void findAmbiguity ()  const throw (Ambiguity);
    // Not important for a deterministic parser
  struct StdParserError : runtime_error  { StdParserError (const string &what_arg) : runtime_error (what_arg) {} };
  void prepare () throw (StdParserError, Ambiguity);
    // Invokes: findParseCycle(), getLeftRecursiveErasable(), setChar2terminal(), findAmbiguity()
    
  double getComplexity () const;
    // Return: >= 0; 0 => Time_parse(s) = O(s.size()), unambiguous grammar

  const Syntagms& parseSentence (Sentence &sentence) const;
    // Top-down parser
    // Return: may be empty() 
    // Update: sentence.seq::Position
    // Invokes: NonTerminal::parse(), Syntagm::setRight()
    // Requires: after prepare() 

  // Analysis
  VectorPtr<NonTerminal> getCutSymbols () const;
  const Terminal* getUniversalDelimiter () const;
    // Return: may be nullptr
  // Invoke: only print() ??
  void splitTerminals () const;  
  void terminals4scc () const;
};




// Analysis

/*
relation:
  on:
    symbols * symbols
    rule occurrences * rule occurrences
    rule occurrences * symbols
    terminal symbols * rules
  property:
    derivation
    first/last symbol derivation
    single symbol derivation
    prev/next symbol
*/



struct SymbolGraph : DiGraph
// nodes: 1-1 with grammar->symbols
// arcs: not parallel
{
  struct Node : DiGraph::Node
  {
    const Symbol& symbol;
    bool reachable;

    Node (SymbolGraph &graph_arg,
          const Symbol &symbol_arg)
      : DiGraph::Node (graph_arg)
      , symbol (symbol_arg)
      , reachable (false)
      {}

    string getName () const
      { return symbol. name; }
    void setReachable ();
      // From parent to children
      // Output: this->reachable = true
    const Terminal* getNextTerminal (bool out) const;
      // Return: !nullptr <=> there is one next Terminal
  };


  const Grammar& grammar;
  map <const Symbol*, Node*> symbol2node;

protected:
  explicit SymbolGraph (const Grammar &grammar_arg)
    : grammar (grammar_arg)
    { // nodes
      for (const Symbol* s : grammar. symbols)
        symbol2node [s] = new Node (*this, *s);
    }
public:

  const Node* getStartSymbol () const
    { return symbol2node. find (grammar. startSymbol) -> second; }
    // Return: !nullptr
  void clearReachable ()
    { for (const DiGraph::Node* node : nodes)
        const_static_cast <SymbolGraph::Node*> (node) -> reachable = false;
    }
  const Node* getNonReachable () const; 
    // Return: nullptr <=> not exist
};



struct DerivedSymbolGraph : SymbolGraph
{
  explicit DerivedSymbolGraph (const Grammar &grammar_arg)
    : SymbolGraph (grammar_arg)
    { initArcs (); }
private:
  void initArcs ();
    // Arcs: {<x,lhs(r)> | x \in rhs(r), r \in R}
};



struct FirstDerivedSymbolGraph : SymbolGraph
{
  explicit FirstDerivedSymbolGraph (const Grammar &grammar_arg)
    : SymbolGraph (grammar_arg)
    { initArcs (); }
private:
  void initArcs ();
};



struct SingleDerivedSymbolGraph : SymbolGraph
{
  explicit SingleDerivedSymbolGraph (const Grammar &grammar_arg)
    : SymbolGraph (grammar_arg)
    { initArcs (); }
private:
  void initArcs ();
    // Arcs: {<x,lhs(r)> | r \in R, rhs(r) = \alpha x \beta, \alpha is erasable, \beta is erasable}
    // Input: Symbol::erasable
};



struct FollowSymbolGraph : SymbolGraph
{
  explicit FollowSymbolGraph (const Grammar &grammar_arg)
    : SymbolGraph (grammar_arg)
    { initArcs (); }
private:
  void initArcs ();
};



struct ROGraph : DiGraph
// nodes: 1-1 with Rule::Occurrence's
// arcs: not parallel
{
  struct Node : DiGraph::Node
  {
    Rule::Occurrence ro;

    Node (ROGraph &graph_arg,
          Rule::Occurrence ro_arg)
      : DiGraph::Node (graph_arg)
      , ro (ro_arg)
      {}

    string getName () const
      { return ro. getName (); }
    Set<const Terminal*> getNextTerminals (bool out) const;
  };


  struct Arc : DiGraph::Arc  
  {
    const Rule& rule;
    size_t rhs_pos1;
    size_t rhs_pos2;
      // rhs_pos1 < rhs_pos2 < rule.rhs.size()
    
    Arc (Node* start,
         Node* end,
         const Rule& rule_arg,
         size_t rhs_pos1_arg,
         size_t rhs_pos2_arg)
      : DiGraph::Arc (start, end)
      , rule (rule_arg)
      , rhs_pos1 (rhs_pos1_arg)
      , rhs_pos2 (rhs_pos2_arg)
      {}
    void qc () const;  
    
    string toString () const
      { return rule. str () + " at positions " +
                         Common_sp::toString (rhs_pos1 + 1)
                 + "-" + Common_sp::toString (rhs_pos2 + 1);
      }
  };


  const Grammar& grammar;
  map <Rule::Occurrence, Node*> ro2node;


protected:
  explicit ROGraph (const Grammar &grammar_arg)
    : grammar (grammar_arg)
    { // nodes
      for (const Symbol* s : grammar. symbols)
        for (const Rule::Occurrence ro : s->ruleOccurrences)
          ro2node [ro] = new Node (*this, ro);
    }
};



struct FollowROGraph : ROGraph
{
  explicit FollowROGraph (const Grammar &grammar_arg)
    : ROGraph (grammar_arg)
    { initArcs (); }
private:
  void initArcs ();
};



}  // namespace



#endif



/* TO DO ??
  first/last --> bool 
*/
