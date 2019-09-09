// grammar_test.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "cfgrammar.hpp"
using namespace Cfgr;



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Test a CF-grammar")
    { addPositional ("grammar", "CF-Grammar file"); }


  void body () const final
  {
    const string grammarFName (getArg ("grammar"));
    ASSERT (! grammarFName. empty ());

    Grammar grammar (grammarFName);
    grammar. qc ();
    grammar. print (cout);
    cout << endl;
    
    grammar. prepare ();    


    const VectorPtr<NonTerminal> cuts (grammar. getCutSymbols ());
    cout << "Cuts:";
    for (const NonTerminal* nt : cuts)
      cout << ' ' << nt->name;
    cout << endl;

    cout << "Split terminals:";
    grammar. splitTerminals ();  // ??
    cout << endl;

    cout << endl;
    cout << "Terminals for SCC:" << endl;
    grammar. terminals4scc ();
    cout << endl;


    // PAR
  //const string s ("(R + R) * R - R / R");  
  //const string s ("123");  // PAR
    const string s ("(1 + 2) * 3e1");

    cout << s << endl;
    Sentence text (grammar, s);
    const Syntagms& syntagms = grammar. parseSentence (text);
    cout << "# Rollbacks: " << text. countRollbacks () << endl;
    cout << "# Wrong terminal syntagms: " << text. countWrongTerminalSyntagms () << endl;  // print ??
    cout << "# Wrong non-terminal syntagms: " << text. countWrongNonTerminalSyntagms () << endl;  // print ??
    if (syntagms. empty ())
      cout << "Not parsed" << endl;
    else
      for (const Syntagm* syntagm : syntagms)
        if (syntagm->right)
        {
          syntagm->print (cout);
          cout << endl;
        }

    cout << endl;
    cout << "Wrong syntagms:" << endl;
    for (const Position& pos : text. seq)
    {
      cout << pos. c << ":";
      for (const auto& it : pos. nonTerminal2syntagms)
        for (const Syntagm* syntagm : * it. second)
          if (! syntagm->right)
          {
            cout << "  ";
            const NonTerminalSyntagm* nts = syntagm->asNonTerminalSyntagm ();
            ASSERT (nts);
            cout << nts->size () << "  ";
            nts->rule.print (cout);
          }
      cout << endl;
    }

  #if 0
    ??
    const Grammar subgr (grammar, cuts);
    subgr. qc ();
    cout << endl;
    cout << "Subgrammar (cut): " << endl;
    subgr. print (cout);
    cout << endl;
  #endif

    if (const Terminal* ud = grammar /*subgr*/. getUniversalDelimiter ())
    {
      cout << "Universal delimiter: " << ud->name << endl;
      const Grammar subgrUd (grammar /*subgr*/, ud);
      subgrUd. qc ();
      cout << endl;
      cout << "Subgrammar (minus universal delimiter): " << endl;
      subgrUd. print (cout);
      cout << endl;
    }
  }
};



int main (int argc,
          const char* argv [])
{
  ThisApplication app;
  return app. run (argc, argv);
}
