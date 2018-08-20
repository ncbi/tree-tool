// NL_test.cpp

#undef NDEBUG
#include "..\common.inc"

#include "..\common.hpp"
using namespace Common_sp;
#include "language.hpp"
using namespace Lang;



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("NL test")
    { addPositional ("grammar", "NL grammar file"); }



  void body () const final
  {
    const string grammarFName (getArg ("grammar"));
    ASSERT (! grammarFName. empty ());


    Language lang (grammarFName);
    lang. qc ();
  }
};



int main (int argc,
          const char* argv [])
{
  ThisApplication app;
  return app. run (argc, argv);
}
