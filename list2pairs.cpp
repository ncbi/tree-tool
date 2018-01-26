// setMinus.cpp

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;


namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print the pairs of words from a list")
  	{
  	  addPositional ("list", "File with words");
  	}



	void body () const final
	{
		const string listFName = getArg ("list");


    StringVector words;  
    {
      LineInput f (listFName);
      words = f. getVector ();
    }
    
    {
    	Progress prog ((uint) words. size ());
	    for (const string& w1 : words)
	    {
	    	prog (w1);
	    	for (const string& w2 : words)
	    		if (& w1 == & w2)
	    			break;
	    	  else
	    	  	cout << w1 << '\t' << w2 << endl;
	    }
	  }
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}
