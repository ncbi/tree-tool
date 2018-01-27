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
    : Application ("Print the pairs of words from a list. In each pair: word1 < word2")
  	{
  	  addPositional ("list", "File with different words");
  	}



	void body () const final
	{
		const string listFName = getArg ("list");


    StringVector words;  
    {
      LineInput f (listFName);
      words = f. getVector ();
    }
    words. sort ();
    ASSERT (words. isUniq ());
    
    {
    	Progress prog ((uint) words. size ());
	    for (const string& w1 : words)
	    {
	    	prog (w1);
	    	for (const string& w2 : words)
	    		if (& w1 == & w2)
	    			break;
	    	  else
	    	  {
	    	  	ASSERT (w1 != w2);
	    	  	const string* p1 = & w1;
	    	  	const string* p2 = & w2;
	    	  	if (*p1 > *p2)
	    	  		swap (p1, p2);
	    	  	cout << *p1 << '\t' << *p2 << endl;
	    	  }
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
