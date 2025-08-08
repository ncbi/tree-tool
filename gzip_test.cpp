// gzip_test.cpp

#undef NDEBUG

#include "common.hpp"
#include "gzip.hpp"  
using namespace Common_sp;


#include "common.inc"



namespace
{
  
  
struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Test")
  	{
  	  addPositional ("arg", "Argument");
  	}
  


 	void body () const final
	{
    const string arg = getArg ("arg");
	  
  
  //size_t n = 0;
    {
      GZip f (arg, 100000);  // PAR
      while (f. nextLine ())
        cout << f. line << '\n';
      //n++;
    }
      
  //cout << n << endl;
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


