// key_test.cpp

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;



namespace
{
    
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("key_test")
  	{
  	  addPositional ("go", "Print key codes after pressing <Enter>. 'q' exits");
  	}
  


 	void body () const final
	{
    char c;
    do
    {
      getChar (cin, c);
      cout << (int) c << '\t';
      if (c <= ' ')
        cout << nonPrintable2str (c);
      else
        cout << c;
      cout << endl;
    } 
    while (c != 'q');   
  }
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


