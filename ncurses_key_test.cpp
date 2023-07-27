// ncurses_key_test.cpp

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;


#define USE_NC  // PAR

#ifdef USE_NC
  #include "ncurses.hpp"
  using namespace NCurses_sp;
#else
  extern "C"
  {
    #include <ncursesw/curses.h>
  }
#endif


namespace
{
    
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("NCurses key test")
  	{
  	  addPositional ("go", "Print key codes. 'q' exits");
  	}
  


 	void body () const final
	{	  
	#ifdef USE_NC
	  NCurses nc (false);
	#else
    initscr (); 
    cbreak ();
    noecho ();
    keypad (stdscr, TRUE);
  #endif

    int key = 0;
    do
    {
      key = 
        #ifdef USE_NC
           NCurses_sp::getKey ()
        #else
           getch ()
        #endif
        ;
      cout << key << '\t';
      if (key > 0 && key < 127)
      {
        if (key <= ' ')
          cout << nonPrintable2str ((char) key);
        else
          cout << key;
      }
      cout << '\t' << keyname (key) << "\n\r";
    } 
    while (key != 'q');   

  #ifndef USE_NC   
    endwin ();
  #endif
  }
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


