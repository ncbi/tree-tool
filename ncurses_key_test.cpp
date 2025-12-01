// ncurses_key_test.cpp

/*===========================================================================
*
*                            PUBLIC DOMAIN NOTICE                          
*               National Center for Biotechnology Information
*                                                                          
*  This software/database is a "United States Government Work" under the   
*  terms of the United States Copyright Act.  It was written as part of    
*  the author's official duties as a United States Government employee and 
*  thus cannot be copyrighted.  This software/database is freely available 
*  to the public for use. The National Library of Medicine and the U.S.    
*  Government have not placed any restriction on its use or reproduction.  
*                                                                          
*  Although all reasonable efforts have been taken to ensure the accuracy  
*  and reliability of the software and data, the NLM and the U.S.          
*  Government do not and cannot warrant the performance or results that    
*  may be obtained by using this software or data. The NLM and the U.S.    
*  Government disclaim all warranties, express or implied, including       
*  warranties of performance, merchantability or fitness for any particular
*  purpose.                                                                
*                                                                          
*  Please cite the author in any work or product based on this material.   
*
* ===========================================================================
*
* Author: Vyacheslav Brover
*
* File Description:
*   NCurses key test
*
*/


#undef NDEBUG

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

#include "common.inc"



namespace
{
    
  
struct ThisApplication final : Application
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
    ::initscr (); 
    ::cbreak ();
    ::noecho ();
    ::keypad (::stdscr, TRUE);
  #endif

    int key = 0;
    do
    {
      key = 
        #ifdef USE_NC
           NCurses_sp::getKey ()
        #else
           ::getch ()
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
      cout << '\t' << keyname (key) << "\r\n";
    } 
    while (key != 'q');   

  #ifndef USE_NC   
    ::endwin ();
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


