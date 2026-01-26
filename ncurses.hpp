// ncurses.hpp

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
*   ncurses classes
*
*/


#ifndef NCURSES_HPP_59734  // random number
#define NCURSES_HPP_59734


#ifdef _MSC_VER
  #error "UNIX is required"
#endif

#include "common.hpp"
using namespace Common_sp;

// ftp://ftp.gnu.org/gnu/ncurses/ncurses-6.2.tar.gz
extern "C"
{
  #include <ncursesw/curses.h>
    // Linking requires: -lncursesw 
    // Test: test/ncurses
    // $TERM info: infocmp
    // man {3|5} terminfo
    // vi does not use $TERM ??
}




namespace NCurses_sp
{
  
 
int getKey ();
  // Invokes: ::getch() (=> ::refresh())

constexpr int ctrl (int x) 
  { return (x) & 0x1f; }



struct NCurses final : Singleton<NCurses>
{
  bool hasColors {false};
  enum Color {white, red, green, yellow, blue, magenta, cyan, colorExtra};
    // <color>-black ncurses color-pairs
  static ::attr_t color2attr (Color color)
    { return COLOR_PAIR (color + 1); }  
  ::chtype background {0};
    
  size_t row_max {0};
  size_t col_max {0};
  
  
  explicit NCurses (bool hideCursor);
    // Invokes: setlocale(LC_ALL,"en_US.UTF-8"), resize()
 ~NCurses ()
    { ::endwin (); }


  static void setColorPair (Color color,
                            short foreColor,
                            short backColor);
  void resize ();
    // Update: row_max, col_max
    // Invokes: ::erase()
  bool print (size_t y,
              size_t x,
              const string &s) const
    { if (col_max <= x)
        return false;
      if (row_max <= y)
        return false;
      ::move ((int) y, (int) x);  
      ::addstr (pad (s, col_max - x, efalse). c_str ());
      return true;
    }
};



struct Attr : Root
// { Attr a (x); Attr b (y); } = Attr c (x | y);
{
  const ::attr_t attr;
  const bool active;


  explicit Attr (::attr_t attr_arg,
                 bool active_arg = true)
		: attr (attr_arg)
		, active (active_arg)
    { if (active)
        ::attron (attr); 
    }
    // A_UNDERLINE can conflict with colors, see infocmp ncv
  ~Attr ()
    { if (active)
        ::attroff (attr); 
    }
};



struct AttrColor final : Attr
{
  explicit AttrColor (NCurses::Color color,
                      bool active_arg = true)
    : Attr (NCurses::color2attr (color), active_arg)
    {}
};



struct Background 
// Used in ::wclear(), ::clrtoeol()
{
  const ::chtype background_old;

  
  explicit Background (::chtype background)
		: background_old (getbkgd (stdscr))
		{ ::bkgdset (background); }
 ~Background ()
    { ::bkgdset (background_old); }
};



struct Window
{
//NCurses &nc;  ??
	// Global
	const size_t global_x;
  const size_t global_y;
  //
  const size_t width;
  const size_t height;
	::WINDOW* win;
	
	
	Window (size_t global_x_arg,
	        size_t global_y_arg,
	        size_t width_arg,
	        size_t height_arg);
	Window (const NCurses &nc,
	        size_t width_arg,
	        size_t height_arg)
	  : Window ( nc. col_max - width_arg  - 1
	  	       , 0 
	  	       , width_arg
	  	       , height_arg
	  	       )
	  {}    
 ~Window ()
    {	//::wborder (win, ' ', ' ', ' ',' ',' ',' ',' ',' ');
    	::werase (win);
    	::wrefresh (win);
			::delwin (win);
    }
    
    
  void print (size_t x,
              size_t y,
              const string &s) 
    { if (s. empty ())
    	  return;
    	::mvwprintw (win, (int) y, (int) x, "%s", s. c_str ());
      ::wrefresh (win);
    }
  void cursor (size_t x,
               size_t y)
    { ::wmove (win, (int) y, int (x));
    	::wrefresh (win);
    }
};



inline void message (const string &text)
	{	Window w (5, 5, text. size () + 4, 3);  // PAR
	  w. print (2, 1, text);
	  getKey ();
	}




struct Field
{
private:
	Window &win;
  const size_t x;
  const size_t y;
	const size_t width;
	  // > 1
	const bool upper;
	  // For ASCII
	size_t val_start {0};
	  // < val.size()
	size_t pos {0};
	  // < width
public:
	StringVector val;  
	  // UTF-8
	
	
	Field (Window &win_arg,
	       size_t x_arg,
         size_t y_arg,
	       size_t width_arg,
	       bool upper_arg,
	       const StringVector &val_arg);
	Field (Window &win_arg,
	       size_t x_arg,
         size_t y_arg,
	       size_t width_arg,
	       bool upper_arg,
	       const string &val_arg)
	  : Field (win_arg, x_arg, y_arg, width_arg, upper_arg, StringVector {{val_arg}})
	  {}
	  
	  
	void print () const;
	enum Exit {fieldDone, fieldCancel, fieldNext, fieldPrev};
	Exit run ();
	  // Update: val, val_start, pos
	string getS () const
	  { return val. toString (); }
};



struct Form : Window
{
  VectorOwn<Field> fields;
    // !empty()
	
	
	Form (const NCurses &nc,
	      size_t width_arg,
	      size_t height_arg)
	  : Window (nc, width_arg, height_arg)
	  {}
	void add (size_t x,
		        size_t y,
			      size_t width_arg,
			      bool upper,
			      const StringVector &val)
	  { fields << new Field (*this, x, y, width_arg, upper, val); }
	void add (size_t x,
		        size_t y,
			      size_t width_arg,
			      bool upper,
			      const string &val)
	  { fields << new Field (*this, x, y, width_arg, upper, val); }


	bool run ();
	  // Return: false <=> cancelled
	StringVector getVal (size_t i) const
	  { return fields [i] -> val; }
	string getS (size_t i) const
	  { return fields [i] -> getS (); }
};




}  // namespace



#endif
