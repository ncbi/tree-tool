// ncurses.cpp

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


#undef NDEBUG
#include "common.inc"

#include "ncurses.hpp"
using namespace Common_sp;



namespace Common_sp
{
  

// NCurses

NCurses::NCurses (bool hideCursor)
{ 
  initscr (); 
  cbreak ();
  noecho ();
  keypad (stdscr, TRUE);
  if (hideCursor)
    curs_set (0);  
  hasColors = has_colors ();
  if (hasColors)
    { EXEC_ASSERT (start_color () == OK); }
  resize ();
  constexpr short bkgdColor = COLOR_BLACK;
  init_pair (1, COLOR_WHITE,   bkgdColor);  // colorNone
  init_pair (2, COLOR_RED,     bkgdColor);
  init_pair (3, COLOR_GREEN,   bkgdColor);
  init_pair (4, COLOR_YELLOW,  bkgdColor);
  init_pair (5, COLOR_BLUE,    bkgdColor);
  init_pair (6, COLOR_MAGENTA, bkgdColor);
  init_pair (7, COLOR_CYAN,    bkgdColor);
  init_pair (8, COLOR_WHITE,   bkgdColor);  
  background = COLOR_PAIR (1);
  bkgdset (background);
  attron (COLOR_PAIR (1)); 
  wclear (stdscr);
}



NCurses::~NCurses ()
{ 
  endwin (); 
}



void NCurses::resize ()
{ 
  int row_max_, col_max_;
  getmaxyx (stdscr, row_max_, col_max_); 
  QC_ASSERT (row_max_ >= 0);
  QC_ASSERT (col_max_ >= 0);
  row_max = (size_t) row_max_;
  col_max = (size_t) col_max_;
}




// NCAttr


NCAttr::NCAttr (attr_t attr_arg,
                bool active_arg)
: attr (attr_arg)
, active (active_arg)
{ 
  if (active)
    attron (attr); 
}



NCAttr::~NCAttr ()
{ 
  if (active)
    attroff (attr); 
}




// NCAttrColor

NCAttrColor::NCAttrColor (NCurses::Color color,
                          bool active_arg)
: NCAttr ((attr_t) COLOR_PAIR (color + 1), active_arg)
{}




// NCBackground 

NCBackground ::NCBackground (chtype background)
: background_old (getbkgd (stdscr))
{ 
  bkgdset (background); 
}



NCBackground ::~NCBackground ()
{ 
  bkgdset (background_old); 
}




}  // namespace


