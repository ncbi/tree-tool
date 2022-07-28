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


#ifdef _MSC_VER
  #error "UNIX is required"
#endif

#include "common.hpp"
using namespace Common_sp;


// ftp://ftp.gnu.org/gnu/ncurses/ncurses-6.2.tar.gz
// Old: ftp://ftp.gnu.org/pub/gnu/ncurses/ncurses.tar.gz
extern "C"
{
  #include <ncurses.h>
}



namespace Common_sp
{
  
  

struct NCurses : Singleton<NCurses>
{
  bool hasColors {false};
  enum Color {colorNone, colorRed, colorGreen, colorYellow, colorBlue, colorMagenta, colorCyan};
  chtype background {0};
  // Change after screen resizing
  size_t row_max {0};
  size_t col_max {0};
  
  explicit NCurses (bool hideCursor);
 ~NCurses ();

  void resize ();
};



struct NCAttr : Root
{
  const int attr;
  const bool active;

  explicit NCAttr (int attr_arg,
                   bool active_arg = true);
 ~NCAttr ();
};



struct NCAttrColor : NCAttr
{
  explicit NCAttrColor (NCurses::Color color,
                        bool active_arg = true);
};



struct NCBackground 
{
  const chtype background_old;
  
  explicit NCBackground (chtype background);
 ~NCBackground ();
};




}  // namespace


