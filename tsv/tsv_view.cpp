// tsv_view.cpp

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
*   tsv-table viewer
*
*/


#undef NDEBUG

#include "../common.hpp"
#include "tsv.hpp"
using namespace Common_sp;
#include "../ncurses.hpp"
using namespace NCurses_sp;
#include "../version.inc"

#include "../common.inc"



#define NUM_P  // ??



namespace
{
  
  
bool consistent = false;
size_t width_max = 0;
  
  
  
bool printString (string s,
                  size_t screen_col_max,
                  size_t &x)
{
  bool truncated = false;
  if (s. size () > screen_col_max - x)  // UTF-8 ??
  {
    s. erase (screen_col_max - x);
    truncated = true;
  }
  printw ("%s", s. c_str ());
  x += s. size ();
  return ! truncated;
}
  
  
  
void printRow (bool is_header,
               const StringVector &values,
               const Vector<bool> &active,
               size_t col_start,
               const Vector<TextTable::Header> &header,
               size_t screen_col_max,
               size_t rows_max)
{
  ASSERT (header. size () == active. size ());
  ASSERT (col_start + values. size () == header. size ());

  size_t x = 0;
  FFOR_START (size_t, col, col_start, header. size ())
  {
    if (! active [col])
      continue;
    ASSERT (x <= screen_col_max);
    const TextTable::Header& h = header [col];
    string value (values [col - col_start]);
    trim (value);
    if (   consistent
        && ! is_header 
        && ! strNull (value)
        && h. numeric 
        && ! h. scientific
        && h. decimals
       )
    {
      bool hasPoint = false;
      streamsize decimals = 0;
      getScientific (value, hasPoint, decimals);
      ASSERT (h. decimals >= decimals);
      // Data modification ??
      if (! hasPoint)
        value += ".";
      value += string ((size_t) (h. decimals - decimals), '0');  
    }
    if (value. size () > h. len_max)
    {
      value. erase (h. len_max - 3);
      value += "...";
    }
    const ebool right = (   rows_max > TextTable::Header::choices_max * 2  // PAR
                         && h. choices. size () < (h. numeric && ! h. scientific ? 3 : TextTable::Header::choices_max)  // PAR
                        )
                          ? enull
                          : h. numeric && ! h. scientific
                            ? etrue
                            : efalse;
    if (! printString (pad (value, h. len_max, right), screen_col_max, x))
      break;
    if (col + 1 < header. size ())
      if (! printString ("  ", screen_col_max, x))
        break;
  }

  FFOR_START (size_t, i, x, screen_col_max)
    printw (" ");
  //clrtoeol ();  // may start erasing before the cursor
}



bool moveRight (size_t &curCol,
                const Vector<bool> &active)
// Return: success
{
  ASSERT (! active. empty ());
  ASSERT (curCol < active. size ());
  
  const size_t curCol_old = curCol;
  curCol++;
  while (curCol < active. size () && ! active [curCol])
    curCol++;
  if (curCol < active. size ())
    return true;
    
  curCol = curCol_old;
  return false;
}
                


bool moveLeft (size_t &curCol,
               const Vector<bool> &active)
// Return: success
{
  ASSERT (! active. empty ());
  
  const size_t curCol_old = curCol;
  do 
  { 
    if (! curCol)
    {
      curCol = curCol_old;
      return false;
    }
    curCol--; 
  }
  while (! active [curCol]);

  return true;
}
  
  
  
struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("View a tsv-table")
  	{
      version = VERSION;
  	  addPositional ("table", "tsv-table");
  	  addFlag ("consistent", "Reformat numbers to make them consistent, e.g., to have the same number of decimals");
  	  addKey ("width", "Max. column width", "30");
  	}



	void body () const final
	{
		const string tableFName = getArg ("table");
		             consistent = getFlag ("consistent");
		             width_max  = str2<size_t> (getArg ("width"));

    if (width_max < 4)
      throw runtime_error ("Too small -width: " + to_string (width_max));
   
    TextTable tt (tableFName);
    tt. qc ();
    if (tt. rows. empty ())
    {
      cout << "No data" << endl;
      return;
    }

    FFOR (size_t, i, tt. header. size ())
    {
      TextTable::Header& h = tt. header [i];      
      minimize (h. len_max, width_max);
      maximize (h. len_max, max (h. name. size (), to_string (i + 1). size ()));
    }
    

    size_t topIndex = 0;
    size_t curIndex = topIndex;
    size_t curCol = 0;
  //size_t curLastCol = curCol;
    string what;  // For search
    Vector<bool> active (tt. header. size (), true);
  #ifdef NUM_P
    constexpr bool numP = true;
  #else
    bool numP = false;
  #endif
    NCurses nc (true);
    bool quit = false;
    Vector<bool> rowFound;
    rowFound. resize (tt. rows. size (), false);
    while (! quit)
    {
      nc. resize ();
      const size_t headerSize = 2 /*file name, header*/ + (size_t) numP; 
      const size_t fieldSize = nc. row_max - (headerSize + 1 /*menu row*/); 
      const size_t pageScroll = fieldSize - 1;
      const size_t bottomIndex_max = topIndex + fieldSize;
      const size_t bottomIndex = min (tt. rows. size (), bottomIndex_max);
      ASSERT (topIndex <= curIndex);
      ASSERT (curIndex < bottomIndex);
      if (   nc. row_max > headerSize + 2
          && nc. col_max > 10  // PAR
         )
      {
        {
          ::move (0, 0);
          const Attr attr (A_BOLD);
          addstr (tableFName. c_str ());
          ::clrtoeol ();
        }
        {
          ::move (1, 0);
          const Attr attr (A_BOLD);
          const Background bkgr (COLOR_PAIR (7) /*nc. background*/ | A_BOLD);
          StringVector values;
          FFOR_START (size_t, j, curCol, tt. header. size ())
            values << tt. header [j]. name;
          /*curLastCol =*/ printRow (true, values, active, curCol, tt. header, nc. col_max, tt. rows. size ());
        }        
        if (numP)
        {
          ::move (2, 0);
          const Attr attr (COLOR_PAIR (4) /*| A_BOLD*/);
          StringVector values;
          FFOR_START (size_t, j, curCol, tt. header. size ())
            values << to_string (j + 1);
          printRow (true, values, active, curCol, tt. header, nc. col_max, tt. rows. size ()); 
        }
        ::move ((int) (fieldSize + headerSize), 0);
        {
          const Attr attr (A_BOLD);
          const Background bkgr (COLOR_PAIR (3) /*nc. background*/ | A_BOLD);
          const string keyS ("Up  Down  Left  Right  PgUp  PgDn  Home  End  F3:Search from cursor  m:(un)mark row  a:(un)mark all rows  -:remove column  +:restore all columns"
                           #ifndef NUM_P
                             "  #:numbers"
                           #endif
                             "  q:Quit"
                            );
            // Non-character keys may be intercepted by the terminal
        #ifdef NUM_P
          const string posS ("  " + to_string (curIndex + 1) + " / " + to_string (tt. rows. size ()));
        #else
          const string posS ("  [Row " + to_string (curIndex + 1) + "/" + to_string (tt. rows. size ()) + "  Col " + to_string (curCol + 1) + "/" + to_string (tt. header. size ()) + "]");
        #endif
          if (nc. col_max > posS. size ())
            addstr ((pad (keyS, nc. col_max - posS. size (), efalse) + posS). c_str ());
          else
            addstr (posS. c_str ());
        }
        FOR_START (size_t, i, topIndex, bottomIndex)
        {
          ::move ((int) (i - topIndex + headerSize), 0);
          const Attr attrCurrent (A_REVERSE, i == curIndex);
          const Attr attrFound (A_BOLD, rowFound [i]);
          StringVector values;
          FFOR_START (size_t, j, curCol, tt. header. size ())
            values << tt. rows [i] [j];
          printRow (false, values, active, curCol, tt. header, nc. col_max, tt. rows. size ());
        }
        FFOR_START (size_t, i, bottomIndex, bottomIndex_max)
        {
          ::move ((int) (i - topIndex + headerSize), 0);
          ::clrtoeol ();
        }
        ::refresh ();
      }

      bool keyAccepted = false;
      while (! keyAccepted)
      {
        const int key = NCurses_sp::getKey ();  
        keyAccepted = true;
        // Show truncated field ??
        switch (key)  
        {
          case 'q':  
            quit = true;
            break;
          case KEY_RESIZE:
            break;
          case KEY_DOWN:
            if (curIndex + 1 < tt. rows. size ())
            {
              curIndex++;
              if (curIndex == bottomIndex)
                topIndex++;
            }
            else
              ::beep ();
            break;
          case ctrl('d'): 
            if (topIndex < curIndex)
              topIndex++;
            else
              ::beep ();
            break;
          case KEY_UP:
            if (curIndex)
            {
              if (curIndex == topIndex)
                topIndex--;
              curIndex--;
            }
            else
              ::beep ();
            break;
          case ctrl('u'): 
            if (   topIndex 
                && curIndex < bottomIndex
               )
              topIndex--;
            else
              ::beep ();
            break;
          case KEY_NPAGE:
            if (curIndex + 1 == tt. rows. size ())
              ::beep ();
            else if (curIndex + 1 < bottomIndex)
              curIndex = bottomIndex - 1;
            else
            {
              curIndex = min (tt. rows. size (), bottomIndex + pageScroll) - 1;
              topIndex = curIndex - pageScroll;
            }
            break;
          case KEY_PPAGE:
            if (! curIndex)
              ::beep ();
            else if (curIndex > topIndex)
              curIndex = topIndex;
            else
            {
              if (topIndex > pageScroll)
                topIndex -= pageScroll;
              else
                topIndex = 0;
              curIndex = topIndex;
            }
            break;
          case KEY_HOME:
            if (! curIndex)
              ::beep ();
            else
            {
              curIndex = 0;
              topIndex = 0;
            }
            break;
          case KEY_END:
            if (curIndex == tt. rows. size () - 1)
              ::beep ();
            else
            {
              curIndex = tt. rows. size () - 1;
              topIndex = tt. rows. size () >= fieldSize ? tt. rows. size () - fieldSize : 0;
            }
            break;
          case KEY_LEFT:
            if (! moveLeft (curCol, active))
              ::beep ();            
            break;
          case KEY_RIGHT:
            if (! moveRight (curCol, active))
              ::beep ();            
            break;
          case KEY_F(3):  // search Form ??!s
            {
              constexpr size_t size = 128;  // PAR
              ASSERT (what. size () <= size);
              char search [size] = "";
              ::echo ();
              ::curs_set (1);
              ::move ((int) (fieldSize + headerSize), 0);
              ::clrtoeol ();
              ::getstr (search);
              ::curs_set (0);
              ::noecho ();
              bool newSearch = false;
              if (*search && what != string (search))
              {
                what = search;
                FFOR (size_t, i, rowFound. size ())
                  rowFound [i] = false;
                newSearch = true;
              }
              if (what. empty ())
              {
                ::beep ();
                continue;
              }
              bool globalFound = false;
              bool curIndexSet = false;
              FFOR (size_t, i, tt. rows. size ())
              {
                bool found = false;
                const StringVector& row = tt. rows [i];
                for (const string& s : row)
                  if (contains (s, what))
                  {
                    found = true;
                    break;
                  }
                if (found)
                {
                  if (   ! curIndexSet
                      && i >= (newSearch ? 0 : (curIndex + rowFound [curIndex]))
                     )
                  {
                    curIndex = i;
                    curIndexSet = true;
                  }
                  rowFound [i] = true;
                  globalFound = true;
                }
              }
              if (! globalFound)
              {
                ::beep ();
                continue;
              }
              if (curIndex < topIndex)
                topIndex = curIndex;
              else if (curIndex >= bottomIndex/*_max*/)
                topIndex = curIndex - pageScroll;
            }
            break;
          case 'm':
            rowFound [curIndex] = ! rowFound [curIndex];
            break;
          case 'a':
            {
              const bool on = rowFound [0];
              FFOR (size_t, i, rowFound. size ())
                rowFound [i] = ! on;
            }
            break;
          case '-':
            {
              size_t actives = 0;
              for (const bool b : active)
                if (b)
                  actives++;
              ASSERT (actives);
              if (actives == 1)
                ::beep ();
              else
              {
                active [curCol] = false;
                if (! moveRight (curCol, active))
                  EXEC_ASSERT (moveLeft (curCol, active));
              }
            }
            break;
          case '+':
            FFOR (size_t, i, active. size ())
              active [i] = true;
            break;
        #ifndef NUM_P
          case '#':   
            toggle (numP);
            break;
        #endif
          default:
            keyAccepted = false;
            ::beep ();
            break;
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
