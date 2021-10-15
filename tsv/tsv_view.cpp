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
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../ncurses.hpp"
#include "../version.inc"



#define NUM_P  // ??



namespace
{
  
  
bool printString (string s,
                  size_t screen_col_max,
                  size_t &x)
{
  bool truncated = false;
  if (s. size () > screen_col_max - x)
  {
    s. erase (screen_col_max - x);
    truncated = true;
  }
  addstr (s. c_str ());
  x += s. size ();
  return ! truncated;
}
  
  
  
size_t printRow (bool is_header,
                 const StringVector &values,
                 size_t col_start,
                 const Vector<TextTable::Header> &header,
                 size_t screen_col_max)
// Return: Last printed column
{
  ASSERT (col_start + values. size () == header. size ());

  size_t x = 0;
  size_t lastCol = col_start;
  FFOR_START (size_t, col, col_start, header. size ())
  {
    ASSERT (x <= screen_col_max);
    const TextTable::Header& h = header [col];
    string value (values [col - col_start]);
    if (   ! is_header 
        && h. numeric 
        && ! h. scientific
        && h. decimals
       )
    {
      bool hasPoint = false;
      streamsize decimals = 0;
      TextTable::getDecimals (value, hasPoint, decimals);
      ASSERT (h. decimals >= decimals);
      // Data modification ??
      if (! hasPoint)
        value += ".";
      value += string ((size_t) (h. decimals - decimals), '0');  
    }
    ASSERT (value. size () <= h. len_max);
    if (! printString (pad (value, h. len_max, ! (h. numeric && ! h. scientific)), screen_col_max, x))
      break;
    lastCol = col;
    if (col + 1 < header. size ())
      if (! printString ("  ", screen_col_max, x))
        break;
  }

  FFOR_START (size_t, i, x, screen_col_max)
    addstr (" ");
  //clrtoeol ();  // may start erasing before the cursor
  
  return lastCol;
}

  
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("View an tsv-table")
  	{
      version = VERSION;
  	  addPositional ("table", "tsv-table");
  	}



	void body () const final
	{
		const string tableFName  = getArg ("table");


    TextTable tt (tableFName);
    tt. qc ();
    FFOR (size_t, i, tt. header. size ())
    {
      TextTable::Header& h = tt. header [i];
      maximize (h. len_max, max (h. name. size (), to_string (i + 1). size ()));
    }
    
    if (tt. rows. empty ())
    {
      cout << "No data" << endl;
      return;
    }


    size_t topIndex = 0;
    size_t curIndex = topIndex;
    size_t curCol = 0;
    size_t curLastCol = curCol;
    string what;  // For search
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
      if (   nc. row_max > headerSize + 2
          && nc. col_max > 10  // PAR
         )
      {
        ASSERT (topIndex <= curIndex);
        ASSERT (topIndex < bottomIndex);
        minimize (curIndex, bottomIndex - 1);
        {
          move (0, 0);
          const NCAttr attr (A_BOLD);
          addstr (tableFName. c_str ());
          clrtoeol ();
        }
        {
          move (1, 0);
          const NCAttr attr (A_BOLD);
          const NCBackground bkgr (COLOR_PAIR (7) /*nc. background*/ | A_BOLD);
          StringVector values;
          FFOR_START (size_t, j, curCol, tt. header. size ())
            values << tt. header [j]. name;
          curLastCol = printRow (true, values, curCol, tt. header, nc. col_max);
        }        
        if (numP)
        {
          move (2, 0);
          const NCAttr attr (COLOR_PAIR (4) /*| A_BOLD*/);
          StringVector values;
          FFOR_START (size_t, j, curCol, tt. header. size ())
            values << to_string (j + 1);
          printRow (true, values, curCol, tt. header, nc. col_max);          
        }
        move ((int) (fieldSize + headerSize), 0);
      #if 0
        // For testing
        addstr (to_string (curLastCol). c_str ());  
      #else
        {
          const NCAttr attr (A_BOLD);
          const NCBackground bkgr (COLOR_PAIR (3) /*nc. background*/ | A_BOLD);
          const string keyS ("Up  Down  Left  Right  PgUp,b  PgDn,f  Home,B  End,F  F3,s:Search from cursor"
                           #ifndef NUM_P
                             "  #:numbers"
                           #endif
                             "  F10,q:Quit"
                            );
            // Non-character keys may be intercepted by the terminal
        #ifdef NUM_P
          const string posS ("  " + to_string (curIndex + 1) + " / " + to_string (tt. rows. size ()));
        #else
          const string posS ("  [Row " + to_string (curIndex + 1) + "/" + to_string (tt. rows. size ()) + "  Col " + to_string (curCol + 1) + "/" + to_string (tt. header. size ()) + "]");
        #endif
          if (nc. col_max > posS. size ())
            addstr ((pad (keyS, nc. col_max - posS. size (), true) + posS). c_str ());
          else
            addstr (posS. c_str ());
        }
      #endif
        FOR_START (size_t, i, topIndex, bottomIndex)
        {
          move ((int) (i - topIndex + headerSize), 0);
          const NCAttr attrCurrent (A_REVERSE, i == curIndex);
          const NCAttr attrFound (A_BOLD, rowFound [i]);
          StringVector values;
          FFOR_START (size_t, j, curCol, tt. header. size ())
            values << tt. rows [i] [j];
          printRow (false, values, curCol, tt. header, nc. col_max);
        }
        FFOR_START (size_t, i, bottomIndex, bottomIndex_max)
        {
          move ((int) (i - topIndex + headerSize), 0);
          clrtoeol ();
        }
        refresh ();
      }

      bool keyAccepted = false;
      while (! keyAccepted)
      {
        const int key = getch ();  // Invokes refresh()
      #if 0
        endwin(); 
        cout << "KEY NAME: " << keyname (key) << " - " << key << endl;
        exit (0);
      #endif
        keyAccepted = true;
        switch (key)
        {
          case 'q':   // ESC
          case KEY_F(10):
            quit = true;
            break;
          case KEY_DOWN:
            if (curIndex + 1 < tt. rows. size ())
            {
              curIndex++;
              if (curIndex == bottomIndex)
                topIndex++;
            }
            else
              beep ();
            break;
          case KEY_UP:
            if (curIndex)
            {
              if (curIndex == topIndex)
                topIndex--;
              curIndex--;
            }
            else
              beep ();
            break;
          case 'f':  
          case KEY_NPAGE:
            if (curIndex + 1 == tt. rows. size ())
              beep ();
            else if (curIndex + 1 < bottomIndex)
              curIndex = bottomIndex - 1;
            else
            {
              curIndex = min (tt. rows. size (), bottomIndex + pageScroll) - 1;
              topIndex = curIndex - pageScroll;
            }
            break;
          case 'b':  
          case KEY_PPAGE:
            if (! curIndex)
              beep ();
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
          case 'B':
          case KEY_HOME:
            if (! curIndex)
              beep ();
            else
            {
              curIndex = 0;
              topIndex = 0;
            }
            break;
          case 'F':
          case KEY_END:
            if (curIndex == tt. rows. size () - 1)
              beep ();
            else
            {
              curIndex = tt. rows. size () - 1;
              topIndex = tt. rows. size () >= fieldSize ? tt. rows. size () - fieldSize : 0;
            }
            break;
          case KEY_LEFT:
            if (curCol)
              curCol--;
            else
              beep ();            
            break;
          case KEY_RIGHT:
            if (curLastCol + 1 < tt. header. size ())
              curCol++;
            else
              beep ();            
            break;
          case 's':
          case KEY_F(3):
            {
              constexpr size_t size = 128;  // PAR
              ASSERT (what. size () <= size);
              char search [size] = "";
              echo ();
              curs_set (1);
              move ((int) (fieldSize + headerSize), 0);
              clrtoeol ();
            #if 0
              char format [32];
              sprintf (format, "%c%lu%s", '%', size, "%s");
            //strcat (format, "%s");
              scanw (format, search);  // does not work
            #else
              getstr (search);
            #endif
              curs_set (0);
              noecho ();
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
                beep ();
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
                beep ();
                continue;
              }
              if (curIndex >= bottomIndex_max)
                topIndex = curIndex - pageScroll;
            }
            break;
        #ifndef NUM_P
          case '#':   
            toggle (numP);
            break;
        #endif
          default:
            keyAccepted = false;
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
