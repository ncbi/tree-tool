// xml_view.cpp

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
*   XML viewer
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "xml.hpp"
#include "../ncurses.hpp"
#include "../version.inc"




namespace
{
  
  
struct Row
{
  const Xml_sp::Data* data {nullptr};
    // !nullptr
  size_t childNum {0};
  size_t nodes {0};
  // Variable
  bool open {false};
  bool found {false};
  NCurses::Color color {NCurses::colorNone};


  Row (const Xml_sp::Data *data_arg,
       size_t childNum_arg)
    : data (data_arg)
    , childNum (childNum_arg)
    , nodes (data->getNodes ())
    { ASSERT (data); }


  size_t getDepth () const
    { return data->getDepth (); }
};




struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("View an XML file.\n\
At line ends: [<# children>|<# nodes in subtree>]\
")
  	{
      version = VERSION;
  	  addPositional ("xml", "XML file");
  	}



	void body () const final
	{
		const string xmlFName  = getArg ("xml");


	  unique_ptr<const Xml_sp::Data> xml;
	  VectorOwn<Xml_sp::Data> markupDeclarations;
	  {
  	  TokenInput ti (xmlFName, '\0', true, false, 1000);  // PAR
      try
      {
        Unverbose unv;
        xml. reset (new Xml_sp::Data (ti, markupDeclarations));
      }
    //catch (const CharInput::Error &e)
      //{ throw e; }
      catch (const exception &e)
        { ti. error (e. what (), false); }
    }
    ASSERT (xml. get ());
    xml->qc ();
    if (verbose ())
    {
      Xml::File f (cout, false, false, "XML");
      xml->saveXml (f);
    }


    Vector<Row> rows;  rows. reserve (10000);  // PAR
    rows << Row (xml. get (), 0);
    size_t topIndex = 0;
    size_t curIndex = topIndex;
    string what;  // For search
    NCurses nc (true);
    bool quit = false;
    while (! quit)
    {
      ASSERT (! rows. empty ());
      if (qc_on)
        FFOR (size_t, i, rows. size () - 1)
        {
          QC_ASSERT (rows [i + 1]. getDepth () <= rows [i]. getDepth () + 1);
          QC_IMPLY (rows [i + 1]. getDepth () == rows [i]. getDepth () + 1, rows [i]. open);
        }
      nc. resize ();
      const size_t fieldSize = nc. row_max - 1;  // Last row is for menu
      const size_t pageScroll = fieldSize - 1;
      const size_t bottomIndex_max = topIndex + fieldSize;
      const size_t bottomIndex = min (rows. size (), bottomIndex_max);
      if (   nc. row_max > 2
          && nc. col_max > 10  // PAR
         )
      {
        ASSERT (topIndex <= curIndex);
        ASSERT (topIndex < bottomIndex);
        minimize (curIndex, bottomIndex - 1);
        move ((int) fieldSize, 0);
        {
          const NCAttr attr (A_REVERSE);
          const NCBackground bkgr (nc. background | A_REVERSE);
          addstr (("[" + getFileName (xmlFName) + "] "). c_str ());
          // Most of keys are intercepted by the terminal
          addstr ("Up   Down   PgUp,b   PgDn,f   Home,B   End,F   Enter:Open/Close   F3,s:Search word from cursor");
          if (nc. hasColors)
            addstr ("   c:color");
          addstr ("   F10,q:Quit");
          clrtoeol ();
            // "F1,h": explain "[%d/%d]" at the end of lines
        }
        FOR_START (size_t, i, topIndex, bottomIndex)
        {
          const Row& row = rows [i];
          move ((int) (i - topIndex), 0);
          const size_t x = row. getDepth () * 2;  // PAR
          FFOR (size_t, j, x)
            addch (' ');
          if (row. data->children. empty ())
            addch (' ');
          else
            addch (row. open ? '-' : '+');
          addch (' ');
          const NCAttrColor attrColor (row. color, row. color != NCurses::colorNone);
          const NCAttr attrCurrent (A_REVERSE, i == curIndex);
          printw ("%lu <%s>", row. childNum + 1, row. data->name. c_str ());
          if (! row. data->token. empty ())
          {
            addch (' ');
            const NCAttr attrFound (A_BOLD, row. found);
            printw ("%s", row. data->token. str (). c_str ());
          }
          if (const size_t n = row. data->children. size ())
            printw (" [%lu/%lu]", n, row. nodes);
        #if 0
          printw (" %s %s %s %s"
                 , to_string (rows. size ()). c_str ()
                 , to_string (topIndex). c_str ()
                 , to_string (bottomIndex). c_str ()
                 , to_string (nc. row_max). c_str ()
                 );
        #endif
          // Too long line ??
          clrtoeol ();
        }
        FFOR_START (size_t, i, bottomIndex, bottomIndex_max)
        {
          move ((int) (i - topIndex), 0);
          clrtoeol ();
        }
        refresh ();
      }

      bool keyAccepted = false;
      while (! keyAccepted)
      {
        const int key = getch ();  // Invokes refresh()
        keyAccepted = true;
        switch (key)
        {
          case 'q':   // ESC
          case KEY_F(10):
            quit = true;
            break;
          case KEY_DOWN:
            if (curIndex + 1 < rows. size ())
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
            if (curIndex + 1 == rows. size ())
              beep ();
            else if (curIndex + 1 < bottomIndex)
              curIndex = bottomIndex - 1;
            else
            {
              curIndex = min (rows. size (), bottomIndex + pageScroll) - 1;
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
            if (curIndex == rows. size () - 1)
              beep ();
            else
            {
              curIndex = rows. size () - 1;
              topIndex = rows. size () >= fieldSize ? rows. size () - fieldSize : 0;
            }
            break;
          case 10:
          case KEY_ENTER:
            if (rows [curIndex]. data->children. empty ())
              beep ();
            else
            {
              toggle (rows [curIndex]. open);
              auto itStart = rows. begin ();
              advance (itStart, curIndex + 1);
              if (rows [curIndex]. open)
              {
                Vector<Row> newRows;
                const VectorOwn<Xml_sp::Data>& children = rows [curIndex]. data->children;
                newRows. reserve (children. size ());
                FFOR (size_t, i, children. size ())
                  newRows << Row (children [i], i);
                rows. insert (itStart, newRows. begin (), newRows. end ());
              }
              else
              {
                size_t i = curIndex + 1;
                const size_t depth = rows [curIndex]. getDepth ();
                for (; i < rows. size (); i++)
                  if (rows [i]. getDepth () <= depth)
                    break;
                auto itEnd = rows. begin ();
                advance (itEnd, i);
                rows. erase (itStart, itEnd);
              }
            }
            break;
          case 's':
          case KEY_F(3):
            {
              constexpr size_t size = 128;  // PAR
              ASSERT (what. size () <= size);
              char search [size] = "";
              echo ();
              curs_set (1);
              move ((int) fieldSize, 0);
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
              if (*search && what != string (search))
              {
                what = search;
                for (Row& row : rows)
                  row. found = false;
              }
              if (what. empty ())
              {
                beep ();
                continue;
              }
              bool equalName = false;  // PAR ??
              bool tokenSubstr = false; // PAR ??
              bool tokenWord = true;  // PAR ??
              size_t i = curIndex + rows [curIndex]. found;
              VectorPtr<Xml_sp::Data> path;
              for (; i < rows. size (); i++)
              {
                const Row& row = rows [i];
                if (row. open)
                {
                  if (row. data->contains (what, equalName, tokenSubstr, tokenWord))
                    break;
                }
                else
                {
                  if (row. data->find (path, what, equalName, tokenSubstr, tokenWord))
                    break;
                }
              }
              if (i == rows. size ())
              {
                beep ();
                continue;
              }
              while (! path. empty ())
              {
                ASSERT (! rows [i]. open);
                rows [i]. open = true;
                // Cf. "Open"
                auto itStart = rows. begin ();
                advance (itStart, i + 1);
                Vector<Row> newRows;
                const VectorOwn<Xml_sp::Data>& children = rows [i]. data->children;
                newRows. reserve (children. size ());
                FFOR (size_t, k, children. size ())
                  newRows << Row (children [k], k);
                rows. insert (itStart, newRows. begin (), newRows. end ());
                const size_t j = children. indexOf (path. back ());
                ASSERT (j != no_index);
                i += 1 + j;
                path. pop ();
              }
              ASSERT (i >= curIndex);
              ASSERT (i < rows. size ());
              curIndex = i;
              rows [curIndex]. found = true;
              if (curIndex >= bottomIndex_max)
                topIndex = curIndex - pageScroll;
            }
            break;
          case 'c':  // Row::color
            {
              // draw a menu ??
              Row& row = rows [curIndex];
              switch (getch ())
              {
                case 'n': row. color = NCurses::colorNone; break;
                case 'r': row. color = NCurses::colorRed; break;
                case 'g': row. color = NCurses::colorGreen; break;
                case 'y': row. color = NCurses::colorYellow; break;
                case 'b': row. color = NCurses::colorBlue; break;
                case 'm': row. color = NCurses::colorMagenta; break;
                case 'c': row. color = NCurses::colorCyan; break;
              }
            }
            break;
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
