// xml_view.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "xml.hpp"

// ftp://ftp.gnu.org/gnu/ncurses/ncurses-6.2.tar.gz
// Old: ftp://ftp.gnu.org/pub/gnu/ncurses/ncurses.tar.gz
extern "C"
{
  #include <ncurses.h>
}



namespace
{


struct Row
{
  const Xml_sp::Data* data {nullptr};
    // !nullptr
  size_t nodes {0};
  // Variable
  bool open {false};
  bool found {false};

  explicit Row (const Xml_sp::Data *data_arg)
    : data (data_arg)
    , nodes (data->getNodes ())
    { ASSERT (data); }

  size_t getDepth () const
    { return data->getDepth (); }
  bool operator== (const Row &other) const
    { return data == other. data; }
};



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("View an XML file")
  	{
  	  addPositional ("xml", "XML file");
  	}



	void body () const final
	{
		const string xmlFName  = getArg ("xml");


	  unique_ptr<const Xml_sp::Data> xml;
	  {
  	  TokenInput ti (xmlFName, '\0', 100 * 1024, 1000);  // PAR
      try
      {
        xml. reset (new Xml_sp::Data (ti));
      }
      catch (const CharInput::Error &e)
        { throw e; }
      catch (const exception &e)
        { ti. error (e. what (), false); }
    }
    ASSERT (xml. get ());
    xml->qc ();


    Vector<Row> rows;  rows. reserve (10000);  // PAR
    rows << Row (xml. get ());
    size_t topIndex = 0;
    size_t curIndex = topIndex;
    int row_max = 0;
    int col_max = 0;
    string what;  // For search
    initscr ();
    cbreak ();
    noecho ();
    keypad (stdscr, TRUE);
    curs_set (0);  // does not work ??
    wclear (stdscr);
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
      getmaxyx (stdscr, row_max, col_max);
      const size_t fieldSize = (size_t) (row_max - 1);  // Last row is for menu
      const size_t pageScroll = fieldSize - 1;
      const size_t bottomIndex_max = topIndex + fieldSize;
      const size_t bottomIndex = min (rows. size (), bottomIndex_max);
      if (   row_max > 2
          && col_max > 10  // PAR
         )
      {
        ASSERT (topIndex <= curIndex);
        ASSERT (topIndex < bottomIndex);
        minimize (curIndex, bottomIndex - 1);
        move ((int) fieldSize, 0);
        attron (A_DIM);
        addstr (("[" + getFileName (xmlFName) + "]   "). c_str ());
        addstr ("Up   Down   PgUp,b   PgDn,f   Home,B   End,F   Enter:Open/Close   F3,s:Search word from cursor   F10,q:Quit");
          // "F1,h": explain "[%d/%d]"
          // Most of keys are intercepted by treminal
        attroff (A_DIM);
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
          if (i == curIndex)
            attron (A_REVERSE);
          printw ("<%s>", row. data->name. c_str ());
          if (! row. data->token. empty ())
          {
            addch (' ');
            if (row. found)
          	  attron (A_BOLD | A_UNDERLINE);
            printw ("%s", row. data->token. str (). c_str ());
            if (row. found)
          	  attroff (A_BOLD | A_UNDERLINE);
          }
          if (const size_t n = row. data->children. size ())
            printw (" [%s/%s]", to_string (n). c_str (), to_string (row. nodes). c_str ());
        #if 0
          printw (" %s %s %s %s"
                 , to_string (rows. size ()). c_str ()
                 , to_string (topIndex). c_str ()
                 , to_string (bottomIndex). c_str ()
                 , to_string (row_max). c_str ()
                 );
        #endif
          // Too long line ??
          clrtoeol ();
          if (i == curIndex)
            attroff (A_REVERSE);
        }
        FFOR_START (size_t, i, bottomIndex, fieldSize)
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
                for (const Xml_sp::Data* child : children)
                  newRows << Row (child);
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
              sprintf (format, ("%" + to_string (size) + "s"). c_str ());
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
              size_t i = curIndex + (rows [curIndex]. found ? 1 : 0);
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
                for (const Xml_sp::Data* child : children)
                  newRows << Row (child);
                const Row back (path. back ());
                const size_t j = newRows. indexOf (back);
                ASSERT (j != no_index);
                rows. insert (itStart, newRows. begin (), newRows. end ());
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
          default:
            keyAccepted = false;
            break;
        }
      }
    }
    curs_set (1);
    endwin ();
	}
};



}  // namespace



int main (int argc,
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}
