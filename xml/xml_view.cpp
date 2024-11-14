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


#define USE_WCHAR 0  
#if USE_WCHAR
  #include <codecvt>
#endif


#undef NDEBUG  

#include "../common.hpp"
using namespace Common_sp;
#include "../ncurses.hpp"
using namespace NCurses_sp;
#include "../version.inc"
#include "xml.hpp"

#include "../common.inc"



#define CTRL(x)     ((x) & 0x1f)



namespace
{
	
	
constexpr size_t token_name_len = 100;   // PAR

  
  
struct Row
{
  const Xml_sp::Data* data {nullptr};
    // !nullptr
  size_t childNum {0};
  size_t nodes {0};
  // Variable
  bool open {false};
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




NCurses::Color mask2color (uchar mask)
{
  if (! mask)
  	return NCurses::colorNone;
	switch (byte2first (mask) % 5)
	{
		case 0: return NCurses::colorRed;
		case 1: return NCurses::colorYellow;
		case 2: return NCurses::colorMagenta;
		case 3: return NCurses::colorCyan;
		case 4: return NCurses::colorWhite;
	}
	NEVER_CALL;
}	



struct Viewer
{
	const string xmlFName;
	const Xml_sp::Data& xml;
  Vector<Row> rows;  
  size_t topIndex {0};
  size_t curIndex {0};
    // rows.size() < curIndex <= topIndex 
  // Search
  StringVector what;  // Of UTF8 characters    
  bool targetNameP {true};  
  bool substrP {true}; 
  bool wordP {true};  
  //
  NCurses nc;


  Viewer (const string& xmlFName_arg,
          const Xml_sp::Data* xml_arg)
  	: xmlFName (xmlFName_arg)
  	, xml (*xml_arg)
  	, nc (true)
  	{ 
  		ASSERT (xml_arg);
  		rows. reserve (10000);  // PAR
		  rows << Row (xml_arg, 0);
  	}
  	
  	
  void run ()
  {
    setTopIndex (nc. row_max - 1, open (curIndex));
	  bool quit = false;
    while (! quit)
    {
      // Input: topIndex, curIndex
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
        ASSERT (topIndex < bottomIndex);
        maximize (curIndex, topIndex);
        minimize (curIndex, bottomIndex - 1);
        drawMenu (fieldSize, "[" + getFileName (xmlFName) + "]  Up  Down  PgUp  PgDn  Home  End  ^F  Enter:Open/Close  a:Open all  F3:Search  n:Next found" + ifS (nc. hasColors, "  c:color") + "  q:Quit");
          // uncolor all ??
          // Complex keys are trapped by the treminal
          // "h": explain [a/b] ??
        FOR_START (size_t, i, topIndex, bottomIndex)
        {
          const Row& row = rows [i];
          move ((int) (i - topIndex), 0);
          {
            const bool current = (i == curIndex);
	          const Attr attr (A_REVERSE, current);
	          const size_t x = row. getDepth ();  
	          FFOR (size_t, j, x)
	          {
	            addch ('|');
	            addch (' ');
	          }
	          if (row. data->children. empty ())
	            addch ('|');
	          else
	            addch (row. open ? '-' : '+');
	        }
          addch (' ');
          {
	          const AttrColor attrColor (NCurses::colorBlue); 
	          size_t width = 0;
	          if (const Xml_sp::Data* parent = row. data->parent)
	          	width = to_string (parent->children. size ()). size ();
            printw ("%*lu", (int) width, row. childNum + 1);
          }
          {
	          const AttrColor attrColor_tag (NCurses::colorGreen); 
	          const uchar mask = row. data->searchFound 
                                ? row. data->searchFound 
                                : row. open 
                               	  ? 0
                                	: row. data->searchFoundAll;
	          const AttrColor attrColor_found (mask2color (mask), mask);  
	          const AttrColor attrColor_color (row. color, row. color != NCurses::colorNone); 
	          printw (" <%s>", row. data->getName (). c_str ());
	        }
          if (! row. data->token. empty ())
          {
          #if USE_WCHAR
          	try
          	{
          		// wchar_t stores character in UTF-16
              std::wstring wstr (converter.from_bytes (row. data->token. name));
	            printw (" %ls", wstr. c_str ());  
            }
            catch (const exception &e)
            {
	            printw (" %s", (row. data->token. name + "  (ERROR: " + string (e. what ()) + ")"). c_str ());
	          #if 0
	            for (const char c : row. data->token. name)
	              printw ("  %x", c & 0xFF);
	          #endif
            }
          #else
            // emphacise current row ?? 
            string s (row. data->token. name);
            size_t j = 0;
            while (j < s. size () && s [j] == ' ')
              j++;
            const string sp ("<SPACES:" + to_string (j) + ">");
            if (j > sp. size ())
              s = sp + s. substr (j);
            if (s. size () >= token_name_len)  // PAR
            {
            	s. erase (token_name_len);
            	s += " ...";
            }
            printw (" %s", s. c_str ());  
          #endif
          }
          if (const size_t n = row. data->children. size ())
          {
            printw (" ");
	          const AttrColor attrColor_suf (NCurses::colorBlue); 
            printw (" %lu/%lu", n, row. nodes);
          }
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
        const int key = getKey ();
        keyAccepted = true;
        switch (key)
        {
          case 'q':  
            quit = true;
            break;
          case KEY_RESIZE:
            break;
          case KEY_DOWN:
            if (curIndex + 1 < rows. size ())
            {
              curIndex++;
              if (curIndex == bottomIndex)
                topIndex++;
            }
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
          case 01016:  // ^KEY_DOWN
          	if (topIndex + 1 < rows. size ())
          		topIndex++;
          	else 
          		::beep ();
          	break;
          case 01067:  // ^KEY_UP
          	if (topIndex)
          		topIndex--;
          	else
          		::beep ();
          	break;
          case KEY_NPAGE:
            if (curIndex + 1 == rows. size ())
              ::beep ();
            else if (curIndex + 1 < bottomIndex)
              curIndex = bottomIndex - 1;
            else
            {
              curIndex = min (rows. size (), bottomIndex + pageScroll) - 1;
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
            if (curIndex == rows. size () - 1)
              ::beep ();
            else
            {
              curIndex = rows. size () - 1;
              topIndex = rows. size () >= fieldSize ? rows. size () - fieldSize : 0;
            }
            break;
          case CTRL('f'):
            do
              curIndex = rows. size () - 1;
            while (open (curIndex));
            topIndex = rows. size () >= fieldSize ? rows. size () - fieldSize : 0;
            break;
          case 10:         // Key above "shift" key
          case KEY_ENTER:  // Bottom-right key
            if (rows [curIndex]. data->children. empty ())
              ::beep ();
            else if (! rows [curIndex]. open)
              setTopIndex (fieldSize, open (curIndex));
            else
            {
              rows [curIndex]. open = false;
              auto itStart = rows. begin ();
              advance (itStart, curIndex + 1);
              size_t i = curIndex + 1;
              const size_t depth = rows [curIndex]. getDepth ();
              for (; i < rows. size (); i++)
                if (rows [i]. getDepth () <= depth)
                  break;
              auto itEnd = rows. begin ();
              advance (itEnd, i);
              rows. erase (itStart, itEnd);
            }
            break;
          case 'a':  // Open _a_ll
            if (   rows [curIndex]. data->children. empty ()
                || rows [curIndex]. open
               )
              ::beep ();
            else
              setTopIndex (fieldSize, openAll ());
            break;
          case KEY_F(3):  
            if (   searchForm ()
            	  && ! what. empty ()
            	 )
            {
              const uchar mask = 1;  // ??
              var_cast (xml). unsetSearchFound (mask);
              // StringMatch::Type ??
              const StringMatch::Type t = wordP 
              	                            ? StringMatch::word 
              	                            : substrP 
              	                            	 ? StringMatch::part 
              	                            	 : StringMatch::whole;
            //message (to_string ((int) t));  
              if (! var_cast (xml). setSearchFound (what. toString (), targetNameP, t, mask))
              	::beep ();
            }
            break;
          case 'n':  // Next found row 
          	{
          		size_t curIndex_new = curIndex;
          		if (rows [curIndex_new]. data->searchFound)
          			curIndex_new++;
          		if (openSearchFound (curIndex_new))
          		{
          			curIndex = curIndex_new;
          		#if 1
          			setTopIndex (fieldSize, 0);
          		#else
       			    if (curIndex >= topIndex + fieldSize)
                   topIndex = curIndex - fieldSize + 1;
              #endif
          		}
          		else
          			::beep ();
          	}
          	break;
          case 'c':  // Row::color
            {
              drawMenu (fieldSize, "n(one) r(ed) g(reen) y(ellow) b(lue) m(agenta) c(yan) w(hite)");
              Row& row = rows [curIndex];
              switch (getKey ())
              {
                case 'n': row. color = NCurses::colorNone; break;
                case 'r': row. color = NCurses::colorRed; break;
                case 'g': row. color = NCurses::colorGreen; break;
                case 'y': row. color = NCurses::colorYellow; break;
                case 'b': row. color = NCurses::colorBlue; break;
                case 'm': row. color = NCurses::colorMagenta; break;
                case 'c': row. color = NCurses::colorCyan; break;
                case 'w': row. color = NCurses::colorWhite; break;  
              }
            }
            break;
          default:
            keyAccepted = false;
            ::beep ();
            break;
        }
      }
    }
  }


  void setTopIndex (size_t fieldSize,
                    size_t opened)
  {
    if (curIndex + opened >= topIndex + fieldSize)
      topIndex = curIndex + opened - fieldSize + 1;
    if (topIndex > curIndex)
      topIndex = curIndex;
  }


	void drawMenu (size_t fieldSize, 
	               const string &s)
	{
	  move ((int) fieldSize, 0);
	  const Attr attr (A_REVERSE);
	  const Background bkgr (nc. background | A_REVERSE);
	  addstr (s. c_str ());
	  clrtoeol ();
	}


	size_t open (size_t index)
	{
	  ASSERT (! rows [index]. open);
	  const VectorOwn<Xml_sp::Data>& children = rows [index]. data->children;
	  if (children. empty ())
	    return 0;
	  rows [index]. open = true;
	  auto itStart = rows. begin ();
	  advance (itStart, index + 1);
	  Vector<Row> newRows;
	  newRows. reserve (children. size ());
	  FFOR (size_t, i, children. size ())
	    newRows << Row (children [i], i);
	  rows. insert (itStart, newRows. begin (), newRows. end ());
	  return newRows. size ();
	}
	
		
	size_t openAll ()
	{
	  size_t opened = 0;
	  FOR_START (size_t, i, curIndex, curIndex + 1 + opened)
	    opened += open (i);
	  return opened;
	}
	
	
	bool openSearchFound (size_t &curIndex_new)
	{
		while (curIndex_new < rows. size ())
		{
			const Row& row = rows [curIndex_new];
			if (row. data->searchFound)
				return true;
			if (row. data->searchFoundAll && ! row. open)
				open (curIndex_new);
			curIndex_new++;
	  }
		return false;
	}
	
	
	bool searchForm () 
	{
		constexpr size_t w0 = 20;
		constexpr size_t w1 = 5;
		constexpr size_t w2 = 5;
		constexpr size_t w3 = 5;
		constexpr size_t c0 = 2;
		constexpr size_t c1 = c0 + w0 + 1;
		constexpr size_t c2 = c1 + w1 + 1;
		constexpr size_t c3 = c2 + w2 + 1;
		constexpr size_t c4 = c3 + w3 + 1;
		Form form (nc, c4, 4);
		form. print (c4 / 2 - 3, 0, " Find ");
		form. print (c0, 1, "What");  // what_xml, what_content ??
		form. print (c1, 1, "Tag");
		form. print (c2, 1, "Whole");
		form. print (c3, 1, "Word");
		form. add (c0, 2, w0, false, what);
		form. add (c1, 2, 2, true, yesNo (targetNameP));
		form. add (c2, 2, 2, true, yesNo (! substrP));
		form. add (c3, 2, 2, true, yesNo (wordP));
		if (! form. run ())
			return false;
		
		// validation ??
			
		what        = form. getVal (0);
		targetNameP = form. getS (1) == "Y";
		// --> StringMatch::Type drop-down list ??
		substrP = form. getS (2) == "N";
		wordP   = form. getS (3) == "Y";
		
	  return true;	
	}
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
  	  addFlag ("bin", "Binary XML");
  	  addKey ("search_tags", "Tag names to be searched, comma-delimited");
  	}



	void body () const final
	{
		const string xmlFName = getArg ("xml");
		const bool   bin      = getFlag ("bin"); 
		const StringVector searchTags (getArg ("search_tags"), ',', true);

    QC_ASSERT (searchTags. size () < 8);  // instruction ??


    EXEC_ASSERT (setlocale (LC_ALL, "en_US.UTF-8"));  
  #if USE_WCHAR
    std::wstring_convert <std::codecvt_utf8_utf16<wchar_t>> converter;
  #endif


    Names names (10000);   // PAR
    Xml_sp::Data* xml = nullptr;
	  unique_ptr<const Xml_sp::Data> xml_;
	  if (bin)
	  	xml = Xml_sp::Data::load (names, xmlFName);
	  else
	  {
		  VectorOwn<Xml_sp::Data> markupDeclarations;
	  	xml = Xml_sp::Data::load (names, xmlFName, markupDeclarations);
	  }
    ASSERT (xml);
    xml_. reset (xml);
    
    // PAR
    xml->text2token ("text");  
    xml->chars2token ("char", token_name_len);
    xml->mergeSingleChildren ();  
    xml->visualizeTokenTrailingSpaces ();
    
    Viewer vw (xmlFName, xml);

    if (! searchTags. empty ())
    {
    	uchar mask = 1;  // NCurses::Color bits
    	for (const string& s : searchTags)
    		if (! s. empty ())
	    	{
	    		// Display s/color in search form window ??
	    		// Use search fields of vw ??
	        xml->setSearchFound (s, true, StringMatch::word, mask);
	        if (vw. what. empty ())
	        	for (const char c : s)
	        	  vw. what << string (1, c);
	    	//mask = uchar (mask * 2);  // White color is bad
	      }
    }
    
    xml->qc ();
    if (verbose ())
    {
      Xml::TextFile f ("xml_view.xml", "XML");
      xml->saveXml (f);
    }

    vw. run ();
	}
};



}  // namespace



int main (int argc,
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}
