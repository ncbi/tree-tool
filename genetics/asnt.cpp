// asnt.cpp

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
*   Parse a text ASN.1 file
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../version.inc"




namespace 
{



struct Value;



struct Item : Root
{
  StringVector names;
    // !empty()
  unique_ptr<Value> value;
    // (bool)get()
    

  Item (TokenInput &in,
        char &followingDelimiter);
    // Output: followingDelimiter: '\0' | ',' | '}'
  void saveText (ostream &os) const final;
};



struct Value : Root
{
  Token t;
    // t.type = eText && t.quote = '\'' then hexadecimal
  Vector<Item> items;
    // !empty() => t.empty()


  explicit Value (const string &s)
    : t (Token (s))
    {}
  explicit Value (Token &&t_arg)
    : t (move (t_arg))
    {}
  Value () = default;
  void saveText (ostream &os) const final
    { if (t. empty ())
      { os << '{';
        bool first = true;
        for (const Item& item : items)
        { if (! first)
            os << ", ";
          item. saveText (os);
          first = false;
        }
        os << '}';
      }
      else
        t. saveText (os);
    }
};



struct Asn : Root
{
private:
	TokenInput in;
	  // ASN.1 text
	Token last;
public:
  string title;    
	unique_ptr<Item> item;
	  //(bool)get()
    

  explicit Asn (const string &fName)
		: in (fName, '\0', true)  		  
		{ const Token titleToken (in. get ());
		  QC_ASSERT (titleToken. type == Token::eName);
    	title = titleToken. name;
    	in. get (':');
    	in. get (':');
    	in. get ('=');    	
    	char followingDelimiter = '\0';
		  item. reset (new Item (in, followingDelimiter));
		}
  void saveText (ostream &os) const final
    { os << title << " ::= ";
      if (item. get ())
        item->saveText (os);
      else
        os << "<no item>";
    }
};




Item::Item (TokenInput &in,
            char &followingDelimiter)
{
  followingDelimiter = '\0';
  bool stop = false;
  while (! stop)
  {
    Token t (in. get ());
    if (t. type == Token::eName)
      names << move (t. name);
    else
    {
    //QC_ASSERT (! names. empty ());
      if (t. type == Token::eDelimiter)
      {
        ASSERT (t. name. size () == 1);
        switch (t. name [0])
        {
          case ',':
          case '}':
            if (! value. get ())
            {
              QC_ASSERT (names. size () >= 2);
              value. reset (new Value (names. pop ()));
            }
            followingDelimiter = t. name [0];
            break;
          case '{':
            QC_ASSERT (! value. get ());
            value. reset (new Value ());
            for (;;)
            { 
              char delimiterChar = '\0';
              Item item (in, delimiterChar);
              value->items << move (item);
              if (delimiterChar == '\0')
              {
                const Token delimiter (in. get ());
                QC_ASSERT (delimiter. type == Token::eDelimiter);
                delimiterChar = delimiter. name [0];
              }
              if (delimiterChar == ',')
                continue;
              if (delimiterChar == '}')
                break;
              throw runtime_error ("Unknown delimiter: " + to_string (delimiterChar));
            }
            break;
          default:
            throw runtime_error ("Unknown delimiter: " + t. name);
        }
        stop = true;
      }
      else
      {
        QC_ASSERT (! value. get ());
        value. reset (new Value (move (t)));
        if (value->t. type == Token::eText && value->t. quote == '\'')
          in. get ("H");
        stop = true;
      }
    }
  }
}



void Item::saveText (ostream &os) const 
{ 
  for (const string& s : names)
    os << s << ' ';
  if (value. get ())
    value->saveText (os);
  else
    os << "<no value>";
}



//

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Parse a text ASN.1 file")
	  {
      version = VERSION;
	  	addPositional ("in", "Text ASN.1 file");
	  }



  void body () const final
  {
    const string inFName = getArg ("in");


    const Asn asn (inFName);
    asn. saveText (cout);
  }  
};


}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



