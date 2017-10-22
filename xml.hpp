// xml.hpp

#ifndef XML_HPP_78429
#define XML_HPP_78429

#include "common.hpp"
using namespace Common_sp;



namespace XML_sp
{


struct XmlText;
struct XmlTags;



struct XmlSentense : Root
{
  virtual const XmlTags* asXmlTags () const 
    { return nullptr; }
  virtual const XmlText* asXmlText () const 
    { return nullptr; }
};



struct XmlText : XmlSentense
{
  string text;
  
  explicit XmlText (const string &text_arg)
    : text (text_arg)
    {}

  const XmlText* asXmlText () const final
    { return this; }
};



struct XmlTag : Root
{
  string name;
    // empty() <=> comment
  string text;
  bool comment {false};
  bool whole {false};
  Common_sp::AutoPtr <const XmlSentense> sentense;


  XmlTag (istream &is,
          string phrase);
  XmlTag (const XmlTag &other)
    : name (other. name)
    , text (other. text)
    , sentense (const_cast <XmlTag&> (other). sentense)
    {}


  bool isEnd (const string& phrase) const
    { if (comment)
        return isRight (phrase, "--");
      return phrase == "</" + name;
    }
  void printSiblings (const string &name1,
                      const string &name2) const;
};



struct XmlTags : XmlSentense
{
  Vector<XmlTag> vec;

  const XmlTags* asXmlTags () const final
    { return this; }
};



XmlTag readXml (istream &is);



}  



#endif
