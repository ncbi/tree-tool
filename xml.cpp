// xml.cpp

#undef NDEBUG
#include "common.inc"

#include "xml.hpp"



namespace XML_sp
{
using namespace std;
using namespace Common_sp;



string readPhrase (istream &is)
{
  char c = '\0';

  while (getChar (is, c) && isspace (c));
  if (isspace (c))
    return string ();

  string s (1, c);
  if (c == '<')
  {
    while (getChar (is, c) && c !='>')
      s += c;
    if (c != '>')
      throw runtime_error ("No '>'");
    if (s. size () == 1 || isspace (s [1]))
      throw runtime_error ("Empty tag");      
  }
  else
  {
    while (getChar (is, c) && c !='<')
    {
      if (c == '>')
        throw runtime_error ("'>' without preceding '<'");
      s += c;
    }
    if (c == '<')
      is. unget ();
  }
  return s;
}




XmlTag::XmlTag (istream &is,
                string phrase)
{ 
  ASSERT (! phrase. empty ());
  ASSERT (phrase [0] == '<');
  ASSERT (phrase. size () > 1);
  ASSERT (! isspace (phrase [1]));
  ASSERT (phrase [1] != '/');
  const bool comment = isLeft (phrase, "<!--");
  const bool whole = comment
                       ? isRight (phrase, "--")  
                       : isRight (phrase, "/");
  name = findSplit (phrase). substr (1);
  text = phrase;
//cout << "Start: " << name << ' ' << comment << ' ' << whole << endl;  
  if (whole)
    return;
  const string end ("</" + name);
  phrase = readPhrase (is);
  if (isLeft (phrase, "<"))
  {
    auto tags = new XmlTags;
    sentense = tags;
    while (   (  comment && ! isRight (phrase, "--"))
           || (! comment && phrase != end)
          )
    {
      if (! isLeft (phrase, "<"))
        throw runtime_error ("'<' expected");
      if (phrase. empty ())
        throw runtime_error ("Closing tag of \"" + name + "\" is not found");
      tags->vec << XmlTag (is, phrase);
      phrase = readPhrase (is);
    }
  } 
  else
  {
    sentense = new XmlText (phrase);
    phrase = readPhrase (is);
    ASSERT (isLeft (phrase, "<"));
    if (comment)
      { ASSERT (isRight (phrase, "--")); }
    else
      { ASSERT (phrase == end); }
  }
  ASSERT (sentense);
//cout << "End: " << name << endl; 
}

    

#if 0
void XmlTag::printSiblings (const string &name1,
                            const string &name2) const
{
  if (! sentense)
    return;
  const XmlTags* tags = sentense->asXmlTags ();
  if (! tags)
    return;
  StringVector vec1;
  StringVector vec2;
  for (const const XmlTag* tag : tags->vec)
    if 
}
#endif

  

XmlTag readXml (istream &is)
{
  string s;
  readLine (is, s);
  ASSERT (s == "<?xml version=\"1.0\"?>");
  return XmlTag (is, readPhrase (is));
}



}
