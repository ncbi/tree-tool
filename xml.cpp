// xml.cpp

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
*   XML utilities
*
*/


#undef NDEBUG
#include "common.inc"

#include "xml.hpp"



namespace XML_sp
{



string readPhrase (istream &is)
// Return: empty() | <[^ ].*[^>] | [^<>].* | >
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
  comment = isLeft (phrase, "<!--");
  whole = comment
            ? isRight (phrase, "--")  
            : isRight (phrase, "/");
  name = findSplit (phrase). substr (1);
  text = phrase;
  if (whole)
    return;
  phrase = readPhrase (is);
  if (isLeft (phrase, "<"))
  {
    auto tags = new XmlTags;
    sentense = tags;
    while (! isEnd (phrase))
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
    ASSERT (isEnd (phrase));
  }
  ASSERT (sentense);
}

    

void XmlTag::printSiblings (const string &name1,
                            const string &name2) const
{
  if (! sentense. get ())
    return;
  const XmlTags* tags = sentense->asXmlTags ();
  if (! tags)
    return;
    
  StringVector vec1;
  StringVector vec2;
  for (const XmlTag& tag : tags->vec)
         if (tag. name == name1)  vec1 << checkPtr (tag. sentense->asXmlText ()) -> text;
    else if (tag. name == name2)  vec2 << checkPtr (tag. sentense->asXmlText ()) -> text;
    else 
      tag. printSiblings (name1, name2);      
  
  for (const string& s1 : vec1)
  for (const string& s2 : vec2)
    cout << s1 << '\t' << s2 << endl;
}

  

XmlTag readXml (istream &is)
{
  string s;
  readLine (is, s);
  ASSERT (s == "<?xml version=\"1.0\"?>");
  return XmlTag (is, readPhrase (is));
}



}
