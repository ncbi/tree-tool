// xml.hpp

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
  Common_sp::AutoPtr<const XmlSentense> sentense;


  XmlTag (istream &is,
          string phrase);
  XmlTag (const XmlTag &other)
    : name (other. name)
    , text (other. text)
    , sentense (var_cast (other). sentense)
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
