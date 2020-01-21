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
#include "../common.inc"

#include "xml.hpp"



namespace Xml
{



// Schema

Schema* Schema::readSchema (const string &fName,
                            string &name)
{
  const StringVector vec (fName, (size_t) 100);  // PAR
  QC_ASSERT (! vec. empty ());
  
  size_t start = 0;
  return new Schema (vec, start, 0, name);
}



Schema::Schema (const StringVector &vec,
                size_t &start,
                size_t depth,
                string &name)
{
  ASSERT (start < vec. size ());
  ASSERT (isLeftBlank (vec [start], 2 * depth));
  
  istringstream iss (vec [start]);
  TokenInput ti (iss);
  
  // name
  Token nameT (ti. get ());
  QC_ASSERT (nameT. type == Token::eName);
  name = move (nameT. name);
  
  // multiple
  Token tableT (ti. get ());
  QC_ASSERT (tableT. type == Token::eName);
  if (tableT. name == "TABLE")
    multiple = true;
  else
    ti. setLast (move (tableT));
  
  // types
  Token typeT (ti. get ());
  QC_ASSERT (typeT. type == Token::eName);
  if (typeT. name == "rows")
    ti. setLast (move (typeT));
  else
  {
    const Token::Type type = Token::str2type (typeT. name);
    types << type;
    if (type == Token::eText)
    {
      // len_max
      ti. get ("size");
      ti. get (':');
      Token lenT (ti. get ());
      QC_ASSERT (lenT. type == Token::eInteger);
      len_max = (size_t) lenT. n;
    }
  }
  
  // rows
  ti. get ("rows");
  ti. get (':');
  Token rowsT (ti. get ());
  QC_ASSERT (rowsT. type == Token::eInteger);
  rows = (size_t) rowsT. n;
  
  // name2schema[]
  start++;
  depth++;
  while (   start < vec. size ()
         && isLeftBlank (vec [start], 2 * depth)
        )
  {
    string name1;
    auto sch = new Schema (vec, start, depth, name1);
    QC_ASSERT (! name1. empty ());
    name2schema [name1] = sch;
  }
}



void Schema::qc () const 
{
  if (! qc_on)
    return;
  Root::qc ();
    
  QC_ASSERT (! types. contains (Token::eName));
  QC_ASSERT (types. size () <= 1);  // Not guaranteed ??
  QC_ASSERT (tokens. size () <= rows);
  if (tokensStored)
    { QC_ASSERT (types. size () <= tokens. size ()); }
  else
    QC_ASSERT (tokens. empty ());
  QC_IMPLY (types. contains (Token::eText), len_max);
  
  for (const auto& it : name2schema)
  {
    QC_ASSERT (isIdentifier (it. first));
    QC_ASSERT (it. second);
    it. second->qc ();
  }
}



void Schema::saveText (ostream &os) const 
{
  if (multiple)
    os << " TABLE";  
  for (const Token::Type type : types)
    os << ' ' << Token::type2str (type);
  if (len_max)
    os << " size:" << len_max;
  if (! tokens. empty ())
    os << " values:" << tokens. size ();
  os << " rows:" << rows;
  {
    const Offset ofs;
    for (const auto& it : name2schema)
    {
      Offset::newLn (os);
      os << it. first;
      it. second->saveText (os);
    }
  }
}



void Schema::merge (Schema& other)
{
  ASSERT (this != & other);
  
  for (auto& it : other. name2schema)
    if (const Schema* sch = findPtr (name2schema, it. first))
      var_cast (sch) -> merge (* var_cast (it. second));
    else
    {
      name2schema [it. first] = it. second;
      it. second = nullptr;
    }
    
  if (other. multiple)
    multiple = true;
  rows += other. rows;
  tokens << other. tokens;
  maximize (len_max, other. len_max);

  // types  
  types << other. types;   
  if (   types. contains (Token::eText)
      && types. size () >= 2
     )
  {
    types. clear ();
    types << Token::eText;
  } 
  else if (   types. contains (Token::eDouble)
           && types. contains (Token::eInteger)
          )
    types. erase (Token::eInteger);
}




// Data

Data::Data (TokenInput &ti)
{
  ti. get ('<');
  ti. get ('?');
  ti. get ("xml");
  for (;;)
  {
    const Token t (ti. get ());
    if (t. isDelimiter ('>'))
      break;
    if (t. empty ())
      throw runtime_error ("XML header is not finished");
  }
  
  readInput (ti);
}



void Data::readInput (TokenInput &ti)
{
  ASSERT (children. empty ());
  ASSERT (token. empty ());
  
  // name, isEnd
  ti. get ('<');
  {
    Token nameT (ti. get ());
    if (nameT. isDelimiter ('/'))
    {
      isEnd = true;
      nameT = move (ti. get ());
    }
    if (nameT. type != Token::eName)
      ti. error ("Name");
    name = move (nameT. name);
  }
  
  // children
  bool finished = false;
  for (;;)
  {
    Token attr (ti. get ());
    if (attr. isDelimiter ('/'))
    {
      finished = true;
      ti. get ('>');
      break;
    }
    if (attr. isDelimiter ('>'))
      break;
    if (attr. type != Token::eName)
      ti. error ("Name");
    Token eq (ti. get ());
    bool xmlAttr = false;
    if (! eq. isDelimiter ('='))
    {
      QC_ASSERT (eq. isDelimiter (':'));
      QC_ASSERT (isLeft (attr. name, "xml"));
      Token attr1 (ti. get ());
      QC_ASSERT (attr1. type == Token::eName);
      attr. name += ":" + attr1. name;
      ti. get ('=');
      xmlAttr = true;      
    }
    Token value (ti. get ());  
    if (value. type != Token::eText)
      ti. error ("Text");
    trim (value. name);
    for (char& c : value. name)
      if ((uchar) c >= 127)
        c = '?';
    value. toNumberDate ();
    if (! xmlAttr)
      children << new Data (this, move (attr. name), move (value));
  }

  if (isEnd)
  {
    if (! children. empty ())
      ti. error ("Ending tag has children");
    return;
  }

  // text, children
  if (! finished)
  {
    if (ti. getNextChar () == '<')
      for (;;)
      {
        auto xml = new Data (this, ti);
        if (xml->isEnd)
        {
          delete xml;
          break;
        }
        children << xml;
      }
    else
    {
      token = move (ti. getXmlText ());
      token. toNumberDate ();
      for (char& c : token. name)
        if ((uchar) c >= 127)
          c = '?';
      ti. get ('/');
      const Token end (ti. get ());
      if (end. type != Token::eName)
        ti. error ("Name");
      if (end. name != name)
        ti. error ("Name " + strQuote (name));
      ti. get ('>');
    }
  }
}



void Data::qc () const
{
  if (! qc_on)
    return;
  Named::qc ();
   
  try 
  {
    QC_ASSERT (isIdentifier (name));
    QC_ASSERT (! isEnd);
  //QC_ASSERT (! children. empty () || ! text. empty ());
    QC_IMPLY (! token. name. empty (), goodName (token. name));
    for (const Data* child : children)
    {
      QC_ASSERT (child);
      child->qc ();
    }
  }
  catch (const exception &e)
  {
    saveText (cout);
    errorExit (e. what ());
  }
}



void Data::saveText (ostream &os) const 
{
  Offset::newLn (os);
  os << '<' << name << '>';
  
  if (children. empty ())
  {
    if (! token. empty ())
    {
      os << ": ";
      token. saveText (os);
    }
  }
  else
  {  
    {
      const Offset ofs;
      for (const Data* child : children)
        child->saveText (os);
      if (! token. empty ())
      {
        Offset::newLn (os);
        token. saveText (os);
      }
    }
    Offset::newLn (os);
    os << "</" << name << '>';
  }
}



Schema* Data::getSchema (bool storeTokens) const
{
  auto sch = new Schema (storeTokens);
  
  for (const Data* child : children)
  {
    Schema* other = child->getSchema (storeTokens);
    if (const Schema* childSch = findPtr (sch->name2schema, child->name))
    {
      var_cast (childSch) -> multiple = true;
      var_cast (childSch) -> merge (*other);
      delete other;
    }
    else
      sch->name2schema [child->name] = other;
  }

  if (! token. empty ())
  {
    sch->types << token. type;
    if (storeTokens)
      sch->tokens << token;
    if (token. type == Token::eText)
      sch->len_max = token. name. size ();
  }
  
  return sch;
}



}
