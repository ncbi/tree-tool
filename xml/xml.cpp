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



namespace Xml_sp
{



// FlatTable

void FlatTable::qc () const
{
  if (! qc_on)
    return;
  Root::qc ();
    
  QC_ASSERT (keys);
  QC_ASSERT (keys <= row. size ());
  QC_ASSERT (keys <= 2);

  FOR (size_t, i, keys)
    QC_ASSERT (row [i]);  
}



void FlatTable::write (size_t xmlNum)
{   
  FOR (size_t, i, keys)
  {
    QC_ASSERT (! row [i] -> empty ());
    QC_ASSERT (row [i] -> n);  
  }

  f << xmlNum;
  FFOR (size_t, i, row. size ())
  { 
    f << '\t';
    if (const Token* token = row [i])
      if (! token->name. empty ())
      {
        if (! goodName (token->name))
          throw runtime_error ("Bad value: " + strQuote (token->name));
        f << token->name;
      }
    if (i >= keys)
      row [i] = nullptr;
  }      
  f << endl;
}




// Schema

Schema* Schema::readSchema (const string &fName,
                            string &name)
{
  const StringVector vec (fName, (size_t) 100);  // PAR
  QC_ASSERT (! vec. empty ());
  
  size_t start = 0;
  return new Schema (nullptr, vec, start, 0, name);
}



Schema::Schema (const Schema* parent_arg,
                const StringVector &vec,
                size_t &start,
                size_t depth,
                string &name)
: parent (parent_arg)
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
  //if (type == Token::eText)
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
    auto sch = new Schema (this, vec, start, depth, name1);
    QC_ASSERT (! name1. empty ());
    name2schema [name1] = sch;
  }
}



void Schema::qc () const 
{
  if (! qc_on)
    return;
  Root::qc ();
    
  QC_IMPLY (! parent, ! multiple);
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
    schema2name (it. second);
    it. second->qc ();
  }
  
  if (flatTable. get ())
  {
    flatTable->qc ();
    QC_ASSERT (multiple);
    QC_ASSERT ((column == no_index) == types. empty ());
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

  if (column != no_index)
    os << " column:" << column;
    
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



void Schema::printTableDdl (ostream &os,
                            const Schema* curTable) const
{ 
  IMPLY (curTable, curTable->multiple);
  
  if (multiple)
  {
    const string table (getColumnName (nullptr));
    // "xml_num_", "id_" ??
    os << "drop   table " + table << ";" << endl 
       << "create table " + table << endl << '(' << endl
       << "  xml_num_ int  not null" << endl
       << ", id_ bigint  not null" << endl;  
    const string refTable (curTable ? curTable->getColumnName (nullptr) : "");
    if (! refTable. empty ())
      os << ", " << refTable << "_id_ bigint  not null" << endl;  
    curTable = this;
    printColumnDdl (os, curTable);  
    os << ");" << endl;
    os << "alter table " << table << " add constraint " << table << "_pk primary key (xml_num_, id_);" << endl;
    os << "grant select on " << table << " to public;" << endl;
    if (! refTable. empty ())
    {
      os << "create index " << table << "_idx on " << table << "(xml_num_," << refTable << "_id_);" << endl;
      os << "alter table " << table << " add constraint " << table << "_fk foreign key (xml_num_, " << refTable << "_id_) references " << refTable << "(xml_num_,id_);" << endl;
    }
    os << endl << endl;
  }

  for (const auto& it : name2schema)
    it. second->printTableDdl (os, curTable);
}



void Schema::printColumnDdl (ostream &os,
                             const Schema* curTable) const
{  
  ASSERT (curTable);
  ASSERT (curTable->multiple);
  
  if (multiple && curTable != this)
    return;

  if (! types. empty ())
  {
    ASSERT (types. size () == 1);
    const Token::Type type = types. front ();
    string sqlType;
    switch (type)
    {
      case Token::eText:     sqlType = "varchar(" + (len_max > 1000 ? "max" : to_string (len_max)) + ")"; break;
      case Token::eInteger:  sqlType = "bigint"; break;
      case Token::eDouble:   sqlType = "float"; break;
      case Token::eDateTime: sqlType = "datetime"; break;
      default: throw runtime_error ("UNknown SQL type");
    }
    ASSERT (! sqlType. empty ());
    string colName (getColumnName (curTable));
    if (colName. empty ())
      colName = "val_";  // ??
    os << ", " << colName << ' ' << sqlType << endl;
  }

  for (const auto& it : name2schema)  
    it. second->printColumnDdl (os, curTable);
}



void Schema::setFlatTables (const string &dirName,
                            const Schema* curTable) 
{ 
  ASSERT (isDirName (dirName));
  IMPLY (curTable, curTable->multiple);
  
  if (multiple)
  {
    const Schema* parentTable = curTable;
    curTable = this;
    const size_t keys = 1/*pk*/ + (bool) parentTable/*fk*/;
    size_t column_max = keys;
    setFlatColumns (curTable, column_max);  
    ASSERT (! flatTable. get ());
    flatTable. reset (new FlatTable (dirName + getColumnName (nullptr), keys, column_max));
    if (parentTable)
    {
      const FlatTable* ft = parentTable->flatTable. get ();
      ASSERT (ft);
      var_cast (flatTable. get ()) -> row [1] = & ft->id;
    }
  }

  for (const auto& it : name2schema)
    var_cast (it. second) -> setFlatTables (dirName, curTable);
}



void Schema::setFlatColumns (const Schema* curTable,
                             size_t &column_max) 
{  
  ASSERT (curTable);
  ASSERT (curTable->multiple);
  ASSERT (column == no_index);
  QC_ASSERT (column_max < no_index);
  
  if (multiple && curTable != this)
    return;

  if (! types. empty ())
  {
    column = column_max;
    column_max++;
  }

  for (const auto& it : name2schema)  
    var_cast (it. second) -> setFlatColumns (curTable, column_max);
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
      string& s = token. name;
      FOR_REV (size_t, i, s. size ())
      {
        const char c = s [i];
        if (   ! printable (c)
            || c == '|'  // Delimiter in bulk-insert
           )
          s = s. substr (0, i) + "%" + uchar2hex ((uchar) c) + s. substr (i + 1);
      }
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
    token. qc ();
    QC_ASSERT (! contains (token. name, '<'));
    QC_ASSERT (! contains (token. name, '>'));
      // Other characters prohibited in XML ??
    for (const Data* child : children)
    {
      QC_ASSERT (child);
      child->qc ();
    }
  }
  catch (const exception &e)
  {
    Xml::File f (cout, "XML");
    saveXml (f);
    errorExit (e. what ());
  }
}



void Data::saveXml (Xml::File &f) const
{
  const Xml::Tag tag (f, name);
  for (const Data* child : children)
    child->saveXml (f);
  f << token. str ();
}



bool Data::find (VectorPtr<Data> &path,
                 const string &what,
                 bool equalName,
                 bool tokenSubstr,
                 bool tokenWord) const
{
  if (   ! equalName
      && ! tokenSubstr
      && ! tokenWord
     )
    return false;
    
  if (equalName && name == what)
    return true;
  if (tokenSubstr || tokenWord)
  {
    const string s (token. str ());
    if (tokenSubstr && contains (s, what))
      return true;
    if (tokenWord && containsWord (s, what))
      return true;
  }
    
  for (const Data* child : children)
    if (child->find (path, what, equalName, tokenSubstr, tokenWord))
    {
      path << child;
      return true;
    }
    
  return false;
}



Schema* Data::createSchema (bool storeTokens) const
{
  auto sch = new Schema (nullptr, storeTokens);
  
  for (const Data* child : children)
  {
    Schema* other = child->createSchema (storeTokens);
    other->parent = sch;
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
  //if (token. type == Token::eText)
      sch->len_max = token. name. size ();
  }
  
  return sch;
}



void Data::writeFiles (size_t xmlNum,
                       const Schema* sch,
                       FlatTable* flatTable) const
{
  QC_ASSERT (sch);
  ASSERT (sch->multiple == (bool) sch->flatTable. get ());

  bool newFlatTable = false;  
  if (const FlatTable* flatTable_ = sch->flatTable. get ())
  {
    flatTable = var_cast (flatTable_);
    newFlatTable = true;
    flatTable->id. n ++;
    flatTable->id. name = to_string (flatTable->id. n);
  }
    
  if (   ! token. empty ()
      && flatTable
     )
  {
    QC_ASSERT (sch->column < no_index);
    QC_ASSERT (! flatTable->row [sch->column]);
    flatTable->row [sch->column] = & token;
  }
  
  for (const Data* child : children)
    if (const Schema* childSch = findPtr (sch->name2schema, child->name))
      child->writeFiles (xmlNum, childSch, flatTable);
    else
      throw runtime_error ("Schema-XML mismatch: " + child->name);
    
  if (newFlatTable)
    flatTable->write (xmlNum);
}


}
