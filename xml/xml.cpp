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

#include "xml.hpp"

#include "../common.inc"



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
      //if (! goodName (token->name))
        //throw runtime_error ("Bad value: " + strQuote (token->name));
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
  const StringVector vec (fName, (size_t) 100, false);  // PAR
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
  TokenInput ti (iss, '\0', true);
  
  // name
  Token nameT (ti. get ());
  QC_ASSERT (nameT. type == Token::eName);
  name = std::move (nameT. name);
  
  // multiple
  Token tableT (ti. get ());
  QC_ASSERT (tableT. type == Token::eName);
  if (tableT. name  == "TABLE")
    multiple = true;
  else
    ti. setLast (std::move (tableT));
  
  // types
  Token typeT (ti. get ());
  QC_ASSERT (typeT. type == Token::eName);
  if (typeT. name == "rows")
    ti. setLast (std::move (typeT));
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
    QC_IMPLY (! contains (it. first, ':'), isIdentifier (it. first, true));
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
      os << "alter table " << table << " add constraint " << table << "_fk foreign key (xml_num_, " << refTable << "_id_) references " << refTable << "(xml_num_,id_) on delete cascade;" << endl;
      os << "-- delete /*top (1000000)*/ from " << table << " with (tablock) where not exists (select null from " << refTable << " (nolock) where " << refTable << ".xml_num_ = " << table << ".xml_num_ and " << refTable << ".id_ = " << table << "." << refTable << "_id_);" << endl;
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
      default: throw runtime_error ("Unknown SQL type");
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

Data::Data (Names &names_arg,
	          TokenInput &ti,
            VectorOwn<Data> &markupDeclarations)
: names (names_arg)
, binary (false)
, attribute (false)
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
  
  for (;;)
  {
    readInput (ti);
    if (   getName (). empty () 
    	  || getName () [0] != '!')
      break;
    markupDeclarations << new Data (*this);
    clear ();
  }
}



void Data::readInput (TokenInput &ti)
{
  ASSERT (children. empty ());
  ASSERT (token. empty ());
  ASSERT (! attribute);
  ASSERT (! colonInName);
  ASSERT (! isEnd);
  ASSERT (! xmlText);
  
  if (ti. getNextChar (true) == '<')
    ti. get ('<');
  else
  {
    nameIndex = var_cast (names). add ("XmlText");   // *this will be deleted
    xmlText = true; 
    token = std::move (ti. getXmlText ());
    token. toNumberDate ();
    return;
  }

  // name, isEnd
  {
    Token nameT (ti. get ());
    if (nameT. isDelimiter ('/'))
    {
      isEnd = true;
      nameT = std::move (ti. get ());
    }
    else if (nameT. isDelimiter ('!'))        
    {
      if (ti. getNextChar (true) == '-')
      {
        nameIndex = var_cast (names). add ("!--");
        token = std::move (ti. getXmlComment ());
      }
      else
      {
        nameT = ti. get ();
        nameIndex = var_cast (names). add ("!" + nameT. name);
        token = std::move (ti. getXmlMarkupDeclaration ());        
      }
      return;
    }
    else if (nameT. isDelimiter ('?'))  
    {
      nameIndex = var_cast (names). add ("ProcessingInstruction");
      token = std::move (ti. getXmlProcessingInstruction ());
      return;
    }
    if (nameT. type != Token::eName)
      ti. error (nameT, "Name");
    nameIndex = var_cast (names). add (nameT. name);
    string name = getName ();
    if (readColonName (ti, name))
    {
    	nameIndex = var_cast (names). add (name);
      colonInName = true;
    }
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
      ti. error (attr, "Name");
    bool colonInName_arg = false;
    if (readColonName (ti, attr. name))
      colonInName_arg = true;
    ti. get ('=');
    Token value (ti. get ());  
    if (value. type != Token::eText)
      ti. error (value, "Text");
    trim (value. name);
    for (char& c : value. name)
      if ((uchar) c >= 127)
        c = '?';
    value. toNumberDate ();
    children << new Data (var_cast (names), this, true, colonInName_arg, attr. name, std::move (value));
  }

  if (isEnd)
  {
    if (! children. empty ())
      ti. error ("Ending tag has children");
    if (finished)
      ti. error ("Ending tag has two slashes");
    return;
  }
  
  if (finished)
    return;

  // text, children
  if (ti. getNextChar (true) == '<')
    for (;;)
    {
      auto xml = new Data (var_cast (names), this, ti);
      if (xml->isEnd)
      {
        if (xml->nameIndex != nameIndex)
          ti. error ("Tag " + strQuote (getName ()) + " is closed by tag " + strQuote (xml->getName ()), false);
        ASSERT (! xml->xmlText);
        delete xml;
        break;
      }
      if (xml->xmlText)
      {
        ASSERT (token. empty ());
        token = std::move (xml->token);
        delete xml;
        trim (token. name);
        ti. get ('/');
        Token end = std::move (ti. get ());
        if (end. type != Token::eName)
          ti. error (end, "Name");
        readColonName (ti, end. name);
        if (end. name != getName ())
          ti. error (end, "Name " + strQuote (getName ()));
        ti. get ('>');
        break;
      }
      children << xml;
    }
  else
  {
    token = std::move (ti. getXmlText ());
    token. toNumberDate ();
    ti. get ('/');
    Token end (ti. get ());
    if (end. type != Token::eName)
      ti. error (end, "Name");
    readColonName (ti, end. name);
    if (end. name != getName ())
      ti. error (end, "Name " + strQuote (getName ()));
    ti. get ('>');
  }
}



bool Data::readColonName (TokenInput &ti,
                          string &name) 
{
  ASSERT (! name. empty ());
  
  if (ti. getNextChar (true) != ':')
    return false;
  ti. get (':');
  
  const Token name1 (ti. get ());
  if (name1. type != Token::eName)
    ti. error (name1, "Name (2nd half)");
  if (name1. name. empty ())
    ti. error (name1, "Name (2nd half) is empty");

  name += ":" + name1. name;
  
  return true;
}



Data::Data (const Names &names_arg,
            ifstream &f,
            bool first)
: names (names_arg)
, binary (true)
, attribute (false)
{
	ASSERT (f. good ());
	ASSERT (! names. size ());
	
	char c = '\0';

	// nameIndex
	{
		EXEC_ASSERT (getChar (f, c));
		const uchar b1 = (uchar) c;
		EXEC_ASSERT (getChar (f, c));
		const uchar b2 = (uchar) c;
		nameIndex = (size_t) b1 * 256 + (size_t) b2;
		ASSERT (nameIndex != no_index);
  }
  if (first)
  {
  	ASSERT (! nameIndex);
  }
  else
		if (! nameIndex)
		{
			isEnd = true;
			return;
		}
	
	// children
	for (;;)
	{
		unique_ptr<Data> child (new Data (names, f, false));
		if (child->isEnd)
			break;
		child->parent = this;
	  children << child. release ();
	}
	
	// token
	for (;;)
	{
		EXEC_ASSERT (getChar (f, c));
		if (! c)
			break;
		token. type = Token::eText;
		if (nonPrintable (c))
			token. name += nonPrintable2str (c);
		else
		  token. name += c;
		xmlText = true;
	}
}



Data* Data::load (Names &names,
	                const string &fName,
                  VectorOwn<Data> &markupDeclarations)
{ 
  unique_ptr<Xml_sp::Data> f;
  { 
    TokenInput ti (fName, '\0', true, false, 10000);  // PAR 
    try 
    { 
		  Unverbose unv;
    	f. reset (new Xml_sp::Data (names, ti, markupDeclarations));	
    }
    catch (const TextPos::Error &e)
      { throw e; }
    catch (const exception &e)
      { ti. error (e. what (), false); }
    const Token t (ti. get ());
    if (! t. empty ())
      ti. error (t, "Token after the XML end: " + strQuote (t. str ()));
  }
  return f. release ();
}



Data* Data::load (Names &names,
	                const string &fName)
{
	checkFile (fName);
	ifstream f (fName, ios_base::binary | ios_base::in);
	QC_ASSERT (f. good ());
	
  unique_ptr<Xml_sp::Data> d;
  { 
	  Unverbose unv;
  	d. reset (new Xml_sp::Data (names, f, true));	
  }
  ASSERT (! names. size ())
  
  char c = '\0';
  string s;
  while (getChar (f, c))
  {
    QC_ASSERT (c);
    if (c == '\n')
    {
      QC_ASSERT (! s. empty ());
    	names. add (s);
    	s. clear ();
    }
    else
      s += c;
  }
  
  return d. release ();
}



void Data::clear ()
{
  ASSERT (! parent);
  ASSERT (children. empty ());
  ASSERT (! attribute);
  
  colonInName = false;
  token. clear ();
  isEnd = false;
  xmlText = false;
//columnTags = 0;
}



void Data::qc () const
{
  if (! qc_on)
    return;
  VirtNamed::qc ();

  try 
  {
    QC_ASSERT (nameIndex < names. size ());
    string name (getName ());
    QC_ASSERT (! name. empty ());
    if (binary)
    {
      QC_ASSERT (names. size () < 256 * 256);
      QC_ASSERT (! Common_sp::contains (name, '\0'));
      QC_ASSERT (! Common_sp::contains (name, '\n'));
    }
    else if (! isLeft (name, "!"))
    {
	  //QC_IMPLY (colonInName, attribute);
	    if (colonInName)
	  	  replace (name, ':', '_');
	  	if (merged)
	  	  replace (name, '/', '_');
	    QC_ASSERT (isIdentifier (name, true));
	  }
    QC_ASSERT (! isEnd);
    token. qc ();
    if (token. type != Token::eText)
    {
    	QC_ASSERT (! Common_sp::contains (token. name, '<'));
      QC_ASSERT (! Common_sp::contains (token. name, '>'));
      // Other characters prohibited in XML ??
    }
    QC_ASSERT (Common_sp::contains (searchFoundAll, searchFound));
    QC_IMPLY (parent, Common_sp::contains (parent->searchFoundAll, searchFoundAll));
    for (const Data* child : children)
    {
      QC_ASSERT (child);
      child->qc ();
    }
  }
  catch (const exception &e)
  {
    PRINT (getName ());
    PRINT (Token::type2str (token. type));
    PRINT (token);
    PRINT (token. name);
    errorExit (e. what ());
  }
}



void Data::saveXml (Xml::File &f) const
{
  const Xml::Tag tag (f, getName ());

  for (const Data* child : children)
    child->saveXml (f);
    
  if (! token. empty ())
    f. print (token. name);
}



string Data::getText () const 
{
  string s;
  
  if (attribute)
    return s;
  
  for (const Data* child : children)
    add (s, "; ", child->getText ());
    
  if (! token. empty ())
    add (s, "; ", token. str ());
    
  return s;
}



bool Data::find (VectorPtr<Data> &path,
                 const string &needle,
                 bool targetNameP,
		             StringMatch::Type matchType) const
{
  if (contains (needle, targetNameP, matchType))
    return true;
    
  for (const Data* child : children)
    if (child->find (path, needle, targetNameP, matchType))
    {
      path << child;
      return true;
    }
    
  return false;
}



size_t Data::setSearchFound (const string &needle,
			  		                 bool targetNameP,
								             StringMatch::Type matchType,
					                   uchar mask)
{
	ASSERT (mask);

  size_t n = 0;

  if (contains (needle, targetNameP, matchType))
  {
    searchFound    |= mask;
    searchFoundAll |= mask;
    n++;
    const Data* parent_ = parent;
    while (   parent_ 
           && ! Common_sp::contains (parent_->searchFoundAll, mask)
          )
    {
    	var_cast (parent_) -> searchFoundAll |= mask;
    	parent_ = parent_->parent;
    }
  }
    
  for (const Data* child : children)
    n += var_cast (child) -> setSearchFound (needle, targetNameP, matchType, mask);
    
  return n;
}



TextTable Data::unify (const Data& query,
                       const string &rowTagName,
                       const string &variableTagName) const
{
  ASSERT (! query. parent);
  
  StringVector header (query. tagName2texts (variableTagName));
  if (header. empty ())  
    throw runtime_error ("No variable tags");
    
  header. sort ();
  header. uniq ();

  TextTable tt;
  tt. pound = true;
  for (const string& s : header)
    tt. header << TextTable::Header (s);
  tt. qc ();
  
  map<string,StringVector> tag2values;
  unify_ (query, rowTagName, variableTagName, tag2values, tt);
  
  return tt;
}



StringVector Data::tagName2texts (const string &tagName) const
{ 
  ASSERT (! tagName. empty ());
  
  StringVector vec;
  if (getName () == tagName)
  {
    QC_ASSERT (parent);
    QC_ASSERT (children. empty ());
    QC_ASSERT (token. type == Token::eText);
    vec << token. name;
  }
  else
    for (const Data* child : children)
      vec << std::move (child->tagName2texts (tagName));

  return vec;
}



bool Data::unify_ (const Data& query,
                   const string &rowTagName,
                   const string &variableTagName,
                   map<string,StringVector> &tag2values,
                   TextTable &tt) const
{
  ASSERT (! rowTagName. empty ());
  ASSERT (! variableTagName. empty ());
  ASSERT (& names == & query. names);  
  ASSERT (query. getName () != variableTagName);
  
  if (nameIndex != query. nameIndex)
    return false;
  if (   ! query. token. empty () 
      && ! (query. token == token)
     )
    return false;
    
  for (const Data* queryChild : query. children)
    if (queryChild->getName () != variableTagName)
    {
      ASSERT (queryChild->parent);
      bool unified = false;
      for (const Data* dataChild : children)
        if (dataChild->unify_ (*queryChild, rowTagName, variableTagName, tag2values, tt))
          unified = true;
      if (   ! unified 
          && ! queryChild->hasDescendantName (variableTagName)  // Otherwise the column is missed
         )
        return false;
    }
  
  if (const Data* columnData = query. name2child (variableTagName))
  {
    const string text (getText ());
    if (! text. empty ())
      tag2values [columnData->token. name] << text;
  }
    
  if (   query. getName () == rowTagName
      && ! tag2values. empty ()
     )
  {
    StringVector row (tt. header. size ());
    for (const auto& it : tag2values)
      row [tt. col2num (it. first)] = it. second. toString ("; ");
    tt. rows << std::move (row);
    tag2values. clear ();
  }
  
  return true;
}



Schema* Data::createSchema (bool storeTokens) const
{
  QC_ASSERT (! colonInName);

  auto sch = new Schema (nullptr, storeTokens);
  
  for (const Data* child : children)
    if (! child->colonInName)
	  {
	    Schema* other = child->createSchema (storeTokens);
	    other->parent = sch;
	    if (const Schema* childSch = findPtr (sch->name2schema, child->getName ()))
	    {
	      var_cast (childSch) -> multiple = true;
	      var_cast (childSch) -> merge (*other);
	      delete other;
	    }
	    else
	      sch->name2schema [child->getName ()] = other;
	  }

  if (! token. empty ())
  {
    sch->types << token. type;
    if (storeTokens)
      sch->tokens << token;
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
  QC_ASSERT (! colonInName);


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
    if (! child->colonInName)
    {
      if (const Schema* childSch = findPtr (sch->name2schema, child->getName ()))
        child->writeFiles (xmlNum, childSch, flatTable);
      else
        throw runtime_error ("Schema-XML mismatch: " + child->getName ());
    }
    
  if (newFlatTable)
    flatTable->write (xmlNum);
}



void Data::text2token (const string &tagName)
{
  if (children. empty ())
    return;
    
  {
  	const Data* child = children. front ();  
  	ASSERT (child);
  	if (   token. empty ()
  		  && child->getName () == tagName
  		  && child->children. empty ()
  		//&& ! child->token. empty ()
  		 )
  	{
  		token = std::move (var_cast (child) -> token);
 		  token. type = Token::eText;
 		  ASSERT (! token. empty ());
  		children. erasePtr (0);
  	}
  }
  
	for (const Data* child : children)
	{
	  ASSERT (child);
    var_cast (child) -> text2token (tagName);
	}
}



bool Data::chars2token (const string &tagName,
	                      size_t len_max)
{
  if (children. empty ())
    return getName () == tagName;

  bool done = false;
  Vector<size_t> toDelete;
  const bool tokenEmpty = token. empty ();
	FFOR (size_t, i, children. size ())
	{
		const Data* child = children [i];
	  ASSERT (child);
    if (var_cast (child) -> chars2token (tagName, len_max))
    	if (tokenEmpty)
    	{
	    	if (len_max > token. name. size ())
	 		    token. name += child->token. name. substr (0, len_max - token. name. size ());
	 		  token. type = Token::eText;
	  		done = true;
	  		if (   child->children. empty ()
	  			  && child->getName () == tagName 
	  			 )
	  			toDelete << i;
	    }
	}
	
	IMPLY (! toDelete. empty (), done);
	FOR_REV (size_t, i, toDelete. size ())
		children. erasePtr (toDelete [i]);
	
  return done;
}



void Data::mergeSingleChildren ()
{
	ASSERT (! merged);
		
  while (children. size () == 1)
  {
  	ASSERT (! xmlText);
  	const Data* child = children. front ();  
  	ASSERT (child);
  	if (   ! token. empty () 
  		  && ! child->token. empty ()
  		  && ! (child->token == token)
  		 )
  		break;
  	nameIndex = var_cast (names). add (getName () + "/" + child->getName ());
  	if (child->colonInName)
  		colonInName = true;
  	if (token. empty ())
	    token = std::move (var_cast (child) -> token);
		children. clear ();
		children = std::move (var_cast (child) -> children);
		ASSERT (child->children. empty ());
		delete child;
		merged = true;
		for (const Data* child_ : children)
		{
		  ASSERT (child_);
	    var_cast (child_) -> parent = this;
		}
  }
  
	for (const Data* child : children)
	{
	  ASSERT (child);
    var_cast (child) -> mergeSingleChildren ();
	}
}



void Data::visualizeTokenTrailingSpaces ()
{
  visualizeTrailingSpaces (token. name);
	for (const Data* child : children)
	{
	  ASSERT (child);
    var_cast (child) -> visualizeTokenTrailingSpaces ();
	}
}



}
