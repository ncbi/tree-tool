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
*   ASCII XML analysis utilities
*
*/


#ifndef XML_HPP_78429
#define XML_HPP_78429

#include "../common.hpp"
#include "../tsv/tsv.hpp"
using namespace Common_sp;




namespace Xml_sp
{

// ASCII XML



// Analysis

struct FlatTable : Root
// name: file name
{
  OFStream f;
  size_t keys {0};
  VectorPtr<Token> row;
  Token id;
  
  
  FlatTable (const string &fileName,
             size_t keys_arg,
             size_t column_max)
    : f (fileName)
    , keys (keys_arg)
    , row (column_max)
    { row [0] = & id; 
      id. type = Token::eInteger;
    }
  void qc () const override;
    
    
  void write (size_t xmlNum);
    // Update: row[]
};



struct Schema : Root
{
  // Tree
  const Schema* parent {nullptr};
  map<string/*name*/,const Schema*> name2schema;
    // !nullptr
    // 1-1
  bool multiple {false};
  Set<Token::Type> types;
  bool tokensStored {false};
  Set<Token> tokens;
    // !tokensStored => empty()
  size_t len_max {0};
  size_t rows {1};
  
  unique_ptr<const FlatTable> flatTable;
    // get() => multiple
  size_t column {no_index};
    // Valid if flatTable.get()
    // no_index <=> types.empty()
  
  
  Schema (const Schema* parent_arg,
          bool tokensStored_arg)
    : parent (parent_arg)
    , tokensStored (tokensStored_arg)
    {}
  static Schema* readSchema (const string &fName,
                             string &name);
    // Input: fName: created by xml2schema
    // Output: name
private:
  Schema (const Schema* parent_arg,
          const StringVector &vec,
          size_t &start,
          size_t depth,
          string &name);
public:
 ~Schema ()
    { for (auto& it : name2schema)
        delete it. second;
    }
  void qc () const override;
  void saveText (ostream &os) const override;
  
  
  const string& schema2name (const Schema* sch) const
    { for (const auto& it : name2schema)
        if (it. second == sch)
          return it. first;
      throw logic_error ("Schema::schema2name()");
    }
  void merge (Schema& other);

  // SQL      
  string getColumnName (const Schema* rootTable) const
    { if (this == rootTable)
        return noString;
      if (! parent)
        return noString;
      string s (parent->getColumnName (rootTable));
      if (! s. empty ())
        s += "_";
      return s + parent->schema2name (this);
    }
  void printTableDdl (ostream &os,
                      bool dataP,
                      bool indexP,
                      const string &sqlSchema) const
    { printTableDdl (os, dataP, indexP, sqlSchema, nullptr); 
      os << "go" << endl;
    }
private:
  void printTableDdl (ostream &os,
                      bool dataP,
                      bool indexP,
                      const string &sqlSchema,
                      const Schema* curTable) const;
  void printColumnDdl (ostream &os,
                       const Schema* curTable) const;
public:
  
  // Output: flatTable, column
  void setFlatTables (const string &dirName,
                      const Schema* curTable);
private:
  void setFlatColumns (const Schema* curTable,
                       size_t &column_max);
public:
};



struct Data : VirtNamed  
{
	const Names &names;  
	size_t nameIndex {no_index};
	  // < names.size()
	
  // Tree
  const Data* parent {nullptr};
  VectorOwn<Data> children;
    // May be empty()

  const bool binary {false};
  const bool attribute {false};
  bool colonInName {false};
  Token token;
    // May be empty()
    // UTF-8
private:
  bool isEnd {false};
  bool xmlText {false};
  bool merged {false};
public:
	
	// Updatable
	uchar searchFound {0};
	uchar searchFoundAll {0};
	  // searchFound => searchFoundAll
	  // searchFoundAll => parent->searchFoundAll
  

private:
  Data (bool headerP,
        Names &names_arg,
        TokenInput &tin,
        VectorOwn<Data> &markupDeclarations);
    // Text
  Data (Names &names_arg,
        Data* parent_arg,
        TokenInput &ti)
    : names (names_arg)
    , parent (parent_arg)
    { readInput (ti); }
  Data (Names &names_arg,
        Data* parent_arg,
        bool attribute_arg,
        bool colonInName_arg,
        const string &attr,
        Token &&value)
    : names (names_arg)
    , nameIndex (names_arg. add (attr))
    , parent (parent_arg)
    , attribute (attribute_arg)
    , colonInName (colonInName_arg)
    , token (std::move (value))
    { if (attribute)
        token. quote = '\0';
    }
  Data (const Names &names_arg,
        ifstream &f,
        bool first);
    // Binary
  Data (const Data &) = default;
  void readInput (TokenInput &ti);
  static bool readColonName (TokenInput &ti,
                             string &name);
    // Update: name
  Data (Names &names_arg,
        Data* parent_arg,
        const string &name_arg,
        const Json& json);
    // name_arg: '-' are replaced by '_'
public:
  static Data* load (bool headerP,
                     Names &names,
                     const string &fName,
                     VectorOwn<Data> &markupDeclarations);
    // Text XML
    // Data ::=   <tag attribute1="value1" attribute2="value2" ... />
    //          | <tag attribute1="value1" attribute2="value2" ... > Data* XmlText </tag>
    //          | <!-- comment -->
    //          | <? ProcessingInstruction ?>
  static Data* load (Names &names,
                     const string &fName);
    // Binary XML, see common.hpp
  Data (Names &names_arg,
        const JsonMap& json)
    : Data (names_arg, nullptr, "Json", json)
    {}
  void clear () override;
  void qc () const override;
  void saveXml (Xml::File &f) const override;
  string getName () const final
    { return nameIndex == no_index ? noString : names. num2elem [nameIndex]; }
  
  
  size_t getDepth () const
    { if (parent)
        return parent->getDepth () + 1;
      return 0;
    }
  size_t getNodes () const
    { size_t n = 1;
      for (const Data* child : children)
        n += child->getNodes ();
      return n;
    }
  size_t getLeaves () const
    { if (children. empty ())
        return 1;
      size_t n = 0;
      for (const Data* child : children)
        n += child->getLeaves ();
      return n;
    }
  StringVector getText () const;
  bool hasDescendantName (const string &name_arg) const
    { if (getName () == name_arg)
        return true;
      for (const Data* child : children)
        if (child->hasDescendantName (name_arg))
          return true;
      return false;
    }
  bool contains (const string &needle,
                 bool targetNameP,
                 StringMatch::Type matchType) const
    { const string s (targetNameP ? getName () : token. name);
    	return matches (s, needle, matchType); 
    }
  bool find (VectorPtr<Data> &path,
             const string &needle,
             bool targetNameP,
             StringMatch::Type matchType) const;
    // DFS
    // Input; !targetNameP <=> target is token
    // Update: path: reverse path to *Data containing <needle>
    //               valid if Return
    // Invokes: contains()
  size_t setSearchFound (const string &needle,
				                 bool targetNameP,
						             StringMatch::Type matchType,
					               uchar mask);
    // Return: number of Data's
    // Output: searchFound, searchFoundAll
    // Invokes: contains()
  void unsetSearchFound (uchar mask)
		{ const uchar b = reverse (mask); 
		  searchFound    &= b; 
			searchFoundAll &= b;
			for (const Data* child : children)
		    var_cast (child) -> unsetSearchFound (mask);
		}
  const Data* name2child (const string &name_arg) const
    { for (const Data* child : children)
        if (child->getName () == name_arg)
          return child;
      return nullptr;
    }
  TextTable unify (const Data& query,
                   const string &rowTagName,
                   const string &variableTagName) const;
    // Unification of query and *this by ancestors' tag names
    // Return: Text of *this unifying with "<" variableTagName ">" column_name "</" variableTagName ">" in query --> column column_name in TextTable
    //         Text is in HTML
    // Requires: !query.parent
private:
  StringVector tagName2texts (const string &tagName) const;
  bool unify_ (const Data& query,
               const string &rowTagName,
               const string &variableTagName,
               map<string,StringVector> &tag2values,
               TextTable &tt) const;
    // Update: tag2values, tt (append)
public:
  Schema* createSchema (bool storeTokens) const;
    // Return: new, !nullptr
    //         Data where colonInName are skipped
    // XML tags which are not descendatnts of any table tags are ignored
  void writeFiles (size_t xmlNum,
                   const Schema* sch,
                   FlatTable* flatTable) const;
  void text2token (const string &tagName);
  bool chars2token (const string &tagName,
                    size_t len_max);
  void mergeSingleChildren ();
  void visualizeTokenTrailingSpaces ();
};




}  



#endif
