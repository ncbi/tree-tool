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

#include "../common.hpp"
using namespace Common_sp;



namespace Xml
{



struct Schema : Root
{
  // Tree
  map<string/*name*/,const Schema*> name2schema;
    // !nullptr
  bool multiple {false};
  Set<Token::Type> types;
  bool tokensStored {false};
  Set<Token> tokens;
    // !tokensStored => empty()
  size_t len_max {0};
  size_t rows {1};
  
  
  explicit Schema (bool tokensStored_arg)
    : tokensStored (tokensStored_arg)
    {}
  static Schema* readSchema (const string &fName,
                             string &name);
    // Input: fName: created by xml2schema
    // Output: name
private:
  Schema (const StringVector &vec,
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
  
  
  void merge (Schema& other);
};



struct Data : Named
{
  // Tree
  const Data* parent {nullptr};
  VectorOwn<Data> children;
    // May be empty()

  Token token;
    // May be empty()
private:
  bool isEnd {false};
public:
  

  explicit Data (TokenInput &ti);
private:
  Data (Data* parent_arg,
        TokenInput &ti)
    : parent (parent_arg)
    { readInput (ti); }
  Data (Data* parent_arg,
        string &&attr,
        Token &&value)
    : Named (move (attr))
    , parent (parent_arg)
    , token (move (value))
    {}
  void readInput (TokenInput &ti);
public:
  void qc () const override;
  void saveText (ostream &os) const override;
  
  
  Schema* getSchema (bool storeTokens) const;
    // Return: new
};




}  



#endif
