// genbank_text.hpp

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
*   GenBank format parsing
*
*/


#ifndef GenBank_Text_HPP
#define GenBank_Text_HPP


#include "../common.hpp"
using namespace Common_sp;



namespace Asn_sp
{


struct GenbankText : Root
{
  struct Record : Named
  {
  	bool isFeature {false};
  private:
  	string value;
  public:

  	Record (const string &name_arg,
  	        bool isFeature_arg);
  	void saveText (ostream &os) const
  	  { os << name << ": " << getValue () << endl; }

  	void addValue (const string &s);
  	string getValue () const;
  };	
  Vector<Record> records;
    // "source" feature is parsed
  string assembler;
  string sequencer;
  double genomeCoverage {0.0};
  

  explicit GenbankText (LineInput &f);
  void saveText (ostream &os) const;
  
  
  string name2value (const string &name,
                     bool isUniq = true) const;
    // Output: May be empty()
    // Requires: if isUniq then name occurs at most once in records
  string keyword2name (const string &keyword,
                       bool caseSensitive,
                       string &content) const;
    // Return: string() <=> not found
    // Output: content
};



}



#endif


