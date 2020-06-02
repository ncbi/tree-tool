// asn.hpp

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
*   NCBI textual ASN.1 parsing
*
*/


#ifndef ASN_HPP
#define ASN_HPP


#include "common.hpp"
using namespace Common_sp;



namespace Asn_sp
{


struct Asn : Root
{
	CharInput in;
	  // ASN.1 text
  string title;    
  StringVector fDict;
    // Index: feature id

  // Current node
  constexpr static const uint no_id {numeric_limits<uint>::max ()};
  uint id {0};
  uint parent {no_id};
  StringVector features;
    // size() = fDict.size()

private:
	Token last;
public:
    

protected:
  explicit Asn (const string &fName)
		: in (fName, 100000)  // PAR
		{}
public:
		
		
  void asnRead ();
    // Invokes: processNode()
  virtual void processNode () = 0;
  bool expectTry (const string& text);
  void expect (const string& text)
    { if (! expectTry (text))
      	in. error (strQuote (text));
    }
	void printFeatures (ostream &os) const;
};



}



#endif


