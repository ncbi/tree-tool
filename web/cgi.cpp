// cgi.cpp

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
*   CGI utilities
*
*/


#undef NDEBUG

#include "cgi.hpp"

#include "../common.hpp"
using namespace Common_sp;

#include "../common.inc"



namespace Cgi_sp
{



KeyValue unCgi ()
{  
  bool post = false;
  {
    const string requestMethod (getEnv ("REQUEST_METHOD"));
    QC_ASSERT (! requestMethod. empty ());
    if      (requestMethod == "GET")
      post = false;
    else if (requestMethod == "POST")
      post = true;
    else
      throw runtime_error (FUNC "Unknown REQUEST_METHOD " + strQuote (requestMethod));
  }

  const string queryString (getEnv ("QUERY_STRING"));

  size_t len = 0;
  if (post)
  {
    const string lenS (getEnv ("CONTENT_LENGTH"));
    if (verbose ())
      cout << lenS << endl;
    len = str2<size_t> (lenS);
    if (verbose ())
      cout << len << endl;
  }
  else
    len = queryString. size ();
  

  KeyValue kv;
  bool isKey = true;
  string key;
  string value;
  bool special = false;
  string specialS;
  FOR (size_t, i, len)
  {
    char c = '\0';
    if (post)
    {
      if (! getChar (cin, c))
        throw runtime_error (FUNC "Cannot get character " + to_string (i + 1));
    }
    else
      c = queryString [i];
    if (verbose ())
      cout << i << ' ' << c << endl;
    QC_ASSERT (c);
  //QC_ASSERT (c != '\'' && c != '\"');

    if (special)
    {
      QC_ASSERT (! isKey);
      ASSERT (specialS. size () <= 2);
      if (! isHex (c))
        throw runtime_error (FUNC "hexadecimial digit exoetced");
      specialS += c;
      if (specialS. size () == 2)
      {
        c = char (hex2uchar (specialS [0]) * 16 + hex2uchar (specialS [1]));
        special = false;
        specialS. clear ();
        if (c != '\x0D')  // Ignore Windows EOL 1st character  // Check that it is followed by '\n' ??
          value += c;
      }
    }
    else
      switch (c)
      {
        case '&':
          QC_ASSERT (! isKey);
          QC_ASSERT (! key. empty ());
          if (contains (kv, key))
            throw runtime_error (FUNC "Duplicate key " + strQuote (key));
          kv [key] = value;
          isKey = true;
          key. clear ();
          value. clear ();
          break;

        case '+':
          ASSERT (! isKey);
          value += ' ';
          break;

        case '=':
          if (isKey)
            isKey = false;
          else
            value += c;
          break;

        case '%':
          ASSERT (! isKey);
          ASSERT (specialS. empty ());
          special = true;
          break;

        default:
          if (isKey)
            key += c;
          else
            value += c;
      }
  }
  QC_ASSERT (! special);
  IMPLY (! kv. empty (), ! isKey);
  if (! key. empty ())
  {
    if (contains (kv, key))
      throw runtime_error (FUNC "Duplicate key " + strQuote (key));
    kv [key] = value;
  }

  
  return kv;
}



}  // namespace

