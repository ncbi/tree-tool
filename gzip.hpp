// gzip.hpp

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
*   GZip utilities
*
*/


#ifndef GZIP_HPP_769476  // random number
#define GZIP_HPP_769476


#include "common.hpp"

extern "C"
{
  #include <zlib.h>
    // Linking requires: -lz 
}



namespace Common_sp
{
  
  
struct GZip final : Root, Nocopy
{
  string line;
private:
  gzFile f;
  static constexpr size_t bufferSize {1024 * 1024};  // PAR
  char buffer [bufferSize + 1];
  size_t bufferSize_real {0};
    // <= bufferSize
  size_t start {0};
    // < bufferSize_real
	Progress prog;  
public:
  
  
  explicit GZip (const string &fName,
                 uint displayPeriod = 0)
    : f (gzopen (fName. c_str (), "rb"))
    , prog (0, displayPeriod)
    { if (! f)
        throw runtime_error ("Cannot open " + strQuote (fName));
      buffer [bufferSize] = '\0';
      read ();
    }
    // Input: fName: gzip'ped UNIX text file
 ~GZip ()
   { gzclose (f); }
  
    
  bool nextLine ();
private:
  void read ();
    // Output: buffer[], bufferSize_real
};



}  
  
  

#endif
