// gzip.cpp

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


#undef NDEBUG

#include "gzip.hpp"

#include "common.inc"



namespace Common_sp
{
 

// GZip

bool GZip::nextLine ()
{
  line. clear ();      
  
  ASSERT (start <= bufferSize_real);
  if (start == bufferSize_real)
    return false;

  for (;;)
  {
    QC_ASSERT (bufferSize_real <= bufferSize);
    ASSERT (start < bufferSize_real);
    size_t stop = start;
    while (stop < bufferSize_real && buffer [stop] != '\n')
      stop++;
    if (stop < bufferSize_real)
    {
      ASSERT (buffer [stop] == '\n');
      buffer [stop] = '\0';
      line += & buffer [start];
      start = stop + 1;
      prog ();
      return true;
    }
    else
    {
      line += & buffer [start];
      if (bufferSize_real < bufferSize)  // Last EOL is missing
      {
        QC_ASSERT (gzeof (f));
        start = bufferSize_real;
        prog ();
        return true;
      }
      start = 0;
      read ();
    }
  }
}



void GZip::read ()
{
  const int bytesRead = gzread (f, buffer, bufferSize);
  if (bytesRead < 0)
    throw runtime_error ("bytesread = " + to_string (bytesRead));
  bufferSize_real = (size_t) bytesRead;
}



}
