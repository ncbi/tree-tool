// evolution.cpp

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
*   Evolution utilities
*
*/


#undef NDEBUG
#include "../common.inc"

#include "evolution.hpp"



namespace DistTree_sp
{


Real intersection2dissim (Real size1,
                          Real size2,
                          Real intersection,
                          Real intersection_min,
	                        Prob sizes_ratio_min,
	                        bool ave_arithP)
{
	ASSERT (size1 >= 0);
	ASSERT (size2 >= 0);
  ASSERT (intersection <= size1);
  ASSERT (intersection <= size2);
  ASSERT (isProb (sizes_ratio_min));
  
  
  if (verbose ())
    cout << '\t' << intersection 
         << '\t' << size1 
         << '\t' << size2 
         << endl;
         
  if (nullReal (max (size1, size2)))
  	return NaN;

  if (  min (size1, size2) 
  	  / max (size1, size2)
  	    < sizes_ratio_min
  	 )
  	return NaN;
         
  if (intersection < intersection_min)
  	return NaN;
  	
  const Real dissim1 = log (size1 / intersection);
  const Real dissim2 = log (size2 / intersection);
  ASSERT (dissim1 >= 0);
  ASSERT (dissim2 >= 0);
  
  if (ave_arithP)
    return ave_arith (dissim1, dissim2);
  //return ave_geom  (dissim1, dissim2);
  else
    return ave_harm  (dissim1, dissim2);
}




// Hashes

Hashes::Hashes (const string &fName)
{
  reserve (10000);  // PAR

  LineInput hf (fName);
  size_t prev = 0;
  while (hf. nextLine ())
  {
    const size_t hash = str2<size_t> (hf. line);
    ASSERT (hash);
    if (hash <= prev)
      throw runtime_error ("Hash " + hf. line + " is not greater than the previous hash " + toString (prev));
    *this << hash;
    prev = hash;
  } 
  searchSorted = true;
/*
  if (empty ())
    throw runtime_error ("No hashes for " + fName);
*/
}



// Read Hashes from a binary file ??




}


