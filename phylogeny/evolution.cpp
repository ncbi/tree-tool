// evolution.cpp

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


