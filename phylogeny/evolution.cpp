// evolution.cpp

#undef NDEBUG
#include "../common.inc"

#include "evolution.hpp"



namespace DistTree_sp
{


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
/*
  if (empty ())
    throw runtime_error ("No hashes for " + fName);
*/
}



// Read Hashes from a binary file ??



Real Hashes::getDissim (const Hashes &other,
	                      size_t intersection_min,
	                      Prob hashes_ratio_min) const
{
  const size_t intersection = getIntersectSize (other);
  ASSERT (intersection <=        size ());
  ASSERT (intersection <= other. size ());
  
  if (verbose ())
    cout << '\t' << intersection 
         << '\t' <<        size ()
         << '\t' << other. size ()
         << endl;
         
  if (  (Real) min (size (), other. size ()) 
  	  / (Real) max (size (), other. size ())
  	    < hashes_ratio_min
  	 )
  	return NAN;
         
  if (intersection < intersection_min)
  	return NAN;
         
  return - 0.5 * (  log ((Real) intersection / (Real)        size ()) 
                  + log ((Real) intersection / (Real) other. size ()) 
                 );
}



}


