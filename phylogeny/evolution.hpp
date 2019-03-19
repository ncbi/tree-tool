// evolution.hpp

#ifndef EVOLUTION_HPP
#define EVOLUTION_HPP

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
using namespace DM_sp;



namespace DistTree_sp
{


Real intersection2dissim (Real size1,
                          Real size2,
                          Real intersection,
                          Real intersection_min,
	                        Prob sizes_ratio_min,
	                        bool ave_arithP);
  // Input: ave_arithP; false <=> ave_harm()
  // Return: >= 0; may be NaN
  // Symmetric



struct Hashes : Vector<size_t>
{
	explicit Hashes (const string &fName);
	Hashes () = default;
	
	Real getDissim (const Hashes &other,
	                size_t intersection_min,
	                Prob hashes_ratio_min) const
		{ return intersection2dissim ( (Real) size ()
			                           , (Real) other. size ()
			                           , (Real) getIntersectSize (other)
			                           , (Real) intersection_min
			                           , hashes_ratio_min
			                           , false  // PAR
			                           ); 
		}
};



}



#endif


