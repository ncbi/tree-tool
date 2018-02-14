// evolution.hpp

#ifndef EVOLUTION_HPP
#define EVOLUTION_HPP

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
using namespace DM_sp;



namespace DistTree_sp
{


struct Hashes : Vector<size_t>
{
	explicit Hashes (const string &fName);
	Hashes ()
	  {}
	
	Real getDissim (const Hashes &other,
	                size_t intersection_min) const;
	  // Symmetric
	  // Return: >= 0; may be NAN
};



}



#endif


