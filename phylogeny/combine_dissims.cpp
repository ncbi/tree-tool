// combine_dissim.cpp

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
*   Combine dissimilarities
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
using namespace DM_sp;
#include "version.inc"



#undef WEIGHT 



namespace 
{
	
	
struct Scale
{
	Real unit {INF};
	Real raw_min {0.0};
	Real raw_max {0.0};
#ifdef WEIGHT
	Real weight {NaN};
#endif
	
	Scale (const string &line,
	       Real coeff)
	  {
	  	static Istringstream iss;
	  	iss. reset (line);
	  	iss >> unit >> raw_min >> raw_max 
        #ifdef WEIGHT
	  	    >> weight
	  	  #endif
	  	  ;
	  	unit *= coeff;
	  	ASSERT (unit > 0.0);
	  //ASSERT (unit <= 1.0); 
	  	ASSERT (raw_min >= 0.0);
	  	ASSERT (raw_min <= raw_max);
    #ifdef WEIGHT
	  	ASSERT (weight > 0.0);
	  #endif
	  }
	Scale () = default;
	  
	Real raw2dissim (Real raw) const
	  { if (isNan (raw))
	  	  return NaN;
	  	ASSERT (raw >= 0.0);
	  	if (raw < raw_min)
	  	  return NaN;
	  	if (raw > raw_max)
	  	  return NaN;
	  	return raw * unit; 
	  }
	Real dissim_max () const
	  { return raw2dissim (raw_max); }
};
	
	
	
struct ObjPair
{
	string name1;
	string name2;
	Real dissim {NaN};
	
	explicit ObjPair (const string &line)
	  { 
	  	static Istringstream iss;
	  	iss. reset (line);
	  	static string dissimS;
	  	dissimS. clear ();
	  	iss >> name1 >> name2 >> dissimS;
	  	ASSERT (iss. eof ());
	  	ASSERT (name1 < name2);
	  	dissim = str2real (dissimS);
	  	IMPLY (! isNan (dissim), dissim >= 0.0);
	  }
};




struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Combine different dissimilarities for the same pairs of objects")
    {
      version = VERSION;
  	  addPositional ("dissims", "File with lines: <obj1> <obj2> <dissimilarity>; # lines = # object pairs times # dissimilarities");
  	  addPositional ("scales", "File with equalizing scales and max. values for each dissimilarity, ordered by scale descending.\n\
Line format: <unit> <raw_min> <raw_max>");
      addKey ("coeff", "Coefficient to multiply all dissimilarities, > 0", "1");
  	}



	void body () const final
	{
		const string dissimFName = getArg ("dissims");
		const string scaleFName  = getArg ("scales");
		const Real coeff         = str2real (getArg ("coeff"));
		ASSERT (coeff > 0.0);
		
		
		Vector<Scale> scales;  scales. reserve (16); // PAR
		{
			LineInput f (scaleFName);
			while (f. nextLine ())
			{
	      scales << Scale (f. line, coeff);
	      const size_t i = scales. size ();
	      if (i >= 2)
	      {
	      	if (  scales [i - 2]. unit 
	      		  > scales [i - 1]. unit
	      		 )
	      	  throw runtime_error ("Scale units must increase");
	      	if (  scales [i - 2]. dissim_max () 
	      		  > scales [i - 1]. dissim_max ()
	      		 )
	      	  throw runtime_error ("dissim_max must increase");
	      }
	    }
		}
		QC_ASSERT (scales. size () >= 2);

      
    Vector<ObjPair> objPairs;  objPairs. reserve (1000);  // PAR
    {
    	LineInput f (dissimFName);
    	while (f. nextLine ())
    		objPairs << move (ObjPair (f. line));
    }

    QC_ASSERT (objPairs. size () % scales. size () == 0);
    const size_t n = objPairs. size () / scales. size ();  // # Pairs
    
    const ONumber on (cout, 6, true);
    Progress prog (n, 1000);  // PAR
	  FOR (size_t, i, n)
	  {
	  	prog ();
	  #ifdef WEIGHT
	    Real dissim_sum = 0.0;
	    Real weight_sum = 0.0;
	  #else
	    Real dissim = NaN;
	  #endif
  	  FFOR (size_t, j, scales. size ())
      {
      	const size_t k = i + j * n;
        if (objPairs [k]. name1 != objPairs [i]. name1)
        	throw runtime_error ("Object " + toString (i + 1) + ", scale " + toString (j + 1) 
        	                     + " has a different name1: " + objPairs [k]. name1 + " vs. " + objPairs [i]. name1
        	                    );
        ASSERT (objPairs [k]. name2 == objPairs [i]. name2);
        const Real raw = objPairs [k]. dissim;  // isNan (objPairs [k]. dissim) ? NaN : max (dissim_max_prev, objPairs [k]. dissim);  // Bad
  	  #ifdef WEIGHT
        const Real dissim = scales [j]. raw2dissim (raw);
        if (isNan (dissim))
        	continue;
        ASSERT (dissim >= 0.0);
        dissim_sum += scales [j]. weight * dissim;
        weight_sum += scales [j]. weight;
      #else
        dissim = scales [j]. raw2dissim (raw);
        if (! isNan (dissim))
        {
          ASSERT (dissim >= 0.0);
          if (! dissim && j)  // Contradiction
            dissim = NaN;
          break;
        }
      #endif
      }
      cout         << objPairs [i]. name1 
           << '\t' << objPairs [i]. name2
           << '\t' << 
		       #ifdef WEIGHT
             dissim_sum / weight_sum
           #else
             dissim 
           #endif
           << endl;
    }
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);

}



