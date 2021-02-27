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
#include "../version.inc"



namespace 
{
	
	
struct Scale
{
	Real unit {NaN};
	Real raw_max {NaN};

	
	Scale (const string &line,
	       Real coeff)
	  {
	  	istringstream iss (line);
	  	string raw_maxS;
	  	iss >> unit >> raw_maxS;
	  	unit *= coeff;
	  	raw_max = str2real (raw_maxS);
	  	QC_ASSERT (unit > 0.0);
	  	QC_ASSERT (raw_max > 0.0);
	  	ASSERT (! isNan (dissim_max ()));
	  }
	Scale () = default;

	  
	Real raw2dissim (Real raw) const
	  { if (isNan (raw))
	  	  return NaN;
	  	QC_ASSERT (raw >= 0.0);
	  	if (raw > raw_max)
	  	  return NaN;
	  	return raw * unit; 
	  }
	Real dissim_max () const
	  { return raw2dissim (raw_max); }
};
	
	
	
	
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Combine different dissimilarities for the same pairs of objects")
    {
      version = VERSION;
  	  addPositional ("dissims", "File with lines: <obj1> <obj2> <dissimilarity 1> <dissimilarity 2> ...");
  	  addPositional ("scales", "File with equalizing scales and max. values for each dissimilarity, ordered by <unit> * <raw_max>.\n\
Line format: <unit> <raw_max>");
      addKey ("coeff", "Coefficient to multiply all dissimilarities, > 0", "1");
  	}



	void body () const final
	{
		const string dissimFName = getArg ("dissims");
		const string scaleFName  = getArg ("scales");
		const Real coeff         = str2real (getArg ("coeff"));
		QC_ASSERT (coeff > 0.0);
		
		
		Vector<Scale> scales;  scales. reserve (16); // PAR
		{
			LineInput f (scaleFName);			
			while (f. nextLine ())
			{
	      scales << Scale (f. line, coeff);
	      const size_t i = scales. size ();
	      if (i >= 2)
	      	if (   scales [i - 2]. dissim_max () 
	      		  >= scales [i - 1]. dissim_max ()
	      		 )
	      	  throw runtime_error ("dissim_max must increase");
	    }
		}
		QC_ASSERT (scales. size () >= 2);

      
  	LineInput f (dissimFName, 1024 * 100, 1);
  	string name1;
  	string name2;
  	string dissimS;
   	Istringstream iss;
  	while (f. nextLine ())
    {
	  	iss. reset (f. line);
	  	iss >> name1 >> name2;
	  	if (name1 > name2)
	  	  swap (name1, name2);
	  	ASSERT (name1 < name2);

	    Real dissim = NaN;
	    Real dissim_prev = NaN;	    
	  	FFOR (size_t, i, scales. size ())
	  	{
	  	  dissimS. clear ();
	  	  iss >> dissimS;
	  	  QC_ASSERT (! dissimS. empty ());
  	  	const Real raw = str2real (dissimS);
  	  	QC_IMPLY (! isNan (raw), raw >= 0.0);
        dissim = scales [i]. raw2dissim (raw);
        if (dissim == 0.0 && i) 
          dissim = NaN;
        if (isNan (dissim))
        {
          if (! isNan (dissim_prev))
          {
            dissim = dissim_prev;
            break;
          }
        }
        else
        {
          ASSERT (dissim >= 0.0);
          if (isNan (dissim_prev))
            dissim_prev = dissim;
          else
          {
            // Merging dissimilarities of adjacent scales
            ASSERT (i);
            const Real prev_max = scales [i - 1]. dissim_max ();
            ASSERT (dissim_prev <= prev_max);
            const Real prev_prev_max = (i - 1 ? scales [i - 2]. dissim_max () : 0.0);
            ASSERT (prev_prev_max < prev_max);
            if (dissim_prev <= prev_prev_max)
              dissim = dissim_prev;
            else
            {
              const Prob weight = (dissim_prev - prev_prev_max) / (prev_max - prev_prev_max);
              dissim = weight * dissim + (1.0 - weight) * dissim_prev;
            }
            break;
          }
        }
      }
      cout         << name1 
           << '\t' << name2
           << '\t' << dissim 
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



