// combine_dissims.cpp

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

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
using namespace DM_sp;
#include "../version.inc"

#include "../common.inc"



namespace 
{
	
	
struct Scale
{
	Real unit {NaN};
	Real raw_max {NaN};
	//	
	Real dissim_min {NaN};
	Real dissim_center {NaN};

	
	Scale (const string &line,
         Real dissim_min_arg,
         Real dissim_center_arg)
    : dissim_min    (nvlReal (dissim_min_arg,    0.0))
    , dissim_center (nvlReal (dissim_center_arg, 0.0))
	  {
	  	istringstream iss (line);
	  	string unitS, raw_maxS;
	  	iss >> unitS >> raw_maxS;
	  	unit    = str2real (unitS);
	  	raw_max = str2real (raw_maxS);
	  	QC_ASSERT (unit > 0.0);
	  	QC_ASSERT (raw_max > 0.0);
	  	QC_ASSERT (dissim_min <= dissim_center);
	  	QC_ASSERT (dissim_center < dissim_max ());
	  	QC_ASSERT (valid ());
	  }
	Scale () = default;
	void saveText (ostream &os) const
	  { os         << unit 
	       << '\t' << raw_max 
	       << '\t' << dissim_min
	       << '\t' << dissim_center
	       << '\t' << dissim_max ()
	       << endl; 
	  }

	  
	bool valid () const 
	  { return ! isNan (unit); }
	Real raw2dissim (Real raw) const
	  { if (isNan (raw))
	  	  return NaN;
	  	QC_ASSERT (raw >= 0.0);
	  	return raw * unit; 
	  }
	Real dissim_max () const
	  { return raw2dissim (raw_max); }
	Prob getWeight (Real dissim,
	                bool first) const
	  { if (dissim < dissim_center)
	    {
	      if (first)
	        return 1.0;
        if (dissim < dissim_min)  // May be: dissim_min = dissim_center = 0
	        return 0.0;
	    	return (dissim - dissim_min) / (dissim_center - dissim_min);
	    }
	    if (dissim >= dissim_max ())
	      return 0.0;
	    if (dissim_max () == inf)
	      return 1.0;
	    return (dissim_max () - dissim) / (dissim_max () - dissim_center);
	  }
};




struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Combine different dissimilarities for the same pairs of objects")
    {
      version = VERSION;
  	  addPositional ("dissims", "File with lines: <obj1> <obj2> <dissimilarity 1> <dissimilarity 2> ...");
  	  addPositional ("scales", "File with equalizing scales and max. values for each dissimilarity, ordered by <unit> * <raw_max>.\n\
Line format: <unit> <raw_max>\n\
unit = 1 / (species_barrier for dissimilarity)\n\
<raw_max> can be inf");
  	  addKey ("coeff", "Coefficient to multiply the dissimilarity by", "1.0");
  	//addKey ("power", "Power to raise the dissimilarity in",  "1.0");
  	  addKey ("barrier", "Scale sequential number, 1-based, 0 - no barrier. If i > barrier then dissim[i] >= dissim[barrier]");
      addFlag ("print_raw", "Print raw dissimilarities");
  	}



	void body () const final
	{
		const string dissimFName = getArg ("dissims");
		const string scaleFName  = getArg ("scales");
	  const Real   coeff       = str2real (getArg ("coeff"));
	//const Real   power       = str2real (getArg ("power"));
	  const size_t barrier    = (size_t) arg2uint ("barrier") - 1;
		const bool   print_raw   = getFlag ("print_raw");		
	  QC_ASSERT (coeff > 0.0);
	//QC_ASSERT (power > 0.0);
	  
	  
		Vector<Scale> scales;  scales. reserve (16); // PAR
		{
			LineInput f (scaleFName);
			Scale prev_prev;
			Scale prev;
			while (f. nextLine ())
			{
			  const Scale scale (f. line, prev_prev. dissim_max (), prev. dissim_max ());
	      scales << scale;
	      if (   prev. valid ()
	      	  && prev. dissim_max () >= scale. dissim_max ()
	      	 )
	      	throw runtime_error ("dissim_max must increase");
        prev_prev = prev;
        prev = scale;
	    }
	    // At least one valid Scale should exist
  		QC_ASSERT (prev. valid ());
  		QC_ASSERT (prev. raw_max == inf);
		}
		QC_ASSERT (scales. size () >= 2);
		QC_IMPLY (barrier != no_index, barrier < scales. size ());

      
  	LineInput f (dissimFName, 1);  // PAR
  	while (f. nextLine ())
    {
     	stringstream iss (f. line);

      {
      	string name1;
      	string name2;
  	  	iss >> name1 >> name2;
  	  	QC_ASSERT (! name2. empty ());
  	  	if (name1 > name2)
  	  	  swap (name1, name2);
  	  	ASSERT (name1 < name2);
        cout         << name1 
             << '\t' << name2;
      }

	    Real dissim_weighted_sum = 0.0;
	    Prob weight_sum = 0.0;
	    Prob weight_whole = 1.0;
	    size_t n = 0;
	  	size_t i = 0;
	  	for (; i < scales. size (); i++)
	  	{
	  	  const Scale& scale = scales [i];
      	string dissimS;
	  	  iss >> dissimS;
	  	  QC_ASSERT (! dissimS. empty ());

  	  	const Real raw = str2real (dissimS);
  	  	if (print_raw)
  	  	  cout << '\t' << raw;
        if (isNan (raw))
        {
          if (n)
            break;
        }
        else
        {
        //if (n >= 2)
          //break;
          const Real dissim_ = scale. raw2dissim (raw);
          const Real dissim = barrier == no_index || i <= barrier ? dissim_ : max (dissim_, scales [barrier]. dissim_max ());
          const Prob weight_raw = scale. getWeight (dissim, ! weight_sum);
          ASSERT (dissim >= 0.0);
          ASSERT (isProb (weight_raw));
          if (weight_raw)
          {
            ASSERT (isProb (weight_whole));
            const Prob weight = weight_whole * weight_raw;
            weight_whole *= 1.0 - weight_raw;
            dissim_weighted_sum += weight * dissim;
            weight_sum          += weight;
            n++;
          }
          if (verbose ())
          {
            cout << endl;
            cout << "scale: ";
            scale. saveText (cout);
            PRINT (raw);
            PRINT (dissim);
            PRINT (weight_raw);
            PRINT (i);
            PRINT (dissim_weighted_sum);
            PRINT (weight_sum);
          }
        }

        if (   ! i
            && dissim_weighted_sum == 0.0
            && weight_sum > 0.0
           )
          break;
      }
      ASSERT (isProb (weight_sum));
      if (   i
          && dissim_weighted_sum == 0.0
          && weight_sum > 0.0
         )
        weight_sum = 0.0;
      if (verbose ())
      {
        cout << endl;
        PRINT (i);
        PRINT (dissim_weighted_sum);
        PRINT (weight_sum);
      }


      const Real d = dissim_weighted_sum / weight_sum;
    //cout << '\t' << coeff * pow (d, power) << endl;
      cout << '\t' << coeff * d << endl;
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



