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
	    if (isNan (unit))
	    {
	      QC_ASSERT (isNan (raw_max));
	    }
	    else
	    {
  	  	QC_ASSERT (unit > 0.0);
  	  	QC_ASSERT (raw_max > 0.0);
  	  	QC_ASSERT (dissim_min <= dissim_center);
  	  	QC_ASSERT (dissim_center < dissim_max ());
  	  }
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
	Real getWeight (Real dissim) const
	  // Return: 0..1
	  { if (dissim < dissim_min)  // May be: dissim_min = dissim_center = 0
	      return 0.0;
	    if (dissim >= dissim_max ())
	      return 0.0;
	    if (dissim < dissim_center)
	      return (dissim - dissim_min) / (dissim_center - dissim_min);
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
Line format: <unit> <raw_max>");
      addFlag ("print_raw", "Print raw dissimilarities");
  	}



	void body () const final
	{
		const string dissimFName = getArg ("dissims");
		const string scaleFName  = getArg ("scales");
		const bool   print_raw   = getFlag ("print_raw");
		
		
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
	          && scale. valid ()
	      	  && prev. dissim_max () >= scale. dissim_max ()
	      	 )
	      	throw runtime_error ("dissim_max must increase");
	      if (scale. valid ())
	      {
	        prev_prev = prev;
	        prev = scale;
	      }
	    }
	    // At least one valid Scale should exist
  		QC_ASSERT (prev. valid ());
  		QC_ASSERT (prev. raw_max == inf);
		}
		QC_ASSERT (scales. size () >= 2);

      
  	LineInput f (dissimFName, 1024 * 100, 1);  // PAE
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
	    Real weight_sum = 0.0;
	    bool first = true;
	    size_t n = 0;
	  	for (const Scale& scale : scales)
	  	{
      	string dissimS;
	  	  iss >> dissimS;
	  	  QC_ASSERT (! dissimS. empty ());

  	  	if (! scale. valid ())
  	  	  continue;

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
          if (n >= 2)
            break;
          const Real dissim = scale. raw2dissim (raw);
          const Real weight = scale. getWeight (dissim);
          ASSERT (dissim >= 0.0);
          ASSERT (weight >= 0.0);
          if (weight)
          {
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
            PRINT (weight);
            PRINT (first);
            PRINT (dissim_weighted_sum);
            PRINT (weight_sum);
          }
        }

        if (   first
            && dissim_weighted_sum == 0.0
            && weight_sum > 0.0
           )
          break;
          
        first = false;
      }
      if (   ! first
          && dissim_weighted_sum == 0.0
          && weight_sum > 0.0
         )
        weight_sum = 0.0;
      if (verbose ())
      {
        PRINT (first);
        PRINT (dissim_weighted_sum);
        PRINT (weight_sum);
      }

      cout << '\t' << dissim_weighted_sum / weight_sum << endl;
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



