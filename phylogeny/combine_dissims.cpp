// combine_dissim.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
using namespace DM_sp;



namespace 
{
	
	
struct ObjPair
{
	string name1;
	string name2;
	Real dissim;
	
	explicit ObjPair (const string &line)
	  { 
	  	istringstream iss (line);
	  	string dissimS;
	  	iss >> name1 >> name2 >> dissimS;
	  	ASSERT (iss. eof ());
	  	ASSERT (name1 < name2);
	  	dissim = str2real (dissimS);
	  	IMPLY (! isNan (dissim), dissim >= 0);
	  }
};



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Combine different dissimilarities for the same pairs of objects. Dissimilarities must have linear-exponential variance")
    {
  	  addPositional ("dissims", "File with lines: <obj1> <obj2> <dissimilarity>; # lines = # object pairs times # dissimilarities");
  	  addPositional ("scales", "File with equalizing scales for each dissimilarity, ordered by scale descending, first scale = 1");
  	}



	void body () const final
	{
		const string dissimFName = getArg ("dissims");
		const string scaleFName  = getArg ("scales");
		
		
		Vector<Real> scales;  scales. reserve (16); // PAR
		{
		  ifstream f (scaleFName);
      Real scale;		
      Real scale_prev = INF;
		  while (f >> scale)
		  {
		  	ASSERT (scale > 0);
		  	ASSERT (scale <= 1);
		    scales << scale;
		    ASSERT (scale_prev > scale);
		    scale_prev = scale;
		  }
		}
		ASSERT (scales. size () >= 2);
		ASSERT (scales [0] == 1);

      
    Vector<ObjPair> objPairs;  objPairs. reserve (1000);  // PAR
    {
    	LineInput f (dissimFName);
    	while (f. nextLine ())
    		objPairs << ObjPair (f. line);
    }
    ASSERT (objPairs. size () % scales. size () == 0);
    const size_t n = objPairs. size () / scales. size ();  // # Pairs
    Vector<Real> dissims (scales. size ());
    Vector<Real> vars    (scales. size ());
    Vector<Real> weights (scales. size ());
    const ONumber on (cout, 6, true);
	  FOR (size_t, i, n)
	  {
	  	dissims. clear ();
	  	vars.    clear ();
	  	weights. clear ();
	  	Real dissim_first = NAN;
  	  FFOR (size_t, j, scales. size ())
      {
      	const size_t k = i + j * n;
        ASSERT (objPairs [k]. name1 == objPairs [i]. name1);
        ASSERT (objPairs [k]. name2 == objPairs [i]. name2);

        const Real coeff = 1 / scales [j];
        ASSERT (coeff >= 1);

        const Real dissim = objPairs [k]. dissim * coeff;
        dissims << dissim;
        
        if (isNan (dissim_first) && ! isNan (dissim))
        	dissim_first = dissim;

        const Real dissim_raw_corrected = dissim_first / coeff;
        Real var = /* N * */ exp (dissim_raw_corrected) - 1;  // Linear-exponential variance   // PAR
        if (isNan (var))
        	var = INF;
        ASSERT (var >= 0);
        vars << var;

        const Real absWeight = isNan (dissim) || dissim == INF
                                  ? 0 
                                  : (1 / (sqr (coeff) * var));
        ASSERT (absWeight >= 0);
        weights << absWeight;
      }

	  	Real dissim = 0;
      Real var = 0;
      Real weight_sum = 0;
      for (const Real w : weights)
    		weight_sum += w;
      ASSERT (weight_sum >= 0);
    	if (weight_sum == 0)
	    {
	    	dissim = NAN;
	    	for (const Real d : dissims)
	    		if (! isNan (d))
	    		{
	    			dissim = d;
	    			break;
	    		}
	    }
	    else if (weight_sum < INF)
    	{
	      // Relative weights
	      for (Real &w : weights)
	      {
	      	w /= weight_sum;
	      	ASSERT (w >= 0);
	      }
	      // dissim
	  	  FFOR (size_t, j, scales. size ())
	  	    if (weights [j])
	     		  dissim += dissims [j] * weights [j];
	     	ASSERT (dissim > 0);
	      // var
	  	  FFOR (size_t, j, scales. size ())
	  	    if (weights [j])
		     		var += sqr (weights [j] / scales [j]) * vars [j];
	      ASSERT (var > 0);
	    }

      if ((isNan (dissim) || dissim == INF || dissim == 0) != (var == 0))
      {
      	cout << objPairs [i]. name1 << ' ' << objPairs [i]. name2 << endl;
      	cout << dissim << ' ' << var << endl;
	  	  FFOR (size_t, j, scales. size ())
	  	    cout << dissims [j] << '\t' << vars [j] << '\t' << weights [j] << endl;
      	ERROR;
      }
      
    #if 0
      cout << dissims;  
      cout << endl;
      cout << weights; 
      cout << endl;
      cout << dissim << endl;
    #endif      
      cout         << objPairs [i]. name1 
           << '\t' << objPairs [i]. name2
           << '\t' << dissim
           << '\t' << var
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



