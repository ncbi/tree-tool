// evolution.hpp

#ifndef EVOLUTION_HPP
#define EVOLUTION_HPP

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
//#include "../dm/matrix.hpp"
#include "../dm/dataset.hpp"
using namespace DM_sp;



namespace DistTree_sp
{


Real intersection2dissim (Real size1,
                          Real size2,
                          Real intersection,
                          Real intersection_min,
	                        Prob sizes_ratio_min,
	                        bool ave_arithP);
  // Return: >= 0; may be NAN
  // Symmetric



struct Hashes : Vector<size_t>
{
	explicit Hashes (const string &fName);
	Hashes ()
	  {}
	
	Real getDissim (const Hashes &other,
	                size_t intersection_min,
	                Prob hashes_ratio_min) const
		{ return intersection2dissim ((Real) size (), (Real) other. size (), (Real) getIntersectSize (other), (Real) intersection_min, hashes_ratio_min, false); }
};



// --> dataset.hpp ??
struct DissimAverage : Root
{
	Real power {NAN};
	
	
	struct DissimAttr : Named
	{
		const DissimAverage& da;
		const PositiveAttr1* attr {nullptr};

		Real center {NAN};
			// > 0
		// For centered value's
		Real var {NAN};
			// >= 0
		Real value {NAN};
			// >= 0
		mutable bool outlier {false};
		WeightedMeanVar outlierMV;  
		  // outlier's statistic
		

		DissimAttr (const DissimAverage &da_arg,
		            const PositiveAttr1* attr_arg,
		            Real center_arg);
		DissimAttr (const DissimAverage &da_arg,
		            const string &line,
		            bool loadStat);
		void qc () const override;
		void saveText (ostream &os) const override
		  { os         << name
 		       << '\t' << center
 		       << '\t' << var
 		       << '\t' << outlierMV. getMean ()
 		       << endl;
 		  }

		  
		static bool goodValue (Real value) 
		  { return ! isNan (value) && value != INF; }
		Real getSD () const
		  { return sqrt (var); }
		bool bad () const
		  { return getSD () >= 0.8; }  // PAR
		Real getWeight () const
		  { return 1 / var; }
		  // Assumption: DissimAttr's are independent, which allows re-weighting of DissimAttr's if some value's are missing
		void setValue (Real value_arg);
	    // Output: value
	private:
		friend struct DissimAverage;
		void setOutlier (Real value_target) const;
	    // Output: outlier
    void setValue (size_t objNum)
			{ setValue ((* checkPtr<const PositiveAttr1> (attr)) [objNum]); }
    void setVar (const PositiveAttr1& averageAttr);
      // Output: var
	};
	Vector<DissimAttr> dissimAttrs;
	
	
	explicit DissimAverage (Real power_arg)
	  : power (power_arg)
	  {}
	DissimAverage (Real power_arg,
	               const string &fName,
	               bool loadStat);
  void qc () const override;
	void saveText (ostream &os) const override
    { for (const DissimAttr& dissimAttr : dissimAttrs)
    	  if (! dissimAttr. bad ())
		      dissimAttr. saveText (os);
		}
   
  Real get () const;
    // Input: DissimAttr::value
    // Output: DissimAttr::outlier
  void calibrate (PositiveAttr1& averageAttr);
    // Output: (averageAttr)[], DissimAttr::{var,outlierMV}
  Real getVar () const
    { Real var = 0;
    	for (const DissimAttr& dissimAttr : dissimAttrs)
    		if (! dissimAttr. bad ())
    		  var += dissimAttr. getWeight ();
    	return 1 / var;
    }    	
private:
	MVector getVars () const;
  void setValues (size_t objNum)
    { for (DissimAttr& dissimAttr : dissimAttrs)
		    dissimAttr. setValue (objNum);
		}
  Real setVars (const PositiveAttr1& averageAttr)
    { MeanVar mv;
    	for (DissimAttr& dissimAttr : dissimAttrs)
    	{ dissimAttr. setVar (averageAttr);
		    if (! dissimAttr. bad ())
		    	mv << dissimAttr. getSD ();
		  }
		  return mv. getMean ();
		}
public:
};
	


}



#endif


