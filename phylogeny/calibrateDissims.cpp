// calibrateDissims.cpp  // --> dm/convex.cpp ??

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
#include "../dm/dataset.hpp"
using namespace DM_sp;
#include "evolution.hpp"
using namespace DistTree_sp;



namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Analyze proportional homoscedastic dissimilarities: centralize to 1, find variances. Print dissimilarity attributes statistics")
    {
      // Input
  	  addPositional ("file", dmSuff + "-file without \"" + dmSuff + "\"");
  	  addPositional ("power", "Power for raw dissimilarities to be raised to");
  	  // Output
  	  addKey ("output_dissim", "Output merged dissimilarity");
  	}



	void body () const final
	{
		const string inFName            = getArg ("file");
		const Real power                = str2real (getArg ("power"));
		const string output_dissimFName = getArg ("output_dissim");


    Dataset ds (inFName);
    ds. qc ();
    

	  DissimAverage dissimAve (power);  dissimAve. dissimAttrs. reserve (ds. attrs. size ());
	  for (const Attr* attr_ : ds. attrs)
	  {
	  	const PositiveAttr1* attr = attr_->asPositiveAttr1 (); 
	  	if (! attr)
	  		throw runtime_error (attr_->name + " should be Positive");
	    const_cast <PositiveAttr1*> (attr) -> inf2missing ();
	  	
	    const Sample sm (ds);
	  #if 0  // PAR
	  	Real center = NAN;
	  	Real scatter = NAN;
	  	attr->getAverageScatter (sm, center, scatter);
	  	// center = INF
	  #else
	    const Real center = attr->getMedian (sm);
	  #endif
	  	dissimAve. dissimAttrs << DissimAverage::DissimAttr (dissimAve, attr, center);
	  }
	  ds. qc ();
    dissimAve. qc ();

	  	  
    auto averageAttr = new PositiveAttr1 ("average", ds, 6);  // PAR
    dissimAve. calibrate (*averageAttr);
	 	dissimAve. print (cout);
	 	
	 	if (! output_dissimFName. empty ())
	 	{
	 		OFStream f (output_dissimFName + dmSuff);
	 		const Sample sm (ds);
	 	  sm. save (VectorPtr<Attr> {averageAttr}, f);
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



