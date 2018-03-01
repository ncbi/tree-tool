// calibrateDissims.cpp

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
    : Application ("Analyze proportional dissimilarities: centralize to 1, find variance")
    {
      // Input
  	  addPositional ("file", dmSuff + "-file without \"" + dmSuff + "\"");
  	}



	void body () const final
	{
		const string inFName = getArg ("file");


    Dataset ds (inFName);
    ds. qc ();
    

	  DissimAverage dissimAve;  dissimAve. attrs. reserve (ds. attrs. size ());
	  for (const Attr* attr_ : ds. attrs)
	  {
	  	const PositiveAttr1* attr = attr_->asPositiveAttr1 (); 
	  	if (! attr)
	  		throw runtime_error (attr_->name + " should be Positive");
	  	
	    const Sample sm (ds);
	  #if 0  // PAR
	    const_cast <PositiveAttr1*> (attr) -> inf2missing ();
	  	Real center = NAN;
	  	Real scatter = NAN;
	  	attr->getAverageScatter (sm, center, scatter);
	  	// center = INF
	  #else
	    const Real center = attr->getMedian (sm);
	  #endif
	  	dissimAve. attrs << DissimAverage::DissimAttr (attr, center);
	  }
	  ds. qc ();
    dissimAve. qc ();

	  	  
    auto averageAttr = new PositiveAttr1 ("average", ds, 6);  // PAR
    dissimAve. calibrate (averageAttr);
	 	dissimAve. print (cout);
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



