// clust.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
using namespace DM_sp;



namespace
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Clustering as the decomposition of a mixture of multivariate Normal distributions")
  	{
  	  addPositional ("file", dmSuff + "-file");	  
  	  addPositional ("clusters_max", "Max. numebr of clusters");	  
  	  addPositional ("sd_min", "Min. SD of each variable in each cluster");	  
  	  addPositional ("prob_min", "probability threshold for the nominal attribute indicating the cluster; 0 - no nominal attribute");
  	}
	
	
	
	void body () const final
	{
		const string inFName      =               getArg ("file");
		const size_t clusters_max = str2<size_t> (getArg ("clusters_max"));
		const Real sd_min         = str2real     (getArg ("sd_min"));
		const Prob prob_min       = str2<Prob>   (getArg ("prob_min"));
		ASSERT (positive (sd_min));


    Dataset ds (inFName);
    const Sample sm (ds);
    const Space1<NumAttr1> sp (ds, true);

	  const Clustering cl (sm, sp, clusters_max, sd_min, false, 0.05);  // PAR
    cl. print (cout);
    
    VectorPtr<Attr1> attrs;
    {
      Common_sp::AutoPtr <const Space1<ProbAttr1>> spOut (cl. createSpace (ds));
      attrs << *spOut;
    }
    if (prob_min)
      attrs << cl. createNominAttr ("Cluster", prob_min, ds);

    cout << endl;
    sm. save (attrs, cout);
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


