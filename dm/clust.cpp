// clust.cpp

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
*   Clustering
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
using namespace DM_sp;
#include "version.inc"



namespace
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Clustering as the decomposition of a mixture of multivariate Normal distributions. Print clusters statsistics")
  	{
  	  version = VERSION;
  	  addPositional ("file", dmSuff + "-file");	  
  	  addPositional ("clusters_max", "Max. numebr of clusters");	  
  	  addPositional ("sd_min", "Min. SD of each variable in each cluster");	  
  	  addPositional ("prob_min", "probability threshold for the nominal attribute indicating the cluster; 0 - no nominal attribute");
  	  addKey ("out", "Output " + dmSuff + "-file with cluster probabilities");
  	  addKey ("threshold_SDs", "Number of SDs to define lower/upper boundaries for unidimensional clustering", "0");
  	}
	
	
	
	void body () const final
	{
		const string inFName      =               getArg ("file");
		const size_t clusters_max = str2<size_t> (getArg ("clusters_max"));
		const Real sd_min         = str2real     (getArg ("sd_min"));
		const Prob prob_min       = str2<Prob>   (getArg ("prob_min"));
		const string outFName     =               getArg ("out");
		const Real threshold_sds  = str2real     (getArg ("threshold_SDs"));
		QC_ASSERT (sd_min > 0.0);
		QC_ASSERT (threshold_sds >= 0.0);


    Dataset ds (inFName);
    const Sample sm (ds);
    const Space1<NumAttr1> sp (ds, true);
   
	  const Clustering cl (sm, sp, clusters_max, sd_min, false);  
    cl. print (cout);
    
    if (cl. mixt. getDim () == 1 && threshold_sds)
    {
      Real boundary_hi_prev = NaN;
      Prob p_prev = 0.0;
      Prob p_next = 1.0;
      FFOR (size_t, i, cl. mixt. components. size ())
      {
        const MultiNormal* mn = cl. getMultiNormal (i);
        ASSERT (mn);
        const Real mean = mn->mu [0];
        const Real sd = sqrt (mn->sigmaInflated. get (false, 0, 0));
        const Real boundary_lo = mean - threshold_sds * sd;
        if (boundary_hi_prev < boundary_lo)
          cout         << "Threshold:" 
               << '\t' << (boundary_hi_prev + boundary_lo) / 2.0 
               << '\t' << p_prev 
               << '\t' << p_next 
               << endl;
        boundary_hi_prev = mean + threshold_sds * sd;
        const Prob p = cl. mixt. components [i] -> prob;
        p_prev += p;
        p_next -= p;
      }
    }


    if (! outFName. empty ())
    {
      OFStream f (outFName);
      VectorPtr<Attr> attrs;
      attrs << cl. createSpace (ds);
      if (prob_min)
        attrs << cl. createNominAttr ("Cluster", prob_min, ds);
      sm. save (attrs, f);
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


