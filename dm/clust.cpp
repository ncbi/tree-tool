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

	  const Clustering cl (sm, sp, clusters_max, sd_min, false);  
    cl. print (cout);
    
    VectorPtr<Attr> attrs;
    attrs << cl. createSpace (ds);
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


