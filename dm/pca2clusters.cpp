// pca2clusters.cpp

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
*   Clusters made by PCA
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "numeric.hpp"
using namespace DM_sp;
#include "../version.inc"



namespace 
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Extract lists of core objects for each cluster made by mds or pca")
  	{
  	  version = VERSION;
  	  addPositional ("in", "Clustering output in Json");
  	  addPositional ("out", "File prefix for lists, each list is sorted");
  	  addPositional ("threshold", "Min. probability to belong to a cluster");
  	}



	void body () const final
	{
		const string fName   = getArg ("in");
		const string outName = getArg ("out");
		const Prob p_min     = str2real (getArg ("threshold"));
		ASSERT (isProb (p_min));


    JsonMap jClust (fName);
  //jClust. print (cout);
    
    Vector<StringVector> clusters (jClust. at ("clusters") -> getSize ());
    const size_t nObjs = jClust. at ("objs") -> getSize ();
    size_t skipped = 0;
    FOR (size_t, i, nObjs)
    {
      const Json* obj = jClust. at ("objs") -> at (i);
      const Json* attr = obj->at ("attr");
      const Real p = attr->at ("Cluster_prob") -> getDouble ();
      if (verbose ())
        cout << p << endl;
      if (p <= p_min)
      {
        skipped++;
        continue;
      }
      string clustS (attr->at ("Cluster") -> getString ());
      ASSERT (clustS [0] == 'C');
      clustS. erase (0, 1);
      const int clust = str2<int> (clustS) - 1;
      ASSERT (clust >= 0);
      clusters [(size_t) clust] << obj->at ("objName") -> getString ();
    }
    
    cout << "# Clusters: " << clusters. size () << endl;
    cout << "# Objects: " << nObjs << endl;
    if (skipped)
      cout << "# Unclustered: " << skipped << endl;
    
    FOR (size_t, i, clusters. size ())
    {
      OFStream ofs (string (), outName + "-" + toString (i + 1), string ());
      sort (clusters [i]);
      save (ofs, clusters [i], '\n');
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

