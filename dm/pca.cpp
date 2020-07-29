// pca.cpp

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
*   Principal components analysis
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "numeric.hpp"
#include "dataset.hpp"
#include "prediction.hpp"
using namespace DM_sp;
#include "../version.inc"



namespace
{

const string clusterAttrName = "Cluster";



struct ThisApplication : Application
{	
  ThisApplication ()
  : Application ("Principal components")
	{
	  version = VERSION;
	  
	  addPositional ("file", dmSuff + "-file without the extension");	  
	  addPositional ("suf", "Files <file>-<suf>.* will be created");	  
	  // Stopping criteria
	  addKey ("maxAttr", "Max. # attributes (0 - to be determined automatically)", "100"); 
	  addKey ("maxTotalExpl", "Min. total explained fraction (0..1)", "1");  // was: 0.99
	  addKey ("minExpl", "Min. explained fraction (0..1)", "0.005");  // was: 0.001
	  // Class attribute
	  addKey ("class", "Class attribute name");  // test ??
	  // Clustering
	  addKey ("maxClusters", "Max. number of clusters; 1 - no clustering, otherwise make a clustering in the MDS space, create an attribute " + strQuote (clusterAttrName), "1");
	  addKey ("minClusteringSDRel", "Min. SD for a cluster / data SD", "0.01"); 
	  addFlag ("mergeClose", "Merge close clusters");
	  addKey ("clusteringDir", "If specified and the number of clusters > 1 then save the original data restricted to each cluster in this directory");
	  //
	  addKey ("attrPValue", "Max. p-value for a cluster-characteristic attribute", "1e-5");  // was : 0.01
	  //
	  addFlag ("mds", "Print the MDS of attributes");
	  addKey ("outlierEValue", "Max. E-value for an outlier", "0.1");  // was: 0.01
	}
	
	
	
	void body () const final
	{
		const string inFName          =             getArg("file");
		const string suf              =             getArg("suf");
		//
		const uint maxAttr            = str2<uint> (getArg("maxAttr"));
		const Prob maxTotalExpl       = str2<Prob> (getArg("maxTotalExpl"));
		const Prob minExpl            = str2<Prob> (getArg("minExpl"));
		//
		const bool mds                =             getFlag("mds");
		const Real outlierEValue      = str2<Real> (getArg("outlierEValue"));
		//
		const string classAttrName    =             getArg("class");
		//
		const size_t maxClusters      = str2<size_t> (getArg ("maxClusters"));
		const string clusteringDir    =             getArg("clusteringDir");
		const Real minClusteringSDRel = str2<Real> (getArg("minClusteringSDRel"));
	  const bool mergeClose         =             getFlag("mergeClose");  
	  //
		const Prob attrPValue         = str2<Prob> (getArg("attrPValue"));
		
		ASSERT (maxAttr);
		ASSERT (isProb (maxTotalExpl));
		ASSERT (isProb (minExpl));
		ASSERT (outlierEValue >= 0);
		ASSERT (maxClusters);
		IMPLY (! clusteringDir. empty (), maxClusters > 1);
		ASSERT (minClusteringSDRel >= 0);
		ASSERT (isProb (attrPValue));
		
    
    const string outFName = inFName + "-" + suf;
    OFStream osPar ("", outFName, "txt");
    OFStream osDat ("", outFName, dmExt);

    Dataset ds (inFName);
    
    if (ds. attrs. size () <= 0)
      throw runtime_error ("No attributes");
    
    const NominAttr1* classAttr = nullptr;
    if (! classAttrName. empty ())
    {
      const Attr* attr = ds. name2attr (classAttrName);
      ASSERT (attr);
      classAttr =  attr->asNominAttr1 ();
      ASSERT (classAttr);
    }

    const Sample sm_orig (ds);
    Sample sm (ds);
    if (sm. mult_sum <= 0.0)
      throw runtime_error ("Too small data size");

    Space1<Attr1> spRaw (ds, true);
    spRaw. removeAttr (*classAttr);

    const VectorPtr<Attr> attrsOrig (spRaw);

    const Space1<RealAttr1> spStnd (spRaw. standardize (sm, ds));
        
    const MultiVariate<NumAttr1> an (sm, spStnd);

    MultiNormal mn; 
    {
      mn. analysis = & an;
      mn. qc ();
      mn. estimate ();
      if (verbose ())
      	mn. print (osPar);    	
    }
    	
    cerr << "PC ..." << endl;
    const PrinComp pc (sm, an. space, mn, maxAttr, maxTotalExpl, minExpl, 1e-4/*PAR*/);
    pc. qc ();
    pc. print (osPar);
    
    const Space1<RealAttr1> spOut (pc. createSpace ("PC_", ds));
    
    const size_t eValueDecimals = 5; // PAR

    cerr << "Outliers ..." << endl;
    RealAttr1* outlierScore = new RealAttr1 ("Outlier_E-Value", ds, eValueDecimals); 
    {
      ExtBoolAttr1* outlier = new ExtBoolAttr1 ("Outlier", ds);  
      for (Iterator it (sm); it ();)  
      {
        const Real score = pc. getChi2 (*it);
        if (score > 0.0)
          (*outlierScore) [*it] = log (score);
      }
      UniVariate<NumAttr1> outlierAn (sm, *outlierScore);     
      // score ~ Exponential ??
      Normal normal;
      normal. analysis = & outlierAn;
	    normal. estimate ();
	    normal. qc ();
	    if (verbose ())
	      normal. print (cout);
      for (Iterator it (sm); it ();)  
      {
        const Prob pValue = 1 - normal. cdf ((*outlierScore) [*it]);
        // Too simplistic criterion ??
        const Real eValue = pValue * sm. mult_sum;
        if (eValue <= outlierEValue * 10)  // PAR
          const_cast <Obj*> (ds. objs [*it]) -> comment += string ("  [E-value = ") + real2str (eValue, eValueDecimals) + "]";
        (*outlierScore) [*it] = eValue;
        (*outlier) [*it] = (eValue <= outlierEValue ? etrue : efalse);  
        if ((*outlier) [*it])
        {
          if (clusteringDir. empty () || ! jRoot)
          {
            osPar << "Outlier: " << ds. objs [*it] -> name << " " << eValue << ": ";
            ONumber on (osPar, 1, false);
            bool first = true;
            for (const RealAttr1* attr : spStnd)
            {
              const Real value = (*attr) [*it];                
              if (fabs (value) >= 4.0)  // PAR
              {
                if (! first)
                  osPar << ';';
                osPar << attr->name << '=' << value;
                first = false;
              }
            }
            osPar << endl;
          }
          sm. mult [*it] = 0;
        }
      }
    }
    sm. finish ();


    Space1<Attr1> sp1 (spOut);


	  if (maxClusters > 1 && pc. getOutDim ())
	  {
      cerr << "Clustering ..." << endl;
      Space1<NumAttr1> spPC (spOut);
      spPC. removeConstants ();
      const Real globalSD = sqrt ((Real) spStnd. size ());
  	  const Clustering cl (sm, spPC, maxClusters, globalSD * minClusteringSDRel, /*true*/ false); 
      cl. qc ();
  	  if (verbose ()) 
        cl. print (cout);
      cout << "# Clusters = " << cl. mixt. components. size () << endl;
      cout << "Entropy = " << cl. mixt. getEntropy_est () << endl;
      JsonArray* jClusts = jRoot ? new JsonArray (jRoot, "clusters") : nullptr;
      cl. processSubclusters (clusterAttrName, mergeClose, attrsOrig, jClusts, attrPValue, clusteringDir, getFileName (inFName), ds, sp1);
    }


    sp1 << outlierScore;
    {
      const VectorPtr<Attr> sm_all (sp1);
      sm_orig. save (sm_all, osDat);
    }


    if (mds)
    {
      cerr << "MDS of attributes ..." << endl;
      Vector<Real> quality;
    	const Dataset dsMds (pc. createAttrMds ("PC_", quality));    	
    	ASSERT (quality. size () == dsMds. attrs. size ());
      OFStream osMds ("", outFName, "mds");
    	dsMds. print (osMds);
    	if (jRoot)
    	{
        auto jMds = new JsonMap (jRoot, "attr_mds");
        const Sample smMds (dsMds);
        const Space1<RealAttr1> spMds (dsMds, true);
        spMds. toJson (smMds, jMds, "objs");
        // Cf. Eigens::toJson()
        {
          auto j = new JsonMap (jMds, "eigens");
          new JsonBoolean (true, j, "psd");
          {
            auto jArr = new JsonArray (j, "explained");  
            FOR (size_t, i, quality. size ())
              new JsonDouble (quality. at (i), 6, jArr);  // PAR
          }
        }
      }
    }


    const EigensLinearRegression elr (pc. eigens);
    elr. qc ();
    
    
    if (classAttr)
      sp1 << classAttr;


    if (jRoot)
    {
      new JsonDouble (outlierEValue, eValueDecimals, jRoot, "outlierEValue");
      sp1. toJson (sm_orig, jRoot, "objs");      
        // Attributes: {PC_<N>,} <processSubclusters()>, outlierScore
      pc. eigens. toJson (jRoot, "eigens");
       // Correlations of original attributes with PC's ??
      elr. toJson (jRoot, "eigenValueLR");  
      new JsonString (classAttrName, jRoot, "class");  // May be ""
      ds. comments2Json (jRoot, "comments");
    }    
    
    
    ds. qc ();
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


