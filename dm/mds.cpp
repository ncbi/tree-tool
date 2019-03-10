// mds.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "numeric.hpp"
#include "matrix.hpp"
#include "dataset.hpp"
#include "prediction.hpp"
using namespace DM_sp;



namespace 
{


const string clusterAttrName = "Cluster";
constexpr uint noAttrType = 3;



struct ThisApplication : Application
{
  ThisApplication ()
  : Application ("Multidimensional scaling (linear)")
  { 
	  // Input
	  addPositional ("file", dmSuff + "-file without the extension");
	  addKey ("attr", "Attribute name of real-valued object-object comparison table in the " + dmSuff + "-file, if empty then the PCA of all attributes except <class> is done");
	  addKey ("attrType", "<attr> type: 0-similarity, 1-dissimilarity, 2-squared dissimilarity", toString (noAttrType));
    // Cf. pca.cpp
	  // Stopping criteria
	  addKey ("maxAttr", "Max. # attributes (0 - to be determined automatically)", "100"); 
	  addKey ("maxTotalExpl", "Min. total explained fraction (0..1)", "1.0");  // was: 0.99
	  addKey ("minExpl", "Min. explained fraction (0..1)", "0.005");  // was: 0.001
	  // Class attribute
	  addKey ("class", "Class attribute name");
	  // Clustering
	  addKey ("maxClusters", "Max. number of clusters; 1 - no clustering, otherwise make a clustering in the MDS space, create an attribute " + strQuote (clusterAttrName), "1");
	  addKey ("minClusteringSDRel", "Min. SD for a cluster / data SD", "0.01"); 
	  addFlag ("mergeClose", "Merge close clusters");
	  addKey ("clusteringDir", "If specified and the number of clusters > 1 then save the original data restricted to each cluster in this directory");
	  //
	  addKey ("attrPValue", "Max. p-value for a cluster-characteristic attribute", "1e-5");  // was : 0.01
  }



	void body () const final
	{
		const string fName            =               getArg ("file");              
		const string attrName         =               getArg ("attr");          
		const uint attrType           = str2<uint>   (getArg ("attrType"));             
    //
		const size_t maxAttr          = str2<size_t> (getArg ("maxAttr"));          
		const Prob maxTotalExpl       = str2<Prob>   (getArg ("maxTotalExpl"));     
		const Prob minExpl            = str2<Prob>   (getArg ("minExpl"));          
		//
		const string classAttrName    =               getArg ("class");             
		//
		const size_t maxClusters      = str2<size_t> (getArg ("maxClusters"));
		const string clusteringDir    =               getArg ("clusteringDir");    
		const Real minClusteringSDRel = str2<Real>   (getArg ("minClusteringSDRel"));
	  const bool mergeClose         =               getFlag ("mergeClose");  
	  //    
		const Prob attrPValue         = str2<Prob>   (getArg ("attrPValue"));      

    IMPLY (attrName. empty (), attrType == noAttrType);
		ASSERT (maxAttr);
		ASSERT (isProb (maxTotalExpl));
		ASSERT (isProb (minExpl));		
		ASSERT (maxClusters);
		IMPLY (! clusteringDir. empty (), maxClusters > 1);
		ASSERT (minClusteringSDRel >= 0);
		ASSERT (isProb (attrPValue));
		
		
    Dataset ds (fName);
    
    const Sample sm (ds);
    if (! positive (sm. mult_sum))
      exit (1);

    const NominAttr1* classAttr = nullptr;
    if (! classAttrName. empty ())
    {
      const Attr* attr = ds. name2attr (classAttrName);
      ASSERT (attr);
      classAttr =  attr->asNominAttr1 ();
      ASSERT (classAttr);
    }
  
    const VectorPtr<Attr> attrsOrig (ds. attrs);
    Space1<Attr1> sp1 (ds, false);

    
    RealAttr2* sim = nullptr;
    if (attrName. empty ())
    {
      cerr << "PCA via MDS ..." << endl;
      Space1<Attr1> spRaw (ds, true);
      if (classAttr)
      {
        sp1 << classAttr;
        spRaw. removeAttr (*classAttr);
      }
      const Space1<RealAttr1> spStnd (spRaw. standardize (sm, ds));
      sim = getSimilarity (spStnd, attrName + "_sim", ds);
    }
    else
    {
      {
        Space1<Attr1> sp1_ (ds, true);
        sp1 << sp1_;
      }
      
      // dist
      const RealAttr2* dist = nullptr;
      {
        const Attr* attr = ds. name2attr (attrName);
        ASSERT (attr);
        dist = attr->asRealAttr2 ();
      }
      ASSERT (dist);
      
      if (const PositiveAttr2* distPos = dist->asPositiveAttr2 ())
      {
        if (const size_t n = distPos->getInfCount ())
          ds. comments << "# Infinities converted to missings: " + toString (n);
        const_cast <PositiveAttr2*> (distPos) -> inf2missing ();
      }
      
  	  // Check dist->matr  
  	  {
        Real maxCorrection;
        size_t row_bad, col_bad;
        const_cast <RealAttr2*> (dist) -> matr. symmetrize (maxCorrection, row_bad, col_bad);
        if (maxCorrection > 2 * pow (10, - (Real) dist->decimals))
          ds. comments << "maxCorrection = " + toString (maxCorrection) + " at " + ds. objs [row_bad] -> name + ", " + ds. objs [col_bad] -> name;
      }
  	  if (   attrType == 1 
  	      || attrType == 2
  	     )
  	  {
  	    size_t row;
  	    if (! dist->matr. zeroDiagonal (row))
        {
         	cout << dist->name << "[" << row + 1 << "] [" << row + 1 << "] != 0" << endl;
         	exit (1);
        }
        if (dist->matr. min () < 0)  // Actually can ??
        {
         	cout << "Matrix cannot have negative values" << endl;
         	exit (1);
        }
      }
  
      sim = new RealAttr2 (attrName + "_sim", ds, *dist); 
      sim->decimals++;
      {
        Matrix& matr = sim->matr;
        if (const size_t missings = matr. undefined2mean ())
          ds. comments << "# Missings converted to row/column means: " + toString (missings);
        Matrix rowMean;
        Real totalMean;
    	  switch (attrType)
        {
          case 0: matr. centerSimilarity (rowMean, totalMean);
                  break;
          case 1: matr. sqrAll ();
                  matr. sqrDistance2centeredSimilarity (rowMean, totalMean);
                  break;
          case 2: matr. sqrDistance2centeredSimilarity (rowMean, totalMean);
          	      break;
          default: throw runtime_error ("No attrType");
        }
      }
    }
    ASSERT (sim);
	  if (verbose ())
	  {
	    cout << endl;
	    cout << "Double-cenetered generalized similarities:" << endl;
	    sim->print (cout);
	  }
	

    WeightedMeanVar globalVarSum;
    for (Iterator it (sm); it ();)  
      globalVarSum. add (sim->get (*it, *it), it. mult);
    const Real globalSD = sqrt (globalVarSum. getMean ());
    
    
    if (false)   
    {
      // Data noise
      Normal norm;
      norm. setParam (0, 50);  // PAR
      Space1<RealAttr1> sp (ds, false);
      FOR (size_t, i, 10000)  // PAR  
      {
        auto attr = new RealAttr1 ("X" + toString (i + 1), ds, 2);  // PAR
        sp << attr;
      //norm. setSeed ((ulong) i + 1);
        FOR (size_t, objNum, ds. objs. size ())
          (*attr) [objNum] = norm. rand ();
      }
      RealAttr2* simExtra = getSimilarity (sp, "simExtra", ds);
      sim->matr. add (false, simExtra->matr, false, 1);
      Matrix rowMean;
      Real totalMean;
      sim->matr. centerSimilarity (rowMean, totalMean);
    }

    if (false)
    {
      // Measurement noise
      Matrix& matr = sim->matr;
      matr. print (cout);
      Normal norm;
      norm. setParam (0, sqr (globalSD) * 0.5);  // PAR
      FOR (size_t, row, matr. rowsSize (false))
        FOR (size_t, col, row)
        {
          matr. putInc (false, row, col, norm. rand ());
          matr. put (false, col, row, matr. get (false, row, col));
        }
      Matrix rowMean;
      Real totalMean;
      matr. centerSimilarity (rowMean, totalMean);
      matr. print (cout);
    }
    

    cerr << "MDS ..." << endl;
    const Mds mds (sm, *sim, maxAttr, maxTotalExpl, minExpl); 
    mds. qc ();
    mds. print (cout);
    cout << endl;
    
    const Space1<RealAttr1> spMds_ (mds. createSpace ("PC_", "IM_", ds));
    ds. qc ();
    
    if (spMds_. empty ())
      return;
        
    // outliers, see pca.cpp ??
        
    sp1 << spMds_;
    

    if (maxClusters > 1 && mds. getOutDim () && positive (globalSD))
	  {
      cerr << "Clustering ..." << endl;
      Space1<NumAttr1> spMds (spMds_);
      spMds. removeConstants ();
  	  const Clustering cl (sm, spMds, maxClusters, globalSD * minClusteringSDRel, /*true*/ false);
      cl. qc ();
  	  if (verbose ())
        cl. print (cout);
      cout << "# Clusters = " << cl. mixt. components. size () << endl;
      cout << "Entropy = " << cl. mixt. getEntropy_est () << endl;
      JsonArray* jClusts = jRoot ? new JsonArray (jRoot, "clusters") : nullptr;
      cl. processSubclusters (clusterAttrName, mergeClose, attrsOrig, jClusts, attrPValue, clusteringDir, getFileName (fName), ds, sp1);
    }
    

    EigensLinearRegression elr (mds. eigens);
    elr. qc ();
      

    IMPLY (classAttr, sp1. contains (classAttr));
    sp1. qc ();
    
    ds. qc ();
    
    
    if (jRoot)
    {
      sp1. toJson (sm, jRoot, "objs");
        // Attributes: all original [with <classAttrName>] (except <attrName>), {PC_<N>,} <processSubclusters()>
      mds. eigens. toJson (jRoot, "eigens");
      elr. toJson (jRoot, "eigenValueLR");

      if (! attrName. empty ())
      {
        new JsonString (attrName, jRoot, "distAttr");  
        new JsonInt ((int) attrType, jRoot, "distType");  
      }

      new JsonString (classAttrName, jRoot, "class");  // May be ""
      ds. comments2Json (jRoot, "comments");
    }
    
    if (! jRoot || verbose ())
    {
      cout << endl;
      sp1. printCsv (sm, cout);
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



