// threshold.cpp

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
*   Find a threshold separating the values of a Boolean attribute
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "numeric.hpp"
#include "dataset.hpp"
using namespace DM_sp;
#include "../version.inc"



namespace
{


#if 0
struct Kernel 
// 1-dimensional
{
  Real x_min, x_max;
    // x_min < x_max
  const Real window;
    // > 0
  Vector<size_t> cells;
private:
  size_t total {0};
public:
  
  
  Kernel (Real x_min_arg, 
          Real x_max_arg,
          const Real window_arg,
          size_t cellsSize)
    : x_min (x_min_arg)
    , x_max (x_max_arg)
    , window (window_arg)
    , cells (cellsSize, 0)
    { ASSERT (x_min < x_max); 
      ASSERT (window > 0.0);
      ASSERT (cellsSize);
    }
  void add (Real x)
    { const long up = min ((long) cells. size () - 1, x2cell (x + window)); 
      for (long i = max ((long) 0, x2cell (x - window)); i <= up; i++)
        cells [(size_t) i] ++; 
    }
  void finish ()
    { total = 0;
      for (const size_t c : cells)
        total += c;
    }
  void saveText (ostream &os) const
    { FOR (size_t, i, cells. size ())
        os << cell2x (i) << '\t' << (Real) cells. at (i) / (Real) total << endl;
    }
    
  
  Real range () const
    { return x_max - x_min; }
  Real step () const
    { return range () / (Real) (cells. size () - 1); }
  Real cell2x (size_t i) const
    { return x_min + (Real) i * step (); }
  long x2cell (Real x) const
    { return DM_sp::round ((x - x_min) / range () * (Real) (cells. size () - 1)); }
};
#endif

  

struct Bin : Root
{
  Real start {NaN};
  Real stop {NaN};
    // start < stop
  size_t classMult {0};
#if 0
  Real mult {0.0};
  Real grayZone_lo {0.0};
  Real grayZone_hi {0.0};
  bool merge {false};
    // With next bin
#endif
  
  explicit Bin (Real stop_arg)
    : stop (stop_arg)
    {}    
  void saveText (ostream &os) const
    { TabDel td;
      td << start << stop << classMult 
       //<< mult << grayZone_lo << grayZone_hi << merge
         ; 
      os << td. str () << endl;
    }

#if 0
  Real range () const
    { return stop - start; }    
  void add (Real score,
            Real objMult)
    { ASSERT (betweenEqualReal (score, start, stop));
      mult += objMult;
      const Real grayLen = range () / 6.0;  // PAR
      ASSERT (grayLen >= 0);
           if (score <= start + grayLen)
        grayZone_lo += objMult;
      else if (score >= stop  - grayLen)
        grayZone_hi += objMult;
    }
#endif
};                   
  


size_t findClass_last (const Vector<Bin> &bins,
                       size_t classMult_max,
                       Prob FP_fraction)
{
  size_t class_last = NO_INDEX;
  FOR (size_t, i, bins. size ())
    if ((Real) bins. at (i). classMult >= FP_fraction * (Real) classMult_max)
      class_last = i;
  ASSERT (class_last != NO_INDEX);
  return class_last;
}

  

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Find a <score-attr> threshold separating the values of Boolean attribute <class_attr>")
    {
      version = VERSION;
      addPositional ("file", dmSuff + "-file, sorted by <score_attr> descending");
      addPositional ("score_attr", "Score number attribute");
      addPositional ("class_attr", "Class Boolean attribute");
      addPositional ("sd_min", "Min. SD of each variable in each cluster");
    }
	
	
	
	void body () const final
	{
  //const Prob confused_max   = args["confused_max"].AsDouble();
		const string inFName    = getArg ("file");
		const string score_attr = getArg ("score_attr");
		const string class_attr = getArg ("class_attr");
		const Real sd_min       = str2real (getArg ("sd_min"));
	//ASSERT (isProb (confused_max));


    Dataset ds (inFName);

    const RealAttr1* scoreAttr = checkPtr (ds. name2attr (score_attr)) -> asRealAttr1 ();
    ASSERT (scoreAttr);
    const BoolAttr1* classAttr = checkPtr (ds. name2attr (class_attr)) -> asBoolAttr1 ();
    ASSERT (classAttr);

    Space1<NumAttr1> sp (ds, false);
    sp << scoreAttr;
    const Sample sm (ds);
        
	  const Clustering cl (sm, sp, 20, sd_min, false);  // PAR  entropyDimensionPrecision = 0: may hang
    cl. qc ();

    NominAttr1* clusterAttr = cl. createNominAttr ("Cluster", 0, ds);
    ASSERT (clusterAttr);
  #if 0
    ProbAttr1* probAttr = cl. createProbAttr ("Cluster_prob", 2, ds);  // PAR
    ASSERT (probAttr);
    cl. mergeClose (*clusterAttr, *probAttr, confused_max);  
  #endif
    
    const size_t clusters = clusterAttr->categories. size ();  
    ASSERT (clusters <= cl. getOutDim ());
	  cout << "# Clusters: " << clusters << endl;
	  if (clusters < 2)
	    return;

  #if 0
    {
      Real score_max = -INF;	    
      for (Iterator it (cl); it ();)  
        if (! scoreAttr->isMissing (*it))
          maximize (score_max, (*scoreAttr) [*it]);
      ASSERT (DM_sp::finite (score_max));
  
  	  Kernel kern (0, score_max, 20, 40);  // PAR
      for (Iterator it (cl); it ();)  
        if (! scoreAttr->isMissing (*it))
          kern. add ((*scoreAttr) [*it]);  // mult ??
      kern. finish ();
      cout << endl;
      kern. saveText (cout);
      cout << endl;
    }
  #endif

    // There can be more bin's than there really are due to background clusters
    Vector<Bin> bins;  bins. reserve (clusters);  // Ordered by score descending
    Real score_prev = NaN;
    size_t cluster_prev = NominAttr1::missing;
    for (Iterator it (sm); it ();)  
      if (! clusterAttr->isMissing (*it))
      {
        const Real score = (*scoreAttr) [*it];
        IMPLY (! isNan (score_prev), geReal (score_prev, score));
        const size_t cluster = (*clusterAttr) [*it];
        if (! isNan (score_prev))
          bins. back (). start = score_prev;
        if (cluster != cluster_prev)
          bins << Bin (score);
        if (classAttr->getBool (*it) == ETRUE)
          bins. back (). classMult++;
        score_prev = score;
        cluster_prev = cluster;
      }
    if (! isNan (score_prev))
      bins. back (). start = score_prev;
  //ASSERT (bins. size () >= clusters);  
     // Clusters with a large variance may be in >= bin's
     // Clusters with a small probabolity may be not coverde by a bin


  #if 0
    auto binIt = bins. begin ();
    for (Iterator it (sm); it ();)  
      if (! clusterAttr->isMissing (*it))
      {
        const Real score = (*scoreAttr) [*it];
        if (score < binIt->start)
          binIt++;
        ASSERT (binIt != bins. end ());
        binIt->add (score, it. mult);
      }
      
    ITER (Vector<Bin>, it, bins)
    {
      const auto &nextIt = next (it);
      if (nextIt == bins. end ())
        break;
      if (it->grayZone_lo + nextIt->grayZone_hi >= min (it->mult, nextIt->mult) * 0.2)  // PAR
        it->merge = true;
    }
  #endif


    if (verbose ())
    {
      {
        TabDel td;
        td << "Start" << "End" << "ClassMult" << "mult" << "grayZone_lo" << "grayZone_hi" << "merge";
        cout << td. str () << endl;
      }
      for (const Bin& bin : bins)
        bin. print (cout);
      cout << endl;

      cl. print (cout);
      
      VectorPtr<Attr> attrs;
      attrs << clusterAttr << scoreAttr << classAttr;
      cout << endl;
      sm. save (attrs, cout);      

      cout << endl;
      Real score_min, score_max;
      scoreAttr->getMinMax (sm, score_min, score_max);
      Histogram h (0, score_max, 5);  // PAR
      for (Iterator it (sm); it ();)  
        h << (*scoreAttr) [*it];
      h. print (cout);
    }

 
    size_t classMult_max = 0;
    for (const Bin& bin : bins)
      maximize (classMult_max, bin. classMult);

    size_t class_last = NO_INDEX;
    for (Prob FP_fraction = 0.01; FP_fraction <= 0.2; FP_fraction += 0.01)  // PAR
    {
      class_last = findClass_last (bins, classMult_max, FP_fraction);  
      if (class_last < bins. size () - 1)
        break;
    }
    ASSERT (class_last != NO_INDEX);
              
    // It is assumed that the class objects are well separated from the others
    const size_t other_first = class_last + 1;
    
    if (other_first == bins. size ())
      cout << "Cannot separate bins" << endl;
    else
    {
      const Real start = bins. at (class_last).  start;
      const Real stop  = bins. at (other_first). stop;
      const Real gap = start - stop;
      cout << "Gap: " << gap << endl;
      cout << "Threshold: " << start - gap / 2 << endl;
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


