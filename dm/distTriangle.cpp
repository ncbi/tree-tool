// distTriangle.cpp

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
*   Find violations of triangle inequality
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "numeric.hpp"
#include "matrix.hpp"
#include "dataset.hpp"
using namespace DM_sp;



namespace 
{


struct ObjNeighbors
{
  size_t objNum {0};
  size_t neighbors {0};
  
  ObjNeighbors (size_t objNum_arg,
                size_t neighbors_arg)
    : objNum (objNum_arg)
    , neighbors (neighbors_arg)
    {}
    
  bool operator< (const ObjNeighbors &other) const
    { LESS_PART (other, *this, neighbors);
      return objNum < other. objNum;
    }
};



struct Violation
// Of triangle inquality
{
  // Dataset::objs index
  size_t x {0};
  size_t y {0};
  size_t z {0};
    // middle object
  Real hybridness {NaN};
  
  Violation () = default;
  Violation (size_t x_arg,
             size_t y_arg, 
             size_t z_arg,
             Real hybridness_arg)
    : x (x_arg)
    , y (y_arg)
    , z (z_arg)
    , hybridness (hybridness_arg)
    { ASSERT (x != y);
      ASSERT (x != z);
      ASSERT (y != z);
      ASSERT (hybridness > 1.0);
    }
  void print (const Dataset &ds) const
    { cout        << ds. objs [z] -> name 
           << ' ' << ds. objs [x] -> name 
           << ' ' << ds. objs [y] -> name;
    }
    
  bool contains (size_t i) const
    { return    i == x
             || i == y
             || i == z;
    }
};



struct Violator : Named
{
  // # Violation's
  size_t end {0};
  size_t middle {0};
  
  Violator (const string &name_arg,
            size_t end_arg,
            size_t middle_arg)
    : Named (name_arg)
    , end (end_arg)
    , middle (middle_arg)
    { ASSERT (num ()); }
  Violator () = default;
  void saveText (ostream &os) const
    { os << name << ": " << end << "+" << middle << " = " << num (); }
    
  size_t num () const
    { return end + middle; }
};




struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Find clusters of objects connected by finite distances.\n\
Remove distance outliers by the triangle inequality analysis.\n\
Save clusters in " + dmSuff + "-files\
")
  	{
  		// Input
  	  addPositional ("file", dmSuff + "-file without the extension");
  	  addPositional ("attrName", "Attribute name of a distance in the <file>");
  	  addFlag("max_cliques", "Clusters must be maximal cliques");
  	  addKey ("hybridness_min", "Min. hybridness to report hybrids, >1", "NAN");
  	  addKey ("hybrid", "Output file with hybrids information. Line format: " + string (PositiveAttr2::hybrid_format));
  	  addKey ("distance_max", "Max. distance to merge into the same cluster; 0 - infinity", "0");
  	  addKey ("objName", "Object name to print violations for");
  	  // Output
  	  addKey ("clustering_dir", "If specified then save the data of each non-singleton cluster in this directory");  	  
  	}



	void body () const
	{
		const string clustering_dir = getArg ("clustering_dir");
		const bool   max_cliques    = getFlag ("max_cliques");
		const string fName          = getArg ("file");
		const string attrName       = getArg ("attrName");
	  const Real   hybridness_min = str2real (getArg ("hybridness_min"));
	  const string hybridFName   = getArg ("hybrid");
		const Real   distance_max   = str2real (getArg ("distance_max"));
		const string objName        = getArg ("objName");
	  if (hybridness_min <= 1.0)
	    throw runtime_error ("-hybridness_min must be > 1.0");
		
		
    Dataset ds (fName);
    
    // dist
    const PositiveAttr2* dist = nullptr;
    {
      const Attr* attr = ds. name2attr (attrName);
      ASSERT (attr);
      dist = attr->asPositiveAttr2 ();
    }
    ASSERT (dist);
    

	  // Check dist->matr  
	  {
      Real maxCorrection;
      size_t row_bad, col_bad;   
      var_cast (dist) -> matr. symmetrize (maxCorrection, row_bad, col_bad);
      if (maxCorrection > 2 * pow (10, - (Real) dist->decimals))
        cout << "maxCorrection = " << maxCorrection << " at " << ds. objs [row_bad] -> name << ", " << ds. objs [col_bad] -> name << endl;
    }
  #if 0
    {
   	  size_t row, col;
  	  if (dist->matr. existsMissing (false, row, col))
      {
      #if 0
      	var_cast (dist->matr). put (false, row, col, INF);
      #else
        cout << dist->name << " [" << row + 1 << "] [" << col + 1 << "] is missing" << endl;
        exit (1);
      #endif
      }
    }
  #endif
	  {
  	  size_t row;
      if (! dist->matr. zeroDiagonal (row))
      {
      #if 0
      	var_cast (dist->matr). put (false, row, row, 0);
      #else
       	cout << dist->name << " [" << row + 1 << "] [" << row + 1 << "] = " << dist->matr. get (false, row, row) << " != 0" << endl;
       	exit (1);
      #endif
      }
    }
    
    size_t objIndex = NO_INDEX;
    if (! objName. empty ())
    { 
      objIndex = ds. getName2objNum (objName);
      ASSERT (objIndex != NO_INDEX);
    }

    // ave
    MeanVar distMV;
    FOR (size_t, row, ds. objs. size ())
      FOR (size_t, col, ds. objs. size ())
      {
        const Real r = dist->get (row, col);
        if (   DM_sp::finite (r)
            && row != col
           )
          distMV << r;
      }
    const Real ave = distMV. getMean ();
    cout << "Average distance = " << ave << endl; 


    unordered_map <const DisjointCluster*, Vector<size_t/*objNum*/>> cluster2objs;  cluster2objs. rehash (ds. objs. size ());
    if (max_cliques)
    {
      // Greedy algorithm
      Vector<bool> clustered (ds. objs. size (), false);
      for (;;)
      {
        Vector<ObjNeighbors> objNeighbors;  objNeighbors. reserve (ds. objs. size());
        FOR (size_t, row, ds. objs. size ())
          if (! clustered [row])
          {
            size_t n = 0;
            FOR (size_t, col, ds. objs. size ())
              if (! clustered [col])
                if (dist->areClose (row, col, distance_max))
                  n++;
            objNeighbors << ObjNeighbors (row, n);
          }
        ASSERT (objNeighbors. size () <= ds. objs. size ());
        objNeighbors. sort ();
        Obj* mainObj = nullptr;
        for (const ObjNeighbors on : objNeighbors)
          if (mainObj)
          {
            ASSERT (! cluster2objs [mainObj]. empty ());
            bool inClique = true;
            for (const size_t objNum : cluster2objs [mainObj])
              if (! dist->areClose (on. objNum, objNum, distance_max))
              {
                inClique = false;
                break;
              }
            if (inClique)
              cluster2objs [mainObj] << on. objNum;
          }
          else
          {
            mainObj = const_cast <Obj*> (ds. objs. at (on. objNum));
            ASSERT (cluster2objs [mainObj]. empty ());
            cluster2objs [mainObj] << on. objNum;
          }
        if (! mainObj)
          break;
        for (const size_t objNum : cluster2objs [mainObj])
          clustered [objNum] = true;
      }
      ASSERT (! clustered. contains (false));
    }
    else
    {
      // Obj::DisjointCluster
     	for (const Obj* obj : ds. objs)
   	    var_cast (obj) -> DisjointCluster::init ();
      FOR (size_t, row, ds. objs. size ())
        FOR (size_t, col, ds. objs. size ())
          if (dist->areClose (row, col, distance_max))
            var_cast (ds. objs [row]) -> merge (* const_cast <Obj*> (ds. objs [col]));
      // cluster2objs
      FOR (size_t, i, ds. objs. size ())
        cluster2objs [var_cast (ds. objs [i]) -> getDisjointCluster ()] << i;
    }
    cout << "# Clusters = " << cluster2objs. size () << endl;
    
    
    unique_ptr<OFStream> hybridF;
    if (! hybridFName. empty ())
      hybridF. reset (new OFStream (hybridFName));
    size_t nClust = 0;
    for (const auto& clustIt : cluster2objs)
    {
      const Vector<size_t>& objs = clustIt. second;
      ASSERT (! objs. empty ());
      if (objs. size () == 1)
      {
        cout << endl << "Singleton: " << ds. objs. at (objs. front () ) -> name << endl;
        continue;  
      }
        
      nClust++;
      cout << endl << "Cluster #" << nClust << ":  size=" << objs. size () << endl;   
      Set<size_t> objSet;
      insertIter (objSet, objs);
      ASSERT (objSet. size () == objs. size ());

      // Distribuiton of distances: exponential ??
  
      // Triangle inequality
      List<Violation> violations;
      {
        Progress prog (objs. size ());
        FOR_REV (size_t, y_, objs. size ())
        {
          prog ();
          FOR (size_t, x_, y_)
            FOR (size_t, z_, objs. size ())
              if (   z_ != x_
                  && z_ != y_
                 )
              {
                const size_t x = objs [x_];
                const size_t y = objs [y_];
                const size_t z = objs [z_];
                const Real d = dist->get (x, y); 
                if (isNan (d))
                  continue;
                ASSERT (d >= 0.0);
                const Real dPair = dist->get (x, z) + dist->get (z, y);
                if (isNan (dPair))
                  continue;
                ASSERT (dPair >= 0.0);
                const Real hybridness = d / dPair;
                if (hybridness < hybridness_min)
                  continue;
                const Violation violation (x, y, z, hybridness);
                violations << violation;
                if (   violation. contains (objIndex)
                    || verbose ()
                   )
                {
                  const string xName (ds. objs [x] -> name);
                  const string yName (ds. objs [y] -> name);
                  const string zName (ds. objs [z] -> name);
                  cout << "d(" << xName << "," << yName << ") = " << d 
                           << " > d(" << xName << "," << zName << ") + d(" << zName << "," << yName << ") = " << dPair
                         //<< " by " << deviation << " (" << fraction * 100.0 << " %)"
                           << " hybridness = " << hybridness
                           << endl;
                }
              }
        }
      }
      cout << "# Violations = " << violations. size () << endl;
          
      // Set-covering problem: approximation alogorithm
      const size_t violationsNum = violations. size ();
      size_t coverage = 0;
      Vector<size_t> endViolations    (ds. objs. size ());
      Vector<size_t> middleViolations (ds. objs. size ());
      while (! violations. empty ())
      {
        size_t i_best = NO_INDEX;
        {
          endViolations.    setAll (0);
          middleViolations. setAll (0);
          for (const Violation& v : violations)
          {
            endViolations    [v. x] ++;
            endViolations    [v. y] ++;
            middleViolations [v. z] ++;
          }
          size_t a = 0;
          FOR (size_t, i, endViolations. size ())
            if (maximize (a, endViolations [i] + middleViolations [i]))
              i_best = i;
        }
        ASSERT (i_best != NO_INDEX);
  
        const Violator v (ds. objs [i_best] -> name, endViolations [i_best], middleViolations [i_best]);
        coverage += v. num ();

        Real hybridness_max = 0.0;
        Violation violation_worst;
        for (Iter <List<Violation> > iter (violations); iter. next (); )
          if (iter->contains (i_best))
          {
            if (maximize (hybridness_max, iter->hybridness))
              violation_worst = *iter;
            iter. erase ();
          }
        ASSERT (hybridness_max > 1.0);

        if (hybridF. get ())
          *hybridF         << ds. objs [violation_worst. z] -> name
                    << '\t' << violation_worst. hybridness
                    << '\t' << ds. objs [violation_worst. x] -> name
                    << '\t' << ds. objs [violation_worst. y] -> name
                    << '\t' << dist->get (violation_worst. z, violation_worst. x)
                    << '\t' << dist->get (violation_worst. z, violation_worst. y)
                    << '\t' << (ds. objs [violation_worst. z] -> name == v. name)
                    << '\t' << (ds. objs [violation_worst. x] -> name == v. name)
                    << '\t' << (ds. objs [violation_worst. y] -> name == v. name)
                    << endl;

        v. saveText (cout);
        const Prob coverageFrac = (Real) coverage / (Real) violationsNum;
        cout << "  coverage = " << coverageFrac * 100.0 << " %";
        cout << "  hybridness = " << hybridness_max;
        cout << "  ";
        violation_worst. print (ds);
        cout << endl;

        objSet. erase (i_best);  // => redo the clustering !??
      }

      if (! clustering_dir. empty ())
      {
        const string dir (clustering_dir + "/" + toString (nClust));
        exec ("mkdir " + dir);
        OFStream f (dir, getFileName (fName), dmExt);

        const VectorPtr<Attr> attrs (ds. attrs);

        Sample sm (ds);
        FOR (size_t, row, ds. objs. size ())
          if (! objSet. contains (row))
            sm. mult [row] = 0.0;
        sm. finish ();    
        sm. save (attrs, f);
      }
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



