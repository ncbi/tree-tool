// feature2dissim.cpp

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
*   Convert feature sets to a Jaccard conservation distance
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/dataset.hpp"
using namespace DM_sp;
#include "version.inc"



namespace 
{


const string dissimName = "cons";



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Convert feature sets to a Jaccard conservation distance " + strQuote (dissimName) + " and print a " + dmSuff + "-file")
    {
      version = VERSION;
  	  addPositional ("objects", "File with a list of object files with features");
  	  addPositional ("objects_dir", "Directory with <object> files containing features");
  	}



	void body () const
	{
		const string objectsFName = getArg ("objects");
		const string objects_dir  = getArg ("objects_dir");
		
		
		typedef  StringVector Features;
		
    Dataset ds;
    Vector<Features> obj2features;  

    size_t features_min = numeric_limits<size_t>::max();
    size_t features_max = 0;
    {
      LineInput objF (objectsFName);
      while (objF. nextLine ())
      {
        trim (objF. line);
        ds. appendObj (objF. line);
        Features features;
        {
          LineInput featF (objects_dir + "/" + objF. line);
          while (featF. nextLine ())
          {
            trim (featF. line);
            features << featF. line;
          }
        }
        features. sort ();
        features. uniq ();
        obj2features << features;
        minimize (features_min, features. size  ());
        maximize (features_max, features. size  ());
      }
    }
    ASSERT (ds. objs. size () == obj2features. size ());
    cerr << "# Objects: " << ds. objs. size () << endl;
    cerr << "Min. features: " << features_min << endl;
    cerr << "Max. features: " << features_max << endl;
    ASSERT (features_min);


    auto attr = new PositiveAttr2 (dissimName, ds, 6);  // PAR
    Progress prog (ds. objs. size ());
    FOR (size_t, i, ds. objs. size ())
    {
      prog ();
      const Features& f1 = obj2features [i];
      FOR_START (size_t, j, i + 1, ds. objs. size ())
      {
        const Features& f2 = obj2features [j];
        const size_t intersection = f1. getIntersectSize (f2);
        Features f (f2);
        f << f1;
        f. sort ();
        f. uniq ();
        const Prob jaccard = (Real) intersection / (Real) f. size ();
        ASSERT (isProb (jaccard));
        const Real dissim = - log (jaccard);
        attr->putSymm (i, j, dissim);
      }
    }
    FOR (size_t, row, ds. objs. size ())
      attr->put (row, row, 0);

    
    ds. qc ();
    ds. saveText (cout);    
	}
};



}  // namespace



int main(int argc, 
         const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



