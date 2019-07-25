// mergeSimilarity.cpp

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
*   Merge similarities
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/matrix.hpp"
#include "../dm/dataset.hpp"
using namespace DM_sp;



namespace 
{


struct ThisApplication : Application
{
  ThisApplication () 
    : Application ("Print a " + dmSuff + "-file with merged <Similarity> attrributes of a list of " + dmSuff + "-files")
    {
  	  addPositional ("file_list", "List of " + dmSuff + "-files without \"" + dmSuff + "\" with a binary <Similarity> attribute");
  	  addKey ("attr", "Similarity attribute", "Similarity");
  	}



	void body () const final
	{
		const string file_list = getArg ("file_list");
		const string attrName  = getArg ("attr");
		
		
		Set<string> objNames;
    {
      // Pass 1
      LineInput li (file_list);
      while (li. nextLine ())
      {
        const Dataset ds (li. line);
        ds. qc ();
        ASSERT (ds. getUnitMult ());
        for (const Obj* obj : ds. objs)
          objNames << obj->name;
      }
    }

		
    Dataset ds_new;
    for (const string& objName : objNames)
      ds_new. appendObj (objName);
    ds_new. setName2objNum ();
    ASSERT (! ds_new. objs. empty ());
    ASSERT (ds_new. objs. size () == objNames. size ());
    

    auto sim_new = new RealAttr2 (attrName, ds_new, 6);  // PAR
    sim_new->setAll (0);
    Matrix count (false, sim_new->matr, false, 0);
  //count. putAll (0);
    {
      // Pass 2
      LineInput li (file_list);
      Progress prog;
      while (li. nextLine ())
      {
        prog (li. line);
        const Dataset ds (li. line);
        
        const Attr* attr = ds. name2attr (attrName);
        ASSERT (attr);
        const RealAttr2* sim = attr->asRealAttr2 ();
        ASSERT (sim);

        Matrix& matr = const_cast <RealAttr2*> (sim) -> matr;
        ASSERT (matr. defined ());
        {
	        Matrix rowMean;
	        Real totalMean;
	        matr. centerSimilarity (rowMean, totalMean);
	      }
        
        // Normalization
        const Real k = matr. getTrace () / (Real) ds. objs. size ();
        ASSERT (positive (k));
        matr. putProdAll (1.0 / k);
        
        FOR (size_t, row, ds. objs. size ())
        {
          const size_t row_new = ds_new. getName2objNum (ds. objs [row] -> name);
          ASSERT (row_new != NO_INDEX);
          FOR (size_t, col, ds. objs. size ())
          {
            const size_t col_new = ds_new. getName2objNum (ds. objs [col] -> name);
            ASSERT (col_new != NO_INDEX);
            const Real val = matr. get (false, row, col);
            ASSERT (! isNan (val));
            sim_new->matr. putInc (false, row_new, col_new, val);  
            count. putInc (false, row_new, col_new, 1.0);  
          }
        }
      }
    }


    // Averaging
    FOR (size_t, row, ds_new. objs. size ())
      FOR (size_t, col, ds_new. objs. size ())
      {
        const Real val = sim_new->matr. get (false, row, col);
        if (! isNan (val))
          sim_new->matr. put (false, row, col, val / count. get (false, row, col));  
      }


    ds_new. print (cout);
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



