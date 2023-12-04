// tsv_cluster.cpp

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
*   Cluster a text column of a tsv-table
*
*/

#undef NDEBUG

#include "../common.hpp"
#include "tsv.hpp"
using namespace Common_sp;
#include "../dissim/evolution.hpp"
#include "../version.inc"

#include "../common.inc"



namespace
{
  

const string suf ("_cluster");


struct TextClust : DisjointCluster
{
  const size_t elem;
  Set<string> keys;
  
  explicit TextClust (size_t elem_arg)
    : elem (elem_arg)
    {}
};
  
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Cluster a text column of a tsv-table base on k-mers. Add numeric column <col>" + suf)
  	{
      version = VERSION;
  	  addPositional ("table", "tsv-table");
  	  addPositional ("col", "Text column");
  	  addKey ("k", "K-mer size", "20");
  	  addKey ("dissim", "Max. dissimlarity for single-linkage clustering", "1.0");  
  	}
  	
  	
 
	void body () const final
	{
		const string fName      =               getArg ("table");
		const string colName    =               getArg ("col");
		const size_t k          = str2<size_t> (getArg ("k"));		
	  const Real   dissim_max = str2real     (getArg ("dissim"));	
	  

    TextTable tt (fName);
    tt. qc ();
    
    const TextTable::ColNum col = tt. col2num (colName);
    if (tt. header [col]. numeric)
      throw runtime_error (strQuote (colName) + " is not a text column");
    
    const string outColName = colName + suf;
    if ( tt. hasColumn (outColName))
      throw runtime_error ("Table already has column " + strQuote (outColName));
    const size_t outCol = tt. header. size ();
    tt. header << TextTable::Header (outColName);
    tt. header [outCol]. numeric = true;
    
    VectorOwn<TextClust> clusters;  clusters. reserve (tt. rows. size ());
    FFOR (size_t, i, tt. rows. size ())
    {
      auto tc = new TextClust (i);
      const StringVector& row = tt. rows [i];
      const string& text = row [col];
      if (! text. empty ())
      {
        if (text. size () < k)
          tc->keys << text;
        else
          FFOR (size_t, j, text. size () - k)
            tc->keys << text. substr (j, k);
      }
      clusters << tc;
      ASSERT (clusters [i] -> elem == i);
    }
    ASSERT (clusters. size () == tt. rows. size ());
    
    {
      Progress prog (clusters. size () * (clusters. size () - 1) / 2, 1000);  // PAR
      for (const auto& it1 : clusters)
        for (const auto& it2 : clusters)
          if (& it1 == & it2)
            break;
          else
          {
            prog ();
            ASSERT (it1->elem > it2->elem);
            const Real dissim = DM_sp::intersection2dissim ( (Real) it1->keys. size ()
                                                           , (Real) it2->keys. size ()
                                                           , (Real) it1->keys. intersectSize (it2->keys)
                                                           , 1.0
                                                           , 0.5
                                                           , true
                                                           );
          //cout << dissim << endl;  // to determine dissim_max ??
            if (dissim <= dissim_max)
              var_cast (it1) -> merge (var_cast (*it2));
          }
    }
    
    FFOR (size_t, i, tt. rows. size ())
    {
      StringVector& row = tt. rows [i];
      row << to_string (static_cast <const TextClust*> (var_cast (clusters [i]) -> getDisjointCluster ()) -> elem + 1);
    }
    tt. qc ();
        
    tt. saveText (cout);  
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



