// tsv2dm.cpp

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
*   Convert a .tsv-table into .dm-format
*
*/

#undef NDEBUG

#include "../../common.hpp"
#include "../../tsv/tsv.hpp"
using namespace Common_sp;
#include "../dataset.hpp"
using namespace DM_sp;
#include "../../version.inc"

#include "../../common.inc"



namespace
{
  


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print a .tsv-table in .dm-format")
  	{
      version = VERSION;
  	  addPositional ("in", "tsv-table");
  	  addPositional ("col", "Object name column");
  	  addKey ("missing", "Replacement for missing values", missingStr);
  	}
  	
  	
 
	void body () const final
	{
		const string inFName    = getArg ("in");
		const string objColName = getArg ("col");
		const string missingVal = getArg ("missing");
	  

    const TextTable tt (inFName);
    tt. qc ();
    
    const TextTable::ColNum objCol = tt. col2num (objColName);
    
    Dataset ds;
    FFOR (TextTable::RowNum, i, tt. rows. size ())
    {
      string name (tt. rows [i] [objCol]);
      replace (name, ' ', '_');
      ds. appendObj (name);
    }
    ds. setName2objNum ();
    
    VectorPtr<Attr1> attrs;  attrs. reserve (tt. header. size ());
    FFOR (TextTable::ColNum, col, tt. header. size ())
    {
      Attr1* attr = nullptr;
      if (col != objCol)
      {
        const TextTable::Header& h = tt. header [col];
        string name (h. name);
        replace (name, ' ', '_');
        if (h. numeric)
          if (h. scientific || h. decimals)
            attr = new RealAttr1 (name, ds, h. decimals);
          else
            attr = new IntAttr1 (name, ds);
        else
          if (h. choices. size () <= TextTable::Header::choices_max)
            attr = new NominAttr1 (name, ds);
      }
      attrs << attr;
    }
    ASSERT (attrs. size () == tt. header. size ());
    ds. qc ();

    FFOR (TextTable::RowNum, i, tt. rows. size ())
    {
      const StringVector& row = tt. rows [i];
    //ASSERT (ds. getName2objNum (row [objCol]) == i);  // ' ' -> '_'
      FFOR (TextTable::ColNum, col, tt. header. size ())
        if (const Attr1* attr = attrs [col])
        {
          if (row [col]. empty ())
          {
            if (missingVal != missingStr)
              var_cast (attr) -> str2value (i, missingVal);
          }
          else
            var_cast (attr) -> str2value (i, row [col]);
        }
    }

    ds. qc ();
    ds. saveText (cout);    
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



