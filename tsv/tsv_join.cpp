// tsv_join.cpp

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
*   Join a tsv-table with a list
*
*/

#undef NDEBUG

#include "../common.hpp"
#include "tsv.hpp"
using namespace Common_sp;
#include "../version.inc"

#include "../common.inc"



namespace
{
  
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Join 2 tsv-tables by identical columns, print the result")
  	{
      version = VERSION;
  	  addPositional ("table1", "tsv-table");
  	  addPositional ("table2", "tsv-table");
  	  addKey ("syn", TextTable::syn_format);
  	  addFlag ("left", "SQL left join");
  	  addFlag ("remove", "Remove the rows from <table1> which are in <table2>");
  	  addKey ("common", "Output file to print common columns");
  	}
  	
  	
 
	void body () const final
	{
		const string tab1FName   = getArg ("table1");
		const string tab2FName   = getArg ("table2");
		const string synFName    = getArg ("syn");
		const bool   leftjoin    = getFlag ("left");
    const bool   remove      = getFlag ("remove");  
    const string commonFName = getArg ("common");
    
    QC_IMPLY (remove, ! leftjoin);
		

    const TextTable t1 (tab1FName, synFName);
    t1. qc ();    
    
    const TextTable t2 (tab2FName, synFName);
    t2. qc ();    
    
    StringVector commonCols;
    {
      for (const TextTable::Header& h : t1. header)
        commonCols << h. name;
      FOR_REV (size_t, i, commonCols. size ())
        if (! t2. hasColumn (commonCols [i]))
          commonCols. eraseAt (i);
      {
        StringVector s (commonCols);
        s. sort ();
        const size_t i = s. findDuplicate ();
        if (i != no_index)
          throw runtime_error ("Duplicate common column in " + strQuote (tab1FName) + ": " + s [i]);
        StringVector s2;
        for (const TextTable::Header& h : t2. header)
          if (s. containsFast (h. name))
            s2 << h. name;
        for (const string& name : s)
          if (s2. countValue (name) > 1)
            throw runtime_error ("Duplicate common column in " + strQuote (tab2FName) + ": " + name);
      }
      if (commonFName. empty ())
      {
        cerr << "Common columns:" << endl;
        save (cerr, commonCols, '\n');
        cerr << endl;
      }
      else
      {
        OFStream f (commonFName);
        save (f, commonCols, '\n');
        f << endl;
      }
    }
    ASSERT (commonCols. size () <= t1. header. size ());
    ASSERT (commonCols. size () <= t2. header. size ());
      
    const TextTable::Index index (t2, commonCols);        


    TextTable tOut (t1. pound, t1. header);
    tOut. name = "Output";
    {
      Vector<TextTable::ColNum> addedColNums;
      if (! remove)
        FOR (TextTable::ColNum, colNum, t2. header. size ())
        {
          const TextTable::Header& h = t2. header [colNum];
          if (! commonCols. contains (h. name))
          {
            tOut. header << h;
            addedColNums << colNum;
          }
        }

      const Vector<TextTable::ColNum> indexColNums1 (t1. columns2nums (commonCols));     
      StringVector indexValues;
      StringVector values2;
      StringVector outValues;
      FOR (TextTable::RowNum, rowNum1, t1. rows. size ())
      {
        t1. colNumsRow2values (indexColNums1, rowNum1, indexValues);
        for (const string& s : indexValues)
          if (s. empty ())
            throw TextTable::Error (t1, "Empty value in index, in row " + to_string (rowNum1 + 1));
        if (const Vector<TextTable::RowNum>* rowNums2 = index. find (indexValues))
        {
          ASSERT (! rowNums2->empty ());
          if (! remove)
            for (TextTable::RowNum rowNum2 : *rowNums2)
            {
              t2. colNumsRow2values (addedColNums, rowNum2, values2);
              outValues. clear ();
              outValues << t1. rows [rowNum1] << values2;
              tOut. rows << std::move (outValues);
            }
        }
        else 
        {
          if (remove || leftjoin)
          {
            tOut. rows << std::move (t1. rows [rowNum1]);
            if (leftjoin)
              tOut. rows. back (). resize (t1. header. size () + t2. header. size () - commonCols. size ());
          }
        }
      }
    }
    tOut. qc ();    
    tOut. saveText (cout);
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



