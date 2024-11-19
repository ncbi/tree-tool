// tsv_comp.cpp  

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
*   Comparison of 2 .tsv-files
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
  

Vector<TextTable::ColNum> getCommonCols (const TextTable &t1,
                                         const TextTable &t2)
// Return: ordered by TextTable::Header::name
// Print: column names of t1 missing in t2
{
  StringVector res;  res. reserve (t1. header. size ());
  for (const TextTable::Header& h1 : t1. header)
    if (t2. col2num_ (h1. name) == no_index)
      cout << h1. name << endl;
    else
      res << h1. name;
  cout << endl;
  
  res. sort ();
  QC_ASSERT (res. isUniq ());

  return t1. columns2nums (res);
}
  
  
typedef  map<string/*key*/, const StringVector* /*row*/>  Key2row;
  
  

Key2row makeKey2row (const TextTable &tab,
                     const Vector<TextTable::ColNum> &keyCols)
{
  Key2row key2row;
  for (const StringVector& row : tab. rows)
  {
    string key;
    bool first = true;
    for (const TextTable::ColNum i : keyCols)
    {
      if (first)
        first = false;
      else
        key += '\t';
      key += row [i];
      ASSERT (! contains (row [i], '\t'));
    }
    if (key2row [key])
      throw runtime_error (tab. name + ": duplicate key: " + key);
    key2row [key] = & row;
  }
  return key2row;
}
  
  

struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Compare two .tsv-files")
  	{
      version = VERSION;
  	  addPositional ("in1", ".tsv-file 1");
  	  addPositional ("in2", ".tsv-file 2");
  	  addPositional ("keys", "Key columns, comma-separated");
  	  addPositional ("diff", "Output file with different column values");
  	}



	void body () const final
  {
	  const string fName1    = getArg ("in1");
	  const string fName2    = getArg ("in2");
	  const string keysS     = getArg ("keys");
	  const string diffFName = getArg ("diff");
	  
	  	  
    const TextTable t1 (fName1, noString, 10000);  // PAR
    cout << "# " << fName1 << " rows: " << t1. rows. size () << endl;
    t1. qc ();

    const TextTable t2 (fName2, noString, 10000);  // PAR
    cout << "# " << fName2 << " rows: " << t2. rows. size () << endl;
    t2. qc ();

    const StringVector keysVec (keysS, ',', true);
    const Vector<TextTable::ColNum> keyCols1 (t1. columns2nums (keysVec));
    const Vector<TextTable::ColNum> keyCols2 (t2. columns2nums (keysVec));
    ASSERT (keyCols1. size () == keyCols2. size ());
    
    cout << "Added columns:";
    const Vector<TextTable::ColNum> commonCols1 (getCommonCols (t1, t2));
    cout << "Deleted columns:";
    const Vector<TextTable::ColNum> commonCols2 (getCommonCols (t2, t1));
    ASSERT (commonCols1. size () == commonCols2. size ());
        
    const Key2row key2row1 (makeKey2row (t1, keyCols1));
    const Key2row key2row2 (makeKey2row (t2, keyCols2));
    ASSERT (key2row1. size () == t1. rows. size ());
    ASSERT (key2row2. size () == t2. rows. size ());
    
    size_t added = 0;
    for (const auto& it : key2row2)
      if (! findPtr (key2row1, it. first))
        added++;
    cout << "# Rows ddded: " << added << endl;

    OFStream diffF (diffFName);    
    size_t deleted = 0;
    for (const auto& it : key2row1)
      if (const StringVector* row2 = findPtr (key2row2, it. first))
      {
        const StringVector* row1 = it. second;
        ASSERT (row1);
        FFOR (size_t, i, commonCols1. size ())
        {
          const string& val1 = (*row1) [commonCols1 [i]];
          const string& val2 = (*row2) [commonCols2 [i]];
          if (val1 != val2)
            diffF << t1. header [commonCols1 [i]]. name << '\t' << val1 << '\t' << val2 << '\n';
        }
      }
      else
        deleted++;
    cout << "# Rows deleted: " << deleted << endl;
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


