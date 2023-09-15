// tsv_comp.cpp  // --> textTab_comp.cpp ??

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
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../version.inc"



struct Tsv : Root
// --> TextTable ??
{
  StringVector header;
    // empty() <=> no file
  bool poundSign {false};
  Vector<size_t> keySchema;
  Vector<StringVector> rows;
    // Element. size(0 = headrer.size()
    
    
  Tsv (const string &fName,
       const StringVector& keyNames)
    {
      LineInput f (fName);
      if (! f. nextLine ())
        return;
      
      if (keyNames. empty ())
        throw runtime_error ("Key is empty");
      
      try 
      {
        // header
        header = move (StringVector (f. line, '\t', true));
        if (header. empty ())
          return;
        FFOR (size_t, i, header. size ())
          if (header [i]. empty ())
            throw runtime_error ("Header name is empty for column " + to_string (i + 1));
        if (header [0] [0] == '#')
        {
          poundSign = true;
          header [0]. erase (0, 1);
        }
        {
          Set<string> headerSet;
          for (const string& s : header)
            if (! headerSet. insert (s). second)
              throw runtime_error ("Duplicate header name: " + strQuote (s));
        }
        
        // keySchema
        for (const string& s : keyNames)
        {
          const size_t index = header. indexOf (s);
          if (index == no_index)
            throw runtime_error ("Key column " + strQuote (s) + " does not exist");
          keySchema << index;
        }
        
        // rows
        while (f. nextLine ())
        {
          StringVector row (f. line, '\t', true);
          if (row. size () != header. size ())
            throw runtime_error ("Row has " + to_string (row. size ()) + " columns whereas header has " + to_string (header. size ()) + " columns");
          rows << move (row);
          ASSERT (row. empty ());
        }
      }
      catch (const exception &e)
      {
        throw runtime_error (".tsv-file " + strQuote (fName) + ", row " + to_string (f. tp. lineNum + 1) + ":\n" + e. what ());
      }
    }
    
    
  StringVector row2key (size_t rowNum) const
    {
      StringVector key;  key. reserve (keySchema. size ());
      for (const size_t i : keySchema)
        key << rows [rowNum] [i];
      return key;
    }
    
  typedef  map<StringVector, size_t/*rowNum*/>  Key2row;
  Key2row getKey2row () const
    { 
      Key2row key2row;
      FFOR (size_t, i, rows. size ())
      {
        const StringVector key (row2key (i));
        if (! key2row. insert (Key2row::value_type (key, i)). second)
          throw runtime_error ("Duplicate key: " + key. toString (","));
      }
      return key2row;
    }
};



namespace 
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Compare two .tsv-files")
  	{
      version = VERSION;
  	  addPositional ("in1", ".tsv-file 1");
  	  addPositional ("in2", ".tsv-file 2");
  	  addPositional ("key", "File with key columns");
  	}



	void body () const final
  {
	  const string fName1   = getArg ("in1");
	  const string fName2   = getArg ("in2");
	  const string keyFName = getArg ("key");
	  
	  
	  const StringVector keyNames (keyFName, (size_t) 10, true);
	  
	  
	  const Tsv f1 (fName1, keyNames);
	  const Tsv f2 (fName2, keyNames);
	  
	  // header
	  if (f1. header != f2. header)
	    throw runtime_error ("Files have diffeernt headers");
	  if (f1. header. empty ())
	    return;
	    
	  // Key2row
	  const Tsv::Key2row key2row1 (f1. getKey2row ());
	  const Tsv::Key2row key2row2 (f2. getKey2row ());
	  	  
	  size_t deleted = 0;
	  map<string/*column name*/,size_t> changes;
	  for (const auto& it : key2row1)
	  {
	    size_t i = no_index;
	    if (find (key2row2, it. first, i))
	    {
	      ASSERT (i != no_index);
	      const StringVector& row1 = f1. rows [it. second];
	      const StringVector& row2 = f2. rows [i];
	      ASSERT (row1. size () == row2. size ());
	      FFOR (size_t, col, row1. size ())
	        if (row1 [col] != row2 [col])
	          changes [f1. header [col]] ++;
	    }
	    else
	      deleted++;
	  }
	  	  
	  size_t added = 0;
	  for (const auto& it : key2row2)
	  {
	    size_t i = no_index;
	    if (! find (key2row2, it. first, i))
	      added++;
	  }
	  
	  cout << "added: " << added << endl;
	  cout << "deleted: " << deleted << endl;
	  for (const auto& it : changes)
	    cout << it. first << ": " << it. second << endl;
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


