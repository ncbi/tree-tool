// index_find.cpp

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
*   Find approximately closest objects ordered by proximity using an index
*
*/

#undef NDEBUG

#include "common.hpp"
using namespace Common_sp;
#include "version.inc"

#include "common.inc"



namespace
{
  
  
struct ObjNum
{
  string obj;
  size_t num {0};
  
  bool operator< (const ObjNum &other) const
    { LESS_PART (other, *this, num);
      LESS_PART (*this, other, obj);
      return false; 
    }
};

  
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Find approximately closest objects ordered by proximity using an index.\n\
An object is a file with a list of attributes.\n\
Average time: O(A S (A (log A + log S) + log N)), where N is the number of objects, and A is the average number of attributes in an object, S = attr_size_max.\n\
")
  	{
      version = VERSION;
      addPositional ("dir", "Directory with objects");
      addPositional ("index", "Directory with an attribute index for <dir>");
      addPositional ("attr_size_max", "Attribute file size limit; objects having only attrbutes whose index files exceed <attr_size_max> are indiscernible");
      addPositional ("target", "Target object");
      addFlag ("large", "Directory <dir> is large: it is subdivided into subdirectories \"0\" .. \"" + to_string (small_hash_class_max - 1) + "\" which are the hashes of file names");
  	}
  	
  	
 
	void body () const final
	{
		const string dir               = getArg ("dir");
		const string index             = getArg ("index");
		const streamsize attr_size_max = (streamsize) arg2uint ("attr_size_max");
		const string target            = getArg ("target");
		const bool   large             = getFlag ("large");


    // size(index) = O(N A)

    StringVector attrObjs;  // Repeated objects
    {      
      StringVector bigAttrs;
      // Time = O(A (log N + log A + S))
      {
        LineInput f (target);
        while (f. nextLine ())
        {
          string& attr = f. line;
          trim (attr);
          const string attrFName (index + "/" + attr);
          if (! fileExists (attrFName))
            continue;
          if (getFileSize (attrFName) >= attr_size_max)
            bigAttrs << attr;
          else
            attrObjs << StringVector (attrFName, (size_t) 1000, true);  // PAR
        }
      }
      // size(attrObjs) = O(A S)
      // size(bigAttrs) = O(A)
      if (verbose ())
      {
        PRINT (attrObjs. size ());
        PRINT (bigAttrs. size ());
      }
      
      // Time: O(A log A)    
      bigAttrs. sort ();
      QC_ASSERT (bigAttrs. isUniq ());
      
      // Time: O(A S (log A + log S))
      StringVector objs (attrObjs);
      objs. sort ();
      objs. uniq ();
      // size(objs) = O(A S)
      if (verbose ())
        PRINT (objs. size ());
      
      // Append: attrObjs
      // Time: O(A S (log N + A log A))
      if (! bigAttrs. empty ())
        for (const string& obj : objs)
        {
          const string objFName (dir + (large ? "/" + to_string (str2hash_class (obj, false)) : "") + "/" + obj);
          if (! fileExists (objFName))
            continue;
          StringVector attrs (objFName, (size_t) 1000, true); // PAR
          attrs. sort ();
          QC_ASSERT (attrs. isUniq ());
          const size_t n = bigAttrs. getIntersectionSize (attrs);   // O(A)
          FOR (size_t, i, n)
            attrObjs << obj;
        }
    }
    // size(attrObjs) = O(A^2 S)
    
    // Time: O(A^2 S)
    unordered_map<string,size_t> obj2num;  obj2num. rehash (attrObjs. size ());
    for (const string& obj : attrObjs)
      obj2num [obj] ++;
    // size(obj2num) = O(A^2 S)

    // Time: O(A^2 S (log A + log S))
    Vector<ObjNum> objNums;  objNums. reserve (obj2num. size ());
    for (const auto& it : obj2num)
      objNums << std::move (ObjNum {it. first, it. second});
    objNums. sort ();
    // size(objNums) = O(A^2 S)
    
    for (const ObjNum& objNum : objNums)
      cout << objNum. obj << endl;
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



