// pairs2attr2.cpp

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
*   Convert pairs of objects to a two-way attribute
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
using namespace DM_sp;



namespace 
{


double parseLine (string line,
                  string &obj1,
                  string &obj2,
                  size_t attr_num)
// Return: attribute value
{
//cout << line << endl;  
  const string line_orig (line);
  try 
  {
	  replace (line, '\t', ' ');
	  replaceStr (line, "  ", " ");
	  trim (line);
	  obj1 = findSplit (line);
	  obj2 = findSplit (line);
	  FOR (size_t, i, attr_num + 1)
	  {
	    const string s (findSplit (line));
	    if (i == attr_num)
	    {
	    //cout << s << endl;
	      return str2real (s);
	    }
	  }
	  throw runtime_error ("Attribute not found for " + obj1 + ", " + obj2);
	}
	catch (const exception &e)
	{
		throw runtime_error ("Line: " + line_orig + "\n" + e. what ());
	}
}




struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Convert two-way attributes to a " + dmSuff + "-file")
    {
  	  addPositional ("pairs", "File with lines: <obj1> <obj2> <attr1> <attr2> ...");
  	  addPositional ("attr_num", "Attribute number (1-based) in each line");
  	  addPositional ("attr_name", "Attribute name");
  	  addPositional ("decimals", "Number of decimals");
  	  addFlag ("distance", "Attribute is a distance");
  	}



	void body () const final
	{
		const string pairsFName = getArg ("pairs");
		      size_t attr_num   = str2<size_t> (getArg ("attr_num"));
		const string attr_name  = getArg ("attr_name");
		const uint   decimals   = str2<uint> (getArg ("decimals"));
		const bool   isDistance = getFlag ("distance");

		
		ASSERT (attr_num);
		attr_num--;
		
		
    Set<string> objNames;
    {
      LineInput f (pairsFName);
      while (f. nextLine ())
      {
        string obj1, obj2;
        parseLine (f. line, obj1, obj2, attr_num);
        objNames << obj1;
        objNames << obj2;
      }
    }

    Dataset ds;
    for (const string& name : objNames)
      ds. appendObj (name);
    ds. setName2objNum ();
    
    RealAttr2* attr = isDistance 
                        ? new PositiveAttr2 (attr_name, ds, decimals)
                        : new RealAttr2     (attr_name, ds, decimals);
    {
      LineInput f (pairsFName);
      while (f. nextLine ())
      {
        string obj1, obj2;
        const Real value = parseLine (f. line, obj1, obj2, attr_num);
        const size_t row = ds. getName2objNum (obj1);
        const size_t col = ds. getName2objNum (obj2);
      #if 0
        if (! attr->isMissing2 (row, col))
        {
          cout << obj1 << ' ' << obj2 << ": " << "duplicate value" << endl;
          ERROR;
        }
      #endif
        if (attr->isMissing2 (col, row))
          attr->putSymm (row, col, value); 
        else
          attr->put (row, col, value); 
      }
    }
    
    if (isDistance)
      FFOR (size_t, row, ds. objs. size ())
        attr->put (row, row, 0.0);
    
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



