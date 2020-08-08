// obj_attr2dm.cpp

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
*   Convert object-attribute pairs to a Data Master file with Boolean attributes
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
using namespace DM_sp;
#include "../version.inc"



namespace 
{


void parseLine (string line,
                string &objName,
                string &attrName)
{
//cout << line << endl;  
  const string line_orig (line);
  try 
  {
	  replace (line, '\t', ' ');
	  replaceStr (line, "  ", " ");
	  trim (line);
	  objName  = findSplit (line);
	  attrName = findSplit (line);
	}
	catch (const exception &e)
	{
		throw runtime_error ("Line: " + line_orig + "\n" + e. what ());
	}
}




struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Convert object-attribute pairs to a " + dmSuff + "-file with Boolean attributes")
    {
      version = VERSION;
  	  addPositional ("pairs", "File with lines: <obj> <attr>");
  	}



	void body () const
	{
		const string pairsFName = getArg ("pairs");

		
    Set<string> objNames;
    Set<string> attrNames;
    { 
      cerr << "Pass 1 ..." << endl;
      LineInput f (pairsFName, 10000000, 1000000);  // PAR
      while (f. nextLine ())
      {
        string objName, attrName;
        parseLine (f. line, objName, attrName);
        objNames  << objName;
        attrNames << attrName;
      }
    }
    cerr << "# Objects: " << objNames. size () << endl;
    cerr << "# Boolean attributes: " << attrNames. size () << endl;
    ASSERT (! objNames. empty ());
    ASSERT (! attrNames. empty ());
    
    Dataset ds;
    for (const string& name : objNames)
      ds. appendObj (name);
    ds. setName2objNum ();

    map<string,CompactBoolAttr1*> name2attr;
    for (const string& name : attrNames)
    {
      auto attr = new CompactBoolAttr1 (name, ds);
      name2attr [name] = attr;
    //attr->setAll (false);
    }
    
    {
      cerr << "Pass 2 ..." << endl;
      LineInput f (pairsFName, 10000000, 1000000);  // PAR
      while (f. nextLine ())
      {
        string objName, attrName;
        parseLine (f. line, objName, attrName);
        const size_t row = ds. getName2objNum (objName);
        CompactBoolAttr1* attr = name2attr [attrName];
        ASSERT (attr);
        attr->setCompactBool (row, true);
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



