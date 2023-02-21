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
*   Convert object-attribute pairs to a Data Master file 
*
*/


#undef NDEBUG
#include "../../common.inc"

#include "../../common.hpp"
using namespace Common_sp;
#include "../dataset.hpp"
using namespace DM_sp;
#include "../../version.inc"



namespace 
{


string parseLine (string line,
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
	  return line;
	}
	catch (const exception &e)
	{
		throw runtime_error ("Line: " + line_orig + "\n" + e. what ());
	}
}




struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Convert object-attribute pairs to a " + dmSuff + "-file with")
    {
      version = VERSION;
  	  addPositional ("pairs", "File with lines: <obj> <attr>");
  	  addPositional ("type", "Boolean|Real|Positive");
  	  addKey ("decimals", "Decimals", "0");
  	}



	void body () const
	{
		const string pairsFName   =                   getArg ("pairs");
		const string type         =                   getArg ("type");
		const streamsize decimals = str2<streamsize> (getArg ("decimals"));
		
		QC_IMPLY (decimals, type == "Real" || type == "Positive");

		
    Set<string> objNames;
    Set<string> attrNames;
    { 
      cerr << "Pass 1 ..." << endl;
      LineInput f (pairsFName, 1000000);  // PAR
      while (f. nextLine ())
      {
        string objName, attrName;
        parseLine (f. line, objName, attrName);
        objNames  << objName;
        attrNames << attrName;
      }
    }
    cerr << "# Objects: " << objNames. size () << endl;
    cerr << "# Attributes: " << attrNames. size () << endl;
    ASSERT (! objNames. empty ());
    ASSERT (! attrNames. empty ());
    
    Dataset ds;
    for (const string& name : objNames)
      ds. appendObj (name);
    ds. setName2objNum ();

    map<string,Attr1*> name2attr;
    for (const string& name : attrNames)
    {
      Attr1* attr = nullptr;
      if (type == "Boolean")
        attr = new CompactBoolAttr1 (name, ds);
      else if (type == "Real")
      {
        auto attr_ = new RealAttr1 (name, ds, decimals);
        attr_->setAll (0.0);
        attr = attr_;
      }
      else if (type == "Positive")
      {
        auto attr_ = new PositiveAttr1 (name, ds, decimals);
        attr_->setAll (0.0);
        attr = attr_;
      }
      else
        throw runtime_error ("Unknown type " + strQuote (type));
      ASSERT (attr);
      name2attr [name] = attr;
    //attr->setAll (false);
    }
    
    {
      cerr << "Pass 2 ..." << endl;
      LineInput f (pairsFName, 1000000);  // PAR
      while (f. nextLine ())
      {
        string objName, attrName;
        const string value (parseLine (f. line, objName, attrName));
        const size_t row = ds. getName2objNum (objName);
        Attr1* attr = name2attr [attrName];
        ASSERT (attr);
        if (type == "Boolean")
          var_cast (attr->asCompactBoolAttr1 ()) -> setCompactBool (row, true);
        else if (type == "Real")
          var_cast (attr->asRealAttr1 ()) -> str2value (row, value);
        else if (type == "Positive")
          var_cast (attr->asPositiveAttr1 ()) -> str2value (row, value);
        else
          ERROR;
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



