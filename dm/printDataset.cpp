// printDatset.cpp

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
*   Print a dataset
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

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Read a dataset and print it")
  	{
  	  version = VERSION;
  	  addPositional ("file", dmSuff + "-file with Number attributes");
  	  addKey ("objs", "List of objects to print");
  	  addKey ("attrs", "List of attributes to print");
  	}
	
	
	
	void body () const final
	{
		const string inFName   = getArg ("file");
		const string objFName  = getArg ("objs");
		const string attrFName = getArg ("attrs");

    
    Dataset ds (inFName);
    ds. qc ();
    
    unique_ptr<Vector<size_t>> objNums;  
    if (! objFName. empty ())
    {
      objNums. reset (new Vector<size_t> ());
      objNums->reserve (ds. objs. size ());
      const StringVector objNames (objFName, ds. objs. size (), true);
      for (const string& s : objNames)
      {
        const size_t objNum = ds. getName2objNum (s);
        if (objNum == no_index)
          cerr << "object not found: " << s << endl;
        else
          *objNums << objNum;
      }
    }

    VectorPtr<Attr> attrs;  attrs. reserve (ds. attrs. size ());
    if (attrFName. empty ())
      insertAll (attrs, ds. attrs);
    else
    {
      const StringVector attrNames (attrFName, ds. attrs. size (), true);
      for (const string& s : attrNames)
        if (const Attr* attr = ds. name2attr (s))
          attrs << attr;
        else
          cerr << "attribute not found: " << s << endl;
    }
    
    const Sample sample (ds);
    sample. qc ();
    sample. save (objNums. get (), attrs, cout);
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


