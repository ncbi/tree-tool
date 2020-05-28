// xml_merge_schema.cpp

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
*   Merge XML schemas
*
*/

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "xml.hpp"
#include "../version.inc"



namespace
{
  
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Merge XML schemas")
  	{
      version = VERSION;
  	  addPositional ("xml_list", "List of XML schemas of the same tag, outputs of xml2schema");
  	}
  	
  	
 
	void body () const final
	{
		const string xmlListFName = getArg ("xml_list");


    string name;
    unique_ptr<Xml::Schema> sch;
    {
      LineInput f (xmlListFName, 100 * 1024, 1);  // PAR
      while (f. nextLine ())
      {
        string name1;	
        unique_ptr<Xml::Schema> sch1 (Xml::Schema::readSchema (f. line, name1));
        sch1->qc ();
        if (name. empty ())
        {
          name = move (name1);
          sch. reset (sch1. release ());
        }
        else
        {
          QC_ASSERT (name == name1);
          sch->merge (*sch1);
        }
      }
    }
    sch->qc ();
    
    cout << name;
    sch->saveText (cout);
    cout << endl;
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



