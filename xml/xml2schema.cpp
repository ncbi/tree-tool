// xml2schema.cpp

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
*   Analyze an XML file and print the derived schema
*
*/

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "xml.hpp"



namespace
{
  
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Analyze an XML file and print the derived schema")
  	{
  	  addPositional ("xml", "XML file");
  	  addFlag ("print", "Print the XML file");
  	  addFlag ("store_values", "Store all field values in schema");
  	}
  	
  	
 
	void body () const final
	{
		const string xmlFName  = getArg ("xml");
		const bool printP      = getFlag ("print");
		const bool storeValues = getFlag ("store_values");
	
	
	  unique_ptr<const Xml::Data> xml;
	  {
  	  TokenInput ti (xmlFName, '\0', 100 * 1024, 1000);  // PAR 
      try
      {	  
        xml. reset (new Xml::Data (ti));	
      }
      catch (const CharInput::Error &e)
        { throw e; }
      catch (const exception &e)
        { ti. error (e. what (), false); }
    }
    
    xml->qc ();
    if (printP)
    {
      xml->saveText (cout);
      cout << endl;
    }

    unique_ptr<Xml::Schema> sch (xml->getSchema (storeValues));
    sch->qc ();
    cout << xml->name;
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



