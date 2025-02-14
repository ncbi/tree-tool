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

#include "../common.hpp"
using namespace Common_sp;
#include "xml.hpp"
#include "../version.inc"

#include "../common.inc"



namespace
{
  
  
struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Analyze an XML file and print the derived schema")
  	{
      version = VERSION;
  	  addPositional ("xml", "XML file");
  	  addFlag ("no_xml_header", "XML file has no header");
  	  addKey ("print", "Output XML file");
  	  addFlag ("store_values", "Store all field values in schema");
  	}
  	
  	
 
	void body () const final
	{
		const string xmlFName   = getArg ("xml");
		const bool headerP      = ! getFlag ("no_xml_header");
		const string printFName = getArg ("print");
		const bool storeValues  = getFlag ("store_values");
	
	
	  Names names (10000);  // PAR
	  VectorOwn<Xml_sp::Data> markupDeclarations;
	  unique_ptr<const Xml_sp::Data> xml (Xml_sp::Data::load (headerP, names, xmlFName, markupDeclarations));
    xml->qc ();
    
    if (! printFName. empty ())
    {
      Xml::TextFile f (printFName,"XML");  // PAR
      xml->saveXml (f);
      cout << endl << endl;
    }

    unique_ptr<Xml_sp::Schema> sch (xml->createSchema (storeValues));
    sch->qc ();
    cout << xml->getName ();
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



