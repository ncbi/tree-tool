// xml_schema2dml.cpp

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
*   Generate tab-delimited files for SQL bulk insert from XML files and their XML schema
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
    : Application ("Generate tab-delimited files for SQL bulk insert from XML files and their XML schema")
  	{
      version = VERSION;
  	  addPositional ("xml", "XML file");
  	  addPositional ("xml_num", "XML file number");
  	  addPositional ("schema", "XML schema file");
  	  addPositional ("out_dir", "Directory for output tab-delimited files");
  	}
  	
  	
 
	void body () const final
	{
		const string xmlFName    =               getArg ("xml");
		const size_t xml_num     = str2<size_t> (getArg ("xml_num"));
		const string schemaFName =               getArg ("schema");
		      string dirName     =               getArg ("out_dir");
		      
		if (! isDirName (dirName))
		  dirName += '/';
	
		
	  string schemaName;
	  unique_ptr<Xml_sp::Schema> sch (Xml_sp::Schema::readSchema (schemaFName, schemaName));
	  sch->qc ();
	  
	  unique_ptr<const Xml_sp::Data> xml;
	  {
  	  TokenInput ti (xmlFName, '\0', false, false, 1000);  // PAR 
      try
      {	  
    	  VectorOwn<Xml_sp::Data> markupDeclarations;
        xml. reset (new Xml_sp::Data (ti, markupDeclarations));	
      }
    //catch (const CharInput::Error &e)
      //{ throw e; }
      catch (const exception &e)
        { ti. error (e. what (), false); }
    }    
    xml->qc ();
    
    sch->setFlatTables (dirName, nullptr);
    sch->qc ();
    
  #if 0
    cout << schemaName;
    sch->saveText (cout);
    cout << endl;
  #endif
    
    xml->writeFiles (xml_num, sch. get (), nullptr);
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



