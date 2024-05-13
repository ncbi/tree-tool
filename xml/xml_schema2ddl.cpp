// xml_schema2ddl.cpp

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
*   Generate SQL "create table" commands from an XML schema
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
  
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Generate SQL \"create table\" commands from an XML schema")
  	{
      version = VERSION;
  	  addPositional ("schema", "XML schema file");
    	addFlag ("data", "'create table' etc. DDL");
    	addFlag ("index", "'create index' etc. DDL");
  	}
  	
  	
 
	void body () const final
	{
		const string schemaFName = getArg ("schema");
		const bool   dataP       = getFlag ("data");
		const bool   indexP      = getFlag ("index");
	
	  QC_ASSERT (dataP || indexP);
		
	  string name;
	  unique_ptr<Xml_sp::Schema> sch (Xml_sp::Schema::readSchema (schemaFName, name));
	  sch->qc ();
	  if (verbose ())
	  {
  	  cout << name;
  	  sch->saveText (cout);
  	  cout << endl << endl << endl;
  	}
	  
	  sch->printTableDdl (cout, dataP, indexP);
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



