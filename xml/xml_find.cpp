// xml_find.cpp

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
*   Find an XML context and print values
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
    : Application ("Find an XML context matching query and print values as a .tsv-file")
  	{
      version = VERSION;
  	  addPositional ("target", "Target XML file");
  	  addPositional ("query", "Query XML file");
  	  addKey ("variable_tag", "Tag name in query indicating .tsv-columns.\n\
Text of query unifying with \"<\" variable_tag \">\" column_name \"</\" variable_tag \">\" goes to the column named column_name in the output .tsv-file", "q");
  	}
  	
  	
 
	void body () const final
	{
		const string targetFName = getArg ("target");
		const string queryFName  = getArg ("query");
		const string variableTag = getArg ("variable_tag");
		
		QC_ASSERT (! variableTag. empty ());
	
	
	  Names names (10000);   // PAR
	  
	  VectorOwn<Xml_sp::Data> targetMarkupDeclarations;
	  unique_ptr<const Xml_sp::Data> target (Xml_sp::Data::load (names, targetFName, targetMarkupDeclarations));
    target->qc ();

	  VectorOwn<Xml_sp::Data> queryMarkupDeclarations;
	  unique_ptr<const Xml_sp::Data> query (Xml_sp::Data::load (names, queryFName, queryMarkupDeclarations));
    query->qc ();
    if (verbose ())
    {
      const string fName ("xml_find.xml");
      Xml::TextFile f (fName, /*false, false,*/ "XML");  // PAR
      query->saveXml (f);
      cerr << "XML file " << strQuote (fName) << " is saved" << endl;
    }
          
    const TextTable tt (target->unify (*query, variableTag));
    tt. qc ();
    tt. saveText (cout);
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



