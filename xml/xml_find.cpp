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
    : Application ("Find an XML context and print values as a tsv-file")
  	{
      version = VERSION;
  	  addPositional ("target", "Target XML file");
  	  addPositional ("query", "Query XML file");
  	  addKey ("variable_tag", "Tag name indicating variables", "q");
  	}
  	
  	
 
	void body () const final
	{
		const string targetFName = getArg ("target");
		const string queryFName  = getArg ("query");
		const string variableTag = getArg ("variable_tag");
	
	
	  unique_ptr<const Xml_sp::Data> target (Xml_sp::Data::load (targetFName));
    target->qc ();

	  unique_ptr<const Xml_sp::Data> query (Xml_sp::Data::load (queryFName));
    query->qc ();
    if (verbose ())
    {
      Xml::File f (cout, true, true, "XML");  // PAR
      query->saveXml (f);
      cout << endl << endl;
    }
    
    const StringVector header (query->tagName2text (variableTag));
    {
      StringVector vec (header);
      vec. sort ();
      if (! vec. isUniq ())
        throw runtime_error ("Variable tag names are not unique");
    }
    cout << '#' << header. toString ("\t") << endl;
      
    Vector<Pair<string>> output;  output. reserve (1000);  // PAR
    target->unify (*query, variableTag, output);
    
    StringVector values (header. size ());
    for (const Pair<string>& p : output)
    {
      const size_t i = header. indexOf (p. first);
      ASSERT (i != no_index);
      values [i] = p. second;
    }
    cout << values. toString ("\t") << endl;
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



