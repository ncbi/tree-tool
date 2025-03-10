// xml_txt2bin.cpp

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
*   Convert text XML to binary XML
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
    : Application ("Convert text XML to binary XML")
  	{
      version = VERSION;
  	  addPositional ("from", "Input text XML file");
  	  addPositional ("to",   "Output binary XML file");
  	}
  	
  	
 
	void body () const final
	{
		const string fromFName = getArg ("from");
		const string toFName   = getArg ("to");
		
		QC_ASSERT (fromFName != toFName);
		QC_ASSERT (! toFName. empty ());
	
	
	  Names names (10000);   // PAR
	  
	  VectorOwn<Xml_sp::Data> targetMarkupDeclarations;
	  unique_ptr<const Xml_sp::Data> in (Xml_sp::Data::load (true /*PAR*/, names, fromFName, targetMarkupDeclarations));
    in->qc ();

    cxml. reset (new Xml::BinFile (toFName, in->getName ()));
    for (const Xml_sp::Data* child : in->children)
      child->saveXml (*cxml);
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



