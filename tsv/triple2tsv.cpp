// triple2tsv.cpp

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
*   Convert a tsv-file of object-attribute-value triples into a tsv-file with attributes as columns
*
*/

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
#include "tsv.hpp"
using namespace Common_sp;
#include "../version.inc"



namespace
{
  
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Convert a tsv-file of object-attribute-value triples into a tsv-file with attributes as columns")
  	{
      version = VERSION;
  	  addPositional ("table", "tsv-table where <<obj>,<attr>> is a key");
  	  addPositional ("obj", "Object column in <table>");
  	  addPositional ("attr", "Attribute column in <table>");
  	  addPositional ("val", "Value column in <table>");  	  
  	}
  	
  	
 
	void body () const final
	{
		const string fName    = getArg ("table");
		const string objName  = getArg ("obj");
		const string attrName = getArg ("attr");
		const string valName  = getArg ("val");


    {
      Set<string> cols;
      cols << objName << attrName << valName;
      if (cols. size () < 3)
        throw runtime_error ("Oject-attribute-value columns must be different");
    }

    const TextTable tt (fName);
    tt. qc ();
    
    const TextTable::ColNum objCol  = tt. col2num (objName);
    const TextTable::ColNum attrCol = tt. col2num (attrName);
    const TextTable::ColNum valCol  = tt. col2num (valName);
      
    if (tt. header [attrCol]. numeric)
      throw runtime_error ("Attribute column is numeric");
      

    TextTable out;
    out. pound = true;

    map<string,size_t> obj2num;
    out. header << tt. header [objCol];
    Set<string> attrs;
    for (const StringVector& row : tt. rows)
    {
      if (   row [objCol].  empty ()
          || row [attrCol]. empty ()
         )
        continue;
      attrs << row [attrCol];
      if (! contains (obj2num, row [objCol]))
        obj2num [row [objCol]] = obj2num. size ();
    }
    
    map<string,size_t> attr2num;
    TextTable::Header h = tt. header [valCol];
    for (const string& attr : attrs)
    {
      h. name = attr; 
      out. header << h;
      attr2num [attr] = attr2num. size () + 1/*object*/;
    }


    for (const StringVector& row : tt. rows)
    {
      if (   row [objCol].  empty ()
          || row [attrCol]. empty ()
         )
        continue;

      ASSERT (contains (obj2num, row [objCol]));
      const size_t objNum = obj2num [row [objCol]];
      ASSERT (objNum <= out. rows. size ());
      if (objNum == out. rows. size ())
      {
        StringVector newRow (attrs. size () + 1);
        newRow [0] = row [objCol];
        out. rows << move (newRow);
      }

      const size_t attrNum = attr2num [row [attrCol]];
      ASSERT (attrNum);

      if (! out. rows [objNum] [attrNum]. empty ())
        throw runtime_error ("Duplicate object/attribute: " + row [objCol] + "/" + row [attrCol]);
      out. rows [objNum] [attrNum] = row [valCol];
    }
    

    out. saveText (cout);
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}




