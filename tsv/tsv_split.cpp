// tsv_split.cpp

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
*   Split a tsv-table into tables with unique key columns
*
*/

#undef NDEBUG

#include "../common.hpp"
#include "tsv.hpp"
using namespace Common_sp;
#include "../version.inc"

#include "../common.inc"



namespace
{
  
  
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Create tables which are a split of a source table by key columns")
  	{
      version = VERSION;
  	  addPositional ("table", "tsv-table");
  	  addPositional ("by", "Column name to split by"); 
  	  addPositional ("dir", "Output directory");
  	}
  	
  	
 
	void body () const final
	{
		const string fName  = getArg ("table");
		const string byS    = getArg ("by");
		const string outDir = getArg ("dir");

    if (byS. empty ())
      throw runtime_error ("-by should not be empty");
   

    TextTable tt (fName);
    tt. qc ();
    if (verbose ())
      tt. printHeader (cout);      
    
    if (tt. rows. empty ())
      throw runtime_error ("No data");
    
    const TextTable::ColNum byIndex = tt. col2num (byS);
    
    Vector<TextTable::Header> header (tt. header);
    header. eraseAt (byIndex);

    {    
      const StringVector by ({byS});
      tt. sort (by);
    }
    
    {
      Progress prog;
      string key;
      unique_ptr<TextTable> out;
      FFOR (TextTable::RowNum, i, tt. rows. size ())
      {
        const StringVector& row (tt. rows [i]);
        string thisKey (row [byIndex]);
        trim (thisKey);
        if (thisKey. empty ())
          throw runtime_error ("No key in row " + to_string (i + 1));
        if (thisKey != key)
        {
          if (out)
          {
            out->qc ();
            ASSERT (! key. empty ());
            OFStream f (outDir, key, "tsv");
            out->saveText (f);
          }
          out. reset (new TextTable (tt. pound, header));
          key = thisKey;
          prog (key);
        }
        ASSERT (out);
        StringVector newRow (row);
        newRow. eraseAt (byIndex);
        out->rows << std::move (newRow);
      }
      {
        ASSERT (out);
        out->qc ();
        ASSERT (! key. empty ());
        OFStream f (outDir, key, "tsv");
        out->saveText (f);
      }
    }  // For Progress
    tt. qc ();
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



