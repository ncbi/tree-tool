// tsv_cat.cpp

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
*   Concatenate tsv-tables
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
    : Application ("Concatenate tsv-tables")
  	{
      version = VERSION;
  	  addPositional ("list", "List of tsv-tables");
  	  addPositional ("file_name", "Name of the added column which will contain the file names in <list>; may be empty");
  	  addKey ("syn", TextTable::syn_format);
  	  addFlag ("skip_bad", "Skip bad tsv-files");
  	}
  	
  	
 
	void body () const final
	{
		const string listFName   = getArg ("list");
		const string fileColName = getArg ("file_name");
		const string synFName    = getArg ("syn");
		const bool   skipBad     = getFlag ("skip_bad");
		
		
		TextTable total;
		total. pound = true;
		{
  		LineInput li (listFName, 1);  // PAR
  		while (li. nextLine ())
  		{
        unique_ptr<TextTable> tab; 
        try 
        { 
          tab. reset (new TextTable (li. line, synFName));
          tab->qc (); 
        }
        catch (const exception &e)
        {
          if (skipBad)
          {
            cerr << li. line << ": " << e. what () << endl;
            continue;
          }
          throw;
        }        
        QC_ASSERT (! tab->header. empty ());
        
        if (! fileColName. empty ())
        {
    		  if (tab->hasColumn (fileColName))
    		    throw runtime_error ("Files already have column " + strQuote (fileColName));
      		tab->header << TextTable::Header (fileColName);
    		  const string fName = getFileName (li. line);
    		  for (StringVector& row : tab->rows)
    		    row << fName;
    		}
    		tab->qc ();

        Vector<TextTable::ColNum> tab2total;  tab2total. reserve (tab->header. size ());
 		    for (const TextTable::Header& h : tab->header)
 		    {
 		      TextTable::ColNum i = total. col2num_ (h. name);
 		      if (i == no_index)
 		      {
 		        i = total. header. size ();
 		        total. header << TextTable::Header (h. name);
      		  for (StringVector& row : total. rows)
      		    row << noString;
 		      }
 		      ASSERT (i != no_index);
 		      tab2total << i;
 		    }
 		    ASSERT (tab2total. size () == tab->header. size ());
 		    total. qc ();

        for (const StringVector& from : tab->rows)
        {
          StringVector to (total. header. size ());
          FFOR (size_t, i, from. size ())
            to [tab2total [i]] = from [i];
          total. rows << std::move (to);
        }
  		}
    }
		total. qc ();


    total. saveText (cout);
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



