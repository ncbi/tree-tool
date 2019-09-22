// selColumn.cpp

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
*   Select columns of a delimited file
*
*/


#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;



namespace 
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Select <Column> out of <File>\nSpaces and tabs are automatically included in <delimiters>")
  	{
  	  addPositional ("in", "Text file name or '-'");
  	  addPositional ("col", "Column number (>=1)");
  	  addKey ("delimiters", "Delimiters");
  	}



	void body () const final
  {
	  const string fName   = getArg ("in");
	  const uint targetCol = str2<uint> (getArg ("col"));
	  string delimiters    = getArg ("delimiters");	  
	  ASSERT (targetCol > 0);
	  
	  
	  delimiters += " \t";

    istream* is = & cin;
    unique_ptr<ifstream> ifs;
    if (fName != "-")
    {
      ifs. reset (new ifstream (fName. c_str ()));  
      is = ifs. get ();
    }
    uint row = 0;
	  for (;;)
	  {
	  	row++;
		  string s; 
		  FOR (uint, col, targetCol)
		  {
		  	s = getColumn (*is, delimiters, "\n");
		  	if (s == "\n")
		  	{
		  		cout << "No column " << targetCol + 1 << " in row " << row << endl;
		  		ERROR;
		  	}
		  }
		  skipLine (*is);  
		  if (s. empty ())
		  	break;
		  cout << s << endl;  
	  }
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


