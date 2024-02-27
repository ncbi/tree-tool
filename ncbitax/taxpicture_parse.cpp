// mutate.cpp

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
*   Print .tsv-file of triples for a \"taxpropdump picture\" output
*
*/
   
   
#undef NDEBUG

#include "../common.hpp"
using namespace Common_sp;

#include "../common.inc"



namespace 
{



Pair<string> get (string &s)
// Return: <attr,value>
// Update: s
{
  QC_ASSERT (s. substr (0, 2) == "\\\"");
  s. erase (0, 2);
  
  const size_t start = s. find ("\\\":\\\"");
  QC_ASSERT (start != string::npos);
  const string attr (s. substr (0, start));
  QC_ASSERT (! attr. empty ());
  QC_ASSERT (! contains (attr, ' '));
  s. erase (0, start + 5);
  
  const size_t end = s. find ("\\\"");
  QC_ASSERT (end != string::npos);
  const string value (s. substr (0, end));
  s. erase (0, end + 2);
    
  return Pair<string> (attr, value);
}



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print .tsv-file of triples for a \"taxpropdump picture\" output")
    {
  	  addPositional ("in", "Output of \"taxpropdump picture\"");
    }


	
	void body () const final
  {
	  const string inFName  = getArg ("in");
  
  
    cout << "#tax\tattr\tval" << endl;
    LineInput f (inFName);
    while (f. nextLine ())
    {
      const uint tax = str2<uint> (findSplit (f. line, '\t'));
      QC_ASSERT (tax);
      try
      {
        const string bar = findSplit (f. line, '\t');
        QC_ASSERT (bar == "|");
        string& s = f. line;
        QC_ASSERT (s. substr (0, 2) == "\"{");
        s. erase (0, 2);
        const size_t len = s. size ();
        QC_ASSERT (len >= 2);
        QC_ASSERT (s. substr (len - 2) == "}\"");
        s. erase (len - 2);
        for (;;)
        {
          const Pair<string> attrVal (get (s));
          cout << tax << '\t' << attrVal. first << '\t' << attrVal. second << endl;
          if (s. empty ())
            break;
          QC_ASSERT (s [0] == ',');
          s. erase (0, 1);
        }
      }
      catch (const exception &e)
      {
        cerr << tax << '\t' << e. what () << endl;
      }
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



