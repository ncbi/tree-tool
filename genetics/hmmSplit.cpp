// splitHmmLib.cpp

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
*   Split an HMM library into HMM files named as the NAME field
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../version.inc"



namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Split an HMM library into HMM files named as the NAME field")
  	{
      version = VERSION;
  	  addPositional ("in", "Input HMM library");
  	  addPositional ("out_dir", "Output directory");
  	}



	void body () const final
  {
		const string in      = getArg ("in");
		const string out_dir = getArg ("out_dir");


    LineInput li (in);
    string s;
    string name;
    Progress prog;
    while (li. nextLine ())
    {
   	  s += li. line + "\n";
   	  replace (li. line, '\t', ' ');
      if (isLeft (li. line, "NAME "))
      {
        ASSERT (name. empty ());
        name = li. line. substr (5);
        trim (name);
      }
      if (! isLeft (li. line, "//"))
        continue;
   //if (contains (s, "hmmbuild -f "))  // for TnpPred_HMM_Profiles.hmm ??
      {
        OFStream ofs (out_dir, name, "HMM");
        ofs << s;
        prog ();
      }
      s. clear ();
      name. clear ();
    }
  }
};



}  // namespace



int main(int argc, 
         const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



