// hmmAddCutoff.cpp

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
*   Add cutoffs to an HMM library
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;




namespace 
{


struct ThisApplication : Application
{
	ThisApplication ()
	  : Application ("Add cutoffs to an HMM library")
		{
		  addPositional ("in", "Input HMM library");
		  addPositional ("cutoffFile", "File with lines: <hmmName> <cuttoff> ...");
		  addPositional ("type", "Cutoff type: TC|GA|NC");
		  addPositional ("out", "Output HMM library");
		}



  void body () const final
	{
		const string in         = getArg ("in");
		const string cutoffFile = getArg ("cutoffFile");
		const string type       = getArg ("type");
		const string out        = getArg ("out");
				
		ASSERT (   type == "TC"
		        || type == "GA"
		        || type == "NC"
		       );


    map<string/*hmmName*/,double> hmm2cutoff;
    {
      LineInput li (cutoffFile);
      while (li. nextLine ())
      {
      	istringstream iss (li. line);
      	string name;
      	double cutoff = 0;
      	iss >> name >> cutoff;
      //ASSERT (cutoff > 0);
      	cutoff = round (cutoff * 100) / 100;
        hmm2cutoff [name] = cutoff;
      }
    }


    OFStream ofs (out);
 	  double cutoff = 0;
    LineInput li (in);
    Progress prog;
    while (li. nextLine ())
    {
      string s (li. line);
   	  replace (s, '\t', ' ');
      if (isLeft (s, "NAME "))
      {
        ASSERT (! cutoff);
        findSplit (s);
        trim (s);
        cutoff = hmm2cutoff [s];
      //ASSERT (cutoff > 0);
        prog (s);
      }
      else if (isLeft (s, "STATS LOCAL MSV "))
      {
      //ASSERT (cutoff > 0);
        ofs << type << "\t" << cutoff << ' ' << 0 << ';' << endl;
        cutoff = 0;
      }
      ofs << li. line << endl;
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



