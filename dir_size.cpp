// dir_size.cpp

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
*   Print .tsv file with columns: directory file_size
*
*/

#undef NDEBUG

#include "../cpp/common.hpp"
#include "../cpp/tsv/tsv.hpp"
using namespace Common_sp;

#include <ftw.h>


#include "../cpp/common.inc"



namespace
{
  
  
//string rootDir;
//const StringVector pipelines {{"Assembly", "Annotation", "SnpClustering", "FTP", "Remapping", "PreSubmission", "KmerAnalysis", "ReadMapping", "Amr", "Browser"}};
unordered_map<string/*directory with files*/,double/*M*/> dir2size;
long long files = 0;  // Not used


int nftw_process (const char* filepath, 
                  const struct stat* info,
                  const int typeflag, 
                  FTW* /*pathinfo*/) 
// Update dir2size, files
{
  ASSERT (info);
//prog (string (filepath));
  if (typeflag == FTW_F)
  {
    const string s (filepath);
    string key (getDirName (s));
    trimSuffix (key, "/");
    dir2size [key] += double (info->st_size) / (1024 * 1024);  
    files++;
  #if 0
    long taskrun_id = 0;
    long request_id = 0;
    ASSERT (isLeft (s, rootDir));
    StringVector v (s, '/', false);
    size_t pipelinePos = 0;
    while (pipelinePos < v. size ())
      if (pipelines. contains (v [pipelinePos]))
        break;
      else
        pipelinePos++;
    if (pipelinePos < v. size ())
    {
      size_t taskPos = pipelinePos + 2;
      bool taxgroup = false;
      if (v [pipelinePos] == "Browser")
      {
        taskPos += 1;
        taxgroup = true;
      }
      else if (v [pipelinePos] == "FTP")
      {
        taskPos += 2;
        taxgroup = true;
      }
      PRINT (v. toString ("/"));  // ??
      if (taskPos < v. size ())
      {
        const StringVector keyDir (v [taskPos], '-', false);
        if (   keyDir. size () == 2
            && str2<long> (keyDir [0], taskrun_id)
            && str2<long> (keyDir [1], request_id)
           )
        {
          ASSERT (taskrun_id > 0);
          ASSERT (request_id > 0);
          v. eraseMany (taskPos + 1);          
          total += info->st_size;  
        }
      }
    }
  #endif
  }
  return 0;
}



struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Print a .tsv-file with columns: " + columns (). toString (" ") + ".\n" + strQuote (columns () [0]) + " is a key.\nTime: as in du")
  	{
  	  addPositional ("dir", "Directory to scan");  // "containing GPipe PathogenDetect output files" ??
  	}
  static StringVector columns ()
    { StringVector v;
      v << "dir" << "size_M";
      return v;
    }



	void body () const final
	{
		const string rootDir = getArg ("dir");
	//addDirSlash (rootDir);
	
	  if (getFiletype (rootDir, false) == Filetype::dir)
    {	
  	  dir2size. rehash (10000000);  // PAR  		
      if (nftw (rootDir. c_str (), nftw_process, 15, FTW_PHYS))  // PAR
        throw runtime_error ("nftw");
    }
  //cout << files << endl;
  //cout << dir2size. size () << endl;
  
    TsvOut tsv (& cout);
    for (const string& s : columns ())
      tsv << s;
    tsv. newLn ();
    for (const auto& it : dir2size)
    {
      tsv << it. first << it. second;
      tsv. newLn ();
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



