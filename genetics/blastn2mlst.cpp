// blastn2mlst.cpp

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
*   7-locus MLST by blastn
*
*/
   
   
#undef NDEBUG 

#include "../common.hpp"
#include "../tsv/tsv.hpp"
using namespace Common_sp;
#include "../version.inc"

#include "../common.inc"



namespace
{
  

struct Gene
{
  string allele;
    // Number
  TextTable::ColNum colNum {no_index};
};
  
  
  
// ThisApplication

struct ThisApplication : Application
{
  static const constexpr char* allFName {"all.fa"};
  static const constexpr char* profilesFName {"profiles.list"};
  static constexpr size_t unknown {0};
  static constexpr size_t novel {999999};
  
  
  ThisApplication ()
    : Application ("Print 7-locus MLST number given the BLASTN result against an MLST scheme. " + to_string (unknown) + " - unknown, " + to_string (novel) + " - novel")
    {
      version = VERSION;
      addPositional ("alleles", "List of MLST alleles (first tab-separated field in each line)");
      addPositional ("mlst_dir", string ("MLST directory with ") + allFName + " (all <locus>.fasta-files) and " + profilesFName);
    }



  void body () const final
  {
    const string blastnFName = getArg ("alleles");
    const string mlstDirName = getArg ("mlst_dir");
    

    map<string,Gene> name2gene;
    {
      LineInput f (mlstDirName + "/" + allFName);
      while (f. nextLine ())
      {
        if (f. line. empty ())
          continue;
        if (f. line [0] != '>')
          continue;
        const size_t pos = f. line. rfind ('_');
        QC_ASSERT (pos != string::npos);
        name2gene [f. line. substr (1, pos - 1)] = Gene ();
      }
    }
    if (verbose ())
    {
      cerr << "Genes:" << endl;
      for (const auto& it : name2gene)
        cerr << it. first << endl;
    }

    
    {
      LineInput f (blastnFName);
      while (f. nextLine ())
      {
        const size_t tab_pos = f. line. find ('\t');
        if (tab_pos != string::npos)
          f. line. erase (tab_pos);
          
        const size_t pos = f. line. rfind ('_');
        QC_ASSERT (pos != string::npos);
          
        const string gene = f. line. substr (0, pos);
        auto it = name2gene. find (gene);
        QC_ASSERT (it != name2gene. end ());
        if (! it->second. allele. empty ())
        {
          cout << unknown << endl;
          return;
        }
        
        string numS (f. line. substr (pos + 1));
        const size_t num = str2<size_t> (numS);
        QC_ASSERT (num > 0);
        QC_ASSERT (num != no_index);
        it->second. allele = std::move (numS);
      }
    }
    for (const auto& it : name2gene)
      if (it. second. allele. empty ())
      {
        cout << unknown << endl;
        return;
      }
    
    
    {
      Progress::disable ();
      TextTable tt (mlstDirName + "/" + profilesFName);
      tt. qc ();
      for (auto& it : name2gene)
      {
        it. second. colNum = tt. col2num (it. first);
        QC_ASSERT (it. second. colNum > 0);
      }
      for (const StringVector& row : tt. rows)
      {
        QC_ASSERT (row [0] != "0");
        bool good = true;
        for (const auto& it : name2gene)
          if (row [it. second. colNum] != it. second. allele)
          {
            good = false;
            break;
          }
        if (good)
        {
          cout << row [0] << endl;
          return;
        }      
      }
    }
    cout << novel << endl;
  }
};



}   // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



