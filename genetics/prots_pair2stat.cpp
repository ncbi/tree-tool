// prots_pair2stat.cpp

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
*   Compute amino acid substitution statistic
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "seq.hpp"
using namespace Seq_sp;
#include "../dm/numeric.hpp"
using namespace DM_sp;
#include "../version.inc"



namespace 
{	


struct Obj 
{
  Vector<Peptide> peptides;  
  

  explicit Obj (const string &fName)
    { 
      peptides. reserve (512);  // PAR
      {
        Multifasta faIn (fName, true);
        while (faIn. next ())
        {
          Peptide pep (faIn, 1000/*PAR*/, true);
          pep. ambig2X ();
          pep. name = pep. getId ();
          pep. qc ();
          peptides << std::move (pep);
        }   
      }   
      peptides. sort (comp);
    }
  Obj () = default;
private:
  static bool comp (const Peptide &p1,
                    const Peptide &p2)
    { return p1. name < p2. name; }
};



static constexpr size_t alphabetSize = 21;



size_t aa2index (char aa)
// Return: < alphabetSize or no_index
{
  static const string alphabet ("ACDEFGHIKLMNPQRSTVWY-");
  ASSERT (alphabet. size () == alphabetSize);
  const size_t pos = alphabet. find (aa);
  if (pos == string::npos)
    return no_index;
  return pos;
}


	
}  // namespace




struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Compute amino acid substitution statistic")
  {
    version = VERSION;
	  addPositional ("pairs_dist", "File with lines: <Alignment file1> <Alignment file2> <evolution distance>");
  }


	
	void body () const final
  {
	  const string pairs_dist = getArg ("pairs_dist");


    MeanVar meanVars [alphabetSize] [alphabetSize];
    {
      map<string/*fName*/,Obj> name2obj;    
      LineInput in (pairs_dist, 1000);  // PAR
      Istringstream iss;
      string fName1, fName2;    
      while (in. nextLine ())
      {
        iss. reset (in. line);
        double dist = NaN; 
        iss >> fName1 >> fName2 >> dist;
        ASSERT (dist >=  0);
        if (! contains (name2obj, fName1))  name2obj [fName1] = std::move (Obj (fName1));
        if (! contains (name2obj, fName2))  name2obj [fName2] = std::move (Obj (fName2));
        const Obj& obj1 = name2obj [fName1];
        const Obj& obj2 = name2obj [fName2];
        size_t i = 0;
        for (const Peptide& p1 : obj1. peptides)
        {
          while (   i < obj2. peptides. size () 
                 && obj2. peptides [i]. name < p1. name
                )
            i++;
          if (i == obj2. peptides. size ())
            break;
          const Peptide& p2 = obj2. peptides [i];
        //cout << p1. name << ' ' << p2. name << endl;  // ??
          if (p2. name == p1. name)
          {
            ASSERT (p1. seq. size () == p2. seq. size ());
            FFOR (size_t, j, p1. seq. size ())
            {
              size_t index1 = aa2index (p1. seq [j]);
              size_t index2 = aa2index (p2. seq [j]);
              if (index1 == no_index)
                continue;
              if (index2 == no_index)
                continue;
              if (index1 > index2)
                swap (index1, index2);
              ASSERT (index1 <= index2);
              ASSERT (index2 < alphabetSize);
              meanVars [index1] [index2] << dist;
            //cout << index1 << ' ' << index2 << ' ' << dist << endl;  // ??
            }
          }
        }
      }
    }
    
    
    FOR (size_t, index2, alphabetSize)
    {
      FOR (size_t, index1, index2 + 1)
      {
        const MeanVar& mv = meanVars [index1] [index2];
        cout << mv. getMean () << ' ' << mv. getVar () << "  ";
      }
      cout << endl;
    }
  }
};



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



