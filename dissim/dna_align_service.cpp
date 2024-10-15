// dna_align_service.cpp

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
*   DNA alignment service
*
*/


#undef NDEBUG

#include <unistd.h>
#include "../common.hpp"
using namespace Common_sp;
#include "../genetics/seq.hpp"
using namespace Seq_sp;
#include "../version.inc"

#include "../common.inc"



namespace 
{
	
	
string queryDirName;
map <string/*name*/, string/*DNA sequence*/> name2seq;
constexpr size_t dissim_inf = numeric_limits<size_t>::max ();

map <string, time_t/*finish time, 0-live*/> fNames;
std::mutex fNamesMtx;
  


size_t seqs2dissim (const string& seq1,
                    const string& seq2,
                    size_t dissim_max)
// Return: < dissim_max or dissim_inf
{
  QC_ASSERT (! seq1. empty ());
  QC_ASSERT (! seq2. empty ());
  QC_ASSERT (seq1. size () == seq2. size ());
  
  size_t diff = 0;
  FFOR (size_t, i, seq1. size ())
    if (   ! isAmbigNucl (seq1 [i])
        && ! isAmbigNucl (seq2 [i])
        && seq1 [i] != seq2 [i]
       )
    {
      diff++;
      if (diff >= dissim_max)
        return dissim_inf;
    }
      
  return diff;
}



void findDissim (const string fName,
                 const StringVector vec)
{
  const string pathName (queryDirName + "/" +  fName);

  if (vec. size () != 3)
  {
    PRINT (pathName);
    PRINT (vec. size ());
    save (cout, vec, '\n');
    cout << endl;
    throw InputError ("dissim command requires: <input file> <output file>");
  }

  {  
    OFStream outF (vec [1]);
    LineInput inF (vec [0]);
    Istringstream iss;
    string name1, name2;
    while (inF. nextLine ())
    {
      iss. reset (inF. line);
      iss >> name1 >> name2;
      const string& seq1 = name2seq [name1];
      const string& seq2 = name2seq [name2];        
      outF         << name1 
           << '\t' << name2 
           << '\t' << seqs2dissim (seq1, seq2, dissim_inf)
           << '\n';
    }
  }

  removeFile (pathName);
  
  {
    const Lock lock (fNamesMtx);
    ASSERT (! fNames [fName]);
    EXEC_ASSERT (fNames [fName] = time (nullptr));
  }
}



struct Neighbor
{
  string name;
  size_t dissim;
  
  
  Neighbor (const string &name_arg,
            size_t dissim_arg)
    : name (name_arg)
    , dissim (dissim_arg)
    {
      ASSERT (! name. empty ());
      ASSERT (dissim < dissim_inf);
    }
};



int compNeighbors (const void* a1,
                   const void* a2)
{
  const Neighbor* n1 = static_cast <const Neighbor*> (a1);
  const Neighbor* n2 = static_cast <const Neighbor*> (a2);
  ASSERT (n1);
  ASSERT (n2);
  if (n1->dissim > n2->dissim)  return  1;
  if (n1->dissim < n2->dissim)  return -1;
  return 0;
}



void findClosest (const string fName,
                  const StringVector vec)
{  
  const string pathName (queryDirName + "/" +  fName);

  if (vec. size () != 4)
  {
    PRINT (pathName);
    PRINT (vec. size ());
    save (cout, vec, '\n');
    cout << endl;
    throw runtime_error ("closest command requires: <object> <search objects> <output file>");
  }


  {
    const string& name1 = vec [0];
    const string& seq1 = name2seq [name1];
    QC_ASSERT (! seq1. empty ());
    
    StringVector scope (vec [1], (size_t) 1000, true);
    scope. sort ();
    QC_ASSERT (scope. isUniq ());

    constexpr size_t heap_size = 100;  // PAR
    // Use index of least frequent nucleotides for the first heap_size items ??
    Heap<const Neighbor> heap (compNeighbors, nullptr, heap_size);  
    size_t dissim_max = 0;
    for (const string& name : scope)
    {
      const string& seq2 = name2seq [name];
      QC_ASSERT (! seq2. empty ());
      ASSERT (heap. size () <= heap_size);  
      const size_t dissim = seqs2dissim (seq1, seq2, heap. size () == heap_size ? dissim_max : dissim_inf);
      if (dissim == dissim_inf)
        continue;
      if (heap. size () == heap_size)
      {
        ASSERT (dissim < dissim_max);
        const Neighbor* n = heap. getMaximum ();
        ASSERT (n);
        ASSERT (n->dissim == dissim_max);
        heap. deleteMaximum ();
        delete n;
        ASSERT (! heap. empty ());
        dissim_max = heap. getMaximum () -> dissim;
      }
      maximize (dissim_max, dissim);
      heap << new Neighbor (name, dissim);
    }

    VectorOwn<Neighbor> neighbors;  neighbors. reserve (heap. size ());
    while (! heap. empty ())
    {
      neighbors << heap. getMaximum ();
      heap. deleteMaximum ();
    }
    
    OFStream outF (vec [2]);
    FOR_REV (size_t, i, neighbors. size ())
      outF         << neighbors [i] -> name 
         //<< '\t' << neighbors [i] -> dissim  ??
           << '\n';
  }
  

  removeFile (pathName);
  
  {
    const Lock lock (fNamesMtx);
    ASSERT (! fNames [fName]);
    EXEC_ASSERT (fNames [fName] = time (nullptr));
  }
}



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("DNA alignment service", true, false, true)
    {
      version = VERSION;
  	  addPositional ("dna_align", "File with DNA alignment: <# sequences> <sequence length><EOL>{<sequence name> <space>+ <DNA sequence>}+");
  	  addPositional ("query_dir", "Directory with queries:\n<file>.dissim: <input file> <nl> <output file> <nl> END\n<file>.closest: <object> <nl> <search objects> <nl> <output file> <nl> END");
  	  addFlag ("quiet", "Quiet");
    }


	
	void body () const final
  {
	  const string alignFName   = getArg ("dna_align");
	               queryDirName = getArg ("query_dir");
	  const bool   quiet        = getFlag ("quiet");
    
    
    if (threads_max <= 1)
      throw runtime_error ("More than 1 thread (including dispatcher) is needed");
    
    
    {
      IFStream is (alignFName);
  	  char* buf = nullptr;
		  if (! is. rdbuf () -> pubsetbuf (buf, 1000000))   // PAR
		  	throw logic_error ("Cannot allocate buffer to " + strQuote (alignFName));      
		  LineInput f (is);
      
      size_t n = 0;
      size_t len = 0;
      {
        EXEC_ASSERT (f. nextLine ());
        istringstream iss (f. line);
        iss >> n >> len;
      }
      QC_ASSERT (n);
      QC_ASSERT (len);
      
      while (f. nextLine ())
      {
        size_t pos = f. line. find (' ');
        QC_ASSERT (pos != string::npos);
        const string name (f. line. substr (0, pos));
        QC_ASSERT (! name. empty ());
        
        while (f. line [pos] == ' ')
          pos++;
        QC_ASSERT (f. line. size () - pos == len);
        strLower (f. line);  // ??
        FFOR_START (size_t, i, pos, f. line. size ())
          QC_ASSERT (strchr (extSparseDnaAlphabet, f. line [i]));
        name2seq [name] = f. line. substr (pos);
      }
      QC_ASSERT (f. lineNum == n + 1);
      QC_ASSERT (name2seq. size () == n);
    }
    

    {
      Progress prog (0, quiet ? 0 : 1);  // PAR
      time_t start = 0;
      bool quit = false;
      while (! quit)
      {
        DirItemGenerator dir (0, queryDirName, false);
        string fName;
        bool found = false;
        while (dir. next (fName))
        {
          ASSERT (! fName. empty ());
          if (isLeft (fName, ".nfs"))  // Techincal temporary Linux file 
            continue;
          const string cmd (getFileExtension (fName));  // doc ??
          if (cmd == "quit")
          {
            quit = true;
            break;
          }
          const string pathName (queryDirName + "/" +  fName);
          StringVector vec;
          try { vec = std::move (StringVector (pathName, (size_t) 10, true)); }
            catch (const exception &)
              { continue; }
          if (   vec. empty () 
              || vec. back () != "END"
             )
            continue;            
          found = true;
          {
            const Lock lock (fNamesMtx);
            if (contains (fNames, fName))
              continue;
          }          
          prog (fName);
          for (;;)
          {
            const Lock lock (fNamesMtx);
            const time_t t = time (nullptr);
            ASSERT (t);
            size_t nLive = 0;
            for (Iter<map<string,time_t>> it (fNames); it. next (); )
              if (it->second)
              {
                if (difftime (t, it->second) > 1/*sec.*/)  // Time sufficient to flush file removal
                  it. erase ();
              }
              else
                nLive++;
            if (nLive < threads_max - 1)
            {
              fNames [fName] = 0;
              break;
            }
          }
          // First create output file, then remove fName!
          thread th;
          if (cmd == "dissim")
            th = std::move (thread (findDissim, fName, vec));  
          else if (cmd == "closest")
            th = std::move (thread (findClosest, fName, vec));
          else
            throw InputError ("Unknown command: " + strQuote (fName));
          th. detach ();
        }
        if (found)
          start = 0;
        else if (start)
        {
          const time_t stop = time (nullptr);
          sleepNano (1000);
          if (difftime (stop, start) > 1/*sec.*/)  // PAR
            sleep (5);  // PAR
          start = stop;
        }
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



