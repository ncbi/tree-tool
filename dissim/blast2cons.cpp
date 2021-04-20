// blast2cons.cpp

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
*   Compute "conservartion" dissimilarity
*
*/


#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm> 
#include <stdexcept>
#include <cmath>
#include <cassert>
using namespace std;



namespace {


struct Hsp
// BLAST HSP
{
  // false - query, true - subject
  string seqid [2/*bool*/];

  // Positions are 0-based
  // End of a segment = start of the next segment
  size_t start [2/*bool*/];
  size_t end [2/*bool*/];
    // start < end
  
  size_t nident;
  size_t length;
  bool strand [2/*bool*/];
  

  Hsp (istream &is,
       bool tblastx)
    { 
      if (tblastx)
      // WU-tblastx: 559197279 595611723 7.8e+02 4 25.48   58  57       18     30   33    31.58   52.63  0     0       2     6       +2     23     193  +3     183192 183344
      //             qid       sid       E-value N Sprime  S   alignlen nident npos nmism pcident pcpos  qgaps qgaplen sgaps sgaplen qframe qstart qend sframe sstart send
      {
        double eValue, Sprime, pcident, pcpos;
        size_t N, score, npos, nmism, qgaps, qgaplen, sgaps, sgaplen;
        int qframe, sframe;
        is >> seqid [false]
           >> seqid [true]
           >> eValue
           >> N
           >> Sprime
           >> score
           >> length
           >> nident
           >> npos
           >> nmism
           >> pcident
           >> pcpos
           >> qgaps
           >> qgaplen
           >> sgaps
           >> sgaplen
           >> qframe
           >> start [false]
           >> end [false]
           >> sframe
           >> start [true]
           >> end [true];
      }
      else
      is >> seqid [false] 
         >> seqid [true] 
         >> start [false] 
         >> end [false] 
         >> start [true] 
         >> end [true] 
         >> nident 
         >> length; 
      if (empty ())
        return;
      if (tblastx)
      { nident *= 3;
        length *= 3;
      }
      for (const bool b : {false, true})
        assert (! seqid [b]. empty ());
      assert (seqid [false] != seqid [true]);
      for (const bool b : {false, true})
      { strand [b] = true;
        if (start [b] > end [b])
        { if (! b && ! tblastx) 
            throw runtime_error ("query positions");
          { strand [b] = false;
            swap (start [b], end [b]);
          }
        }
      }
      for (const bool b : {false, true})
      { assert (start [b] >= 1);
        start [b] --;
      }
      assert (nident > 0);
      for (const bool b : {false, true})
      { assert (nident <= size (b));
        if (size (b) > length)
        {
          cout << (int) b << endl;
          print (cout);
          throw runtime_error ("ERROR");
        }
      }
    } 
  
    
  bool empty () const
    { return seqid [false]. empty (); }  
  void print (ostream &os) const
    { os        << seqid [false]
         << ' ' << seqid [true]
         << ' ' << start [false] + 1
         << ' ' << end [false]
         << ' ' << strand [false]
         << ' ' << start [true] + 1
         << ' ' << end [true]
         << ' ' << strand [true]
         << ' ' << nident
         << ' ' << length
         << ' ' << size (false)
         << ' ' << size (true)
         << endl; 
    }
  size_t size (bool b) const
    { return end [b] - start [b]; }
};



template <const bool b>
  bool compareStart (const Hsp &h1,
                     const Hsp &h2)
  {
    if (h1. seqid [b] < h2. seqid [b])  return true;
    if (h1. seqid [b] > h2. seqid [b])  return false;
    return h1. start [b] < h2. start [b];
  }


}  // namespace




int main (int argc,
          const char* argv [])
{
  size_t lineNum = 0;
  try 
  {
    if (argc != 4)
      throw runtime_error ("\
Console: file with BLAST output: qseqid sseqid qstart qend sstart send nident length\n\
#1: 0 - BLASTN, 1 - TBLASTX\n\
#2: q-genome length\n\
#3: s-genome length\n\
Print: q-coverage s-coverage\
");
    
    const bool tblastx = (string (argv [1]) == "1");
    size_t genomeLen [2/*bool*/];
    genomeLen [false] = (size_t) atol (argv [2]);
    genomeLen [true]  = (size_t) atol (argv [3]);
    

    vector<Hsp> hsps;
    for (;;)
    {
      const Hsp hsp (cin, tblastx);
      lineNum++;
      if (cin. eof ())
        break;
      if (hsp. empty ())
        break;
      if (min ( hsp. size (false)
              , hsp. size (true)
              ) < (tblastx ? 40 * 3 : 60)  // PAR
         )
        continue;
      hsps. push_back (hsp);
    }

        
  //const double accessory = 0.05;  // PAR
  //double dist = 0;
    size_t coverage_ave = 0;
    size_t genomeLen_total = 0;
    for (const bool b : {false, true})
    {
      sort (hsps. begin (), hsps. end (), b ? compareStart<true> : compareStart<false>);
      const Hsp* prev = nullptr;
      size_t start = 0;
      size_t end = 0;
      size_t coverage = 0;
      for (const Hsp& hsp : hsps)
      {
      //hsp. print (cout);
        if (   prev
            && prev->seqid [b] == hsp. seqid [b]
            && hsp. start [b] <= end + 400   // PAR
            // --> hsp.fullStart[b] = start[b] - 200 ??
            //     hsp.fullEnd  [b] = end  [b] + 200
           )
        {
          if (end < hsp. end [b])
              end = hsp. end [b];
        }
        else
        {
          coverage += end - start;
          start = hsp. start [b];
          end   = hsp. end   [b];
        }
        prev = & hsp;
      }    
      coverage += end - start;
      assert (coverage <= genomeLen [b]);
      coverage_ave += coverage;
      genomeLen_total += genomeLen [b];
    #if 0
      const double genomeCore = (double) genomeLen [b] * (1 - 0 /*accessory*/);
      const double coverageRel = min (1.0, (double) coverage / genomeCore);
      dist += log (coverageRel);
    #endif
    }
    coverage_ave /= 2;
    const double jaccard = (double) coverage_ave / (double) (genomeLen_total - coverage_ave);
    assert (jaccard >= 0);
    assert (jaccard <= 1);
    const double dist = - log (jaccard);  // PAR  // was: 10 *
    cout << dist << endl;
  }
  catch (const exception &e)
  {
    cout << e. what () << endl;
    cout << "Line: " << lineNum << endl;
    return 1;
  }
  

  return 0;
}
