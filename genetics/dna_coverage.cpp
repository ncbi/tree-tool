// dna_coverage.cpp

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
*   Coverage length of DNA
*
*/
   
   
#undef NDEBUG 
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../version.inc"



#define FORMAT  sseqid >> length >> nident >> qstart >> qend >> qlen >> sstart >> send >> slen >> stitle




struct Subject
{
  string sseqid;
  size_t qstart {numeric_limits<size_t>::max()};
    // Minimum
  size_t qend {0};
    // Maximum
  size_t sstart {numeric_limits<size_t>::max()};
    // Minimum
  size_t send {0};
    // Maximum
  size_t slen {0};
  size_t length {0};
    // Alignment length, cumulative
  size_t nident {0};
    // Cumulative
  string stitle;
  string plasmidName;
  size_t coverage {0};
    // Query coverage
  
  Subject () = default;
  void saveText (ostream &os) const
    { string stitle_ (stitle);
      trimPrefix (stitle_, sseqid + "_");
      replace (stitle_, '_', '#');
      os         << sseqid                  // 1
         << '\t' << qstart                  // 2
         << '\t' << qend                    // 3
         << '\t' << sstart                  // 4
         << '\t' << send                    // 5
         << '\t' << slen                    // 6
         << '\t' << nvl (plasmidName, "NA") // 7
         << '\t' << coverage                // 8
         << '\t' << double (coverage) / double (slen)  // 9
         << '\t' << double (nident) / double (length)  // 10
         << '\t' << length                  // 11
         << '\t' << stitle_                 // 12
         ;
    }
    
  bool operator< (const Subject& other) const
    { LESS_PART (other, *this, coverage);
      LESS_PART (*this, other, sseqid);
      return false;
    }
};



void add (Vector<Subject> &subjects, 
          Subject &subj,
          const Vector<uint> &chars)
{
  if (subj. sseqid. empty ())
    return;
    
  subj. coverage = 0;
  for (const uint c : chars)
    subj. coverage += (bool) c;    

  subjects << move (subj);
  subj = Subject ();
}




// ThisApplication

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Find coverage of DNA")
    {
      version = VERSION;
      addPositional ("in", "BLASTN output in format: " XSTR(FORMAT) ", sorted by sseqid");
      addFlag ("plasmid", "Report plasmid name");
    }



  void body () const final
  {
    const string inFName = getArg ("in");
    const bool   plasmid = getFlag ("plasmid");
    
    
    Vector<Subject> subjects;  subjects. reserve (10000);  // PAR
    Subject subj;  // previous
    size_t qlen_prev = 0;  // >0 => constant
    {
      LineInput f (inFName);
      string qseqid, sseqid;
      size_t length, nident, qstart, qend, qlen, sstart, send, slen;
      string stitle;
      string plasmidName_prev;
      Vector<uint> chars;    
      while (f. nextLine ())
      {
        replace (f. line, ' ', '_');
        istringstream iss (f. line);
        slen = 0;
        iss >> FORMAT;
        const size_t sstart_ = min (sstart, send);
        const size_t send_   = max (sstart, send);
        QC_ASSERT (nident);
        QC_ASSERT (nident <= length);
        QC_ASSERT (qstart);
        QC_ASSERT (qstart < qend);
        QC_ASSERT (qend <= qlen);
        QC_ASSERT (sstart_);
        QC_ASSERT (send_ <= slen);
        QC_ASSERT (nident <= min (qlen, slen));
        QC_ASSERT (qlen);
        QC_ASSERT (slen);
        QC_ASSERT (subj. sseqid <= sseqid);
        
        string plasmidName;
        if (plasmid)
        {
          string s = stitle;
          do 
          { 
            trimSuffix (s, ",");
            if (   s. size () >= 2 
                && s [0] == 'p'
                && isUpper (s [1])
               )
              plasmidName = s;
            s. clear ();
            iss >> s; 
          }
          while (plasmidName. empty () && ! s. empty ());
        }
        
        if (qlen_prev)
          { QC_ASSERT (qlen_prev == qlen); }
        else  // First time
        {
          qlen_prev = qlen;
          chars. resize (qlen, 0);
        }
        
        if (subj. sseqid != sseqid)
        {
          add (subjects, subj, chars);
          chars. setAll (0);
        }
        FOR_START (size_t, i, qstart - 1, qend)  // 1-based
          chars [i] ++;
          
        // Subject attributes
        subj. sseqid = move (sseqid);
        minimize (subj. qstart, qstart);
        maximize (subj. qend, qend);
        minimize (subj. sstart, sstart);
        maximize (subj. send, send);
        subj. slen = slen;
        subj. length += length;
        subj. nident += nident;
        subj. stitle = stitle;
        subj. plasmidName = plasmidName;
      }
      add (subjects, subj, chars);
    }

      
    subjects. sort ();
    if (! subjects. empty ())
    {
      // Best Subject
      subjects [0]. saveText (cout);      
      cout << '\t' << qlen_prev << endl;
    }
  }
};




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



