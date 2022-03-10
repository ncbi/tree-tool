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
*   Coverage of DNA
*
*/
   
   
#undef NDEBUG 
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../version.inc"



#define FORMAT  qseqid >> sseqid >> length >> nident >> qstart >> qend >> qlen >> sstart >> send >> slen >> stitle



struct Hsp
{
  size_t qstart {0};
  size_t qend {0};
  size_t sstart {0};
  size_t send {0};
  bool sstrand {false};
  // Alignment
  size_t length {0};
  size_t nident {0};
  
  
  Hsp (size_t qstart_arg,
       size_t qend_arg,
       size_t sstart_arg,
       size_t send_arg,
       size_t length_arg,
       size_t nident_arg)
    // 1-based, end = last position
    : qstart (qstart_arg)
    , qend (qend_arg)
    , sstart (min (sstart_arg, send_arg))
    , send   (max (sstart_arg, send_arg))
    , sstrand (sstart_arg <= send_arg)
    , length (length_arg)
    , nident (nident_arg)
    {
      QC_ASSERT (qstart);
      qstart--;
      QC_ASSERT (qstart < qend);
      
      QC_ASSERT (sstart);
      sstart--;
      QC_ASSERT (sstart < send);
      
      QC_ASSERT (nident);
      QC_ASSERT (nident <= length);
      QC_ASSERT (nident <= qLen ());
      QC_ASSERT (nident <= sLen ());
      QC_ASSERT (length < qLen () + sLen ());
    }
  Hsp () = default;
    
    
  size_t qLen () const
    { return qend - qstart; }
  size_t sLen () const
    { return send - sstart; }
  void saveText (ostream &os) const
    { os         << qstart + 1
         << '\t' << qend
         << '\t' << qLen ()
         << '\t' << sstart + 1
         << '\t' << send
         << '\t' << sLen ()
         << '\t' << sstrand
         << '\t' << length
         << '\t' << nident
         << '\t' << (double) nident / (double) length * 100.0;
    }
  static void saveHeader (ostream &os)
    { os         << "qstart"
         << '\t' << "qend"
         << '\t' << "qcoverage"
         << '\t' << "sstart"
         << '\t' << "send"
         << '\t' << "scoverage"
         << '\t' << "sstrand"
         << '\t' << "length"
         << '\t' << "nident"
         << '\t' << "pident";
    }
};
      


struct Subject
{
  string stitle;
  Vector<Hsp> hsps;
    // !empty()
  // Function of hsps
  Vector<uint> schars;    
    // slen = schars.size()
  size_t length_sum {0};
  size_t nident_sum {0}; 
  size_t scoverage {0};
  
  
  Subject (string &&stitle_arg,
           size_t slen)
    : stitle (move (stitle_arg))
    , schars (slen)
    { QC_ASSERT (! schars. empty ()); }
  Subject () = default;
  Subject (Subject &&) = default;
  Subject& operator= (Subject &&) = default;
  void addHsp (const Hsp& hsp)
    {
      hsps << hsp;
      QC_ASSERT (hsp. send <= schars. size ());
      FOR_START (size_t, i, hsp. sstart, hsp. send) 
        schars [i] ++;
      length_sum += hsp. length;
      nident_sum += hsp. nident;
    }
  void finish ()
    {
      ASSERT (! hsps. empty ());
      ASSERT (scoverage == 0);
      for (const uint c : schars)
        scoverage += (bool) c;   
      
    }
};



static const string na ("NA");
string queryName;
string subjectName;
enum class Mode {none, combine, /*missed,*/ all};
Mode mode {Mode::none};
Vector<uint> qchars;    
map<string/*sseqid*/,Subject> sseqid2subject; 



void reportSubjects (const string &qseqid,
                     bool force)
{
  ASSERT (mode != Mode::none);
  ASSERT (qseqid. empty () == qchars. empty ());
  IMPLY (qseqid. empty (), sseqid2subject. empty ());

//if (mode == Mode::missed)
  //return;
    
  const size_t qlen = qchars. size ();

  size_t qcoverage = 0;
  for (const uint c : qchars)
    qcoverage += (bool) c;   
  ASSERT (qcoverage <= qlen);

  if (sseqid2subject. empty ())
  {
    if (! force)
      return;
  }
  else
    ASSERT (qcoverage);
    
  // Subject::scoverage
  for (auto& it : sseqid2subject)
    it. second. finish ();

  // blast has no qtitle
  string qtitle (qseqid);
  const string qseqid_ (findSplit (qtitle, '|'));
  IMPLY (qseqid_. empty (), force && mode == Mode::combine);
  replace (qtitle, '_', ' ');
  trim (qtitle);
  

  // cout
  switch (mode)
  {
    case Mode::combine:
      {
        StringVector sseqids;
        size_t slen_sum = 0;
        size_t scoverage_sum = 0;
        size_t length_sum = 0;
        size_t nident_sum = 0;        
        for (const auto& it : sseqid2subject)
        {
          QC_ASSERT (! contains (it. first, ' '));
          sseqids << it. first;
          const Subject& subj = it. second;
          slen_sum      += subj. schars. size ();
          scoverage_sum += subj. scoverage;
          length_sum    += subj. length_sum;
          nident_sum    += subj. nident_sum;
        }
        ASSERT (scoverage_sum <= slen_sum);
        ASSERT (nident_sum <= length_sum);
        sseqids. sort ();
        ASSERT (sseqids. isUniq ());

        if (! queryName. empty ())
          /* 0 */ cout << queryName << '\t';
        cout << nvl (qseqid_, na) << '\t' << nvl (qtitle, na);
          //    1                  2
        const ONumber on (cout, 2, false);  // PAR
        cout         
          /* 3*/ << '\t' << qlen
          /* 4*/ << '\t' << qcoverage                
          /* 5*/ << '\t' << double (qcoverage) / double (qlen) * 100.0;
      //if (verbose ())
          cout 
            /* 6*/ << '\t' << length_sum                  
            /* 7*/ << '\t' << nident_sum;
        cout
          /* 8*/ << '\t' << double (nident_sum) / double (length_sum) * 100.0;
        if (! subjectName. empty ())
          cout /* 9 */ << '\t' << subjectName;
      //if (verbose ())
          cout 
            /*10*/ << '\t' << slen_sum
            /*11*/ << '\t' << scoverage_sum;
        cout
          /*12*/ << '\t' << double (scoverage_sum) / double (slen_sum) * 100.0
          /*13*/ << '\t' << sseqids. toString (" ");        
        cout << endl;
      }
      break;
      
    case Mode::all:
      {
        for (const auto& it : sseqid2subject)
        {
          const string& sseqid = it. first;
          const Subject& subj  = it. second;
          const size_t slen = subj. schars. size ();
          for (const Hsp& hsp : subj. hsps)
          {
            if (! queryName. empty ())
              /* 0 */ cout << queryName << '\t';
            cout 
              /* 1 */         << qseqid_ 
              /* 2 */ << '\t' << nvl (qtitle, na)
              /* 3 */ << '\t' << qlen;
            if (! subjectName. empty ())
              cout /* 4 */ << '\t' << subjectName;
            cout << /* 5 */ '\t' << sseqid 
                 << /* 6 */ '\t' << subj. stitle
                 << /* 7 */ '\t' << slen
                 << /* 8 */ '\t' << subj. scoverage
                 << '\t';
            hsp. saveText (cout);
            cout << '\t' << (double) hsp. qLen () / (double) qlen * 100.0
                 << '\t' << (double) hsp. sLen () / (double) slen * 100.0;
            cout << endl;
          }
        }
      }
      break;
      
    default:
      NEVER_CALL;
  }
}



#if 0
void reportMissed (const string &qseqid,
                   const Vector<uint> &qchars)
{
  if (qseqid. empty ())
    return;
    
  ASSERT (! qchars. empty ());

  string qtitle (qseqid);
  const string qseqid_ (findSplit (qtitle, '|'));
  replace (qtitle, '_', ' ');
  trim (qtitle);

  size_t start = 0;    
  FFOR (size_t, i, qchars. size () + 1)
    if (i == qchars. size () || qchars [i])
    {
      if (start < i)
      {
        if (! queryName. empty ())
          /* 0 */ cout << queryName << '\t';
        cout         << qseqid_ 
             << '\t' << nvl (qtitle, na) 
             << '\t' << qchars. size ()
             << '\t' << start + 1
             << '\t' << i
             << '\t' << i - start
             << '\t' << double (i - start) / double (qchars. size ()) * 100.0
             << endl;
      }
      start = i + 1;
    }
}
#endif




// ThisApplication

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print the coverage of a query DNA by subject DNAs, sorted by qseqid")
    {
      version = VERSION;
      string format (XSTR(FORMAT));
      replaceStr (format, " >>", "");
      addPositional ("in", "BLASTN output in format: " + format  + ", sorted by qseqid and nident descending.\n\
qseqid can have a suffix \"|qtitle\"");
      addKey ("mode", "\
combine: combine the coverages by subject DNAs\n\
all: report all covered segments", "all");
// best: print only the best coverage by a subject DNA
// missed: report non-covered query DNA segments
      addKey ("query", "Query DNA name");
      addKey ("subject", "Subject DNA name, used if mode = combine");
      addKey ("pident_min", "Min. percent of identity", "0");
      addKey ("align_min", "Min. alignment length", "0");
      addFlag ("force", "Force one-line report if there is no match");
    }



  void body () const final
  {
    const string inFName     = getArg ("in");
    const string modeS       = getArg ("mode");
                 queryName   = getArg ("query");
                 subjectName = getArg ("subject");
    const double pident_min  = str2<double> (getArg ("pident_min"));
    const size_t align_min   = str2<size_t> (getArg ("align_min"));
    const bool force         = getFlag ("force");
    
    // mode
    ASSERT (mode == Mode::none);
  /*if (modeS == "best")
      mode = Mode::best;
    else*/ if (modeS == "combine")
      mode = Mode::combine;
  //else if (modeS == "missed")
    //mode = Mode::missed;
    else if (modeS == "all")
      mode = Mode::all;
    else
      throw runtime_error ("Unknown mode: " + strQuote (modeS));
    ASSERT (mode != Mode::none);
      
    QC_IMPLY (force, mode == Mode::combine);
              
    
    // Header
    cout << '#';
    if (! queryName. empty ())
      cout << "query\t";
        //     0
    switch (mode)
    {
      case Mode::combine:
        {    
          cout << "qseqid\tqtitle\tqlen\tqcoverage\tpqcoverage";
            //     1       2       3     4          5        
        //if (verbose ())
            cout << "\talign_length\tnident";
              //       6             7  
          cout  << "\tpident";  // 8
          if (! subjectName. empty ())
            cout << "\tsubject";
              //       9
        //if (verbose ())
            cout << "\tslen\tscoverage";
              //       10     11
          cout << "\tpscoverage\tscontigs";  
              //     12          13
        }
        break;
      case Mode::all:
        {    
          cout << "qseqid\tqtitle\tqlen";
            //     1       2       3    
          if (! subjectName. empty ())
            cout << "\tsubject";
              //       4
          cout << "\tsseqid\tstitle\tslen\tscoverage_sum\t";
            //       5       6        7    8
          Hsp::saveHeader (cout);
          cout << "\tpqcoverage\tpscoverage";
        }
        break;
      default: 
        NEVER_CALL;
    }
    cout << endl;
    
    
    if (inFName. empty ())
      return;


    string qseqid_prev;
    size_t qlen_prev = 0;  // >0 => constant
    LineInput f (inFName);
    string qseqid, sseqid;
    string stitle;
    size_t nident_prev = 0;
    while (f. nextLine ())
      try
      {
        replace (f. line, ' ', '_');  // For qtitle, stitle processing
        if (verbose (-1))
          cerr << f. line << endl;
        istringstream iss (f. line);
        size_t length, nident, qstart, qend, qlen, sstart, send, slen;
        slen = 0;
        iss >> FORMAT;
        QC_ASSERT (! sseqid. empty ());
        QC_ASSERT (! qseqid. empty ());
        QC_ASSERT (qseqid_prev <= qseqid);
        QC_ASSERT (qlen);
        QC_ASSERT (slen);
        QC_IMPLY (qseqid_prev == qseqid, qlen == qlen_prev); 
        QC_IMPLY (qseqid_prev == qseqid, nident <= nident_prev); 
        nident_prev = nident;
        
        const Hsp hsp (qstart, qend, sstart, send, length, nident);
          // Invokes: QC_ASSERT
        QC_ASSERT (hsp. qend <= qlen);
        QC_ASSERT (hsp. send <= slen);

        if (length < align_min)
          continue;
        if ((double) nident / (double) length * 100.0 < pident_min)
          continue;

        if (qseqid_prev != qseqid)
        {
          reportSubjects (qseqid_prev, false);  
          sseqid2subject. clear ();
          qseqid_prev = qseqid;
          qlen_prev = qlen;
          qchars. resize (qlen, 0);
          qchars. setAll (0);
        }

      #if 0
        if (   mode == Mode::missed 
            && qseqid_prev != qseqid
           )
          reportMissed (qseqid_prev, qchars);
      #endif
      
        bool novel = false;
        ASSERT (hsp. qend <= qchars. size ());
        FOR_START (size_t, i, hsp. qstart, hsp. qend) 
        {
          if (! qchars [i])
            novel = true;
          qchars [i] ++;
        }
        if (! novel)
        {
          ASSERT (! sseqid2subject. empty ());
          continue;
        }
      
        const Subject* subj = findPtr (sseqid2subject, sseqid);
        if (! subj)
        {
          stitle += "_";
          trimPrefix (stitle, sseqid + "_");  // stitle has sseqid as a prefix
          replace (stitle, '_', ' ');
          trim (stitle);
          sseqid2subject [sseqid] = move (Subject (move (stitle), slen));
          ASSERT (stitle. empty ());
          subj = & sseqid2subject [sseqid];
        }
        ASSERT (subj);
      
        var_cast (subj) -> addHsp (hsp);
      }
      catch (const exception &e)
      {
        throw runtime_error (string (e. what ()) + "\nat line " + to_string (f. lineNum));
      }
  #if 0
    if (mode == Mode::missed)
      reportMissed (qseqid_prev, qchars);
  #endif      
    reportSubjects (qseqid_prev, force);
  }
};




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



