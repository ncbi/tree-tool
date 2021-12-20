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



#define FORMAT  qseqid >> sseqid >> length >> nident >> qstart >> qend >> qlen >> sstart >> send >> slen >> stitle



static const string na ("NA");
string query;
string subject;
      


struct Subject
{
  // 1-based, end - last character position
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
//string plasmidName;
  size_t qcoverage {0};
  size_t scoverage {0};
  
  Subject () = default;
  void saveText (ostream &os,
                 size_t qlen) const
    { QC_ASSERT (qend <= qlen);
      const ONumber on (os, 2, false);  // PAR
      os         
        /* 3*/         << qlen
        /* 4*/ << '\t' << qstart                  
        /* 5*/ << '\t' << qend                    
        /* 6*/ << '\t' << qcoverage                
        /* 7*/ << '\t' << double (qcoverage) / double (qlen) * 100.0
        /* 8*/ << '\t' << length                  
        /* 9*/ << '\t' << nident
        /*10*/ << '\t' << double (nident) / double (length) * 100.0
        /*11*/ << '\t' << (sseqid. empty () ? nvl(subject,"NA") : sseqid)
        /*12*/ << '\t' << slen;
      if (! sseqid. empty ())
      {
        os /*13*/ << '\t' << sstart                  
           /*14*/ << '\t' << send;
      }
      os
        /*15*/ << '\t' << scoverage                
        /*16*/ << '\t' << double (scoverage) / double (slen) * 100.0;
      if (! sseqid. empty ())
      {
        string stitle_ (stitle + "_");
        trimPrefix (stitle_, sseqid + "_");
        replace (stitle_, '_', ' ');
        trim (stitle_);
        os /*17*/ << '\t' << nvl (stitle_, na);
      }
    }
  void qc () const
    { if (! qc_on)
        return;
      QC_ASSERT (! sseqid. empty ());
      QC_ASSERT (qstart < qend);
      QC_ASSERT (sstart < send);
      QC_ASSERT (send <= slen);
      QC_ASSERT (nident <= length);
      QC_ASSERT (qcoverage); 
      QC_ASSERT (scoverage); 
      QC_ASSERT (qcoverage <= qend - qstart + 1);
      QC_ASSERT (scoverage <= send - sstart + 1);
    }
    
  bool operator< (const Subject& other) const
    { LESS_PART (other, *this, qcoverage);
      LESS_PART (other, *this, scoverage);
      LESS_PART (*this, other, sseqid);  // Tie resolution
      return false;
    }
};



void addPrevSubject (Vector<Subject> &subjects, 
                     Subject &subj,
                     const Vector<uint> &qchars,
                     const Vector<uint> &schars)
// Update: subj
{
  if (subj. sseqid. empty ())
    return;
    
  ASSERT (schars. size () == subj. slen);
    
  subj. qcoverage = 0;
  for (const uint c : qchars)
    subj. qcoverage += (bool) c;   

  subj. scoverage = 0;
  for (const uint c : schars)
    subj. scoverage += (bool) c;   
    
  subj. qc ();

  subjects << move (subj);
  subj = Subject ();
}



enum class Mode {best, combine, missed, all};



void processSubjects (const string &qseqid,
                      size_t qlen,
                      Vector<Subject> &subjects,
                      Mode mode,
                      bool force)
{
  IMPLY (qseqid. empty (), subjects. empty ());

  if (mode == Mode::missed)
    return;
    
  if (subjects. empty ())
  {
    if (force)
    {
      Subject s;
      s. qstart = 0;
      s. sstart = 0;
      subjects << move (s);
    }
    else
      return;
  }
  else
    ASSERT (qlen);
    
  if (mode == Mode::combine)
  {
    Subject last (subjects. pop ());  // Has right qcoverage
    last. sseqid. clear ();
    for (const Subject& subj : subjects)
    {
      minimize (last. qstart, subj. qstart);
      maximize (last. qend,   subj. qend);
      last. slen      += subj. slen;
      last. length    += subj. length;
      last. nident    += subj. nident;
      last. scoverage += subj. scoverage;
    }
    subjects. clear ();
    subjects << last;
    ASSERT (subjects. size () == 1);
  }

  subjects. sort ();
  string qtitle (qseqid);
  const string qseqid_ (findSplit (qtitle, '|'));
  replace (qtitle, '_', ' ');
  trim (qtitle);
  for (const Subject& subj : subjects)
  {
    if (! query. empty ())
      /* 0 */ cout << query << '\t';
    cout << qseqid_ << '\t' << nvl (qtitle, na) << '\t';
      //    1                  2
    subj. saveText (cout, qlen);      
    cout << endl;
    if (mode == Mode::best)
      break;
  }
  
  subjects. clear ();
}



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
        if (! query. empty ())
          /* 0 */ cout << query << '\t';
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




// ThisApplication

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print the coverage of a query DNA by subject DNAs, sorted by qcoverage descending")
    {
      version = VERSION;
      string format (XSTR(FORMAT));
      replaceStr (format, " >>", "");
      addPositional ("in", "BLASTN output in format: " + format  + ", sorted by qseqid, sseqid.\n\
qseqid can have a suffix \"|qtitle\"");
      addKey ("mode", "\
best: print only the best coverage by a subject DNA\n\
combine: combine the coverages by subject DNAs\n\
missed: report non-covered query DNA segments\n\
all: report all covered segments", "all");
      addKey ("query", "Query DNA name");
      addKey ("subject", "Subject DNA name, used if mode = combine");
      addKey ("pident_min", "Min. percent of identity", "0");
      addFlag ("force", "Force one-line report if there is no match");
    }



  void body () const final
  {
    const string inFName = getArg ("in");
    const string modeS   = getArg ("mode");
                 query   = getArg ("query");
                 subject = getArg ("subject");
    const double pident_min = str2<double> (getArg ("pident_min"));
    const bool   force   = getFlag ("force");
    
    Mode mode;
    if (modeS == "best")
      mode = Mode::best;
    else if (modeS == "combine")
      mode = Mode::combine;
    else if (modeS == "missed")
      mode = Mode::missed;
    else if (modeS == "all")
      mode = Mode::all;
    else
      throw runtime_error ("Unknown mode: " + strQuote (modeS));
              
    
    cout << '#';
    if (! query. empty ())
      cout << "query\t";
        //     0
    cout << "qseqid\tqtitle\tqlen\tqstart\tqend\tqcoverage\tpqcoverage";
      //     1       2       3     4       5     6          7           
    if (mode != Mode::missed)
    {
      cout << "\talign_length\tnident\tpident\tsseqid\tslen";
        //       8             9       10      11      12
      if (mode != Mode::combine)
        cout << "\tsstart\tsend";
          //       13      14
      cout << "\tscoverage\tpscoverage";
        //       15         16         
      if (mode != Mode::combine)
        cout << "\tstitle";
          //       17
    }
    cout << endl;
    
    
    if (inFName. empty ())
      return;


    string qseqid_prev;
    size_t qlen_prev = 0;  // >0 => constant
    Vector<Subject> subjects;  subjects. reserve (10000);  // PAR
    Subject subj;  // previous
    {
      LineInput f (inFName);
      string qseqid, sseqid;
      string stitle;
      Vector<uint> qchars;    
      Vector<uint> schars;    
      while (f. nextLine ())
        try
        {
          replace (f. line, ' ', '_');  // For stitle processing
          istringstream iss (f. line);
          size_t length, nident, qstart, qend, qlen, sstart, send, slen;
          slen = 0;
          iss >> FORMAT;
          const size_t sstart_ = min (sstart, send);
          const size_t send_   = max (sstart, send);
          QC_ASSERT (! qseqid. empty ());
          QC_ASSERT (! sseqid. empty ());
          QC_ASSERT (qseqid_prev <= qseqid);
          QC_IMPLY (qseqid_prev == qseqid, subj. sseqid <= sseqid);
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
          QC_IMPLY (qseqid_prev == qseqid, qlen_prev == qlen); 
          
          if ((double) nident / (double) length * 100.0 < pident_min)
            continue;

          if (   mode == Mode::missed 
              && qseqid_prev != qseqid
             )
            reportMissed (qseqid_prev, qchars);

          if (! (         qseqid_prev == qseqid
                 && subj. sseqid      == sseqid
                )
             )
          {
            addPrevSubject (subjects, subj, qchars, schars);
            if (! ((   mode == Mode::combine 
                    || mode == Mode::missed
                   )
                   && qseqid_prev == qseqid
                  )
               )
            {
              qchars. resize (qlen, 0);
              qchars. setAll (0);
            }
            schars. resize (slen, 0);
            schars. setAll (0);
          }
          
          if (qseqid_prev != qseqid)
          {
            processSubjects (qseqid_prev, qlen_prev, subjects, mode, false);
            qseqid_prev = qseqid;
            qlen_prev = qlen;
          }
          
          FOR_START (size_t, i, qstart - 1, qend)  // 1-based
            qchars [i] ++;
          FOR_START (size_t, i, sstart_ - 1, send_)  // 1-based
            schars [i] ++;
            
          // Subject attributes
          subj. sseqid = move (sseqid);
          minimize (subj. qstart, qstart);
          maximize (subj. qend, qend);
          minimize (subj. sstart, sstart_);
          maximize (subj. send, send_);
          subj. slen = slen;
          subj. length += length;
          subj. nident += nident;
          subj. stitle = stitle;
        //subj. plasmidName = plasmidName;
        }
        catch (const exception &e)
        {
          throw runtime_error (string (e. what ()) + "\nat line " + to_string (f. lineNum));
        }
      addPrevSubject (subjects, subj, qchars, schars);
      if (mode == Mode::missed)
        reportMissed (qseqid_prev, qchars);
        
      processSubjects (qseqid_prev, qlen_prev, subjects, mode, force);
    }
  }
};




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



