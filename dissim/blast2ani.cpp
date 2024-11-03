// blast2ani.cpp

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
*   Compute ANI-based Jukes-Cantor dissimilarity
*
*/


#undef NDEBUG

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
using namespace DM_sp;
#include "../version.inc"

#include "../common.inc"



#undef INTERVAL_STAT



namespace 
{


inline Real jukesCantor (Real ani,
                         Real p_chance) 
  { ASSERT (ani <= 1);
    ASSERT (p_chance >= 0);
    if (ani <= p_chance)
      return NaN;
    return - log ((ani - p_chance) / (1 - p_chance));  // was: 100 *
  }



inline bool is_nt (char c)
  { return    c == 'A'
           || c == 'C'
           || c == 'G'
           || c == 'T';
  }
  
  
  
inline bool is_nt_M (char c)
  { return    c == 'A'
           || c == 'C';
  }



bool isAmbiguousTriplet (const char triplet [3])
{
  for (size_t i = 0; i < 3; i++)
    if (! is_nt (triplet [i]))
      return true;
  return false;
}



bool isStopCodon (const char triplet [3])
{
  switch (triplet [0])
  {
    case 'T':
      switch (triplet [1])
      {
        case 'A':
          switch (triplet [2])
          {
            case 'A': return true;  // TAA
            case 'G': return true;  // TAG
            default:  return false;
          }
        case 'G':
          switch (triplet [2])
          {
            case 'A': return true;  // TGA
            default:  return false;
          }
        default: return false;
      }
    default: return false;
  } 
}



#if 0  
// not used
bool hasSynonymy (const char triplet [3])
{
  switch (triplet [0])
  {
    case 'A':
      switch (triplet [1])
      {
        case 'T': return false;  // ATA Ile, Met
        default:  return true;
      }
    case 'T':
      switch (triplet [1])
      {
        case 'A':
          switch (triplet [2])
          {
            case 'A': return false;  // stop codon
            case 'G': return false;  // stop codon
            default:  return true;
          }
        case 'G':
          switch (triplet [2])
          {
            case 'A': return false;  // srop codon
            case 'G': return false;  // Trp
            default:  return true;
          }
        default: return true;
      }
    default: return true;
  } 
}
#endif



bool hasRedundant3 (const char triplet [3])
{
  switch (triplet [0])
  {
    case 'A':
      switch (triplet [1])
      {
        case 'C': return true;  // Thr
        default:  return false;
      }
    case 'C':
      switch (triplet [1])
      {
        case 'C': return true;  // Pro
        default:  return false;
      }
    case 'G':
      switch (triplet [1])
      {
        case 'C': return true;  // Ala
        case 'G': return true;  // Gly
        case 'T': return true;  // Val
        default:  return false;
      }
    default: return false;
  } 
}



void complementNt (char &c)
{
  switch (c)
  {
    case 'A': c = 'T'; break;
    case 'C': c = 'G'; break;
    case 'G': c = 'C'; break;
    case 'T': c = 'A'; break;
  }
}



string reverseDna (const string &seq)
{
  string rev (seq);
  for (size_t i = 0; i < rev. size () / 2; i++)
    swap (rev [i], rev [rev. size () - 1 - i]);
  for (char& c : rev)
    complementNt (c);
  return rev;
}



struct SimpleAnnotation
{
  struct Codon
  {
    size_t pos3;
      // 3rd codon position in seq_arg
    bool toCount;
      // non-ambiguous
    Codon (size_t pos3_arg,
           bool toCount_arg)
      : pos3 (pos3_arg)
      , toCount (toCount_arg)
      {}
    bool operator< (const Codon &other) const
      { return pos3 < other. pos3; }
  };
  typedef  vector<Codon>  Orf;  
  
  struct Codon3
  {
    size_t orfLen;
    bool strand;
    bool toCount;
      // non-ambiguous
    Codon3 (size_t orfLen_arg,
            bool strand_arg,
            bool toCount_arg)
      : orfLen (orfLen_arg)
      , strand (strand_arg)
      , toCount (toCount_arg)
      { ASSERT (orfLen); }
    Codon3 ()
      : orfLen (0)
      , strand (false)
      , toCount (false)
      {}
  };
  vector<Codon3> codon3_vec;
    // Index: position in seq
    // Only 1 Codon3 per position in seq: best strand
  
  
  explicit SimpleAnnotation (const string &seq_arg)
    : codon3_vec (seq_arg. size ())
    {
      ASSERT (! seq_arg. empty ());
      for (const bool strand : {false, true})
      {
        const string seq (strand ? seq_arg : reverseDna (seq_arg));
        for (size_t frame = 0; frame < 3; frame++)
        {
          size_t codonPos = 0;
          char triplet [3] = {' ', ' ', ' '};
          Orf orf;
          for (size_t i = frame; i < seq. size (); i++)
          {
            const char c = seq. at (i);
            if (c == '-')
              continue;
            triplet [codonPos] = c;
            if (codonPos == 2)  // 3rd codon position, triplet[] is filled
            {
              if (isStopCodon (triplet) || isAmbiguousTriplet (triplet))
              {
                orf2vec (orf, seq_arg, strand);
                orf. clear ();
              }
              else
                orf. push_back (Codon (strand ? i : (seq. size () - 1 - i), ! isAmbiguousTriplet (triplet) && hasRedundant3 (triplet)));
            }
            codonPos++;
            codonPos %= 3;
          }
          orf2vec (orf, seq_arg, strand);
        }
      }
    }
private:
  void orf2vec (const Orf &orf,
                const string &seq,
                bool strand)
    // Update: codon3_vec
    {
      ASSERT (! seq. empty ());
      if (orf. size () < 120)  // PAR
        return;
      for (const Codon& codon : orf)
      {
        ASSERT (seq [codon. pos3] != '-');
        //
        if (codon. pos3 == 0)
          continue;
        size_t pos3_prev = codon. pos3 - 1;
        while (pos3_prev > 0 && seq [pos3_prev] == '-')
          pos3_prev--;
        ASSERT (seq [pos3_prev] != '-');
        //
        if (codon. pos3 == seq. size () - 1)
          continue;
        size_t pos3_next = codon. pos3 + 1;
        while (pos3_next < seq. size () - 1 && seq [pos3_next] == '-')
          pos3_next++;
        ASSERT (seq [pos3_next] != '-');
        // ORFs in all frames around codon. pos3
        if (orf. size () <= max ( max ( codon3_vec [codon. pos3]. orfLen
                                      , codon3_vec [pos3_prev].   orfLen
                                      )
                                ,       codon3_vec [pos3_next].   orfLen
                                )
           )
          continue;
        // Codons in wrong frames are suppressed
        codon3_vec [pos3_prev] = Codon3 ();
        codon3_vec [pos3_next] = Codon3 ();
        // Codon in the right frame
        codon3_vec [codon. pos3] = Codon3 (orf. size (), strand, codon. toCount);
      }
    }
public:
  void print () const
    {
      for (size_t i = 0; i < codon3_vec. size (); i++)
      {
        const Codon3& c = codon3_vec [i];
        if (c. orfLen && c. toCount)
          cout << ' ' << i + 1;
      }
      cout << endl;
    }
};



struct Alignment
{
  const string& seq1;
  const string& seq2;
  SimpleAnnotation ann1;
  SimpleAnnotation ann2;
  vector<bool> dense;

  
  Alignment (const string &seq1_arg,
             const string &seq2_arg)
    : seq1 (seq1_arg)
    , seq2 (seq2_arg)
    , ann1 (seq1)
    , ann2 (seq2)
    , dense (seq1_arg. size (), false)
    {
      ASSERT (! seq1. empty ());
      ASSERT (seq1. size () == seq2. size ());
      setDense (25);  // PAR  
    }
private:
  void setDense (size_t windowSize)
  // Density filtering
  // Update: dense[]
  {
    ASSERT (dense. size () == seq1. size ());
    const size_t window_mismatches_max = getMaxMismatches (windowSize);
    if (window_mismatches_max > windowSize)
      return;
    if (window_mismatches_max <= 1)
    {
      dense. assign (dense. size (), true);
      return;
    }
    size_t windowMismatches = 0;  
    vector<bool> window (windowSize, false);  // true <=> mismatch
    size_t windowPos = 0;
    for (size_t i = 0; i < seq1. size (); i++)
    {
      const bool mismatch = (seq1 [i] != seq2 [i]);
      if (mismatch)
        windowMismatches++;
      if (window [windowPos])  // Outside the current window
      {
        ASSERT (windowMismatches);
        windowMismatches--;
      }
      window [windowPos] = mismatch;
      windowPos++;
      windowPos %= windowSize;
      ASSERT (windowMismatches <= windowSize);
      if (windowMismatches < window_mismatches_max)
        continue;
    //cout << "density at: " << i + 1 << endl;
      for (size_t j = i + windowSize; j + 2 * windowSize > i; j--)
      {
        if (j < dense. size ())
        {
          if (dense [j])
            break;
          dense [j] = true;
        }
        if (j == 0)
          break;
      }
    }
  }
  size_t getMaxMismatches (size_t windowSize) const
  // Return: 0 .. windowSize+1
  {
    ASSERT (windowSize);
    size_t mismatches = 0;  // global mismatches
    for (size_t i = 0; i < seq1. size (); i++)
      if (seq1 [i] != seq2 [i])
        mismatches++;
    ASSERT (mismatches <= seq1. size ());
    if (! mismatches)
      return windowSize + 1;
    const Real pValue = 1e-6;  // PAR
    // Poisson distribution
    const Real lambda = (Real) windowSize * (Real) mismatches / (Real) seq1. size ();
    ASSERT (lambda > 0);    
    size_t k = 0;
    Real p = exp (- lambda);
    while (p > pValue && k <= windowSize)
    {
      k++;
      p *= lambda / (Real) k;
    //cout << "Poisson[" << k << "] = " << p << endl; 
    }
    ASSERT (p >= 0);
    ASSERT (p < 1);
    ASSERT (k <= windowSize + 1);
  //cout << lambda << ' ' << k << endl;  
    return k;
  }
public:
    
    
  bool getCodonPos3 (size_t i,
                     char &c1,
                     char &c2) const
    // Return: success
    // Output: c1, c2; if Return
    {
      ASSERT (i < seq1. size ());
      if (dense [i])
        return false;  
      const SimpleAnnotation::Codon3 codon1 = ann1. codon3_vec [i];
      const SimpleAnnotation::Codon3 codon2 = ann2. codon3_vec [i];
      if (   codon1. orfLen && codon1. toCount 
          && codon2. orfLen && codon2. toCount
          && codon1. strand == codon2. strand
         )
      {
        c1 = seq1 [i];
        c2 = seq2 [i];
        ASSERT (is_nt (c1));
        ASSERT (is_nt (c2));
        if (! codon1. strand)
        {
          complementNt (c1);
          complementNt (c2);
        }
        return true;
      }
      return false;
    }
  size_t getDensitySize () const
    {
      size_t n = 0;
      for (const bool c : dense)
        if (c)
          n++;
      return n;
    }
};



struct Segment
{
  string seqid;
  size_t start;
  size_t end;
  bool strand;
  string seq;
    // DNA alphabet, may contain '-'
    
    
  Segment ()
    : start (0)
    , end (0)
    , strand (false)
    {}
  Segment (const Segment &other)
    : seqid  (other. seqid)
    , start  (other. start)
    , end    (other. end)
    , strand (other. strand)
    , seq    (other. seq)
    {}
  void finish ()
    {
      ASSERT (start);
      ASSERT (end);
      ASSERT (start != end);
      strand = (start < end);
      if (strand)
        start--;
      else
        end--;
    }

  void qc () const
    {
      ASSERT (! seqid. empty ());
      ASSERT (! contains (seqid, ' '));
      ASSERT (start != end);
      ASSERT (strand == (start < end));
      IMPLY (! seq. empty (), size () <= seq. size ());
      for (const char c : seq)
        ASSERT (c == (char) toupper (c));
    }


  void saveName (ostream &os) const
    { os << seqid << ":" << start + 1 << '-' << end; }
  size_t size () const
    { return max (start, end) - min (start, end); }
  void setRelEnd (size_t len)
    { 
      if (strand)
        end = start + len;
      else
      {
        ASSERT (start >= len);
        end = start - len;
      }
    }
};



struct Hsp
{
  Vector<Segment> segs;
  Real p_chance;
  
  size_t len;
    // # aligned nucleotides
  size_t nident;
  
  size_t len_sum;
  Real jc_sum;
  Real jc2_sum;
  
  
  Hsp (const Vector<Segment> &segs_arg,
       Real p_chance_arg,
       size_t len_arg,
       size_t nident_arg)
    : segs (segs_arg)
    , p_chance (p_chance_arg)
    , len (len_arg)
    , nident (nident_arg)
    , len_sum (0)
    , jc_sum (0)
    , jc2_sum (0)
    { ASSERT (p_chance <= 1);
      ASSERT (p_chance >= 0);
      ASSERT (len);
      ASSERT (nident <= len);
    }
  void print (ostream &os) const
    { segs [0]. saveName (os);
      os << '/';
      segs [1]. saveName (os);
      os << ' ' << len 
         << ' ' << nident
         << ' ' << ani ()
         << ' ' << jukesCantor ()
         << ' ' << jc_ave ()
         << ' ' << jc_sd ()
         << endl;
    }


  Real ani () const
    // "Average Nucleotide Identity"
    { return (Real) nident / (Real) len; }
  Real jukesCantor () const
    { return ::jukesCantor (ani (), p_chance); }

  void setSum (const Hsp* prev)
    // Requires: sort()'ed
    { 
      const Real jc = jukesCantor ();
      len_sum = len;
      jc_sum =       jc  * (Real) len;
      jc2_sum = sqr (jc) * (Real) len;
      if (prev)
      {
        len_sum += prev->len_sum;
        jc_sum  += prev->jc_sum;
        jc2_sum += prev->jc2_sum;
      } 
    }
  Real jc_ave () const
    { return jc_sum / (Real) len_sum; }
  Real jc_sd () const
    { return sqrt (jc2_sum / (Real) len_sum - sqr (jc_ave ())); }
  bool isOutlier (Real sds) const
    { const Real jc = jukesCantor ();
      return isNan (jc) || jc > jc_ave () + sds * jc_sd (); 
    }
};



bool compareMatches (const Hsp &h1,
                     const Hsp &h2)
{
  return h1. ani () > h2. ani ();
}




struct ThisApplication final : Application
{
  ThisApplication ()
  : Application ("Compute ANI-based Jukes-Cantor dissimilarity.\n\
Console input: BLASTN HSP: 'qseqid sseqid qstart qend sstart send nident length qseq sseq' (capital letters)\n\
Print: <jc_whole> [<jc_ave> <jc_sd/jc_ave>] <jc_cluster_ave>")
  { 
    version = VERSION;
	  addPositional ("piece_len", "Minimum piece length");
	  addFlag ("cut", "cut into pieces of <piece_len>: if true then 'qseq sseq' should be absent, only 3rd codon positions in reduced alphabet (M/not M) are used");
	  addFlag ("print_matches", "print matching nt pairs (if <cut>)");
	  addKey ("pieces_file", "Save data of HSP pieces into a file", "");
	}



	void body () const final
	{
    const size_t pieceLen    = str2<size_t> (getArg  ("piece_len"));   
    const bool cut           = getFlag ("cut");  
    const bool printMatches  = getFlag ("print_matches");  
    const string pieces_file = getArg ("pieces_file");  
    ASSERT (pieceLen);
    ASSERT (printMatches <= cut);
        
    
    const Real p_chance = cut ? 0.5 : 0.25;  // PAR


    vector<Hsp> hsps;
    map <pair<char,char>, size_t> pair2num;
  #ifdef INTERVAL_STAT
    const size_t intervalLen = 25;  
    vector<size_t> intervalMatches (intervalLen + 1, 0); 
  #endif
    while (! cin. eof ())
    {
      // BLASTN HSP
      Vector<Segment> segs (2);
        // 0 - query, 1 - subject
      size_t nident;
      size_t length;
  
      cin >> segs [0]. seqid
          >> segs [1]. seqid
          >> segs [0]. start
          >> segs [0]. end
          >> segs [1]. start
          >> segs [1]. end
          >> nident 
          >> length;
      if (cut)
        cin >> segs [0]. seq
            >> segs [1]. seq; 
      if (segs [0]. seqid. empty ())
        break;          
        
      if (length < pieceLen)
        continue;
        
      for (const bool b : {false, true})
      {
        segs [b]. finish ();
        segs [b]. qc ();
        IMPLY (cut, segs [b]. seq. size () == length);
      }
      ASSERT (segs [0]. strand);
      ASSERT (nident > 0);


      if (cut)
      {
        const Alignment al (segs [0]. seq, segs [1]. seq);
        Vector<Segment> pieces (segs);
        for (const bool b : {false, true})
          pieces [b]. seq. clear ();
      #ifdef INTERVAL_STAT
        size_t intervalPos = 0;
        size_t intervalMismatches = 0;
      #endif
        for (size_t i = 0; i + pieceLen <= length; i += pieceLen)  // tail is ignored; different results if pieceLen = 10000 ??
        {
        #if 0
          cout << endl;
          for (size_t j = i; j < i + pieceLen && j < length; j++) 
            cout << seg [0]. seq [j];
          cout << endl;
          for (size_t j = i; j < i + pieceLen && j < length; j++) 
            cout << seg [1]. seq [j];
          cout << endl;
        #endif
          size_t nidentAl = 0;
          size_t len = 0;
          size_t segLen [2/*bool*/] = {0, 0};
          for (size_t j = i; j < i + pieceLen && j < length; j++) 
          {
            for (const bool b : {false, true})
              if (segs [b]. seq [j] != '-')
                segLen [b] ++;
                
            char c1 = ' ';
            char c2 = ' ';
            if (al. getCodonPos3 (j, c1, c2))
            {
              if (c1 > c2)
                swap (c1, c2);
              ASSERT (is_nt (c1));
              ASSERT (is_nt (c2));
              pair2num [pair<char,char> (c1, c2)] ++;
              len++;
              if (is_nt_M (c1) == is_nt_M (c2))
              {
                nidentAl++;
              //cout << 1;  
              }
            #if 0
              else
                cout << 0; 
            #endif
            }
          #if 0
            else
              cout << ' ';
          #endif
          #ifdef INTERVAL_STAT
            if (seg [0]. seq [j] != seg [1]. seq [j])
              intervalMismatches++;
            if (intervalPos == intervalLen - 1)
            {
              ASSERT (intervalMismatches <= intervalLen);
              intervalMatches [intervalMismatches]++;
              intervalPos = 0;
              intervalMismatches = 0;
            }
            else
              intervalPos++;
          #endif
          }
        #if 0
          cout << endl;  
          exit (2); 
        #endif
          for (const bool b : {false, true})
          {
          #if 0
            cout << "Before: ";
            pieces [b]. saveName (cout);
            cout << "  by " << segLen [b] << "  "; 
          #endif
            pieces [b]. setRelEnd (segLen [b]);
          #if 0
            cout << "  After: ";
            pieces [b]. saveName (cout);
            cout << "   ";
          #endif
            pieces [b]. qc ();
          }
          if (len >= 100)  // PAR
            hsps. push_back (Hsp (pieces, p_chance, len, nidentAl));
       // al. getDensitySize (): compute conservation ??!
          for (const bool b : {false, true})
            pieces [b]. start = pieces [b]. end;
        }
      }
      else
      {
        size_t gaps [2/*bool*/];
        for (const bool b : {false, true})
        {
          const size_t seg_size = segs [b]. size ();
          ASSERT (length >= seg_size);
          gaps [b] = length - seg_size; 
        }
        hsps. push_back (Hsp (segs, p_chance, length - (gaps [false] + gaps [true]), nident));
      }
    }
    ASSERT (! hsps. empty ());


    size_t len = 0;
    {
      // Average, SD
      size_t nident = 0;
      size_t jc_len = 0;
      Real s = 0;
      Real s2 = 0;
      for (const Hsp& hsp : hsps)
      {
        nident += hsp. nident;
        len    += hsp. len;
        const Real jc = hsp. jukesCantor ();
        if (! isNan (jc))
        {
          jc_len +=                      hsp. len;
          s      +=  jc       * (Real) hsp. len;
          s2     += (jc * jc) * (Real) hsp. len;
        }
      //cout << hsp. nident << ' ' << hsp. len << endl;  
      }
      const Real ani = (Real) nident / (Real) len;
      const Real jc_whole = jukesCantor (ani, p_chance);
      cout << jc_whole;
      if (cut)
      {
        const Real jc_ave = s / (Real) jc_len;
        const Real jc_sd = sqrt (s2 / (Real) jc_len - jc_ave * jc_ave);
        cout << ' ' << jc_ave << ' ' << jc_sd / jc_ave;
      }
    }


    // Cluster
  //if (cut)
    {
      sort (hsps. begin (), hsps. end (), compareMatches);
    #if 0  // A little better for Mycobacterium tuberculosis group 
      size_t cluster_start = 0;
      size_t cluster_stop = 0;
      {
        Real jc_range_min = (Real) len;
        size_t i = 0;
        size_t clusterLen = 0;
        for (size_t j = 0; j < hsps. size (); j++)
        {
          const Hsp& hsp2 = hsps. at (j);
          clusterLen += hsp2. len;
          while ((Real) clusterLen / (Real) len >= 0.75)  // PAR
          {
            const Hsp& hsp1 = hsps. at (i);
            const Real jc_range =   hsp2. jukesCantor () 
                                    - hsp1. jukesCantor ();
            if (jc_range_min > jc_range)
            {
              jc_range_min = jc_range;
              cluster_start = i;
              cluster_stop  = j;
            }
            clusterLen -= hsp1. len;
            i++;
          }
        }
      }
      // Cluster average
      size_t clusterLen = 0;
      Real s = 0;
      for (size_t i = cluster_start; i <= cluster_stop; i++)
      {
        const Hsp& hsp = hsps. at (i);
        clusterLen += hsp. len;
        s += hsp. jukesCantor () * (Real) hsp. len;
      }
      const Real jc_cluster_ave = s / (Real) clusterLen;
    #else
      const Hsp* prev = nullptr;
      for (Hsp& hsp : hsps)
      {
        hsp. setSum (prev);
        prev = & hsp;
      }
      if (! pieces_file. empty ())
      {
        OFStream f (pieces_file);
        ONumber on (f, 6, false);
        for (const Hsp& hsp : hsps)
          hsp. print (f);  
      }
      // Trimming outliers
      size_t stop = hsps. size () - 1; 
      while (stop > 0 && hsps [stop]. isOutlier (4.0))  // PAR
        stop--;
      if (verbose ())
        cout << endl << "jc_max = " << stop << ' ' << hsps [stop]. jukesCantor () << endl; 
      const Real jc_cluster_ave = hsps [stop]. jc_ave ();  
      ASSERT (jc_cluster_ave >= 0);   
    #endif
      cout << ' ' << jc_cluster_ave;
    }
    cout << endl;


    if (cut && printMatches)
      for (const auto& it : pair2num)
        cout << it. first. first << ' ' << it. first. second << ' ' << it. second << endl;
        
  #ifdef INTERVAL_STAT
    for (size_t i = 0; i <= intervalLen; i++)
      cout << i << ' ' << intervalMatches [i] << endl;
  #endif
  }
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



