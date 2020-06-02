// hmm.cpp

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
*   HMM utilities
*
*/


#undef NDEBUG
#include "../common.inc"

#include "hmm.hpp"




namespace Hmm_sp
{



constexpr const Real delta = 1e-4;  // PAR



// Profile

Profile::Profile (const string &line) 
: transition (3, 0.0)
, state (3, 0.0)
{
  istringstream iss (line);
  string val;
  Prob sum = 0;
  ASSERT (aa_size - 1 == 20);
  FOR (size_t, i, aa_size - 1)
  {
    ASSERT (! iss. eof ());
    iss >> val;
    const Prob p = str2prob (val);
    sum += p;
    prob [i] = p;
  }
  ASSERT_EQ (sum, 1, delta); 
}



Prob Profile::str2prob (const string &s)
{
  ASSERT (! s. empty ());
  if (s == "*")
    return 0;
  const Prob p = exp (- atof (s. c_str ()));
  ASSERT (isProb (p)); 
  return p;
}



void Profile::finish (const Profile* prev)
{
  FOR (size_t, i, 3)
    ASSERT_EQ (transition. sumRow (false, i), 1.0, delta);
  
  // state  
  if (prev)
    state. multiply ( true
                    , prev->state, true
                    , prev->transition, false
                    );
  else
    state. copyRow ( true, 0
                   , transition, false, 0
                   );
  ASSERT (state. min () >= 0.0);
  ASSERT_EQ (state. sum (), 1, delta);
  
  // prob[]
  prob [aa_size - 1] = state [2];  // State "d", aa = '-'
  const Prob match = 1 - state [2];
  FOR (size_t, i, aa_size - 1)
    prob [i] *= match;

  array<Real,aa_size> vars;
  Real var_sum = 0;
  Prob prob_sum = 0;    
  FOR (size_t, i, aa_size)
  {
    prob_sum += prob [i];
    const Real var = prob [i] * (1 - prob [i]);
    vars [i] = var;
    var_sum += var;
  }
  ASSERT (var_sum > 0.0);
  ASSERT_EQ (prob_sum, 1.0, delta);  

  // dist2[]
  FOR (size_t, i, aa_size)
  {
    const Real var = vars [i] /* / var_sum*/;  // ??
    // dist = (1-mean)/SD - (0-mean)/SD = 1/SD.
    // dist2 = 1/var
    dist2 [i] = 1 / var;  
  }
}



size_t Profile::aa2index (char aa)
{
  ASSERT (aa > 0);
  
  static constexpr const size_t aa2index_size = 128;
  static array<size_t,aa2index_size /*amino acid ASCII*/> aa2index_;
  static bool initialized = false;
  if (! initialized)
  {
    FOR (size_t, i, aa2index_size)
      aa2index_ [i] = NO_INDEX;
    const string aaVec ("ACDEFGHIKLMNPQRSTVWY-");
    ASSERT (aaVec. size () == aa_size);
    FFOR (size_t, i, aa_size)
    {
      const unsigned char c = (unsigned char) aaVec [i];
      ASSERT (c < aa2index_size);
      aa2index_ [c] = i;
    }
    initialized = true;
  }

  const size_t i = (size_t) aa;
  if (i >= aa2index_size)
    return NO_INDEX;
  return aa2index_ [i];
}



Real Profile::getDist2 (AminoAcid c1,
                        AminoAcid c2) const
{
  const size_t i1 = aa2index (c1);
  const size_t i2 = aa2index (c2);
  if (   i1 == NO_INDEX
      || i2 == NO_INDEX
     )
    return NaN;
  
  if (c1 == c2)
    return 0;

  ASSERT (i1 != i2);    
  return   dist2 [i1]  
         + dist2 [i2];
}




// Hmm

Hmm::Hmm (LineInput &li)
{
  static const string eoh ("//");
  string field;
  
  while (li. nextLine ()) 
  {
    replace (li. line, '\t', ' ');
    field = findSplit (li. line, ' ');
    trim (li. line);
    if (field == "HMM") 
      break;
    meta [field] = li. line;
  }
  
  name = meta ["NAME"];
  ASSERT (goodName (name));
  
  const size_t len = (size_t) atoi (meta ["LENG"]. c_str ());
  ASSERT (len);
  
  ASSERT (meta ["ALPH"] == "amino");
  
//cerr << "Loading " << name << endl;  

  // Model started
  // compo, profiles
  profiles. reserve (len);
  EXEC_ASSERT (li. nextLine ());  // m->m     m->i     m->d     i->m     i->i     d->m     d->d  
  size_t phase = 0;  
  int num = 0;
  Istringstream iss;
  string m2m, m2i, m2d, i2m, i2i, d2m, d2d;
  bool compoExists = false;
  while (li. nextLine () && field != eoh) 
  {
  //cerr << "Line " << num + 1 << endl;  
    replace (li. line, '\t', ' ');
    trim (li. line);
    Profile* prev = nullptr;
    switch (phase)  
    {
      case 0:  // Match emissions
        {
          field = findSplit (li. line, ' ');
          if (field == eoh)
            break;
          if (field == "COMPO")
          {
            ASSERT (! compoExists);
            compoExists = true;
            compo = Profile (li. line);
          }
          else
          {
            const int num1 = atoi (field. c_str ());
            ASSERT (num1 == num + 1);
            num = num1;
            profiles << Profile (li. line);
          }
        }      
        break;
      case 1:  // Insert emissions
        break;  
      case 2:  // State transitions
        {
          iss. reset (li. line);
          iss >> m2m >> m2i >> m2d >> i2m >> i2i >> d2m >> d2d;
          ASSERT (iss. eof ());
          if (! compoExists)
            throw runtime_error ("HMM COMPO line is required");
          Profile& p = profiles. empty () ? compo : profiles. back ();
          p. transition. put (false, 0, 0, Profile::str2prob (m2m));
          p. transition. put (false, 0, 1, Profile::str2prob (m2i));
          p. transition. put (false, 0, 2, Profile::str2prob (m2d));
          p. transition. put (false, 1, 0, Profile::str2prob (i2m));
          p. transition. put (false, 1, 1, Profile::str2prob (i2i));
          p. transition. put (false, 2, 0, Profile::str2prob (d2m));
          p. transition. put (false, 2, 2, Profile::str2prob (d2d));
          p. finish (prev);
          prev = & p;
        }
        break;
      default:
        ERROR;
    }
    phase++;
    phase %= 3;
  }  
  ASSERT (field == eoh);
  ASSERT (profiles. size () == len);
  ASSERT ((size_t) num == len);
}



Real Hmm::getDist2 (const Peptide &pep1,
                    const Peptide &pep2) const
{
  ASSERT (pep1. sparse);
  ASSERT (pep2. sparse);
  ASSERT (pep1. seq. size () == pep2. seq. size ());
  ASSERT (pep1. seq. size () == profiles. size ());
  
  Real dist2 = 0;
  FFOR (size_t, i, profiles. size ())
  {
    const Real d = profiles [i]. getDist2 ( pep1. seq [i]
                                          , pep2. seq [i]
                                          );
    if (! isNan (d))
      dist2 += d;
  }
  return dist2;
}



//

HmmLib loadHmmLib (const string &fName)
{
  HmmLib hmms;
  LineInput li (fName);
  Progress prog;
  while (! li. eof)
  {
    Hmm hmm (li);
    prog (hmm. name);
    hmms [hmm. name] = move (hmm);
  }
  
  return hmms;
}




// HmmSearch

bool Hmmsearch::next ()
{
  while (li. nextLine ())
  {
    if (li. line. empty () || li. line [0] == '#')
      continue;
      
    iss. reset (li. line);
    //     prot name         accession  query name        accession   E-value  score     bias     E-value  score  bias   exp reg clu  ov env dom rep inc description of prot
    iss >> prot_name >> prot_accession 
        >> hmm_name  >> hmm_accession 
        >> eValue1 >> score1 >> bias1
        >> eValue2 >> score2 >> bias2;

  //ASSERT (score1 > 0);
  //ASSERT (score2 > 0);

    ASSERT (! prot_name. empty ());
    ASSERT (! contains (prot_name, ','));

    ASSERT (! hmm_accession. empty ());
    ASSERT (! contains (hmm_accession, ','));

    if (prot_accession == "-")  prot_accession. clear ();
    if (hmm_name == "-")        hmm_name.       clear ();

    return true;
  }
  return false;
}



// HmmDom

bool HmmDom::next ()
{
  while (li. nextLine ())
  {
    if (li. line. empty () || li. line [0] == '#')
      continue;
      
    iss. reset (li. line);
    //      target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
    iss >> prot_name >> prot_accession >> prot_len >> hmm_name >> hmm_accession >> hmm_len 
        >> eValue >> full_score >> full_bias 
        >> n >> of >> cValue >> i_eValue >> domain_score >> domain_bias 
        >> hmm_from >> hmm_to >> ali_from >> ali_to >> env_from >> env_to
        >> accuracy;
    hmm_from--;
    ali_from--;
    env_from--;

    ASSERT (! prot_name. empty ());
    ASSERT (! contains (prot_name, ','));

    ASSERT (! hmm_name. empty ());
    ASSERT (! contains (hmm_name, ','));

    if (prot_accession == "-")  prot_accession. clear ();
    if (hmm_name == "-")        hmm_name.       clear ();

    return true;
  }
  return false;
}




}  // namespace
