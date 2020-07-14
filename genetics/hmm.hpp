// hmm.hpp

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


#ifndef HMM_HPP
#define HMM_HPP


#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
#include "../dm/matrix.hpp"
using namespace DM_sp;
#include "seq.hpp"
using namespace Seq_sp;



namespace Hmm_sp
{
// Protein HMM
  
  
  
typedef  char  AminoAcid;
  // In peptideAlphabet or '-'
  
  
struct Hmm;

  
  
struct Profile
// One position in protein sequence
{
  friend Hmm;

  static constexpr const size_t aa_size = 20 + 1;
  array<Prob,aa_size> prob;  
    // AminoAcid distribution
private:
  Matrix transition;
    // State transition probabilities
  MVector state;
    // Probabilities of states
    // States: m, i, d
public:
  array<Real,aa_size> dist2;

  
  explicit Profile (const string &line);
  Profile () = default;
private:
  static Prob str2prob (const string &s);
    // Return: !isNan()
  void finish (const Profile* prev);
    // Output: state, dist2
public:
    

private:    
  static size_t aa2index (char aa);
    // Return: < aa_size; no_index <=> bad aa
public:
  Real getDist2 (AminoAcid c1,
                 AminoAcid c2) const;
    // Return: NaN <=> bad c1 or c2
};



struct Hmm : Named
// COMPO line is required
{
  map<string/*field*/,string/*content*/> meta;
  Profile compo;
  Vector<Profile> profiles;
    // For each position
  
  
  explicit Hmm (LineInput &li);
  Hmm () = default;
    

  Real getDist2 (const Peptide &pep1,
                 const Peptide &pep2) const;
    // Input: pep1, pep2: sparse, size() = profiles.size()
};



typedef  map<string/*Hmm::name*/,Hmm>  HmmLib;
HmmLib loadHmmLib (const string &fName);



struct Hmmsearch : Root
{
private:
  LineInput li;
  Istringstream iss;
public:	
  string prot_name, prot_accession
       , hmm_name,  hmm_accession;
  double eValue1, score1, bias1
       , eValue2, score2, bias2;
  

  explicit Hmmsearch (const string &fName)
    // Input: fName: hmmsearch --tblout result
    : li (fName)
    {}
  

  bool next ();
    // Return: !li.eof
};



struct HmmDom : Root
{
private:
  LineInput li;
  Istringstream iss;
public:
  string prot_name, prot_accession, hmm_name, hmm_accession;
  size_t prot_len, hmm_len, n, of, hmm_from, hmm_to, ali_from, ali_to, env_from, env_to;
  double eValue, full_score, full_bias, cValue, i_eValue, domain_score, domain_bias, accuracy;
  

  explicit HmmDom (const string &fName)
    // Input: fName: hmmsearch --domtblout result
    : li (fName)
    {}
  

  bool next ();
    // Return: !li.eof
};



}  // namespace




#endif
