// combinatorics.hpp

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
*   Combinatorics utilities
*
*/


#ifndef COMBINATORICS_HPP_5843  // random number
#define COMBINATORICS_HPP_5843


#include "common.hpp"
using namespace Common_sp;



namespace Cmb_sp
{


// Usage:  Obj obj; do <use obj>; while (obj. next ());

bool next (vector<bool> &v);
  // Return: false <=> search is finished
  // Update: v
  // Search size = 2^|v|

bool next (vector<size_t> &indexes,
           const vector<size_t> &indexes_max);
  // Return: false <=> search is finished
  // Update: indexes: 0 <= indexes[i] <= indexes_max[i]
  // Search size = \prod_i (indexes_max[i] + 1)



struct SubsetSearch
// Search size = choice(whole_size,subset_size)
{
  Vector<size_t> subset;
    // Indexes of elements
    // size() = subset_size + 1
    // subset.last() = whole_size: fictitious element
    // i <= subset[i] < subset[i+1]

  SubsetSearch (size_t whole_size,
                size_t subset_size);

  bool next ();
    // Update: subset
    // Return: false <=> search is finished
  void complement (Vector<size_t> &disjoint) const;
    // Output: disjoint = complement of subset
};



struct MinSubsetSearch
// Search time = O(2^n)
// Usage: MinSubsetSearch mss; ASSERT(!mss.complete()); while (mss.next()) <use mss.subset>; ASSERT(mss.subset.empty());
{
  const size_t n;
    // Set size
  Vector<size_t> subset;
    // subset[i] < subset[i+1] < n


  explicit MinSubsetSearch (size_t n_arg)
    : n (n_arg)
    { subset. reserve (n); }
  void saveText (ostream &os) const
    { for (const size_t i : subset)
        os << ' ' << i;
    }


  virtual bool complete (size_t exclude = no_index) const = 0;
    // Input: subset
    //        exclude != no_index => subset element exclude must be excluded
    // Requires: complete(), elements are added to subset => complete()
    //           subset.empty() => !complete()
  bool next ();
    // DFS
    // Return: false <=> search is finished
    //         true <=> subset is a minimal complete subset
private:
  bool advance ();
  bool advance_wide ();
public:
  static void test (ulong seed);
};



struct Permute
// Search size = n!
{
  Vector<size_t> vec;
    // Indexes of elements
    // Init: 0, 1, 2, ...

  explicit Permute (size_t n);

  bool next ()
    { return next_permutation (vec. begin (), vec. end ()); }
    // Update: vec (in lexicographic order)
    // Return: false <=> search is finished
};



}



#endif
