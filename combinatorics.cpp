// combinatorics.cpp

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


#undef NDEBUG

#include "combinatorics.hpp"

#include "common.inc"



namespace Cmb_sp
{


bool next (vector<bool> &v)
{
  const size_t s = v. size ();

  size_t i = 0;
  while (i < s && v [i])
  {
    v [i] = false;
    i++;
  }
  if (i == s)
    return false;

  v [i] = true;

  return true;
}



bool next (vector<size_t> &indexes,
           const vector<size_t> &indexes_max)
{
  ASSERT (indexes. size () == indexes_max. size ());

  FFOR (size_t, i, indexes. size ())
  {
    ASSERT (indexes [i] <= indexes_max [i]);
    if (indexes [i] < indexes_max [i])
    {
      indexes [i] ++;
      return true;
    }
    indexes [i] = 0;
  }

  return false;
}




// SubsetSearch

SubsetSearch::SubsetSearch (size_t whole_size,
                            size_t subset_size)
: subset (subset_size + 1)
{
  ASSERT (subset_size <= whole_size);
  FOR (size_t, i, subset_size)
    subset [i] = i;
  subset [subset_size] = whole_size;
}



bool SubsetSearch::next ()
{
  FFOR (size_t, i, subset. size () - 1)
    if (subset [i] + 1 == subset [i + 1])
      subset [i] = i;
    else
    {
      subset [i] ++;
      return true;
    }
  return false;
}



void SubsetSearch::complement (Vector<size_t> &disjoint) const
{
  ASSERT (! subset. empty ());
  ASSERT (subset. back () >= subset. size () - 1);
  const size_t size = subset. back () - (subset. size () - 1);
  disjoint. reserve (size);
  disjoint. clear ();
  size_t j = 0;
  for (size_t i : subset)
  {
    for (; j < i; j++)
      disjoint << j;
    j++;
  }
  ASSERT (disjoint. size () == size);
}




// MinSubsetSearch

bool MinSubsetSearch::next ()
{
  for (;;)
  {
    if (! subset. empty ())  // Previous invocation of next() has returned true
      if (! advance_wide ())
        return false;
    while (! complete ())
      if (! advance ())
        return false;
    // Check minimality
    size_t i = 0;
    while (i + 1 < subset. size ())
      if (complete (subset [i]))
        break;
      else
        i++;
    if (i + 1 == subset. size ())
      return true;
    ASSERT (! subset. empty ());
  }
}



bool MinSubsetSearch::advance ()
{
  if (verbose (-1))
  {
    cout << "advance:";
    saveText (cout);
    cout << endl;
  }
  if (! n)
    return false;
  if (subset. empty ())
    subset << 0;
  else if (subset. back () + 1 < n)
  {
    ASSERT (subset. size () < n);
    subset << subset. back () + 1;
  }
  else
    return advance_wide ();
  return true;
}



bool MinSubsetSearch::advance_wide ()
{
  size_t limit = n;
  while (! subset. empty ())
  {
    ASSERT (subset. back () + 1 <= limit);
    if (subset. back () + 1 == limit)
    {
      subset. pop ();
      limit--;
    }
    else
    {
      subset. back () ++;
      return true;
    }
  }
  return false;
}


namespace
{

struct ThisSubsetSearch : MinSubsetSearch
{
  Vector<Vector<bool>> arr;
  size_t whole {0};

  ThisSubsetSearch (size_t n_arg,
                    ulong seed)
    : MinSubsetSearch (n_arg)
    , arr (n)
    , whole (5)  // PAR
    { Rand rand (seed);
      for (Vector<bool>& vec : arr)
      {
        FOR (size_t, i, whole)
          vec << (rand. getProb () > 0.2);
        for (const bool b : vec)
          cout << b;
        cout << endl;
      }
    }

  bool complete (size_t exclude = no_index) const final
    {
      FOR (size_t, i, whole)
      {
        bool intersect = true;
        for (const size_t j : subset)
          if (j != exclude)
            if (! arr [j] [i])
            {
              intersect = false;
              break;
            }
        if (intersect)
          return false;
      }
      return true;
    }
};

}



void MinSubsetSearch::test (ulong seed)
{
  ThisSubsetSearch ss (5, seed);  // PAR
  ASSERT (! ss. complete ());
  while (ss. next ())
  {
    cout << "Found:";
    ss. saveText (cout);
    cout << endl;
  }
  ASSERT (ss. subset. empty ());
}




// Permute

Permute::Permute (size_t n)
: vec (n)
{
  FOR (size_t, i, n)
    vec [i] = i;
}



}  // namespace
