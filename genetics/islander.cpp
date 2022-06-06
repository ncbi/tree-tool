// islander.cpp

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
*   Find annotation islands
*
*/
   
   
#undef NDEBUG 
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../version.inc"



// ThisApplication


namespace
{  
  
struct Island
{
  // 0-based
  size_t start;
  size_t stop;
  // start < stop
};  



void process (const string &accession,
              size_t len,
              size_t gap_max,
              Vector<Island> &islands,
              ostream* nonIslandCoor)
{
  IMPLY (accession. empty (), islands. empty ());
  if (accession. empty ())
    return;
  
  size_t islands_len = 0;  
  size_t dropped = 0;
  size_t dropped_len = 0;
  if (! islands. empty ())
  {
    const Island* prev = nullptr;
    for (const Island& island : islands)
    {
      ASSERT (island. start <= island. stop);
      IMPLY (prev, prev->stop < island. start);
      islands_len += island. stop - island. start;
      prev = & island;
      if (verbose ())
        cout << accession << '\t' << island. start << '\t' << island. stop << endl;
    }
    Island& front = islands. front ();
    if (front. start && front. start <= gap_max)
    {
      dropped++;
      dropped_len += front. start;
      front. start = 0;
    }
    Island& back = islands. back ();
    ASSERT (back. stop <= len);
    const size_t tail = len - back. stop;
    if (tail && tail <= gap_max)
    {
      dropped++;
      dropped_len += tail;
      back. stop = len;
    }
  }
  
  if (nonIslandCoor)
  {
    size_t start = 0;
    for (const Island& island : islands)
    {
      if (island. start > start)
        *nonIslandCoor << accession << '\t' << start << '\t' << island. start << endl;
      start = island. stop;
    }
    if (len > start)
      *nonIslandCoor << accession << '\t' << start << '\t' << len << endl;
  }

  cout         << accession 
       << '\t' << islands_len
       << '\t' << islands. size ()
       << '\t' << len
       << '\t' << dropped_len
       << '\t' << dropped
       << '\t' << (double) islands_len / (double) len * 100.0
       << endl;
  
  islands. clear ();
}
  
}



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Find annotation islands, print summary statistics")
    {
      version = VERSION;
      // Input
      addPositional ("annot", "Annotation tsv-file with columns: plasmid\tcontig_len\tstart\tstop (1-based), sorted by plasmid (accession), start, stop desc"); 
      addPositional ("gap_max", "Max. gap to merge annotations into an island");
      addKey ("non_island_coor", "Coordinates of non-islands");
    }



  void body () const final
  {
    const string annotFName =               getArg ("annot");
    const size_t gap_max    = str2<size_t> (getArg ("gap_max"));    
    const string non_island_coor = getArg ("non_island_coor");
    

    const TextTable annot (annotFName);
    annot. qc ();
    QC_ASSERT (annot. header. size () >= 4);
    QC_ASSERT (annot. header [1]. numeric);
    QC_ASSERT (annot. header [2]. numeric);
    QC_ASSERT (annot. header [3]. numeric);
    
    unique_ptr<OFStream> non_island_coorF;
    if (! non_island_coor. empty ())
      non_island_coorF. reset (new OFStream (non_island_coor));
    
    cout         << "#plasmid"
         << '\t' << "non_core_length"
         << '\t' << "non_core_num"
         << '\t' << "contig_len"
         << '\t' << "dropped_len"
         << '\t' << "dropped_num"
         << '\t' << "percent_non_core_length"
         << endl;
    string accession_prev;
    size_t len_prev = 0;
    size_t start_prev = 0;
    size_t stop_prev = 0;
    Vector<Island> islands;
    for (const StringVector& row : annot. rows)
      try
      {      
        const string& accession = row [0];
        QC_ASSERT (! accession. empty ());
        QC_ASSERT (accession_prev <= accession);
        const size_t len = (size_t) stol (row [1]);
        QC_ASSERT (len);
        QC_IMPLY (accession_prev == accession, len_prev == len);
        QC_ASSERT (row [2]. empty () == row [3]. empty ());
        if (row [2]. empty ())
        {
          QC_ASSERT (accession_prev < accession);
          process (accession_prev, len_prev, gap_max, islands, non_island_coorF. get ());
          start_prev = 0;
          stop_prev = 0;
        }
        else
        {
                size_t start = (size_t) stol (row [2]);
          const size_t stop  = (size_t) stol (row [3]);
          QC_ASSERT (start >= 1);      
          QC_ASSERT (start <= stop);
          QC_ASSERT (stop <= len);
          start--;
          QC_IMPLY (accession_prev == accession, start_prev <= start);
          QC_IMPLY (accession_prev == accession && start_prev == start, stop_prev >= stop);
          bool newIsland = true;
          if (accession_prev == accession)
          {
            if (! islands. empty ())
            {
              Island& island = islands. back (); 
              if ((int) start - (int) island. stop <= (int) gap_max)
              {
                maximize (island. stop, stop);
                newIsland = false;
              }
            }
          }
          else
            process (accession_prev, len_prev, gap_max, islands, non_island_coorF. get ());
          if (newIsland)
            islands << Island {start, stop};
          start_prev = start;
          stop_prev = stop;
        }
        accession_prev = accession;
        len_prev = len;
      }
      catch (...)
      {
        save (cout, row, ',');
        throw;
      }
    process (accession_prev, len_prev, gap_max, islands, non_island_coorF. get ());
  }
};




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



