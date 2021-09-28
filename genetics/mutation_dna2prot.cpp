// mutation_dna2prot.cpp

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
*   Print protein mutations given DNA mutations and CDS annotations
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "seq.hpp"
using namespace Seq_sp;
#include "../version.inc"



namespace 
{
  
  
struct Annot : Named
{
  size_t start {no_index};
  size_t stop {no_index};
  bool strand {true};
  
  
  Annot (const string &name_arg,
         size_t start_arg,
         size_t stop_arg)
    : Named (name_arg)
    , start (start_arg)
    , stop (stop_arg)    
    {
      QC_ASSERT (start != stop);
      if (start > stop)
      {
        swap (start, stop);
        strand = false;
      }
      QC_ASSERT (start);
      start--;
    }
  Annot (const Annot &other) = default;
  void qc () const override
    {
      if (! qc_on)
        return;
      Named::qc ();
      QC_ASSERT (start <= stop);
      QC_ASSERT (getLen () % 3 == 0);
    }
  void saveText (ostream &os) const override
    { os << name << ' ' << start + 1 << ".." << stop << ' ' << (strand ? '+' : '-'); }
    
    
  size_t getLen () const
    { return stop - start; }
  Dna getOrf (const Dna &ref) const
    { ASSERT (stop <= ref. seq. size ());
      Dna orf (name, ref. seq. substr (start, getLen ()), false);
      if (! strand)
        orf. reverse ();
      return orf;
    }
  bool overlap (const Mutation &mut) const
    { ASSERT (! mut. prot);
      return    mut. stop () > start
             && mut. pos     < stop;
    }
  bool contains (const Mutation &mut) const
    { ASSERT (! mut. prot);
      return    mut. pos     >= start
             && mut. stop () <= stop;
    }    
  size_t getOffset5 (const Mutation &mut) const
    { ASSERT (contains (mut));
      if (strand)
        return mut. pos - start;
      return stop - mut. stop ();
    }
  Mutation trimMutation (const Mutation &mut) const
    { ASSERT (overlap (mut));
      if (contains (mut))
        return mut;
      Mutation mut1 (mut);
      if (mut1. pos < start)
      { const size_t utr = start - mut1. pos;
        ASSERT (mut1. ref. size () > utr);
        mut1. ref. erase (0, utr);
        if (mut1. allele. size () > mut1. ref. size ())
          mut1. allele. erase (0, mut1. allele. size () - mut1. ref. size ());        
        mut1. pos = start;
      }
      if (mut1. stop () > stop)
      { ASSERT (mut1. pos < stop);
        const size_t cds = stop - mut1. pos;
        ASSERT (mut1. ref. size () > cds);
        mut1. ref. erase (cds);
        if (mut1. allele. size () > cds)
          mut1. allele. erase (cds);
      }
      ASSERT (mut1. allele. size () <= mut1. ref. size ());
      ASSERT (contains (mut1));
      return mut1;
    }
  Annot trimAnnot (const Mutation &mut) const
    { ASSERT (contains (mut));
      Annot an (*this);
      if (mut. ref. size () == mut. allele. size ())
        return an;
      if (mut. allele. size () > mut. ref. size ())
        an. stop += mut. allele. size () - mut. ref. size ();
      else
        an. stop -= mut. ref. size () - mut. allele. size ();
      return an;
    }
};




struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print protein mutations given DNA mutations and CDS annotations")
  	{
      version = VERSION;
  	  addPositional ("mut", "Input list of DNA mutations");
  	  addPositional ("ref", "Reference DNA sequence file");
  	  addPositional ("annot", "File with CDS annotations in <ref>, line format: <CDS name> <start> <stop>; 1-based, <stop> is the last nucleotide position");
  	  addKey ("gencode", "NCBI genetic code", "1");
  	  addFlag ("dna_mut", "print DNA mutation");
  	}



	void body () const final
  {
		const string  mutFName    = getArg ("mut");
		const string  refFName    = getArg ("ref");
		const string  annotFName  = getArg ("annot");
		const Gencode gencode     = (Gencode) arg2uint ("gencode");
		const bool    printDnaMut = getFlag ("dna_mut");


    unique_ptr<const Dna> ref;
    {
      LineInput f (refFName);
      EXEC_ASSERT (f. nextLine ());
      ref. reset (new Dna (f, 100000, false));
      QC_ASSERT (! ref->seq. empty ());
    }

    
    Vector<Annot> annots;
    {
      LineInput f (annotFName);
      Istringstream iss;
      while (f. nextLine ())
      {
        iss. reset (f. line);
        string name;
        size_t start = 0;
        size_t stop = 0; 
        iss >> name >> start >> stop;
        QC_ASSERT (start || stop);
        const Annot annot (name, start, stop);
        annot. qc ();
        QC_ASSERT (annot. getLen ());
        QC_ASSERT (annot. stop <= ref->seq. size ());
        annots << move (annot);
      }
    }


    // 2 DNA mutations in the same codon produce 2 protein mutations: for a'b'c' a'bc -> aa1 and abc' -> aa2 ??
    {
      LineInput f (mutFName);
      while (f. nextLine ())
      {
        const Mutation dnaMut (false, f. line);
        dnaMut. qc ();
        QC_ASSERT (dnaMut. geneName. empty ());
        QC_ASSERT (dnaMut. stop () <= ref->seq. size ());
        for (const Annot& annot : annots)
        {
          if (! annot. overlap (dnaMut))
            continue;
          if (verbose ())
            cout << annot << endl
                 << dnaMut << endl;
          const Mutation dnaMutTrimmed (annot. trimMutation (dnaMut));
          if (dnaMutTrimmed. ref == dnaMutTrimmed. allele)
            continue;
          dnaMutTrimmed. qc ();
          ASSERT (annot. contains (dnaMutTrimmed));
          Dna mutDna (*ref);
          dnaMutTrimmed. replace (mutDna);
          size_t translationStart = 0;
          const Peptide refProt (annot. getOrf (*ref). makePeptide (1/*frame*/, gencode, false, false, translationStart));
          refProt. qc ();
          ASSERT (! refProt. seq. empty ());
          Mutation protMut;
          Mutation protMut_prev;
          if (dnaMutTrimmed. isFrameshift ())
          {
            const size_t prot_start = annot. getOffset5 (dnaMutTrimmed) / 3;
            const string prot_ref (refProt. seq. substr (prot_start, dna2codons_len (dnaMutTrimmed. ref. size ())));
            if (dnaMutTrimmed. ambig && ! prot_ref. empty ())
              protMut = Mutation (annot. name, prot_start, prot_ref, string (prot_ref. size (), 'X'), false); 
            else 
              protMut = Mutation (annot. name, prot_start, prot_ref, string (), true); 
          }
          else
          {
            const Annot mutAnnot (annot. trimAnnot (dnaMutTrimmed));
            mutAnnot. qc ();
            const Peptide mutProt (mutAnnot. getOrf (mutDna). makePeptide (1/*frame*/, gencode, false, false, translationStart));
            if (refProt. seq == mutProt. seq)
              continue;
            if (mutProt. seq. empty ())
              protMut = Mutation (annot. name, 0, refProt. seq, string (), false);
            else
            {
              if (verbose ())
              {
                PRINT (refProt);
                PRINT (mutProt);
              }
              mutProt. qc ();
              size_t pos = 0;
              while (refProt. seq [pos] == mutProt. seq [pos])
                pos++;
              size_t refEnd = refProt. seq. size ();
              size_t mutEnd = mutProt. seq. size ();
              while (   refEnd > pos
                     && mutEnd > pos
                     && refProt. seq [refEnd - 1] == mutProt. seq [mutEnd - 1]
                    )
              {
                ASSERT (refEnd);
                ASSERT (mutEnd);
                refEnd--;
                mutEnd--;
              }
              ASSERT (pos <= refEnd);
              ASSERT (pos <= mutEnd);
              protMut = Mutation ( annot. name
                                 , pos
                                 , refProt. seq. substr (pos, refEnd - pos)
                                 , mutProt. seq. substr (pos, mutEnd - pos)
                                 , false
                                 );
            }
          }
          protMut. qc ();
          ASSERT (protMut. prot); 
          if (! (protMut == protMut_prev))
          {
            if (printDnaMut)
              cout << dnaMut << '\t';
            cout << protMut << endl;
          }
          protMut_prev = protMut;
        }
      }
    }
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);  
}



