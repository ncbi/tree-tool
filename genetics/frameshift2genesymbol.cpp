// frameshift2genesymbol.cpp

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
*   Convert a raw frameshift symbol to gene symbol
*
*/


#undef NDEBUG

#include "../common.hpp"
using namespace Common_sp;
#include "seq.hpp"
using namespace Seq_sp;
#include "../version.inc"

#include "../common.inc"



namespace 
{
  
  
struct FS final : Root
{
  // Input
  string contig;
  string prot;
  // 0-based
  size_t protPos;    // aa
  size_t contigPos;  // bp
  Strand strand;
  // Output
  char ref {'\0'};
  char allele {'\0'};
  size_t termLen {no_index};
    // Distance from contigPos to a stop codon
    // allele = '*' => 0
    // no_index <=> stop codon is not found


  explicit FS (const string &line)
    {
      istringstream iss (line);
      string fs;
      iss >> prot >> contig >> fs;
      QC_ASSERT (isLeft (fs, "fs_"));
      ASSERT (! prot. empty ());
      ASSERT (! contig. empty ());
      
      // strand
      const string strandS = rfindSplit (fs, '_');
      if (strandS == "0")
        strand = -1;
      else if (strandS == "1")
        strand = 1;
      else
        throw runtime_error ("Unknown strand: " + strQuote (strandS));
          
      contigPos = str2<size_t> (rfindSplit (fs, '_'));
      protPos   = str2<size_t> (rfindSplit (fs, '_'));
      
      QC_ASSERT (fs == "fs");
    }  
  void saveText (ostream &os) const final
    {
      os         << contig 
         << '\t' << prot
         << '\t' << ref << protPos + 1 << allele << "fs";
      if (termLen != no_index)
        os << "Ter" << termLen;
      os << '\t' << "fs_" << protPos << '_' << contigPos << '_' << (strand == 1 ? 1 : 0);
    }
    
    
  char contig2aa (const Dna &dna,
                  size_t offset,
                  Gencode gencode) const
    // Input: offset: from contigPos
    // Return: '\0' <=> offset is outside dna
    {
      QC_ASSERT (contigPos + 3 <= dna. seq. size ());
      
      if (strand == 1)
      {
        const size_t i = contigPos + offset * 3;
        if (i + 3 > dna. seq. size ())
          return '\0';
        return codon2aa (& dna. seq [i], gencode, false);
      }
      
      ASSERT (strand == -1);
      if (contigPos < (offset + 1) * 3)
        return '\0';
      const size_t i = contigPos - (offset + 1) * 3;
      if (i + 3 > dna. seq. size ())
        return '\0';
      string s (dna. seq. substr (i, 3));
      reverseDna (s);
      return codon2aa (s. c_str (), gencode, false);
    }
};

	
	
struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Convert raw frame shift symbols to standard gene symbols according to https://hgvs-nomenclature.org/stable/recommendations/protein/frameshift/")
	  {
      version = VERSION;
		  addPositional ("nucl", "Input nucleotide FASTA file");
		  addPositional ("prot", "Input protein FASTA file");
		  addPositional ("fs_tab", "Table with lines: <protein identifier in <prot>>  <contig identifier in <nucl>>  <raw frameshift symbol: fs_<qpos>_<spos>_<sstrand(1/0)>>");
		  addKey ("gencode", "NCBI genetic code for translated BLAST", "11");
	  }

  
  
	void body () const final
  {
	  const string nuclFName = getArg ("nucl");
	  const string protFName = getArg ("prot");
	  const string fsFName   = getArg ("fs_tab");
    const Gencode gencode  = (Gencode) arg2uint ("gencode"); 

	  
	  Vector<FS> fss;
	  {
	    LineInput f (fsFName);
	    while (f. nextLine ())
	      fss << std::move (FS (f. line));
	  }
	  if (fss. empty ())
	    return;
	  
	  // FS::{allele,termLen}
		{
		  Multifasta fa (nuclFName, false);
		  while (fa. next ())
		  {
	      const Dna dna (fa, 100000/*PAR*/, true);
		    dna. qc ();	
		    const string id (dna. getId ());
		    for (FS& fs : fss)
		      if (fs. contig == id)
		      {
  		      fs. allele = fs. contig2aa (dna, 0, gencode);
  		      QC_ASSERT (fs. allele);
  		      fs. termLen = 0;
  		      for (size_t offset = 0; ; offset++)
  		      {
  		        const char aa = fs. contig2aa (dna, offset, gencode);
  		        if (aa == '\0')
  		        {
  		          fs. termLen = no_index;
  		          break;
  		        }
  		        if (aa == '*')
  		          break;
  		        fs. termLen++;
  		      }
		        break;
  		    }
		  }
		}

    // FS::ref
		{
		  Multifasta fa (protFName, true);
		  while (fa. next ())
		  {
		    const Peptide pep (fa, 1000/*PAR*/, true);
		    pep. qc ();	
		    const string id (pep. getId ());
		    for (FS& fs : fss)
		      if (fs. prot == id)
		      {
		        QC_ASSERT (fs. protPos < pep. seq. size ());
		        fs. ref = pep. seq [fs. protPos];
		        break;
		      }
		  }
		}
		
		// fss
		for (const FS& fs : fss)
		{
		  fs. saveText (cout);
		  cout << '\n';
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



