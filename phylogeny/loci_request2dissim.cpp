// loci_request2dissim.cpp

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
*   Compute loci dissimilarities
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../version.inc"



namespace 
{


constexpr size_t locus_max = 30000;  // PAR
typedef array<size_t,locus_max> Loci;  



Loci getLoci (const string &fName)
{
  static_assert (sizeof (size_t) == 8, "sizeof (size_t) == 8");
  LineInput f (fName);
  Loci loci;  
  loci. fill (0);
  string locus;
  size_t n = 0;
  Istringstream iss;
  while (f. nextLine ())
  {
    iss. reset (f. line);
    size_t locusNum, allele;
    iss >> locus >> allele;
    ASSERT (iss. eof ());
    ASSERT (locus [0] == 'G');
    locusNum = str2<size_t> (locus. substr (1));
    ASSERT (locusNum);
    ASSERT (allele);
    if (locusNum > locus_max)
      throw runtime_error ("Too large locus number: " + toString (locusNum));
    loci [locusNum] = allele;
    n++;
  } 
  if (! n)
    throw runtime_error ("No loci for " + fName);
    
  return loci;
}
    



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Compute loci dissimilarities for requested pairs")
    {
      version = VERSION;
  	  addPositional ("pairs", "File with pairs of objects");
  	  addPositional ("loci_dir", "Directory with a file of loci for each object. File line format: G<locus number> <allele number>");
  	  addPositional ("out", "Output file");
  	//check min. number of loci ??
  	}



	void body () const final
	{
		const string pairsFName = getArg  ("pairs");
		const string loci_dir   = getArg  ("loci_dir");
		const string out        = getArg  ("out");
		ASSERT (! loci_dir. empty ());
		ASSERT (! out. empty ());
		
		
    OFStream output (out);
    ONumber on (output, 6, true);  // PAR
    PairFile input (pairsFName, false, false);
    while (input. next ())
    {
      const Loci loci1 (getLoci (loci_dir + "/" + input. name1));
      const Loci loci2 (getLoci (loci_dir + "/" + input. name2));
      size_t sameLoci = 0;
      size_t diffLoci = 0;
      FOR (size_t, i, locus_max)
        if (   loci1 [i]
            && loci2 [i]
           )
        {
          if (loci1 [i] == loci2 [i])
    	      sameLoci++;
    	    else
            diffLoci++;
        }
      ASSERT (sameLoci + diffLoci);
      const double dissim = sameLoci 
                              ? log ((double) (sameLoci + diffLoci) / (double) sameLoci)
                              : numeric_limits<double>::infinity ();        
      ASSERT (dissim >= 0);        
      output << input. name1 << '\t' << input. name2 << '\t' << dissim;
      if (verbose ())
        output << '\t' << sameLoci << '\t' << diffLoci;
      output << endl;
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



