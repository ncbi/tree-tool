// snp2dm.cpp

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
*   Make the SNP-organism Data Master file and the organism-SNP Data Master file with the distance matrix on organisms
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/matrix.hpp"
#include "../dm/dataset.hpp"
using namespace DM_sp;
#include "../version.inc"



namespace 
{
	
	
struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Make the SNP-organism " + dmSuff + "-file and the organism-SNP "+ dmSuff + "-file with the distance matrix on organisms")
    {
  	  version = VERSION;
  	  addKey ("snps", "File", "File with selected SNPs (no file - all SNPs)");
  	  addPositional ("HGDP", "HGDP SNP file");
  	  addPositional ("dataset", "Output " + dmSuff + "-file with a distance matrix in KJukes");  
    }


	
	void body () const final
  {
	//const string snps    = args["snps"].AsString();
	  const string inFName = getArg ("HGDP");
	  const string dsFName = getArg ("dataset");
    QC_ASSERT (! dsFName. empty ());


  #if 0
    Set<string> snpSet;
    if (! snps. empty ())
    {
      ifstream is (snps);
      string s;
      while (! is. eof ())
      {
        is >> s;
        snpSet << s;
      }
      cout << "# Selected SNPs: " << snpSet. size () << endl;
    }
  #endif
    

    Dataset ds;  
  //Space1<ExtBoolAttr1> spSnp (ds, false);
        
    
    auto place = new NominAttr1 ("Place", ds);
    size_t allSnps = 0;
    size_t missings_max = 0;
    size_t redundantSnps = 0;
    size_t singleSnps = 0;
    {
      LineInput li (inFName, 10 * 1024 * 1024, 1000);  // PAR
      
      EXEC_ASSERT (li. nextLine ());
      {
        istringstream iss (li. line);
        string name;
        while (! iss. eof ())
        {
          iss >> name;  // individual
          if (name. empty ())
            break;
  		    ds. appendObj (name);
  		  //spSnp << new ExtBoolAttr1 (name, ds);
        }
      }
    //ASSERT (dsSnp. attrs. size () == dsOrg. objs. size ());
      
      // *place
      EXEC_ASSERT (li. nextLine ());
      {
        istringstream iss (li. line);
        string name;
        size_t objNum = 0;
        while (! iss. eof ())
        {
          iss >> name;  // place_name
          if (name. empty ())
            break;
  		    (*place) [objNum] = place->category2index (name);
  		    const_cast <Obj*> (ds. objs [objNum]) -> comment = name;
  		    objNum++;
        }
        ASSERT (objNum == ds. objs. size ());
      }
      
      while (li. nextLine ())
      {
        allSnps++;
        // SNP
        istringstream iss (li. line);
        size_t chromosome;
        string name;
        size_t position;
        char nuc1;
        string nuc2;
        //
        iss >> chromosome >> name >> position >> nuc1 >> nuc2;
        ASSERT (chromosome >= 1);
        ASSERT (chromosome <= 22);
        const string nucs ("ACGT");
        ASSERT (charInSet (nuc1, nucs));
        if (nuc2 != "-1")
        {
          ASSERT (nuc2. size () == 1);
          ASSERT (charInSet (nuc2 [0], nucs));
          ASSERT (nuc1 != nuc2 [0]);
        }
      //if (! snpSet. empty () && ! snpSet. contains (name))
        //continue;
        auto snp = new ExtBoolAttr1 (name, ds);
        size_t missings = 0;
        size_t count_1 = 0;
        size_t count_0 = 0;
        size_t objNum = 0;
        int value;
        while (! iss. eof ())
        {
          iss >> value;
          switch (value)
          {
            case -9: missings++; break;
            case  0: (*snp) [objNum] = EFALSE; count_0++; break;
            case  1: (*snp) [objNum] = ETRUE;  count_1++; break;
            default: ERROR_MSG ("Wrong value " + toString (value));
          }
          objNum++;
        }
        ASSERT (objNum == ds. objs. size ());
        maximize (missings_max, missings);
        if (   count_0 == 0
            || count_1 == 0
           )
          redundantSnps++;
        if (   count_0 == 1
            || count_1 == 1
           )
          singleSnps++;
      }
      ds. setName2objNum ();
      ds. qc ();
		}
		cout << "# All SNPs: " << allSnps << endl;
	//IMPLY (! snpSet. empty (), snpSet. size () == dsSnp. objs. size ());
		cout << "# Objects: " << ds. objs. size () << endl;
	//cout << "# Selected SNPs: " << dsSnp. objs. size () << endl;
		cout << "# Max. missings in a SNP: " << missings_max << endl;
		cout << "# Redundant SNPs: " << redundantSnps << endl;
		cout << "# Singleton SNPs: " << singleSnps << endl;
		

  #if 0
		auto flipped = new ExtBoolAttr1 ("Flipped", dsSnp);
		FOR (size_t, objNum, dsSnp. objs. size ())
		{
		  size_t n = 0;
		  size_t m = 0;
		  for (const ExtBoolAttr1* attr : spSnp)
		    if (! attr->isMissing (objNum))
		    {
		      n++;
		      if ((*attr) [objNum] == ETRUE)
		        m++;
		    }
		  (*flipped) [objNum] = toEbool (m > n / 2);
		  if ((*flipped) [objNum])
  		  for (const ExtBoolAttr1* attr : spSnp)
  		    if (! attr->isMissing (objNum))
  		      toggle ((* const_cast <ExtBoolAttr1*> (attr)) [objNum]);
		}
	#endif
		

    {
      OFStream os ("", dsFName, dmExt); 
      ds. saveText (os);
    }

  
  #if 0
	  auto dist = new PositiveAttr2 ("KJukes", dsOrg, 2);  // PAR
	  {
		  Progress prog (dsOrg. objs. size ());
		  FOR (size_t, row, dsOrg. objs. size ())
		  {
		  	prog ();
		  	dist->matr. putDiag (row, 0);
        ExtBoolAttr1* rowAttr = const_cast <ExtBoolAttr1*> (dsSnp. attrs. at (row) -> asExtBoolAttr1 ());
        ASSERT (rowAttr);
		    FOR_START (size_t, col, row + 1, dsOrg. objs. size ())
		    {
          ExtBoolAttr1* colAttr = const_cast <ExtBoolAttr1*> (dsSnp. attrs. at (col) -> asExtBoolAttr1 ());
          ASSERT (colAttr);
		      size_t n = 0;
		      size_t m = 0;
		      FOR (size_t, objNum, dsSnp. objs. size ())
		        if (   (*rowAttr) [objNum] != UBOOL
		            && (*colAttr) [objNum] != UBOOL
		           )
		        {
		          n++;
		          if ((*rowAttr) [objNum] == (*colAttr) [objNum])
		            m++;
		        }
		      const Prob p = (Real) m / (Real) n;
		      ASSERT (p > 0.25);  // p = 0.25 <=> evolutionally independent objects
		      const Real kJukes = - 100 * log (4.0/3.0 * (p - 0.25));
		      dist->matr. putSymmetric (row, col, kJukes);
		    }
		  }
		}
    dsOrg. setName2objNum ();
    dsOrg. qc ();
				
		
		{
      OFStream os ("", dsFName, dmExt); 
      dsOrg. saveText (os);
    }
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



