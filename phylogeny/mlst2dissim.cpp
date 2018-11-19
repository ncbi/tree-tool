// mlst2dissim.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../dm/numeric.hpp"
#include "../dm/dataset.hpp"
using namespace DM_sp;



namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Convert MLST data to a dissimilarity and print a " + dmSuff + "-file")
    {
  	  addPositional ("biosamples", "File with a list of BioSamples");
  	  addPositional ("pairs", "File with tab-delimited lines: <obj1> <obj2> <X1> <X2> <X3> <X4> <Result> <number of loci in common with same allele> <number of loci in common with different allele> ...");
  	  addFlag ("ani", "Compute density-filtered ANI-based Jukes-Cantor distance, otherwise MLST DNA Jaccard evolution distance");
  	  addFlag ("filtering", "Use density filtering?");
  	}



	void body () const
	{
		const string biosamplesFName = getArg  ("biosamples");
		const string pairsFName      = getArg ("pairs");
		const bool aniP              = getFlag ("ani");
		const bool filtering         = getFlag ("filtering");
		IMPLY (aniP, filtering);

		
		
    Set<string> objNames;
    {
      LineInput f (biosamplesFName);
      while (f. nextLine ())
      {
        trim (f. line);
        objNames << f. line;
      }
    }
    cerr << "# Objects: " << objNames. size () << endl;

    Dataset ds;
    for (const string& name : objNames)
      ds. appendObj (name);
    ds. setName2objNum ();
    
    auto attr = new PositiveAttr2 ("dist", ds, 6);  // PAR
    {
      LineInput f (pairsFName);
      while (f. nextLine ())
      {
        const string obj1 = findSplit (f. line, '\t');
        const string obj2 = findSplit (f. line, '\t');
        if (! objNames. contains (obj1))
          continue;
        if (! objNames. contains (obj2))
          continue;

/*          
Columns 1 and 2: genome pair
3: Number of loci after filtering common to the pair and with same allele
4: Number of loci after filtering common to the pair and with different allele
5: Total length of alleles in loci retained after filtering (shorter of two used for each locus)
6: Total length of alleles in loci removed during filtering
7: Number of differences in alleles in column 4 <- this is the column David asked me to add for getting at ANI estimate
8: Dummy ‘Result’ tag
9-13 are results before density filtering
9: number of loci in common with same allele
10: number of loci in common with different allele
11: number of loci seen only in the first genome
12: number of loci seen only in the second genome
13: number of differences in alleles in column 10
*/        

        const size_t sameAllelesFiltered = str2<size_t> (findSplit (f. line, '\t'));  // X3
        const size_t diffAllelesFiltered = str2<size_t> (findSplit (f. line, '\t'));  // X4
        const size_t lociDnaLen          = str2<size_t> (findSplit (f. line, '\t'));  // X5
        findSplit (f. line, '\t');  // X6
        const size_t diffDnaLen          = str2<size_t> (findSplit (f. line, '\t'));  // X7
        EXEC_ASSERT (findSplit (f. line, '\t') == "Result");                          // x8
        const size_t sameAlleles = str2<size_t> (findSplit (f. line, '\t'));          // x9
        const size_t diffAlleles = str2<size_t> (findSplit (f. line, '\t'));          // x10
        ASSERT (sameAllelesFiltered <= sameAlleles);
        ASSERT (diffAllelesFiltered <= diffAlleles);
        ASSERT (diffDnaLen <= lociDnaLen);

        Real dist = NaN;
        if (aniP)
        {
          const Prob ani = (Real) (lociDnaLen - diffDnaLen) / (Real) lociDnaLen;
          ASSERT (isProb (ani));
          dist = - log ((ani - 0.25) / 0.75);   // was: 100 *
        }
        else
        {
          const Prob cons = filtering
                              ? (Real) sameAllelesFiltered / (Real) (sameAllelesFiltered + diffAllelesFiltered)
                              : (Real) sameAlleles / (Real) (sameAlleles + diffAlleles);
          ASSERT (isProb (cons));
          dist = - log (cons);   // was: 10 *
        }
          
        const size_t row = ds. getName2objNum (obj1);
        const size_t col = ds. getName2objNum (obj2);
        ASSERT (row != NO_INDEX);
        ASSERT (col != NO_INDEX);
        if (! attr->isMissing2 (row, col))
          cerr << obj1 << ' ' << obj2 << ": " << "duplicate value" << endl;
        attr->putSymm (row, col, dist);
      }
    }    
    FOR (size_t, row, ds. objs. size ())
      attr->put (row, row, 0);

    
    ds. qc ();
    ds. saveText (cout);    
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);

}



