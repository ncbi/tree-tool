// attr2_2attr1.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "dataset.hpp"
using namespace DM_sp;



namespace 
{


struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Print a " + dmSuff + "-file where all two-way attributes are converted into one-way attributes.\n\
Objects will be named: <objName1>-<objName2>")
    {
  	  addPositional ("file", dmSuff + "-file");
  	  addFlag ("symmetrize", "Make all two-way attributes symmetric by averaging");
  	  addFlag ("zero_diagonal", "Skip diagonal values because they are 0");
  	}



	void body () const final
	{
		const string fName       = getArg ("file");
		const bool symmetrize    = getFlag ("symmetrize");
		const bool zero_diagonal = getFlag ("zero_diagonal");
		
		
    const Dataset ds (fName);
    
    Dataset ds1;
    FFOR (size_t, i, ds. objs. size ())
      FOR (size_t, j, symmetrize ? i + 1 : ds. objs. size ())
        if (! zero_diagonal || i != j)
        {
          const size_t objNum = ds1. appendObj (ds. pair2name (i, j, symmetrize));
          var_cast (ds1. objs [objNum]) -> mult =   ds. objs [i] -> mult
                                                  * ds. objs [j] -> mult;
        }
    ds1. setName2objNum ();
    
    for (const Attr* attr : ds. attrs)
      if (const Attr2* attr2 = attr->asAttr2 ())
      {
        Attr1* attr1 = attr2->createAttr1 (ds1);
        if (symmetrize)
          var_cast (attr2) -> symmetrize ();
        FFOR (size_t, i, ds. objs. size ())
          FFOR (size_t, j, symmetrize ? i + 1 : ds. objs. size ())
            if (! zero_diagonal || i != j)
              {
                const size_t objNum = ds1. getName2objNum (ds. pair2name (i, j, symmetrize));
                ASSERT (objNum != NO_INDEX);
                attr1->str2value (objNum, attr2->value2str (i, j));
              }
      }
      
    ds1. qc ();
    
    ds1. saveText (cout);    
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



