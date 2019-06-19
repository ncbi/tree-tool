// testDistTree.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "distTree.hpp"
using namespace DistTree_sp;



namespace 
{



struct ThisApplication : Application
{
	ThisApplication ()
	: Application ("Test a tree made by makeDistTree")
	{
    // Input
	  addKey ("input_tree", "Tree file");
	  addPositional ("data", dmSuff + "-file without " + strQuote (dmSuff) + " to read object comments");
	  addKey ("dissim_attr", "Dissimilarity attribute name in the <data> file");
	  addKey ("variance", "Dissimilarity variance: " + varianceTypeNames. toString (" | "), varianceTypeNames [varianceType]); 
	  addKey ("variance_power", "Power for -variance pow; > 0", "NaN");
	}



	void body () const final
  {
	  const string input_tree     = getArg ("input_tree");
	  const string dataFName      = getArg ("data");
	  const string dissimAttrName = getArg ("dissim_attr");
	               varianceType   = str2varianceType (getArg ("variance"));  // Global    
	               variancePower  = str2real (getArg ("variance_power"));    // Global

  //if (input_tree. empty ())
    //throw runtime_error ("-input_tree must be present");
		if (! isNan (variancePower) && varianceType != varianceType_pow)
		  throw runtime_error ("-variance_power requires -variance pow");
		if (isNan (variancePower) && varianceType == varianceType_pow)
		  throw runtime_error ("-variance_power is needed by -variance pow");
		if (variancePower <= 0.0)
		  throw runtime_error ("-variance_power must be positive");
    

    unique_ptr<DistTree> tree;
    tree. reset (isDirName (dataFName)
                   ? new DistTree (dataFName, input_tree, true, true, false)
                   : input_tree. empty ()
                     ? new DistTree (            dataFName, dissimAttrName, string())
                     : new DistTree (input_tree, dataFName, dissimAttrName, string())
                );
  //tree->multFixed = true;
  //tree->setDissimMult (false);  
    tree->qc ();
    
    
    // DistTree::optimizeReinsert()
    size_t improvements = 0;
    {
      Progress prog (tree->nodes. size ());
      for (const DiGraph::Node* node_ : tree->nodes)
      {
        prog ();
        const DTNode* from = static_cast <const DTNode*> (node_);
        if (   from->inDiscernible ()
            || from == tree->root
           )
          continue;
        Real nodeAbsCriterion_old = NaN;
        const NewLeaf nl (from, tree->name2leaf. size (), nodeAbsCriterion_old);
        ASSERT (nodeAbsCriterion_old >= 0.0);
        nl. qc ();
        const DTNode* to = nl. location. anchor;
        ASSERT (to);
        const Real improvement = nodeAbsCriterion_old - nl. location. absCriterion_leaf;
        if (improvement > 1e-6)  // PAR
          improvements++;
        else
          try 
          {
            QC_ASSERT (! negative (improvement));
            QC_ASSERT (from->getParent () == to->getParent () || from->getParent () == to);
            QC_ASSERT (fabs (from->len - nl. location. leafLen) < 1e-4);  // PAR
          }
          catch (const exception &e)
          {
            if (verbose ())
              cout         << from->getLcaName ()
                   << '\t' << from->len 
                   << '\t' << nl. location. leafLen 
                   << '\t' << (from->getParent () == to->getParent () || from->getParent () == to)
                   << '\t' << improvement
                   << '\t' << to->getLcaName ()
                   << '\t' << e. what ()
                   << endl;
            else
              throw;
          }
      } 
    }
    PRINT (improvements);
	}
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


