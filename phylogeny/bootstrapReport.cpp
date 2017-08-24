// bootstrapReport.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
//#include "distTree.hpp"
//using namespace DistTree_sp;



namespace 
{
  
  
struct ThisApplication : Application
{
	ThisApplication ()
	: Application ("Analyze bootstrap result")
	{
		// Input
	//addPositional ("input_tree", "Tree");
	  addPositional ("bootstrap", "Bootstrap result: file with lines: <N> match[+-]: <node LCA-name>");
	  addKey ("replicas", "# replicas (match+ + match-) - for QC", "0");
	  addKey ("support_min", "Min. match+ / (match+ + match-) to report a node", "0.333");   // =1.0/3.0
	  addFlag ("print_nodes", "Print support for each interior arc identified by a child node");
	}



	void body () const
  {
  //const string input_tree        = getArg ("input_tree");
		const string bootstrap         = getArg ("bootstrap");
		const size_t replicas_expected = str2<size_t> (getArg ("replicas"));
		const double support_min       = str2<double> (getArg ("support_min"));
		const bool print_nodes         = getFlag ("print_nodes");
		ASSERT (replicas_expected);
		ASSERT (support_min >= 0);
    

  #if 0
    DistTree tree (input_tree, string (), string (), false);  
    if (verbose ())
      tree. qc ();     
  #endif
      
    typedef  Vector<size_t>  Support;  // size() = 2
    map <string/*node name*/, Support>  node2support;
    {
      LineInput f (bootstrap);
      while (f. nextLine ())
      {
        istringstream oss (f. line);
        size_t n;
        string match, name;
        oss >> n >> match >> name;
        ASSERT (n);
        ASSERT (match == "match+" || match == "match-");
        ASSERT (! name. empty ());
        const bool matchP = (match == "match+");
      #if 1
        Support& support = node2support [name];
        if (support. empty ())
          support. resize (2, 0);
        ASSERT (! support [matchP]);
        support [matchP] = n;
      #else
        const Steiner* node = tree. lcaName2node (name) -> asSteiner ();
        ASSERT (node);
        ASSERT (! node->bootstrap [matchP]);
        const_cast <Steiner*> (node) -> bootstrap [matchP] = n;
      #endif
      }
    }
    
    
  #if 1 
    double supportSum = 0;
    size_t supported = 0;
    double varSum = 0;
    size_t n = 0;
    for (const auto it : node2support)
    {
 	    const size_t replicas = it. second [0] + it. second [1];
 	    ASSERT (replicas);
 	    if (replicas != replicas_expected)
 	    {
 	      cout << it. first << ' ' << replicas << endl;
 	      ERROR;
 	    }
 	    // Probability of a Bernoulli random variable
 	    const double support = (double) it. second [1] / (double) replicas;
 	    if (support > 0.5)
 	      supported++;
 	    varSum += support * (1 - support);
 	    n++;
 	    if (support < support_min)
 	      continue;
 	    supportSum += support;
 	    if (print_nodes)
        cout        << it. first 
             << ' ' << it. second [0] 
             << ' ' << it. second [1] 
             << ' ' << support
             << endl;
    }
    ASSERT (supported <= n);
    ONumber on (cout, 3, false);
    cout << "Sum: " << supportSum << endl;
    cout << "Var: " << varSum / (double) n << endl;
    cout << "Supported fraction: " << (double) supported / (double) n << endl;
  #else
   	for (const DiGraph::Node* node : tree. nodes)  
   	  if (const Steiner* st = static_cast <const DTNode*> (node) -> asSteiner ())
   	  {
   	    const size_t replicas = st->bootstrap [0] + st->bootstrap [1];
   	    if (replicas >= replicas_min)
   	      cout << st->getLcaName () << ' ' << st->bootstrap [0] << ' ' << st->bootstrap [1] << endl;
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


