// compareTrees.cpp 

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "distTree.hpp"
using namespace DistTree_sp;
#include "featureTree.hpp"
using namespace FeatureTree_sp;



namespace 
{
  
  
string frequencyS;

  

map <const Tree::TreeNode* /*isInteriorType()*/, string/*lcaName*/>  node2name;
  // For all Tree's
  

void tree2names (const Tree &tree)
// Update: node2name
{
 	for (const DiGraph::Node* node_ : tree. nodes)  
 	{
 	  const Tree::TreeNode* node = static_cast <const Tree::TreeNode*> (node_);
    if (node->isInteriorType ())
      node2name [node] = node->getLcaName ();
  }
}


string getNodeName (const Tree::TreeNode* node)
// Input: node2name
{
  ASSERT (node);
  const string s (node2name [node]);
  ASSERT (! s. empty ());
  return s;
}

  
  
typedef  StringVector  Leaves;
  // searchSorted



Leaves tree2leaves (const Tree &tree,
                    bool frequentOnly)
{
  Leaves leaves;  leaves. reserve (tree. nodes. size ());
 	for (const DiGraph::Node* node : tree. nodes)  
 	{
 	  const Tree::TreeNode* tn = static_cast <const Tree::TreeNode*> (node);
    if (tn->isLeafType ())
      if (! frequentOnly || (    frequencyS == "none"
                             || (frequencyS ==   "directed" && tn->frequentChild)
                             || (frequencyS == "undirected" && tn->frequentDegree == 1)
                            )
         )
        leaves << node->getName ();
  }
  leaves. sort ();
  ASSERT (leaves. isUniq ());
  return leaves;
}

  
  
map <const Tree::TreeNode*, Leaves> node2leaves;
  // For all Tree's



void setNode2leaves (const Tree::TreeNode* node)
// Output: node2leaves
{
  ASSERT (node);
  ASSERT (node2leaves [node]. empty ());
  if (node->isLeafType ())
    node2leaves [node] << node->getName ();
  else 
  	for (const DiGraph::Arc* arc : node->arcs [false])
  	{
  	  const Tree::TreeNode* child = static_cast <const Tree::TreeNode*> (arc->node [false]);
  	  setNode2leaves (child);
  	  node2leaves [node] << node2leaves [child];
  	}
  node2leaves [node]. sort ();
  ASSERT (node2leaves [node]. isUniq ());
}



inline bool even (size_t x)
{
  return x % 2 == 0;
}
  
  
  
void adjustNode2leaves (const Tree &tree,
                        const Leaves &allLeaves)
// Update: node2leaves
{
  ASSERT (allLeaves. searchSorted);
  
 	const size_t all = allLeaves. size ();
  for (auto& it : node2leaves)
  {
    const Tree::TreeNode* node = it. first;
    ASSERT (node->graph);
    if (& node->getTree () != & tree)
      continue;
    Leaves& leaves = it. second;
  	const size_t size = leaves. size ();
  	ASSERT (size <= all);
  	if (   size > all / 2
  	    || (even (all) && size == all / 2 && ! leaves. containsFast (allLeaves. front ()))
  	   )
  	{
  	  Leaves s (allLeaves);
      ASSERT (s. searchSorted);  
      ASSERT (leaves. searchSorted);  
  	  s. setMinus (leaves);
  	  leaves = s;
      ASSERT (leaves. searchSorted);  
  	}
  	ASSERT (leaves. size () <= all / 2);
  	IMPLY (even (all) && leaves. size () == all / 2, leaves. containsFast (allLeaves. front ()));
  }
}



struct Signature
{
  size_t leaves {0};
  string front;
  string back;
  
  Signature () = default;
  Signature (size_t leaves_arg,
             const string &front_arg,
             const string &back_arg)
    : leaves (leaves_arg)
    , front (front_arg)
    , back (back_arg)
    { ASSERT (leaves);
      ASSERT (! front. empty ());
      ASSERT (! back.  empty ()); 
    }    
};



bool operator< (const Signature &a,
                const Signature &b)
{ 
  LESS_PART (a, b, leaves);
  LESS_PART (a, b, front);
  LESS_PART (a, b, back);
  return false;
}



Signature leaves2signature (const Leaves &leaves)
{
  ASSERT (leaves. searchSorted);
  return leaves. empty () 
           ? Signature ()
           : Signature (leaves. size (), leaves. front (), leaves. back ());
}




//

struct ThisApplication : Application
{
	ThisApplication ()
	: Application ("Remove " + real2str (DistTree_sp::rareProb, 2) + "-infrequent leaves from tree 1 depending on -frequency, compare two trees by Robinson-Foulds method, print matching and mismatching interior nodes for Tree 1.\n\
A node with an empty set of leaves matches.")
	{
		// Input
	  addPositional ("input_tree1", "Tree 1");
	  addPositional ("input_tree2", "Tree 2");
	  addKey ("type", "Tree type: dist|feature", "dist");
	  addKey ("frequency", "Node frequency is computed for directed|undirected tree; 'none' - not used", "none");
	    // "local"|"global" ??
	//addFlag ("skip_empty", "Skip arcs with empty sets of leaves");
	  addFlag ("arc_info", "Print arc information: length, depth, # leaves");
	}



	void body () const final
  {
		const string input_tree1 = getArg ("input_tree1");
		const string input_tree2 = getArg ("input_tree2");
		const string treeType    = getArg ("type");
		             frequencyS  = getArg ("frequency");
	//const bool   skip_empty  = getFlag ("skip_empty");
	  const bool   arc_info    = getFlag ("arc_info");
		             
		if (! (   treeType == "dist" 
		       || treeType == "feature" 
		      )
		   )
		  throw runtime_error ("Wrong tree type");
		if (! (   frequencyS == "none"       
		       || frequencyS == "directed"
		       || frequencyS == "undirected"  
		      )
		   )
		  throw runtime_error ("Wrong frequency");
		       
		       
    Common_sp::AutoPtr<Tree> tree1;
    if (treeType == "dist")
      tree1 = new DistTree (input_tree1, string (), string (), false);
    else
      tree1 = new FeatureTree (input_tree1, string (), string (), false, true, false);
    if (verbose ())
      tree1->qc ();     
      
    Common_sp::AutoPtr<Tree> tree2;
    if (treeType == "dist")
      tree2 = new DistTree (input_tree2, string (), string (), false);
    else
      tree2 = new FeatureTree (input_tree2, string (), string (), false, true, false);
    if (verbose ())
      tree2->qc (); 
      
    tree2names (*tree1);
    tree2names (*tree2);


    VectorPtr<Tree::TreeNode> interiorArcNodes1;  interiorArcNodes1. reserve (tree1->nodes. size ());
      // Does not depend on *tree2
   	for (const DiGraph::Node* node_ : tree1->nodes)  
   	{
   	  const Tree::TreeNode* node = static_cast <const Tree::TreeNode*> (node_);
   	  if (node == tree1->root)
   	    continue;
      if (! node->isInteriorType ())
        continue;
      const Tree::TreeNode* parent = node->getParent ();
      ASSERT (parent);
      if (parent == tree1->root)
      {
        const VectorPtr<DiGraph::Node> children (parent->getChildren ());
        ASSERT (children. size () >= 2);
        if (children. size () == 2)  // root is transient
        {
          VectorPtr<Tree::TreeNode> interiors;  interiors. reserve (2);
          for (const DiGraph::Node* child : children)
          {
            const Tree::TreeNode* treeChild = static_cast <const Tree::TreeNode*> (child);
            if (treeChild->isInteriorType ())
              interiors << treeChild;
          }
          if (interiors. size () < 2)
            continue;
          ASSERT (interiors. size () == 2);
          ASSERT (interiors. contains (node));
          if (interiors [0] == node)
            continue;
        }
      }
   	  interiorArcNodes1 << node;
   	}

    
    {  
      tree1->setFrequentChild  (DistTree_sp::rareProb); 
      tree1->setFrequentDegree (DistTree_sp::rareProb); 

      const Leaves leaves1 (tree2leaves (*tree1, true));
      const Leaves leaves2 (tree2leaves (*tree2, false));
      cout << "# Leaves in " << input_tree1 << ": " << leaves1. size () << endl;
      cout << "# Leaves in " << input_tree2 << ": " << leaves2. size () << endl;
      cout << endl;
      
      cout << "Deleting from " << input_tree1 << endl;
      const size_t n1 =   tree1->restrictLeaves (leaves1, false)   // Keep only frequent (stable) leaves
                        + tree1->restrictLeaves (leaves2, false);
      cout << "# Deleted: " << n1 << endl;
      cout << endl;
      
      cout << "Deleting from " << input_tree2 << endl;
      const size_t n2 = tree2->restrictLeaves (leaves1, true);
      cout << "# Deleted: " << n2 << endl;
      cout << endl;
      
    //tree1->qc ();
      tree2->qc ();
    }
    {
    	const auto pred = [] (const DiGraph::Node* n) { return static_cast <const Tree::TreeNode*> (n) -> isLeafType (); };
      const auto leaf_num_1 = Common_sp::count_if (tree1->nodes, pred);
      const auto leaf_num_2 = Common_sp::count_if (tree2->nodes, pred);
      ASSERT (leaf_num_1 == leaf_num_2);
    }
        
    setNode2leaves (tree1->root);
    setNode2leaves (tree2->root); 
    {
      const Leaves allLeaves (tree2leaves (*tree1, false));  // Same for *tree2
      adjustNode2leaves (*tree1, allLeaves);   
      adjustNode2leaves (*tree2, allLeaves);   
    }

    map <Signature, VectorPtr<Tree::TreeNode>> signature2nodes2;  // tree2
   	for (const DiGraph::Node* node2_ : tree2->nodes)  
   	{
   	  const Tree::TreeNode* node2 = static_cast <const Tree::TreeNode*> (node2_);
   	  signature2nodes2 [leaves2signature (node2leaves [node2])] << node2;
   	}


    map <const Tree::TreeNode* /*tree1*/, const Tree::TreeNode* /*tree2*/> node2node;  
      // key: !nullptr
      // value: not deterministic <= Leaf's in tree2 are removed 
   	for (const Tree::TreeNode* node1 : interiorArcNodes1)  
   	{
   		const Signature sig (leaves2signature (node2leaves [node1]));
      if (sig. leaves)
      {
	   	  const VectorPtr<Tree::TreeNode> nodes2 (signature2nodes2 [sig]);
	   	  for (const Tree::TreeNode* node2 : nodes2)
	     	  if (node2leaves [node1] == node2leaves [node2])
	   	    {
	 	        node2node [node1] = node2;
	   	      break;  
	   	    }
	   	}
	   	else /*if (! skip_empty)*/
	   		node2node [node1] = nullptr;
   	}

   	
   	MeanVar mv;
   	for (const Tree::TreeNode* node1 : interiorArcNodes1)  
   	{
      cout << "match" << (contains (node2node, node1) ? '+' : '-') 
      	   << '\t' << getNodeName (node1);
      if (arc_info)
      	cout 
      	   << '\t' << node2leaves [node1]. size ()
      	   << '\t' << node1->getParentDistance ()
      	   << '\t' << node1->getRootDistance ();
      cout << endl;
      mv << node1->getRootDistance ();
    }
    cout << endl << "Depth" << '\t';
    mv. saveText (cout);
    cout << endl;
        
      
    // ??
    // Correlation between arc lengths of node2node keys and values
    // Distribution of arc lengths of non-matching node2node keys
    // Distribution of arc lengths of non-matching node2node values
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


