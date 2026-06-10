// assignment.cpp

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
*   Assignment problem
*
*/


#undef NDEBUG

#include "common.hpp"
#include "graph.hpp"
using namespace Common_sp;
#include "version.inc"

#include "common.inc"



namespace 
{
  

// https://en.wikipedia.org/wiki/Hungarian_algorithm

template <typename T, typename U>
  using Pair = std::pair<T, U>;
template <typename T>
  using Vector = std::vector<T>;

template <typename T>
  using NumericLimits = std::numeric_limits<T>;

/**
 * @brief Checks if b < a
 *
 * Sets a = min(a, b)
 * @param a The first parameter to check
 * @param b The second parameter to check
 * @tparam The type to perform the check on
 * @return true if b < a
 */
template <typename T> 
  constexpr bool ckmin(T& a, const T& b) { 
    return b < a ? a = b, true : false; 
}

/**
 * @brief Performs the Hungarian algorithm.
 *
 * Given J jobs and W workers (J <= W), computes the minimum cost to assign each
 * prefix of jobs to distinct workers.
 *
 * @tparam T a type large enough to represent integers on the order of J *
 * max(|C|)
 * @param C a matrix of dimensions JxW such that C[j][w] = cost to assign j-th
 * job to w-th worker (possibly negative)
 *
 * @return a vector of length J, with the j-th entry equaling the minimum cost
 * to assign the first (j+1) jobs to distinct workers
 */
template <typename T> 
Vector<T> hungarian(const Vector<Vector<T>>& C) {
    const size_t/*int*/ J = /*static_cast<int>*/(C.size());
    const size_t/*int*/ W = /*static_cast<int>*/(C[0].size());
    assert(J <= W);
    // job[w] = job assigned to w-th worker, or -1 if no job assigned
    // note: a W-th worker was added for convenience
    Vector<size_t/*int*/> job(W + 1, no_index /*-1*/);
    Vector<T> ys(J); 
    Vector<T> yt(W + 1);  // potentials
    // -yt[W] will equal the sum of all deltas
    Vector<T> answers;
    const T inf = NumericLimits<T>::max();
    for (size_t/*int*/ jCur = 0; jCur < J; ++jCur) 
    {  
        // assign jCur-th job
        size_t/*int*/ wCur = W;
        job[wCur] = jCur;
        // min reduced cost over edges from Z to worker w
        Vector<T> minTo(W + 1, inf);
        Vector<size_t/*int*/> prev(W + 1, no_index/*-1*/);  // previous worker on alternating path
        Vector<bool> inZ(W + 1);     // whether worker is in Z
        while (job[wCur] != no_index /*-1*/) {    // runs at most jCur + 1 times
            inZ[wCur] = true;
            const size_t/*int*/ j = job[wCur];
            T delta = inf;
            size_t/*int*/ wNext = 0;
            for (size_t/*int*/ w = 0; w < W; ++w) {
                if (!inZ[w]) {
                    if (ckmin(minTo[w], C[j][w] - ys[j] - yt[w]))
                        prev[w] = wCur;
                    if (ckmin(delta, minTo[w])) 
                        wNext = w;
                }
            }
            // delta will always be nonnegative,
            // except possibly during the first time this loop runs
            // if any entries of C[jCur] are negative
            for (size_t/*int*/ w = 0; w <= W; ++w) {
                if (inZ[w]) {
                    ys[job[w]] += delta;
                    yt[w] -= delta;
                } else {
                    minTo[w] -= delta;
                }
            }
            wCur = wNext;
        }
        // update assignments along alternating path
        for (size_t /*int*/ w; wCur != W; wCur = w) 
            job[wCur] = job[w = prev[wCur]];
        answers.push_back(-yt[W]);
    }
    return answers;
}

#if 0
/**
 * @brief Performs a sanity check for the Hungarian algorithm.
 *
 * Sanity check: https://en.wikipedia.org/wiki/Hungarian_algorithm#Example
 * First job (5):
 *   clean bathroom: Bob -> 5
 * First + second jobs (9):
 *   clean bathroom: Bob -> 5
 *   sweep floors: Alice -> 4
 * First + second + third jobs (15):
 *   clean bathroom: Alice -> 8
 *   sweep floors: Carol -> 4
 *   wash windows: Bob -> 3
 */
void sanityCheckHungarian() {
    Vector<Vector<int>> costs{{8, 5, 9}, {4, 2, 4}, {7, 3, 8}};
    assert((hungarian(costs) == Vector<int>{5, 9, 15}));
    cerr << "Sanity check passed." << endl;
}
#endif


  
struct ThisApplication final : Application
{
  ThisApplication ()
    : Application ("Assignment problem. Print minimized or maximized total cost")
    {
      version = VERSION;
      addPositional ("in", "File with bipartite arcs: <object1>\\t<object2>\\t<cost>");
      addFlag ("max", "Maximize total cost");
    }


  struct AssignmentNode final : DiGraph::Node
  {
    const string name;
    bool matched {false};
    
    AssignmentNode (DiGraph &graph_arg,
                    const string &name_arg)
      : DiGraph::Node (graph_arg)
      , name (name_arg)
      {}
  };


  struct CostArc final : DiGraph::Arc
  {
    double cost {NaN};
    
    CostArc (AssignmentNode* start,
             AssignmentNode* end,
             double cost_arg)
      : DiGraph::Arc (start, end)
      , cost (cost_arg)
      {}    
      
    bool operator< (const CostArc& other) const
      { return cost < other. cost; }
    bool available () const
      { for (const bool out : {false, true})
          if (static_cast <const AssignmentNode*> (node [out]) -> matched)
            return false;
        return true;
      }
    void use ()
      { for (const bool out : {false, true})
          var_cast (static_cast <const AssignmentNode*> (node [out])) -> matched = true;
      }
      // Postcondition: !available()
  };



	void body () const final
	{
	  const string fName    = getArg ("in");
	  const double maxCoeff = getFlag ("max") ? -1.0 : 1.0;
	  
	  
	  // Naive algorithm ??
	  DiGraph gr;  // Bipartite: from obj1 to obj2
	  {
  	  LineInput f (fName);
  	  map<string,AssignmentNode*> obj1_2node;
  	  map<string,AssignmentNode*> obj2_2node;
  	  while (f. nextLine ())
  	  {
  	    string obj1 (findSplit (f. line, '\t'));
  	    string obj2 (findSplit (f. line, '\t'));
  	    const double cost = maxCoeff * str2<double> (f. line);
  	    QC_ASSERT (! isNan (cost));
  	    trim (obj1);
  	    QC_ASSERT (! obj1. empty ());
  	    trim (obj2);
  	    QC_ASSERT (! obj2. empty ());
  	    AssignmentNode* &n1 = obj1_2node [obj1];
  	    if (! n1)
  	      n1 = new AssignmentNode (gr, obj1);
  	    AssignmentNode* &n2 = obj2_2node [obj2];
  	    if (! n2)
  	      n2 = new AssignmentNode (gr, obj2);
  	    if (const DiGraph::Arc* arc = n1->incident (n2, true))
  	      var_cast (static_cast <const CostArc*> (arc)) -> cost += cost;
  	    else
  	      new CostArc (n1, n2, cost);
  	  }
  	}
  	
  	VectorPtr<CostArc> arcs;  arcs. reserve (obj1_2node. size ());
  	for (const auto& it : obj1_2node)
  	  for (const DiGraph::Arc* arc : it. second->arcs [true])
  	  {
  	    ASSERT (arc);
  	    arcs << static_cast <const CostArc*> (arc);
  	  }  	  
  	arcs. sortPtr ();  	
  	  
    double total = 0.0;
    for (const CostArc* ca : arcs)
    {
      ASSERT (ca);
      if (ca->available ())
      {
        total += ca->cost;
        var_cast (ca) -> use ();
      }
      ASSERT (! ca->available ());
    }
    cout << maxCoeff * total << endl;  
  }  
};



}  // namespace




int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}



