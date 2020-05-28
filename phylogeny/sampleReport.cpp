// sampleReport.cpp

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
*   Analyze tree sampling result
*
*/


#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "../version.inc"



namespace 
{
  
  
struct ThisApplication : Application
{
	ThisApplication ()
	: Application ("Analyze tree sampling result")
	{
	  version = VERSION;
		// Input
	  addPositional ("sampling", "Tree sampling result: file with lines: <N> match[+-]: <node LCA-name>");
	  addKey ("replicas", "# replicas (match+ + match-) - for QC", "0");
	//addKey ("support_min", "Min. match+ / (match+ + match-) to report a node", "0.3333333");   // =1.0/3.0
	  addKey ("stable_min", "Min. match+ / (match+ + match-) for a node to be stable", "0.92");   // PAR
	  addFlag ("print_nodes", "Print support for each interior arc identified by a child node");
	}



	void body () const final
  {
		const string sampling          = getArg ("sampling");
		const size_t replicas_expected = str2<size_t> (getArg ("replicas"));
  //const double support_min       = str2<double> (getArg ("support_min"));
		const double stable_min        = str2<double> (getArg ("stable_min"));
		const bool print_nodes         = getFlag ("print_nodes");
		ASSERT (replicas_expected);
	//ASSERT (support_min >= 0);
	//ASSERT (support_min <= 1);
		ASSERT (stable_min >= 0);
		ASSERT (stable_min <= 1);
    

    typedef  Vector<size_t>  Support;  // size() = 2
    map <string/*node name*/, Support>  node2support;
    {
      LineInput f (sampling);
      Istringstream iss;
      while (f. nextLine ())
      {
        iss. reset (f. line);
        size_t n;
        string match, name;
        iss >> n >> match >> name;
        ASSERT (n);
        ASSERT (match == "match+" || match == "match-");
        ASSERT (! name. empty ());
        const bool matchP = (match == "match+");
        Support& support = node2support [name];
        if (support. empty ())
          support. resize (2, 0);
        ASSERT (! support [matchP]);
        support [matchP] = n;
      }
    }
    
    
    double supportSum = 0;
    size_t stable = 0;
    double varSum = 0;
    size_t n = 0;
    for (const auto& it : node2support)
    {
 	    const size_t replicas = it. second [0] + it. second [1];
 	    ASSERT (replicas);
 	    if (replicas != replicas_expected)
 	      throw runtime_error (it. first + ": # replicas = " +  toString (replicas) + ", # expected = "  + toString (replicas_expected));
 	    // Probability of a Bernoulli random variable
 	    const double support = (double) it. second [1] / (double) replicas;
 	    if (support >= stable_min)  
 	      stable++;
 	    varSum += support * (1 - support);
 	    n++;
 	    if (support < 1.0/3.0 /*support_min*/)
 	      continue;
 	    supportSum += support;
 	    if (print_nodes)
        cout        << it. first 
             << ' ' << it. second [0] 
             << ' ' << it. second [1] 
             << ' ' << support
             << endl;
    }
    ASSERT (stable <= n);
    const ONumber on (cout, 3, false);
    cout << "RF sum: " << supportSum << endl;
    cout << "Var: " << varSum / (double) n << endl;
    cout << stable_min * 100 << "%-Stable: " << stable << endl;
	}
};



}  // namespace



int main (int argc, 
          const char* argv[])
{
  ThisApplication app;
  return app. run (argc, argv);
}


