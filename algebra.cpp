// algebra.cpp

#undef NDEBUG
#include "common.inc"

#include "common.hpp"
using namespace Common_sp;
#include "expr.hpp"
using namespace Algebra_sp;



namespace 
{


void test (Expr* expr)
{
  ASSERT (expr);

  cout << endl;  
  expr->qc ();
  cout << expr->str () << endl;
  expr = expr->stndize ();
  cout << expr->str () << endl;
  
  delete expr;
}



struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Algebraic manipulations")
    {}



	void body () const
	{
	//const CArgs& args = GetArgs();
		
		{
	    Expr* expr = new OperExpr 
	                 ( Algebra_sp::times
	                 , new NatExpr (2)
	                 , new OperExpr ( Algebra_sp::plus
	                                , new OperExpr ( Algebra_sp::times
	                                	             , new VarExpr ("x")
	                                	             , new OperExpr (over, new NatExpr (2))
	                                	             )
	                                , new OperExpr (Algebra_sp::minus, new NatExpr (5))
	                                )
	                 );
	    test (expr);
	  }
    
    {
      Expr* expr = new OperExpr ( Algebra_sp::plus
      	                        , new OperExpr ( Algebra_sp::times
      	                        	             , new NatExpr (1)
      	                        	             , new OperExpr (over, new NatExpr (2))
      	                        	             )
      	                        , new OperExpr ( Algebra_sp::times
      	                        	             , new NatExpr (1)
      	                        	             , new OperExpr (over, new NatExpr (6))
      	                        	             )
      	                        );
      test (expr);      	                        	             
    }

    {
      Expr* expr = new OperExpr ( Algebra_sp::plus
      	                        , new NatExpr (1)
      	                        , new OperExpr (over, new NatExpr (2))
      	                        );
      test (expr);      	                        	             
    }

    {
      Expr* expr = new OperExpr ( Algebra_sp::plus
      	                        , new VarExpr ("x")
      	                        , new VarExpr ("x")
      	                        );
      test (expr);      	                        	             
    }

    {
      Expr* expr = new OperExpr ( Algebra_sp::plus
      	                        , new VarExpr ("x")
      	                        , new OperExpr ( Algebra_sp::plus
      	                        	             , new NatExpr (2)
      	                        	             , new VarExpr ("x")
      	                        	             )
      	                        );
      test (expr);      	                        	             
    }

    {
      Expr* expr = new OperExpr ( Algebra_sp::power
      	                        , new NatExpr (2)
      	                        , new OperExpr ( Algebra_sp::power
      	                        	             , new NatExpr (2)
      	                        	             , new NatExpr (3)
      	                        	             )
      	                        );
      test (expr);      	                        	             
    }

    {
      Expr* expr = new OperExpr ( Algebra_sp::power
      	                        , new OperExpr ( Algebra_sp::power
      	                        	             , new NatExpr (2)
      	                        	             , new NatExpr (2)
      	                        	             )
      	                        , new NatExpr (3)
      	                        );
      test (expr);      	                        	             
    }

    {
      Expr* expr = new OperExpr ( Algebra_sp::power
      	                        , new OperExpr ( Algebra_sp::power
      	                        	             , new NatExpr (2)
      	                        	             , new NatExpr (2)
      	                        	             )
      	                        , new OperExpr ( Algebra_sp::power
      	                        	             , new NatExpr (2)
      	                        	             , new NatExpr (3)
      	                        	             )
      	                        );
      test (expr);      	                        	             
    }

    {
      Expr* expr = new OperExpr ( Algebra_sp::power
      	                        , new NatExpr (2)
      	                        , new OperExpr ( Algebra_sp::power
								      	                       , new OperExpr ( Algebra_sp::power
								      	                       	              , new NatExpr (2)
								      	                       	              , new NatExpr (1)
								      	                        	            )
      	                        	             , new NatExpr (3)
      	                        	             )
      	                        );
      test (expr);      	                        	             
    }

    {
      // Not complete simplification ??
      Expr* expr = new OperExpr ( Algebra_sp::plus
                                , new OperExpr ( Algebra_sp::times
                                                , new OperExpr ( Algebra_sp::times
                                      	                       , new OperExpr ( Algebra_sp::plus
                                      	                                      , new VarExpr ("b")
                                      	                                      , new VarExpr ("n")
                                      	                                      )
                                      	                       , new OperExpr ( Algebra_sp::plus
                                      	                                      , new VarExpr ("b")
                                      	                                      , new OperExpr ( Algebra_sp::plus
                                                      	                                     , new VarExpr ("n")
                                                      	                                     , new NatExpr (1)
                                                      	                                     )
                                                      	                      )
                                      	                       )
                                      	        , new OperExpr ( Algebra_sp::plus
                                      	                       , new OperExpr ( Algebra_sp::plus
                                        	                                     , new OperExpr ( Algebra_sp::times
                                                        	                                    , new NatExpr (2)
                                                        	                                    , new VarExpr ("b")
                                                        	                                    )
                                                        	                     , new OperExpr ( Algebra_sp::times
                                                        	                                    , new NatExpr (2)
                                                        	                                    , new VarExpr ("n")
                                                        	                                    )
                                                                        	     )
                                                               , new NatExpr (1)
                                                               )
                                        	      )
                                , new OperExpr ( Algebra_sp::minus
                                               , new OperExpr ( Algebra_sp::times
                                                              , new OperExpr ( Algebra_sp::times
                                                    	                       , new VarExpr ("b")
                                                    	                       , new OperExpr ( Algebra_sp::plus
                                                    	                                      , new VarExpr ("b")
                                                    	                                      , new NatExpr (1)
                                                    	                                      )
                                                    	                       )
                                                    	        , new OperExpr ( Algebra_sp::plus
                                                    	                       , new OperExpr ( Algebra_sp::times
                                                      	                                    , new NatExpr (2)
                                                      	                                    , new VarExpr ("b")
                                                      	                                    )
                                                                             , new NatExpr (1)
                                                                             )
                                                      	      )
                                                )
                                );
                                        	      
      test (expr);      	                        	             
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
