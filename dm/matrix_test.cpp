// matrix_test.cpp

#undef NDEBUG
#include "../common.inc"

#include "../common.hpp"
using namespace Common_sp;
#include "matrix.hpp"
using namespace DM_sp;




namespace
{

struct ThisApplication : Application
{
  ThisApplication ()
    : Application ("Test matrix.cpp")
    {
      addPositional ("go", "Go");
      addFlag ("eigens", "test Eigens");
    }

	

	void body () const final
	{
    const bool eigensP = getFlag ("eigens");
    

	  Matrix pd (4);
	  {
		  Matrix m (false, pd, false);
		  m. putAll (0);
		  m. put (false, 0, 0, 5);
		  m. put (false, 0, 1, -3);
		  m. put (false, 0, 2, 7);
		  m. put (false, 0, 3, -2);
		  m. put (false, 1, 0, 0);
		  m. put (false, 1, 1, 2);
		  m. put (false, 1, 2, 8);
		  m. put (false, 1, 3, 10);
		  m. put (false, 2, 0, -7);
		  m. put (false, 2, 1, -5);
		  m. put (false, 2, 2, 9);
		  m. put (false, 2, 3, 2);
		  m. put (false, 3, 0, 4);
		  m. put (false, 3, 1, -9);
		  m. put (false, 3, 2, 2);
		  m. put (false, 3, 3, 1);
		  m. qc ();
		  pd. multiply (false, m, false, m, true);
		  pd. qc ();
  	  if (verbose ())
		  {
		  	cout << endl << "PD:" << endl;
	  	  pd. print (cout);
	  	}
		}
	  
    {
		  const Matrix ch (pd. getCholesky (false));
		  ASSERT (ch. defined ());
		  ch. qc ();
		  if (verbose ())
		  {
		  	cout << endl << "Choleski:" << endl;
		    ch. print (cout);
		  }
	  
		  Matrix m1 (false, pd, false);
		  m1. multiply (false, ch, false, ch, true);
		  m1. qc ();
		  if (verbose ())
		  {
		  	cout << endl << "Back to PD:" << endl;
	  	  m1. print (cout);
	  	}
		  ASSERT (pd. maxAbsDiff (false, m1, false) < 1e-3);
		}
			  
	  {
		  if (verbose ())
		  	cout << endl << "Eigens1:" << endl;
	    ASSERT (pd. psd);
		  const Eigens eigens (pd, 4, 1, 0, 1e-6, 1000/*, true*/);
		  eigens. qc ();
		  Matrix back (false, pd, false);
		  eigens. restore (back);
		  back. qc ();
		  if (verbose ())
		  {
		  	cout << "Back form eigens to PD: " << endl;
		    back. print (cout);
		  }
		  ASSERT (pd. maxAbsDiff (false, back, false) < 1e-3);
		}

    {
      if (verbose ())
        cout << endl << "eigens2" << endl;
  	  Matrix pd2 (4);  // Rank = 2
  	  {
  		  Matrix m (false, 2, 4);
  		  m. putAll (0);
  		  m. put (false, 0, 0, 5);
  		  m. put (false, 0, 1, -3);
  		  m. put (false, 0, 2, 7);
  		  m. put (false, 0, 3, -2);
  		  m. put (false, 1, 0, 0);
  		  m. put (false, 1, 1, 2);
  		  m. put (false, 1, 2, 8);
  		  m. put (false, 1, 3, 10);
  		  pd2. multiply (false, m, true, m, false);
  		}
  		ASSERT (pd2. psd);
  		pd2. qc ();
  		Eigens eigens (pd2, 4, 1, 0, 1e-6, 1000/*, true*/);
  		eigens. qc ();
      ASSERT (eigens. getDim () == 2);
  	  if (verbose ())
  	  {
  			eigens. print (cout);
  			cout << endl;
  			eigens. basis. print (cout);
  		}
    }
		
		{
      if (verbose ())
        cout << endl << "eigens3" << endl;
  	  Matrix m (3);  
		  m. put (false, 0, 0, 5);
		  m. put (false, 0, 1, -3);
		  m. put (false, 0, 2, 7);
		  m. put (false, 1, 0, 0);
		  m. put (false, 1, 1, 2);
		  m. put (false, 1, 2, 8);
		  m. put (false, 2, 0, 3);
		  m. put (false, 2, 1, -2);
		  m. put (false, 2, 2, 3);
      Real maxCorrection;
      size_t row_bad;
      size_t col_bad;	  
		  m. symmetrize (maxCorrection, row_bad, col_bad);
		  if (verbose ())
		    m. print (cout);
		  ASSERT (! m. psd);
		  m. qc ();
  		Eigens eigens (m, 3, 1, 0, 1e-6, 1000/*, false*/);
  		eigens. qc ();
  	  if (verbose ())
  	  {
  			eigens. print (cout);
  			cout << endl;
  			eigens. basis. print (cout);
  		}
      ASSERT (eigens. getDim () == 3);
      ASSERT_EQ (eigens. totalExplainedFrac (), 1, 1e-4);
      ASSERT_EQ (eigens. values. sumSqr (), m. sumSqr (), 1e-4);
		}
		
		{
			// http://en.wikipedia.org/wiki/Fisher's_exact_test
			Matrix mat (2);
			mat. put (false, 0, 0,  1);
			mat. put (false, 0, 1,  9);
			mat. put (false, 1, 0, 11);
			mat. put (false, 1, 1,  3);
			mat. qc ();
			const Prob pValue = exp (mat. getLnFisherExact (true));
			ASSERT_EQ (pValue, 0.00137973, 1e-8);
		}

		{
			// J.G.Upton, Analysis of Cross-Tabluated Data, table 2.3
			Matrix mat (2);
			mat. put (false, 1, 0, 10);
			mat. put (false, 1, 1, 20);
			mat. put (false, 0, 0,  5);
			mat. put (false, 0, 1, 25);
			mat. qc ();
			const Prob pValue = exp (mat. getLnFisherExact (true));
			ASSERT_EQ (pValue, 0.1163, 1e-4);
			const Real chi2 = mat. getChi2 ();
			ASSERT_EQ (chi2, 2.2222, 1e-4);
		}

		{
			Matrix mat (false, 5, 3);
			mat. put (false, 0, 0, 31);
			mat. put (false, 0, 1, 14);
			mat. put (false, 0, 2, 52);
			mat. put (false, 1, 0,  2);
			mat. put (false, 1, 1,  1);
			mat. put (false, 1, 2, 12);
			mat. put (false, 2, 0,  6 );
			mat. put (false, 2, 1,  8);
			mat. put (false, 2, 2,  9);
			mat. put (false, 3, 0, 16);
			mat. put (false, 3, 1, 20);
			mat. put (false, 3, 2, 24);
			mat. put (false, 4, 0, 19);
			mat. put (false, 4, 1, 17);
			mat. put (false, 4, 2,  2);
			mat. qc ();
			const Real chi2 = mat. getChi2 ();
			ASSERT_EQ (chi2, 40.9748, 1e-4);
		}
		
		
		// MATLAB: ~1000 times faster (with 8 cores)
		{
			ifstream f ("data/masten_60.mat");
			if (! f. good ())
			  throw runtime_error ("File data/masten_60.mat does not exist");
		  Matrix mat (false, f, true);
		  mat. psd = true;
		  mat. qc ();
		 
		  { 
			  const Matrix res (mat. getSqrt ());
	      if (false)
			  {
				  OFStream out ("res.mat");
				  res. saveText (out);
				}
				
	      Matrix test (mat. rowsSize ());
	      test. multiply (false, res, true, res, false);
	      ASSERT (mat. maxAbsDiff (false, test, false) < 1e-5);  // PAR
	      if (false)
	      {
				  OFStream out ("test.mat");
				  test. saveText (out);      	
	      }
	    }

      if (Chronometer::enabled)
      {
	      const Chronometer_OnePass cop ("sqrt(matrix)");  
	      FOR (size_t, i, 100)
				  mat. getSqrt ();
		  }
		  

		  if (eigensP)
		  {
			  Eigens ei (mat, mat. rowsSize (), 1, 0, 1e-5, 10000);   // PAR
			  ei. values. print (cout);
			  cout << endl;
		  			  	
		  	// Choleski iterations
	      const Chronometer_OnePass cop ("Choleski iterations");  
	      FOR (size_t, tries, 100)
	      {
	        Matrix prod (mat);
	        size_t iter = 0;
	        for (;;)
	        {
	        	iter++;
				  	const Matrix ch (prod. getCholesky (false));
				  	ASSERT (ch. defined ());
		        prod. multiply (false, ch, true, ch, false);
		        // Sum of off-diagonal elements
		        Real s = 0;
		        FOR (size_t, i, prod. rowsSize (false))
			        FOR (size_t, j, prod. rowsSize (false))
			          if (i != j)
			          	s += fabs (prod. get (false, i, j));
			      if (s / sqr ((Real) prod. rowsSize (false)) < 1e-6)  // PAR   
			      	break;
			    }
			    cout << iter << endl;  // ??
			    Vector<Real> vec;  vec. reserve (prod. rowsSize (false)); 
	        FOR (size_t, i, prod. rowsSize (false))
	          vec << prod. getDiag (i);
	        vec. sort ();
	        if (false)
 	        {
		        for (const Real x : vec)
		          cout << ' ' << x;
		        cout << endl;
		      }
	      }
		  }
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


