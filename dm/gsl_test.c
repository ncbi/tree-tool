#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>



int main ()
{
  const size_t size = 60;  // PAR

  gsl_matrix* mat_init = gsl_matrix_alloc (size, size);
  assert (mat_init);
  {
	  FILE* f = fopen ("data/masten_60.txt", "r");  // PAR
	  assert (f);
	  size_t i, j;
	  for (i = 0; i < size; i++)
		  for (j = 0; j < size; j++)
		  {
			  double x;
		  	assert (fscanf (f, "%lf", & x) == 1);
		  	gsl_matrix_set (mat_init, i, j, x);
		  }
		fclose (f);
  }
  

  gsl_matrix* res = gsl_matrix_alloc (size, size);
  assert (res);

  size_t iter;
  const clock_t start = clock ();
  for (iter = 0; iter < 10000; iter++)  // PAR
  {
	  // sqrtm
	  size_t i, j;

	  // To be changed
	  gsl_matrix* mat = gsl_matrix_alloc (size, size);
	  assert (mat);
	  for (i = 0; i < size; i++)
		  for (j = 0; j < size; j++)
		    gsl_matrix_set (mat, i, j, gsl_matrix_get (mat_init, i, j));
	
	  gsl_eigen_symmv_workspace* ws = gsl_eigen_symmv_alloc (size);
	  assert (ws);	  
	  gsl_vector* eVec = gsl_vector_alloc (size);
	  assert (eVec);
	  gsl_matrix* eMat = gsl_matrix_alloc (size, size);
	  assert (eMat);	  
	  assert (gsl_eigen_symmv (mat, eVec, eMat, ws) == 0);
	  gsl_matrix_free (mat);
	  gsl_eigen_symmv_free (ws);
	  
	  for (i = 0; i < size; i++)
	  {
	  	const double x = gsl_vector_get (eVec, i);
	  	assert (x >= 0);
	  	gsl_vector_set (eVec, i, sqrt (x));
	  }  
	  for (i = 0; i < size; i++)
		  for (j = 0; j < size; j++)
		  {
		  	double s = 0;
		  	size_t k;
			  for (k = 0; k < size; k++)
		  	  s += gsl_matrix_get (eMat, i, k) * gsl_vector_get (eVec, k) * gsl_matrix_get (eMat, j, k);
		  	gsl_matrix_set (res, i, j, s);
		  }

    if (iter == 0)
    {
    	double diff2_sum = 0;
		  for (i = 0; i < size; i++)
			  for (j = 0; j < size; j++)
			  {
			  	double s = 0;
			  	size_t k;
				  for (k = 0; k < size; k++)
			  	  s += gsl_matrix_get (res, i, k) * gsl_matrix_get (res, k, j);
			  	const double diff = s - gsl_matrix_get (mat_init, i, j);
			  	diff2_sum += diff * diff;
				}
			printf ("diff2_sum = %g\n", diff2_sum);
		}	


	  gsl_matrix_free (eMat);
	  gsl_vector_free (eVec);
  }
  printf ("%lf sec.\n", (double) (clock () - start) / CLOCKS_PER_SEC);


  gsl_matrix_free (res);
  gsl_matrix_free (mat_init);


  return 0;
}
