// matrix.hpp
// C++, ISO/EIC 14882:1998(E), hosted implementation


#ifndef MATRIX_HPP_59376
#define MATRIX_HPP_59376



#include "../common.hpp"
#include "numeric.hpp"



namespace DM_sp
{
using namespace std;
using namespace Common_sp;

// In scanning operations NAN-elements are skipped
// "bool t" means whether the operation is applied to the transposed matrix, i.e. rows <-> columns
// <Function>Row means <Function> applied to a row
// Object is a row and column of a square matrix



struct Matrix;



Real multiplyVec (const Matrix &m1,
                  bool         t1,
                  size_t       row1,
                  const Matrix &m2,
                  bool         t2,
                  size_t       row2);
  // Multiplying row1 of m1 and row2 of m2
  // Requires: m1.rowsSize(!t1) == m2.rowsSize(!t2)

Real sumSqrDifferenceVec (const Matrix &m1,
                          bool         t1,
                          size_t       row1,
                          const Matrix &m2,
                          bool         t2,
                          size_t       row2);
  // Return: Squared distance in metric L_2
  // Requires: m1.rowsSize(!t1) == m2.rowsSize(!t2)

Real sumAbsDifferenceVec (const Matrix &m1,
                          bool         t1,
                          size_t       row1,
                          const Matrix &m2,
                          bool         t2,
                          size_t       row2);
  // Return: Distance in metric L_1
  // Requires: m1.rowsSize(!t1) = m2.rowsSize(!t2)

Real maxAbsDifferenceVec (const Matrix &m1,
                          bool         t1,
                          size_t       row1,
                          const Matrix &m2,
                          bool         t2,
                          size_t       row2);
  // Return: Distance in metric L_inf
  // Requires: m1.rowsSize(!t1) = m2.rowsSize(!t2)

Real getCovarianceVec (const Matrix &m1,
                       bool         t1,
                       size_t       row1,
                       const Matrix &m2,
                       bool         t2,
                       size_t       row2,
                       bool         Biased);
  // Covariance of two vectors, corrected
  // Requires: m1.rowsSize(!t1) = m2.rowsSize(!t2)

Real getCorrelationVec (const Matrix &m1,
                        bool         t1,
                        size_t       row1,
                        const Matrix &m2,
                        bool         t2,
                        size_t       row2);
  // Correlation of two vectors
  // Requires: m1.rowsSize(!t1) = m2.rowsSize(!t2)


void distribution2MeanVar (const Matrix &distribution,
                           size_t       DistributionRow,
                           bool         distributionT,
                           const Matrix &Values,
                           size_t       ValuesRow,
                           bool         ValuesT,
                           Real         &mean, 
                           Real         &Variance);
  // Requires: distribution.rowsSize(!t) = Values.rowsSize(!t)
  //           distribution.definedRow ()
  //           distribution must be normalized



struct Determinant : LogReal
{
  size_t size;
  
  explicit Determinant (size_t size_arg)
    : size (size_arg)
    {}
  Determinant (size_t size_arg,
               Real r)
    : LogReal (r)
    , size (size_arg)
    {}
  void qc () const;
  
  Real getAverage () const
    { return exp (n / (Real) size) * (sign ? 1 : -1); }
};



struct Eigen;  



struct Matrix : Root
// Matrix of Real
{
protected:
  Vector <Vector <Real> /*columns*/> rows;
  size_t colsSize;
public:  
  bool psd;

 
  Matrix (bool t,
          size_t maxRow,
          size_t maxCol) 
    { init (t, maxRow, maxCol); }          
  Matrix (bool t,
          const Matrix &source,
          bool sourceT) 
    { init (t, source. rowsSize (sourceT), source. rowsSize (! sourceT)); }
    // Rows are not copied
  Matrix () 
    { init (false, 0, 0); }
  Matrix (const Matrix &source) 
    : colsSize (0)
    { init (false, source. rowsSize (false), source. rowsSize (true));
      copyDataCheck (false, source, false);
    }
  void operator= (const Matrix &source)
    { copyData (source); }
//Matrix (Matrix&&) noexcept = default;
//Matrix& operator= (Matrix&&) = default;
  Matrix* copy () const override
    { return new Matrix (*this); }
  Matrix (bool t,
          istream &f,
          bool square);
    // Invokes: init(), loadSize(), loadData()
    // See formatInstruction()
private:
  void init (bool t,
             size_t maxRow,
             size_t maxCol);
    // rows[][] = 0
    // psd = isSquare()
  void loadSize (istream &f,
                 bool square,
                 size_t &maxRow,
                 size_t &maxCol);
    // Output: maxRow, maxCol
  void loadData (bool t,
                 istream &f);
public:
  void load (bool t,
             istream &f,
             bool square);
    // Invokes: loadSize(), loadData()
  void reserve (bool t,
                size_t maxRow,
                size_t maxCol);
	void resize (bool t,
               size_t maxRow,
               size_t maxCol)
    { // delData ();
    	init (t, maxRow, maxCol);
    }
  static string formatInstruction (const string &fName) 
	  { return   "Format of <" + fName + ">:\n" 
	           + "  <# Rows> [<# Cols>]\n" 
	           + "  <X_0_0> <X_0_1> ... <X_0_Cols>\n"
	           + "  <X_1_0> <X_1_1> ... <X_1_Cols>\n" 
	           + "  ...\n"
	           + "  <X_Rows_0> <X_Rows_1> ... <X_Rows_Cols>\n"
	           + "If <# Cols> is omitted then # Cols = # Rows\n";
	  }
private:
  void delData ();
public:
  void qc () const override;
  void saveText (ostream &os) const override
    { saveFile_ (false, false, os); }


protected:
  static void swapRowCol (bool t,
                          size_t &row, 
                          size_t &col) 
    { if (t) std::swap (row, col); }
public:
 
  // Output
  void printRow (bool  t,
                 size_t row,
                 ostream &f) const;
    // "%f". No '\n'
  void saveFile_ (bool t,
                  bool square,
                  ostream &f) const;

  // Basic access
  Real get (bool t,
            size_t row,
            size_t col) const
    { swapRowCol (t, row, col);
      return rows [row] [col]; 
    }
  bool get (bool t,
            size_t row, 
            size_t col,
            Real &a) const
    { a = get (t, row, col); return ! isNan (a); }
  Real getDiag (size_t row) const
    { return get (false, row, row); }
  bool getDiag (size_t row,
                Real &a) const
    { return get (false, row, row, a); }
  void put (bool t,
            size_t row,
            size_t col,
            Real a)
    { swapRowCol (t, row, col);
      rows [row] [col] = a;
      psd = false;
    }
  void putDiag (size_t row,
                Real a)
    { Keep<bool> kp (psd);
      put (false, row, row, a); 
    }
  void putSymmetric (size_t row, 
                     size_t col, 
                     Real a)
    { Keep<bool> kp (psd);
      put (false, row, col, a);
      put (false, col, row, a);
    }
      
  // Change rowsSize()
  void insertRows (bool t,
		               size_t firstRow,
		               size_t rowInc);
  void deleteRows (bool t,
                   size_t firstRow,
                   size_t rowDec);
  void trimRows (bool t,
                 size_t newMaxRow);
    // Requires: newMaxRow <= rowsSize(t)
  void moveDeleteRow (bool t,
                      size_t sourceRow,
                      size_t destRow)
    { if (sourceRow != destRow)
        copyRow (t, destRow, *this, t, sourceRow);
      deleteRows (t, sourceRow, 1);
    }
  void deleteObject (size_t objNum);
    // Requires: isSquare()

  // processRow2Col
  typedef Real (Matrix::*ProcessRowFunc) (bool /*t*/,
                                          size_t /*row*/) const;
  Matrix* processRow2Col (bool t,
                          ProcessRowFunc processRow) const;
    // Return: !nullptr
    //         Result->rowsSize(t)  == rowsSize(t)
    //         Result->rowsSize(!t) == 1
    // To avoid warning invoke processRow as & Matrix::processRow

  // Basic info
  size_t rowsSize (bool t) const  
    { return t ? colsSize : rows. size (); }
  size_t minSize () const
    { return std::min (rowsSize (false), rowsSize (true)); }
  bool equalValueRow (bool t,
                      size_t row,
                      Real value) const;
    // Return: True if each element in row is undefined or equal to value
  bool equalValue (Real value) const;
    // Invokes: equalValueRow()
  size_t maxValueLenRow (bool t,
                         size_t row) const;
  size_t maxValueLen () const;
    // Invokes: maxValueLenRow()
  size_t* getMaxValueLenCol (bool t) const;
    // Size == rowsSize(t)
    // Invokes: maxValueLenRow()
  bool isInteger () const;
  // Square matrix
  bool isSquare () const
    { return rowsSize (false) == rowsSize (true); }
  bool isSymmetric (size_t &row,
             		    size_t &col) const;
    // Output: row, col - a non-symmetric element
    //           valid if result is false
  bool isSymmetric () const
    { size_t row, col; return isSymmetric (row, col); }
  bool zeroDiagonal (size_t &row) const;
    // Output: row - a non-0 diagonal element; 
    //           valid if result is false
    // Return: true if all elements of the diagonal are 0
    // Requires: isSquare()
  bool greaterDiagonal (bool t,
                        size_t &row,
                        size_t &col) const;
    // Output: row, col - element which is > getDiag(row); 
    //           valid if result is false
    // Return: true if all elements of the diagonal are >= the corresponding non-diagonal elements
    // Requires: isSquare()

  // defined (= !isNan())
  size_t getDefinedNumRow (bool t,
                           size_t row) const;
  size_t getDefinedNum () const;
  bool definedRow (bool t,
                   size_t row) const
		{ if (row >= rowsSize (t))
		    return false;
	    return getDefinedNumRow (t, row) == rowsSize (! t);
		}
  bool defined () const;
    // Return: true if all elements are defined
  bool existsMissing (bool t,
                      size_t &MissingRow,
	                    size_t &MissingCol) const;
		// Output: MissingRow|col
  bool emptyRow (bool t,
                 size_t row) const
    { return getDefinedNumRow (t, row) == 0; }
  size_t firstEmptyRow (bool t) const;
    // Return: If no empty rows then NO_INDEX
  bool noEmptyRows (bool t) const
    { return firstEmptyRow (t) == NO_INDEX; }
  bool empty () const;
  size_t firstColDefined (bool t,
                          size_t row) const
    { Real r; return nextColDefined (t, row, 0, r); }
  size_t lastColDefined (bool t,
                         size_t row) const
    { Real r; return prevColDefined (t, row, rowsSize (! t) - 1, r); }
  size_t nextColDefined (bool  t,
	                       size_t  row,
	                       size_t  StartCol,
	                       Real &R) const;
    // Return: column >= StartCol with a defined element;
    //         NO_INDEX if not found
    // Output: R if result != NO_INDEX
  size_t prevColDefined (bool  t,
	                       size_t  row,
	                       size_t  StartCol,
	                       Real &R) const;
    // Return: column <= StartCol with a defined defined; 
    //         NO_INDEX if not found
    // Output: R if result != NO_INDEX
  Matrix* deleteUndefinedRows (bool t) const;
    // Return: matrix with defined rows

  // Sum
  Real sumRow (bool t,
               size_t row) const;
  Real sum () const;
  void posNegSumRow (bool t,
                     size_t row,
                     Real &PosSum, 
                     Real &NegSum) const;
    // Sums of positive and negative elements
    // Output: Pos|NegSum
  void getPosNegSumCols (bool   t,
                         Matrix &tagret,
                         bool   tagretT) const;
    // Output: tagret
    //           col = 0: sum of positive values
    //           col = 1: sum of negative values
    // Requires: tagret.rowsSize(tagretT) = rowsSize(t)
    //           tagret.rowsSize(!tagretT] = 2
  Real getTrace () const;
    // Requires: isSquare()

  // Max|Min
  // If empty then -INF|INF
  Real getMaxRow (bool t,
                  size_t row) const;
  Real getMinRow (bool t,
                  size_t row) const;
  Real max () const;
  Real min () const;
  Real maxAbsRow (bool t,
                  size_t row) const;
  Real minAbsRow (bool t,
                  size_t row) const;
  Real maxAbs () const;
  Real minAbs () const;
  void posNegSumCols2MaxMin (bool t,
                             Real &PosMax, 
                             Real &NegMin) const;
    // Output: PosMax of PosSumCol, NegMin of NegSumCol
  Real closestZeroRow (bool t,
                       size_t row) const;
  Real firstNonZeroRow (bool t,
                        size_t row) const;
  size_t argMaxRow (bool t,
                    size_t row) const;
    // Return: the greatest column of row
    //         if no result then NO_INDEX
  size_t argMinRow (bool t,
                  size_t row) const;
    // Return: the smallest column of row;
    //         if no result then NO_INDEX
  size_t greatestColLess (bool  t,
	                        size_t  row,
	                        Real value) const;
    // Return: the greatest column col s.t. get (t, row, col) < value
    // Requires: elements of row increase
  bool finite () const
    { return    max () <  INF 
    	       && min () > -INF;
    }

  // Center
  void meanDefinedNumRow (bool  t,
                          size_t row,
                          Real &mean,
                          size_t &DefinedNum) const;
    // Output: mean, DefinedNum
    //         If DefinedNum = 0 then mean = NAN
  Real meanRow (bool t,
                size_t row) const
    { Real mean; size_t definedNum;
      meanDefinedNumRow (t, row, mean, definedNum);
		  return mean;
		}

  // Deviation
  Real sumSqrDevRow (bool t,
                     size_t row,
                     Real CenterValue) const;
    // Return: Sum of squares
  Real sumSqrDev (Real CenterValue) const;
  Real sumAbsDevRow (bool t,
                     size_t row,
                     Real CenterValue) const;
    // Return: Sum of absolute differences
  Real sumAbsDev (Real CenterValue) const;
  size_t meanVarianceRow (bool t,
	                        size_t row,
	                        bool  Biased,
	                        Real &mean, 
	                        Real &Variance) const;
    // Variance: denominator is (N - Correction); If too few data then INF
    // Return: number of defined
    // Invokes: meanDefinedNumRow ()

  // Deviation from 0
  Real sumSqrRow (bool t,
                  size_t row) const
    { return sumSqrDevRow (t, row, 0.0); }
  Real getNorm (bool t,
                size_t row) const
    { return sqrt (sumSqrRow (t, row)); }
  Real sumSqr () const
    { return sumSqrDev (0.0); }
  Real sumAbsRow (bool t,
                  size_t row) const
    { return sumAbsDevRow (t, row, 0.0); }
  Real sumAbs () const
    { return sumAbsDev (0.0); }

  // Entropy
  Real entropyRow (bool t,
                   size_t row) const;
    // Requires: the row is a distribution

  // Symmetric
  Real getCovariance (bool t,
                      size_t row1, 
                      size_t row2,
                      bool biased) const
    { return getCovarianceVec (*this, t, row1,
		                           *this, t, row2, 
		                           biased);
    }
  Real getCorrelation (bool t,
                       size_t row1, 
                       size_t row2) const
    { return getCorrelationVec (*this, t, row1,
                                *this, t, row2);
    }
  void makeCovarianceMatrix (bool   t,
                             Matrix &tagret,
                             bool   biased) const;
    // Covariance of rows
    // Output: tagret
    // Invokes: getCovariance ()
    // Requires: tagret.isSquare()
    //           tagret.rowsSize() == rowsSize(t)
  void covariance2correlation ();
    // Requires: a covariance matrix
    //           isSymmetric(), defined()
  bool isSimilarity () const;
  Real getBilinear (bool t,
	                  const Matrix &m1,
	                  bool t1,
	                  size_t row1,
	                  const Matrix &m2,
	                  bool t2,
	                  size_t row2) const;
  Real getMahalanobis (bool t,
                       const Matrix &m1,
		                   bool t1,
		                   size_t row1,
		                   const Matrix &m2,
		                   bool t2,
		                   size_t row2) const;
  bool checkInverse (const Matrix &inverseM,
                     bool         inverseMT) const;
    // Return: false if inverseM is not the inverse of *this
  Matrix* getDiagVector (bool VecT) const;
    // Requires: isSquare()
    // Invokes: diag2Row()
//Find minimal subsets of linearly dependent rows ??
    // Requires: isSquare()
  Determinant getDeterminant () const;
    // Requires: isSquare(), defined()
    // Invokes: source.gaussElimination ()
  bool getEigen (Eigen &eigen,
                 Real error,
				 	       size_t maxIter) const;
    // Return: Success
    // Input: error: max. average absolute difference between eigen.vec and (*this)*eigen.vec
    // Update: eigen.vec: normalized
    // Output: eigen.value - maximum absolute eigenvalue
    // Requires: eigen.vec is normalized
    //           eigen.vec.rowsSize(eigen.t) = rowsSize(false)
    //           isSymmetric(), defined()
  Matrix* getCholesky (bool t) const;
    // Return: Lower-triangle matrix
    //         (*Return)^t * (*Return)^(!t) = *this
    //         nullptr <=> numerical problem
    // Requires: psd

  // Dependence
  Real getChi2 () const;
  Real getLnFisherExact (bool oneTail) const;
    // If oneTail then the rejection region is low values of get(false,0,0)
    // Return: ln p-value
    // Requires: isInteger()
    //           2*2 or 2*3 matrix 
  
  // Comparison
  bool equalLen (bool t,
                 const Matrix &source,
                 bool sourceT) const
    { return    rowsSize (t)   == source. rowsSize (sourceT) 
             && rowsSize (! t) == source. rowsSize (! sourceT);
    }
  // Requires: this->equalLen(source,sourceT) 
  Real maxAbsDiff (bool         t,
                   const Matrix &source,
                   bool         sourceT) const;
  Real meanRelativeError (bool         t,
                          const Matrix &source,
                          bool         sourceT) const;

  // Update
  void copyDataCheck (bool         t,
                      const Matrix &source,
                      bool         sourceT);
    // Requires: equalLen()
  void copyData (const Matrix &source)
    { resize (false, source. rowsSize (false), source. rowsSize (true));
      copyDataCheck (false, source, false);
    }
  void copyShift (bool         t,
                  size_t       firstRow, 
                  size_t       firstCol,
                  const Matrix &source,
                  bool         sourceT); 
    // Copy source to *this s.t. source[0,0] -> *this[firstRow,firstCol]
    // Input: firstRow|col; may be < 0
  void copyRow (bool         t,
                size_t       row,
                const Matrix &source,
                bool         sourceT,
                size_t       sourceRow);
    // Requires: rowsSize(!t) = source.rowsSize(!sourceT)
  void compressRow (bool t,
                    size_t row);
  void copyDefined (bool         t,
                    const Matrix &source,
                    bool         sourceT);
    // Copy only defined elements
    // Requires: EqualLen ()
  void swap (bool t,
             size_t row1,
             size_t col1,
             size_t row2,
             size_t col2)
    { const Real f = get (t, row1, col1);
		  put (t, row1, col1, get (t, row2, col2));
		  put (t, row2, col2, f);
		}
  void swapRows (bool t,
                 size_t row1, 
                 size_t row2);
  void putRow (bool t,
               size_t row,
               Real r);
  void putAll (Real R);
  void clear ()
    { putAll (NAN); }
  size_t undefined2ValueRow (bool t,
                             size_t row,
                             Real R);
    // Return: # old undefined values in row
  size_t undefined2ValueAll (Real R);
    // Return: # old undefined values
  size_t undefined2mean ();
    // Return: # old undefined values
    // Replace undefined values by the average of row-wise and column-wise means
  void reverseOrderRow (bool t,
                        size_t row);
    // Reverse the order of columns in row
  void reverseOrder (bool t);
  void putInc (bool t,
               size_t row, 
               size_t col,
               Real delta);
    // put (get () + delta)
  void putIncRow (bool t,
                  size_t row,
                  Real delta);
  void putIncAll (Real delta);
  void putIncSeries (bool t,
                     size_t row,
                     Real InitValue, 
                     Real delta);
    // put series: InitValue, InitValue + delta, ...
  Matrix* centerRows (bool t);
    // Metric L_2
  void negateRow (bool t,
                  size_t row);
  void negate ();
  void absRow (bool t,
               size_t row);
  void absAll ();
  void roundRow (bool t,
                 size_t row);
  void roundAll ();
  void add (bool         t,
            const Matrix &source,
            bool         sourceT,
            Real         factor);
    // *this += factor * source
  void addRow (bool         t,
               size_t       row,
		           const Matrix &source,
		           bool         sourceT,
		           size_t       sourceRow,
		           Real         factor);
    // (*this)_row += factor * source_sourceRow
  void putProd (bool t,
                size_t row, 
                size_t col,
                Real factor)
    { Real a;
		  if (get (t, row, col, a))
		    put (t, row, col, a * factor);
    }
  void putProdRow (bool t,
                   size_t row,
                   Real factor);
  void putProdAll (Real factor);
  void multiply (bool         t,
                 const Matrix &m1, 
                 bool         t1,
                 const Matrix &m2,
                 bool         t2);
    // *this = m1 * m2
  void multiplyBilinear (bool         t,
                         const Matrix &m1, 
                         bool         t1,
                         const Matrix &m2,
                         bool         t2)
    { Matrix m (false, m1. rowsSize (t1), m2. rowsSize (! t2));
      m. multiply (false, m1, t1, m2, t2);
      multiply (t, m2, ! t2, m, false);
      psd = m1. psd;
    }
    // *this = m2' * m1 * m2
  void expRow (bool t,
               size_t row);
  void expAll ();
  void logRow (bool t,
               size_t row);
  void logAll ();
  bool maximize (bool t,
                 size_t row, 
                 size_t col,
                 Real value);
    // Return: true if rows^t[row][col] is maximized
  bool minimize (bool t,
                 size_t row, 
                 size_t col,
                 Real value);
    // Return: true if rows^t[row][col] is minimized
  size_t maximizeRow (bool t,
	                    size_t row, 
	                    Real value);
  size_t minimizeRow (bool t,
	                    size_t  row, 
	                    Real value);
  size_t maximizeAll (Real value);
  size_t minimizeAll (Real value);
  void putRandomRow (bool t,
                     size_t row,
                     Rand &rand);
    // Uniformly distributed on [-0.5, 0.5]
  bool subtractProjectionRow (bool         t,
			                        size_t       row,
			                        const Matrix &source,
			                        bool         sourceT,
			                        size_t       sourceRow);
	  // Return: success
	  // Requries: vector of sourceRow != 0
	  // Post-condition: (*this)_row is orthogonal to source_sourceRow
  void cumulate (bool t);
    // putInc(t,row,col,get(t,row,col-1))
    // Undefined elements are skipped
  Real balanceRow (bool t,
                   size_t row,
                   Real rowSum);
    // Make sumRow() = rowSum
    //   If nullReal(sumRow()) then put 0 into the row
    // Return: rowSum-sumRow(t,row)  (old sumRow() -??)
  Real balance (bool  t,
                Real rowSum);
    // Return: sum of sqr (balanceRow())
    // Invokes: balanceRow (t) for each row
  Real balanceVec (bool         t,
                   const Matrix &sumColVector,
                   bool         sumColT);
    // Return: sum of sqr (balanceRow())
    // Requires: rowsSize(t) <= sumColVector.rowsSize(sumColT)
  void expBalanceLogRow (bool t,
                         size_t row,
                         Real expRowSum);
    // Requires: expRowSum >= 0.0
    //           At least 1 element is defined
  void expBalanceLog (bool t,
                      Real expRowSum);
    // Requires: expRowSum >= 0.0
    //           At least 1 element in each row is defined
    // Invokes: expBalanceLogRow ()
  bool normalizeRow (bool t,
                     size_t row,
                     Real &sqrNorma);
    // Makes row*row^t = 1
    // Return: success (<=> sqrNorma > 0)
    // Output: srqNorma = old row*row^t
  bool normalizeRow (bool t,
                     size_t row)
    { Real sqrNorma;
    	return normalizeRow (t, row, sqrNorma);
    }
  bool normalize (bool t,
                  Real &sumSqrNorma);
    // Output: sumSqrNorma = sum of sqrNorma of each row
    // Return: true if all rows are normalized
    // Invokes: normalizeRow()
  void sqrRow (bool t,
               size_t row);
  void sqrAll ();
  void sqrtRow (bool t,
                size_t row);
  void sqrtAll ();
  Determinant gaussElimination (bool t);
    // Makes an upper-triangle matrix using elementary row operations
    // Return: determinant if isSqaure()
    // Requires: defined()
  Real ipfp (bool         t,
             const Matrix &sumRowVector, 
             bool         sumRowT,
             const Matrix &sumColVector,
             bool         sumColT,
             Real         minError,
             size_t       maxIter);
    // See Kiruta, Shevyakov
    // Update: *this: balance to sumColVector and sumRowVector
    // Return: Error: > 0
    // Requires: defined ()
    //             rowsSize(true)  = sumColVector.rowsSize(true)
    //             rowsSize(false) = sumRowVector.rowsSize(false)
    //             sumRow/ColVector.defined()
    // If sumRowVector.sum() != sumColVector.sum() then minError is corrected
  // Square matrix
  void diag2row (bool         t,
                 size_t       row,
                 const Matrix &source);
    // Requires: source.isSquare()
    //           source.rowsSize(false) = rowsSize(!t)
  void row2diag (const Matrix &source,
                 bool         sourceT,
                 size_t       sourceRow,
                 Real         factor);
    // put sourceRow * factor to *this as a diagonal
    // *this may be not square
    // Requires: source.rowsSize(!sourceT) = minSize()
  void addVecVecT (bool         t,
                   const Matrix &source,
                   bool         sourceT,
                   size_t       sourceRow,
                   Real         factor);
    // *this += factor * sourceRow * sourceRow^t
    // Requires: isQuare(), rowsSize(false) = source.rowsSize(!sourceT)
  bool solveSystem (bool   t,
                    size_t col,
                    Matrix &source,
                    bool   sourceT);
    // Solve a system of linear equations - rows of source
    // Output: column col - solution
    // Return: true if success
    // Invokes: source.gaussElimination ()
    // Requires: rowsSize(t) = source.rowsSize(sourceT)
    // If *this is p.d. then use getCholesky() ??
  Determinant inverse ();
    // Generalized elimination method with the selection of the leading element
    // Return: getDeterminant()
    // If result is 0.0 then no inversion is made
    // Requries: isSquare(), defined()
  Real symmetrize (Real &maxCorrection,
                   size_t &row_bad,
                   size_t &col_bad);
    // rows[i][j] = (rows[i][j] + rows[j][i]) / 2
    // Output: maxCorrection, row_bad, col_bad
    // Return: sum of corrections
    // Requires: isSquare()
  void lower2upper (bool t);
    // Requires: isSquare()
  // Symmetric matrix
  void putSquaredDistances (const Matrix& source,
                            bool          sourceT);
    // Output: *this contains squared distances in metric L_2 between rows of source
  void sqrDistance2centeredSimilarity (Matrix* &rowMean,
                                       Real    &totalMean)
    // |x_i - x_j|^2  ->  x_i^t x_j
    // Similarities (attributes) are centered
    { centerSimilarity (rowMean, totalMean);
      putProdAll (-0.5);
    }
  void similarity2sqrDistance ();
    // x_i^t x_j -> |x_i - x_j|^2
    // Requires: isSymmetric(), defined()
  void centerSimilarity (Matrix* &rowMean,
                         Real    &totalMean); 
    // Output: rowMean, totalMean - sufficient statistics
    //         rowMean rowsSize = rowsSize(false) * 1
    //         totalMean = mean of rowMean[]
    // Requires: isSymmetric(), defined()
};



struct MVector : Matrix
{
  explicit MVector (size_t maxRow = 0) 
    : Matrix (false, maxRow, 1)
    {}
//MVector (const MVector &vec) = default;
//MVector& operator= (const MVector&) = default;
//MVector (MVector&& other) noexcept = default;
//MVector& operator= (MVector&&) = default;
  explicit MVector (const Vector<Real> &vec)
    : Matrix (false, vec. size (), 1)
    { operator= (vec); }
  MVector& operator= (const Vector<Real> &vec);
  MVector* copy () const override
    { return new MVector (*this); }

 
  size_t size () const
    { return rowsSize (false); }
  void resize/*Vec*/ (size_t maxRow)
    { Matrix::resize (false, maxRow, 1); }
  Real operator[] (size_t row) const
    { return get (false, row, 0); }
  Real& operator[] (size_t row)
    { return rows [row] [0]; }
};



// Eigen

struct Eigen : Root
{
  Real value;
    // Init: 0
  MVector vec;
    // Normalized

  explicit Eigen (size_t vectorLen) 
    : value (0)
    , vec (vectorLen)
    {}
//Eigen (Eigen&&) noexcept = default;
//Eigen& operator= (Eigen&&) = default;
  void qc () const override;
  void saveText (ostream& os) const override;

  Real getNorm2 () const
    { return vec. sumSqr (); }  
  bool operator> (const Eigen &other) const
    { return abs (value) > abs (other. value); }
  bool makePositiveVec ();
    // Return: vec is negateRow()'d
    // Update: vec
};



// Eigens

struct Eigens : Root 
{
  // INPUT
  bool psd;
  Real error;
	Real totalExplained_max;
	  // > 0
	  
	// OUTPUT
  Matrix basis;
  MVector values;
    // size() == basis.rowsSize(true)
  Real explainedVarianceFrac;
    // May be: isNan()
  // psd => isProb(explainedVarianceFrac)
  Prob explainedFrac_next;
    // May be isNan()
  bool orthogonal;
    // false <=> the next Eigen::vector is not orthogonal to basis

 
  Eigens (const Matrix &matr,
          size_t dim_max,
          Prob totalExplainedFrac_max,
          Prob explainedFrac_min,
          Real relError,
          size_t maxIter);
    // Randomized algorithm
    // Input: matr: defined(), isSymmetric()
    // Output: totalExplained_max = psd ? matr.getTrace() : matr.sumSqr()
    //         basis: is maximum s.t.:
    //                  rowsSize(true) <= dim_max
    //                  totalExplainedFrac() <= totalExplainedFrac_max 
    //                  explainedFrac() > explainedFrac_min 
    // Invokes: matr.getEigen(error,maxIter)
//Eigens (const Eigens&) = default;
//Eigens& operator= (const Eigens&) = default;    
  Eigens* copy () const final
    { return new Eigens (*this); }
//Eigens (Eigens&&) noexcept = default;
//Eigens& operator= (Eigens&&) = default;    
  void qc () const override;
  void saveText (ostream& os) const override;
  JsonMap* toJson (JsonContainer* parent,
                   const string& name) const override;
    

  size_t getInitSize () const
    { return basis. rowsSize (false); }
  size_t getDim () const
    { return values. size (); }
    // Return: reason can be derived from the constructor parameters and OUTPUT
  Prob explainedFrac (size_t eigenNum) const
    { Real explained = values [eigenNum];
      if (! psd)  
        explained = sqr (explained);
      return explained / totalExplained_max; 
    }
  Prob cumulativeExplainedFrac (size_t eigenNum) const;
  Prob totalExplainedFrac () const
    { return values. empty () ? 0 : cumulativeExplainedFrac (values. size () - 1); }
  Real cumulativeValues (size_t eigenNum) const;
  Real cumulativeValues () const
    { return values. empty () ? 0 : cumulativeValues (values. size () - 1); }
  Prob cumulativeValuesFrac (size_t eigenNum) const
    { return values. empty () ? 0 : cumulativeValues (eigenNum) / cumulativeValues (values. size () - 1); }
  bool imaginaryMdsAttr (size_t attrNum) const
    { return values [attrNum] < 0; }
  Real getMds (size_t objNum,
				       size_t attrNum,
				       bool imaginary) const
    { return imaginaryMdsAttr (attrNum) == imaginary
               ? basis. get (false, objNum, attrNum) * sqrt (abs (values [attrNum])) 
               : 0;  
    }
	size_t makeSquare ();
	  // Return: # basis's (= # values's) added
	  //         > 0 => basis.getDeterminant() = 0
	void restore (Matrix &matr) const;
		// Output: matr = basis values^d basis'
		// Requires: matr: isSquare(), rowsSize() = getInitSize()
};



}



#endif

