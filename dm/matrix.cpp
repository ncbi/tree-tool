// matrix.cpp

#undef NDEBUG
#include "../common.inc"

#include "matrix.hpp"
using namespace std;
using namespace Common_sp;

#include <string>
#include <fstream>
#include <cmath>

 

namespace DM_sp
{



Real multiplyVec (const Matrix &m1,
                  bool t1,
                  size_t row1,
                  const Matrix &m2,
                  bool t2,
                  size_t row2)
{
  ASSERT (m1. rowsSize (! t1) == m2. rowsSize (! t2));

  Real s = 0.0;
  Real r1, r2;
  FOR (size_t, col, m1. rowsSize (! t1))
    if (m1. get (t1, row1, col, r1) &&
        m2. get (t2, row2, col, r2)
       )
      s += r1 * r2;
  return s;
}



Real sumSqrDifferenceVec (const Matrix &m1,
                          bool t1,
                          size_t row1,
                          const Matrix &m2,
                          bool t2,
                          size_t row2)
{
  ASSERT (m1. rowsSize (! t1) == m2. rowsSize (! t2));

  Real s = 0;
  Real r1, r2;
  FOR (size_t, col, m1. rowsSize (! t1))
    if (m1. get (t1, row1, col, r1) &&
        m2. get (t2, row2, col, r2)
       )
     s += sqr (r1 - r2);
  return s;
}



Real sumAbsDifferenceVec (const Matrix &m1,
                          bool         t1,
                          size_t       row1,
                          const Matrix &m2,
                          bool         t2,
                          size_t       row2)
{
  ASSERT (m1. rowsSize (! t1) == m2. rowsSize (! t2));

  Real s = 0;
  Real r1, r2;
  FOR (size_t, col, m1. rowsSize (! t1))
    if (m1. get (t1, row1, col, r1) &&
        m2. get (t2, row2, col, r2)
       )
     s += abs (r1 - r2);
  return s;
}



Real maxAbsDifferenceVec (const Matrix &m1,
                          bool         t1,
                          size_t       row1,
                          const Matrix &m2,
                          bool         t2,
                          size_t       row2)
{
  ASSERT (m1. rowsSize (! t1) == m2. rowsSize (! t2));

  Real s = 0;
  Real r1, r2;
  FOR (size_t, col, m1. rowsSize (! t1))
    if (m1. get (t1, row1, col, r1) &&
        m2. get (t2, row2, col, r2)
       )
     Common_sp::maximize (s, abs (r1 - r2));
  return s;
}



Real getCovarianceVec (const Matrix &m1,
                            bool      t1,
                            size_t         row1,
                            const Matrix &m2,
                            bool      t2,
                            size_t         row2,
                            bool      Biased)
{
  ASSERT (m1. rowsSize (! t1) == m2. rowsSize (! t2));

  size_t n = 0;
  Real Sxy = 0.0;
  Real Sx  = 0.0;
  Real Sy  = 0.0;
  Real r1, r2;
  FOR (size_t, col, m1. rowsSize (! t1))
    if (m1. get (t1, row1, col, r1) &&
        m2. get (t2, row2, col, r2))
    {
      n++;
      Sxy += r1 * r2;
      Sx  += r1;
      Sy  += r2;
    }

  const size_t Correction = ! Biased;
  if (n <= Correction)
    return NAN;
  const Real a = Sxy - Sx * Sy / (Real) n;
  return a / ((int) n - (int) Correction);
}



Real getCorrelationVec (const Matrix  &m1,
                        bool         t1,
                        size_t       row1,
                        const Matrix &m2,
                        bool         t2,
                        size_t       row2)
{
  ASSERT (m1. rowsSize (! t1) == m2. rowsSize (! t2));


  size_t n = 0;
  Real Sxy = 0.0;
  Real Sx  = 0.0;
  Real Sx2 = 0.0;
  Real Sy  = 0.0;
  Real Sy2 = 0.0;
  Real r1, r2;
  FOR (size_t, col, m1. rowsSize (! t1))
    if (m1. get (t1, row1, col, r1) &&
        m2. get (t2, row2, col, r2))
    {
      n++;
      Sxy += r1 * r2;
      Sx  += r1;
      Sx2 += sqr (r1);
      Sy  += r2;
      Sy2 += sqr (r2);
    }

  if (n == 0)
    return NAN;

  const Real a     = Sxy - Sx * Sy  / (Real) n;
  const Real Var1N = Sx2 - sqr (Sx) / (Real) n;
  const Real Var2N = Sy2 - sqr (Sy) / (Real) n;
  const Real b = Var1N * Var2N;
  if (leReal (b, 0.0))
    return 1.0;
  return a / sqrt (b);
}



void distribution2MeanVar (const Matrix &distribution,
                           size_t       distributionRow,
                           bool         distributionT,
                           const Matrix &Values,
                           size_t       ValuesRow,
                           bool         ValuesT,
                           Real         &mean, 
                           Real         &Variance)
{
  ASSERT (distribution. rowsSize (! distributionT) == Values. rowsSize (! ValuesT));
  ASSERT (distribution. definedRow (distributionT, distributionRow));
//ASSERT (Values.       definedRow (ValuesT, ValuesRow));  -??

  // mean
  mean = 0.0;
  FOR (size_t, col, distribution. rowsSize (! distributionT))
  {
    const Real r = distribution. get (distributionT, distributionRow, col);
    if (! nullReal (r))
      mean += r * Values. get (ValuesT, ValuesRow, col);
  }

  // StDev
  Variance = 0.0;
  FOR (size_t, col, distribution. rowsSize (! distributionT))
  {
    const Real r = distribution. get (distributionT, distributionRow, col);
    if (! nullReal (r))
      Variance += r * sqr (Values. get (ValuesT, ValuesRow, col) - mean);
  }
}




// Determinant

void Determinant::qc () const
{
  if (! qc_on)
    return;
  LogReal::qc ();
    
  ASSERT (size);
}




///////////////////////////////// Matrix ///////////////////////////////

void Matrix::init (bool t,
                   size_t maxRow,
                   size_t maxCol)
{
//ASSERT (rows. empty ());

  swapRowCol (t, maxRow, maxCol);
  
  rows. resize (maxRow);
  FOR (size_t, row, maxRow)
    rows [row]. resize (maxCol, NAN);
    
  colsSize = maxCol;
  
  psd = isSquare ();  
}



Matrix::Matrix (bool t,
                istream &f,
                bool square) 
{
  ASSERT (rows. empty ());

  size_t maxRow, maxCol;
  loadSize (f, square, maxRow, maxCol); 

  init (t, maxRow, maxCol);

  loadData (t, f);
}



void Matrix::loadSize (istream &f,
                       bool square,
                       size_t &maxRow,
                       size_t &maxCol)
{
  ASSERT (! f. bad ());

  f >> maxRow;
  if (! f. good ())
    ERROR;
  if (square)
    maxCol = maxRow;
  else
    if (! (f >> maxCol). good ())
      ERROR;
      
  psd = false;  // store in file ??
}



void Matrix::loadData (bool t,
                       istream &f)
{
  ASSERT (! f. bad ());

  FOR (size_t, row, rowsSize (false))
  FOR (size_t, col, rowsSize (true))
  {
    Real r;
    if (! (f >> r). good ())
      throw ios_base::failure (string("row ") + toString (row) + ", col " + toString (col));
    put (t, row, col, r);
  }
}



void Matrix::load (bool t,
                   istream &f,
                   bool square) 
{
  size_t maxRow, maxCol;
  loadSize (f, square, maxRow, maxCol); 

  swapRowCol (t, maxRow, maxCol);
  ASSERT (maxRow == rowsSize (false));
  ASSERT (maxCol == rowsSize (true));

  loadData (t, f);
}



void Matrix::reserve (bool t,
                      size_t maxRow,
                      size_t maxCol)
{
  swapRowCol (t, maxRow, maxCol);
  
  rows. reserve (maxRow);
  FOR (size_t, row, rowsSize (false))
    rows [row]. reserve (maxCol);
}



void Matrix::delData ()
{
  FOR (size_t, row, rowsSize (false))
    rows [row]. clear ();
  rows. clear ();
  
  colsSize = 0;
  psd = true;
}



void Matrix::qc () const
{
  if (! qc_on)
    return;
  IMPLY (psd, isSymmetric ());
}




// Output


namespace 
{
  const char* undefStr = "?";
}



void Matrix::printRow (bool t,
                       size_t row,
                       ostream &f) const
{
  Real r;
  FOR (size_t, col, rowsSize (! t))
  {
    if (col)
      f << '\t';

    if (get (t, row, col, r))
    {
      f << scientific; f. precision (3); f << r;  // PAR
    }
    else
      f << undefStr;
  }
}



void Matrix::saveFile_ (bool t,
                        bool square,
                        ostream &f) const
{
  ASSERT (! f. bad ());
  
  
  if (square)
  {
    ASSERT (isSquare ());
    f << rowsSize (false) << endl;
  }
  else
    f << rowsSize (false) << ' ' << rowsSize (true) << endl;

  FOR (size_t, row, rowsSize (t))
  {
    printRow (t, row, f);
    f << endl;
  }
}



void Matrix::insertRows (bool t,
	                       size_t firstRow,
	                       size_t rowInc)
{
  if (rowInc == 0)
    return;

  const size_t oldSize = rowsSize (t);
  ASSERT (firstRow <= oldSize);
  if (t)
  {
    FOR (size_t, row, rowsSize (false))
    {
    //const size_t oldRowSize = rowsSize (false); ??
      rows [row]. resize (oldSize + rowInc, NAN);  
      FOR_REV_END (size_t, col, firstRow, oldSize)
        rows [row] [col + rowInc] = rows [row] [col];
      FOR_START (size_t, col, firstRow, firstRow + rowInc)
        rows [row] [col] = NAN;
    }
    colsSize += rowInc;
  }
  else
  {
    const size_t colSize_ = rowsSize (true);
    rows. resize (oldSize + rowInc);  
    FOR_REV_END (size_t, row, firstRow, oldSize)
      rows [row + rowInc] = rows [row];      
    FOR_START (size_t, row, firstRow, firstRow + rowInc)
      if (row >= oldSize)
        rows [row]. resize (colSize_, NAN);
      else
        putRow (false, row, NAN);
  }
  
  psd = false;
}



void Matrix::deleteRows (bool t,
                         size_t firstRow,
                         size_t rowDec)
{
  if (rowDec == 0)
    return;

  ASSERT (firstRow + rowDec <= rowsSize (t));
  ASSERT (rowDec < rowsSize (t));
  
  if (t)
    FOR (size_t, row, rowsSize (false))
      rows [row]. eraseMany (firstRow, firstRow + rowDec);
  else
    rows. eraseMany (firstRow, firstRow + rowDec);

  psd = false;
}



void Matrix::trimRows (bool t,
                       size_t NewMaxRow)
{
  ASSERT (NewMaxRow <= rowsSize (t));
  
  deleteRows (t, NewMaxRow, rowsSize (t) - NewMaxRow);
  ASSERT (rowsSize (t) == NewMaxRow);

  psd = false;
}



void Matrix::deleteObject (size_t objNum)
{
  ASSERT (isSquare ());
  ASSERT (objNum < rowsSize (false));

  Keep<bool> kp (psd);
  for (const bool b : Bool)
    deleteRows (b, objNum, 1);
}




// processRow2Col

Matrix* Matrix::processRow2Col (bool           t,
                                ProcessRowFunc ProcessRow) const
{
  Matrix* target = new Matrix (t, rowsSize (t), 1);
  FOR (size_t, row, rowsSize (t))
  {
    const Real a = (this->*ProcessRow) (t, row);
    target->put (t, row, 0, a);
  }
    
  return target;
}




// Basic info

bool Matrix::equalValueRow (bool t,
                            size_t row,
                            Real value) const
{
  Real r;
  FOR (size_t, col, rowsSize (! t))
    if (   get (t, row, col, r) 
    	  && ! eqReal (r, value)
    	 )
      return false;
  return true;
}



bool Matrix::equalValue (Real value) const
{
  FOR (size_t, row, rowsSize (false))
    if (! equalValueRow (false, row, value))
      return false;
  return true;
}



size_t Matrix::maxValueLenRow (bool t,
                               size_t row) const
{
  size_t ValueLen = string (undefStr). size ();
  Real r;
  FOR (size_t, col, rowsSize (! t))
    if (get (t, row, col, r))
    {
      ostringstream oss;
      oss << r;
      Common_sp::maximize (ValueLen, oss. str (). size ());
    }
  return ValueLen;
}



size_t Matrix::maxValueLen () const
{
  size_t ValueLen = 0;
  FOR (size_t, row, rowsSize (false))
    Common_sp::maximize (ValueLen, maxValueLenRow (false, row));
  return ValueLen;
}



size_t* Matrix::getMaxValueLenCol (bool t) const
{
  size_t* maxValueLenCol = new size_t [rowsSize (t)];
  FOR (size_t, row, rowsSize (t))
    maxValueLenCol [row] = maxValueLenRow (t, row);
  return maxValueLenCol;
}



bool Matrix::isInteger () const
{
	Real x;
  FOR (size_t, row, rowsSize (false))
  FOR (size_t, col, rowsSize (true))
    if (get (false, row, col, x))
      if (! DM_sp::isInteger (x))
    	  return false;
  return true;
}



bool Matrix::isSymmetric (size_t &row,
				                  size_t &col) const
{
  if (! isSquare ())
  	return false;

  const Real delta = std::max (1.0, maxAbs ()) * DM_sp::epsilon; 
  Real a, b;
  for (row = 0;       row < rowsSize (false); row++)
  for (col = row + 1; col < rowsSize (false); col++)
    if (   get (false, row, col, a) 
        && get (false, col, row, b) 
        && ! eqReal (a, b, delta)
       )
    {
    //if (verbose ())
      //cout << row << ' ' << col << ' ' << a << ' ' << b << ' ' << delta << endl;
     	return false;
    }
  return true;
}



bool Matrix::zeroDiagonal (size_t &row) const
{
  ASSERT (isSquare ());

  Real a;
  for (row = 0; row < rowsSize (false); row++)
    if (  getDiag (row, a) &&
        ! nullReal (a)
       )
      return false;
  return true;
}



bool Matrix::greaterDiagonal (bool t,
                              size_t    &row,
                         		  size_t    &col) const
{
  ASSERT (isSquare ());

  Real a, b;
  for (row = 0; row < rowsSize (false); row++)
  for (col = 0; col < rowsSize (false); col++)
    if (getDiag (row, a)     &&
        get (t, row, col, b) &&
        lessReal (a, b)
       )
      return false;
  return true;
}



size_t Matrix::getDefinedNumRow (bool t,
                                 size_t row) const
{
  size_t s = 0;
  Real r;
  FOR (size_t, col, rowsSize (! t))
    if (get (t, row, col, r))
      s++;
  return s;
}



size_t Matrix::getDefinedNum () const
{
  size_t s = 0;
  FOR (size_t, row, rowsSize (false))
    s += getDefinedNumRow (false, row);
  return s;
}



bool Matrix::defined () const
{
  FOR (size_t, row, rowsSize (false))
    if (! definedRow (false, row))
      return false;
  return true;
}



bool Matrix::existsMissing (bool t,
                            size_t &MissingRow,
				                    size_t &MissingCol) const
{
  Real r;
  FOR (size_t, row, rowsSize (t))
  FOR (size_t, col, rowsSize (! t))
    if (! get (t, row, col, r))
    {
      MissingRow = row;
      MissingCol = col;
      return true;
    }
  return false;
}



size_t Matrix::firstEmptyRow (bool t) const
{
  FOR (size_t, row, rowsSize (t))
    if (emptyRow (t, row))
      return row;

  return NO_INDEX;
}



bool Matrix::empty () const
{
  FOR (size_t, row, rowsSize (false))
    if (! emptyRow (false, row))
      return false;
  return true;
}



size_t Matrix::nextColDefined (bool  t,
	                             size_t  row,
	                             size_t  StartCol,
	                             Real &r) const
{
  FOR_START (size_t, col, StartCol, rowsSize (! t))
    if (get (t, row, col, r))
      return col; 
  return NO_INDEX;
}



size_t Matrix::prevColDefined (bool  t,
	                             size_t  row,
	                             size_t  StartCol,
	                             Real &r) const
{
  FOR_REV (size_t, col, StartCol + 1)
    if (get (t, row, col, r))
      return col; 
  return NO_INDEX;
}



Matrix* Matrix::deleteUndefinedRows (bool t) const
{
  Matrix* M = 0;
  {
    size_t s = 0;
    FOR (size_t, row, rowsSize (t))
      if (definedRow (t, row))
        s++;
    M = new Matrix (false, s, rowsSize (! t));
  }

  {
    size_t s = 0;
    FOR (size_t, row, rowsSize (t))
      if (definedRow (t, row))
      {
        M->copyRow (false, s, *this, t, row);
        s++;
      }
  }
    
  return M;
}



Real Matrix::sumRow (bool t,
                     size_t row) const
{
  Real s = 0.0;
  Real r;
  FOR (size_t, col, rowsSize (! t))
    if (get (t, row, col, r))
      s += r;
  return s;
}



Real Matrix::sum () const
{
  Real s = 0.0;
  FOR (size_t, row, rowsSize (false))
    s += sumRow (false, row);
  return s;
}



void Matrix::posNegSumRow (bool  t,
                           size_t row,
                           Real &PosSum, 
                           Real &NegSum) const
{
  PosSum = 0.0;
  NegSum = 0.0;
  Real r;
  FOR (size_t, col, rowsSize (! t))
    if (get (t, row, col, r))
    {
      if (r > 0.0)
        PosSum += r;
      else
        NegSum += r;
    }
}



void Matrix::getPosNegSumCols (bool   t,
                               Matrix &tagret,
                               bool   tagretT) const
{
  ASSERT (tagret. rowsSize (  tagretT) == rowsSize (t));
  ASSERT (tagret. rowsSize (! tagretT) == 2);
  
  FOR (size_t, row, rowsSize (t))
  {
    Real PosSum, NegSum;
    posNegSumRow (t, row, PosSum, NegSum);
    tagret. put (tagretT, row, 0, PosSum);
    tagret. put (tagretT, row, 1, NegSum);
  }
}



Real Matrix::getTrace () const
{
  ASSERT (isSquare ());
  
  Real s = 0.0;
  FOR (size_t, row, rowsSize (false))
    s += get (false, row, row);
  return s;
}



Real Matrix::getMaxRow (bool t,
                        size_t row) const
{
  Real r = - INF;
  Real a;
  FOR (size_t, col, rowsSize (! t))
    if (get (t, row, col, a))
      Common_sp::maximize (r, a);
  return r;
}



Real Matrix::getMinRow (bool t,
                        size_t row) const
{
  Real r = INF;
  Real a;
  FOR (size_t, col, rowsSize (! t))
    if (get (t, row, col, a))
      Common_sp::minimize (r, a);
  return r;
}



Real Matrix::max () const
{
  Real maxValue = - INF;
  FOR (size_t, row, rowsSize (false))
    Common_sp::maximize (maxValue, getMaxRow (false, row));
  return maxValue;
}



Real Matrix::min () const
{
  Real minValue = INF;
  FOR (size_t, row, rowsSize (false))
    Common_sp::minimize (minValue, getMinRow (false, row));
  return minValue;
}



Real Matrix::maxAbsRow (bool t,
                        size_t row) const
{
  Real r = - INF;
  Real a;
  FOR (size_t, col, rowsSize (! t))
    if (get (t, row, col, a))
      Common_sp::maximize (r, abs (a));
  return r;
}



Real Matrix::minAbsRow (bool t,
                        size_t row) const
{
  Real r = INF;
  Real a;
  FOR (size_t, col, rowsSize (! t))
    if (get (t, row, col, a))
      Common_sp::minimize (r, abs (a));
  return r;
}



Real Matrix::maxAbs () const
{
  Real maxValue = - INF;
  FOR (size_t, row, rowsSize (false))
    Common_sp::maximize (maxValue, maxAbsRow (false, row));
  return maxValue;
}



Real Matrix::minAbs () const
{
  Real minValue = INF;
  FOR (size_t, row, rowsSize (false))
    Common_sp::minimize (minValue, minAbsRow (false, row));
  return minValue;
}



void Matrix::posNegSumCols2MaxMin (bool t,
                                   Real   &PosMax, 
                                   Real   &NegMin) const
{
  Matrix PosNegSumCols (t, rowsSize (t), 2);
  getPosNegSumCols (t, PosNegSumCols, t);
  PosMax = PosNegSumCols. getMaxRow (! t, 0);
  NegMin = PosNegSumCols. getMinRow (! t, 1);
}



Real Matrix::closestZeroRow (bool t,
                             size_t row) const
{
  Real r = INF;
  Real a;
  FOR (size_t, col, rowsSize (! t))
    if (get (t, row, col, a))
      Common_sp::minimize (r, abs (a));
  return r;
}



Real Matrix::firstNonZeroRow (bool t,
                              size_t row) const
{
  Real a;
  FOR (size_t, col, rowsSize (! t))
    if (   get (t, row, col, a) 
        && ! nullReal (a)
       )
      return a;
  return 0.0;
}



size_t Matrix::argMaxRow (bool t,
                          size_t row) const
{
  size_t max_ = NO_INDEX;
  Real r = - INF;
  Real a;
  FOR (size_t, col, rowsSize (! t))
    if (   get (t, row, col, a) 
        && Common_sp::maximize (r, a)
       )
      max_ = col;
  return max_;
}



size_t Matrix::argMinRow (bool t,
                          size_t row) const
{
  size_t min_ = NO_INDEX;
  Real r = INF;
  Real a;
  FOR (size_t, col, rowsSize (! t))
    if (   get (t, row, col, a) 
        && Common_sp::minimize (r, a)
       )
      min_ = col;
  return min_;
}



size_t Matrix::greatestColLess (bool  t,
	                              size_t row,
	                              Real value) const
{
  Real r;
  FOR (size_t, col, rowsSize (! t))
    if (   get (t, row, col, r) 
        && r >= value
       )
      return col - 1;
  return rowsSize (! t) - 1;
}



void Matrix::meanDefinedNumRow (bool t,
                                size_t row,
                                Real &mean,
                                size_t &definedNum) const
{
  definedNum = 0;
  Real s = 0.0;
  Real r;
  FOR (size_t, col, rowsSize (! t))
    if (get (t, row, col, r))
    {
      definedNum++;
      s += r;
    }

  mean = definedNum ? s / (Real) definedNum : NAN;
}



Real Matrix::sumSqrDevRow (bool  t,
                           size_t  row,
                           Real centerValue) const
{
  Real s = 0.0;
  Real a;
  FOR (size_t, col, rowsSize (! t))
    if (get (t, row, col, a))
      s += sqr (a - centerValue);
  return s;
}



Real Matrix::sumSqrDev (Real centerValue) const
{
  Real s = 0.0;
  FOR (size_t, row, rowsSize (false))
    s += sumSqrDevRow (false, row, centerValue);
  return s;
}



Real Matrix::sumAbsDevRow (bool t,
                           size_t row,
                           Real centerValue) const
{
  Real s = 0.0;
  Real a;
  FOR (size_t, col, rowsSize (! t))
    if (get (t, row, col, a))
      s += abs (a - centerValue);
  return s;
}



Real Matrix::sumAbsDev (Real centerValue) const
{
  Real s = 0.0;
  FOR (size_t, row, rowsSize (false))
    s += sumAbsDevRow (false, row, centerValue);
  return s;
}



size_t Matrix::meanVarianceRow (bool  t,
	                              size_t  row,
	                              bool  Biased,
	                              Real &mean, 
	                              Real &Variance) const
{
  size_t DefinedNum;
  meanDefinedNumRow (t, row, mean, DefinedNum);

  const size_t Correction = ! Biased;
  if (DefinedNum <= Correction)
    Variance = INF;
  else 
    Variance = sumSqrDevRow (t, row, mean) / ((int) DefinedNum - (int) Correction);
    
  return DefinedNum;
}



Real Matrix::entropyRow (bool t,
                         size_t row) const
{
  Real s = 0.0;
  Real r;
  FOR (size_t, col, rowsSize (! t))
    if (   get (t, row, col, r) 
        && ! nullReal (r)
       )
      s += r * prob2info (r);
  return s;
}



void Matrix::makeCovarianceMatrix (bool   t,
                                   Matrix &tagret,
                                   bool   Biased) const
{
  ASSERT (tagret. isSquare ());
  ASSERT (rowsSize (t) == tagret. rowsSize (false));

  FOR (size_t, row1, rowsSize (t))
    FOR_START (size_t, row2, row1, rowsSize (t))
      tagret. putSymmetric (row1, row2, getCovariance (t, row1, row2, Biased));
}



void Matrix::covariance2correlation ()
{
  ASSERT (isSquare ()); 

  Matrix StDev (false, 1, rowsSize (true));
  FOR (size_t, row, rowsSize (false))
    StDev. put (false, 0, row, sqrt (get (false, row, row)));

  FOR (size_t, row, rowsSize (false))
  {
    putDiag (row, 1);
    FOR_START (size_t, col, row + 1, rowsSize (true))
    {
      const Real a = StDev. get (false, 0, row) *
                     StDev. get (false, 0, col);
      putSymmetric (row, col, nullReal (a) ? 1 : (get (false, row, col) / a));
    }
  }
}



bool Matrix::isSimilarity () const
{
	if (! isSquare ())
		return false;
	if (! isSymmetric ())
		return false;
		
	FOR (size_t, row, rowsSize (false))
  	FOR (size_t, col, rowsSize (true))
  	  if (greaterReal (sqr (get (false, row, col)), getDiag (row) * getDiag (col)))
  	  {
  	  //if (verbose ())
  	    //cout << row << ' ' << col << ' ' << sqr (get (false, row, col)) << ' ' << getDiag (row) * getDiag (col) << endl;
  	    return false;
  	  }
	  	
	return true;
}



Real Matrix::getBilinear (bool t,
		                      const Matrix &m1,
				                  bool t1,
				                  size_t row1,
				                  const Matrix &m2,
				                  bool t2,
				                  size_t row2) const
{
  ASSERT (m1. rowsSize (! t1) == rowsSize (t));
  ASSERT (m2. rowsSize (! t2) == rowsSize (! t));
  
  Real s = 0;
  FOR (size_t, row, rowsSize (t)) 
  FOR (size_t, col, rowsSize (! t))
    s +=       get (t, row, col) 
         * m1. get (t1, row1, row) 
         * m2. get (t2, row2, col);
  return s;
}



Real Matrix::getMahalanobis (bool t,
			                       const Matrix &m1,
					                   bool t1,
					                   size_t row1,
					                   const Matrix &m2,
					                   bool t2,
					                   size_t row2) const
{
	ASSERT (m1. rowsSize (! t1) == m2. rowsSize (! t2));
	
  MVector x (m1. rowsSize (! t1));
  FOR (size_t, col, m1. rowsSize (! t1))
    x [col] =   m1. get (t1, row1, col) 
              - m2. get (t2, row2, col);
	
	return getBilinear ( t
	                   , x, true, 0
	                   , x, true, 0
	                   );
}
	


bool Matrix::checkInverse (const Matrix &inverseM,
                           bool         inverseMT) const
{
  ASSERT (isSquare ());
  ASSERT (defined ());
  ASSERT (inverseM. isSquare ());
  ASSERT (inverseM. defined ());
  ASSERT (rowsSize (false) == inverseM. rowsSize (false));


  Matrix Ident (false, rowsSize (false), rowsSize (false));
  
  Ident. multiply (false, *this, false, inverseM, inverseMT);
  
  const Real error = 0.1;  // PAR ??
  FOR (size_t, row, rowsSize (false))
    FOR (size_t, col, rowsSize (false))
      if (row == col)
      {
        if (abs (Ident. get (false, row, col) - 1.0) >= error)
          return false;
      }
      else
        if (abs (Ident. get (false, row, col)) >= error)
          return false;
          
  return true;
}



Matrix* Matrix::getDiagVector (bool VecT) const
{
  ASSERT (isSquare ());

  Matrix* diagVec = new Matrix (VecT, rowsSize (false), 1);
  diagVec->diag2row (! VecT, 0, *this);
  return diagVec;
}



Determinant Matrix::getDeterminant () const
{
  ASSERT (isSquare ());

  Matrix matr (*this);
  return matr. gaussElimination (false);
}



bool Matrix::getEigen (Eigen &eigen,
                       Real error,
				 	             size_t maxIter) const
{
	ASSERT (defined ());
	{
  	size_t row = NO_INDEX;
  	size_t col = NO_INDEX;
    if (! isSymmetric (row, col))
    {
      cout << row << " " << col << endl;
      print (cout);
      ERROR;
    }
  }
  ASSERT (eigen. vec. size () == rowsSize (false));
	ASSERT (eqReal (eigen. vec. sumSqr (), 1));
  ASSERT (! negative (error));
  ASSERT (maxIter >= 1);


  MVector vec (rowsSize (false));
  Real diff = INF;  
  FOR (size_t, iteration, maxIter)
  {                                        
    vec. multiply (false, 
                   *this, false, 
                   eigen. vec, false);
      
    // eigen.value
    WeightedMeanVar mv;
    FOR (size_t, row, vec. size ())
    {
      const Real a =        vec [row];
      const Real b = eigen. vec [row];
      if (! nullReal (b))
        mv. add (a / b, abs (a) + abs (b));
    }
    eigen. value = mv. getMean ();
      
    Real sqrNorma;  // = eigenvalue^2
    if (! vec. normalizeRow (true, 0, sqrNorma))
    {
      eigen. vec. putAll (0); 
      eigen. value = 0;
      return false;
    }
    
    diff = maxAbsDifferenceVec (       vec, true, 0,
                                eigen. vec, true, 0);
                      
    eigen. vec. copyDataCheck (false, vec, false);
    
    if (leReal (diff, error))
    {
      if (psd)
      {
        if (negative (eigen. value, 1e-3))  // PAR
        {
          cout << "eigen.value = " << eigen. value << endl;
          ERROR;
        }
        Common_sp::maximize (eigen. value, 0.0);
      }
      eigen. makePositiveVec ();
      return true;
    }
  }
  

  if (verbose ())
  {
	  eigen. vec. print (cout);         
	  cout << "diff = " << diff << "  error = " << error << endl;
	}
  return false;
}



Real Matrix::getChi2 () const
{
  ASSERT (minSize () >= 2); 
	ASSERT (defined ());
	ASSERT (! negative (min ()));
	
	Real chi2 = 0;
  Common_sp::AutoPtr <Matrix> rowSum (processRow2Col (false, & Matrix::sumRow));
  Common_sp::AutoPtr <Matrix> colSum (processRow2Col (true,  & Matrix::sumRow));
  const Real s = rowSum->sumRow (true, 0);
  ASSERT (s > 0);
  FOR (size_t, row, rowsSize (false))
    FOR (size_t, col, rowsSize (true))
    	if (const Real expect =   rowSum->get (false, row, 0) 
    		                      * colSum->get (false, 0, col)
    		                      / s
    		  )
    	  chi2 += sqr (get (false, row, col) - expect) / expect;
    	  
  return chi2;
}



Real Matrix::getLnFisherExact (bool oneTail) const
// Table 2*x ==> Multivariate hypergeometric distribution
{
  ASSERT (minSize () >= 2); 
  ASSERT (defined ());
  ASSERT (! negative (min ())); 
  ASSERT (isInteger ());
  
  
  Real threshold = 0;
  if (! oneTail)
	  FOR (size_t, row, rowsSize (false))
	  FOR (size_t, col, rowsSize (true))
	    threshold += lnFactorial ((uint) round (get (false, row, col)));
  
  
  const uint sum00 = (uint) round (sum ());


  if (rowsSize (false) == 2 &&
      rowsSize (true)  == 2)
    {
      SumLn sumLn;
      const uint sum01 = (uint) round (sumRow (true,  0));
      const uint sum02 = (uint) round (sumRow (true,  1));
      const uint sum10 = (uint) round (sumRow (false, 0));
      const uint sum20 = (uint) round (sumRow (false, 1));
      ASSERT (sum01 + sum02 == sum00);
      ASSERT (sum10 + sum20 == sum00);
      const uint sum11_max = oneTail ? (uint) round (get (false, 0, 0)) : std::min (sum01, sum10);
      FOR_START (  uint, sum11
                , sum10 <= sum02 ? 0 : sum10 - sum02
                , sum11_max + 1
                )
      {
        const Real f =   lnFactorial (sum11)
                       + lnFactorial (sum01 - sum11)
                       + lnFactorial (sum10 - sum11)
                       + lnFactorial (sum02 - (sum10 - sum11));
        if (oneTail || geReal (f, threshold))
          sumLn. addExp (- f);
      }
      return
          lnFactorial (sum01)
        + lnFactorial (sum02)
        + lnFactorial (sum10)
        + lnFactorial (sum20)
        - lnFactorial (sum00)
        + sumLn. getLn ();  
    }
  

  if (   ! oneTail 
  	  && rowsSize (false) == 2 
      && rowsSize (true)  == 3
     )
    {
      SumLn sumLn;
      const uint sum01 = (uint) round (sumRow (true,  0));
      const uint sum02 = (uint) round (sumRow (true,  1));
      const uint sum03 = (uint) round (sumRow (true,  2));
      const uint sum10 = (uint) round (sumRow (false, 0));
      const uint sum20 = (uint) round (sumRow (false, 1));
      ASSERT (sum01 + sum02 + sum03 == sum00);
      ASSERT (sum10 + sum20 == sum00);
      FOR (uint, sum11, std::min (sum01, sum10) + 1)
      FOR_START (uint, sum12, 
                 sum10 <= (sum11 + sum03) ? 0 : sum10 - (sum11 + sum03),
                 std::min (sum02, sum10 - sum11) + 1
                )
      {
        const uint sum23 = sum03 - (sum10 - (sum11 + sum12));
        const Real f =   lnFactorial (sum11)
                        + lnFactorial (sum12)
                        + lnFactorial (sum10 - (sum11 + sum12))
                        + lnFactorial (sum01 - sum11)
                        + lnFactorial (sum02 - sum12)
                        + lnFactorial (sum23);
        if (geReal (f, threshold))
          sumLn. addExp (- f);
      }
      return
          lnFactorial (sum01)
        + lnFactorial (sum02)
        + lnFactorial (sum03)
        + lnFactorial (sum10)
        + lnFactorial (sum20)
        - lnFactorial (sum00)
        + sumLn. getLn ();  
    }

    
  NOT_IMPLEMENTED;
  return NAN;
}



Real Matrix::maxAbsDiff (bool         t,
                         const Matrix &source,
                         bool         sourceT) const
{
  ASSERT (equalLen (t, source, sourceT));

  Real diff = 0;
  Real a, b;
  FOR (size_t, i, rowsSize (t)) 
  FOR (size_t, j, rowsSize (! t))
    if (   get (t, i, j, a)
    	  && source. get (sourceT, i, j, b)
    	 )
      Common_sp::maximize (diff, abs (a - b));
  return diff;
}



Real Matrix::meanRelativeError (bool         t,
                                const Matrix &source,
                                bool         sourceT) const
{
  ASSERT (equalLen (t, source, sourceT));

  size_t C = 0;
  Real s = 0.0;
  Real r;
  FOR (size_t, i, rowsSize (t)) 
	  FOR (size_t, j, rowsSize (! t))
	    if (get (t, i, j, r))
	    {
	      if (r == 0.0) 
	        return INF;
	      Real X;
	      if (source. get (sourceT, i, j, X))
	      {
	        C++;
	        s += abs ((X - r) / r);
	      }
	    }
  if (C == 0)
    return INF;

  return s / (Real) C;
}



void Matrix::copyDataCheck (bool t,
                            const Matrix &source,
                            bool sourceT)
{
  ASSERT (equalLen (t, source, sourceT));
  copyShift (t, 0, 0, source, sourceT);
  psd = source. psd;
}



void Matrix::copyShift (bool         t,
                        size_t       firstRow, 
                        size_t       firstCol,
                        const Matrix &source,
                        bool         sourceT)
{
  FOR (size_t, sourceRow, source. rowsSize (sourceT))
  {
    const size_t row = firstRow + sourceRow;
    if (between (row, (size_t) 0, rowsSize (t)))
	    FOR (size_t, sourceCol, source. rowsSize (! sourceT))
	    {
	      const size_t col = firstCol + sourceCol;
	      if (between (col, (size_t) 0, rowsSize (! t)))
  	      put (t, row, col, source. get (sourceT, sourceRow, sourceCol));
	    }
  }
}



void Matrix::copyRow (bool         t,
                      size_t       row,
                      const Matrix &source,
                      bool         sourceT,
                      size_t       sourceRow)
{
  ASSERT (rowsSize (! t) == source. rowsSize (! sourceT));
  FOR (size_t, col, rowsSize (! t))
    put (t, row, col, source. get (sourceT, sourceRow, col));
}



void Matrix::compressRow (bool t,
                          size_t row)
{
  FOR_START (size_t, i, row + 1, rowsSize (t))
    copyRow (t, i - 1, *this, t, i);
}



void Matrix::copyDefined (bool      t,
                          const Matrix &source,
                          bool      sourceT)
{
  ASSERT (equalLen (t, source, sourceT));

 
  Real r;
  FOR (size_t, row, rowsSize (t))
  FOR (size_t, col, rowsSize (! t))
    if (source. get (sourceT, row, col, r))
      put (t, row, col, r);
}



void Matrix::swapRows (bool t,
                           size_t    row1, 
                           size_t    row2)
{
  if (row1 == row2)
    return;

  FOR (size_t, col, rowsSize (! t))
  {
    const Real r = get (t, row1, col);
    put (t, row1, col, get (t, row2, col));
    put (t, row2, col, r);
  }
}



void Matrix::putRow (bool t,
                     size_t row,
                     Real r)
{
  FOR (size_t, col, rowsSize (! t))
    put (t, row, col, r);
}



void Matrix::putAll (Real r)
{
  FOR (size_t, row, rowsSize (false))
    putRow (false, row, r);
  psd = isSquare () && (isNan (r) || r >= 0);
}



size_t Matrix::undefined2ValueRow (bool t,
	                                 size_t row,
	                                 Real r)
{
  ASSERT (! isNan (r));

  size_t n = 0;
  FOR (size_t, col, rowsSize (! t))
    if (isNan (get (t, row, col)))
    {
      put (t, row, col, r);
      n++;
    }

  return n;  
}



size_t Matrix::undefined2ValueAll (Real r)
{ 
  size_t n = 0;
  FOR (size_t, row, rowsSize (false))
    n += undefined2ValueRow (false, row, r);
  return n;
}



size_t Matrix::undefined2mean ()
{ 
  size_t n = 0;
  Common_sp::AutoPtr <Matrix> rowMean (processRow2Col (false, & Matrix::meanRow));
  Common_sp::AutoPtr <Matrix> colMean (processRow2Col (true,  & Matrix::meanRow));
  FOR (size_t, row, rowsSize (false))
    FOR (size_t, col, rowsSize (true))
      if (isNan (get (false, row, col)))
      {
        const Real r = (  rowMean->get (false, row, 0) 
                        + colMean->get (true,  col, 0)
                       ) / 2;
        put (false, row, col, r);
        n++;
      }
  return n;
}



void Matrix::reverseOrderRow (bool t,
                              size_t row)
{
  FOR (size_t, col, rowsSize (! t) / 2)
    swap ( t
         , row, col
         , row, rowsSize (! t) - 1 - col
         );
}



void Matrix::reverseOrder (bool t)
{
  FOR (size_t, row, rowsSize (t))
    reverseOrderRow (t, row);
}



void Matrix::putInc (bool  t,
                     size_t  row, 
                     size_t  col,
                     Real delta)
{
  Real a;
  if (get (t, row, col, a))
    put (t, row, col, a + delta);
}



void Matrix::putIncRow (bool  t,
                        size_t  row,
                        Real delta)
{
  FOR (size_t, col, rowsSize (! t))
    putInc (t, row, col, delta);
}



void Matrix::putIncAll (Real delta)
{
  Keep<bool> kp (psd);
  FOR (size_t, row, rowsSize (false))
    putIncRow (false, row, delta);
}



void Matrix::putIncSeries (bool  t,
                           size_t  row,
                           Real InitValue, 
                           Real delta)
{
  Real r = InitValue;
  FOR (size_t, col, rowsSize (! t))
  {
    put (t, row, col, r);
    r += delta;
  }
}



Matrix* Matrix::centerRows (bool t)
{
  Matrix* rowMean = processRow2Col (t, & Matrix::meanRow);
  FOR (size_t, row, rowsSize (t))
    putIncRow (t, row, - rowMean->get (t, row, 0));
  return rowMean;
}



void Matrix::negateRow (bool t,
                        size_t row)  
{
  Real a;
  FOR (size_t, col, rowsSize (! t))
    if (get (t, row, col, a))
      put (t, row, col, - a);
}



void Matrix::negate ()
{
  FOR (size_t, row, rowsSize (false))
    negateRow (false, row);
}



void Matrix::absRow (bool t,
                     size_t row)  
{
  Real a;
  FOR (size_t, col, rowsSize (! t))
    if (get (t, row, col, a))
      put (t, row, col, abs (a));
}



void Matrix::absAll ()
{
  FOR (size_t, row, rowsSize (false))
    absRow (false, row);
}



void Matrix::roundRow (bool t,
                       size_t row)  
{
  Real a;
  FOR (size_t, col, rowsSize (! t))
    if (get (t, row, col, a))
      put (t, row, col, std::round (a));
}



void Matrix::roundAll ()
{
  FOR (size_t, row, rowsSize (false))
    roundRow (false, row);
}



void Matrix::add (bool         t,
                  const Matrix &source,
                  bool         sourceT,
                  Real         factor)
{
  ASSERT (equalLen (t, source, sourceT));
  
  const bool psd_new = psd && source. psd && factor >= 0;

  FOR (size_t, row, rowsSize (t)) 
    addRow (         t,       row
           , source, sourceT, row
           , factor
           );

  psd = psd_new;
}



void Matrix::addRow (bool         t,
			               size_t       row,
					           const Matrix &source,
					           bool         sourceT,
					           size_t       sourceRow,
					           Real         factor)
{
  ASSERT (equalLen (t, source, sourceT));

  FOR (size_t, col, rowsSize (! t))
    putInc (t, row, col, factor * source. get (sourceT, sourceRow, col));  	
}



void Matrix::putProdRow (bool  t,
                         size_t  row,
                         Real factor)
{
  FOR (size_t, col, rowsSize (! t))
    putProd (t, row, col, factor);
}



void Matrix::putProdAll (Real factor)
{
  Keep<bool> kp (psd);
  FOR (size_t, row, rowsSize (false))
    putProdRow (false, row, factor);
}



void Matrix::multiply (bool         t,
                       const Matrix &m1, 
                       bool         t1,
                       const Matrix &m2,
                       bool         t2)
{
  ASSERT (this != & m1);
  ASSERT (this != & m2);
  ASSERT (    rowsSize (t)    == m1. rowsSize (t1));
  ASSERT (    rowsSize (! t)  == m2. rowsSize (! t2));
  ASSERT (m1. rowsSize (! t1) == m2. rowsSize (t2));

  FOR (size_t, row, rowsSize (t)) 
  FOR (size_t, col, rowsSize (! t))
    put (t, row, col, multiplyVec (m1, t1,   row,
                                   m2, ! t2, col));
  psd = & m1 == & m2 && t1 != t2;
}



void Matrix::expRow (bool t,
                     size_t row)
{
  Real r;
  FOR (size_t, col, rowsSize (! t))
    if (get (t, row, col, r))
      put (t, row, col, exp (r));
}



void Matrix::expAll ()
{
  FOR (size_t, row, rowsSize (false))
    expRow (false, row);
}



void Matrix::logRow (bool t,
                     size_t row)
{
  Real r;
  FOR (size_t, col, rowsSize (! t))
    if (get (t, row, col, r))
      put (t, row, col, log (r));
}



void Matrix::logAll ()
{
  FOR (size_t, row, rowsSize (false))
    logRow (false, row);
}



bool Matrix::maximize (bool  t,
                       size_t  row, 
                       size_t  col,
                       Real value)
{
  Real a;
  if (! get (t, row, col, a))
    return false;
  if (! Common_sp::maximize (a, value))
    return false;

  put (t, row, col, a);
  return true;
}



bool Matrix::minimize (bool  t,
                       size_t  row, 
                       size_t  col,
                       Real value)
{
  Real a;
  if (! get (t, row, col, a))
    return false;
  if (! Common_sp::minimize (a, value))
    return false;

  put (t, row, col, a);
  return true;
}



size_t Matrix::maximizeRow (bool t,
	                          size_t row, 
	                          Real value)
{
  size_t n = 0;
  FOR (size_t, col, rowsSize (! t))
    if (maximize (t, row, col, value))
      n++;
  return n;
}



size_t Matrix::minimizeRow (bool t,
	                          size_t row, 
	                          Real value)
{
  size_t n = 0;
  FOR (size_t, col, rowsSize (! t))
    if (minimize (t, row, col, value))
      n++;
  return n;
}



size_t Matrix::maximizeAll (Real value)
{
  size_t n = 0;
  FOR (size_t, row, rowsSize (false))
    n += maximizeRow (false, row, value);
  return n;
}



size_t Matrix::minimizeAll (Real value)
{
  size_t n = 0;
  FOR (size_t, row, rowsSize (false))
    n += minimizeRow (false, row, value);
  return n;
}



void Matrix::putRandomRow (bool t,
                           size_t row,
                           Rand &rand)
{
  FOR (size_t, col, rowsSize (! t))
    put (t, row, col, rand. getProb () - 0.5);
}



bool Matrix::subtractProjectionRow (bool         t,
						                        size_t       row,
						                        const Matrix &source,
						                        bool         sourceT,
						                        size_t       sourceRow)
{
  ASSERT (rowsSize (! t) == source. rowsSize (! sourceT));
	
	const Real s = source. sumSqrRow (sourceT, sourceRow);
	if (! positive (s))
	  return false;
	  
	const Real prod = multiplyVec (*this, t, row, source, sourceT, sourceRow);
	addRow (t, row, source, sourceT, sourceRow, - prod / s);
	
	return true;
}



void Matrix::cumulate (bool t)
{
  FOR (size_t, row, rowsSize (t))
  {
    Real currentSum = 0.0;
    Real a;
    FOR (size_t, col, rowsSize (! t))
      if (get (t, row, col, a))
      {
        currentSum += a;
        put (t, row, col, currentSum);
      }
  }
}



Real Matrix::balanceRow (bool  t,
                         size_t  row,
                         Real rowSum)
{
  const Real s = sumRow (t, row);
  if (nullReal (s))
    putRow (t, row, 0.0);
  else
    putProdRow (t, row, rowSum / s);
  return rowSum - s;
}



Real Matrix::balance (bool  t,
                      Real rowSum)
{
  Real s = 0.0;
  FOR (size_t, row, rowsSize (t))
    s += sqr (balanceRow (t, row, rowSum));
  return s;
}



Real Matrix::balanceVec (bool         t,
                         const Matrix &sumColVector,
                         bool         sumColT)
{
  ASSERT (sumColVector. rowsSize (! sumColT) == 1);
  ASSERT (rowsSize (t) <= sumColVector. rowsSize (sumColT));

  Real sumSqrDeviation = 0.0;
  FOR (size_t, row, rowsSize (t))
  {
    Real r;
    EXEC_ASSERT (sumColVector. get (sumColT, row, 0, r));
    sumSqrDeviation += sqr (balanceRow (t, row, r));
  }
  return sumSqrDeviation;
}



void Matrix::expBalanceLogRow (bool  t,
                               size_t  row,
                               Real expRowSum)
{
  ASSERT (expRowSum >= 0);
  
  const Real maxRow = getMaxRow (t, row);
  putIncRow (t, row, - maxRow);

  // s
  Real s = 0.0;
  Real r;
  FOR (size_t, col, rowsSize (! t))
    if (get (t, row, col, r))
      s += exp (r);
  if (s < 1.0)
  {
    cout << s << endl;
    ERROR;
  }
  
  const Real delta = log (expRowSum) - log (s);
  putIncRow (t, row, delta);
}



void Matrix::expBalanceLog (bool  t,
                            Real expRowSum)
{
  FOR (size_t, row, rowsSize (t))
    expBalanceLogRow (t, row, expRowSum);
}



bool Matrix::normalizeRow (bool  t,
                           size_t row,
                           Real &sqrNorma)
{
  sqrNorma = sumSqrRow (t, row);
  
  const Real norm = sqrt (sqrNorma);
  if (nullReal (norm))
    return false;

  putProdRow (t, row, 1 / norm);
  return true;
}



bool Matrix::normalize (bool  t,
                        Real &sumSqrNorma)
{     
  bool Ok = true;
  sumSqrNorma = 0.0;
  FOR (size_t, row, rowsSize (t))
  {
    Real sqrNorma;
    if (! normalizeRow (t, row, sqrNorma))
      Ok = false;
    sumSqrNorma += sqrNorma;
  }
    
  return Ok;
}



void Matrix::sqrRow (bool t,
                     size_t row)
{
  Real a;
  FOR (size_t, col, rowsSize (! t))
    if (get (t, row, col, a))
      put (t, row, col, sqr (a));
}



void Matrix::sqrAll ()
{
  FOR (size_t, row, rowsSize (false))
    sqrRow (false, row);
}



void Matrix::sqrtRow (bool t,
                      size_t row)
{
  Real a;
  FOR (size_t, col, rowsSize (! t))
    if (get (t, row, col, a))
      put (t, row, col, sqrt (a));
}



void Matrix::sqrtAll ()
{
  FOR (size_t, row, rowsSize (false))
    sqrtRow (false, row);
}



Determinant Matrix::gaussElimination (bool t)
{
  ASSERT (defined ());

  Determinant d (std::min (rowsSize (t), rowsSize (! t)));
  // Use partial pivoting ??
  FOR (size_t, row, d. size)
  {
    size_t row1 = row;
    while (row1 < rowsSize (t) &&
           nullReal (get (t, row1, row)))
      row1++;
    if (row1 < rowsSize (t))
    {
      // eor 1
      if (row != row1)
      {
        swapRows (t, row, row1);
        d *= -1;
      }
      // eor 2
      ASSERT (! nullReal (get (t, row, row)));      
      FOR_START (size_t, row2, row + 1, rowsSize (t))
      {
        const Real k = get (t, row2, row) / get (t, row, row);
        if (! nullReal (k))
	        FOR_START (size_t, col, row, rowsSize (! t))
	          putInc (t, row2, col, - get (t, row, col) * k);
      }
    }
      
    d *= get (t, row, row);
  }
  return d;
}



Real Matrix::ipfp (bool         t,
                   const Matrix &sumRowVector, 
                   bool         sumRowT,
                   const Matrix &sumColVector,
                   bool         sumColT,
                   Real         minError,
                   size_t       maxIter)
{
  ASSERT (rowsSize (t)   == sumColVector. rowsSize (sumColT));
  ASSERT (rowsSize (! t) == sumRowVector. rowsSize (! sumRowT));
  ASSERT (sumRowVector. rowsSize (sumRowT)   == 1);
  ASSERT (sumColVector. rowsSize (! sumColT) == 1);
  ASSERT (defined ());

  Common_sp::maximize (minError, abs (sumRowVector. sum () - sumColVector. sum ()));

  Real error = NAN;
  FOR (size_t, i, maxIter)
  {
    const Real sumDeviation = balanceVec (  t, sumColVector,   sumColT) +
                              balanceVec (! t, sumRowVector, ! sumRowT);
    error = sqrt (sumDeviation);
    if (error <= minError)
      break;
  }
  return error;
}



void Matrix::diag2row (bool t,
                       size_t row,
                       const Matrix &source)
{
  ASSERT (source. isSquare ());
  ASSERT (source. rowsSize (false) == rowsSize (! t));

  FOR (size_t, col, rowsSize (! t))
    put (t, row, col, source. getDiag (col));
}



void Matrix::row2diag (const Matrix &source,
                       bool sourceT,
                       size_t sourceRow,
                       Real factor)
{
  ASSERT (source. rowsSize (! sourceT) == minSize ());
  
  FOR (size_t, col, source. rowsSize (! sourceT))
    putDiag (col, source. get (sourceT, sourceRow, col) * factor);
}



void Matrix::addVecVecT (bool         t,
                         const Matrix &source,
                         bool         sourceT,
                         size_t       sourceRow,
                         Real         factor)
{
  ASSERT (isSquare ());
  ASSERT (rowsSize (false) == source. rowsSize (! sourceT));
  
  Keep<bool> kp (psd);
  FOR (size_t, row, rowsSize (false))
  {
    const Real a = factor * source. get (sourceT, sourceRow, row);
    FOR (size_t, col, rowsSize (false))
      putInc (t, row, col, a * source. get (sourceT, sourceRow, col));
  }
}



bool Matrix::solveSystem (bool   t,
                          size_t col,
                          Matrix &source,
                          bool   sourceT)
{
  ASSERT (rowsSize (t) == source. rowsSize (sourceT));
  ASSERT (source. rowsSize (! sourceT) == source. rowsSize (sourceT) + 1);

  source. gaussElimination (sourceT);
  FOR_REV (size_t, row, rowsSize (t))
  {
    Real s = 0.0;
    FOR_START (size_t, i, row + 1, rowsSize (t))
      s += source. get (sourceT, row, i) * get (t, i, col);
    const Real r = source. get (sourceT, row, source. rowsSize (! sourceT) - 1) - s;
    const Real a = source. get (sourceT, row, row);
    if (nullReal (a))
      if (nullReal (r))
        put (t, row, col, 0.0);  // Solution is not unique
      else 
        return false;
    else 
      put (t, row, col, r / a);
  }

  return true;
}



Matrix* Matrix::getCholesky (bool t) const
// Numerical Recipes in C, p. 97
{
	ASSERT (isSymmetric ());
	ASSERT (defined ());
	ASSERT (psd);
	
	Common_sp::AutoPtr <Matrix> m (new Matrix (false, *this, false));
	
	const Real abs_max = maxAbs ();
	FOR (size_t, row, rowsSize (false))
	  FOR (size_t, col, rowsSize (false))
	    if (col < row)
	    	m->put (t, col, row, 0);
	    else
		  {
		  	Real s = get (false, row, col);
		  	FOR (size_t, k, row) 
		  	  s -=   m->get (t, row, k) 
		  	       * m->get (t, col, k);
		  	if (row == col)
		  		if (negative (s, abs_max * epsilon))
		  		{
		  		#if 1
		  		  return nullptr;  // Numerical problem
		  		#else
		  		  cout << s << endl; 
		  		  print (cout); 
		  			ERROR;
		  		#endif
		  	  }		  		
		  		else
		  		  m->putDiag (row, sqrt (std::max (s, 0.0)));
		  	else
		  	{
    	    const Real d = m->getDiag (row);
		  	  if (nullReal (s, abs_max * epsilon) || nullReal (d, abs_max * epsilon))  // nullReal(d) => nullReal(s)
		  		  m->put (t, col, row, 0);
		  	  else
		  		  m->put (t, col, row, s / d);
		    }
		  }
		  
  ASSERT (m. get ());
  return m. release ();
}



Determinant Matrix::inverse ()
{
  ASSERT (isSquare ());
  ASSERT (defined ());

  // Inversion may work incorrectly if there is linear dependence -??
  // Determinant
  const Determinant determinant = getDeterminant ();
//if (determinant. get () == 0)
  if (nullReal (determinant. getAverage ()))
    return determinant; 


  Keep<bool> kp (psd);

  // M
  const size_t M = rowsSize (false);

  // U, V
  Vector <size_t> U (M);
  Vector <size_t> V (M);
  FOR (size_t, i, M)
  {
    U [i] = NO_INDEX;
    V [i] = NO_INDEX;
  }


  FOR (size_t, k, M)
  {
    // r, Q, coeff
   	// abs(get(false,r,Q)) -> max, s.t. U[r] = NO_INDEX, V[Q] = NO_INDEX
    Real XMax = -INF;
    size_t r = NO_INDEX;
    size_t q = NO_INDEX;
    FOR (size_t, j, M)
      if (U [j] == NO_INDEX)
        FOR (size_t, i, M)
          if (V [i] == NO_INDEX)
            if (Common_sp::maximize (XMax, abs (get (false, i, j))))
            {
              r = i;
              q = j;
            }
    if (q == NO_INDEX || 
        r == NO_INDEX) 
      return Determinant (rowsSize (false), 0);
    //throw runtime_error ("Inverse: no r or q");
    const Real f = get (false, r, q);
    // XMax = abs (f)
    U [q] = r;
    V [r] = q;
    if (nullReal (f)) 
      return Determinant (rowsSize (false), 0);
    //throw runtime_error ("Inverse: f=0");
    const Real coeff = 1.0 / f;

    // data
    FOR (size_t, j, M)
      if (j != q)
        {
          putProd (false, r, j, coeff);
          FOR (size_t, i, M)
            if (i != r)
              putInc (false, i, j, - get (false, r, j) * get (false, i, q));
        }
    put (false, r, q, coeff);
    FOR (size_t, i, M)
      if (i != r)
        putProd (false, i, q, - coeff);
  }

  FOR (size_t, k, M)
    if (k != U [k])
    {
      FOR (size_t, j, M)
        swap (false, k, j, U [k], j);

      FOR (size_t, i, M)
        swap (false, i, k, i, V [k]);

      U [V [k]] = U [k];
      V [U [k]] = V [k];
      U [k] = k;
      V [k] = k;
    }

  return determinant;
}



/*
//Inversion by Cholesky decomposition

  // a
  Matrix a (false, "a.mat");
  ASSERT (a. isSquare ());
  // Requires: a is symmentric
  
  
  // U
  Matrix U (false, "U.mat");
  ASSERT (U. isSquare ());
  // Requires: a = U^t * U
  // Cholesky decomposition


  // C = U^{-1}
  // Make in time O(n^2) -??
  Matrix C (false, U, false);
  FOR (size_t, col, C. rowsSize (true))
  {
    Matrix b (false, U. rowsSize (false), U. rowsSize (true) + 1);
    
    b. copyShift (false, 0, 0, U, false);
    b. putRow (true, b. rowsSize (true) - 1, 0.0);
    b. put (false, col, b. rowsSize (true) - 1, 1.0);
    
    C. SolveSystem (false, col, b, false);
  }
  

  // D = a^{-1}
  Matrix D (false, U, false);
  D. multiply (false, C, false, C, true);
*/



Real Matrix::symmetrize (Real &maxCorrection,
                         size_t &row_bad,
                         size_t &col_bad)
{
  ASSERT (isSquare ());


  Keep<bool> kp (psd);
  
  maxCorrection  = 0;
  Real corrected = 0;
  Real a, b;
  FOR (size_t, row, rowsSize (false))
    FOR_START (size_t, col, row + 1, rowsSize (false))
      if (get (false, row, col, a) &&
          get (false, col, row, b)
         )
      {
       	const Real c = (a + b) / 2;
       	const Real correction = abs (a - c);
       	if (Common_sp::maximize (maxCorrection, correction))
       	{
       	  row_bad = row;
       	  col_bad = col;
       	}
       	corrected += 2 * correction;
       	putSymmetric (row, col, c);
      }
      else
        if      (get (false, row, col, a))
          put (false, col, row, a);
        else if (get (false, col, row, a))
          put (false, row, col, a);

  return corrected;
}



void Matrix::lower2upper (bool t)
{
	ASSERT (isSquare ());
	
	FOR (size_t, row, rowsSize (false))
  	FOR (size_t, col, row)
  	  put (t, col, row, get (t, row, col));
}



void Matrix::putSquaredDistances (const Matrix& source,
                                  bool          sourceT)
{
  ASSERT (isSquare ());
  ASSERT (rowsSize (false) == source. rowsSize (false));
  
  FOR (size_t, row1, source. rowsSize (false))
  {
    putDiag (row1, 0.0);
    FOR_START (size_t, row2, row1 + 1, source. rowsSize (false))
      putSymmetric (row1, row2, sumSqrDifferenceVec (source, sourceT, row1,
                                                     source, sourceT, row2));
  }
}



#if 0
Real Matrix::sqrDistance2centeredSimilarity ()
{
  ASSERT (isSymmetric ());
	ASSERT (defined ());

  // sum_i x_i = 0
  
  Common_sp::AutoPtr <Matrix> rowSum (processRow2Col (false, & Matrix::sumRow));
  
  // = sum_i |x_i|^2
  const Real normSum = rowSum->sum () / (2 * (Real) rowsSize (false));
  rowSum->putIncAll (- normSum);
  rowSum->putProdAll (1.0 / (Real) rowsSize (false));
  // rowSum: |x_i|^2

  // x_i^t * x_j
  FOR (size_t, row, rowsSize (false))
  FOR (size_t, col, rowsSize (false))
    put (false, row, col, (  rowSum->get (false, row, 0) 
                           + rowSum->get (false, col, 0) 
                           - get (false, row, col)
                          ) / 2
        );

  Real defect = 0;
  Real diagSum = 0;
  FOR (size_t, row, rowsSize (false))
  {
    Real x = getDiag (row);    
    if (x < 0)
    {
      defect -= x;
      x = 0;
      putDiag (row, x);
    }
    diagSum += x;
  }
  
  return defect / diagSum;
}
#endif



void Matrix::similarity2sqrDistance ()
{
  ASSERT (isSymmetric ());
	ASSERT (defined ());

  Common_sp::AutoPtr <Matrix> DiagVec (getDiagVector (false));
  FOR (size_t, row, rowsSize (false))
  FOR (size_t, col, rowsSize (false))
    put (false, row, col,   DiagVec->get (false, row, 0) 
                          + DiagVec->get (false, col, 0) 
                          - 2 * get (false, row, col)
        );
}



void Matrix::centerSimilarity (Matrix* &rowMean,
                               Real &totalMean)
{
  ASSERT (isSymmetric ());
	ASSERT (defined ());
      
  rowMean = processRow2Col (false, & Matrix::sumRow);
  rowMean->putProdAll (1.0 / (Real) rowsSize (false));

  totalMean = rowMean->sum () / (Real) rowsSize (false);

  FOR (size_t, row, rowsSize (false))
  FOR (size_t, col, rowsSize (false))
    putInc (false, row, col, (  totalMean
                              - rowMean->get (false, row, 0) 
                              - rowMean->get (false, col, 0)));
}




// MVector

MVector& MVector::operator= (const Vector<Real> &vec)
{
  ASSERT (size () == vec. size ());
  FOR (size_t, i, size ())
    (*this) [i] = vec [i];
  return *this;
}




// Eigen

void Eigen::qc () const
{
  if (! qc_on)
    return;
  Root::qc ();
    
  ASSERT (eqReal (getNorm2 (), 1));
}



void Eigen::saveText (ostream& os) const
{
  os << "value = " << fixed; os. precision (4); os << value << endl;  
  os << "Vector: ";
  vec. saveFile_ (true, false, os);
}



bool Eigen::makePositiveVec ()
{
  FOR (size_t, i, vec. size ())
  {
    const Real a = vec [i];
    if (nullReal (a))
      continue;
    if (a < 0)
    {
      vec. negateRow (true, 0);
      return true;
    }
    else
      return false;
  }
  
  return false;
}




// Eigens

Eigens::Eigens (const Matrix &matr,
                size_t dim_max,
                Prob totalExplainedFrac_max,
                Prob explainedFrac_min,
                Real relError,
                size_t maxIter)
: psd (matr. psd)
, error (relError / sqrt ((Real) matr. rowsSize (false)))
, totalExplained_max (matr. psd ? matr. getTrace () : matr. sumSqr ())
, basis (false, matr. rowsSize (false), 0)
, values (0)
, explainedVarianceFrac (NAN)
, explainedFrac_next (NAN)
, orthogonal (true)
{
  if (verbose ())
    matr. print (cout);
  ASSERT (matr. defined ());
  ASSERT (matr. isSymmetric ());
  ASSERT (isProb (totalExplainedFrac_max));  
  ASSERT (isProb (explainedFrac_min));
  ASSERT (relError >= 0.0);
  
  
  const Real maxAbs = matr. maxAbs ();
  if (verbose ())
    cout << "maxAbs = " << maxAbs 
         << "  totalExplained_max = " << totalExplained_max
         << endl;
  if (nullReal (maxAbs))
    return;  
  if (   nullReal (totalExplained_max, 1e-6)  // PAR
      && maxAbs <= 1e-3  // PAR
     )
    return;

  ASSERT (positive (totalExplained_max));
  
  
  Unverbose unv;
  

  VectorOwn<Eigen> vecs;  
    // Eigen::vector's are orthogonal
  Rand rand;
  Matrix work (matr);
  Real totalExplained = 0;
  Common_sp::AutoPtr<Eigen> eigen;
  const size_t len = matr. rowsSize (false);
  if (verbose ())
    cout << "dim_max = " << dim_max << "  len = " << len << endl;
  Progress prog ((uint) min (dim_max, len), len >= 300 ? 1 : 0);  // PAR
  while (vecs. size () < min (dim_max, len))
  {
    prog ();

    Eigen* prevEigen = eigen. get ();  

    // work: removing the previous eigenvector
    if (prevEigen)
    {
      work. addVecVecT (false, 
                        prevEigen->vec, true, 0, 
                        - prevEigen->value);
      Real maxCorrection;
      size_t row_bad, col_bad;
      work. symmetrize (maxCorrection, row_bad, col_bad);
    }

    // eigen
    eigen. release ();
    eigen. reset (new Eigen (len));
    
    // eigen->vec: init
    // P (# Itertaions = 1) = 1
    do
    {
	    eigen->vec. putRandomRow (true, 0, rand);
	    for (const Eigen* other : vecs)
	      EXEC_ASSERT (eigen->vec. subtractProjectionRow (            true, 0,
	                                                      other->vec, true, 0));
    }
    while (! eigen->vec. normalizeRow (true, 0));  
    
    {
      Unverbose unv1;
      if (verbose ())
        work. saveText (cout);
    }

    ASSERT (work. psd == matr. psd);
    if (   ! work. getEigen (*eigen, error, maxIter)
        && ! eigen->getNorm2 ()
       )
    {
      if (verbose ())
      {
        eigen->print (cout);
        cout << "getEigen failed" << endl;
      }
      break;      
    }

    if (psd)
    { 
      if (negative (eigen->value / totalExplained_max, 1e-3))  // PAR
      {
        cout << eigen->value << ' ' << totalExplained_max << endl;
        ERROR;
      }
      maximize (eigen->value, 0.0);
    }
    
    eigen->qc ();

    const Real explained = psd ? eigen->value : sqr (eigen->value);
    ASSERT (explained >= 0);
    explainedFrac_next = explained / totalExplained_max;

    if (lessReal (explainedFrac_next, explainedFrac_min))
    {
      if (verbose ())
        cout << endl 
             << "  explainedFrac_next = " << explainedFrac_next
             << "  explainedFrac_min = " << explainedFrac_min 
             << endl;
      break;
    }
    
    totalExplained += explained;
    if (totalExplainedFrac_max < 1 && greaterReal (totalExplained, totalExplainedFrac_max * totalExplained_max, 1e-5))  // PAR
    {
      if (verbose ())
        cout << endl 
             << "  totalExplained = " << totalExplained
             << "  totalExplainedFrac_max * totalExplained_max = " << totalExplainedFrac_max * totalExplained_max
             << endl;
      break;
    }
    
  //IMPLY (prevEigen, leReal (eigen->value, prevEigen->value));  // May not hold

  #if 0  
    if (removeNegativeEigenValues && ! positive (eigen->value))
    {
      eigenValue_negative = true;
      break;
    }
  #endif

    // orthogonal
    for (const Eigen* other : vecs)
    {
    	const Real prod = abs (multiplyVec (other->vec, true, 0,
                                          eigen->vec, true, 0) 
									           / (Real) len
									          );
      if (greaterReal (prod, 0.05))  // PAR  
      {
        if (verbose ())
        	cout << "prod = " << scientific << prod 
        	         << "  vec # = " << vecs. size () 
        	         << "  dim_max = " << dim_max
        	         << "  len = " << len
        	         << endl;
    	  orthogonal = false;
    	  break;
    	}
    }
    if (! orthogonal)  // Try a different random normalized vector ??
      break;
  
    explainedFrac_next = NAN;
    vecs << eigen. get ();
  }
  eigen. release ();
  

  vecs. sortBubble ();    

  basis. resize (false, len, vecs. size ());
  values. resize (vecs. size ());
  explainedVarianceFrac = 0;
  FOR (size_t, col, vecs. size ())
  {
    FOR (size_t, row, len)
      basis. put (false, row, col, vecs [col] -> vec [row]);
    values [col] = vecs [col] -> value;
    explainedVarianceFrac += values [col];
  }   
  const Real tr = matr. getTrace ();
  ASSERT (! nullReal (tr));
//cout << vecs. size () << ' ' << tr << ' ' << explainedVarianceFrac << endl;  
  explainedVarianceFrac /= tr;
  if (psd)
  {
    ASSERT (positive (tr));
    minimize (explainedVarianceFrac, 1.0);
  }
}



void Eigens::qc () const
{
  if (! qc_on)
    return;
	ASSERT (basis. rowsSize (true) == getDim ());
	ASSERT (values. rowsSize (false) == getDim ());
	FOR (size_t, i, basis. rowsSize (true))
	{
	  ASSERT (eqReal (basis. sumSqrRow (true, i), 1));
	  IMPLY (psd, values [i] >= 0);
	}
  if (! leReal (totalExplainedFrac (), 1 + 1e-2))  // PAR
  {
    cout << totalExplainedFrac () << endl;
    ERROR;
  }
  if (! isNan (explainedVarianceFrac))
    if (psd && ! isProb (explainedVarianceFrac))
    {
      cout << explainedVarianceFrac << endl;
      ERROR;
    }
}



void Eigens::saveText (ostream& os) const
{
	os << fixed;
  os << "Matrix size: " << getInitSize () << endl;
  os << "Total explained fraction = ";
  os << fixed; os. precision (3); os << totalExplainedFrac () * 100 << " % (";
  os. width (8); os << fixed; os. precision (2); os << totalExplainedFrac () * totalExplained_max << " / ";
  os. width (8); os << fixed; os. precision (2); os << totalExplained_max << ")" << endl;

  os << "Total explained variance fraction = ";
  os << fixed; os. precision (3); os << explainedVarianceFrac * 100 << " %" << endl;
  
  FOR (size_t, i, getDim ())
  { 
    os << endl << "Eigen " << i + 1 << ':' << endl;
    os << "Explained fraction            = " << fixed; os. precision (3); os << explainedFrac (i) * 100 << " %" << endl;
    os << "Cumulative explained fraction = " << fixed; os. precision (3); os << cumulativeExplainedFrac (i) * 100 << " %" << endl;
    os << "Value = " << fixed; os. precision (4); os << values [i] << endl;
  }
  if (! isNan (explainedFrac_next))
  {
    os << endl << "Next explained fraction = " << fixed; os. precision (3); os << explainedFrac_next * 100 << " %" << endl;
  }
  os << "Orthogonal: " << orthogonal << endl;
//os << "Negative eigenvalue: " << eigenValue_negative << endl;
  
  os << endl;
  basis. saveText (os);
}



JsonMap* Eigens::toJson (JsonContainer* parent,
                         const string& name) const
{
  JsonMap* j = new JsonMap (parent, name);
  
  const size_t decimals = 6; // PAR

  new JsonBoolean (psd, j, "psd");

  {
    JsonArray* jArr = new JsonArray (j, "explained");  
    FOR (size_t, i, getDim ())
      new JsonDouble (explainedFrac (i) * (values [i] >= 0 ? 1 : -1), decimals, jArr); 
  }

  new JsonDouble (explainedFrac_next, decimals, j, "explainedFrac_next");
  new JsonBoolean (orthogonal, j, "orthogonal");
//new JsonBoolean (eigenValue_negative, j, "eigenValue_negative");
  new JsonDouble (explainedVarianceFrac, decimals, j, "explainedVarianceFrac");
  
  return j;
}



Prob Eigens::cumulativeExplainedFrac (size_t eigenNum) const
{
  Real s = 0;
  FOR (size_t, i, eigenNum + 1)
    s += explainedFrac (i);
  return s;
}



Real Eigens::cumulativeValues (size_t eigenNum) const
{
  Real s = 0;
  FOR (size_t, i, eigenNum + 1)
    s += values [i];
  return s;
}



size_t Eigens::makeSquare ()
{
  const size_t size_old = basis. rowsSize (true);
	ASSERT (basis. rowsSize (false) >= size_old);
	const size_t toAdd = basis. rowsSize (false) - size_old;

  // basis
  basis. insertRows (true, size_old, toAdd);
  Rand rand;
  FOR_START (size_t, i, size_old, basis. rowsSize (false))
    // Cf. Eigens::Eigens()
    // P (# Itertaions = 1) = 1
    do
    {
	    basis. putRandomRow (true, i, rand);
	    FOR (size_t, j, i)
	      EXEC_ASSERT (basis. subtractProjectionRow (       true, i,
	                                                 basis, true, j));
    }
    while (! basis. normalizeRow (true, i));  

    
  // values
  values. insertRows (false, size_old, toAdd);
  FOR (size_t, i, toAdd) 
    values [size_old + i] = 0;

	ASSERT (basis. isSquare ());
	
	return toAdd;
}



void Eigens::restore (Matrix &matr) const
{
	ASSERT (matr. isSquare ());
	ASSERT (matr. rowsSize (false) == getInitSize ());

	Matrix a (basis);
	FOR (size_t, col, values. size ())
	  a. putProdRow (true, col, values [col]);
	  
  matr. multiply (false, a, false, basis, true);
  matr. psd = true;
  
  Real maxCorrection;
  size_t row_bad, col_bad;
  matr. symmetrize (maxCorrection, row_bad, col_bad);
}



}



