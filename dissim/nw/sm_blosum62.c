/*  $Id: sm_blosum62.c 90506 2006-09-25 19:30:59Z madden $
* ===========================================================================
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
* Author:  Aaron Ucko, Mike Gertz
*
* Slava Brover: removed NCBI C ToolKit
*
* File Description:
*   Protein alignment score matrices; shared between the two toolkits.
*
* ===========================================================================
*/

#include "raw_scoremat.h"

/** Entries for the BLOSUM62 matrix at a scale of ln(2)/2.0. */

static const TNCBIScore s_Blosum62PSM[25 * 25] = {
    /*       A,  R,  N,  D,  C,  Q,  E,  G,  H,  I,  L,  K,  M,
             F,  P,  S,  T,  W,  Y,  V,  B,  J,  Z,  X,  *        */ 
    /*A*/    4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1,
            -2, -1,  1,  0, -3, -2,  0, -2, -1, -1, -1, -4,
    /*R*/   -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1,
            -3, -2, -1, -1, -3, -2, -3, -1, -2,  0, -1, -4,
    /*N*/   -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2,
            -3, -2,  1,  0, -4, -2, -3,  4, -3,  0, -1, -4,
    /*D*/   -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3,
            -3, -1,  0, -1, -4, -3, -3,  4, -3,  1, -1, -4,
    /*C*/    0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1,
            -2, -3, -1, -1, -2, -2, -1, -3, -1, -3, -1, -4,
    /*Q*/   -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0,
            -3, -1,  0, -1, -2, -1, -2,  0, -2,  4, -1, -4,
    /*E*/   -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2,
            -3, -1,  0, -1, -3, -2, -2,  1, -3,  4, -1, -4,
    /*G*/    0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3,
            -3, -2,  0, -2, -2, -3, -3, -1, -4, -2, -1, -4,
    /*H*/   -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2,
            -1, -2, -1, -2, -2,  2, -3,  0, -3,  0, -1, -4,
    /*I*/   -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,
             0, -3, -2, -1, -3, -1,  3, -3,  3, -3, -1, -4,
    /*L*/   -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,
             0, -3, -2, -1, -2, -1,  1, -4,  3, -3, -1, -4,
    /*K*/   -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1,
            -3, -1,  0, -1, -3, -2, -2,  0, -3,  1, -1, -4,
    /*M*/   -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,
             0, -2, -1, -1, -1, -1,  1, -3,  2, -1, -1, -4,
    /*F*/   -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,
             6, -4, -2, -2,  1,  3, -1, -3,  0, -3, -1, -4,
    /*P*/   -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2,
            -4,  7, -1, -1, -4, -3, -2, -2, -3, -1, -1, -4,
    /*S*/    1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1,
            -2, -1,  4,  1, -3, -2, -2,  0, -2,  0, -1, -4,
    /*T*/    0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1,
            -2, -1,  1,  5, -2, -2,  0, -1, -1, -1, -1, -4,
    /*W*/   -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,
             1, -4, -3, -2, 11,  2, -3, -4, -2, -2, -1, -4,
    /*Y*/   -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,
             3, -3, -2, -2,  2,  7, -1, -3, -1, -2, -1, -4,
    /*V*/    0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1,
            -1, -2, -2,  0, -3, -1,  4, -3,  2, -2, -1, -4,
    /*B*/   -2, -1,  4,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3,
            -3, -2,  0, -1, -4, -3, -3,  4, -3,  0, -1, -4,
    /*J*/   -1, -2, -3, -3, -1, -2, -3, -4, -3,  3,  3, -3,  2,
             0, -3, -2, -1, -2, -1,  2, -3,  3, -3, -1, -4,
    /*Z*/   -1,  0,  0,  1, -3,  4,  4, -2,  0, -3, -3,  1, -1,
            -3, -1,  0, -1, -2, -2, -2,  0, -3,  4, -1, -4,
    /*X*/   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
            -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -4,
    /***/   -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,
            -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1
};
const SNCBIPackedScoreMatrix NCBISM_Blosum62 = {
    "ARNDCQEGHILKMFPSTWYVBJZX*",
    s_Blosum62PSM,
    -4
};

