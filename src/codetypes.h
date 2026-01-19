// Copyright TIOGA Developers. See COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD 3-Clause)

#ifndef CODETYPES_H
#define CODETYPES_H

#define MPICH_SKIP_MPICXX
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// Use CMake-provided MPI includes
#include <mpi.h>
/*====================================================================*/
/*  Floating point definition                                         */
/*====================================================================*/
# define REAL double
/*====================================================================*/
/*  Base for indexing (0 or 1) */
/*====================================================================*/
# define BASE 1

/*====================================================================*/
/*  Define arithmetic constants                                       */
/*====================================================================*/
#define BIGVALUE           1.0e+15
#define BIGINT             2147483647
#define TOL                1.0e-10
#define HOLEMAPSIZE        192
// #define NFRINGE            3
// #define NVAR               6
/*==================================================================*/
/* inline debugging tools                                             */
/*==================================================================*/
# define TRACEI(x)  printf("#tioga:\t"#x" =%d\n",x);
# define TRACED(x)  printf("#tioga:\t"#x" =%.16e\n",x);
# define TIOGA_MIN(x,y)  (x) < (y) ? (x) : (y)
# define TIOGA_MAX(x,y)  (x) > (y) ? (x) : (y)
# define TIOGA_FREE(a1)  {free(a1);a1=NULL;}
/*====================================================================*/
/*  Numerical Tools                                                   */
/*====================================================================*/
#define TIOGA_Max(a1,a2) (((a1) >= (a2))? (a1): (a2))
#define TIOGA_Min(a1,a2) (((a1) <= (a2))? (a1): (a2))
/*********************************************************************/
/* Code specific types */
/*********************************************************************/
typedef struct HOLEMAP
{
  int existWall;
  int nx[3];
  int *samLocal;
  int *sam;
  double extents[6];
} HOLEMAP;


typedef struct OBB
{
  double xc[3];
  double dxc[3];
  double vec[3][3];

  int comm_idx;    /* Index in comm map for this OBB                       */
  int iblk_local;  /* Index of this mesh block                             */
  int iblk_remote; /* Index of the remote mesh block (intersecting  pair)  */
  int tag_remote;
  int send_tag;
  int recv_tag;
}OBB;
typedef struct DONORLIST
{
  int donorData[4];
  double donorRes;
  double receptorRes;
  int cancel;
  struct DONORLIST *next;
} DONORLIST;

typedef struct PACKET
{
  int nints;
  int nreals;
  int *intData;
  REAL *realData;
} PACKET;

typedef struct INTERPLIST
{
  int cancel;
  int nweights;
  int receptorInfo[3];
  double xtmp[3];
  int *inode;
  double *weights;
} INTERPLIST;

typedef struct INTERPLIST2
{
  int cancel;
  int nweights;
  int receptorInfo[3];
  int *inode;
  double *weights;
  struct INTERPLIST2 *next;
} INTERPLIST2;

typedef struct INTEGERLIST
{
  int inode;
  struct INTEGERLIST *next;
} INTEGERLIST;

typedef struct INTEGERLIST2
{
  int intDataSize,realDataSize;
  int *intData;
  double *realData;
  struct INTEGERLIST2 *next;
}INTEGERLIST2; 

#endif /* CODETYPES_H */
