// FORTRAN names called from C use "_"
#define __UNDERSCORE__
// ILUPACK drivers use by default 64 bit long int as "integer"
//#define _LONG_INTEGER_
// we only use the real double precision version as "FLOAT"
#define _DOUBLE_REAL_

#include "ilupack.h"
#include "ilupackmacros.h"
#include "blas.h"
#include "ilupack/include/ilupack.h" /* has to be FIRST to avoid warnings about _Complex */
#include "pcildl.h"
#include <petsc/private/pcimpl.h>   /*I "petscpc.h" I*/

// elbow factor for the expected relative fill-in computed by the incomplete LDL^T factorization
#define ELBOW    3.0

/* define this to do some extra timing 
#define PCILDL_MORE_INFO
*/

/* define this to print the fill factors (on each rank)
#define PCILDL_PRINT_FILL
*/
#define PCILDL_PRINT_FILL

const char* const PCILDLOrderingTypes[] = {"METISN","METISE","RCM","AMD","ILDLOrderingType","ILDL_ORDERING_",0}; 

typedef struct {
  integer          *ia,*ja,*p,*invq,*jlu,nB;
  FLOAT            *pcolscale, *prowscale, *alu;
  PetscBool        matching;
  PetscReal        droptol;
  ILDLOrderingType ordering;
} PC_ILDL;

static PetscErrorCode ILDLSetUp(PC_ILDL *ildl,PetscInt nA,PetscScalar val[],PetscInt row_ptr[],PetscInt col_ind[])
{
  // ILDL drivers are written in FORTRAN 77, indexing is done in FORTRAN style (1,...,nA) and
  // symmetric matrices only consist of half of the matrix, say upper triangular part
  integer       i,j,k,l,nnz=0, *ia, *ja, *p, *invq, nB, ierr, *ibuff, *jlu, posiwk, PILUCparam;
  FLOAT         *a, *pcolscale, *prowscale, *dbuff, *alu;
  LONG_INT      imem, myimem;
  size_t        mem;
  DILUPACKparam options;
  Dmat          B;
#ifdef PCILDL_MORE_INFO
  double wtime;
  PetscMPIInt rank;
#endif

  PetscFunctionBegin;

  /* We assume that these types are the same. We have lazily left code referring to all types, but have removed conversions
     that would be required to use different types */
  if (sizeof(FLOAT) != sizeof(PetscScalar)) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"ILDL must use the same floating point type as PETSc!");
  if (sizeof(integer) != sizeof(PetscInt)) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"ILDL must use the same integer type as PETSc!");
  
  nB = (integer)nA;

#ifdef PCILDL_MORE_INFO
  wtime = MPI_Wtime(); 
  MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
#endif
  
  // logical pointer array, initialized with 0
  ia=(integer *)calloc((nB+1),sizeof(integer));
  
  // count number of nonzeros in the upper triangular part
  for (i=0; i<nB; i++) {
	  for (j=row_ptr[i]; j<row_ptr[i+1]; j++) {
      k=col_ind[j];
      // upper triangular part only
      if (k>=i) {
        // count nnz by row, do not use ia[0] (see below, why)
        ia[i+1]++;
        nnz++;
      } // end if
	  } // end for j
  } // end for i
  // switch from nz to pointer. Now this refers to the physical beginning of each (empty) row
  // this is why we used ia[1],...,ia[nB] previously and left ia[0] untouched
  for (i=0; i<nB; i++)
	  ia[i+1]+=ia[i];
  // array of column indices
  ja=(integer *)malloc(nnz*sizeof(integer));
  // array of numerical values
  a =(FLOAT *)  malloc(nnz*sizeof(FLOAT));
  // extract upper triangular part
  for (i=0; i<nB; i++) {
	  for (j=row_ptr[i]; j<row_ptr[i+1]; j++) {
      k=col_ind[j];
      // upper triangular part only
      if (k>=i) {
        // current start of the unoccupied part of row i
        l=ia[i];
        // FORTRAN starts indexing with 1
        ja[l]=k+1;
        // copy numerical value
        a[l]=val[j];
        // free part of row i incremented
        ia[i]=l+1;
      } // end if
	  } // end for j
  } // end for i
  // now the entries ia[0],...,ia[nB-1] of the pointer array must have those values
  // that were previously stored as pointers in ia[1],...,ia[nB]
  // shift pointers back, FORTRAN starts indexing with 1
  for (i=nB; i>0; i--)
	  ia[i]=ia[i-1]+1;
  ia[0]=1;
  
  
  
  // now ia,ja,a contain the upper triangular part of the matrix in FORTRAN format
  B.nr = B.nc=nB;
  B.nnz = nnz;
  B.ia = ia;
  B.ja = ja;
  B.a = a;

#ifdef PCILDL_MORE_INFO
  wtime = MPI_Wtime() - wtime;
  ierr = PetscPrintf(PETSC_COMM_SELF,"[%D] PCILDL setup phase 1 : %g s\n",rank,wtime);
  wtime = MPI_Wtime();
#endif


  // initialize ILUPACK options structure to its default options
  DSYMAMGinit(&B,&options);

  // threshold for ILU
  options.droptol = (FLOAT) ildl->droptol;

  // turn on matching
  options.matching = ildl->matching ? 1 : 0;

  switch (ildl->ordering) {
    case ILDL_ORDERING_METISE:
      options.ordering = "metise";
      break;
    case ILDL_ORDERING_METISN:
      options.ordering = "metisn";
      break;
    case ILDL_ORDERING_RCM:
      options.ordering = "rcm";
      break;
    case ILDL_ORDERING_AMD:
      options.ordering = "amd";
      break;
    default :
      SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unrecognized ILDL ordering type");
  }
  
  // -----------------------------------------------------------------
  // ----- preprocessing part: scaling and reordering the system -----
  // memory for scaling and permutation
  p   =(integer *)MAlloc((size_t)B.nc*sizeof(integer),"CSRInterface:p");
  invq=(integer *)MAlloc((size_t)B.nc*sizeof(integer),"CSRInterface:invq");
  pcolscale=(FLOAT *)MAlloc((size_t)B.nc*sizeof(FLOAT),"CSRInterface:pcolscale");
  prowscale=pcolscale;
  // matching turned on?
  if (options.matching) {
    if (!strcmp("metise",options.ordering))
	    ierr=DSYMperm_mc64_metis_e(B, prowscale, pcolscale, p, invq, &nB, &options);
    else if (!strcmp("amd",options.ordering))
	    ierr=DSYMperm_mc64_amd(B, prowscale, pcolscale, p, invq, &nB, &options);
    else if (!strcmp("rcm",options.ordering))
	    ierr=DSYMperm_mc64_rcm(B, prowscale, pcolscale, p, invq, &nB, &options);
    else // METIS nested dissection by nodes used otherwise
	    ierr=DSYMperm_mc64_metis_n(B, prowscale, pcolscale, p, invq, &nB, &options);
  }
  else { // no matching
    if (!strcmp("metisn",options.ordering))
	    ierr=DSYMperm_metis_n(B, prowscale, pcolscale, p, invq, &nB, &options);
    else if (!strcmp("metise",options.ordering))
	    ierr=DSYMperm_metis_e(B, prowscale, pcolscale, p, invq, &nB, &options);
    else if (!strcmp("amd",options.ordering))
	    ierr=DSYMperm_amd(B, prowscale, pcolscale, p, invq, &nB, &options);
    else if (!strcmp("rcm",options.ordering))
	    ierr=DSYMperm_rcm(B, prowscale, pcolscale, p, invq, &nB, &options);
    else  // none
	    ierr=DSYMperm_null(B, prowscale, pcolscale, p, invq, &nB, &options);
  }
  // -------------------- END preprocessing part ---------------------
  // -----------------------------------------------------------------

  /* Free some data that we won't use anymore */
  FREE(options.dbuff);
  FREE(options.ibuff);
  
#ifdef PCILDL_MORE_INFO
  wtime = MPI_Wtime() - wtime;
  ierr = PetscPrintf(PETSC_COMM_SELF,"[%D] PCILDL setup phase 2 : %g s\n",rank,wtime);
  wtime = MPI_Wtime();
#endif
  
  // -----------------------------------------------------------------
  // ------------ compute incomplete LDL^T factorization -------------
  PILUCparam=0;
  // bit 0: simple dual threshold dropping strategy
  // bit 1: a zero pivot at step k is replaced by a small number
  PILUCparam|=2;
  // bit 2: simple Schur complement
  // bit 3: ILU is computed for the first time
  // bit 4: simple inverse-based dropping would have been used if it applies
  // bit 5: diagonal compensation is possibly done
  
  // auxiliary buffers
  ibuff=(integer *)MAlloc((size_t)14*B.nc*sizeof(integer),"CSRInterface:ibuff");
  dbuff=(FLOAT *)  MAlloc((size_t)4*B.nc*sizeof(FLOAT),"CSRInterface:dbuff");
  
  // ILU buffers
  // just a simple guess at least 2nB+1!
  mem=ELBOW*B.ia[B.nr]+1;
  jlu=(integer *)MAlloc((size_t)mem*sizeof(integer),"DSYMildlfactor");
  alu=(FLOAT *)  MAlloc((size_t)mem*sizeof(FLOAT),"DSYMildlfactor");
  
  // since FORTRAN77 does not have memory allocation we need to provide
  // enough memory in advance (mem is a guess). If it turns out that the memory
  // is not sufficient, then the FORTRAN77 routine DSYMiluc will stop, the
  // external C-wrapper has to do a "realloc" and then calls the FORTRAN77 code
  // again. Note that DSYMiluc will not restart from scratch but return to
  // the position where it stopped previously (reverse communincation, static
  // variables) and continue the incomplete LDL^T factorization
  posiwk=0;
  ierr=0;
  do {
    imem=(LONG_INT)mem;
    DSYMiluc(&(B.nc),B.a,B.ja,B.ia,
             &(options.lfil),&(options.droptol),&PILUCparam,
             p,invq,alu,jlu,&imem,
             dbuff,ibuff,&posiwk,&ierr);
    
    // not enough memory for ILDL?
    if (posiwk>0) {
	    // total amount of memory requested by the parameters
	    myimem=ELBOW*(size_t)B.ia[B.nr];
	    mem+=myimem;
	    jlu=(integer *)ReAlloc(jlu,mem*sizeof(integer),"CSRInterface:jlu");
	    alu=(FLOAT *)  ReAlloc(alu,mem*sizeof(FLOAT),  "CSRInterface:alu");
    } // end if
  } while (posiwk>0);
  if (ierr) {
    printf("ILDL terminated with error code %d\n",ierr);
    exit (ierr);
  }
  free(ibuff);
  free(dbuff);

#ifdef PCILDL_MORE_INFO
  wtime = MPI_Wtime() - wtime;
  ierr = PetscPrintf(PETSC_COMM_SELF,"[%D] PCILDL setup phase 3 : %g s\n",rank,wtime);
  wtime = MPI_Wtime();
#endif

#ifdef PCILDL_PRINT_FILL
  ierr = PetscPrintf(PETSC_COMM_SELF,"relative fill ILDL/A: %8.1le (wrt %D nz)\n",((double)jlu[nB])/B.ia[nB],B.ia[nB]);
#endif
  // ---------------------- END incomplete LDL^T ---------------------
  // -----------------------------------------------------------------

  free(a);
  free(ia);
  free(ja);

  /* Retain this data to apply the factors */
  ildl->p = p;
  ildl->invq = invq;
  ildl->jlu = jlu;
  ildl->nB = nB;
  ildl->pcolscale = pcolscale;
  ildl->prowscale = prowscale;
  ildl->alu = alu;

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCSetUp_ILDL"
static PetscErrorCode PCSetUp_ILDL(PC pc)
{
  PC_ILDL      *ildl = (PC_ILDL*)pc->data;
  PetscErrorCode ierr;
  Mat A;
  PetscBool issame = PETSC_FALSE;
  
  PetscFunctionBegin;
  PCGetOperators(pc,&A,NULL);
  PetscObjectTypeCompare((PetscObject)A,MATSEQAIJ,&issame);
  if (issame) {
    PetscInt  nA;
    const PetscInt  *col_ind;
    const PetscInt  *row_ptr;
    PetscScalar     *val;
    const PetscBool shift=PETSC_FALSE, symmetric=PETSC_FALSE, inodecompressed=PETSC_FALSE;
    PetscBool       done;
    
    /* Get Raw CSR data */
    ierr = MatSeqAIJGetArray(A,&val);CHKERRQ(ierr);
    ierr = MatGetRowIJ(A,shift,symmetric,inodecompressed,&nA,&row_ptr,&col_ind,&done);CHKERRQ(ierr);

    /* ILDL setup */
    ierr = ILDLSetUp(ildl,nA,(PetscScalar*)val,(PetscInt*)row_ptr,(PetscInt*)col_ind);CHKERRQ(ierr);
    
    ierr = MatRestoreRowIJ(A,shift,symmetric,inodecompressed,&nA,&row_ptr,&col_ind,&done);CHKERRQ(ierr);
    ierr = MatSeqAIJRestoreArray(A,&val);CHKERRQ(ierr);
    
  } else SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Only valid for MatSEQAIJ");

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCApply_ILDL"
static PetscErrorCode PCApply_ILDL(PC pc,Vec vecB,Vec vecX)
{
  PC_ILDL           *ildl = (PC_ILDL*)pc->data;
  PetscInt          i;
  const PetscScalar *_vecB;
  PetscScalar       *_vecX;
  FLOAT             *x,*b;
  PetscErrorCode    ierr;
  
  PetscFunctionBegin;
  ierr = VecGetArrayRead(vecB,&_vecB);CHKERRQ(ierr);
  ierr = VecGetArray(vecX,&_vecX);CHKERRQ(ierr);
  x = _vecX;
  b = (FLOAT*) malloc(ildl->nB * sizeof(FLOAT));
  for (i=0; i<ildl->nB; i++) {
	  b[i] = (FLOAT)_vecB[i];
  }
  ierr = VecRestoreArrayRead(vecB,&_vecB);CHKERRQ(ierr);

  /* row scaling */
  for (i=0; i<ildl->nB; i++) {
	  b[i] *= ildl->prowscale[i];
  }
  /* permutation p (FORTRAN indexing) */
  for (i=0; i<ildl->nB; i++) {
	  x[i] = b[ ildl->p[i] - 1 ];
  }
  
  /* forward/backward solve LDL^Tb=x
     This accounts for almost all of the time for this routine */
  integer nB = ildl->nB;
  DSYMpilucsol(&nB, x, b, ildl->alu, ildl->jlu);
  
  /* inverse permutation invq=p^{-1} (FORTRAN indexing) */
  for (i=0; i<ildl->nB; i++) {
	  x[i] = b[ ildl->invq[i] - 1 ];
  }

  /* column scaling */
  for (i=0; i<ildl->nB; i++) {
	  x[i] *= ildl->pcolscale[i];
  }
  
  ierr = VecRestoreArray(vecX,&_vecX);CHKERRQ(ierr);
  free(b);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCReset_ILDL"
static PetscErrorCode PCReset_ILDL(PC pc)
{
  PC_ILDL      *ildl = (PC_ILDL*)pc->data;

  PetscFunctionBegin;

  if (ildl->p) { free(ildl->p); ildl->p = NULL; }
  if (ildl->invq) { free(ildl->invq); ildl->invq = NULL; }
  if (ildl->jlu) { free(ildl->jlu); ildl->jlu = NULL; }

  if (ildl->pcolscale == ildl->prowscale) {
    if (ildl->pcolscale) { free(ildl->pcolscale); ildl->pcolscale = NULL; ildl->prowscale = NULL; }
  } else {
    if (ildl->pcolscale) { free(ildl->pcolscale); ildl->pcolscale = NULL; }
    if (ildl->prowscale) { free(ildl->prowscale); ildl->prowscale = NULL; }
  }
  if (ildl->alu) { free(ildl->alu); ildl->alu = NULL; }

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCDestroy_ILDL"
static PetscErrorCode PCDestroy_ILDL(PC pc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PCReset_ILDL(pc);CHKERRQ(ierr);
  ierr = PetscFree(pc->data);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCSetFromOptions_ILDL"
static PetscErrorCode PCSetFromOptions_ILDL(PetscOptionItems *PetscOptionsObject,PC pc)
{
  PC_ILDL      *ildl = (PC_ILDL*)pc->data;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"ILDL options");CHKERRQ(ierr);
  ierr = PetscOptionsBool("-pc_ildl_matching","Use matching",NULL,ildl->matching,&ildl->matching,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-pc_ildl_droptol","Set drop tol",NULL,ildl->droptol,&ildl->droptol,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnum("-pc_ildl_ordering","Choose ordering (rcm, amd, metise, metisn)",NULL,PCILDLOrderingTypes,(PetscEnum)ildl->ordering,(PetscEnum*)&ildl->ordering,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCView_ILDL"
static PetscErrorCode PCView_ILDL(PC pc,PetscViewer viewer)
{
  PC_ILDL        *ildl = (PC_ILDL*)pc->data;
  PetscErrorCode ierr;
  PetscBool      iascii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii);CHKERRQ(ierr);
  if (iascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"  ILDL: matching : %d\n",ildl->matching);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  ILDL: droptol : %g\n",ildl->droptol);CHKERRQ(ierr);
    switch (ildl->ordering) {
      case ILDL_ORDERING_METISE:
        ierr = PetscViewerASCIIPrintf(viewer,"  ILDL: ordering : metise\n");CHKERRQ(ierr);
        break;
      case ILDL_ORDERING_METISN:
        ierr = PetscViewerASCIIPrintf(viewer,"  ILDL: ordering : metisn\n");CHKERRQ(ierr);
        break;
      case ILDL_ORDERING_RCM:
        ierr = PetscViewerASCIIPrintf(viewer,"  ILDL: ordering : rcm\n");CHKERRQ(ierr);
        break;
      case ILDL_ORDERING_AMD:
        ierr = PetscViewerASCIIPrintf(viewer,"  ILDL: ordering : amd\n");CHKERRQ(ierr);
        break;
      default :
        SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"Unrecognized ILDL ordering type");
    }
  } 
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCCreate_ILDL"
PETSC_EXTERN PetscErrorCode PCCreate_ILDL(PC pc)
{
  PC_ILDL      *ildl;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr     = PetscNewLog(pc,&ildl);CHKERRQ(ierr);
  pc->data = (void*)ildl;

  pc->ops->apply               = PCApply_ILDL;
  pc->ops->applytranspose      = 0;
  pc->ops->setup               = PCSetUp_ILDL;
  pc->ops->reset               = PCReset_ILDL;
  pc->ops->destroy             = PCDestroy_ILDL;
  pc->ops->setfromoptions      = PCSetFromOptions_ILDL;
  pc->ops->view                = PCView_ILDL;
  pc->ops->applyrichardson     = 0;
  pc->ops->applysymmetricleft  = 0;
  pc->ops->applysymmetricright = 0;

  ildl->matching  = PETSC_TRUE;
  ildl->droptol   = 0.01;
  ildl->ordering = ILDL_ORDERING_METISN;

  PetscFunctionReturn(0);
}

#undef PCILDL_MORE_INFO
#undef PCILDL_PRINT_FILL

static PetscBool PCILDLPackageInitialized = PETSC_FALSE;

/* Either call this function from your code, or load the shared library
 * dynamically by adding "-dll_append /path/to/libpcildl.so" to your
 * command line argument */
#undef __FUNCT__
#define __FUNCT__ "PCILDLInitializePackage"
PetscErrorCode PCILDLInitializePackage(void)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (PCILDLPackageInitialized) PetscFunctionReturn(0);
  PCILDLPackageInitialized = PETSC_TRUE;
  ierr = PCRegister(PCILDL,PCCreate_ILDL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#if defined(PETSC_HAVE_DYNAMIC_LIBRARIES)
#undef __FUNCT__
#define __FUNCT__ "PetscDLLibraryRegister_pcildl"
PETSC_EXTERN PetscErrorCode PetscDLLibraryRegister_pcildl(void)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = PCILDLInitializePackage();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#endif
