#include "ilupack.h" /* has to be FIRST to avoid warnings about _Complex */
#include "pcilupack.h"
#include <petsc/private/pcimpl.h>   /*I "petscpc.h" I*/

/* Hard-coded - this will only work in double */
#define FLOAT double

/*
#define PRINT_INFO
*/
#define PCILUPACK_FILL_INFO

typedef struct {
  ILUPACKparam  param;
  SPARSEmat     B; 
  AMGlevelmat   PRE;

  integer       *ia, *ja;
  FLOAT         *a;

  PetscReal     droptol;
  PetscReal     droptolS;
  PetscReal     condest;
   
  PetscBool     ilupackinit;
} PC_ILUPACK;

/* A local function with the ILUpack setup calls */
static PetscErrorCode ILUPACKSetUp(PC_ILUPACK *ilupack,PetscInt nA,PetscScalar val[],PetscInt row_ptr[],PetscInt col_ind[])
{
  integer       i,j,k,l,nnz=0;
  int           ilupackerr;

  integer       *ia,*ja,nB;
  FLOAT         *a; /* convenience */

  PetscFunctionBegin;

  /* We assume that these types are the same. We have lazily left code referring to all types, but have removed conversions
     that would be required to use different types */
  if (sizeof(FLOAT) != sizeof(PetscScalar)) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"ILUPACK must use the same floating point type as PETSc!");
  if (sizeof(double) != sizeof(PetscScalar)) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"PCILUPACK only supports double for now!");
  if (sizeof(integer) != sizeof(PetscInt)) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_SUP,"ILUPACK must use the same integer type as PETSc!");
  
  nB = (integer)nA;

#ifdef PRINT_INFO

  for (i=0; i<nB; i++) {
    printf("%d : ",i);
	  for (j=row_ptr[i]; j<row_ptr[i+1]; j++) {
      k=col_ind[j];
        printf("(%d,%f) ",k,val[j]);
    }
    printf("\n");
  }
#endif
  
  /* logical pointer array, initialized with 0 */
  ilupack->ia=(integer *)calloc((nB+1),sizeof(integer));
  ia=ilupack->ia;
  
  /* count number of nonzeros in the upper triangular part */
  for (i=0; i<nB; i++) {
	  for (j=row_ptr[i]; j<row_ptr[i+1]; j++) {
      k=col_ind[j];
      /* upper triangular part only */
      if (k>=i) {
        /* count nnz by row, do not use ia[0] (see below, why) */
        ia[i+1]++;
        nnz++;
      }
	  }
  } 
  /* switch from nz to pointer. Now this refers to the physical beginning of each (empty) row */
  /* this is why we used ia[1],...,ia[nB] previously and left ia[0] untouched */
  for (i=0; i<nB; i++)
	  ia[i+1]+=ia[i];
  /* array of column indices */
  ilupack->ja=(integer *)malloc(nnz*sizeof(integer));
  ja=ilupack->ja;
  /* array of numerical values */
  ilupack->a =(FLOAT *)  malloc(nnz*sizeof(FLOAT));
  a=ilupack->a;
  /* extract upper triangular part */
  for (i=0; i<nB; i++) {
	  for (j=row_ptr[i]; j<row_ptr[i+1]; j++) {
      k=col_ind[j];
      /* upper triangular part only */
      if (k>=i) {
        /* current start of the unoccupied part of row i */
        l=ia[i];
        /* FORTRAN starts indexing with 1 */
        ja[l]=k+1;
        /* copy numerical value */
        a[l]=val[j];
        /* free part of row i incremented */
        ia[i]=l+1;
      }
	  }
  }
  /* now the entries ia[0],...,ia[nB-1] of the pointer array must have those values */
  /* that were previously stored as pointers in ia[1],...,ia[nB] */
  /* shift pointers back, FORTRAN starts indexing with 1 */
  for (i=nB; i>0; i--) ia[i]=ia[i-1]+1;
  ia[0]=1;
  
  /* now ia,ja,a contain the upper triangular part of the matrix in FORTRAN format */

  /* Hard-coded! */
  ilupack->B.isreal=1;      /* real matrix */
  ilupack->B.issymmetric=1; /* symmetric matrix */
  ilupack->B.isdefinite=0;  /* not positive definite */
  ilupack->B.issingle=0;    /* use double precision */

  ilupack->B.nr = ilupack->B.nc=nB;
  ilupack->B.nnz = nnz;
  ilupack->B.ia = ilupack->ia;
  ilupack->B.ja = ilupack->ja;
  ilupack->B.a = (void*) ilupack->a;

#ifdef PRINT_INFO
  printf("nr : %d\n",ilupack->B.nr);
  printf("nc : %d\n",ilupack->B.nc);
  printf("nnz : %d\n",ilupack->B.nnz);
  printf("ia : ");
  for (i=0; i<ilupack->B.nr+1; i++){
    printf("%d ",ilupack->B.ia[i]);
  }
  printf("\n");
  printf("ja : ");
  for (i=0; i<ilupack->B.nnz; i++){
    printf("%d ",ilupack->B.ja[i]);
  }
  printf("\n");
  printf("a : ");
  for (i=0; i<ilupack->B.nnz; i++){
    printf("%f ",((FLOAT*)ilupack->B.a)[i]);
  }
  printf("\n");
#endif

  /* initialize ILUPACK options structure to its default options */
  AMGinit(&ilupack->B,&ilupack->param);

  /* thresholds for ILU */
  ilupack->param.droptol = (FLOAT) ilupack->droptol;
  ilupack->param.droptolS= (FLOAT) ilupack->droptolS;

  /* turn on matching */
  ilupack->param.matching = 1;

  /* METIS-by-nodes ordering */
  ilupack->param.ordering = "metisn";

  /* Parameter which controls how many levels. */
  ilupack->param.condest=(FLOAT) ilupack->condest;

  /* Call ILUPack Setup routine */
#if 1
  ilupackerr=AMGfactor(&ilupack->B, &ilupack->PRE, &ilupack->param);
  if (ilupackerr) {
    SETERRQ1(PETSC_COMM_WORLD,PETSC_ERR_LIB,"ILUPACK encountered an error during factorization: %d",ilupackerr);
  } else {
#ifdef PCILUPACK_FILL_INFO 
    PetscErrorCode ierr;
    ierr=PetscPrintf(PETSC_COMM_WORLD,"factorization successful with %d levels completed\n", ilupack->PRE.nlev);CHKERRQ(ierr);
    ierr=PetscPrintf(PETSC_COMM_WORLD,"final elbow space factor=%8.2f\n",ilupack->param.elbow+0.005);CHKERRQ(ierr);
#endif
  }
  ilupack->ilupackinit=PETSC_TRUE;
#endif

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCSetUp_ILUPACK"
static PetscErrorCode PCSetUp_ILUPACK(PC pc)
{
  PC_ILUPACK     *ilupack = (PC_ILUPACK*)pc->data;
  PetscErrorCode ierr;
  Mat            A;
  PetscBool      issame;
  
  PetscFunctionBegin;
  ierr = PCGetOperators(pc,&A,NULL);CHKERRQ(ierr);
  ierr = PetscObjectTypeCompare((PetscObject)A,MATSEQAIJ,&issame);CHKERRQ(ierr);
  
  if (issame) {
    PetscInt          nA;
    const PetscInt    *col_ind,*row_ptr;
    PetscScalar       *val;
    const PetscBool   shift=PETSC_FALSE,symmetric=PETSC_FALSE,inodecompressed=PETSC_FALSE;
    PetscBool         done;
    
    /* Get Raw CSR data (Note that we are being wasteful and not using sym types) */
    ierr = MatSeqAIJGetArray(A,&val);CHKERRQ(ierr); /* PETSc doesn't have a "Read" version of this, hence val isn't marked const */
    ierr = MatGetRowIJ(A,shift,symmetric,inodecompressed,&nA,&row_ptr,&col_ind,&done);CHKERRQ(ierr);

    /* ILUPACK setup (local function) */
    ierr = ILUPACKSetUp(ilupack,nA,(PetscScalar*)val,(PetscInt*)row_ptr,(PetscInt*)col_ind);CHKERRQ(ierr);
    
    ierr = MatRestoreRowIJ(A,shift,symmetric,inodecompressed,&nA,&row_ptr,&col_ind,&done);CHKERRQ(ierr);
    ierr = MatSeqAIJRestoreArray(A,&val);CHKERRQ(ierr);
    
  } else SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_USER,"Only valid for MATSEQAIJ"); /* Yes, this sucks - should be SBAIJ but more work */

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCApply_ILUPACK"
static PetscErrorCode PCApply_ILUPACK(PC pc,Vec vecB,Vec vecX)
{
  PC_ILUPACK        *ilupack = (PC_ILUPACK*)pc->data;
  PetscInt          nB;
  FLOAT             *x,*b;
  PetscErrorCode    ierr;
  
  PetscFunctionBegin;
  
  ierr = VecGetSize(vecB,&nB);CHKERRQ(ierr); /* note this is uniprocessor only */
  ierr = VecGetArrayRead(vecB,&b);CHKERRQ(ierr); 
  ierr = VecGetArray(vecX,&x);CHKERRQ(ierr);

  AMGsol(&ilupack->PRE,&ilupack->param,b,x);
  
  ierr = VecRestoreArray(vecX,&x);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(vecB,&b);CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCReset_ILUPACK"
static PetscErrorCode PCReset_ILUPACK(PC pc)
{
  PC_ILUPACK      *ilupack = (PC_ILUPACK*)pc->data;

  PetscFunctionBegin;
  if (ilupack->ilupackinit){
    ilupack->ilupackinit=PETSC_FALSE;
    AMGdelete(&ilupack->B,&ilupack->PRE,&ilupack->param); 
  }
  if(ilupack->ia) {free(ilupack->ia); ilupack->ia=NULL;}
  if(ilupack->ja) {free(ilupack->ja); ilupack->ja=NULL;}
  if(ilupack->a)  {free(ilupack->a);  ilupack->a=NULL;}
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCDestroy_ILUPACK"
static PetscErrorCode PCDestroy_ILUPACK(PC pc)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PCReset_ILUPACK(pc);CHKERRQ(ierr);
  ierr = PetscFree(pc->data);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCSetFromOptions_ILUPACK"
static PetscErrorCode PCSetFromOptions_ILUPACK(PetscOptionItems *PetscOptionsObject,PC pc)
{
  PC_ILUPACK      *ilupack = (PC_ILUPACK*)pc->data;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscOptionsHead(PetscOptionsObject,"ILUPACK options");CHKERRQ(ierr);
  ierr = PetscOptionsReal("-pc_ilupack_droptol","Set drop tol",NULL,ilupack->droptol,&ilupack->droptol,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-pc_ilupack_droptolS","Set Schur complement drop tol",NULL,ilupack->droptolS,&ilupack->droptolS,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsReal("-pc_ilupack_condest","Set condition number used to choose level structure",NULL,ilupack->condest,&ilupack->condest,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsTail();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCView_ILUPACK"
static PetscErrorCode PCView_ILUPACK(PC pc,PetscViewer viewer)
{
  PC_ILUPACK     *ilupack = (PC_ILUPACK*)pc->data;
  PetscErrorCode ierr;
  PetscBool      iascii;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii);CHKERRQ(ierr);
  if (iascii) {
    ierr = PetscViewerASCIIPrintf(viewer,"  ILUPACK: droptol  : %g\n",ilupack->droptol );CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  ILUPACK: droptolS : %g\n",ilupack->droptolS);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer,"  ILUPACK: condest  : %g\n",ilupack->condest );CHKERRQ(ierr);
  } 
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PCCreate_ILUPACK"
PETSC_EXTERN PetscErrorCode PCCreate_ILUPACK(PC pc)
{
  PC_ILUPACK     *ilupack;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr     = PetscNewLog(pc,&ilupack);CHKERRQ(ierr);
  pc->data = (void*)ilupack;

  pc->ops->apply               = PCApply_ILUPACK;
  pc->ops->applytranspose      = 0;
  pc->ops->setup               = PCSetUp_ILUPACK;
  pc->ops->reset               = PCReset_ILUPACK;
  pc->ops->destroy             = PCDestroy_ILUPACK;
  pc->ops->setfromoptions      = PCSetFromOptions_ILUPACK;
  pc->ops->view                = PCView_ILUPACK;
  pc->ops->applyrichardson     = 0;
  pc->ops->applysymmetricleft  = 0;
  pc->ops->applysymmetricright = 0;

  ilupack->ilupackinit = PETSC_FALSE;

  ilupack->droptol   = 0.01;
  ilupack->droptolS  = 0.01;
  ilupack->condest   = 10.0;
  PetscFunctionReturn(0);
}

static PetscBool PCILUPACKPackageInitialized = PETSC_FALSE;

/* Either call this function from your code, or load the shared library
 * dynamically by adding "-dll_append /path/to/libpcilupack.so" to your
 * command line argument */
#undef __FUNCT__
#define __FUNCT__ "PCILUPACKInitializePackage"
PetscErrorCode PCILUPACKInitializePackage(void)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (PCILUPACKPackageInitialized) PetscFunctionReturn(0);
  PCILUPACKPackageInitialized = PETSC_TRUE;
  ierr = PCRegister(PCILUPACK,PCCreate_ILUPACK);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#if defined(PETSC_HAVE_DYNAMIC_LIBRARIES)
#undef __FUNCT__
#define __FUNCT__ "PetscDLLibraryRegister_pcilupack"
PETSC_EXTERN PetscErrorCode PetscDLLibraryRegister_pcilupack(void)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = PCILUPACKInitializePackage();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#endif
