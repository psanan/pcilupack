#ifndef __PCILUPACK_h__
#define __PCILUPACK_h__

#include <petscpc.h>

#define PCILUPACK "ilupack"

PetscErrorCode PCCreate_ILUPACK(PC pc);

#endif
