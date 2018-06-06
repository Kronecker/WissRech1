//
// Created by grabiger on 05.12.2017.
//

#ifndef WISSRECH1_MEMALIGNOS_H
#include <stdlib.h>
#define WISSRECH1_MEMALIGNOS_H
void memAlignOS(void** ptr,size_t alignment,  size_t size);
void freeAlignedMemOS(void* ptr);
#endif //WISSRECH1_MEMALIGNOS_H
