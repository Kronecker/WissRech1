//
// Created by grabiger on 05.12.2017.
//

#include "memAlignOS.h"
#include <stdlib.h>
//
// Created by grabiger on 05.12.2017.
//


#ifdef _WIN32
void memAlignOS(void** ptr, size_t alignment, size_t size) {
    *ptr=_aligned_malloc(size,alignment);
}
void freeAlignedMemOS(void* ptr) {
    _aligned_free(ptr);
}

#else
void memAlignOS(void** ptr, size_t alignment, size_t size) {
    posix_memalign(ptr, alignment, size);
}
void freeAlignedMemOS(void* ptr) {
    free(ptr);
}



#endif













































