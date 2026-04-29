#ifndef CUDALAND_H
#define CUDALAND_H

// On AiMOS, CUDA functions only appear to work when they are compiled in a .cu
// file and with only nvcc This is compiled as a separate library and linked to
// MAKAsolve (like assignment 4)
namespace cudaland {

void setDevice(int dev);

}

#endif // CUDALAND_H