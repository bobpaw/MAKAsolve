#include <MAKAsolve/CUDALand.h>
#include <cuda_runtime.h>
#include <stdio.h>

// On AiMOS, CUDA functions only appear to work when they are compiled in a .cu
// file and with only nvcc This is compiled as a separate library and linked to
// MAKAsolve (like assignment 4)
namespace cudaland {

void setDevice(int dev) {
	cudaError_t err = cudaSetDevice(dev);
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to set to GPU %d\n", dev, cudaGetErrorString(err));
	}
}

} // namespace cudaland