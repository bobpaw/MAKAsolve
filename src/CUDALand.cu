#include "CUDALand.h"
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

int getDeviceCount() {
	int count = 0;
	cudaError_t err = cudaGetDeviceCount(&count);
	if (err != cudaSuccess) {
		fprintf(stderr, "Failed to get device count: %s\n",
						cudaGetErrorString(err));
		return 0;
	}
	return count;
}

} // namespace cudaland