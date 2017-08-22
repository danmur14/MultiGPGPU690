//this header file is used so that computeContourEdgesFor can call these functions located in the .cu file that are still cpu code
typedef float vec2[2];

int expectedEdgesKernel(float* vertexes, int nRows, int nCols, float level); //runs the expected edges kernel and does all the CUDA mallocs and memory management, returns int of expected edges
int actualEdgesKernel(float* vertexes, int nRows, int nCols, float level, int numExpectedPoints, vec2* buffer); //runs the actual edges kernel and does all the CUDA mallocs and memory management, returns int of actual edges present