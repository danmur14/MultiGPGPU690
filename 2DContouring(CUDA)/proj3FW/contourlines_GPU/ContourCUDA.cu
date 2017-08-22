#include <stdio.h>
#include <iostream>
#include <cmath>
#include "ContourCUDA.h"

__global__ void countEdges(float *vertexes, int nRows, int nCols, int *numExpectedPoints, float level)
{
  int x = (blockIdx.x * blockDim.x) + threadIdx.x; //found using the dimensions, based on the file given on the in-class materials page
  int y = (blockIdx.y * blockDim.y) + threadIdx.y;
  int index = y * gridDim.x * blockDim.x + x;

  int bound = ((nRows * nCols) - 1) - nCols; 
  if (index <= bound && ((index + 1) % nCols != 0)) //check if row below top row, and column before last column
  {
    //local values to determine how many edges to expect from the grid
    int count = 0;
    int nAbove = 0;
    int nBelow = 0;

    //each point to check
    float bottomL = vertexes[index];
    float bottomR = vertexes[index + 1];
    float topL = vertexes[index + nCols];
    float topR = vertexes[index + nCols + 1];

    //check if values are above or below the level, add accordingly
    if (bottomL > level)
    {
      nAbove++;
    }
    else
    {
      nBelow++;
    }
    if (bottomR > level)
    {
      nAbove++;
    }
    else
    {
      nBelow++;
    }
    if (topL > level)
    {
      nAbove++;
    }
    else
    {
      nBelow++;
    }
    if (topR > level)
    {
      nAbove++;
    }
    else
    {
      nBelow++;
    }

    //calculate number of expected edges based on how many vertices were below or above the desired level
    if (nAbove == 3 && nBelow == 1)
    {
      count = 1;
    }
    else if (nAbove == 1 && nBelow == 3)
    {
      count = 1;
    }
    else if (nAbove == 2 && nBelow == 2)
    {
      count = 2;
    }
    else
    {
      count = 0;
    }

    atomicAdd(numExpectedPoints, count); //add to the number of expected edges total
  }
}

__global__ void computeKernel(float *vertexes, int nRows, int nCols, int level, int *edgeCount, vec2 *actualEdgePoints, int *buf_location)
{
  int x = (blockIdx.x * blockDim.x) + threadIdx.x; 
  int y = (blockIdx.y * blockDim.y) + threadIdx.y;
  int index = y * gridDim.x * blockDim.x + x; //the index in the vertex array, acquired by multiplying the dimensions of each block and the threads within those blocks

  int bound = ((nRows * nCols) - 1) - nCols; //we do not want to check the top row of the grid or the index at the farthest right
  if (index <= bound && ((index % nCols) != nCols - 1)) //check if row below top row, and column before last column
  {
    //each point to check
    float bottomL = vertexes[index];
    float bottomR = vertexes[index + 1];
    float topL = vertexes[index + nCols];
    float topR = vertexes[index + nCols + 1];

    int loc; //the location of our index in the actualEdgePoints array so that we do not overlap edge points
    bool vertfound = false; //if we have found one vertex already
    int count = 0; //the number of vertexes we have found so far pertaining to one edge
    float x_coord = -1.0;
    float y_coord = -1.0;

    //check for missing data and return if found missing
    if (bottomL == -9999 || bottomR == -9999 || topL == -9999 || topR == -9999)
    {
      return; //do not check
    }

    //check every corner of the square starting from the bottom line, to right vertical, top horizontal, then left vertical
    if ((bottomL <= level && level <= bottomR) || (bottomL > level && level > bottomR)) //if the level is between the two points, not dependent on which corner is greater
    {
      if (bottomL <= level && level <= bottomR) //if the bottom right is greater
      {
        float f = (level - bottomL) / (bottomR - bottomL); //using the function given to find the coordinate between points
        x_coord = (1.0 - f) * (index % nCols) + f * ((index + 1) % nCols); //use that percentage and attribute it to x and y values, depending on which part of the square we are checking
        y_coord = (float)(index / nCols); //use the normal y coordinate
      }
      else if (bottomL > level && level > bottomR) //bottom left is greater, so the function is switched backwards
      {
        float f = (level - bottomR) / (bottomL - bottomR);
        x_coord = (1.0 - f) * ((index + 1) % nCols) + f * (index % nCols);
        y_coord = (float)(index / nCols);
      }
      if (!vertfound) //we have not found a vertice already, this is the first point of our edge
      {
        loc = atomicAdd(buf_location, 2); //get the index to add this vertex coordinate set to the actualEdgePoint array
        vertfound = true; //set to true so that we know we are on our second vertex of a certain edge
      }
      actualEdgePoints[loc + count][0] = x_coord; //set the coordinates of the vertex
      actualEdgePoints[loc + count][1] = y_coord;
      count++; //add to know how many vertices we have added so far
      if (count == 2) //checks if we have completed our edge with 2 vertices, reset the edge count
      {
        vertfound = 0;
        count = 0;
        atomicAdd(edgeCount, 1); //add to the total number of edges that we have
      }
    }

    //repeat
    if ((bottomL <= level && level <= topL) || (bottomL > level && level > topL))
    {
      if (bottomL <= level && level <= topL)
      {
        float f = (level - bottomL) / (topL - bottomL);
        x_coord = (float)(index % nCols);
        y_coord = (1.0 - f) * (index / nCols) + f * ((index + nCols) / nCols);
      }
      else if (bottomL > level && level > topL)
      {
        float f = (level - topL) / (bottomL - topL);
        x_coord = (float)(index % nCols);
        y_coord = (1.0 - f) * ((index + nCols) / nCols) + f * (index / nCols);
      }
      if (!vertfound)
      {
        loc = atomicAdd(buf_location, 2);
        vertfound = true;
      }
      actualEdgePoints[loc + count][0] = x_coord;
      actualEdgePoints[loc + count][1] = y_coord;
      count++;
      if (count == 2)
      {
        vertfound = 0;
        count = 0;
        atomicAdd(edgeCount, 1);
      }
    }

    if ((topR <= level && level <= topL) || (topR > level && level > topL))
    {
      if (topR <= level && level <= topL)
      {
        float f = (level - topR) / (topL - topR);
      x_coord = (1.0 - f) * ((index + nCols + 1) % nCols) + f * ((index + nCols) % nCols);
      y_coord = (float)((index + nCols) / nCols);
      }
      else if (topR > level && level > topL)
      {
        float f = (level - topL) / (topR - topL);
      x_coord = (1.0 - f) * ((index + nCols) % nCols) + f * ((index + nCols + 1) % nCols);
      y_coord = (float)((index + nCols) / nCols);
      }
      if (!vertfound)
      {
        loc = atomicAdd(buf_location, 2);
        vertfound = true;
      }
      actualEdgePoints[loc + count][0] = x_coord;
      actualEdgePoints[loc + count][1] = y_coord;
      count++;
      if (count == 2)
      {
        vertfound = 0;
        count = 0;
        atomicAdd(edgeCount, 1);
      }
    }

    if ((topR <= level && level <= bottomR) || (topR > level && level > bottomR))
    {
      if (topR <= level && level <= bottomR)
      {
        float f = (level - topR) / (bottomR - topR);
      x_coord = (float)((index + 1) % nCols);
      y_coord = (1.0 - f) * ((index + nCols + 1) / nCols) + f * ((index + 1) / nCols);
      }
      else if (topR > level && level > bottomR)
      {
         float f = (level - bottomR) / (topR - bottomR);
      x_coord = (float)((index + 1) % nCols);
      y_coord = (1.0 - f) * ((index + 1) / nCols) + f * ((index + nCols + 1) / nCols);
      }

      if (!vertfound)
      {
        loc = atomicAdd(buf_location, 2);
        vertfound = true;
      }
      actualEdgePoints[loc + count][0] = x_coord;
      actualEdgePoints[loc + count][1] = y_coord;
      count++;
      if (count == 2)
      {
        vertfound = 0;
        count = 0;
        atomicAdd(edgeCount, 1);
      }
    }
  }
}

int expectedEdgesKernel(float *vertexes, int nRows, int nCols, float level)
{
  float *dev_varray;                               //device vertex array buffer to copy
  int vert_size = (nRows * nCols) * sizeof(float); //size of vertex array to copy to gpu

  int *dev_count;          //expected edge device count variable to copy to gpu
  int zero = 0;            //start the device count at 0
  int *host_count = &zero; //host count to copy gpu value back to cpu

  cudaMalloc((void**)&dev_varray, vert_size);  //allocate size to hold the vertex array in gpu
  cudaMalloc((void**)&dev_count, sizeof(int)); //allocate one int variable on gpu to hold edge count

  cudaMemcpy(dev_varray, vertexes, vert_size, cudaMemcpyHostToDevice);    //copy vertexValues to the gpu in dev_varray
  cudaMemcpy(dev_count, host_count, sizeof(int), cudaMemcpyHostToDevice); //copy edge count to gpu starting at 0

  dim3 block(16, 16);                                                          //placeholder size for blocks only optimized for warps
  dim3 grid((nRows + block.x - 1) / block.x, (nCols + block.y - 1) / block.y); //launch grid based on size of vertexValues divided by block thread size

  countEdges<<<grid, block>>>(dev_varray, nRows, nCols, dev_count, level); //call kernel to count expected edges
  cudaThreadSynchronize();                                                 //barrier

  cudaMemcpy(host_count, dev_count, sizeof(int), cudaMemcpyDeviceToHost); //copy device count back to host count to pass back

  cudaFree(dev_varray); //free gpu vertex array
  cudaFree(dev_count);  //free device count

  return *host_count;
}

int actualEdgesKernel(float *vertexes, int nRows, int nCols, float level, int numExpectedPoints, vec2 *buffer)
{
  float *dev_varray;                               //device vertex array buffer to copy
  int vert_size = (nRows * nCols) * sizeof(float); //size of vertex array to copy to gpu

  int *dev_count;          //actual edges device count variable to copy to gpu
  int zero = 0;            //start the device count at 0
  int *host_count = &zero; //host count to copy gpu value back to cpu
  int *buf_location;       //index of the buffer that the coordinates should be placed at so edges are correct

  cudaMalloc(&dev_varray, vert_size);     //allocate size to hold the vertex array in gpu
  cudaMalloc(&dev_count, sizeof(int));    //allocate one int variable on gpu to hold actual edge count
  cudaMalloc(&buf_location, sizeof(int)); //allocate index of buffer we are writing coordinates to

  vec2 *dev_buffer;                                               //allocate buffer to hold points of actual edges calculated
  cudaMalloc(&dev_buffer, 2 * numExpectedPoints * sizeof(float)); //two points for each edge calculated

  cudaMemcpy(dev_varray, vertexes, vert_size, cudaMemcpyHostToDevice);       //copy vertexValues to the gpu in dev_varray
  cudaMemcpy(dev_count, host_count, sizeof(int), cudaMemcpyHostToDevice);    //copy edge count to gpu starting at 0
  cudaMemcpy(buf_location, host_count, sizeof(int), cudaMemcpyHostToDevice); //copy buffer index location to gpu
  cudaMemcpy(dev_buffer, buffer, 2 * numExpectedPoints * sizeof(float), cudaMemcpyHostToDevice);

  dim3 block(16, 16);                                                          //placeholder size for blocks only optimized for warps
  dim3 grid((nRows + block.x - 1) / block.x, (nCols + block.y - 1) / block.y); //launch grid based on size of vertexValues divided by block thread size

  computeKernel<<<grid, block>>>(dev_varray, nRows, nCols, level, dev_count, dev_buffer, buf_location); //compute actual number of edges in vertex array
  cudaThreadSynchronize();                                                                              //barrier

  cudaMemcpy(host_count, dev_count, sizeof(int), cudaMemcpyDeviceToHost);                        //copy back actual number of edges calculated
  cudaMemcpy(buffer, dev_buffer, 2 * numExpectedPoints * sizeof(float), cudaMemcpyDeviceToHost); //copy the actual edges from the gpu to the actual_edge_buffer on the cpu to then fill lines

  cudaFree(dev_varray); //free gpu vertex array
  cudaFree(dev_count);  //free device count
  cudaFree(dev_buffer); //free gpu actual edge buffer

  return *host_count;
}
