//
// ToyPlaneFDTDnaanan
//
// This is planar of 2D Finite-Difference Time-Domain coded in C with numerous comments
// as a working tool for instruction and further deeper discussion of time-domain local-
// operator codes.
//
// Copyright (C) 1999 Paul Robert Hayes
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// When reporting or displaying any results or animations created
// using this code or modification of this code, make the appropriate
// citation referencing ToyFDTD by name and including the version
// number.  You should know better. 8)
//   
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// http://cemtach.com
//


#include <stdio.h>
#include <float.h>
#include <math.h>


#include "toyPlane.h"





void
main(int argc, char **argv)
  {
  //
  // Current simulation time
  //
  double currentSimulationTime;
  //
  // The current iteration which will be used to track output data files
  //
  int currentIteration;
  //
  // Time step of the simulation
  //
  double deltaTime;
  //
  // Estimation variables which help locate an appropriate deltaTime
  //
  double deltaTimeEstimate;
  double velocity;
  //
  // Begin simulating as if this were time zero
  //
  double startSimulationTime;
  //
  // Simulate no further than the following
  //
  double stopSimulationTime;
  //
  // Plotting variables where plotTime is the next time at which plotting will occur
  //
  double plotTime;
  double plotDelta;
  //
  // Actual data arrays
  //
#ifdef E_NORMAL
  double **ez, *ezBase;
  double **hx, *hxBase;
  double **hy, *hyBase;
#endif
#ifdef H_NORMAL
  double **hz, *hzBase;
  double **ex, *exBase;
  double **ey, *eyBase;
#endif
  Material ***material, *materialBase;
  //
  // Size of the mesh in terms of the `out of the screen' quantity (ez or hz)
  //
  int nx, ny;
  //
  // Spatial step in each direction, which need not be equal values
  //
  double dx, dy;
  //
  // Actual positioning of the 2D mesh as a rectangle of space
  //
  double xLow, xHigh;
  double yLow, yHigh;
  //
  // material definition arrays
  //
  Material *freeSpace, *newMaterial, *copper;
  //
  // Sourcing variables
  //
  double frequency;
  //
  // Generic temporary variables commonly used in multiple places
  //
  FILE *singleFilePointer;
  char filename[1024];
  int i,j;
  double etmp;
  int iLow, iHigh;
  int jLow, jHigh;
#ifdef SINGLE_POINT
  int iSingle, jSingle;
#endif
  double eps, mu, one, two;
  double value;
  //
  // Temporary values for the non-compact higher order stencil calculations
  //
  double d2x, d2y, d4x, d4y;

  //
  // Print out some identifying information
  //
  fprintf(stdout, "\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "ToyAirportFDTDnaanan\n");
  fprintf(stdout, "Copyright (C) 1999 Paul Hayes\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "ToyAirportFDTDnaanan is free software published under the terms\n"); 
  fprintf(stdout, "of the GNU General Public License as published by the\n"); 
  fprintf(stdout, "Free Software Foundation.\n");  
  fprintf(stdout, "\n");



  ///////////////////////////////////////////////////////////////////////////////
  // Section one of problem definition where initial values will
  // determine the size of the mesh elements.  A second section
  // will be required after memory allocation in order to fill the mesh
  // with relevant values.
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  //
  //
  // Establish physical limits for the 2D mesh based off a WG-16 waveguide with
  // the larger dimension as 0.90 inch and the smaller dimension as 0.40 inch.
  //
  // I chose these values on a whim, the Low values could start at 0.0 and the High
  // quantities could be 0.90*0.0254, for example.  The code doesn't force Low to be
  // 0.0, in other words.
  //
//   xLow = -0.5*(0.90*0.0254);
//   xHigh = 0.5*(0.90*0.0254);
//   yLow = -0.5*(0.40*0.0254);
//   yHigh = 0.5*(0.40*0.0254);
  //
  // Establish boundaries comparable to xfdtd example
  //
  xLow = -0.5*7.68;
  xHigh = 0.5*7.68;
  yLow = -0.5*7.68;
  yHigh = 0.5*7.68;
  //
  // Experiment #2 code
  //
//   xLow = -0.5*15.36;
//   xHigh = 0.5*15.36;
//   yLow = -0.5*15.36;
//   yHigh = 0.5*15.36;
#ifdef CALCULATE_MESH_SIZE
  //
  // Establish the spatial deltas for each coordinate direction
  //
  dx = 0.03;
  dy = 0.03;
  //
  // Experiment #1 code
  //
//   dx = 0.015;
//   dy = 0.015;
//   dx = 0.01;
//   dy = 0.01;
  //
  // Given the two above, the mesh size is determined
  //
  nx = (int)ceil((xHigh - xLow)/dx);
  ny = (int)ceil((yHigh - yLow)/dy);
#endif
#ifdef CALCULATE_SPATIAL_DELTAS
  //
  // Establish the mesh size
  //
  nx = 10;
  ny = 5;
  //
  // Calculate the spatial deltas based on the above data
  //
  dx = (xHigh - xLow)/(double)nx;
  dy = (yHigh - yLow)/(double)ny;
#endif


  //
  // Output information to the user
  //
  fprintf(stdout, "xLow: %lg\n", xLow);
  fprintf(stdout, "xHigh: %lg\n", xHigh);
  fprintf(stdout, "yLow: %lg\n", yLow);
  fprintf(stdout, "yHigh: %lg\n", yHigh);
  fprintf(stdout, "\n");
  //
  fprintf(stdout, "dx: %lg\n", dx);
  fprintf(stdout, "dy: %lg\n", dy);
  fprintf(stdout, "\n");
  //
  fprintf(stdout, "nx: %d\n", nx);
  fprintf(stdout, "ny: %d\n", ny);
  fprintf(stdout, "\n");
  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  ///////////////////////////////////////////////////////////////////////////////









  //
  // Establish the output array dimensions based on the user's definitions.  This can't
  // be performed until now because the nx, ny, et cetera variables wree just set in
  // a user input section.
  //
#ifdef OUTPUT_X_FIELD
  iLow = 0;
  iHigh = nx-1;
  jLow = 0;
  jHigh = ny;
#endif
#ifdef OUTPUT_Y_FIELD
  iLow = 0;
  iHigh = nx;
  jLow = 0;
  jHigh = ny-1;
#endif
#ifdef OUTPUT_Z_FIELD
  iLow = 0;
  iHigh = nx+1;
  jLow = 0;
  jHigh = ny+1;
#endif

#ifdef SINGLE_POINT
  iSingle = (int)ceil(SINGLE_POINT_X_OFFSET*(double)(iHigh - iLow) + (double)iLow);
  jSingle = (int)ceil(SINGLE_POINT_Y_OFFSET*(double)(jHigh - jLow) + (double)jLow);
  fprintf(stdout, "Single output at (x,y):(%lg, %lg) (i,j):(%d, %d)\n",
	  xLow + dx*(double)iSingle, yLow + dy*(double)jSingle,
	  iSingle, jSingle);
#endif











  //
  // Define freespace material constants
  //
  freeSpace = (struct material *)malloc(sizeof(struct material));
  freeSpace->permittivity = 1.0*EPS_0;
  freeSpace->permeability = 1.0*MU_0;
  freeSpace->conductivity = 0.0;
  //
  // Allocate memory as one contiguous array.  Ez contains one more row/column
  // of data on order to allow for enforcement of boundary conditions.  The extra
  // layer is not necessary on the planar field arrays.  However, the planar field
  // arrays require one less on column/row in one of the directions in order to match
  // the mesh of the normal (Ez) dimensions.
  //
  // The valid inner region of the normal field is (1 -> nx) with 0 and nx+1 as boundary
  // zones where absorbing boundaries might be implemented.  Similarly, (1 -> ny) are the
  // valid inner values with 0 and nx+1 as boundary zones where absorbing boundaries might
  // be implemented.  Both variables are, of course, allowed to vary at the same time
  // making up this two-dimensional mesh.
  //
#ifdef E_NORMAL
  ezBase = (double *)malloc(((nx)+2)*((ny)+2)*sizeof(double));
  ez = (double **)malloc(((nx)+2)*sizeof(double *));
  for(i=0; i<=(nx)+1; i++)
    {
    ez[i] = &(ezBase[i*((ny)+2)]);
    for(j=0; j<=ny+1; j++)
      {
      ez[i][j] = 0.0;
      }
    }
  hxBase = (double *)malloc((nx)*(ny+1)*sizeof(double));
  hx = (double **)malloc((nx)*sizeof(double *));
  for(i=0; i<nx; i++)
    {
    hx[i] = &(hxBase[i*((ny)+1)]);
    for(j=0; j<ny+1; j++)
      {
      hx[i][j] = 0.0;
      }
    }
  hyBase = (double *)malloc((nx+1)*(ny)*sizeof(double));
  hy = (double **)malloc((nx+1)*sizeof(double *));
  for(i=0; i<(nx)+1; i++)
    {
    hy[i] = &(hyBase[i*(ny)]);
    for(j=0; j<ny; j++)
      {
      hy[i][j] = 0.0;
      }
    }
#endif
  //
  // While at this level the fields could have been labeled h1, h2, h3 as generic names
  // resulting in only one set of allocation loops, this may have caused confusion
  // later in the code where the actual field updates occur.  Thus, this extra code allows
  // clarity later when presenting the update equations.  This is, after all, an example code
  // for FDTD.
  //
#ifdef H_NORMAL
  hzBase = (double *)malloc(((nx)+2)*((ny)+2)*sizeof(double));
  hz = (double **)malloc(((nx)+2)*sizeof(double *));
  for(i=0; i<=(nx)+1; i++)
    {
    hz[i] = &(hzBase[i*((ny)+2)]);
    for(j=0; j<=ny+1; j++)
      {
      hz[i][j] = 0.0;
      }
    }
  exBase = (double *)malloc((nx)*(ny+1)*sizeof(double));
  ex = (double **)malloc((nx)*sizeof(double *));
  for(i=0; i<nx; i++)
    {
    ex[i] = &(exBase[i*((ny)+1)]);
    for(j=0; j<ny+1; j++)
      {
      ex[i][j] = 0.0;
      }
    }
  eyBase = (double *)malloc((nx+1)*(ny)*sizeof(double));
  ey = (double **)malloc((nx+1)*sizeof(double *));
  for(i=0; i<(nx)+1; i++)
    {
    ey[i] = &(eyBase[i*(ny)]);
    for(j=0; j<ny; j++)
      {
      ey[i][j] = 0.0;
      }
    }
#endif
  material = (struct material ***)malloc(((nx)+2)*sizeof(struct material **));
  for(i=0; i<=(nx)+1; i++)
    {
    material[i] = (struct material **)malloc(((ny)+2)*sizeof(struct material *));
    }
  //
  // Start out all point as free space
  //
  for(i=0; i<=nx+1; i++)
    {
    for(j=0; j<=ny+1; j++)
      {
      material[i][j] = freeSpace;
      }
    }





  ///////////////////////////////////////////////////////////////////////////////
  // Section two of the problem specification where the allocated memory
  // may now be filled according to the actual geometry.
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  //


  //
  // Establish start and stop times for this simulation
  //
  startSimulationTime = 0.0;
  stopSimulationTime = 3.5e-8;
  




  

  //
  // Establish plotting values
  //
  plotTime = startSimulationTime;
  //
  // A specific plot delta
  //
  plotDelta = stopSimulationTime/200.0;;


  //
  // Define a material as copper
  //
  copper = (struct material *)malloc(sizeof(struct material));
  copper->permittivity = 1.0*EPS_0;
  copper->permeability = 1.0*MU_0;
  copper->conductivity = 5.8e7;

  //
  // Create a plate of copper within the mesh.  Actually, this is an
  // infinitely long strip coming out and into the screen. Position it at -3.0
  // in the x direction with a thickness of 0.07.
  //
  for(i=(int)ceil((-2.0-xLow)/(double)dx); i<=(int)ceil((-2.0+0.17-xLow)/(double)dx); i++)
    {
    for(j=(int)ceil((-1.5-yLow)/(double)dy); j<=(int)ceil((1.5-yLow)/(double)dy); j++)
      {
      material[i][j] = copper;
      }
    }







  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  ///////////////////////////////////////////////////////////////////////////////




  //
  // Calculate the delta time from the previous data.  One must wait
  // until this late because, the local velocity of light within ALL
  // regions of the mesh has not been set.
  //
  deltaTime = FLT_MAX;
  for(i=0; i<=nx+1; i++)
    {
    for(j=0; j<=ny+1; j++)
      {
      velocity = 1.0/sqrt(material[i][j]->permittivity*material[i][j]->permeability);
      deltaTimeEstimate = (6.0/7.0) / (velocity * sqrt( 1.0/(dx*dx) + 1.0/(dy*dy) ));
      //
      // Reducing the deltaTime is generally good and increasing it is generally not.
      // FDTD and most other time-domain local-operator algorithms will support smaller
      // temporal deltas than required.
      //
      if (deltaTimeEstimate < deltaTime)
	{
	deltaTime = deltaTimeEstimate;
	}
      }
    }




  /////////////////////////////////////////////////////////////////////////
  // Sourcing arena where all static time-dependent sources would be applied
  // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv















  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  /////////////////////////////////////////////////////////////////////////

  //
  // Output information to the user
  //
  fprintf(stdout, "startSimulationTime: %lg\n", startSimulationTime);
  fprintf(stdout, "calculated deltaTime: %lg\n", deltaTime);
  fprintf(stdout, "stopSimulationTime: %lg\n", stopSimulationTime);
  fprintf(stdout, "plotDelta: %lg\n", plotDelta);
  fprintf(stdout, "\n");



#ifdef SINGLE_POINT
  //
  // Open the single point file, if requested.
  //
  sprintf(filename, "c_%f_%f.dat", SINGLE_POINT_X_OFFSET, SINGLE_POINT_Y_OFFSET);
  singleFilePointer = fopen(filename, "w");
#endif


  //
  // Perform the time advances utilizing a time based counter sequence
  //
  currentSimulationTime = startSimulationTime;
  currentIteration = 0;
  while(currentSimulationTime < stopSimulationTime)
    {
    //
    // Print to standard output the iteration number and current simulated time
    // This yields the mapping between currentIteration and currentSimulationTime
    // in case the deltaTime variable changes during the simulation.  The
    // real quantity is currentSimulationTime, however it is painful to name a file
    // with a floating-point number.  Thus, the currentSimulationTime maps into
    // a currentIteration which yields the appropriate file fora particular point
    // in time.
    //
    fprintf(stdout, "#%d %lgsec", currentIteration, currentSimulationTime);
    //
    // Plotting only at the requested times
    //
    if ( fabs(currentSimulationTime - plotTime) <= plotDelta)
      {
      //
      // E_NORMAL selected, so pick from amongst the output options
      //
#ifdef E_NORMAL
#ifdef OUTPUT_X_FIELD
      moduloOutput(hx, "hx", iLow, iHigh, jLow, jHigh,
		   xLow, xHigh, yLow, yHigh, currentIteration);
#endif
#ifdef OUTPUT_Y_FIELD
      moduloOutput(hy, "hy", iLow, iHigh, jLow, jHigh,
		   xLow, xHigh, yLow, yHigh, currentIteration);
#endif
#ifdef OUTPUT_Z_FIELD
      moduloOutput(ez, "ez", iLow, iHigh, jLow, jHigh,
		   xLow, xHigh, yLow, yHigh, currentIteration);
#endif
#endif
      //
      // H_NORMAL selected, so pick from amongst the output options
      //
#ifdef H_NORMAL
#ifdef OUTPUT_X_FIELD
      moduloOutput(ex, "ex", iLow, iHigh, jLow, jHigh,
		   xLow, xHigh, yLow, yHigh, currentIteration);
#endif
#ifdef OUTPUT_Y_FIELD
      moduloOutput(ey, "ey", iLow, iHigh, jLow, jHigh,
		   xLow, xHigh, yLow, yHigh, currentIteration);
#endif
#ifdef OUTPUT_Z_FIELD
      moduloOutput(hz, "hz", iLow, iHigh, jLow, jHigh,
		   xLow, xHigh, yLow, yHigh, currentIteration);
#endif
#endif
      //
      // Upgrade the plotTime to the next requested output time
      //
      plotTime += plotDelta;
      }





    //
    // Output single point information, if requested.
    //
#ifdef SINGLE_POINT
    //
    // E_NORMAL selected, so pick from amongst the output options
    //
#ifdef E_NORMAL
#ifdef OUTPUT_X_FIELD
    value = hx[iSingle][jSingle];  
#endif
#ifdef OUTPUT_Y_FIELD
    value = hy[iSingle][jSingle];  
#endif
#ifdef OUTPUT_Z_FIELD
    value = ez[iSingle][jSingle];
#endif
#endif
    //
    // H_NORMAL selected, so pick from amongst the output options
    //
#ifdef H_NORMAL
#ifdef OUTPUT_X_FIELD
    value = ex[iSingle][jSingle];  
#endif
#ifdef OUTPUT_Y_FIELD
    value = ey[iSingle][jSingle];  
#endif
#ifdef OUTPUT_Z_FIELD
    value = hz[iSingle][jSingle];
#endif
#endif
    fprintf(singleFilePointer, "%lg %lg ", currentSimulationTime, value);
#endif





    /////////////////////////////////////////////////////////////////////////
    // Sourcing arena where all time-dependent sources would be applied in
    // either hard or soft forms
    // vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

    //
    // Apply a source of a plane wave at the yHigh portion of the mesh with
    // a frequency of 8.3 GHz.   This can easily be a hard source since it is
    // at the outer boundary and no field updates are performed on this actual
    // location.
    //
//     for(i=0; i<=nx+1; i++)
//       {
//       ez[i][ny+1] = 1.0*sin(2.0*M_PI*8.3e9*currentSimulationTime);
//       }
    //
    // A point source sinusoid.  Hails from xfdtd/xyfdtd examples.
    //
//     frequency = 1.0/(50.0*deltaTime);
//     ez[(int)((double)nx*4.0/10.0)][(int)((double)ny*3.0/10.0)] +=
//       1.0*sin(2.0*M_PI*frequency*currentSimulationTime);
    //
    // A shaped pulse.  Hails from xfdtd/xyfdtd examples.
    //
    etmp = currentSimulationTime/(20.0*deltaTime);
    if (etmp <= 1.0)
      {
      ez[(int)((double)nx*6.0/10.0)][(int)((double)ny*6.5/10.0)] +=
  	0.5*(sin(2.0*M_PI*etmp - M_PI/2.0) + 1.0);
//       hz[(int)((double)nx*4.0/10.0)][(int)((double)ny*3.0/10.0)] +=
//  	0.5*(sin(2.0*M_PI*etmp - M_PI/2.0) + 1.0);
#ifdef SINGLE_POINT
    fprintf(singleFilePointer, " %lg ", 0.5*(sin(2.0*M_PI*etmp - M_PI/2.0) + 1.0));
#endif
      }
    else
      {
#ifdef SINGLE_POINT
      fprintf(singleFilePointer, " %lg ", 0.0);
#endif
      }
      

    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    /////////////////////////////////////////////////////////////////////////


    //
    // Truncate the `single' data file
    //
    fprintf(singleFilePointer, "\n");




#ifdef E_NORMAL
    //
    // Advance the normal field in time via the surrounding planar fields.
    // Remember that the 0 and nx+1 normal field quantities are reserved
    // for boundary condition specification and thus are NOT updated here.
    //
    for(i=2; i<nx; i++)
      {
      for(j=2; j<ny; j++)
	{
	eps = material[i][j]->permittivity;
// 	etmp = (eps + material[i][j]->conductivity*deltaTime);
// 	one = eps/etmp;
// 	two = deltaTime/etmp;
// 	ez[i][j] = ( one*ez[i][j] + two*( ( hy[i][j-1] - hy[i-1][j-1] )/dx -
// 					  ( hx[i-1][j] - hx[i-1][j-1] )/dy ) );
	//
	d2x = ( hy[i][j-1] - hy[i-1][j-1] )/dx;
	d2y = ( hx[i-1][j] - hx[i-1][j-1] )/dy;
	//
	d4x = d2x + (3.0*hy[i][j-1] - 3.0*hy[i-1][j-1] - hy[i+1][j-1] + hy[i-2][j-1])/(24.0*dx);
	d4y = d2y + (3.0*hx[i-1][j] - 3.0*hx[i-1][j-1] - hx[i-1][j+1] + hx[i-1][j-2])/(24.0*dy);
	//
	ez[i][j] = ez[i][j] + (deltaTime/eps)*d4x - (deltaTime/eps)*d4y;
	}
      }
    

    //
    // Advance the planar fields via the surrounding normal fields.  Remember
    // that these fields are spatially offset from the normal fields which
    // accounts for the array indices not matching to identical physical
    // locations. 
    //
    for(i=0; i<nx; i++)
      {
      for(j=1; j<ny; j++)
	{
	mu = material[i][j]->permeability;
//   	one = deltaTime/mu;
//   	hx[i][j] -= one*( ez[i+1][j+1] - ez[i+1][j] )/dy;
	//
	d2y = ( ez[i+1][j+1] - ez[i+1][j] )/dy;
	//
	d4y = d2y + (3.0*ez[i+1][j+1] - 3.0*ez[i+1][j] - ez[i+1][j+2] + ez[i+1][j-1])/(24.0*dy);
	//
	hx[i][j] = hx[i][j] - (deltaTime/mu)*d4y;
	}
      }
    for(i=1; i<nx; i++)
      {
      for(j=0; j<ny; j++)
	{
	mu = material[i][j]->permeability;
//  	one = deltaTime/mu;
//  	hy[i][j] += one*( ez[i+1][j+1] - ez[i][j+1] )/dx;
	//
	d2x = ( ez[i+1][j+1] - ez[i][j+1] )/dx;
	//
	d4x = d2x + (3.0*ez[i+1][j+1] - 3.0*ez[i][j+1] - ez[i+2][j+1] + ez[i-1][j+1])/(24.0*dx);
	//
	hy[i][j] = hy[i][j] + (deltaTime/mu)*d4x;
	}
      }
#endif


#ifdef H_NORMAL
    //
    // Advance the normal field in time via the surrounding planar fields.
    // Remember that the 0 and nx+1 normal field quantities are reserved
    // for boundary condition specification and thus are NOT updated here.
    //
    for(i=1; i<=nx; i++)
      {
      for(j=1; j<=ny; j++)
	{
	mu = material[i][j]->permeability;
	two = deltaTime/mu;;
	hz[i][j] += two*( ( ex[i-1][j] - ex[i-1][j-1] )/dy -
			  ( ey[i][j-1] - ey[i-1][j-1] )/dx );
	}
      }


    //
    // Advance the planar fields via the surrounding normal fields.  Remember
    // that these fields are spatially offset from the normal fields which
    // accounts for the array indices not matching to identical physical
    // locations. 
    //
    for(i=0; i<nx; i++)
      {
      for(j=0; j<=ny; j++)
	{
	eps = material[i][j]->permittivity;
	etmp = (eps + material[i][j]->conductivity*deltaTime);
	one = eps/etmp;
	two = deltaTime/etmp;
	ex[i][j] = one*ex[i][j] + two*( hz[i+1][j+1] - hz[i+1][j] )/dy;
	}
      }
    for(i=0; i<=nx; i++)
      {
      for(j=0; j<ny; j++)
	{
	eps = material[i][j]->permittivity;
	etmp = (eps + material[i][j]->conductivity*deltaTime);
	one = eps/etmp;
	two = deltaTime/etmp;
	ey[i][j] = one*ey[i][j] - two*( hz[i+1][j+1] - hz[i][j+1] )/dx;
	}
      }
#endif


    //
    // Apply boundary conditions roughly here for the next iteration.
    // Without anything in this section, the outer boundary may act as
    // a PEC in the case where the normal field component is magnetic.
    // This is because the outer field values which were reserved for
    // boundary specifications are the normal components and are set to
    // 0.0 by default.  This would conform to the PEC for electric normal
    // fields and to the PMC for magnetic normal fields.
    //



    fprintf(stdout, "\n");
    //
    // Advance time.
    //
    currentSimulationTime += deltaTime;
    currentIteration++;
    }
  


  //
  // Output the last iteration regardless of the modulo
  //



  //
  // Print to standard output the iteration number and current simulated time
  // This yields the mapping between currentIteration and currentSimulationTime
  // in case the deltaTime variable changes during the simulation.  The
  // real quantity is currentSimulationTime, however it is painful to name a file
  // with a floating-point number.  Thus, the currentSimulationTime maps into
  // a currentIteration which yields the appropriate file fora particular point
  // in time.
  //
  fprintf(stdout, "#%d %lgsec", currentIteration, currentSimulationTime);
  //
  // Output single point information, if requested.
  //
#ifdef SINGLE_POINT
  //
  // E_NORMAL selected, so pick from amongst the output options
  //
#ifdef E_NORMAL
#ifdef OUTPUT_X_FIELD
  value = hx[iSingle][jSingle];  
#endif
#ifdef OUTPUT_Y_FIELD
  value = hy[iSingle][jSingle];  
#endif
#ifdef OUTPUT_Z_FIELD
  value = ez[iSingle][jSingle];
#endif
#endif
  //
  // H_NORMAL selected, so pick from amongst the output options
  //
#ifdef H_NORMAL
#ifdef OUTPUT_X_FIELD
  value = ex[iSingle][jSingle];  
#endif
#ifdef OUTPUT_Y_FIELD
  value = ey[iSingle][jSingle];  
#endif
#ifdef OUTPUT_Z_FIELD
  value = hz[iSingle][jSingle];
#endif
#endif
  fprintf(singleFilePointer, "%lg %lg\n", currentSimulationTime, value);
#endif






  //
  // E_NORMAL selected, so pick from amongst the output options
  //
#ifdef E_NORMAL
#ifdef OUTPUT_X_FIELD
  moduloOutput(hx, "hx", iLow, iHigh, jLow, jHigh,
	       xLow, xHigh, yLow, yHigh, currentIteration);
#endif
#ifdef OUTPUT_Y_FIELD
  moduloOutput(hy, "hy", iLow, iHigh, jLow, jHigh,
	       xLow, xHigh, yLow, yHigh, currentIteration);
#endif
#ifdef OUTPUT_Z_FIELD
  moduloOutput(ez, "ez", iLow, iHigh, jLow, jHigh,
	       xLow, xHigh, yLow, yHigh, currentIteration);
#endif
#endif
  //
  // H_NORMAL selected, so pick from amongst the output options
  //
#ifdef H_NORMAL
#ifdef OUTPUT_X_FIELD
  moduloOutput(ex, "ex", iLow, iHigh, jLow, jHigh,
	       xLow, xHigh, yLow, yHigh, currentIteration);
#endif
#ifdef OUTPUT_Y_FIELD
  moduloOutput(ey, "ey", iLow, iHigh, jLow, jHigh,
	       xLow, xHigh, yLow, yHigh, currentIteration);
#endif
#ifdef OUTPUT_Z_FIELD
  moduloOutput(hz, "hz", iLow, iHigh, jLow, jHigh,
	       xLow, xHigh, yLow, yHigh, currentIteration);
#endif
#endif





  fprintf(stdout, "\n");
  //
  // Close any remaining open files
  //
  fclose(singleFilePointer);

  }
