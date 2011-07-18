

//
// The common output routine from ToyPlaneFDTD. Since the code is
// called twice from the code, it is actually easier to maintain in a
// separate routine and also lends clarity to the numerical algorithm
// in ToyPlaneFDTD.c.
//
#include <stdio.h>
#include <float.h>
#include <math.h>

#include "toyPlane.h"

void
moduloOutput(double **fieldData,
	     char *label,
	     int iLow, int iHigh,
	     int jLow, int jHigh,
	     double xLow, double xHigh,
	     double yLow, double yHigh,
	     int currentIteration)
{
  int i,j;
  double min, max, norm, drange, value;
  FILE *imageFilePointer;
  char filename[1024];
  double etmp;
#ifdef PNM
  //
  // Creation of PNM image files variables
  //
  int ir,ib,ig;
  double h,s,v;
  char red, green, blue;
#endif

#ifdef PNM
  //
  // Reset the min and max values for a new search of the array
  //
  min = FLT_MAX;
  max = -FLT_MAX;
  //
  // Autoscale over the array ranges just picked using the requested array
  //
  for(i=iLow; i<=iHigh; i++)
    {
    for(j=jLow; j<=jHigh; j++)
      {
      value = fieldData[i][j];
      //
      // Actually autoscale according to the chosen array value
      //
      if (value < min) {min = value;}
      if (value > max) {max = value;}
      }
    }
  //
  // The norm is the largest value of either max or min quantities.
  //
  if (fabs(max) > fabs(min))
    {
    norm = fabs(max);
    }
  else
    {
    norm = fabs(min);
    }
  //
  // Enforce a smallest value to avoid devide by zero failures.
  //
  if (norm < FLT_EPSILON)
    {
    norm = 1.0;
    }
  //
  // Dynamic range calculation
  //
  drange = max - min;
  //
  // Echo min and max to output since the image cannot display this information
  //
  fprintf(stdout, "\t%lg < %s < %lg", min, label, max);


#endif






#ifdef MTV
  //
  // Create an mtv datafile
  //
  sprintf(filename, "c_%05d.mtv", currentIteration);
  imageFilePointer = fopen(filename, "w");
  //
  // Label the MTV plot
  //
  fprintf(imageFilePointer,"$ DATA=CONTOUR NAME=\"%s\"\n", label );
  //
  // Output the MTV header information describing the data
  //
  fprintf(imageFilePointer,"%% contfill\n");
  fprintf(imageFilePointer,"%% nsteps=25\n");
  fprintf(imageFilePointer,"%% xmin = %lg  xmax = %lg\n", xLow, xHigh);
  fprintf(imageFilePointer,"%% ymin = %lg  ymax = %lg\n", yLow, yHigh);
  fprintf(imageFilePointer,"%% nx   = %d\n",iHigh - iLow + 1);
  fprintf(imageFilePointer,"%% ny   = %d\n",jHigh - jLow + 1);
  //
  // MTV allows either a ASCII or binary dump of the array data.
  // ASCII is more portable across systems, yet looses precision.
  //
#ifdef MTV_ASCII
  for(j=jLow; j<=jHigh; j++)
    {
    for(i=iLow; i<=iHigh; i++)
      {
      value = fieldData[i][j];
      //
      // Okay, so write the value to the file!
      //
      fprintf(imageFilePointer,"%lf\n",value);
      }
    }
#endif
  //
  // MTV allows either a ASCII or binary dump of the array data.
  // Binary maintains all precision, yet is less portable across systems.
  //
#ifdef MTV_BINARY
  fprintf(imageFilePointer,"%% binary=True\n");
  //
  // Simply write the entire array directly
  //
  fwrite((char*)fieldData, sizeof(float), (iHigh - iLow + 1)*(jHigh - jLow + 1), imageFilePointer);
#endif
  //
  // Truncate the MTV file
  //
  fprintf(imageFilePointer,"$ end\n");
  fclose(imageFilePointer);
#endif









#ifdef PNM
  //
  // Create a PNM datafile.  Note that if the dx and dy values are not identical
  // this image plot may appear squished on one of the coordinate directions.
  //
  sprintf(filename, "c_%05d.pnm", currentIteration);
  imageFilePointer = fopen(filename, "wb");
  //
  // Use P3 for ASCII and P6 for binary data in color
  // Add 10 to the mesh size for color bar, if present.
  //
  fprintf(imageFilePointer,"P6\n");
#ifdef PNM_COLORBAR
  fprintf(imageFilePointer,"%d %d\n", (iHigh-iLow+1)+10, (jHigh-jLow+1));
#else
  fprintf(imageFilePointer,"%d %d\n", (iHigh-iLow+1), (jHigh-jLow+1));
#endif
  //
  // Color 3 byte values each limited to 0->255 in range
  //
  fprintf(imageFilePointer,"255\n");
  //
  // Okay, write out the data in character form
  //
  for(j=jHigh; j>=jLow; j--)
    {
#ifdef PNM_COLORBAR
    //
    // Create the color bar segment on this line of the image
    //
    etmp = 2.0*(double)((j-(jHigh - jLow + 1)/2))/(jHigh - jLow + 1);
    if (etmp > 1.0) {etmp = 1.0;}
    if (etmp < -1.0) {etmp = -1.0;}
    for(i=-10; i<iLow; i++)
      {
      //
      // Locate the appropriate color within the ranged color map
      // and write the RGB values into the image file.
      //
      cmap2(&h,&s,&v,etmp);
      hsvrgb(h,s,v,&ir,&ig,&ib);
      red = ir; green = ig; blue = ib;
      fwrite(&red, sizeof(char), 1, imageFilePointer);
      fwrite(&green, sizeof(char), 1, imageFilePointer);
      fwrite(&blue, sizeof(char), 1, imageFilePointer);
      }
#endif
    for(i=iLow; i<=iHigh; i++)
      {
      value = fieldData[i][j];
      //
      // Locate the appropriate color within the ranged color map
      // and write the RGB values into the image file.
      //
      cmap2(&h,&s,&v,value/norm);
      //
      // Convert from HSV color map into RGB values for PNM format
      //
      hsvrgb(h,s,v,&ir,&ig,&ib);
      red = ir; green = ig; blue = ib;
      fwrite(&red, sizeof(char), 1, imageFilePointer);
      fwrite(&green, sizeof(char), 1, imageFilePointer);
      fwrite(&blue, sizeof(char), 1, imageFilePointer);
      }
    }
  fclose(imageFilePointer);
#endif




}
