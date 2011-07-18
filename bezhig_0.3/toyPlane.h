


//
// Calculate the mesh size based on the user specified spatial deltas
//
#define CALCULATE_MESH_SIZE
//
// Calculate the spatial deltas based on the user specified mesh size
//
//#define CALCULATE_SPATIAL_DELTAS
//
// Use the PNM image format to output data
//
// PNM_COLORBAR produces a colorbar along the left side of the 2D PNM image
//
#define PNM
#define PNM_COLORBAR
//
// Use the MTV output format which connects to the PlotMTV plotting package
// This has nothing to do with MTV on the TV or the MTV raytracing software
// formatted files.  Both PNM and MTV output formats may be used simultaneously.
//
// MTV_ASCII uses ASCII numbers in the file
// MTV_BINARY uses binary data values in the file
//
// Only use one if these options.  Binary uses less space and yields more precision yet
// is MUCH less readable than ASCII for humans.
//
#define MTV
#define MTV_ASCII
//#define MTV_BINARY
//
// Choose either E_NORMAL or H_NORMAL to pick between (Ez, Hx, Hy) and (Hz, Ex, Ey)
// two-dimensional propagation, respectively.
//
// Only one may be chosen at a time.
//
#define E_NORMAL
//#define H_NORMAL


//
// Define output options, note that pending the E_NORMAL or H_NORMAL choice
// that the OUTPUT_X_FIELD will output E or H values.
//
// Only one may be chosen at a time.
//
//#define OUTPUT_X_FIELD
//#define OUTPUT_Y_FIELD
#define OUTPUT_Z_FIELD


//
// Allow a single data point to be output throughout the entire simulation.
// Because things are getting messy already, this will be specified as a percentage
// offset from the lower mesh boundaries.  So, 0.40 would calculated 40% from the
// Low boundary towards the High boundary.  The field output is defined by the
// aforementioned #defines.  Again, I could make the code more flexible, yet it already
// hides much of the algorithm.  This output is 55% from the xLow to the xHigh and 45% 
// from the yLow to the yHigh.
//
#define SINGLE_POINT
#define SINGLE_POINT_X_OFFSET 0.55
#define SINGLE_POINT_Y_OFFSET 0.45








//
// Speed of light in free-space, c, in meter/second
//
#define LIGHT_SPEED	299792458.0
#define INVERSE_LIGHT_SPEED 0.0000000033356409519815204957557671447491851179258151984597290969874899254470237540131846812503868926
#define LIGHT_SPEED_SQUARED 89875517873681764.0
//
// Permeability of free-space, mu  , in Newton/Ampere^2
//                               0                    
// Approximates to 12.566370614e-7 N/(A^2)
//
//#define MU_0		(M_PI*4.0e-7)
//#define MU_0 0.00000125663706143591729538505735331180115367886775
#define MU_0 1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6
//
// Permittivity of free-space, e  , in Farad/meter
//                              0    
// Approximates to 8.854187817e-12 F/m
//
//#define EPS_0		(1.0/(MU_0*LIGHT_SPEED_SQUARED))
//#define EPS_0 0.00000000000885418781762038985053656303171075026060
#define EPS_0 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12








//
// A simple materials framework allowing definition of a material once
// and reference to it many times.
//
typedef struct material Material;

struct material
{
  float conductivity;
  float permittivity;
  float permeability;
};









//
// Two very handy functions from the `Schneider' plotting code
// samples and links found on the CEMTACH web page.  cmap2 defines
// a specific color mapping to floating point values and hsvrgb
// converts that output into rgb which is used to create PNM formatted
// 2d image files.
//
//
void
cmap2(double *h,
      double *s,
      double *v,
      double xpnt);
void
hsvrgb(double h,
       double s,
       double v,
       int *ir,
       int *ig,
       int *ib);

void
moduloOutput(double **fieldData,
	     char *label,
	     int iLow, int iHigh,
	     int jLow, int jHigh,
	     double xLow, double xHigh,
	     double yLow, double yHigh,
	     int currentIteration);
