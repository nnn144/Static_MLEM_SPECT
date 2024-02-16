#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

// !!!!!!!!!!!!! most important pre-compile options !!!!!!!!!!!

//#define ATTENUATION_YES //-- uncomment to add attenuation, then, though, we need attenuation files
#define BLURRING_YES // uncomment to allow blurring due to geom correction and stuff
                     // If defined AXIS_DETECTOR_DISTANCE is needed and
                     // kernel function and its parameters, (convolution) function,
                     // and function computing distance detector-voxel need to be
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define POSITION_OF_CENTER_OF_ROTATION 63.5 //63.5 //43.5 //31.5 // pixels
#define N_VOLUME_XY 128  //  VOlume size in pixel units
#define N_VOLUME_Z  1 // N_VOLUME_XY x N_VOLUME_XY x N_VOLUME_Z
// -- voxel size in cm !!!
#define VOXEL_SIZE 1 // in units of pixel size!!!  At so far only 1 is supported
#define VOXEL_SUBDIVISION   3 // Means that each voxel is V_S * V_S * V_S subpoints
#define AXIS_DETECTOR_DISTANCE 170.0 // in cm, we will convert it to voxel sizes later
#define NO_DIVERGENT_RAYS 7 // can be 1, 7, 13, 19, 31, 37, 43

#define N_detector_r 128  // detectror
#define N_detector_z 1  // pixelsx
#define N_views      64  // no-views
// detector pixel size in cm -- needed for attenuation purposes and to determine axis-detector-distance
#define d_pixel 0.47952

#define ANGULAR_SPAN -M_PI // negative or positive for rotation direction, initial negative
#define STARTING_ANGLE 0 //-0.7853981633974483 // == -M_PI / 4, mostly, we use zero or pi

#define EPSILON 0.000001
#define ONE_OVER_EPSILON 1000000.
static int iswap_arg;
#define ISWAP(a, b) iswap_arg = a;  a = b;  b = iswap_arg
#define ISWAP3(a, b, c) iswap_arg = a;  a = b;  b = c;  c = iswap_arg


#ifdef BLURRING_YES  // --------- BLURRING-related stuff -------
// required definitions
#define point_response_a 0.0256685  // inverse cm
#define point_response_b  0.208467  // cm
// with these two defined, point response on the detector is
//    ~exp(-r^2/s^2), where r is pixel_
#define MAX_ARGUMENT_TO_EXP 3. // this is good enough
// required variables
int length_kernel_lookup_table = 2097152; // this length is smaller than
double one_over_kernel_lookup_scaling_factor;
double *kernel_lookup_table;
// required functions
void initialize_kernel_lookup_table();
double kernel_sigma_squared(double sinus, double cosinus, int nx, int ny);
void convolve_with_blurring_kernel(float **a, float **a_tmp, float **h, int N, int K, float Nh);
#endif    // -------------------- end blurring stuff  ----------


int main(int argc, char **argv);
#ifdef ATTENUATION_YES
void project_volume(int Nx, int Ny, int Nz, int Nr, int NZ, int Nang, char *sysmat_fname, float ***attn_map);
float integrate_one_ray(float ***volume, int N[], double R0[], double directions[]);
#else
void project_volume(int Nx, int Ny, int Nz, int Nr, int NZ, int Nang, char *sysmat_fname);
#endif
float **DivergentRays(int N, float max_distance);
#ifdef NO_RESPIRATORY_MOTION_FRAMES
void deformation_calculator(float *x, float *y, float *z, int nang);
#endif

void Read_volume(float ***Data, int x, int y, int z, char *fname);
// // // // // // // // // // // // // // // // // // // // // // // // // // //
 // // // // // // // // // // // // // // // // // // // // // // // // // // //
// // // // // // // // // // // // // // // // // // // // // // // // // // //
int main(int argc, char **argv)
{
  int Nx, Ny, Nz;      // these things should be changeable
  int Nr, NZ, Nang;  // same here
  FILE *fid;
  float *params;
#ifdef ATTENUATION_YES  // attn tensor declaration
  float ***attn_tensor;

  float *tempMap;
  int n, m, k;
#endif


#ifdef ATTENUATION_YES
// --- initial assignments ---
  Nx = Ny = N_VOLUME_XY;  Nz = N_VOLUME_Z;
  Nr = N_detector_r; NZ = N_detector_z;  Nang = N_views;

// ----------------------------
  if(argc != 3)
  { fprintf(stderr, "Usage: %c a.out sysmat.fname attn_matrix.fname\n", 37);
    return 1;  }
  if( (fid=fopen(argv[2], "rb"))==NULL) // read attn map
  {
    fprintf(stderr, "Failed to open attenuation matrix file %s\n", argv[3]);
    exit(1);
  }
  else
  {
    attn_tensor =  f3tensor(0, Nz-1, 0, Ny-1, 0, Nx-1);
    fread(&attn_tensor[0][0][0], sizeof(short), Nx*Ny*Nz, fid);
    fclose(fid);

    project_volume(Nx, Ny, Nz, Nr, NZ, Nang, argv[1], attn_tensor);
    free_f3tensor(attn_tensor, 0, Nx-1, 0, Ny-1, 0, Nz-1);
  }
#else
// --- initial assignments ---
  Nx = Ny = N_VOLUME_XY;  Nz = N_VOLUME_Z;
  Nr = N_detector_r;  NZ = N_detector_z;  Nang = N_views;

// --- code
  if(argc != 2)
  { fprintf(stderr, " Usage: \n %c a.out sysmat.fname \n", 37);
    return 1;}
  project_volume(Nx, Ny, Nz, Nr, NZ, Nang, argv[1]);
#endif
  return 1;
}

// // // // // // // // // // // // // // // // // // // // // // // // // // //
 // // // // // // //  function that creates the system matrix
// // // // // // // // // // // // // // // // // // // // // // // // // // //
#ifdef ATTENUATION_YES
void project_volume(int Nx, int Ny, int Nz, int Nr, int NZ, int Nang, char *sysmat_fname, float ***attn_map)
#else
void project_volume(int Nx, int Ny, int Nz, int Nr, int NZ, int Nang, char *sysmat_fname)
#endif
{
  int voxel_subsamples, volume_index, sinogram_index, index,  no_entries, *indices;
  int Ndet, Nsub, nx, ny, nz, nZ, nr, nang, n,m,k, n_div, NDivergGeom =  NO_DIVERGENT_RAYS;
  size_t size_float, size_int;
  float *entries;
  float x, y, z;
  float **Rsub, weight, *cos_theta, *sin_theta, r_shift, voxel_to_det, center_to_det, **RDivergGeom;
  FILE *fid;
#ifdef BLURRING_YES
  float **detector_tmp, **detector_matrix, **kernel;
  int N_kernel;
#endif

#ifdef ATTENUATION_YES
  double directions[3], R0[3], attn_factor;
  int N3[3];
  N3[0] = Nx;  N3[1] = Ny;  N3[2] = Nz;
  // scale attenuation tensor by 0.001 * d_pixel
  // for(k=0;k<Nz;k++)   for(m=0;m<Ny;m++)   for(n=0;n<Nx;n++) attn_map[k][m][n] *= d_pixel;
  float temp = 0.001 * d_pixel * 0.15; // Comment this will be the original code (HAORAN)
  for(k=0;k<Nz;k++)   for(m=0;m<Ny;m++)   for(n=0;n<Nx;n++) attn_map[k][m][n] *= temp;

#endif
  // ---- end declacations -----
  fprintf(stderr, "Nz = %d\n", Nz);

// --- check for supported options ---
//  if(VOXEL_SIZE != 1) nrerror("Voxel size has to be the same as detector pixel size");

  // --- helper variables ---
  size_float = sizeof(float);
  size_int   = sizeof(int);
  Ndet = Nr * NZ;
  r_shift = POSITION_OF_CENTER_OF_ROTATION;
  center_to_det = AXIS_DETECTOR_DISTANCE / d_pixel;

  // --- precompute cos and sin of theta ---
  cos_theta = (float *) malloc(size_float * Nang);
  sin_theta = (float *) malloc(size_float * Nang);
  weight = ANGULAR_SPAN / ( (double) Nang ); // here, we recycle variable "weight", it will be changed later
  for(n=0; n<Nang; n++){
    cos_theta[n] = cos(STARTING_ANGLE + weight * n);
    sin_theta[n] = sin(STARTING_ANGLE + weight * n);  }

  // --- declare indexes and entries vectors
  indices = (int   *) malloc(size_int   * Nr * NZ * Nang);
  entries = (float *) malloc(size_float * Nr * NZ * Nang);

  // --- assign voxel subdivision parameters ---
  Nsub = VOXEL_SUBDIVISION * VOXEL_SUBDIVISION * VOXEL_SUBDIVISION;
  Rsub = matrix(0, Nsub-1, 0, 2);
  weight = 1. / ( (double) VOXEL_SUBDIVISION);
  for(k=0; k<VOXEL_SUBDIVISION; k++)
    for(m=0; m<VOXEL_SUBDIVISION; m++)
      for(n=0; n<VOXEL_SUBDIVISION; n++)
      {
        volume_index = n + VOXEL_SUBDIVISION * (m + VOXEL_SUBDIVISION * k); // recucle variable volume_index
        Rsub[volume_index][0] = (n - (VOXEL_SUBDIVISION-1)*0.5) * weight - (Nx-1.)/2.;
        Rsub[volume_index][1] = (m - (VOXEL_SUBDIVISION-1)*0.5) * weight - (Ny-1.)/2.;
        Rsub[volume_index][2] = (k - (VOXEL_SUBDIVISION-1)*0.5) * weight;
      }
  weight = weight * weight * weight; // this is the right value of weight

  // compute divergent ray matrix : we do not need it if we have blurring kernel
  RDivergGeom = DivergentRays(NDivergGeom, center_to_det + Nx / sqrt(2.));

  // --- open file for writing sys mat ---
  if( (fid=fopen(sysmat_fname, "wb"))==NULL) // read attn map
  {
    fprintf(stderr, "Failed to open system matrix file %s for writing, exiting...\n", sysmat_fname);
    exit(1);
  }
  // write Nsinogram
  n = Ndet * Nang;   fwrite(&n, size_int, 1, fid); fprintf(stderr, "Nsinogram = %d\n", n);
  // write Nvolume
  n = Nx * Ny * Nz;   fwrite(&n, size_int, 1, fid); fprintf(stderr, "Nvolume = %d\n", n);

  fprintf(stderr, "Loop, slices\n");

  // main loop
  for(nz=0; nz<Nz; nz++)
  { fprintf(stderr, "%d|", nz);
    for(ny=0; ny<Ny; ny++)
      for(nx=0; nx<Nx; nx++)
      {
        volume_index = nx + Nx*(ny + Ny*nz);
        no_entries = n = 0;// n is the lowest entry used for given angle
        for(nang=0; nang<Nang; nang++)
        {
          x = nx + Rsub[0][0];
          y = ny + Rsub[0][1];
          z = nz + Rsub[0][2];
          voxel_to_det = -sin_theta[nang]*x + cos_theta[nang]*y + center_to_det;
          nr = (int) ( cos_theta[nang]*x + sin_theta[nang]*y + r_shift + voxel_to_det*RDivergGeom[0][0]);
          nZ = (int) (z + voxel_to_det*RDivergGeom[0][1]);
          index = nr + Nr*nZ;
          if(nr>=0 && nr < Nr && nZ>=0 && nZ<NZ)
          {
            indices[no_entries] = index + nang * Ndet;
            entries[no_entries] = weight * RDivergGeom[0][2];
            no_entries = n+1;
          }
          for(k=0; k<Nsub; k++)
          {
           x = nx + Rsub[0][0];
           y = ny + Rsub[0][1];
           z = nz + Rsub[0][2];
            for(n_div =0; n_div<NDivergGeom; n_div++)
            if(! (n_div==0 && k==0))
            {
              nr = (int) ( cos_theta[nang]*x + sin_theta[nang]*y + r_shift + voxel_to_det*RDivergGeom[n_div][0]);
              nZ = (int) (z + voxel_to_det*RDivergGeom[n_div][1]);
              if(nr<0 || nr >= Nr || nZ<0 || nZ >=NZ) continue; // if outside of sinogram
              index = nr + Nr*nZ;
              sinogram_index = index + nang * Ndet;
              for(m=n; m<no_entries; m++)
                if(indices[m] == sinogram_index)
                {
                  entries[m] += weight * RDivergGeom[n_div][2];
                  break;
                }
              if(m>=no_entries)
              {
                indices[no_entries] = sinogram_index;
                entries[no_entries] = weight * RDivergGeom[n_div][2];
                no_entries++;
              }
            } // end s_subdivision
          }
#ifdef ATTENUATION_YES
          R0[0] = nx+0.5;                     R0[1] = ny+0.5;                    R0[2] = nz+0.5;
          directions[0] = -sin_theta[nang];   directions[1] =  cos_theta[nang];  directions[2] = 0.;
          attn_factor = exp(-integrate_one_ray(attn_map, N3, R0, directions)); // remember, attn_map is scaled so this formula is right!
          for(m=n; m<no_entries; m++)  entries[m] *= attn_factor;
#endif
          // before finishing loop for given value of nang, reassign n
          n = no_entries;
        } // end angle loop
        // at this point, all sys matindeces/entries for given voxel are known, now all we do is save them
        fwrite(&volume_index, size_int, 1, fid);
        fwrite(&no_entries, size_int, 1, fid);
        if(no_entries == 0) continue;
        fwrite(indices,  size_int,   no_entries, fid);
        fwrite(entries, size_float, no_entries, fid);
      }
  }
  // --- done main cycle --
  fclose(fid);
  fprintf(stderr, "\nDone, cleaning up memory...\n");
  // --- free memory
  free(indices);
  free(entries);
  free(cos_theta);
  free(sin_theta);
  free_matrix(Rsub, 0, Nsub-1, 0, 2);
  return;
}


// ------ volume voxel positions ------
// ---- this function is experiment-dependant, it calculates motion

// ------ divergent rays subdivision ------
// Divergent rays, distribution follows scheme in point_response.pdf
float **DivergentRays(int N, float max_distance) // max_distance in voxels
{
  int n;
  float **R;
  double nearest_neighbor_distance, gaussian_width, weight, sqrt3_2;

  if(! (N==1 || N==7 || N==19 || N==31 || N==37 || N==43))     nrerror("Bad no of divergent rays");
  nearest_neighbor_distance = 1./max_distance; // chosen so points are never more than one pixel apart
  // remmeber that Gaussian widths are 0.0256685* D + 0. 2 CM, where D -- cm distance to detector
  gaussian_width = max_distance * VOXEL_SIZE *  d_pixel * 0.0256685 + 0.2; // sigma in terms of dimensionless pixel cooords
  //
  R = matrix(0, N-1, 0, 2);  // first two components -- coordinates, third -- weight
  sqrt3_2 = 0.5 * sqrt(3.);
  switch(N)
  {
    case 43:
      R[37][0] = -3.; R[37][1] =  2.*sqrt3_2;
      R[38][0] = -3.; R[38][1] = -2.*sqrt3_2;
      R[39][0] =  0.; R[39][1] =  4.*sqrt3_2;
      R[40][0] =  0.; R[40][1] = -4.*sqrt3_2;
      R[41][0] =  3.; R[41][1] =  2.*sqrt3_2;
      R[42][0] =  3.; R[42][1] = -2.*sqrt3_2;
      weight = exp(- 12./SQR(gaussian_width)); // Distance is sqrt(12)
      for(n=37; n<43; n++) R[n][2] = weight;
    case 37:
      R[31][0] = -3.0; R[31][1] =  0.;
      R[32][0] = -1.5; R[32][1] =  3.*sqrt3_2;
      R[33][0] = -1.5; R[33][1] = -3.*sqrt3_2;
      R[34][0] =  1.5; R[34][1] =  3.*sqrt3_2;
      R[35][0] =  1.5; R[35][1] = -3.*sqrt3_2;
      R[36][0] =  3.0; R[36][1] =  0.;
      weight = exp(- 9./SQR(gaussian_width)); // Distance is sqrt(7)
      for(n=31; n<37; n++) R[n][2] = weight;
    case 31:
      R[19][0] = -2.5; R[19][1] =     sqrt3_2;
      R[20][0] = -2.0; R[20][1] =  2.*sqrt3_2;
      R[21][0] = -2.5; R[21][1] = -   sqrt3_2;
      R[22][0] = -0.5; R[22][1] =  3.*sqrt3_2;
      R[23][0] = -2.0; R[23][1] = -2.*sqrt3_2;
      R[24][0] =  0.5; R[24][1] =  3.*sqrt3_2;
      R[25][0] = -0.5; R[25][1] = -3.*sqrt3_2;
      R[26][0] =  2.0; R[26][1] =  2.*sqrt3_2;
      R[27][0] =  0.5; R[27][1] = -3.*sqrt3_2;
      R[28][0] =  2.5; R[28][1] =     sqrt3_2;
      R[29][0] =  2.0; R[29][1] = -2.*sqrt3_2;
      R[30][0] =  2.5; R[30][1] = -   sqrt3_2;
      weight = exp(- 7./SQR(gaussian_width)); // Distance is sqrt(7)
      for(n=19; n<31; n++) R[n][2] = weight;
    case 19:
      R[13][0] = R[17][0] =  1.;
      R[14][0] = R[16][0] = -1.;
      R[18][0] =  2.;
      R[15][0] = -2.;
      R[13][1] = R[14][1] =  sqrt3_2*2.;
      R[15][1] = R[18][1] =  0.0;
      R[16][1] = R[17][1] = -sqrt3_2*2.;
      weight = exp(- 4./SQR(gaussian_width));
      for(n=13;n<19; n++) R[n][2] = weight; // Distance is 2
    case 13:
      R[7][0] =  R[12][0] =  1.5;
      R[8][0] =  R[11][0] =  0.0;
      R[9][0] =  R[10][0] = -1.5;
      R[ 7][1] =  R[ 9][1] =  sqrt3_2;
      R[10][1] =  R[12][1] = -sqrt3_2;
      R[ 8][1] =  sqrt3_2*2.;
      R[11][1] = -sqrt3_2*2.;
      weight = exp(- 3./SQR(gaussian_width));
      for(n=7;n<13; n++) R[n][2] = weight; // Distance is sqrt(3)
    case 7:
      R[1][0] =  R[5][0] =  0.5;
      R[2][0] =  R[4][0] = -0.5;
      R[6][0] =  1.0;
      R[3][0] = -1.0;
      R[3][1] =  R[6][1] =  0.0;
      R[1][1] =  R[2][1] =  sqrt3_2;
      R[4][1] =  R[5][1] = -sqrt3_2;
      weight = exp(- 1./SQR(gaussian_width));// Distance is 1 (pixel)
      for(n=1;n<7; n++) R[n][2] = weight;
    case 1:
      R[0][0] = R[0][1] = 0.;
      R[0][2] = 1.;
      break;
    default:
      nrerror("Cannot reach this point in DivergentRays");
  }
  // -- scale by nearest_neighbor_distance and scale weights
  if(N>0)
  {
    weight=1.; // automatically take center into account
    for(n=0; n<N; n++)
    {
      R[n][0] *= nearest_neighbor_distance;
      R[n][1] *= nearest_neighbor_distance;
      weight += R[n][2];
    }
    weight = 1./weight;
    // now, scale weights so they add up to one
    for(n=0; n<N; n++) R[n][2] *= weight;
  }
  return(R);
}



#ifdef ATTENUATION_YES
// ------------------------ integrate through a volume -- forattenuation -------------------------------------------
// -----------------------------------------------------------------------------------------------------------------
float integrate_one_ray(float ***volume, int N[], double R0[], double directions[])
{   // integrate ray starting at R0 and moving towards directions throught volume
// N[] == [Nx Ny Nz] == volume dimensions
// this backprojector is designed for the case if A is numbered from 0 to Nx*Ny*Nz-1
// propagation happens according to the rule:
//   x = X0 + s cosX             s = (x - X0) / cosX == (x-X0) * dx
//   y = Y0 + s cosY             s = (y - Y0) / cosY == (y-Y0) * dy
//   z = Z0 + s cosZ             s = (z - Z0) / cosZ == (z-Z0) * dz
 int coord[3] = {0, 1, 2};
 int n[3], nz, ny, nx;
 int k, tmp_index;
 double s_in, s_out, pixel_in[3], pixel_out[3];
 double path_through_pixel;
 int dn[3];
 double ds[3];
 float  ray_value;

 for(k=0; k<3; k++)
   if(directions[k]>0.)
   {
     dn[k] = 1;
     ds[k] = 1./directions[k];
   }
   else if(directions[k]<0.)
   {
     dn[k] = -1;
     ds[k] = -1./directions[k];
   }
   else
   {
     dn[k] = 0;
     ds[k] = 1.0e+10;
   }

// We move through the volume using linear variable s
// for each cube there is s_in and s_out

// Initial assignment of s_in and s_out.
 s_in = 0.;   // value at R = R0
 s_out = 8388606.;
// fprintf(stderr, "(%f,%f,%f)/(%f,%f,%f) by (%f, %f, %f)\n", R0[0], R0[1], R0[2], directions[0], directions[1], directions[2], ds[0]*dn[0], ds[1]*dn[1], ds[2]*dn[2] );

// ---------------------- all components -----------------
 for(k=0; k<3; k++)
 {
   if(dn[k] > 0) // positive propagation in the direction k
   {
     if(R0[k] >= N[k]) return(0.); // we know that we miss the volume
     n[k] = 0;   // This means that the entrance pixel <for dimension k> is zero
     pixel_in[k]  = -R0[k] * ds[k]; // The minus sign is OK, we expect R0 to be negative unless ray enters from the side
     pixel_out[k] = (N[k] - R0[k]) * ds[k];
     if(pixel_in[k] > s_in) // assign this as entrance coordinate -- positive propagation case
     {
       s_in = pixel_in[k];
       ISWAP(coord[0], coord[k]);
       n[coord[k]] = 0;
     }
   }
   else                 // negative propagation
   {
     if(R0[k] <= 0.) return 0.; // we know that we miss the volume
     n[k]         = N[k]-1;
     pixel_in[k]  =  -(N[k] - R0[k]) * ds[k];
     pixel_out[k] = R0[k] * ds[k]; /// === "-" sign is there because we also divide by dn, which is negative here
     if(pixel_in[k] > s_in) // assign this as entrance coordinate -- negative propagation case
     {
       s_in = pixel_in[k];
       ISWAP(coord[0], coord[k]);
       n[coord[k]] = N[coord[k]]-1;
     }
   }
// compute s_in and s_out for the whole volume
   if(pixel_out[k] < s_out)  s_out = pixel_out[k];
 }
  if(s_in >= s_out)
    return(0.); // ------------------ return zero if in is larger than out... however, means, we missed the volume


//  -------- now, we need to assign the right values to the indeces and to the pixel_outs AND sort them
//  --------- this process is different if we are starting outside the ROI (s_in > 0) or if we are starting inside the ROI (s_in ==0)
 if(s_in > 0.0001) // we are computing entry parameters for all rays but one, which enters the ROI last
 {
   pixel_out[coord[0]] = pixel_in[coord[0]] + ds[coord[0]];
   for(k=1; k<3;k++)
   {
     pixel_in[coord[k]] = s_in;
     n[coord[k]] = (int) floor(R0[coord[k]] + s_in * directions[coord[k]]);
     if(dn[coord[k]]==1) pixel_out[coord[k]] =  (n[coord[k]] + 1 - R0[coord[k]]) * ds[coord[k]];
     else                pixel_out[coord[k]] = -(n[coord[k]]     - R0[coord[k]]) * ds[coord[k]];
   }
 }
 else /// in this case there is no preferred entry ray, everyone starts inside the ROI
   for(k=0; k<3;k++)
   {
     pixel_in[k] = 0;
     n[k] = (int) floor(R0[k]);
     if(dn[k]==1) pixel_out[k] =  (n[k] + 1 - R0[k]) * ds[k];
     else         pixel_out[k] = -(n[k]     - R0[k]) * ds[k];
 }

 if(pixel_out[coord[1]] > pixel_out[coord[2]])
 { ISWAP(coord[1], coord[2]); }

// fprintf(stderr, "(pixel_in=%f %f %f  /// pixel_out = %f  %f %f\n ", pixel_in[0] ,pixel_in[1] ,pixel_in[2], pixel_out[0],  pixel_out[1],  pixel_out[2]);

 if(pixel_out[coord[1]] < pixel_in[coord[0]])   return 0;


 //fprintf(stderr, "\nIn BP one ray, dn = (%d, %d, %d), ds = (%g, %g, %g)\n", dn[0], dn[1], dn[2], ds[0], ds[1], ds[2]);
// ----------------- trace the ray from s_in to exit:  Quadrant independent --------------------
  ray_value = 0.;
  while((n[0] >= 0)  &&  (n[0] < N[0])  && (n[1] >= 0)  &&  (n[1] < N[1]) && (n[2] >= 0)  &&  (n[2] < N[2]))
  {
    nz = n[2];     ny = n[1];     nx = n[0];
    if(pixel_out[coord[0]] < pixel_out[coord[1]])  //  out0 < out1 =< out2  :: only change coord0 stuff
    {
   //   fprintf(stderr, "@1@");
      path_through_pixel = pixel_out[coord[0]] - pixel_in[coord[0]];
      if(path_through_pixel < 0) nrerror("integration error 0");
      pixel_in[coord[0]]  = pixel_out[coord[0]];
      pixel_out[coord[0]] += ds[coord[0]];
      n[coord[0]] += dn[coord[0]];
    }
    else   // out1 >= out0 :: have to reshuffle possibly
    {

      path_through_pixel = pixel_out[coord[1]] - pixel_in[coord[0]]; // this is always true if out[coord[0]] is not the smallest, then out[coord[1]] is the smallest!
 //     fprintf(stderr, "@2@ path = %f\n",path_through_pixel);
      if(path_through_pixel < 0)        nrerror("integration error 1");
// here, we need to reshuffle coordinates so pixel_out[coord[2]] >= pixel+out[coord[1]]
      if(pixel_out[coord[1]] < pixel_out[coord[0]]  && pixel_out[coord[0]] <= pixel_out[coord[2]])  // out1 < out0 <= out2: single swap, only coord1 gets updated
      {
 //     fprintf(stderr, "@3@");
        n[coord[1]] += dn[coord[1]];
        pixel_in[coord[1]] = pixel_out[coord[1]];
        pixel_out[coord[1]] += ds[coord[1]];
        ISWAP(coord[0], coord[1]);
      }
      else if(pixel_out[coord[1]] < pixel_out[coord[0]]  && pixel_out[coord[2]] < pixel_out[coord[0]])  // out1 < out2 < out0: permute swap, only coord1 gets updated
      {
 //     fprintf(stderr, "@4@");
        n[coord[1]] += dn[coord[1]];
        pixel_in[coord[1]] = pixel_out[coord[1]];
        pixel_out[coord[1]] += ds[coord[1]];
        ISWAP3(coord[0], coord[1], coord[2]);
      }
      else if(pixel_out[coord[1]] == pixel_out[coord[0]]  && pixel_out[coord[1]] < pixel_out[coord[2]])  // out1 = out0 < out2: coord0 and coord1 get updated, swap 2 and 1 if needed after update
      {
 //       fprintf(stderr, "@5@");
        n[coord[0]] += dn[coord[0]];
        n[coord[1]] += dn[coord[1]];
        pixel_in[coord[0]] = pixel_out[coord[0]];
        pixel_in[coord[1]] = pixel_out[coord[1]];
        pixel_out[coord[0]] += ds[coord[0]];
        pixel_out[coord[1]] += ds[coord[1]];
        if(pixel_out[coord[2]] < pixel_out[coord[1]])  // shuffle is needed only if pixel_out[coord[1]] BECAME larger than pixel_out[coord[2]]
        {ISWAP(coord[1], coord[2]);}

      }
      else if(pixel_out[coord[1]] == pixel_out[coord[0]]  && pixel_out[coord[0]] == pixel_out[coord[2]])  // out1 = out0 == out2: no swap, all coord-s get updated
      {
//        fprintf(stderr, "@6@");
        n[coord[0]] += dn[coord[0]];
        n[coord[1]] += dn[coord[1]];
        n[coord[2]] += dn[coord[2]];
        pixel_in[coord[0]] = pixel_out[coord[0]];
        pixel_in[coord[1]] = pixel_out[coord[1]];
        pixel_in[coord[2]] = pixel_out[coord[2]];
        pixel_out[coord[0]] += ds[coord[0]];
        pixel_out[coord[1]] += ds[coord[1]];
        pixel_out[coord[2]] += ds[coord[2]];
        if(pixel_out[coord[2]] < pixel_out[coord[1]])  // shuffle is needed only if pixel_out[coord[1]] BECAME larger than pixel_out[coord[2]]
        {ISWAP(coord[1], coord[2]);}
      }
      else
      {
        fprintf(stderr, "k = %d, n = (%d  %d  %d), coord = (%d %d %d)\n", k, n[0], n[1], n[2], coord[0], coord[1], coord[2]);
        fprintf(stderr, "pixel_in = (%f  %f  %f)\n", pixel_in[coord[0]], pixel_in[coord[1]], pixel_in[coord[2]]);
        fprintf(stderr, "pixel_out = (%f  %f  %f)\n", pixel_out[coord[0]], pixel_out[coord[1]], pixel_out[coord[2]]);
        fprintf(stderr, "ds =        (%f, %f, %f)\n", ds[coord[0]]*dn[coord[0]], ds[coord[1]]*dn[coord[1]], ds[coord[2]]*dn[coord[2]]);
        fprintf(stderr, "path = %f\n", path_through_pixel);
        nrerror("Error!");
      }
    }
// fprintf(stdout, "%d %d %d   %g \n", nx, ny, nz, path_through_pixel);
    ray_value +=  volume[nz][ny][nx]* path_through_pixel;
  } // end while loop
  return(ray_value);
}
// --- done integrating
#endif
