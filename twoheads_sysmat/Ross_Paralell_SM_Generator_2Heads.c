#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nrutil.h"

#define PARM_NUM 13
#define MAX_STR_LEN		1024
#define TWO_HEADS
// #define ATTENUATION_YES //-- uncomment to add attenuation, then, though, we need attenuation files

#define EPSILON 0.00001
//#define EPSILON 0.00005  //changed by Haoran on 2016/8/29
static int iswap_arg;
#define ISWAP(a, b) iswap_arg = a;  a = b;  b = iswap_arg
#define ISWAP3(a, b, c) iswap_arg = a;  a = b;  b = c;  c = iswap_arg

int read_paramaters(int parm_num, char escape, float *params, FILE *fptr);
int main();
void BP_one_ray_new_attn(  double backprojection_value, float ***volume, int Nxyz[], double R0[], double directions[], float ***attn_map);
void BP_one_ray_new_noattn(double ray_value, float ***volume, int Nxyz[], double R0[], double directions[]);
int extract_nonzero_entries(float ***volume, int Nx, int Ny, int Nz, int *indices, float *entries);
void point_response_weights_parallel02(int Ndir, float *xdir, float *ydir, float *zdir, float *weight);

// // // // // // // // // // // // // // // // // // // // // // // // // // //
 // // // // // // // // // // // // // // // // // // // // // // // // // // //
// // // // // // // // // // // // // // // // // // // // // // // // // // //
int read_paramaters(int parm_num, char escape, float *params, FILE *fptr)
{
  int i=0;
  char buffer[MAX_STR_LEN],parm[MAX_STR_LEN];
  char *str_pos;

  while (i<parm_num)
  {
    if(fgets (buffer , MAX_STR_LEN , fptr)==NULL)
    	return i;

    str_pos = strchr(buffer,escape);
    /*printf ("string:  %s",buffer);*/
    if((str_pos-buffer+1 >= 0) && (str_pos-buffer+1 <= MAX_STR_LEN))
    {
	    strncpy (parm,buffer,str_pos-buffer+1);
	    parm[str_pos-buffer+1]='\0';
	    params[i] = atof(parm);
	    /*printf ("From: %s   we get: %f    ",parm,params[i]);*/
	    i++;
    }
    /*printf ("--> found at %d \n",str_pos-buffer+1);*/
  }
  return i;
}

// // // // // // // // // // // // // // // // // // // // // // // // // // //
 // // // // // // // // // // // // // // // // // // // // // // // // // // //
// // // // // // // // // // // // // // // // // // // // // // // // // // //
int main(int argc, char **argv)
{
  int Nx, Ny, Nz, Nr, NZ, Nang; // volume and sinogram dimensions
  int Ndir; // no of
  int Nvolume, n, m, k;
  int N_pixel_subdivision, N3[3];
  int nz, nr, nang;
  int Ndet, Nsino;
  int index;
  size_t sizeof_float, sizeof_int;
  double d_angle, one_over_dx;
  double directions[3], R0[3];
  float *xdir, *ydir, *zdir, *weight, *xdir_tmp, *ydir_tmp;
  float *r_pixel_subdivision, *z_pixel_subdivision;
  double backprojection_value, pixel_subdivision_factor, t, r_center, x_center, y_center;
  float ***volume;

  int sinogram_index, no_entries, *indices;
  float *entries;
  float *params;
  double theta, cos_theta, sin_theta;

  FILE *fid;  // file id for reading attn map and writing system matrix

   if( (fid=fopen(argv[2], "rb"))==NULL)
   {
     fprintf(stderr, "Failed to open parameters file %s\n", argv[2]);
     exit(1);
   }
   else
   {
     params = (float *) malloc(sizeof(float)*PARM_NUM);

 	if(index=read_paramaters(PARM_NUM, ':', params, fid) != PARM_NUM)
 	{
 		fprintf(stderr, "Expected %d paramaters in file %s, only read %d paramaters\n",PARM_NUM, argv[2],index);
 		exit(1);
 	}
 	//convert degrees to radians
    fprintf(stdout, "Number of view: %f. Starting angle: %f. Rotation: %f", params[6], params[7], params[9]);
    params[7] = params[7] * M_PI / 180;
 	params[8] = params[8] * M_PI / 180;
 	params[9] = params[9] * M_PI / 180;
    fclose(fid);

   }

   N_pixel_subdivision = (int)params[12];
   Ndir = (int)params[11];

#ifdef ATTENUATION_YES
  float ***attn_tensor;
#endif

  // ---- end declacations -----


// --- initial assignments ---
  // Reconstructed volume dimensions
  Nx = (int)params[2];
  Ny = (int)params[2];
  Nz = (int)params[3];
  Nvolume = Nx * Ny * Nz;
  // sinogram/detector dimensions
  Nr = (int)params[0];
  NZ = (int)params[1];
  Nang = (int)params[6];
  Ndet = Nr * NZ;
  Nsino = Nang * Ndet * 2;
  
  // variables to be used often
  sizeof_float = sizeof(float);
  sizeof_int = sizeof(int);
  // --- more helper variable
  N3[0] = Nx;  N3[1] = Ny;  N3[2] = Nz;
  one_over_dx = 1./ params[4];
  d_angle = params[9] / ((double) Nang);
  r_center = Nr * 0.5;
  x_center = Nx * 0.5;
  y_center = Ny * 0.5;

  //////////////////////////////////////////////
  //////////////////////////////////////////////
  
  // ---Check input variables----
#ifdef ATTENUATION_YES
  if(argc != 4)
  { fprintf(stderr, "Usage: %c a.out output_sysmat.fname params.fname attn_matrix.fname\n", 37); return 1;  }
  // read in attenuation map array, has to be Nx x Ny x Nz, 4-byte floats, in 1/1000 of attn coeff 1/cm
  if( (fid=fopen(argv[3], "rb"))==NULL)
  { fprintf(stderr, "Failed to open attenuation matrix file %s\n", argv[3]); exit(1);}
  //
  attn_tensor =  f3tensor(0, Nz-1, 0, Ny-1, 0, Nx-1);
  fread(&attn_tensor[0][0][0], sizeof(float), Nvolume, fid);
  fclose(fid);
  // scale attn_tensor, scaling factor depends on the format of the data
  // scaling includes the pixel length. This scaling is here because attn map
  // is typically saved using unsigned shorts or whatever int variables which are
  // round(1000 * what_we_need)
  t = 0.001 * params[5] * 0.15; // times 0.15 if the \mu units are in Hounsfeld
  for(n=0; n<Nz; n++) for(m=0; m<Ny; m++) for(k=0; k<Nx; k++) attn_tensor[n][m][k] *= t;
#else
  if(argc != 3)
  { fprintf(stderr, "Usage: %c a.out output_sysmat.fname params.fname\n", 37); return 1;  }
#endif

  // --- Pixel subdivision ---
  if(params[12] < 1 || params[12]> 5){
    fprintf(stderr, "N_pixel_subdivisions should be 1, 2x2, 3x3, 4x4, or 5x5, intstead its %d^2\n",  params[12]);
    fprintf(stderr, "Change N_pixel_subdivisions in the #define section\n");   exit(1);}
  n = N_pixel_subdivision*N_pixel_subdivision;
  r_pixel_subdivision = (float *) malloc(n * sizeof_float);
  z_pixel_subdivision = (float *) malloc(n * sizeof_float);
  t = 1./((double) N_pixel_subdivision);
  for(n=0; n<N_pixel_subdivision; n++)
    for(k=0; k<N_pixel_subdivision; k++){
      m = k + N_pixel_subdivision*n;
      r_pixel_subdivision[m] = t*(k + 0.5);
      z_pixel_subdivision[m] = t*(n + 0.5);
    }

  //
  N_pixel_subdivision*=N_pixel_subdivision;
  pixel_subdivision_factor = 1. / ((double) N_pixel_subdivision);
  fprintf(stderr, "NO_PIXEL_SUBDIVISION = %d, N_subdivisions = %d, weight = %f\n", (int)params[12], (int)N_pixel_subdivision, pixel_subdivision_factor);
  //

  // --- Ray propagation directions ---
  xdir   = (float *) malloc(sizeof_float * Ndir);
  ydir   = (float *) malloc(sizeof_float * Ndir);
  zdir   = (float *) malloc(sizeof_float * Ndir);
  weight   = (float *) malloc(sizeof_float * Ndir);

  xdir_tmp   = (float *) malloc(sizeof_float * Ndir);
  ydir_tmp   = (float *) malloc(sizeof_float * Ndir);
  point_response_weights_parallel02(Ndir, xdir, ydir, zdir, weight); // gives an error if Ndir is wrong
//fprintf(stderr, "%d ray weights:\n", Ndir);  for(n=0; n<Ndir; n++)    fprintf(stderr, " %g   %g   %g  / %g\n", xdir[n], ydir[n], zdir[n], weight[n]);  exit(23);



  // allocate memory for the pieces that should get written into the file
  // the largest possible size they can get is the size of the volume
  indices = (int   *) malloc(sizeof_int   * Nvolume);
  entries = (float *) malloc(sizeof_float * Nvolume);

  // declare volume and sinogram mask, compute sinogram mask
  volume = f3tensor(0, Nz-1, 0, Ny-1, 0, Nx-1);

  // open system matrix file and GET READY TO ROCK ;-)
  if( (fid=fopen(argv[1], "wb")) == NULL){
    fprintf(stderr, "Failed to open file %s for writing system matrix!\n", argv[1]);
    exit(1);}
  else fprintf(stderr, "Opened file %s for writing system matrix\n", argv[1]);

  // write sinogram length and volume size
  //Nsino = Ndet * Nang;
  fwrite(&Nsino, sizeof_int,   1, fid);
  fwrite(&Nvolume, sizeof_int, 1, fid);

  //db("2");

// --------Assume one head for now, second head --------
  fprintf(stderr, "\nHead 1, %d views:\n", Nang);
  // --- angular loop ---
  for(nang=0; nang<Nang; nang++)
  {
    fprintf(stderr, "<%d>", nang+1);
	// compute sin and cos of the current rotation angle and the propagations angles that take into account geometric correction
    theta = (nang * d_angle + params[7]); // for  head 2, add M_PI or M_PI*0.5, for H and L mode respectively
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    for(k=0; k<Ndir; k++){
	    xdir_tmp[k] =  xdir[k]*cos_theta - ydir[k]*sin_theta;
	    ydir_tmp[k] =  xdir[k]*sin_theta + ydir[k]*cos_theta;}
 	//
    for(nz=0; nz<NZ; nz++)
      for(nr=0; nr<Nr; nr++)
      {
        sinogram_index = nr + Nr*(nz+NZ*nang); // head no 1, for head 2, add Nr*NZ*Nang
	    // set volume index to zeros
        memset(&volume[0][0][0], 0, Nvolume*sizeof_float);
        for(m=0; m<N_pixel_subdivision; m++)
        {
          // compute the starting point of the ray (
          R0[0] = params[10]*cos_theta - (nr-r_center+ r_pixel_subdivision[m])*sin_theta + x_center;
  		    R0[1] = params[10]*sin_theta + (nr-r_center+ r_pixel_subdivision[m])*cos_theta + y_center;
		      R0[2] = nz + z_pixel_subdivision[m];
          //if(nz == 23 && nr == 0 && m==0) fprintf(stderr, "%g     %g\n", R0[0], R0[1]);
		  //
  		    for(k=0; k<Ndir; k++)
		      {
//fprintf(stderr, "[%d, %d] ", m, k);
		        directions[0] = xdir_tmp[k];
//fprintf(stderr, " , ");
            directions[1] = ydir_tmp[k];
//fprintf(stderr, " ; ");
			      directions[2] = zdir[k];
//fprintf(stderr, " ^ ");
            backprojection_value = pixel_subdivision_factor*weight[k];
//fprintf(stderr, " a ");
            #ifdef ATTENUATION_YES
            BP_one_ray_new_attn(backprojection_value,   volume, N3, R0, directions, attn_tensor);
            #else
            BP_one_ray_new_noattn(backprojection_value, volume, N3, R0, directions);
            #endif
//fprintf(stderr, " b ");
          }
        }
//fprintf(stderr, ">\n");
		// at this point, we are done with backprojecting detector pixel (nang, nr, nz) into the volume
        no_entries = extract_nonzero_entries(volume, Nx, Ny, Nz, indices, entries);
        if(no_entries<0) 	          nrerror("no_entries should be non-negative!");
        else if(no_entries>Nvolume) nrerror("no entries larger than Nvolume -- this should not happen");
        fwrite(&sinogram_index, sizeof_int, 1, fid);
        fwrite(&no_entries,     sizeof_int, 1, fid);
        if(no_entries==0) continue;
        fwrite(indices, sizeof_int,   no_entries, fid);
        fwrite(entries, sizeof_float, no_entries, fid);
      } // end of Nr loop
  } // end angular loop

  ///============================= head2 =======================
#ifdef TWO_HEADS

  fprintf(stderr, "\nHead 2, %d views:\n", Nang);
  // --- angular loop ---
  for(nang=0; nang<Nang; nang++)
  {
    fprintf(stderr, "<%d>", nang+1);
	// compute sin and cos of the current rotation angle and the propagations angles that take into account geometric correction
    theta = nang * d_angle + params[8];
    cos_theta = cos(theta);
    sin_theta = sin(theta);
    for(k=0; k<Ndir; k++){
	    xdir_tmp[k] =  xdir[k]*cos_theta - ydir[k]*sin_theta;
	    ydir_tmp[k] =  xdir[k]*sin_theta + ydir[k]*cos_theta;}
 	//
    for(nz=0; nz<NZ; nz++)
      for(nr=0; nr<Nr; nr++)
      {
        //sinogram_index = nr + Nr*(nz+NZ*nang); // head no 1, for head 2, add Nr*NZ*Nang
        sinogram_index = nr + Nr*(nz+NZ*nang) + Nr*NZ*Nang; // Edit 2022/4/27 by Haoran
	    // set volume index to zeros
        memset(&volume[0][0][0], 0, Nvolume*sizeof_float);
        for(m=0; m<N_pixel_subdivision; m++)
        {
          // compute the starting point of the ray (
          R0[0] = params[10]*cos_theta - (nr-r_center+ r_pixel_subdivision[m])*sin_theta + x_center;
  		    R0[1] = params[10]*sin_theta + (nr-r_center+ r_pixel_subdivision[m])*cos_theta + y_center;
		      R0[2] = nz + z_pixel_subdivision[m];
  		    for(k=0; k<Ndir; k++)
		      {
		        directions[0] = xdir_tmp[k];
            directions[1] = ydir_tmp[k];
			      directions[2] = zdir[k];
            backprojection_value = pixel_subdivision_factor*weight[k];
            #ifdef ATTENUATION_YES
            BP_one_ray_new_attn(backprojection_value,   volume, N3, R0, directions, attn_tensor);
            #else
            BP_one_ray_new_noattn(backprojection_value, volume, N3, R0, directions);
            #endif
          }
        }
		// at this point, we are done with backprojecting detector pixel (nang, nr, nz) into the volume
        no_entries = extract_nonzero_entries(volume, Nx, Ny, Nz, indices, entries);
        if(no_entries<0) 	          nrerror("no_entries should be non-negative!");
        else if(no_entries>Nvolume) nrerror("no entries larger than Nvolume -- this should not happen");
        fwrite(&sinogram_index, sizeof_int, 1, fid);
        fwrite(&no_entries,     sizeof_int, 1, fid);
        if(no_entries==0) continue;
        fwrite(indices, sizeof_int,   no_entries, fid);
        fwrite(entries, sizeof_float, no_entries, fid);
      } // end of Nr loop
  } // end angular loop
#endif
  fclose(fid); // done creating system matrix file
  fprintf(stderr, "\nDone, cleaning up memory...\n");
  //
  free_f3tensor(volume, 0, Nz-1, 0, Ny-1, 0, Nx-1);
  free(indices);
  free(entries);
  free(xdir);
  free(ydir);
  free(xdir_tmp);
  free(ydir_tmp);
  free(zdir);
  free(weight);
  free(r_pixel_subdivision);
  free(z_pixel_subdivision);
  //
  return 0;
}


// // // // // // // // // // // // // // // // // // // // // // // // // // //
 // // // // // // // // // // // // // // // // // // // // // // // // // // //
// // // // // // // // // // // // // // // // // // // // // // // // // // //
int extract_nonzero_entries(float ***volume, int Nx, int Ny, int Nz, int *indices, float *entries)
{
  int nx, ny, nz, index, total=0;

  for(nz=0; nz<Nz; nz++)
    for(ny=0; ny<Ny; ny++)
      for(nx=0; nx<Nx; nx++)
        if(volume[nz][ny][nx] > 1.0e-013)
        {
          indices[total] = nx + Nx*(ny+Ny*nz);
          entries[total] = volume[nz][ny][nx];
          total++;
        }
  return(total);
}

// ------------------------ No attenuation,  backproject  ----------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------
void BP_one_ray_new_noattn(double ray_value, float ***volume, int N[], double R0[], double directions[])
// backproject_one_ray_noattn_nomask(float ray_value, float *A, const int *N, const double *R0, const double *directions, const double *ds, const int *dn)
{
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
     dn[k] = 1;//0; ------------------- !!!!!!!!!!!!!!!!!
     ds[k] = 1.0e+10;
   }

// Initial assignment of s_in and s_out -- should change as we go through the thing
 s_in = 0.;  s_out = 8388606.;
// fprintf(stderr, "(%f,%f,%f)/(%f,%f,%f) by (%f, %f, %f)\n", R0[0], R0[1], R0[2], directions[0], directions[1], directions[2], ds[0]*dn[0], ds[1]*dn[1], ds[2]*dn[2] );

// ---------------------- all components -----------------
 for(k=0; k<3; k++)
 {
   if(dn[k] > 0) // positive propagation in the direction k
   {
     if(R0[k] >= N[k]) return; // we know that we miss the volume
     n[k] = 0;
     pixel_in[k]  = -R0[k] * ds[k]; // The minus sign is OK, we expect R0 to be negative unless ray enters from the side
     pixel_out[k] = (N[k] - R0[k]) * ds[k];
     if(pixel_in[k] > s_in) // assign this as entrance coordinate -- negative propagation case
     {
       s_in = pixel_in[k];
       ISWAP(coord[0], coord[k]);
       n[coord[k]] = 0;
     }
   }
   else                 // negative propagation
   {
     if(R0[k] <= 0.) return; // we know that we miss the volume
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
  if(s_in >= s_out) return; // ------------------ return zero if in is larger than out

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

 if(pixel_out[coord[1]] < pixel_in[coord[0]])   return;


 //fprintf(stderr, "\nIn BP one ray, dn = (%d, %d, %d), ds = (%g, %g, %g)\n", dn[0], dn[1], dn[2], ds[0], ds[1], ds[2]);
// ----------------- trace the ray from s_in to exit:  Quadrant independent --------------------
  while((n[0] >= 0)  &&  (n[0] < N[0])  && (n[1] >= 0)  &&  (n[1] < N[1]) && (n[2] >= 0)  &&  (n[2] < N[2]))
  {

    nz = n[2];     ny = n[1];     nx = n[0];
    if(pixel_out[coord[0]] < pixel_out[coord[1]])  //  out0 < out1 =< out2  :: only change coord0 stuff
    {
   //   fprintf(stderr, "@1@");
      path_through_pixel = pixel_out[coord[0]] - pixel_in[coord[0]];
      if(path_through_pixel < 0) nrerror("system matrix attn 0");
      pixel_in[coord[0]]  = pixel_out[coord[0]];
      pixel_out[coord[0]] += ds[coord[0]];
      n[coord[0]] += dn[coord[0]];
    }
    else   // out1 >= out0 :: have to reshuffle possibly
    {

      path_through_pixel = pixel_out[coord[1]] - pixel_in[coord[0]]; // this is always true if out[coord[0]] is not the smallest, then out[coord[1]] is the smallest!
 //     fprintf(stderr, "@2@ path = %f\n",path_through_pixel);
      if(path_through_pixel < 0)        nrerror("system matrix attn 1");
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
    } // end "else": shuffling and updating
///*-------------------------------------------------------------------*/ fprintf(stdout, "%d %d %d   %g \n", nx, ny, nz, path_through_pixel);
    volume[nz][ny][nx] += ray_value * path_through_pixel;
  } // end while loop
  return;
}




// ------------------------ Yes attenuation,  backproject  ----------------------------------------------------------
// -----------------------------------------------------------------------------------------------------------------
void BP_one_ray_new_attn(double ray_value, float ***volume, int N[], double R0[], double directions[], float ***attn_map)
// backproject_one_ray_noattn_nomask(float ray_value, float *A, const int *N, const double *R0, const double *directions, const double *ds, const int *dn)
{
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
 double path_through_pixel, exp_mu_ds;
 int dn[3];
 double ds[3];

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
     dn[k] = 1;//0; ------------------- !!!!!!!!!!!!!!!!!
     ds[k] = 1.0e+10;
   }

// Initial assignment of s_in and s_out -- should change as we go through the thing
 s_in = 0.;  s_out = 8388606.;
// fprintf(stderr, "(%f,%f,%f)/(%f,%f,%f) by (%f, %f, %f)\n", R0[0], R0[1], R0[2], directions[0], directions[1], directions[2], ds[0]*dn[0], ds[1]*dn[1], ds[2]*dn[2] );

// ---------------------- all components -----------------
 for(k=0; k<3; k++)
 {
   if(dn[k] > 0) // positive propagation in the direction k
   {
     if(R0[k] >= N[k]) return; // we know that we miss the volume
     n[k] = 0;
     pixel_in[k]  = -R0[k] * ds[k]; // The minus sign is OK, we expect R0 to be negative unless ray enters from the side
     pixel_out[k] = (N[k] - R0[k]) * ds[k];
     if(pixel_in[k] > s_in) // assign this as entrance coordinate -- negative propagation case
     {
       s_in = pixel_in[k];
       ISWAP(coord[0], coord[k]);
       n[coord[k]] = 0;
     }
   }
   else                 // negative propagation
   {
     if(R0[k] <= 0.) return; // we know that we miss the volume
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
  if(s_in >= s_out) return; // ------------------ return zero if in is larger than out

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

 if(pixel_out[coord[1]] < pixel_in[coord[0]])   return;


 //fprintf(stderr, "\nIn BP one ray, dn = (%d, %d, %d), ds = (%g, %g, %g)\n", dn[0], dn[1], dn[2], ds[0], ds[1], ds[2]);
// ----------------- trace the ray from s_in to exit:  Quadrant independent --------------------
  while((n[0] >= 0)  &&  (n[0] < N[0])  && (n[1] >= 0)  &&  (n[1] < N[1]) && (n[2] >= 0)  &&  (n[2] < N[2]))
  {

    nz = n[2];     ny = n[1];     nx = n[0];
    if(pixel_out[coord[0]] < pixel_out[coord[1]])  //  out0 < out1 =< out2  :: only change coord0 stuff
    {
   //   fprintf(stderr, "@1@");
      path_through_pixel = pixel_out[coord[0]] - pixel_in[coord[0]];
      if(path_through_pixel < 0) nrerror("system matrix attn 0");
      pixel_in[coord[0]]  = pixel_out[coord[0]];
      pixel_out[coord[0]] += ds[coord[0]];
      n[coord[0]] += dn[coord[0]];
    }
    else   // out1 >= out0 :: have to reshuffle possibly
    {

      path_through_pixel = pixel_out[coord[1]] - pixel_in[coord[0]]; // this is always true if out[coord[0]] is not the smallest, then out[coord[1]] is the smallest!
 //     fprintf(stderr, "@2@ path = %f\n",path_through_pixel);
      if(path_through_pixel < 0)        nrerror("system matrix attn 1");
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
    } // end "else": shuffling and updating
    if(attn_map[nz][ny][nx] > EPSILON)
    {
       exp_mu_ds = exp(-attn_map[nz][ny][nx] * path_through_pixel);
       volume[nz][ny][nx] += ray_value * (1.-exp_mu_ds) /  attn_map[nz][ny][nx]; // update the pixel
//       ray_value /= exp_mu_ds; // de-attenuate the value
       ray_value *= exp_mu_ds; // de-attenuate the value
       if(attn_map[nz][ny][nx] > 0.5){ fprintf(stderr, "Unusially large attn map[%d][%d][%d] = %g\n", nz, ny, nx, attn_map[nz][ny][nx]);exit(1);}
    }
    else // no attenuation in the pixel
       volume[nz][ny][nx] += ray_value * path_through_pixel;

  } // end while loop
  return;
}


// ---------------------------------------------------------------------
// ----- computes point response weights: phenomenon of non-parallel ---
// ----- propagation of rays in case of real (non-ideal) collimation ---
// ---------- the parameters here were computed empiricallly -----------
// ---------------------------------------------------------------------
void point_response_weights_parallel02(int Ndir, float *xdir, float *ydir, float *zdir, float *weight) // tau = 0.15 sigma ?
{
  int n;
  double a = 0.0256685, b = 0.208467, dist=1000., tau, tau3_2;

  if(Ndir == 1) // ------------ Zeroth order,          N = 0   ------------
  {
    ydir[0] = zdir[0] = 0.;
    weight[0] = 1.;
  }
  else if(Ndir ==7)  // ------ first order nearest neighbors. N = 1 ------
  {
    tau = a * dist + b;  // width of the Gaussian in plane distance dist away
    tau *= 0.2;// distance to the nearest point
    tau3_2 = tau * sqrt(3.) * 0.5; // this is what I denoted as SIXTY times tau
// central point
    ydir[0] = zdir[0] = 0.;
// six surrounding points
/**/ydir[1] =  ydir[6] =  tau3_2;
    ydir[3] =  ydir[4] = -tau3_2;
    ydir[2] =  ydir[5] =  0.;
/**/zdir[1] =  zdir[3] =  tau * 0.5;
    zdir[4] =  zdir[6] = -tau * 0.5;
    zdir[2] =  tau;
    zdir[5] = -tau;
/**/weight[0] = 0.14782541493764;
    for(n=1; n<7; n++) weight[n] = 0.14202909751039;
  }
  else if(Ndir == 13) // ------ second order nearest neighbors. N = 2 ------
  {
    tau = a * dist + b;  // width of the Gaussian in plane distance dist away
    tau *= 0.2;// distance to the nearest point
    tau3_2 = tau * sqrt(3.) * 0.5; // this is what I denoted as SIXTY times tau
// central point
    ydir[0] = zdir[0] = 0.;
// six surrounding points
/**/ydir[1] =  ydir[6] =  tau3_2;
    ydir[3] =  ydir[4] = -tau3_2;
    ydir[2] =  ydir[5] =  0.;
/**/zdir[1] =  zdir[3] =  tau * 0.5;
    zdir[4] =  zdir[6] = -tau * 0.5;
    zdir[2] =  tau;
    zdir[5] = -tau;
// middles of the second hexagon
/**/ydir[7] = 2. * tau3_2;
    ydir[8] = ydir[12] =  tau3_2;
    ydir[9] = ydir[11] = -tau3_2;
    ydir[10] = -2. * tau3_2;
/**/zdir[7] = zdir[10] = 0.;
    zdir[8] = zdir[9] = 1.5 * tau;
    zdir[11] = zdir[12] = -1.5 * tau;
/**/weight[0] = 0.08273858593331;
    for(n=1; n<7; n++) weight[n] = 0.07949435957512;
    for(n=7; n<13; n++) weight[n] = 0.07338254276933;
  }
  else if(Ndir == 19) // ------ third order nearest neighbors. N = 3 ------
  {
    tau = a * dist + b;  // width of the Gaussian in plane distance dist away
    tau *= 0.2;// distance to the nearest point
    tau3_2 = tau * sqrt(3.) * 0.5; // this is what I denoted as SIXTY times tau
// central point
    ydir[0] = zdir[0] = 0.;
// six surrounding points
/**/ydir[1] =  ydir[6] =  tau3_2;
    ydir[3] =  ydir[4] = -tau3_2;
    ydir[2] =  ydir[5] =  0.;
/**/zdir[1] =  zdir[3] =  tau * 0.5;
    zdir[4] =  zdir[6] = -tau * 0.5;
    zdir[2] =  tau;
    zdir[5] = -tau;
// middles of the second hexagon
/**/ydir[7] = 2. * tau3_2;
    ydir[8] = ydir[12] =  tau3_2;
    ydir[9] = ydir[11] = -tau3_2;
    ydir[10] = -2. * tau3_2;
/**/zdir[7] = zdir[10] = 0.;
    zdir[8] = zdir[9] = 1.5 * tau;
    zdir[11] = zdir[12] = -1.5 * tau;
// vertices of the second hexagon
/**/ydir[13] = ydir[18] =  2. * tau3_2;
    ydir[14] = ydir[17] = 0.;
    ydir[15] = ydir[16] = -2. * tau3_2;
/**/zdir[13] = zdir[15] =  tau;
    zdir[16] = zdir[18] = -tau;
    zdir[14] =  2. * tau;
    zdir[17] = -2. * tau;
/**/weight[0] = 0.05814250289281;
    for(n=1; n<7; n++) weight[n] =   0.05586270274529;
    for(n=7; n<13; n++) weight[n] =   0.05156777405752;
    for(n=13; n<19; n++) weight[n] =   0.04954577271506;
  }
  else if(Ndir == 31) // ------ fourth order nearest neighbors. N = 4 ------
  {
    tau = a * dist + b;  // width of the Gaussian in plane distance dist away
    tau *= 0.2;// distance to the nearest point
    tau3_2 = tau * sqrt(3.) * 0.5; // this is what I denoted as SIXTY times tau
// central point
    ydir[0] = zdir[0] = 0.;
// six surrounding points
/**/ydir[1] =  ydir[6] =  tau3_2;
    ydir[3] =  ydir[4] = -tau3_2;
    ydir[2] =  ydir[5] =  0.;
/**/zdir[1] =  zdir[3] =  tau * 0.5;
    zdir[4] =  zdir[6] = -tau * 0.5;
    zdir[2] =  tau;
    zdir[5] = -tau;
// middles of the second hexagon
/**/ydir[7] = 2. * tau3_2;
    ydir[8] = ydir[12] =  tau3_2;
    ydir[9] = ydir[11] = -tau3_2;
    ydir[10] = -2. * tau3_2;
/**/zdir[7] = zdir[10] = 0.;
    zdir[8] = zdir[9] = 1.5 * tau;
    zdir[11] = zdir[12] = -1.5 * tau;
// vertices of the second hexagon
/**/ydir[13] = ydir[18] =  2. * tau3_2;
    ydir[14] = ydir[17] = 0.;
    ydir[15] = ydir[16] = -2. * tau3_2;
/**/zdir[13] = zdir[15] =  tau;
    zdir[16] = zdir[18] = -tau;
    zdir[14] =  2. * tau;
    zdir[17] = -2. * tau;
// sides of the third hexagon
/**/ydir[19] = ydir[30] = 3. * tau3_2;
    ydir[20] = ydir[29] = 2. * tau3_2;
    ydir[21] = ydir[28] = tau3_2;
    ydir[22] = ydir[27] = -tau3_2;
    ydir[23] = ydir[26] = -2. * tau3_2;
    ydir[24] = ydir[25] = -3. * tau3_2;
/**/zdir[19] = zdir[24] = tau * 0.5;
    zdir[20] = zdir[23] = tau * 2.;
    zdir[21] = zdir[22] = tau * 2.5;
    zdir[30] = zdir[25] = -tau * 0.5;
    zdir[29] = zdir[26] = -tau * 2.;
    zdir[27] = zdir[28] = -tau * 2.5;

/**/weight[0] = 0.03806836996690;
    for(n=1; n<7; n++) weight[n] =   0.03657568782994;
    for(n=7; n<13; n++) weight[n] =   0.03376361531615;
    for(n=13; n<19; n++) weight[n] =   0.03243972502336;
    for(n=19; n<31; n++) weight[n] =   0.02877145508470;
  }
  else if(Ndir == 37) // ------ fourth order nearest neighbors. N = 5 ------
  {
    tau = a * dist + b;  // width of the Gaussian in plane distance dist away
    tau *= 0.2;// distance to the nearest point
    tau3_2 = tau * sqrt(3.) * 0.5; // this is what I denoted as SIXTY times tau
// central point
    ydir[0] = zdir[0] = 0.;
// six surrounding points
/**/ydir[1] =  ydir[6] =  tau3_2;
    ydir[3] =  ydir[4] = -tau3_2;
    ydir[2] =  ydir[5] =  0.;
/**/zdir[1] =  zdir[3] =  tau * 0.5;
    zdir[4] =  zdir[6] = -tau * 0.5;
    zdir[2] =  tau;
    zdir[5] = -tau;
// middles of the second hexagon
/**/ydir[7] = 2. * tau3_2;
    ydir[8] = ydir[12] =  tau3_2;
    ydir[9] = ydir[11] = -tau3_2;
    ydir[10] = -2. * tau3_2;
/**/zdir[7] = zdir[10] = 0.;
    zdir[8] = zdir[9] = 1.5 * tau;
    zdir[11] = zdir[12] = -1.5 * tau;
// vertices of the second hexagon
/**/ydir[13] = ydir[18] =  2. * tau3_2;
    ydir[14] = ydir[17] = 0.;
    ydir[15] = ydir[16] = -2. * tau3_2;
/**/zdir[13] = zdir[15] =  tau;
    zdir[16] = zdir[18] = -tau;
    zdir[14] =  2. * tau;
    zdir[17] = -2. * tau;
// sides of the third hexagon
/**/ydir[19] = ydir[30] = 3. * tau3_2;
    ydir[20] = ydir[29] = 2. * tau3_2;
    ydir[21] = ydir[28] = tau3_2;
    ydir[22] = ydir[27] = -tau3_2;
    ydir[23] = ydir[26] = -2. * tau3_2;
    ydir[24] = ydir[25] = -3. * tau3_2;
/**/zdir[19] = zdir[24] = tau * 0.5;
    zdir[20] = zdir[23] = tau * 2.;
    zdir[21] = zdir[22] = tau * 2.5;
    zdir[30] = zdir[25] = -tau * 0.5;
    zdir[29] = zdir[26] = -tau * 2.;
    zdir[27] = zdir[28] = -tau * 2.5;
// vertices of the third hexagon
/**/ydir[31] = ydir[36] =  3. * tau3_2;
    ydir[32] = ydir[25] = 0.;
    ydir[33] = ydir[34] = -2. * tau3_2;
/**/zdir[31] = zdir[33] =  1.5 * tau;
    zdir[36] = zdir[34] = -1.5 * tau;
    zdir[32] =  3. * tau;
    zdir[35] = -3. * tau;
/**/weight[0] = 0.03283577842739;
    for(n=1; n<7; n++) weight[n] = 0.03154826913938;
    for(n=7; n<13; n++) weight[n] = 0.02912272294277;
    for(n=13; n<19; n++) weight[n] = 0.02798080464277;
    for(n=19; n<31; n++) weight[n] = 0.02481674747347;
    for(n=31; n<37; n++) weight[n] = 0.02290874525691;
  }
  else
  {
    fprintf(stderr, "\nNdir = %d, we did not expect this number !!!\n", Ndir);
    exit(1);
  }

 for(n=0; n<Ndir; n++)
 {
  a = sqrt(dist*dist +  ydir[n]*ydir[n] + zdir[n]*zdir[n]);
  xdir[n] = -dist / a;
  ydir[n] /= a;
  zdir[n] /= a;
 }
}

