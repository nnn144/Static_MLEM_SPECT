#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#pragma warning(disable:4996)

#define RECON_3D


// ----- This block defines implementation-specific parameters -----

//=======================================================================================================
// ==============================     Select sinogram and volume Data types    ==========================
//=======================================================================================================
// input sinogram can be either 16-bit unsigned integer (unsigned short) or float
#define FLOAT_INPUT_SINOGRAM  // uncomment this, if input sinogram is float

// Volume is best output as 4-bit float
#define FLOAT_INPUT_VOLUME

//=======================================================================================================
// ==============================     Select reconstruction algorithm: MLEM or CG    ====================
//=======================================================================================================

// Uncomment this if you want to reconstruct with MLEM Otherwise reconstruction with Conjugate Gradient
// MLEM is much faster than CG and provides better results
#define RECON_MLEM

//=======================================================================================================
// ======================================     Select number of iterations    ============================
//=======================================================================================================

// No of iterations
// Typically MLEM needs 20 iterations, while CG needs 70 or more
int N_ITERATIONS_MLEM;        // MLEM max iterations 20
#define N_ITERATIONS_CG 30          // CG max iterations 50

//=======================================================================================================
// ==============================     Set regularization parameters    ==================================
//=======================================================================================================

// Uncomment this if you want to reconstruct with TV smoothness (this is for static)
//#define STATIC_RECON_WITH_TV


// Ucomment this if you want to select the regularizarion parameter dynamically
//#define DYNAMIC_LAMBDA_SELECTION

// Provide the needed values for the regularization parameter.
// If the DYNAMIC_LAMBDA_SELECTION is commented, these value will be used instead
#define DEFAULT_SMOOTHNESS_LAMBDA 0.0//.7       // Total variation lambda
#define DEFAULT_MIX_LAMBDA 0.0//10                // coefficient mix lambda
#define DEFAULT_TB_SMOOTHNESS_LAMBDA 0.0//20      // Curve smoothing lambda

//=======================================================================================================
// ======================================     Do not change these parameters    =========================
//=======================================================================================================
// This defines how small zero is, and what is one divided by zero
#define ONE_OVER_EPSILON 3.0e+13
#define EPSILON 3.0e-13

size_t size_float= sizeof(float);
size_t size_int= sizeof(int);

// a structure to store all passed parameters. this is just for simplisty.
struct PARAMETERS
{
    float *sinogram;           // Sinogram array
    float *basis;              // time baiss functions, bsplines,factors, TACs, or any temporal curves
    float *coefficients;       // coefficients array
    char *sysmat_fname;        // name of system matrix
    float *StaticMask;         // mask created out of segmented static reconstruction
    float *Mask;
    
    int singleProjectionSize;  // size of projection/number of pins on head
    int basisNum;              // number of time basis/factors/TACs
    int basisLen;              // length of time basis/factors/TACs
    int Nvolume;               // size of imaged volume
    int Ncoefficients;         // size of coefficients = Nvolume * basisNum
    int Nsinogram;             // size of sinogram
    int No_Heads;              // number of detectors/heads use for acquiring the data
    int Nx;                    // dimenrsion of imaged volume in x direction
    int Ny;                    // dimenrsion of imaged volume in y direction
    int Nz;                    // dimenrsion of imaged volume in z direction
    
    int flag;                   // to turn on/off regularization
    
    float TB_SMOOTHNESS_LAMBDA; // Time basis/factors/TACs Smoothness parameter
    float SMOOTHNESS_LAMBDA;    // coefficient nearest niebghour (total variation)
    float MIX_LAMBDA;           // coefficient mix
    
    float gamma_0;              // parameter for dynamic lambda estimation
    
    float chi;                  // container to store the current value of least squares term in the objective function
    float mix;                  // container to store the current value of coefficients mix function in the objective function
    float smooth;               // container to store the current value of coefficients smoothness (TV) function in the objective function
    float tbsmooth;             // container to store the current value of curve smoothness function in the objective function
};

int main(int argc, char **argv);
void Save_Volume(float *volume, int Save_size,char Save_Name[128],char Save_Extension[128]);
void Read_Data(float *Data, int read_size,int Data_Type,  char *fname);
//===========================
void ST_MLEM_recon(int Nsinogram, float *sinogram, int Nvolume, float *volume,int Nx,int Ny, int Nz,float SMOOTHNESS_LAMBDA, char *fname );

void RayDrBP_SM(int Nsinogram, float *sinogram, int Nvolume, float *volume,char *fname);
void RayDrFP_SM(int Nsinogram, float *sinogram, int Nvolume, float *volume,
                char *fname);

// ------------ END HEADER, BEGIN FUNCTIONS ---------

// --- main code.
//     command line input arguments are explained in error message printed when a.out is started proj_weights
int main(int argc, char **argv) {
    
    FILE *fid;
    int i,j;
    int Nvolume, Nsinogram, Nx,Ny,Nz,Ncoefficients=0,No_Heads;
    int tbFlag = 0, tbNum, tbLen, singleProjectionSize=0, numbSegments;
    float temp;
    
    float *basis,*Tvolumes;
    float *sinogram, *volume, *coefficients=0;
    float *Mask,*StaticMask;
    struct PARAMETERS func_parameters;
    
    char sino_fname[128], image_fname[128], sysmat_fname[128],Segmented_Volume_fname[128],
    basis_fname[128], Action, *opArg;
    
    
    // to estimate running time
    clock_t t1, t2;
    float diff;
    t1 = clock();
    
    if  ((argc != 5) && (argc != 9) && (argc != 13))
    {
        for (i = 0; i < argc; i++)
            printf("%s\n", argv[i]);
        printf("#argc %d\n",argc);
        fprintf(
                stderr,
                "\nStatic Usage:\n sysmat.file sinogram.file image.file operation Volume_x Volume_y Volume_z \n");
        fprintf(
                stderr,
                "Where \"operation\" can be either R for reconstruction, or F or B for forward or backprojection. Volume_x, Volume_y, and Volume_z are the volume dimensions (Dimensions are only needed in case of TV is used in the static reconstruction)\n");
        
        fprintf(
                stderr,
                "\nDynamic Usage:\n sysmat.file sinogram.file image.file operation timebasis.file SingleProjectionSize Volume_x Volume_y Volume_z Segmented.Static.Volume NumnerOfHeads NumnerOfSegments \n");
        fprintf(
                stderr,
                "Where \"operation\" can be either TR for reconstruction, or TF or TB for forward or backprojection. Volume_x, Volume_y, and Volume_z are the static volume dimensions (Dimensions are always needed in dynamic reconstruction). NumnerOfHeads can be 1 or 2. NumnerOfSegments is the number of segments in the Segmented.Static.Volume\n");
        exit(1);
    }
    
    // Copy files' names
    strcpy(sysmat_fname, argv[1]);
    strcpy(sino_fname, argv[2]);
    strcpy(image_fname, argv[3]);
    
    // find the entered action
    opArg =strpbrk(argv[4], "rRfFbB");
    
    if (opArg == NULL)
    {
        fprintf(
                stderr,
                "Action setting should be: \"R\" or \"r\"  for reconstruction, \"F\" or \"f\" for forwardprojection,\n");
        fprintf(
                stderr,
                "or \"B\" or \"b\" for backprojection, unrecongnized action \"%s\"\n",
                argv[4]);
        exit(1);
    }
    else
    {
        switch (*opArg) {
            case 'r':
            case 'R':
                Action = 'r';
                break;
                
            case 'f':
            case 'F':
                Action = 'f';
                break;
                
            case 'b':
            case 'B':
                Action = 'b';
                break;
        }
    }
    
    // Check if the action is for dynamic reconstruction
    opArg = strpbrk(argv[4], "Tt");
    if (opArg != NULL) {
        tbFlag = 1;
    }
    
    // -- read sinogram and image size from the sysmat file.
    if ((fid = fopen(sysmat_fname, "rb")) == NULL) {
        fprintf(stderr, "Could not open sysmat file \"%s\" for reading\n",
                sysmat_fname);
        exit(1);
    }
    fread(&Nsinogram, sizeof(int), 1, fid);
    fread(&Nvolume, sizeof(int), 1, fid);
    fclose(fid);
    
    // copy the volume dimensions in case of static reconstruction with TV or dynamic reconstruction
    if(argc == 9)
    {
        Nx = atoi(argv[5]);
        Ny =atoi(argv[6]) ;
        Nz = atoi(argv[7]);
        N_ITERATIONS_MLEM = atoi(argv[8]);
    }
    else if(argc == 13)
    {
        Nx = atoi(argv[7]);
        Ny =atoi(argv[8]) ;
        Nz = atoi(argv[9]);
    }
    
    // check if the size of volume is consistent with the input dimensions
    if(Nvolume != (Nx*Ny*Nz))
    {
        fprintf(
                stderr,
                "The input volume dimensions are not consistent with the volume size in the system matrix \n");
        exit(1);
    }
    
    
    if (tbFlag) // Dynamic -- if using time basis file, than read it and define all needed parameters
    {

    }
    else // Staic -- define all needed parameters
    {
        // allocate sinogram and volume
        sinogram = (float *) malloc(size_float * Nsinogram);
        volume = (float *) malloc(size_float * Nvolume);
        
        // check if the size of volume is consistent with the input dimensions
        if(Nvolume != (Nx*Ny*Nz))
        {
            fprintf(
                    stderr,
                    "The input volume dimensions are not consistent with the volume size in the system matrix \n");
            exit(1);
        }
    }
       
    // read input data
    if (Action == 'f') // read volume
    {
        
#ifndef FLOAT_INPUT_VOLUME
        
        // 1 = 16-bit volume, 2 = 32-bit volume (float)
        if(tbFlag)
            Read_Data(coefficients,  Ncoefficients, 1, image_fname);
        else
            Read_Data(volume,  Nvolume, 1, image_fname);
#else
        if(tbFlag)
            Read_Data(coefficients,  Ncoefficients, 2, image_fname);
        else
            Read_Data(volume,  Nvolume, 2, image_fname);
#endif
    }
    else
    {  //-------------- read sinogram
#ifndef FLOAT_INPUT_SINOGRAM
        
        // 1 = 16-bit singram, 2 = 32-bit sinogram (float)
            Read_Data(sinogram,  Nsinogram, 1, sino_fname);
#else
            Read_Data(sinogram,  Nsinogram, 2, sino_fname);
#endif
    }
    //  --- Main action ---
    switch (Action) {
        case 'r':
            if (tbFlag)  // Dynamic reconstruction
            {
             }
            else // Static reconstruction
            {
#ifdef RECON_MLEM

                fprintf( stderr,"Static MLEM reconstruction volume of length %d from sinogram of length %d\n",Nvolume, Nsinogram);
                ST_MLEM_recon(Nsinogram, sinogram, Nvolume, volume, Nx, Ny, Nz,DEFAULT_SMOOTHNESS_LAMBDA,sysmat_fname);

#endif
                // Save final volume
                Save_Volume(volume,Nvolume,image_fname,"Vol");
            }
            break;
        case 'f':
            if (tbFlag)
            {

            }
            else
            {
                fprintf( stderr,"Forward project volume of length %d into sinogram of length %d\n",Nvolume, Nsinogram);
                RayDrFP_SM(Nsinogram, sinogram, Nvolume, volume, sysmat_fname);
            }
            
            // Save final sinogram
            Save_Volume(sinogram,Nsinogram,sino_fname,"sino");
            break;
        case 'b':
            if (tbFlag)
            {

            }
            else
            {
                fprintf(stderr,"Backproject volume of length %d from sinogram of length %d\n",Nvolume, Nsinogram);
                RayDrBP_SM(Nsinogram, sinogram, Nvolume, volume, sysmat_fname);
                
                // Save final volume
                Save_Volume(volume,Nvolume,image_fname,"Vol");
            }
            break;
    }
    
    //free memory
    if (tbFlag)
    {
     }
    else
    {
        free(sinogram);
        free(volume);
    }
    
    t2 = clock();
    diff = ((float)t2 - (float)t1) / (CLOCKS_PER_SEC);
    fprintf(stderr, "Program running time: %f \n",diff);
    return 1;
}
// reads volume or sinogram from HDD to container volume in memory
void Read_Data(float *Data, int read_size,int Data_Type,  char *fname)
{
    FILE *fid;
    int n;
    unsigned short tmp;

    // open the file for reading
    if ((fid = fopen(fname, "rb")) == NULL) {
        fprintf(stderr,
                "Could not open volume data file \"%s\" for reading\n",
                fname);
        exit(1);
    }
    
    // read the file
    if(Data_Type == 1) // 16-bit data
    {
        for (n = 0; n < read_size; n++) {
            fread(&tmp, 2, 1, fid); // 2 = sizeof(unsigned short) == 16 bit
            Data[n] = tmp;
        }
    }
    else if (Data_Type == 2) // 32-bit
    {
        fread(Data, size_float, read_size, fid);
    }
    else
    {
        fprintf(stderr,"Can't read data. Data type must either 16-bit or 32-bit\n");
        exit(1);
    }
    fclose(fid);
    
    return;
}

// Saves volume or singram to HDD
void Save_Volume(float *volume, int Save_size,char Save_Name[128],char Save_Extension[128])
{
    FILE *fin;
    char buf[128];
    
    // construct the save name
    strcpy(buf, Save_Name);
    strcat(buf, ".");
    strcat(buf, Save_Extension);
    
    //open file for writing
    if ((fin = fopen(buf, "wb")) == NULL) {
        fprintf(stderr, "Could not open data file %s for writing\n",buf);
        exit(1);
    }
    // write the data
    fwrite(volume, size_float, Save_size, fin);
    
    //close the file
    fclose(fin);
    
    return;
}


/////////////////////////////////
/////////////////////////////////
// recostructs a static volume with MLEM algrithm + TV smmoothing regularization
void ST_MLEM_recon(int Nsinogram, float *sinogram, int Nvolume, float *volume,int Nx,int Ny, int Nz,float SMOOTHNESS_LAMBDA, char *fname ) {
    int n, niter;
    float *tmp_volume, *tmp_sinogram, *constant_denominator;
    float denominator,no_reg_denominator;
    float reg_denominator,norm0,xm_square;
    float *mask, *grad;
    
    //these are used for dynamic lambda estimation
    float NN,T, gamma_1,gamma_0;
    NN=T=0;
    gamma_1=1.;
    gamma_0=200;
    
    // ------ create ml_em variables
    mask = (float *) malloc(size_float * Nvolume);
    grad = (float *) malloc(size_float * Nvolume);
    
    tmp_sinogram = (float *) malloc( (size_float * Nsinogram));
    constant_denominator = (float *) malloc( (size_float * Nvolume));
    tmp_volume = (float *) malloc( (size_float * Nvolume));
    
    // --- initial volume assignment: all pixels are one
    for (n = 0; n < Nvolume; n++)
    {
        volume[n] = 1.;
        grad[n]=0.0;
        mask[n]=1.0;
    }
    
    // --- compute  element-by-element inverse of efficiency matrix
    for (n = 0; n < Nsinogram; n++)
        tmp_sinogram[n] = 1.;
    
    RayDrBP_SM(Nsinogram, tmp_sinogram, Nvolume, constant_denominator, fname);
    
    // Estimate the initial noise level
    norm0=0.0;
    for (n = 0; n < Nsinogram; n++)
    {
        norm0  += pow(sinogram[n],2.0);
    }
    
    //  -------- ITERATION LOOP --------
	for (niter = 1; niter <= N_ITERATIONS_MLEM; niter++) {
        fprintf(stderr, "Iteration No [%d] of [%d] ", niter, N_ITERATIONS_MLEM);
        
        // compute the reprojection through the n-1 version of the file into tmp_sinogram
        RayDrFP_SM(Nsinogram, tmp_sinogram, Nvolume, volume, fname);
        
        // the error value (difference between sinogram and estimated sinogram)
        xm_square = 0.0;
        // divide the sinogram by the tmp_sinogram
        for (n = 0; n < Nsinogram; n++)
        {
            // Calculate the error value (difference between sinogram and estimated sinogram)
            xm_square  += pow(sinogram[n] - tmp_sinogram[n],2.0);
            
            // Divide the sinogram by the estimated sinogram
            if (sinogram[n] == 0.)
                tmp_sinogram[n] = 0.;
            else if (sinogram[n] < 0.)
                fprintf(stderr, "sinogram in MLEM smaller than zero");
            else if (tmp_sinogram[n] > EPSILON)
                tmp_sinogram[n] = sinogram[n] / tmp_sinogram[n];
            else
                tmp_sinogram[n] = sinogram[n] * ONE_OVER_EPSILON;
        }
        
        
// #ifdef STATIC_RECON_WITH_TV
        fprintf( stderr,"Error value=:%f: ",xm_square);
        
//         // Nearest neighbour constraint
//         NN = NN_3D(grad, volume, mask,Nx, Ny, Nz,100);
//         fprintf( stderr,"Smoothness value=:%f: ",NN);
        
// #ifdef DYNAMIC_LAMBDA_SELECTION
//         //Estimate lambdas value
//         gamma_1 = GAMMA(xm_square,norm0,gamma_0);
        
//         // find a new lambda value
//         if(NN>=1.0)
//             SMOOTHNESS_LAMBDA =(xm_square/Nsinogram)/(gamma_1 * NN/Nvolume);
        
//         fprintf( stderr,"GAMMA=:%f:",gamma_1);
        
// #endif
//         fprintf( stderr,"TV_LAMBDA=:%f:",SMOOTHNESS_LAMBDA);
// #endif
        
        // backproject the result into tmp_volume
        RayDrBP_SM(Nsinogram, tmp_sinogram, Nvolume, tmp_volume, fname);
        

        float xc=0;
        for (n = 0; n < Nvolume; n++) {
            xc += pow(tmp_volume[n]-constant_denominator[n],2.0); //1.0/constant_denominator[n]
        }
        //printf("convergence :\t%f",xc);


        // multiply by the constant denominator
        for (n = 0; n < Nvolume; n++) {
            
            no_reg_denominator= (constant_denominator[n]>EPSILON? constant_denominator[n]:EPSILON);
            
            reg_denominator= (constant_denominator[n] + SMOOTHNESS_LAMBDA * grad[n]);
            
            denominator = 1. / ((reg_denominator>1.)? reg_denominator:no_reg_denominator);
            volume[n] *=  denominator * tmp_volume[n];
            
            grad[n] =0.0;
        }
        
        fprintf(stderr, "\n");
    }
    // end: free memory up
    free(mask);
    free(grad);
    free(tmp_sinogram);
    free(tmp_volume);
    free(constant_denominator);
    
    return;
}

/////////////// Sinogram driven system-matrix based backprojection
void RayDrBP_SM(int Nsinogram, float *sinogram, int Nvolume, float *volume,
                char *fname) {
    #ifndef RECON_3D
    RayDrBP_SM_2D(Nsinogram, sinogram, Nvolume, volume, fname);
    return;
    #endif
    int ns, n, Nchunk;
    float *SMtemp;
    int *Itemp;
    FILE *fid;
    
    if ((fid = fopen(fname, "rb")) == NULL) {
        fprintf(stderr, "Could not open sysmat file %s\n", fname);
        exit(1);
    }
	
    // read in Nsinogram, compare to given
    fread(&n, size_int, 1, fid);
    if (n != Nsinogram) {
        fprintf(stderr, "Read in Nvsinogram %d not equal to expected %d\n", n,
                Nsinogram);
        exit(1);
    }
	
    // read in Nvolume, compare to given
    fread(&n, size_int, 1, fid);
    if (n != Nvolume) {
        fprintf(stderr, "Read in Nvolume %d not equal to expected %d\n", n,
                Nvolume);
        exit(1);
    }
    
    // set volume values to zero
    memset(volume, 0, size_float*Nvolume);
    // Declare sm chunk and volume index chunk
    Itemp = (int *) malloc(size_int * Nvolume);
    SMtemp = (float *) malloc(size_float * Nvolume);
    
    // Main loop
    for (ns = 0; ns < Nsinogram; ns++) {
        // read volume index and compare to expected
        fread(&n, size_int, 1, fid);
        if (n != ns) {
            fprintf(stderr,
                    "Read in sinogram index %d not equal to expected %d\n", n,
                    ns);
            exit(1);
        }
        // fread chunk size and indices and SM chunks
        fread(&Nchunk, size_int, 1, fid);
        if (Nchunk > Nvolume) {
            fprintf(stderr, "Sinogram chunk %d is longer than Nvolume %d\n",
                    Nchunk, Nvolume);
            exit(1);
        }
        fread(Itemp, size_int, Nchunk, fid);
        fread(SMtemp, size_float, Nchunk, fid);
        // is sinogram pixel is non-zero, do loop muptiplication
        if (fabs(sinogram[ns]) > 1.0e-14)
            for (n = 0; n < Nchunk; n++)
                volume[Itemp[n]] += sinogram[ns] * SMtemp[n];
    }
    // done, free memory, close file and return
    fclose(fid);
    free(Itemp);
    free(SMtemp);
    return;
}

/////////////// Sinogram driven system-matrix based forward projection
void RayDrFP_SM(int Nsinogram, float *sinogram, int Nvolume, float *volume,
                char *fname) {

    #ifndef RECON_3D
    RayDrFP_SM_2D(Nsinogram, sinogram, Nvolume, volume, fname);
    return;
    #endif
    int ns, n, Nchunk;
    float *SMtemp;
    int *Itemp;
    FILE *fid;
    
    if ((fid = fopen(fname, "rb")) == NULL) {
        fprintf(stderr, "Could not open sysmat file %s\n", fname);
        exit(1);
    }
    // read in Nsinogram, compare to given
    fread(&n, size_int, 1, fid);
    if (n != Nsinogram) {
        fprintf(stderr, "Read in Nvsinogram %d not equal to expected %d\n", n,
                Nsinogram);
        exit(1);
    }
    // read in Nvolume, compare to given
    fread(&n, size_int, 1, fid);
    if (n != Nvolume) {
        fprintf(stderr, "Read in Nvolume %d not equal to expected %d\n", n,
                Nvolume);
        exit(1);
    }
    // set sinogram values to zero
    // memset(sinogram, 0, size_float*Nsinogram); // don't need to do it here, done in the loop
    // Declare sm chunk and volume index chunk
    Itemp = (int *) malloc(size_int * Nvolume);
    SMtemp = (float *) malloc(size_float * Nvolume);
    
    // Main loop
    for (ns = 0; ns < Nsinogram; ns++) {
        // read volume index and compare to expected
        fread(&n, size_int, 1, fid);
        if (n != ns) {
            fprintf(stderr,
                    "Read in sinogram index %d not equal to expected %d\n", n,
                    ns);
            exit(1);
        }
        // fread chunk size and indices and SM chunks
        fread(&Nchunk, size_int, 1, fid);
        if (Nchunk > Nvolume) {
            fprintf(stderr, "ns = %d, nchunk = %d\n", ns, Nchunk);
            fprintf(
                    stderr,
                    "Forward proj: System matrix chunk length %d is longer than Nvolume %d\n",
                    Nchunk, Nvolume);
            exit(1);
        }
        fread(Itemp, size_int, Nchunk, fid);
        fread(SMtemp, size_float, Nchunk, fid);
        sinogram[ns] = 0.;
        // do loop muptiplication
        for (n = 0; n < Nchunk; n++)
            if (fabs(volume[Itemp[n]]) > 1.0e-14)
                sinogram[ns] += volume[Itemp[n]] * SMtemp[n];
    }
    // done, free memory, close file and return
    fclose(fid);
    free(Itemp);
    free(SMtemp);
    return;
}