#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define max_array_size 134217728 //50331648

int main(int argc, char ** argv);
float **matrix_my_onepiece(long nrow, long ncol);
int **imatrix_my_onepiece(long nrow, long ncol);

// transpose_system_matrix(fname_in, fname_out, transpose_direction)
//transposes ray-driven to volume-driven sm if transpose_direction = 0
// or volume driven to ray driven sm if transpose_direction = 1

        
        
int main(int argc, char ** argv)
{
  FILE *fid1, *fid2;
  int j, k,m,n, Nin, Nout, Kout,Kout1, Nchunk, Nchunk_max, No_filereads, transpose_direction;
  int base_out_index, max_out_index, no_fileread;
  int *Nchunks, *Nchunkses, **indeces, *tmp_indeces, *nchunkses;
  float **entries, *tmp_entries;
  size_t sizeof_int, sizeof_float;
  
  //
  sizeof_int = sizeof(int);
  sizeof_float = sizeof(float);
  //
  
  if(argc != 4){
    fprintf(stderr, "Usage:\n >> a.out sm.in sm.out transpose_direction (0 for ray-> vox, 1 otherwise)\n");
    exit(1);}
  //
  if(atoi(argv[3]) == 0)
    transpose_direction = 0;
  else if(atoi(argv[3]) == 1)
    transpose_direction = 1;
  else{
    fprintf(stderr, "Fourth input argument should be 0 or 1, not %s\n", argv[3]);
    exit(1);}

// - open sm files
  if( (fid1=fopen(argv[1], "rb")) == NULL){
    fprintf(stderr, "Could not open system matrix file %s for reading\n", argv[1]); exit(1);} 
  if( (fid2=fopen(argv[2], "wb")) == NULL){
    fprintf(stderr, "Could not open system matrix file %s for writing\n", argv[2]); exit(1);} 
  
// read Ninput and Noutput from the input system matrix
  fread(&k, sizeof_int, 1, fid1);
  fread(&n, sizeof_int, 1, fid1);
  fwrite(&k, sizeof_int, 1, fid2);
  fwrite(&n, sizeof_int, 1, fid2);
  
  if(transpose_direction == 0){Nin = k; Nout = n;}
  else{Nin = n; Nout = k;}
  
  // some allocations
  Nchunks = (int *) malloc(sizeof_int * Nout);
  tmp_indeces = (int *) malloc(sizeof_int * Nin);
  tmp_entries = (float *) malloc(sizeof_float * Nin);
  
  fprintf(stderr, "Nin = %d, Nout = %d\n", Nin, Nout);
  
  // read sm.in once to compute Nchunks for each voxel
  for(n=0; n<Nin; n++){
     fread(&k, sizeof_int, 1, fid1);
     if(k!=n){ fprintf(stderr, "sm index expected %d, read %d\n", n,k); exit(1);}
     // read nchunk
     fread(&k, sizeof_int, 1, fid1);
     if(k<0){ 
      fclose(fid1); fclose(fid2);
      fprintf(stderr, "negative Nchunk = %d for voxel %d\n", k, n); exit(1);}
     else if(k==0) {printf("k=%d\n", k); continue;};
     //      
     fread(tmp_indeces, sizeof_int, k, fid1);
     fread(tmp_entries , sizeof_float, k, fid1);
     for(m=0; m<k; m++)
     {
         Nchunks[tmp_indeces[m]]++;
     }   
  }
  fclose(fid1);

  
  // determine size of running variables
  Nchunk_max = Nchunks[0];
  for(n=1; n<Nout; n++) if(Nchunks[n]>Nchunk_max) Nchunk_max = Nchunks[n];
  // - - - - - - - - - - 
  if(Nchunk_max > max_array_size || Nchunk_max <= 0){
    fprintf(stderr, "Nchunk max = %d, weird value, max aray size %d\n", Nchunk_max, max_array_size);
    exit(1);}
  else
    Kout = max_array_size / Nchunk_max;
  // Kout -- max number of output columns that we can hold in memory
  No_filereads = Nout / Kout + 1;
  fprintf(stderr,  "We will have to read the system matrix %d times\n", No_filereads);
  // allocate [BIG] output arrays
  indeces = imatrix_my_onepiece(Kout, Nchunk_max);
  entries =  matrix_my_onepiece(Kout, Nchunk_max);
  nchunkses = (int *) malloc(Kout * sizeof_int);
   
  // - - - - - - -main loop - - - - - - - 
  for(no_fileread=0; no_fileread<No_filereads; no_fileread++){
    fprintf(stderr, "fileread %d out of %d, nout", no_fileread, No_filereads);
    // establish brackets of reading
    base_out_index = no_fileread * Kout;
    if(base_out_index >= Nout) break;
    //
    max_out_index = base_out_index + Kout;
    if(max_out_index > Nout)   max_out_index = Nout; 
    Kout1 = max_out_index - base_out_index;    
    fprintf(stderr, " %d to %d\n", base_out_index, max_out_index);
    //
    for(n=0; n<Kout1; n++) nchunkses[n] = 0;
    //
    if( (fid1=fopen(argv[1], "rb")) == NULL){
      fprintf(stderr, "cannot open file %s for reading\n", argv[1]);exit(1);}
    fread(&k, sizeof_int, 1, fid1); // Nsinogram, checked
    fread(&k, sizeof_int, 1, fid1); // Nsvolume, checked
    // - begin main reading cycle -
    for(n=0; n<Nin; n++){
      fread(&k, sizeof_int, 1, fid1); // Nsinogram, checked
      if(k != n) {fprintf(stderr, "expect index %d, read %d\n", n,k);exit(1);}
      // read Nchunk
      fread(&k, sizeof_int, 1, fid1); // Nchunk
      if(k == 0) continue;
      fread(&tmp_indeces[0], sizeof_int,   k, fid1); // indeces
      fread(&tmp_entries[0], sizeof_float, k, fid1); // entries
      //
      for(m=0; m<k; m++)
        if(tmp_indeces[m] >= base_out_index && tmp_indeces[m] < max_out_index){
          j = tmp_indeces[m]-base_out_index;
          indeces[j][ nchunkses[j] ] = n;
          entries[j][ nchunkses[j] ] = tmp_entries[m];
          nchunkses[j]++;
          if(nchunkses[j]>Nchunks[tmp_indeces[m]]){
            fprintf(stderr, "Exceeded expected Nchunk[%d] = %d not %d\n", tmp_indeces[m], nchunkses[j], Nchunks[tmp_indeces[m]]);exit(1);}
        }        
    } // end intermediate sm.in file read, index variable n
    fclose(fid1);
    // fwrite portions of 
    for(n=base_out_index; n<max_out_index; n++){
      j = n-base_out_index;
      k = nchunkses[j];
      //
      if(k != Nchunks[n]){
        fprintf(stderr, "column %d, nchunks = %d, nchunkses = %d\n", n, k, Nchunks[n]);exit(1);}
      fwrite(&n, sizeof_int, 1, fid2);      
      fwrite(&k, sizeof_int, 1, fid2);
      //
      if(k == 0) continue;
      fwrite(indeces[j], sizeof_int,   k, fid2);
      fwrite(entries[j], sizeof_float, k, fid2);}
  }
  fclose(fid2);
  free(Nchunks);
  free(nchunkses);
  free(tmp_indeces);
  free(tmp_entries);
  free(indeces[0]); free(indeces);
  free(entries[0]); free(entries);
  fprintf(stderr, "Done !!!\n");
  return 1;
}
        
// - - - -
 
// 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8
/* allocate a float matrix m[Nrows x Ncolumns] such that all the memory is at one piece */
float **matrix_my_onepiece(long nrow, long ncol)
{
  long i;
  float **m;
  size_t size_float;
  size_float = sizeof(float);
  
  /* allocate pointers to rows */
  m=(float **) malloc((size_t)(nrow*sizeof(float*)));
  if (!m) {fprintf(stderr, "allocation failure 1 in matrix_my_onepiece()");exit(1);}

  /* allocate rows and set pointers to them */
  m[0]=(float *) malloc((size_t)(nrow*ncol*size_float));
  if (!m[0]){fprintf(stderr, "allocation failure 2 in matrix_my_onepiece()");exit(1);}

  for(i=1;i<nrow;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

// 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8 8
/* allocate a int matrix m[Nrows x Ncolumns] such that all the memory is at one piece */
int **imatrix_my_onepiece(long nrow, long ncol)
{
  long i;
  int **m;
  size_t size_int;
  size_int = sizeof(int);
  /* allocate pointers to rows */
  m=(int **) malloc((size_t)(nrow*sizeof(int*)));
  if (!m) {fprintf(stderr, "allocation failure 1 in matrix_my_onepiece()"); exit(1);}

  /* allocate rows and set pointers to them */
  m[0]=(int *) malloc((size_t)(nrow*ncol*size_int));
  if (!m[0]) {fprintf(stderr, "allocation failure 2 in matrix_my_onepiece()");exit(1);}

  for(i=1;i<nrow;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

