#include "mex.h"
#include "math.h"
/*
 * oasis.c - implement the OASIS algorithm 
 *
 *
 * This is a MEX-file for MATLAB.
 * Copyright Gal Chechik 2009
 *
 *
 * I use a sparse representation as in pamir:
 * [num_paires, index,value, index,value, index,value... ] 
 *
 *
 */

int CountNonzerosOfFullVec(double *full, int dim) {
  int count=0;
  int i;
  for (i=0; i<dim; i++) {
    count += (full[i] != 0);
  }
  return count;
}

void FillASparseVec(double *full, double *sparse, int dim, int count) {
  sparse[0] = count;
  /* printf("\t\t\t FillASparseVec: count=%d\n", count); */

  int sparse_entry=1;
  int i;
  for (i=0; i<dim; i++) { 
    if (full[i] != 0) {
      sparse[sparse_entry] = i;
      sparse_entry++;
      sparse[sparse_entry] = full[i];
      sparse_entry++;
      /* printf("\t\t\t\t%d: %d->%f\n", sparse_entry, i, full[i]);*/
    }
  }
  return;
}

/* Translate a given full vector into a sparse one */
double* SparsifyAFullVec(double *full, int dim) {  
  int count = CountNonzerosOfFullVec(full, dim);
  int i;
 /* printf("\t\t\t dim = %d ", dim);*/
  /* 
     printf("\t\t\t SparsifyAFullVec: Full=");
     for (i=0; i<dim; i++) {
     printf("%f ", full[i]); 
     }
  */
  /*   printf(" count=%d\n", count);*/
    double *sparse = malloc(sizeof(double)*(1+2*count));
  
  FillASparseVec(full, sparse, dim, count);
  return sparse;
}

double** SparsifyImages(double *images_full, int num_images, int dim) {
  int i;
  double **images_sparse = malloc(sizeof(double*)*num_images);
  

  for (i=0; i<num_images; i++) {
    images_sparse[i] = SparsifyAFullVec(&images_full[i*dim], dim);
  }
  return images_sparse;
}




/* Sample an integer in [0,..,num-1] */
int Sample(int num) {
  int s;
  s = (int)((double)num * rand()/(1.0+RAND_MAX));
  return s;
}


/* Sample an image from a class 
 * Note that class_start is in "matlab counting" from 1..n
 * Output is in 0..n-1 indices
 */
int SampleImageForClass(int class_start, int class_size) { 
  int s = (class_start - 1) + Sample(class_size);
  return s;
}

/* dot product between two vectors */
double DotProduct(double* image1, double *image2) {
  int count1 = image1[0];
  int count2 = image2[0];
  int index1, index2;
  double norm = 0.0;
  int position1 = 0;
  int position2 = 0;
  while (position1 < count1 && position2 < count2) {
    index1  = (int)(image1[1+position1*2]+0.1);
    index2  = (int)(image2[1+position2*2]+0.1);
    if (index1==index2) { 
      norm += image1[1+position1*2+1] * image2[1+position2*2+1];
      position1++;
      position2++;
    } else{
      if (index1 < index2) {
        position1++;
      } else {
        position2++;
      }
    }
  }
  return norm;
}


double SumSquares(double* image) {
  int i;
  int count = image[0];
  double sum_of_squares = 0.0;  
  for (i=0; i<count; i++) {
    sum_of_squares += image[1+i*2+1] * image[1+i*2+1];
    /* printf("%i,%f ",i, image[1+i*2+1]); */
  }
  return sum_of_squares;
}

double ComputeSquareNormGradW(double *query, double* image_pos, double *image_neg) {
  /* For query q, positive image p and negative image n, the squared norm
   * |q * (p - n)|^2 is implemented as q^2 * (p^2 - 2 p n + n^2).
   */
  double query_norm, plus_norm, minus_norm, images_norm;
  query_norm = SumSquares(query);
  plus_norm = SumSquares(image_pos);
  minus_norm = SumSquares(image_neg);
  
  images_norm = plus_norm + minus_norm - 2.0 * DotProduct(image_pos, image_neg);

  /*  printf("\t\t ComputeSquareNormGradW: images_norm=%f, query_norm=%f\n",
      images_norm, query_norm);*/
  return query_norm * images_norm;
}

double ComputeScore(double* image1, double *W, double *image2, int dim) {
  int i,j;
  int count1 = image1[0];
  int count2 = image2[0];
  double score = 0.0;
  int index1, index2;
  double value1, value2;  

  /* 
  printf("\t\t ComputScore: ");
  printf("count(image1)=%d, ", count1);  
  printf("count(image2)=%d\n", count2);
  */
  
  for (i=1; i<1+count1*2; i+=2) {
    index1 = (int)image1[i];
    value1 = image1[i+1];
    /* printf("\t\t\t i=%d, index1=%d, value1=%f\n", i, index1, value1);*/
    
    for (j=1; j<1+count2*2; j+=2) {      
      index2 = (int)image2[j];
      value2 = image2[j+1];
      /* printf("\t\t\t\t j=%d index1=%d, value1=%f\n", j, index2, value2);*/
      score += value1 * W[index1+index2*dim] * value2;
    }
  }
  /* printf("\t\t Computed score = %f\n", score);*/
  return score;
}


void UpdateW(double *W, double tau, int dim, 
             double *image_query,
             double *image_pos,
             double *image_neg) {
  int i, t, j;
  int count_query = image_query[0];
  int count_pos = image_pos[0];
  int count_neg = image_neg[0];
  int index_query, index_pos, index_neg;
  double value_query, value_pos, value_neg;
  for (t=1; t < 1+count_query*2; t+=2) {
    index_query = (int)image_query[t];
    value_query = image_query[t+1];
  
    for (i=1; i < 1+count_pos*2; i+=2) {
      index_pos = (int)image_pos[i];
      value_pos = image_pos[i+1];    
      W[index_query + index_pos*dim] += value_query * tau * value_pos;
    }

    for (i=1; i < 1+count_neg*2; i+=2) {
      index_neg = (int)image_neg[i];
      value_neg = image_neg[i+1];    
      W[index_query + index_neg*dim] -= value_query * tau * value_neg;
    }
  }
  return;
}


/* Check arguments */
void CheckInputs(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

  mwSize mrows_W, ncols_W;             /* Size of first input */
  mwSize mrows_images, ncols_images;   /* Size of second input */
  mwSize mrows_class_labels, ncols_class_labels; /* Size of third input */
  mwSize mrows_class_start, ncols_class_start; /* Size of fourth input */
  mwSize mrows_class_size, ncols_class_size; /* Size of fifth input */
  mwSize mrows_numsteps, ncols_numsteps; /* Size of sixth input */
  mwSize mrows_aggress, ncols_aggress; /* Size of seventh input */
  mwSize mrows_rseed, ncols_rseed; /* Size of eighth input */
  mwSize mrows_call_no, ncols_call_no; /*Size of ninth input */
  
  static char msgbuf[256];
  int expected = 9;
    
  /* Check for proper number of arguments. */
  if(nrhs != expected) {
    printf("\n\n\tUsage:\n\t\t");    
    printf("W=oasis(W0, images, class_start, class_size, num_steps, aggress)");
    printf("\n\n");
    sprintf(msgbuf, "%d inputs required, %d found", expected, nrhs);
    mexErrMsgTxt(msgbuf);
  } else if(nlhs > 2) {
    mexErrMsgTxt("Too many output arguments");
  } else if (2 > nlhs) {
    mexErrMsgTxt("Too few output arguments");
  }

  /* W */  
  if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])){
    mexErrMsgTxt("Oasis: Input #1 must be a real vector.");
  }
  mrows_W = mxGetM(prhs[0]);
  ncols_W = mxGetN(prhs[0]);
  if (mrows_W != ncols_W ) {
    sprintf(msgbuf, "Oasis: Input #1 must be a square matrix, but dim=[%d %d]",
            mrows_W, ncols_W);
    mexErrMsgTxt(msgbuf);
  }
  printf("\t\t Good size W %d x %d\n", mrows_W, ncols_W);

  
  /* images */
  if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])){
    mexErrMsgTxt( "Oasis: Input #2 must be a real vector.");
  }
  mrows_images = mxGetM(prhs[1]);
  ncols_images = mxGetN(prhs[1]);
  if (mrows_images!=mrows_W) {
    sprintf(msgbuf, "Oasis: Input #2 and #1 must have same n-rows, but %d!=%d",
            ncols_images, ncols_W);
    mexErrMsgTxt(msgbuf);
  }
  printf("\t\t Good size images %d features x %d images\n", mrows_images,
         ncols_images);
  
  
  /* class_labels */
  if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])){
    mexErrMsgTxt( "Oasis: Input #2 must be a real vector.");
  }
  mrows_class_labels = mxGetM(prhs[2]);
  ncols_class_labels = mxGetN(prhs[2]);
  if (mrows_class_labels!=ncols_images) {
    sprintf(msgbuf, "Oasis: Input #3 and #2 must have same n-rows, but %d!=%d",
            mrows_class_labels,ncols_images);
    mexErrMsgTxt(msgbuf);
  }
  
  if (ncols_class_labels!=1) {
    sprintf(msgbuf, "Oasis: Input #3 is not a column vector, size = %d x %d",
            mrows_class_labels, ncols_class_labels);
    mexErrMsgTxt(msgbuf);
  }
  printf("\t\t Good size class_labels: %d x %d\n", mrows_class_labels,
         ncols_class_labels);  
  
  
  /* class_starts */
  if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])){
    mexErrMsgTxt( "Oasis: Input #4 must be a real vector.");
  }
  mrows_class_start = mxGetM(prhs[3]);
  ncols_class_start = mxGetN(prhs[3]);
  if (ncols_class_start!=1) {
    sprintf(msgbuf, "Oasis: Input #4 is not a column vector, size = %d x %d",
            mrows_class_start, ncols_class_start);
    mexErrMsgTxt(msgbuf);
  }
  printf("\t\t Good size class_starts %d x %d \n", mrows_class_start, ncols_class_start);
  
  /* class_sizes */
  if (!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4])){
    mexErrMsgTxt( "Oasis: Input #5 must be a real vector.");
  }
  mrows_class_size = mxGetM(prhs[4]);
  ncols_class_size = mxGetN(prhs[4]);
  if (ncols_class_size!=1) {
    sprintf(msgbuf, "Oasis: Input #5 is not a column vector, size = %d x %d",
            mrows_class_size, ncols_class_size);
    mexErrMsgTxt(msgbuf);
  }
  if (mrows_class_size!= mrows_class_start) {
    sprintf(msgbuf, "Oasis: Input #4,#5 should have same number of rows, %d,%d",
            mrows_class_size, ncols_class_size);
    mexErrMsgTxt(msgbuf);
  }
  printf("\t\t Good size class_sizes %d x %d \n", mrows_class_size, ncols_class_size);

  
  /* num_steps */
  if (!mxIsDouble(prhs[5]) || mxIsComplex(prhs[5])){
    mexErrMsgTxt( "Oasis: Input #6 must be real.");
  }
  mrows_numsteps = mxGetM(prhs[5]);
  ncols_numsteps = mxGetN(prhs[5]);
  if (mrows_numsteps!=1 || ncols_numsteps!=1) {
    sprintf(msgbuf, "Oasis: Input #6 (num_steps) is not a scalar, size = %dx%d",
            mrows_numsteps, ncols_numsteps);
    mexErrMsgTxt(msgbuf);
  }
  printf("\t\t Good size num_steps %d x %d \n", mrows_numsteps, ncols_numsteps);
  
  
  /* aggress */
  if (!mxIsDouble(prhs[6]) || mxIsComplex(prhs[6])){
    mexErrMsgTxt( "Oasis: Input #7 must be real.");
  }
  mrows_aggress = mxGetM(prhs[6]);
  ncols_aggress = mxGetN(prhs[6]);
  if (mrows_aggress!=1 || ncols_aggress!=1) {
    sprintf(msgbuf, "Oasis: Input #7 (aggres) is not a scalar, size = %dx%d",
            mrows_aggress, ncols_aggress);
    mexErrMsgTxt(msgbuf);
  }
  printf("\t\t Good size aggress %d x %d \n", mrows_aggress, ncols_aggress);
  
  
  /* rseed */
  if (!mxIsDouble(prhs[7]) || mxIsComplex(prhs[7])){
    mexErrMsgTxt( "Oasis: Input #8 must be real.");
  }
  mrows_rseed = mxGetM(prhs[7]);
  ncols_rseed = mxGetN(prhs[7]);
  if (mrows_rseed!=1 || ncols_rseed!=1) {
    sprintf(msgbuf, "Oasis: Input #8 (seed) is not a scalar, size = %dx%d",
            mrows_rseed, ncols_rseed);
    mexErrMsgTxt(msgbuf);
  }
  
  
   /* call_no */
  if (!mxIsDouble(prhs[8]) || mxIsComplex(prhs[8])){
    mexErrMsgTxt( "Oasis: Input #9 must be real.");
  }
  mrows_call_no = mxGetM(prhs[8]);
  ncols_call_no = mxGetN(prhs[8]);
  if (mrows_call_no!=1 || ncols_call_no!=1) {
    sprintf(msgbuf, "Oasis: Input #9 (call_no) is not a scalar, size = %dx%d",
            mrows_call_no, ncols_call_no);
    mexErrMsgTxt(msgbuf);
  }
  /*printf("\t\t Good size call_no %d x %d \n", mrows_seed, ncols_seed); */
  
  
  printf("\t\t Finished checking parameters\n");
  
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  /*
   * Asumes a call of the form:
   *  [W,loss_vec]=oasis(W0, images, class_labels, class_start, class_size, num_steps, aggress, rseed, call_no);
   *
   * Inputs:
   * 0  W0      - initial weight matrix
   * 1  images  - a matrix with images in each COLUMN!
   * 2  class_labels - a column vector, for each image what class it belongs to
   * 3  class_start - a column vector, which image is first on each class
   * 4  class_size - a column vector. Number of images per class   
   * 5  n_steps - number of traiing steps
   * 6  aggres  - the aggrersssivness parameter  
   * 7  rseed    - the seed for the random sampling
   * 8  call_no  - a flag that indicates whether to initialize the random seed 
   
   * Outputs:
   *   W       - the updated matrix 
   *   loss_vec - binary vector of loss events
   */

  double *images_vec; /* the images */
  double *class_labels; /* the labels of the images */
  double **images; /* A sparse matrix with the images */  
  double *W0, *W; /* the weight matrices (input and output) */
  double *loss_vec; /* vector of loss occurences */
  double *class_start; /* Class labels, vector of starting points (0,...n-1)*/
  double *class_size; /* Class labels, vector of class sizes */  

  
  int i, j;
  int dim; /* vocabulary size (dimension of W) */
  static char msgbuf[256];
  int image_query, image_pos, image_neg;
  int num_images;
  int i_step, num_steps;
  int rseed;
  int call_no;
  double aggres; /* aggresiveness parameter */
  double loss; /* loss per triplet */
  double norm_grad_w; /* norm of V = dloss/dW */
  double tau; /* learning rate */    
  int query_label;
  
  /* Check input dimensionality */
  CheckInputs(nlhs, plhs, nrhs,prhs);

  /* Extract inputs */
  W0 = mxGetPr(prhs[0]);  
  dim = (int)mxGetM(prhs[0]);
  images_vec = mxGetPr(prhs[1]);
  num_images = mxGetN(prhs[1]);
  class_labels = mxGetPr(prhs[2]);
  class_start = mxGetPr(prhs[3]);
  class_size = mxGetPr(prhs[4]);
  num_steps =*mxGetPr(prhs[5]);
  aggres = *mxGetPr(prhs[6]);
  rseed = *mxGetPr(prhs[7]);
  call_no = *mxGetPr(prhs[8]);
  
  
  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix(dim, dim, mxREAL);
  W = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(1, num_steps, mxREAL);
  loss_vec = mxGetPr(plhs[1]);
/*  printf("\t\t dim=%d, num_images=%d, num_steps=%d, aggres=%g, seed=%d\n",
         dim,  num_images, num_steps, aggres, rseed); */
  
  /* Copy W0 to W */
  for (i=0; i<dim*dim; i++) {
    W[i] = W0[i];
  }
  printf("\t\t Copied W0\n");
  
  /* Sparsify the images data */
  
  images = SparsifyImages(images_vec, num_images, dim);
  printf("\t\t Sparsified images\n");

  

  if (call_no == 1) {
      srandom(rseed);
  }
  
  /* Train */
  for (i_step=0; i_step<num_steps; i_step++) {

    if (i_step % 1000 == 0) {
      printf(".");
    }
    
    /*
    label = Sample(num_classes);       
    label_neg = Sample(num_classes-1);
    if (label_neg>(label-1)) {
        label_neg = label_neg+1;
    }
     */
    
    image_query = Sample(num_images);

    query_label = class_labels[image_query];
    /* query labels are in matlab notation 1..n */

    image_neg = Sample(num_images-class_size[query_label-1]);
    /* image_neg is sampled from all the images without those of the same class as the query
     * so the block of images from the query class should be yanked out
     */
    if (image_neg >= class_start[query_label-1]) {
        image_neg = image_neg+class_size[query_label-1];
	}
  
    image_pos = SampleImageForClass(class_start[query_label-1], class_size[query_label-1]);
    
    /* printf("Main: Sampled images qry=%d pos=%d neg=%d (num_images=%d)\n",
       image_query, image_pos, image_neg, num_images); */
    
    loss = 1.0
        - ComputeScore(images[image_query], W, images[image_pos], dim)
        + ComputeScore(images[image_query], W, images[image_neg], dim);
    
    /* double dd =  ComputeScore(images[image_query], W, images[image_query], dim);    
    printf("\t\t loss: %f \n", dd); */
    if (loss > 0) {
      /* printf("Main: loss = %f>0: Update W\n", loss);*/
      
      norm_grad_w = ComputeSquareNormGradW(images[image_query],
                                           images[image_neg], 
                                           images[image_pos]);
      tau = loss/norm_grad_w;
      if (tau> aggres) { 
        tau = aggres;
      }

      UpdateW(W, tau, dim,
              images[image_query],
              images[image_pos],
              images[image_neg]);
      loss_vec[i_step] = 1;
    }
    else {
        loss_vec[i_step] = 0;
    }
    
  }/* Loop over steps */

  
  /* Cleanup memory */
  for (i=0; i<num_images; i++) {
    free(images[i]);
  }
  free(images);

  return;
  
}/* End of main */
