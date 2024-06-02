#include <stdio.h>

#include "ct-cbn.h"

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Function to convert char* to SEXP
SEXP char_to_sexp(const char* input) {
    // Create a character vector of length 1
    SEXP result = PROTECT(allocVector(STRSXP, 1));

    // Set the value of the character vector
    SET_STRING_ELT(result, 0, mkChar(input));

    // Unprotect the SEXP object to manage R's garbage collection
    UNPROTECT(1);

    // Return the SEXP object
    return result;
}

char* read_file(const char *filename) {
    FILE *file = fopen(filename, "rb");
    if (file == NULL) {
        perror("Error opening file");
        return NULL;
    }

    // Determine the file size
    fseek(file, 0, SEEK_END);
    long file_size = ftell(file);
    rewind(file);

    // Allocate memory for the file content (+1 for the null-terminator)
    char *buffer = (char *)malloc((file_size + 1) * sizeof(char));
    if (buffer == NULL) {
        perror("Memory allocation error");
        fclose(file);
        return NULL;
    }

    // Read the file into the buffer
    size_t read_size = fread(buffer, sizeof(char), file_size, file);
    if (read_size != file_size) {
        perror("Error reading file");
        free(buffer);
        fclose(file);
        return NULL;
    }

    // Null-terminate the string
    buffer[file_size] = '\0';

    // Close the file
    fclose(file);

    return buffer;
}

SEXP ctcbn(SEXP fs, SEXP mb, SEXP bs, SEXP rs, SEXP sr, SEXP epsilon, SEXP nd, SEXP emr)
{
  const char *filestem = CHAR(STRING_ELT(fs, 0));

  // defaults:
  double eps = REAL(epsilon)[0];                         // e
  int R = INTEGER(emr)[0];                                // # of EM runs
  double S = REAL(sr)[0];                           // sampling rate \lambda_s
  unsigned int seed = (unsigned) INTEGER(rs)[0]; // r, random seed
  int verbose = 0;
  int N_draw = INTEGER(nd)[0]; // # of samples to draw
  int error_flag = 0;
  int f_flag = 1;
  int e_flag = 0;
  int GPS = 0;
  int mode = LEARN_PARAM;
  int B = INTEGER(bs)[0];  // bootstrap samples
  int BM = INTEGER(mb)[0]; // bootstrap mode
  int c = 0;

  // no epsilon provided from R
  if (eps >= 1.0) {
    e_flag = 0;
  } else {
    e_flag = 1;
  }

  if ((outputFileObj = fopen(outputFile, "w")) == NULL) {
    fprintf(stderr, "ERROR: Could not create file, %s\n", outputFile);

    return fs;
  }

  if (B > 0)
    if ((eps < 0.0) || (N_draw > 0) || ((BM == LEARN_POSET) && (!e_flag)))
      error_flag++;

  if (error_flag)
  {
    fprintf(stderr, "Error: Bad parameter values!\n");
    // exit(1);
  }

  if (!f_flag)
  {
    fprintf(stderr, "Error: No input file specified!\n");
    // exit(1);
  }

  if (B > 0)
  {
    mode = BM;
  }
  else
  {
    if (e_flag)
      mode = LEARN_BOTH;
    else
      mode = LEARN_PARAM;
  }

  fprintf(outputFileObj, "mode = %d,   %d\n", mode, e_flag);

  srand(seed);
  RNG = gsl_rng_alloc(gsl_rng_taus); // global variable
  gsl_rng_set(RNG, seed);            // seed rng

  int i, j, k;

  model M;
  read_poset(filestem, &M);

  // precompute binary expansions
  precompute_binary(M.n + 1);

  M.lin_ext = get_int_array(M.n);             // a linear extension of the poset
  double *lambda = get_double_array(M.n + 1); // Exp rates
  lambda[0] = S;

  if (GPS) // Calculate genetic progression score (GPS)
  {
    // Read data
    int N, N_u;
    int **pat = read_patterns(filestem, &N, M.n);
    int *pat_idx = get_int_array(N);

    data *D = make_data_set(pat, N, M.n, &N_u, pat_idx);

    for (k = 0; k < N; k++)
      free(pat[k]);
    free(pat);

    read_lambda(filestem, lambda, M.n);

    // Inititalize data structures
    M.J_P = bfs_order_ideals(M.P, M.n + 1, &(M.m), M.lin_ext);
    parents(&M);
    children(&M);
    construct_sublattices(D, N_u, &M);
    compatibility(D, N_u, &M);

    double *lambda_exit = get_double_array(M.m);
    compute_lambda_exit(lambda, &M, lambda_exit);

    double *Prob = get_double_array(M.m);
    double *Exp = get_double_array(M.m);

    int c = 0;
    double gps;

    // Calculate GPS for all observations and print to std
    for (k = 0; k < N; k++)
    {
      c = pat_idx[k];
      if (D[c].is_compatible)
      {
        compute_prob(lambda, &M, D[c].J_Q, lambda_exit, Prob);
        gps = compute_GPS(lambda, &M, D[c].J_Q, lambda_exit, Prob, Exp);
        for (i = 0; i < M.n; i++)
          fprintf(outputFileObj, "%i ", D[c].g[i]);
        fprintf(outputFileObj, " %f\n", gps);
      }
      else
      {
        for (i = 0; i < M.n; i++)
          fprintf(outputFileObj, "%i ", D[c].g[i]);
        fprintf(outputFileObj, " NA\n");
      }
    }
  }
  else
  {

    if (N_draw == 0) // learn model
    {
      int N, N_u;
      int **pat = read_patterns(filestem, &N, M.n);
      int *pat_idx = get_int_array(N);
      data *D = make_data_set(pat, N, M.n, &N_u, pat_idx);
      for (k = 0; k < N; k++)
        free(pat[k]);
      free(pat);

      fprintf(outputFileObj, "Poset\tEps\tAlpha\tLoglike\tlambda_s");
      for (i = 1; i <= M.n; i++)
        fprintf(outputFileObj, "\tlambda_%d", i);
      fprintf(outputFileObj, "\n");

      if (eps >= 0.0) // fixed epsilon
      {
        int b;
        eps = MIN(eps, 1.0);

        // single run:
        select_poset(0, eps, &M, lambda, D, N_u, R, mode, 1);
        if (e_flag)
          write_poset(0, filestem, M.P, M.n, -1);
        write_lambda(filestem, lambda, M.n);

        // bootstrap runs:
        double *p_orig = get_double_array(N_u); // frequencies of original data
        for (k = 0; k < N_u; k++)
          p_orig[k] = (double)D[k].count;

        int **bootstrap_count = get_int_matrix(M.n + 1, M.n + 1);       // all relations
        int **bootstrap_cover_count = get_int_matrix(M.n + 1, M.n + 1); // cover relations only
        for (i = 0; i <= M.n; i++)
          for (j = 0; j <= M.n; j++)
            bootstrap_count[i][j] = bootstrap_cover_count[i][j] = 0;

        int **T = get_int_matrix(M.n + 1, M.n + 1);
        for (b = 1; b <= B; b++)
        {
          resample(D, p_orig, N_u);
          select_poset(b, eps, &M, lambda, D, N_u, R, mode, 1); // e_flag

          int_matrix_sum(bootstrap_cover_count, M.P, bootstrap_cover_count, M.n + 1);
          transitive_closure(M.P, T, M.n + 1);
          int_matrix_sum(bootstrap_count, T, bootstrap_count, M.n + 1);

          if (e_flag)
            write_poset(0, filestem, M.P, M.n, b);
        }

        if ((B > 0) && (BM == 0))
        {
          fprintf(outputFileObj, "\nbootstrap counts: matrix entry (i,j) counts edge i-->j\n");
          fprintf(outputFileObj, "\nall poset relations =\n");
          print_int_matrix(bootstrap_count, M.n + 1, M.n + 1);
          fprintf(outputFileObj, "\ncover relations =\n");
          print_int_matrix(bootstrap_cover_count, M.n + 1, M.n + 1);
        }

        for (i = 0; i <= M.n; i++)
        {
          free(T[i]);
          free(bootstrap_count[i]);
        }
        free(T);
        free(bootstrap_count);
        free(p_orig);
      }
      else // sample epsilon
      {
        double epsilon;
        int k;
        for (k = 0; k < N; k++)
        {
          epsilon = (double)k / (2.0 * (double)N);
          select_poset(k, epsilon, &M, lambda, D, N_u, R, LEARN_BOTH, 1);
          if (e_flag)
            write_poset(k, filestem, M.P, M.n, -1);
        }
      }
      free_data(D, N_u, M.n);
      free_poset(&M);
    }
    else // sample from given model:
    {
      // construct model:
      M.lin_ext = get_int_array(M.n); // a linear extension of the poset
      M.J_P = bfs_order_ideals(M.P, M.n + 1, &(M.m), M.lin_ext);
      parents(&M);
      children(&M);
      if (verbose)
        print_model(&M);

      int **pat_draw = get_int_matrix(N_draw, M.n + 1);
      double **t_draw = get_double_matrix(N_draw, M.n + 1);

      read_lambda(filestem, lambda, M.n);

      draw_samples(M.P, lambda, M.lin_ext, M.n, pat_draw, t_draw, N_draw);
      write_patterns(filestem, pat_draw, N_draw, M.n);
      write_times(filestem, t_draw, N_draw, M.n);

      for (k = 0; k < N_draw; k++)
      {
        free(pat_draw[k]);
        free(t_draw[k]);
      }
      free(pat_draw);
      free(t_draw);

      free_model(&M);
    }
  }

  free(lambda);
  fclose(outputFileObj);

  return char_to_sexp(read_file(outputFile));
}

R_CallMethodDef callMethods[] =
    {
        // C functions extended by using .Call interface
        {"ctcbn", (DL_FUNC)&ctcbn, 2},
        {NULL, NULL, 0}};
