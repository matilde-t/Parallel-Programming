#include "dgemm.h"
#include <cstdio>
#include <cstdlib>
#include <immintrin.h>

void dmv(const float *mat, const float *in_vec, float *out_vec, size_t mat_size,
         float alpha, float beta) {
  float sum = 0;
  __m256 alpha_vec = _mm256_set1_ps(alpha);
  for (size_t row = 0; row < mat_size; row++) {
    out_vec[row] = beta;
    __m256 partial_sum = _mm256_set1_ps(0);
    float partial_sum_array[8] = {0, 0, 0, 0, 0, 0, 0, 0};
    __m256 partial_row;
    __m256 partial_vec;
    __m256 partial_product;

    for (size_t col = 0; col < mat_size; col += 8) {
      // TODO: define a __m256 type and load 8 float values in the matrix row
      // into it
      partial_row = _mm256_load_ps(mat + row * mat_size + col);

      // TODO: define a __m256 type and load 8 float values in the vector into
      // it
      partial_vec = _mm256_load_ps(in_vec + col);

      // TODO: perform element-wise product between the above two __m256 type
      // and store it in a new __m256 type
      partial_product = _mm256_mul_ps(partial_row, partial_vec);
      partial_product = _mm256_mul_ps(partial_product, alpha_vec);

      // TODO: add partial_sum and the product result, and assign the result to
      // partial_sum
      partial_sum = _mm256_add_ps(partial_sum, partial_product);
    }

    // TODO: store the partial_sum into partial_sum_array
    _mm256_store_ps(partial_sum_array, partial_sum);

    for (int i = 0; i < 8; i++) {
      out_vec[row] += partial_sum_array[i];
    }
    sum += out_vec[row];
  }
  for (size_t row = 0; row < mat_size; row++) {
    out_vec[row] /= sum;
  }
}

// void dgemm(float alpha, const float *a, const float *b, float beta, float *c)
// {
//     for (int i = 0; i < MATRIX_SIZE; i++) { // row of c
//         for (int j = 0; j < MATRIX_SIZE; j++) { // col of c
//             c[i * MATRIX_SIZE + j] *= beta;h
//             for (int k = 0; k < MATRIX_SIZE; k++) { // i-th row of a, j-th
//             row of b
//                 c[i * MATRIX_SIZE + j] += alpha * a[i * MATRIX_SIZE + k] *
//                 b[j * MATRIX_SIZE + k];
//             }
//         }
//     }
// }

// void dgemm(float alpha, const float *a, const float *b, float beta, float *c) {
//   for (int i = 0; i < MATRIX_SIZE; i++) {
//     dmv(a, b + i * MATRIX_SIZE, c + i * MATRIX_SIZE, MATRIX_SIZE, alpha, beta);
//   }
// }

int main(int, char **) {
  float alpha, beta;

  // mem allocations
  int mem_size = MATRIX_SIZE * MATRIX_SIZE * sizeof(float);
  auto a = (float *)malloc(mem_size);
  auto b = (float *)malloc(mem_size);
  auto c = (float *)malloc(mem_size);

  // check if allocated
  if (nullptr == a || nullptr == b || nullptr == c) {
    printf("Memory allocation failed\n");
    if (nullptr != a)
      free(a);
    if (nullptr != b)
      free(b);
    if (nullptr != c)
      free(c);
    return 0;
  }

  generateProblemFromInput(alpha, a, b, beta, c);

  std::cerr << "Launching dgemm step." << std::endl;
  // matrix-multiplication
  // dgemm(alpha, a, b, beta, c);
  dmv(a, b, c, MATRIX_SIZE, alpha, beta);

  outputSolution(c);

  free(a);
  free(b);
  free(c);
  return 0;
}
