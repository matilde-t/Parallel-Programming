#include "dgemm.h"
#include <cstdio>
#include <cstdlib>
#include <immintrin.h>
#include <iostream>

void dgemm(float alpha, const float *a, const float *b, float beta, float *c) {
  __m256 alpha_vec = _mm256_set1_ps(alpha);
  for (int row = 0; row < MATRIX_SIZE; row++) { // row of c
    __m256 partial_mul;
    for (int col = 0; col < MATRIX_SIZE; col++) { // column of c
      __m256 partial_a;
      __m256 partial_b;
      float partial_sum_array[8] = {0, 0, 0, 0, 0, 0, 0, 0};

      __m256 partial_sum = _mm256_set1_ps(0);
      int align = 0;
      for (int j = 0; j < MATRIX_SIZE - 8; j += 8) {
        partial_a = _mm256_loadu_ps(a + row * MATRIX_SIZE + j);
        partial_b = _mm256_loadu_ps(b + col * MATRIX_SIZE + j);
        partial_mul = _mm256_mul_ps(partial_a, partial_b);
        partial_mul = _mm256_mul_ps(partial_mul, alpha_vec);
        partial_sum = _mm256_add_ps(partial_sum, partial_mul);
        align = j + 8;
      }
      _mm256_store_ps(partial_sum_array, partial_sum);
      float sum = 0;
      for (int k = 0; k < MATRIX_SIZE % 8; k++) {
        sum += alpha * a[row * MATRIX_SIZE + align + k] *
               b[col * MATRIX_SIZE + align + k];
      }
      for (int k = 0; k < 8; k++) {
        sum += partial_sum_array[k];
      }

      c[row * MATRIX_SIZE + col] *= beta;
      c[row * MATRIX_SIZE + col] += sum;
    }
  }
}

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
  dgemm(alpha, a, b, beta, c);

  outputSolution(c);

  free(a);
  free(b);
  free(c);
  return 0;
}
