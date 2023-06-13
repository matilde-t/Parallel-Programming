#include "Utility.h"
#include <cstdlib>
#include <math.h>
#include <string.h>
#include <unistd.h>

// include library for intel intrinsics
#include <immintrin.h>

// uncomment this line to print used time
// comment this line before submission
#define PRINT_TIME 1

/*
 * Dense matrix vector multiplication
 */
void dmv(float *mat, float *in_vec, float *out_vec, size_t mat_size) {
  float sum = 0;
  for (size_t row = 0; row < mat_size; row++) {
    out_vec[row] = 0;
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

int main() {
  unsigned int seed = readInput();

  // allocate aligned memory
  float *mat = (float *)_mm_malloc(sizeof(float) * MAT_SIZE * MAT_SIZE, 32);
  float *in_vec = (float *)_mm_malloc(sizeof(float) * MAT_SIZE, 32);
  float *out_vec = (float *)_mm_malloc(sizeof(float) * MAT_SIZE, 32);

#ifdef PRINT_TIME
  TicToc total_time;
#endif
  generate_test(seed, mat, in_vec);

  for (int i = 0; i < ITER_NUM; i++) {
    dmv(mat, in_vec, out_vec, MAT_SIZE);
    out_vec[i] = 1.0 / double(i + 1);
    memcpy(in_vec, out_vec, sizeof(float) * MAT_SIZE);
  }

  outputResult(seed, out_vec);
#ifdef PRINT_TIME
  std::cerr << "time used: " << total_time.toc() << "ms.\n";
#endif
}
