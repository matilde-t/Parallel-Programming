#include <cstring>
#include <iostream>

#include "vv-aes.h"

/**
 * This function takes the characters stored in the 7x7 message array and
 * substitutes each character for the corresponding replacement as specified in
 * the originalCharacter and substitutedCharacter array. This corresponds to
 * step 2.1 in the VV-AES explanation.
 */

uint8_t dict[UNIQUE_CHARACTERS];

int powers[UNIQUE_CHARACTERS][BLOCK_SIZE + 1];

/*
 * This function shifts (rotates) a row in the message array by one place to the
 * left.
 * @param row The row which to shift.
 */

/*
 * This function shifts each row by the number of places it is meant to be
 * shifted according to the AES specification. Row zero is shifted by zero
 * places. Row one by one, etc. This corresponds to step 2.2 in the VV-AES
 * explanation.
 */
void shift_rows(int row, int times) {
  // Shift each row, where the row index corresponds to how many columns the
  // data is shifted.
  for (int shifts = 0; shifts < times; ++shifts) {
    auto first_char = message[row][0];
    for (unsigned int i = 0; i < BLOCK_SIZE - 1; i++) {
      message[row][i] = message[row][i + 1];
    }
    message[row][BLOCK_SIZE - 1] = first_char;
  }
}

void substitute_bytes_and_shift() {
  // For each byte in the message
  for (int row = 0; row < BLOCK_SIZE; row++) {
    for (int column = 0; column < BLOCK_SIZE; column++) {
      // Search for the byte in the original character list
      // and replace it with corresponding the element in the substituted
      // character list
      message[row][column] = dict[message[row][column]];
    }
    shift_rows(row, row);
  }
}

/*
 * This function calculates x^n for polynomial evaluation.
 */
void precalculate_powers() {
  for (int i = 0; i < UNIQUE_CHARACTERS; ++i) {
    powers[i][0] = 1;
    for (int j = 1; j < BLOCK_SIZE + 1; ++j) {
      powers[i][j] = powers[i][j - 1] * i;
    }
  }
}

/*
 * This function evaluates four different polynomials, one for each row in the
 * column. Each polynomial evaluated is of the form m'[row, column] = c[r][3]
 * m[3][column]^4 + c[r][2] m[2][column]^3 + c[r][1] m[1][column]^2 +
 * c[r][0]m[0][column]^1 where m' is the new message value, c[r] is an array of
 * polynomial coefficients for the current result row (each result row gets a
 * different polynomial), and m is the current message value.
 *
 */

void multiply_with_polynomial(int row) {
  for (int column = 0; column < BLOCK_SIZE; ++column) {
    int result = 0;
    for (int degree = 0; degree < BLOCK_SIZE; degree++) {
      result += polynomialCoefficients[row][degree] *
                powers[message[degree][column]][degree + 1];
    }
    message[row][column] = result;
  }
}

/*
 * For each column, mix the values by evaluating them as parameters of multiple
 * polynomials. This corresponds to step 2.3 in the VV-AES explanation.
 */
void mix_columns() {
  for (int column = 0; column < BLOCK_SIZE; ++column) {
    multiply_with_polynomial(column);
  }
}

/*
 * Add the current key to the message using the XOR operation.
 */
void add_key() {
  for (int row = 0; row < BLOCK_SIZE; ++row) {
    for (int column = 0; column < BLOCK_SIZE; column++) {
      // ^ == XOR
      message[row][column] = message[row][column] ^ key[row][column];
    }
  }
}

void create_dict() {
  // Create dictionary in the form substitutedCharacter =
  // dict[originalCharacter]
  for (unsigned int i = 0; i < UNIQUE_CHARACTERS; i++) {
    auto key = originalCharacter[i];
    dict[key] = substitutedCharacter[i];
  }
}

/*
 * Your main encryption routine.
 */
int main() {
  create_dict();
  precalculate_powers();
  // Receive the problem from the system.
  readInput();

  // For extra security (and because Varys wasn't able to find enough test
  // messages to keep you occupied) each message is put through VV-AES lots of
  // times. If we can't stop the adverse Maesters from decrypting our highly
  // secure encryption scheme, we can at least slow them down.
  for (int i = 0; i < ITERATIONS; i++) {
    // For each message, we use a predetermined key (e.g. the password). In our
    // case, its just pseudo random.
    set_next_key();

    // First, we add the key to the message once:
    add_key();

    // We do 9+1 rounds for 128 bit keys
    for (int round = 0; round < ROUNDS; round++) {
      // In each round, we use a different key derived from the original (refer
      // to the key schedule).
      set_next_key();

      // These are the four steps described in the slides.
      substitute_bytes_and_shift();
      mix_columns();
      add_key();
    }
    // Final round
    substitute_bytes_and_shift();
    add_key();
  }

  // Submit our solution back to the system.
  writeOutput();
  return 0;
}
