#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "ntt.h"
#include "reduce.h"
#include "cbd.h"
#include "symmetric.h"

#if defined(KYBER_N_ENCODING_DIMENSIONS) && (KYBER_N_ENCODING_DIMENSIONS == 2)

///////////// GENERATED WITH ../../scripts/generate_decoding_lines.py /////////////////////////////
#if KYBER_K == 2

#define CODE_ALPHA (KYBER_Q/2)
#define CODE_BETA 422

static int32_t L[5][3] = {
    {-211, 832, 1438953},
    {844, 3330, -2950309},
    {-832, 211, 736745},
    {-1, 1, -1243},
    {3330, 844, 139789},
};

#define ABOVE_L1(x, y) lower_than_mask(L[0][0]*x + L[0][2], L[0][1]*y)
#define ABOVE_L2(x, y) lower_than_mask(L[1][0]*x + L[1][2], L[1][1]*y)
#define ABOVE_L3(x, y) lower_than_mask(L[2][0]*x + L[2][2], L[2][1]*y)
#define ABOVE_L4(x, y) lower_than_mask(L[3][0]*x + L[3][2], L[3][1]*y)
#define ABOVE_L5(x, y) lower_than_mask(L[4][0]*x + L[4][2], L[4][1]*y)

#elif KYBER_K == 3

#define CODE_ALPHA (KYBER_Q/2)
#define CODE_BETA 422

static int32_t L[5][3] = {
    {-211, 832, 1438953},
    {844, 3330, -2950309},
    {-832, 211, 736745},
    {-1, 1, -1243},
    {3330, 844, 139789},
};

#define ABOVE_L1(x, y) lower_than_mask(L[0][0]*x + L[0][2], L[0][1]*y)
#define ABOVE_L2(x, y) lower_than_mask(L[1][0]*x + L[1][2], L[1][1]*y)
#define ABOVE_L3(x, y) lower_than_mask(L[2][0]*x + L[2][2], L[2][1]*y)
#define ABOVE_L4(x, y) lower_than_mask(L[3][0]*x + L[3][2], L[3][1]*y)
#define ABOVE_L5(x, y) lower_than_mask(L[4][0]*x + L[4][2], L[4][1]*y)

#elif KYBER_K == 4

#define CODE_ALPHA (KYBER_Q/2)
#define CODE_BETA 436

static int32_t L[5][3] = {
    {-109, 416, 732626},
    {872, 3330, -2962321},
    {-416, 109, 369874},
    {-1, 1, -1229},
    {3330, 872, 58561},
};

#define ABOVE_L1(x, y) lower_than_mask(L[0][0]*x + L[0][2], L[0][1]*y)
#define ABOVE_L2(x, y) lower_than_mask(L[1][0]*x + L[1][2], L[1][1]*y)
#define ABOVE_L3(x, y) lower_than_mask(L[2][0]*x + L[2][2], L[2][1]*y)
#define ABOVE_L4(x, y) lower_than_mask(L[3][0]*x + L[3][2], L[3][1]*y)
#define ABOVE_L5(x, y) lower_than_mask(L[4][0]*x + L[4][2], L[4][1]*y)

#endif
///////////////////////////////////////////////////////////////////////////////////////////////////

#endif

/*************************************************
* Name:        poly_compress
*
* Description: Compression and subsequent serialization of a polynomial
*
* Arguments:   - uint8_t *r: pointer to output byte array
*                            (of length KYBER_POLYCOMPRESSEDBYTES)
*              - const poly *a: pointer to input polynomial
**************************************************/
void poly_compress(uint8_t r[KYBER_POLYCOMPRESSEDBYTES], const poly *a)
{
  unsigned int i,j;
  int16_t u;
  uint32_t d0;
  uint8_t t[8];

#if (KYBER_POLYCOMPRESSEDBYTES == 128)
  for(i=0;i<KYBER_N/8;i++) {
    for(j=0;j<8;j++) {
      // map to positive standard representatives
      u  = a->coeffs[8*i+j];
      u += (u >> 15) & KYBER_Q;
/*    t[j] = ((((uint16_t)u << 4) + KYBER_Q/2)/KYBER_Q) & 15; */
      d0 = u << 4;
      d0 += 1665;
      d0 *= 80635;
      d0 >>= 28;
      t[j] = d0 & 0xf;
    }

    r[0] = t[0] | (t[1] << 4);
    r[1] = t[2] | (t[3] << 4);
    r[2] = t[4] | (t[5] << 4);
    r[3] = t[6] | (t[7] << 4);
    r += 4;
  }
#elif (KYBER_POLYCOMPRESSEDBYTES == 160)
  for(i=0;i<KYBER_N/8;i++) {
    for(j=0;j<8;j++) {
      // map to positive standard representatives
      u  = a->coeffs[8*i+j];
      u += (u >> 15) & KYBER_Q;
/*      t[j] = ((((uint32_t)u << 5) + KYBER_Q/2)/KYBER_Q) & 31; */
      d0 = u << 5;
      d0 += 1664;
      d0 *= 40318;
      d0 >>= 27;
      t[j] = d0 & 0x1f;
    }

    r[0] = (t[0] >> 0) | (t[1] << 5);
    r[1] = (t[1] >> 3) | (t[2] << 2) | (t[3] << 7);
    r[2] = (t[3] >> 1) | (t[4] << 4);
    r[3] = (t[4] >> 4) | (t[5] << 1) | (t[6] << 6);
    r[4] = (t[6] >> 2) | (t[7] << 3);
    r += 5;
  }
#else
#error "KYBER_POLYCOMPRESSEDBYTES needs to be in {128, 160}"
#endif
}

/*************************************************
* Name:        poly_decompress
*
* Description: De-serialization and subsequent decompression of a polynomial;
*              approximate inverse of poly_compress
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: pointer to input byte array
*                                  (of length KYBER_POLYCOMPRESSEDBYTES bytes)
**************************************************/
void poly_decompress(poly *r, const uint8_t a[KYBER_POLYCOMPRESSEDBYTES])
{
  unsigned int i;

#if (KYBER_POLYCOMPRESSEDBYTES == 128)
  for(i=0;i<KYBER_N/2;i++) {
    r->coeffs[2*i+0] = (((uint16_t)(a[0] & 15)*KYBER_Q) + 8) >> 4;
    r->coeffs[2*i+1] = (((uint16_t)(a[0] >> 4)*KYBER_Q) + 8) >> 4;
    a += 1;
  }
#elif (KYBER_POLYCOMPRESSEDBYTES == 160)
  unsigned int j;
  uint8_t t[8];
  for(i=0;i<KYBER_N/8;i++) {
    t[0] = (a[0] >> 0);
    t[1] = (a[0] >> 5) | (a[1] << 3);
    t[2] = (a[1] >> 2);
    t[3] = (a[1] >> 7) | (a[2] << 1);
    t[4] = (a[2] >> 4) | (a[3] << 4);
    t[5] = (a[3] >> 1);
    t[6] = (a[3] >> 6) | (a[4] << 2);
    t[7] = (a[4] >> 3);
    a += 5;

    for(j=0;j<8;j++)
      r->coeffs[8*i+j] = ((uint32_t)(t[j] & 31)*KYBER_Q + 16) >> 5;
  }
#else
#error "KYBER_POLYCOMPRESSEDBYTES needs to be in {128, 160}"
#endif
}

/*************************************************
* Name:        poly_tobytes
*
* Description: Serialization of a polynomial
*
* Arguments:   - uint8_t *r: pointer to output byte array
*                            (needs space for KYBER_POLYBYTES bytes)
*              - const poly *a: pointer to input polynomial
**************************************************/
void poly_tobytes(uint8_t r[KYBER_POLYBYTES], const poly *a)
{
  unsigned int i;
  uint16_t t0, t1;

  for(i=0;i<KYBER_N/2;i++) {
    // map to positive standard representatives
    t0  = a->coeffs[2*i];
    t0 += ((int16_t)t0 >> 15) & KYBER_Q;
    t1 = a->coeffs[2*i+1];
    t1 += ((int16_t)t1 >> 15) & KYBER_Q;
    r[3*i+0] = (t0 >> 0);
    r[3*i+1] = (t0 >> 8) | (t1 << 4);
    r[3*i+2] = (t1 >> 4);
  }
}

/*************************************************
* Name:        poly_frombytes
*
* Description: De-serialization of a polynomial;
*              inverse of poly_tobytes
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: pointer to input byte array
*                                  (of KYBER_POLYBYTES bytes)
**************************************************/
void poly_frombytes(poly *r, const uint8_t a[KYBER_POLYBYTES])
{
  unsigned int i;
  for(i=0;i<KYBER_N/2;i++) {
    r->coeffs[2*i]   = ((a[3*i+0] >> 0) | ((uint16_t)a[3*i+1] << 8)) & 0xFFF;
    r->coeffs[2*i+1] = ((a[3*i+1] >> 4) | ((uint16_t)a[3*i+2] << 4)) & 0xFFF;
  }
}


#if !defined(KYBER_N_ENCODING_DIMENSIONS) || (KYBER_N_ENCODING_DIMENSIONS == 1)

/*************************************************
* Name:        poly_frommsg
*
* Description: Convert 32-byte message to polynomial
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *msg: pointer to input message
**************************************************/
void poly_frommsg(poly *r, const uint8_t msg[KYBER_INDCPA_MSGBYTES])
{
  unsigned int i,j;
  int16_t mask;

#if (KYBER_INDCPA_MSGBYTES != KYBER_N/8)
#error "KYBER_INDCPA_MSGBYTES must be equal to KYBER_N/8 bytes!"
#endif

  for(i=0;i<KYBER_N/8;i++) {
    for(j=0;j<8;j++) {
      mask = -(int16_t)((msg[i] >> j)&1);
      r->coeffs[8*i+j] = mask & ((KYBER_Q+1)/2);
    }
  }
}

/*************************************************
* Name:        poly_tomsg
*
* Description: Convert polynomial to 32-byte message
*
* Arguments:   - uint8_t *msg: pointer to output message
*              - const poly *a: pointer to input polynomial
**************************************************/
void poly_tomsg(uint8_t msg[KYBER_INDCPA_MSGBYTES], const poly *a)
{
  unsigned int i,j;
  uint32_t t;

  for(i=0;i<KYBER_N/8;i++) {
    msg[i] = 0;
    for(j=0;j<8;j++) {
      t  = a->coeffs[8*i+j];
      // t += ((int16_t)t >> 15) & KYBER_Q;
      // t  = (((t << 1) + KYBER_Q/2)/KYBER_Q) & 1;
      t <<= 1;
      t += 1665;
      t *= 80635;
      t >>= 28;
      t &= 1;
      msg[i] |= t << j;
    }
  }
}

#elif defined(KYBER_N_ENCODING_DIMENSIONS) && (KYBER_N_ENCODING_DIMENSIONS == 2)

// Returns 0xffffffff if (v1 < v2) and 0x00000000 otherwise
static __inline__ uint32_t lower_than_mask(const uint32_t v1, const uint32_t v2) {
  return -((v1 - v2) >> 31);
}

void poly_frommsg(poly *r, const uint8_t msg[KYBER_INDCPA_MSGBYTES])
{
  unsigned int i,j;

#if (KYBER_INDCPA_MSGBYTES != KYBER_N/8)
#error "KYBER_INDCPA_MSGBYTES must be equal to KYBER_N/8 bytes!"
#endif

  size_t base = 0;
  for(i=0;i<KYBER_N/8;i++) {
    for(j=0;j<8;j+=2) {
      uint16_t bit0_mask = -((msg[i] >> j) & 1);
      uint16_t bit1_mask = -((msg[i] >> (j + 1)) & 1);
      r->coeffs[base] = (CODE_ALPHA & bit1_mask) + (CODE_BETA & bit0_mask);
      r->coeffs[base + KYBER_N/2] = (CODE_BETA & bit1_mask) + (CODE_ALPHA & bit0_mask);

      base++;
    }
  }
}

static __inline__ int decode_msg_pair(int32_t x, int32_t y) {
  uint32_t reflect_mask = lower_than_mask(x, y);

  int32_t x_prime = (x & ~reflect_mask) | (y & reflect_mask);
  int32_t y_prime = (y & ~reflect_mask) | (x & reflect_mask);

  uint8_t above_l1 = ABOVE_L1(x_prime, y_prime);
  uint8_t above_l2 = ABOVE_L2(x_prime, y_prime);
  uint8_t above_l3 = ABOVE_L3(x_prime, y_prime);
  uint8_t above_l4 = ABOVE_L4(x_prime, y_prime);
  uint8_t above_l5 = ABOVE_L5(x_prime, y_prime);

  // Not needed, but conceptually: uint8_t c00 = (~above_l3 & above_l2 & above_l4);
  uint8_t c01 = ~above_l2 & ~above_l5 & ~above_l3;
  uint8_t c10 = above_l2 & above_l3 & ~above_l1;
  uint8_t c11 = above_l1 | (above_l3 & ~above_l2) | (above_l5 & ~above_l4);

  c01 &= (1 ^ reflect_mask);
  c10 &= (2 ^ reflect_mask);

  uint8_t bits = c01 | c10 | c11;

  return bits & 3;
}

/*************************************************
* Name:        poly_tomsg
*
* Description: Convert polynomial to 32-byte message
*
* Arguments:   - uint8_t *msg: pointer to output message
*              - const poly *a: pointer to input polynomial
**************************************************/
void poly_tomsg(uint8_t msg[KYBER_INDCPA_MSGBYTES], const poly *a)
{
  unsigned int i,j;
  uint32_t base = 0;
  for(i=0;i<KYBER_N/8;i++) {
    msg[i] = 0;
    for(j=0;j<8;j+=2) {
      int x = a->coeffs[base];
      int y = a->coeffs[base + KYBER_N/2];
      base++;

      msg[i] |= (decode_msg_pair(x, y) << j);
    }
  }
}

#elif defined(KYBER_N_ENCODING_DIMENSIONS) && (KYBER_N_ENCODING_DIMENSIONS == 4)

// FishCode of dimension 4 and generator matrix:
//     [   0    0 1664  752]
//     [ 752    0    0 1664]
//     [1664  752    0    0]
//     [   0 1664  752    0]

#define CODE_ALPHA 1664
#define CODE_BETA  752

static int16_t CODE_VALS[4] = {0, CODE_BETA, CODE_ALPHA, CODE_BETA + CODE_ALPHA - KYBER_Q};

static int8_t CODEWORDS[16][4] = {
   {0, 0, 0, 0},
   {1, 2, 0, 0},
   {2, 0, 0, 1},
   {3, 2, 0, 1},
   {0, 0, 1, 2},
   {1, 2, 1, 2},
   {2, 0, 1, 3},
   {3, 2, 1, 3},
   {0, 1, 2, 0},
   {1, 3, 2, 0},
   {2, 1, 2, 1},
   {3, 3, 2, 1},
   {0, 1, 3, 2},
   {1, 3, 3, 2},
   {2, 1, 3, 3},
   {3, 3, 3, 3},
};

// Returns 0xffffffff if (v1 < v2) and 0x00000000 otherwise
static __inline__ uint32_t lower_than_mask(const uint32_t v1, const uint32_t v2) {
  return -((v1 - v2) >> 31);
}

void poly_frommsg(poly *r, const uint8_t msg[KYBER_INDCPA_MSGBYTES])
{
  unsigned int i,j;

#if (KYBER_INDCPA_MSGBYTES != KYBER_N/8)
#error "KYBER_INDCPA_MSGBYTES must be equal to KYBER_N/8 bytes!"
#endif

  size_t base = 0;
  for(i=0;i<KYBER_N/8;i++) {
    for(j=0;j<8;j+=4) {
      uint16_t b3 = -((msg[i] >> (j + 0)) & 1);
      uint16_t b2 = -((msg[i] >> (j + 1)) & 1);
      uint16_t b1 = -((msg[i] >> (j + 2)) & 1);
      uint16_t b0 = -((msg[i] >> (j + 3)) & 1);

      r->coeffs[base]               = (CODE_ALPHA & b2) + (CODE_BETA & b3);
      r->coeffs[base + KYBER_N/4]   = (CODE_ALPHA & b3) + (CODE_BETA & b0);
      r->coeffs[base + 2*KYBER_N/4] = (CODE_ALPHA & b0) + (CODE_BETA & b1);
      r->coeffs[base + 3*KYBER_N/4] = (CODE_ALPHA & b1) + (CODE_BETA & b2);

      base++;
    }
  }
}

// Returns `abs(centered_mod(value, KYBER_Q)` assuming `-KYBER_Q <= value <= KYBER_Q`
static __inline__ int32_t abs_center_mod_of_2q_centered_value(int32_t value) {
  // value = abs(value):
  uint32_t mask_sign = value >> 31;
  value ^= mask_sign;
  value += mask_sign & 1;
  value -= KYBER_Q & lower_than_mask(KYBER_Q/2, value);
  return value;
}

static __inline__ uint32_t get_distance_sqr_to_codeword(uint16_t idx, uint32_t distsqr_matrix[4][4]) {
  int32_t a0 = distsqr_matrix[0][CODEWORDS[idx][0]];
  int32_t a1 = distsqr_matrix[1][CODEWORDS[idx][1]];
  int32_t a2 = distsqr_matrix[2][CODEWORDS[idx][2]];
  int32_t a3 = distsqr_matrix[3][CODEWORDS[idx][3]];

  // Returns `distance_sqr | codeword_index`
  return (a0 + a1 + a2 + a3) << 4 | idx;
}

// Computes min(v1, v2) in constant time
static __inline__ uint32_t secure_min(int32_t v1, int32_t v2) {
  uint32_t mask_min_v1 = lower_than_mask(v1, v2);
  return (mask_min_v1 & v1) | (~mask_min_v1 & v2);
}

static __inline__ int decode_msg_4d(int16_t target[4]) {
  // This function assumes that:
  //    * -q/2 < target[i] <= q/2
  //    * -q/2 < CODEWORD[idx][i] <= q/2

  // Build matrix with square distances to target coordinates
  uint32_t distsqr_matrix[4][4] = {0};
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      distsqr_matrix[i][j] = abs_center_mod_of_2q_centered_value(target[i] - CODE_VALS[j]);
      distsqr_matrix[i][j] *= distsqr_matrix[i][j];
    }
  }

  uint32_t min_dist_codeword = get_distance_sqr_to_codeword(0, distsqr_matrix);
  for (size_t i = 1; i < 16; i++) {
    min_dist_codeword = secure_min(get_distance_sqr_to_codeword(i, distsqr_matrix), min_dist_codeword);
  }

  // Retrieves `codeword_index` from `distance_sqr | codeword_index`
  return min_dist_codeword & 0xF;
}

/*************************************************
* Name:        poly_tomsg
*
* Description: Convert polynomial to 32-byte message
*
* Arguments:   - uint8_t *msg: pointer to output message
*              - const poly *a: pointer to input polynomial
**************************************************/
void poly_tomsg(uint8_t msg[KYBER_INDCPA_MSGBYTES], const poly *a)
{
  unsigned int i;
  uint32_t base = 0;
  for(i=0;i<KYBER_N/8;i++) {
    msg[i] = 0;

    int16_t target[4] = {0};
    target[0] = a->coeffs[base];
    target[1] = a->coeffs[base + 1*KYBER_N/4];
    target[2] = a->coeffs[base + 2*KYBER_N/4];
    target[3] = a->coeffs[base + 3*KYBER_N/4];
    msg[i] |= (decode_msg_4d(target) << 0);

    target[0] = a->coeffs[base + 1];
    target[1] = a->coeffs[base + 1 + 1*KYBER_N/4];
    target[2] = a->coeffs[base + 1 + 2*KYBER_N/4];
    target[3] = a->coeffs[base + 1 + 3*KYBER_N/4];
    msg[i] |= (decode_msg_4d(target) << 4);

    base += 2;
  }
}

#elif defined(KYBER_N_ENCODING_DIMENSIONS) && (KYBER_N_ENCODING_DIMENSIONS == 8)

#define CODE_ALPHA 1664
#define CODE_BETA  752

static int16_t CODE_VALS[4] = {0, CODE_BETA, CODE_ALPHA, CODE_BETA + CODE_ALPHA - KYBER_Q};

static int8_t CODEWORDS[256][8] = {
   {0, 0, 0, 0, 0, 0, 0, 0},
   {1, 2, 0, 0, 0, 0, 0, 0},
   {2, 0, 0, 0, 0, 0, 0, 1},
   {3, 2, 0, 0, 0, 0, 0, 1},
   {0, 0, 0, 0, 0, 0, 1, 2},
   {1, 2, 0, 0, 0, 0, 1, 2},
   {2, 0, 0, 0, 0, 0, 1, 3},
   {3, 2, 0, 0, 0, 0, 1, 3},
   {0, 0, 0, 0, 0, 1, 2, 0},
   {1, 2, 0, 0, 0, 1, 2, 0},
   {2, 0, 0, 0, 0, 1, 2, 1},
   {3, 2, 0, 0, 0, 1, 2, 1},
   {0, 0, 0, 0, 0, 1, 3, 2},
   {1, 2, 0, 0, 0, 1, 3, 2},
   {2, 0, 0, 0, 0, 1, 3, 3},
   {3, 2, 0, 0, 0, 1, 3, 3},
   {0, 0, 0, 0, 1, 2, 0, 0},
   {1, 2, 0, 0, 1, 2, 0, 0},
   {2, 0, 0, 0, 1, 2, 0, 1},
   {3, 2, 0, 0, 1, 2, 0, 1},
   {0, 0, 0, 0, 1, 2, 1, 2},
   {1, 2, 0, 0, 1, 2, 1, 2},
   {2, 0, 0, 0, 1, 2, 1, 3},
   {3, 2, 0, 0, 1, 2, 1, 3},
   {0, 0, 0, 0, 1, 3, 2, 0},
   {1, 2, 0, 0, 1, 3, 2, 0},
   {2, 0, 0, 0, 1, 3, 2, 1},
   {3, 2, 0, 0, 1, 3, 2, 1},
   {0, 0, 0, 0, 1, 3, 3, 2},
   {1, 2, 0, 0, 1, 3, 3, 2},
   {2, 0, 0, 0, 1, 3, 3, 3},
   {3, 2, 0, 0, 1, 3, 3, 3},
   {0, 0, 0, 1, 2, 0, 0, 0},
   {1, 2, 0, 1, 2, 0, 0, 0},
   {2, 0, 0, 1, 2, 0, 0, 1},
   {3, 2, 0, 1, 2, 0, 0, 1},
   {0, 0, 0, 1, 2, 0, 1, 2},
   {1, 2, 0, 1, 2, 0, 1, 2},
   {2, 0, 0, 1, 2, 0, 1, 3},
   {3, 2, 0, 1, 2, 0, 1, 3},
   {0, 0, 0, 1, 2, 1, 2, 0},
   {1, 2, 0, 1, 2, 1, 2, 0},
   {2, 0, 0, 1, 2, 1, 2, 1},
   {3, 2, 0, 1, 2, 1, 2, 1},
   {0, 0, 0, 1, 2, 1, 3, 2},
   {1, 2, 0, 1, 2, 1, 3, 2},
   {2, 0, 0, 1, 2, 1, 3, 3},
   {3, 2, 0, 1, 2, 1, 3, 3},
   {0, 0, 0, 1, 3, 2, 0, 0},
   {1, 2, 0, 1, 3, 2, 0, 0},
   {2, 0, 0, 1, 3, 2, 0, 1},
   {3, 2, 0, 1, 3, 2, 0, 1},
   {0, 0, 0, 1, 3, 2, 1, 2},
   {1, 2, 0, 1, 3, 2, 1, 2},
   {2, 0, 0, 1, 3, 2, 1, 3},
   {3, 2, 0, 1, 3, 2, 1, 3},
   {0, 0, 0, 1, 3, 3, 2, 0},
   {1, 2, 0, 1, 3, 3, 2, 0},
   {2, 0, 0, 1, 3, 3, 2, 1},
   {3, 2, 0, 1, 3, 3, 2, 1},
   {0, 0, 0, 1, 3, 3, 3, 2},
   {1, 2, 0, 1, 3, 3, 3, 2},
   {2, 0, 0, 1, 3, 3, 3, 3},
   {3, 2, 0, 1, 3, 3, 3, 3},
   {0, 0, 1, 2, 0, 0, 0, 0},
   {1, 2, 1, 2, 0, 0, 0, 0},
   {2, 0, 1, 2, 0, 0, 0, 1},
   {3, 2, 1, 2, 0, 0, 0, 1},
   {0, 0, 1, 2, 0, 0, 1, 2},
   {1, 2, 1, 2, 0, 0, 1, 2},
   {2, 0, 1, 2, 0, 0, 1, 3},
   {3, 2, 1, 2, 0, 0, 1, 3},
   {0, 0, 1, 2, 0, 1, 2, 0},
   {1, 2, 1, 2, 0, 1, 2, 0},
   {2, 0, 1, 2, 0, 1, 2, 1},
   {3, 2, 1, 2, 0, 1, 2, 1},
   {0, 0, 1, 2, 0, 1, 3, 2},
   {1, 2, 1, 2, 0, 1, 3, 2},
   {2, 0, 1, 2, 0, 1, 3, 3},
   {3, 2, 1, 2, 0, 1, 3, 3},
   {0, 0, 1, 2, 1, 2, 0, 0},
   {1, 2, 1, 2, 1, 2, 0, 0},
   {2, 0, 1, 2, 1, 2, 0, 1},
   {3, 2, 1, 2, 1, 2, 0, 1},
   {0, 0, 1, 2, 1, 2, 1, 2},
   {1, 2, 1, 2, 1, 2, 1, 2},
   {2, 0, 1, 2, 1, 2, 1, 3},
   {3, 2, 1, 2, 1, 2, 1, 3},
   {0, 0, 1, 2, 1, 3, 2, 0},
   {1, 2, 1, 2, 1, 3, 2, 0},
   {2, 0, 1, 2, 1, 3, 2, 1},
   {3, 2, 1, 2, 1, 3, 2, 1},
   {0, 0, 1, 2, 1, 3, 3, 2},
   {1, 2, 1, 2, 1, 3, 3, 2},
   {2, 0, 1, 2, 1, 3, 3, 3},
   {3, 2, 1, 2, 1, 3, 3, 3},
   {0, 0, 1, 3, 2, 0, 0, 0},
   {1, 2, 1, 3, 2, 0, 0, 0},
   {2, 0, 1, 3, 2, 0, 0, 1},
   {3, 2, 1, 3, 2, 0, 0, 1},
   {0, 0, 1, 3, 2, 0, 1, 2},
   {1, 2, 1, 3, 2, 0, 1, 2},
   {2, 0, 1, 3, 2, 0, 1, 3},
   {3, 2, 1, 3, 2, 0, 1, 3},
   {0, 0, 1, 3, 2, 1, 2, 0},
   {1, 2, 1, 3, 2, 1, 2, 0},
   {2, 0, 1, 3, 2, 1, 2, 1},
   {3, 2, 1, 3, 2, 1, 2, 1},
   {0, 0, 1, 3, 2, 1, 3, 2},
   {1, 2, 1, 3, 2, 1, 3, 2},
   {2, 0, 1, 3, 2, 1, 3, 3},
   {3, 2, 1, 3, 2, 1, 3, 3},
   {0, 0, 1, 3, 3, 2, 0, 0},
   {1, 2, 1, 3, 3, 2, 0, 0},
   {2, 0, 1, 3, 3, 2, 0, 1},
   {3, 2, 1, 3, 3, 2, 0, 1},
   {0, 0, 1, 3, 3, 2, 1, 2},
   {1, 2, 1, 3, 3, 2, 1, 2},
   {2, 0, 1, 3, 3, 2, 1, 3},
   {3, 2, 1, 3, 3, 2, 1, 3},
   {0, 0, 1, 3, 3, 3, 2, 0},
   {1, 2, 1, 3, 3, 3, 2, 0},
   {2, 0, 1, 3, 3, 3, 2, 1},
   {3, 2, 1, 3, 3, 3, 2, 1},
   {0, 0, 1, 3, 3, 3, 3, 2},
   {1, 2, 1, 3, 3, 3, 3, 2},
   {2, 0, 1, 3, 3, 3, 3, 3},
   {3, 2, 1, 3, 3, 3, 3, 3},
   {0, 1, 2, 0, 0, 0, 0, 0},
   {1, 3, 2, 0, 0, 0, 0, 0},
   {2, 1, 2, 0, 0, 0, 0, 1},
   {3, 3, 2, 0, 0, 0, 0, 1},
   {0, 1, 2, 0, 0, 0, 1, 2},
   {1, 3, 2, 0, 0, 0, 1, 2},
   {2, 1, 2, 0, 0, 0, 1, 3},
   {3, 3, 2, 0, 0, 0, 1, 3},
   {0, 1, 2, 0, 0, 1, 2, 0},
   {1, 3, 2, 0, 0, 1, 2, 0},
   {2, 1, 2, 0, 0, 1, 2, 1},
   {3, 3, 2, 0, 0, 1, 2, 1},
   {0, 1, 2, 0, 0, 1, 3, 2},
   {1, 3, 2, 0, 0, 1, 3, 2},
   {2, 1, 2, 0, 0, 1, 3, 3},
   {3, 3, 2, 0, 0, 1, 3, 3},
   {0, 1, 2, 0, 1, 2, 0, 0},
   {1, 3, 2, 0, 1, 2, 0, 0},
   {2, 1, 2, 0, 1, 2, 0, 1},
   {3, 3, 2, 0, 1, 2, 0, 1},
   {0, 1, 2, 0, 1, 2, 1, 2},
   {1, 3, 2, 0, 1, 2, 1, 2},
   {2, 1, 2, 0, 1, 2, 1, 3},
   {3, 3, 2, 0, 1, 2, 1, 3},
   {0, 1, 2, 0, 1, 3, 2, 0},
   {1, 3, 2, 0, 1, 3, 2, 0},
   {2, 1, 2, 0, 1, 3, 2, 1},
   {3, 3, 2, 0, 1, 3, 2, 1},
   {0, 1, 2, 0, 1, 3, 3, 2},
   {1, 3, 2, 0, 1, 3, 3, 2},
   {2, 1, 2, 0, 1, 3, 3, 3},
   {3, 3, 2, 0, 1, 3, 3, 3},
   {0, 1, 2, 1, 2, 0, 0, 0},
   {1, 3, 2, 1, 2, 0, 0, 0},
   {2, 1, 2, 1, 2, 0, 0, 1},
   {3, 3, 2, 1, 2, 0, 0, 1},
   {0, 1, 2, 1, 2, 0, 1, 2},
   {1, 3, 2, 1, 2, 0, 1, 2},
   {2, 1, 2, 1, 2, 0, 1, 3},
   {3, 3, 2, 1, 2, 0, 1, 3},
   {0, 1, 2, 1, 2, 1, 2, 0},
   {1, 3, 2, 1, 2, 1, 2, 0},
   {2, 1, 2, 1, 2, 1, 2, 1},
   {3, 3, 2, 1, 2, 1, 2, 1},
   {0, 1, 2, 1, 2, 1, 3, 2},
   {1, 3, 2, 1, 2, 1, 3, 2},
   {2, 1, 2, 1, 2, 1, 3, 3},
   {3, 3, 2, 1, 2, 1, 3, 3},
   {0, 1, 2, 1, 3, 2, 0, 0},
   {1, 3, 2, 1, 3, 2, 0, 0},
   {2, 1, 2, 1, 3, 2, 0, 1},
   {3, 3, 2, 1, 3, 2, 0, 1},
   {0, 1, 2, 1, 3, 2, 1, 2},
   {1, 3, 2, 1, 3, 2, 1, 2},
   {2, 1, 2, 1, 3, 2, 1, 3},
   {3, 3, 2, 1, 3, 2, 1, 3},
   {0, 1, 2, 1, 3, 3, 2, 0},
   {1, 3, 2, 1, 3, 3, 2, 0},
   {2, 1, 2, 1, 3, 3, 2, 1},
   {3, 3, 2, 1, 3, 3, 2, 1},
   {0, 1, 2, 1, 3, 3, 3, 2},
   {1, 3, 2, 1, 3, 3, 3, 2},
   {2, 1, 2, 1, 3, 3, 3, 3},
   {3, 3, 2, 1, 3, 3, 3, 3},
   {0, 1, 3, 2, 0, 0, 0, 0},
   {1, 3, 3, 2, 0, 0, 0, 0},
   {2, 1, 3, 2, 0, 0, 0, 1},
   {3, 3, 3, 2, 0, 0, 0, 1},
   {0, 1, 3, 2, 0, 0, 1, 2},
   {1, 3, 3, 2, 0, 0, 1, 2},
   {2, 1, 3, 2, 0, 0, 1, 3},
   {3, 3, 3, 2, 0, 0, 1, 3},
   {0, 1, 3, 2, 0, 1, 2, 0},
   {1, 3, 3, 2, 0, 1, 2, 0},
   {2, 1, 3, 2, 0, 1, 2, 1},
   {3, 3, 3, 2, 0, 1, 2, 1},
   {0, 1, 3, 2, 0, 1, 3, 2},
   {1, 3, 3, 2, 0, 1, 3, 2},
   {2, 1, 3, 2, 0, 1, 3, 3},
   {3, 3, 3, 2, 0, 1, 3, 3},
   {0, 1, 3, 2, 1, 2, 0, 0},
   {1, 3, 3, 2, 1, 2, 0, 0},
   {2, 1, 3, 2, 1, 2, 0, 1},
   {3, 3, 3, 2, 1, 2, 0, 1},
   {0, 1, 3, 2, 1, 2, 1, 2},
   {1, 3, 3, 2, 1, 2, 1, 2},
   {2, 1, 3, 2, 1, 2, 1, 3},
   {3, 3, 3, 2, 1, 2, 1, 3},
   {0, 1, 3, 2, 1, 3, 2, 0},
   {1, 3, 3, 2, 1, 3, 2, 0},
   {2, 1, 3, 2, 1, 3, 2, 1},
   {3, 3, 3, 2, 1, 3, 2, 1},
   {0, 1, 3, 2, 1, 3, 3, 2},
   {1, 3, 3, 2, 1, 3, 3, 2},
   {2, 1, 3, 2, 1, 3, 3, 3},
   {3, 3, 3, 2, 1, 3, 3, 3},
   {0, 1, 3, 3, 2, 0, 0, 0},
   {1, 3, 3, 3, 2, 0, 0, 0},
   {2, 1, 3, 3, 2, 0, 0, 1},
   {3, 3, 3, 3, 2, 0, 0, 1},
   {0, 1, 3, 3, 2, 0, 1, 2},
   {1, 3, 3, 3, 2, 0, 1, 2},
   {2, 1, 3, 3, 2, 0, 1, 3},
   {3, 3, 3, 3, 2, 0, 1, 3},
   {0, 1, 3, 3, 2, 1, 2, 0},
   {1, 3, 3, 3, 2, 1, 2, 0},
   {2, 1, 3, 3, 2, 1, 2, 1},
   {3, 3, 3, 3, 2, 1, 2, 1},
   {0, 1, 3, 3, 2, 1, 3, 2},
   {1, 3, 3, 3, 2, 1, 3, 2},
   {2, 1, 3, 3, 2, 1, 3, 3},
   {3, 3, 3, 3, 2, 1, 3, 3},
   {0, 1, 3, 3, 3, 2, 0, 0},
   {1, 3, 3, 3, 3, 2, 0, 0},
   {2, 1, 3, 3, 3, 2, 0, 1},
   {3, 3, 3, 3, 3, 2, 0, 1},
   {0, 1, 3, 3, 3, 2, 1, 2},
   {1, 3, 3, 3, 3, 2, 1, 2},
   {2, 1, 3, 3, 3, 2, 1, 3},
   {3, 3, 3, 3, 3, 2, 1, 3},
   {0, 1, 3, 3, 3, 3, 2, 0},
   {1, 3, 3, 3, 3, 3, 2, 0},
   {2, 1, 3, 3, 3, 3, 2, 1},
   {3, 3, 3, 3, 3, 3, 2, 1},
   {0, 1, 3, 3, 3, 3, 3, 2},
   {1, 3, 3, 3, 3, 3, 3, 2},
   {2, 1, 3, 3, 3, 3, 3, 3},
   {3, 3, 3, 3, 3, 3, 3, 3},
};

// Returns 0xffffffff if (v1 < v2) and 0x00000000 otherwise
static __inline__ uint64_t lower_than_mask(const uint64_t v1, const uint64_t v2) {
  return -((v1 - v2) >> 63);
}

void poly_frommsg(poly *r, const uint8_t msg[KYBER_INDCPA_MSGBYTES])
{
  unsigned int i;

#if (KYBER_INDCPA_MSGBYTES != KYBER_N/8)
#error "KYBER_INDCPA_MSGBYTES must be equal to KYBER_N/8 bytes!"
#endif

  size_t base = 0;
  for(i=0;i<KYBER_N/8;i++) {
    uint16_t b7 = -((msg[i] >> 0) & 1);
    uint16_t b6 = -((msg[i] >> 1) & 1);
    uint16_t b5 = -((msg[i] >> 2) & 1);
    uint16_t b4 = -((msg[i] >> 3) & 1);
    uint16_t b3 = -((msg[i] >> 4) & 1);
    uint16_t b2 = -((msg[i] >> 5) & 1);
    uint16_t b1 = -((msg[i] >> 6) & 1);
    uint16_t b0 = -((msg[i] >> 7) & 1);

    r->coeffs[base]               = (CODE_ALPHA & b6) + (CODE_BETA & b7);
    r->coeffs[base + 1*KYBER_N/8] = (CODE_ALPHA & b7) + (CODE_BETA & b0);
    r->coeffs[base + 2*KYBER_N/8] = (CODE_ALPHA & b0) + (CODE_BETA & b1);
    r->coeffs[base + 3*KYBER_N/8] = (CODE_ALPHA & b1) + (CODE_BETA & b2);

    r->coeffs[base + 4*KYBER_N/8] = (CODE_ALPHA & b2) + (CODE_BETA & b3);
    r->coeffs[base + 5*KYBER_N/8] = (CODE_ALPHA & b3) + (CODE_BETA & b4);
    r->coeffs[base + 6*KYBER_N/8] = (CODE_ALPHA & b4) + (CODE_BETA & b5);
    r->coeffs[base + 7*KYBER_N/8] = (CODE_ALPHA & b5) + (CODE_BETA & b6);

    base++;
  }
}

// Returns `abs(centered_mod(value, KYBER_Q)` assuming `-KYBER_Q <= value <= KYBER_Q`
static __inline__ uint64_t abs_center_mod_of_2q_centered_value(int64_t value) {
  // value = abs(value):
  uint64_t mask_sign = value >> 63;
  value ^= mask_sign;
  value += mask_sign & 1;
  value -= KYBER_Q & lower_than_mask(KYBER_Q/2, value);
  return value;
}

static __inline__ uint64_t get_distance_sqr_to_codeword(uint16_t idx, uint64_t distsqr_matrix[8][4]) {
  uint64_t a0 = distsqr_matrix[0][CODEWORDS[idx][0]];
  uint64_t a1 = distsqr_matrix[1][CODEWORDS[idx][1]];
  uint64_t a2 = distsqr_matrix[2][CODEWORDS[idx][2]];
  uint64_t a3 = distsqr_matrix[3][CODEWORDS[idx][3]];
  uint64_t a4 = distsqr_matrix[4][CODEWORDS[idx][4]];
  uint64_t a5 = distsqr_matrix[5][CODEWORDS[idx][5]];
  uint64_t a6 = distsqr_matrix[6][CODEWORDS[idx][6]];
  uint64_t a7 = distsqr_matrix[7][CODEWORDS[idx][7]];

  // Returns `distance_sqr | codeword_index`
  return (a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7) << 8 | idx;
}

// Computes min(v1, v2) in constant time
static __inline__ uint32_t secure_min(int64_t v1, int64_t v2) {
  uint64_t mask_min_v1 = lower_than_mask(v1, v2);
  return (mask_min_v1 & v1) | (~mask_min_v1 & v2);
}

static __inline__ int decode_msg_8d(int16_t target[4]) {
  // This function assumes that:
  //    * -q/2 < target[i] <= q/2
  //    * -q/2 < CODEWORD[idx][i] <= q/2

  // Build matrix with square distances to target coordinates
  uint64_t distsqr_matrix[8][4] = {0};
  for (int i = 0; i < 8; i++) {
    for (int j = 0; j < 4; j++) {
      distsqr_matrix[i][j] = abs_center_mod_of_2q_centered_value(target[i] - CODE_VALS[j]);
      distsqr_matrix[i][j] *= distsqr_matrix[i][j];
    }
  }

  uint64_t min_dist_codeword = get_distance_sqr_to_codeword(0, distsqr_matrix);
  for (size_t i = 1; i < 256; i++) {
    min_dist_codeword = secure_min(get_distance_sqr_to_codeword(i, distsqr_matrix), min_dist_codeword);
  }

  // Retrieves `codeword_index` from `distance_sqr | codeword_index`
  return min_dist_codeword & 0xFF;
}

/*************************************************
* Name:        poly_tomsg
*
* Description: Convert polynomial to 32-byte message
*
* Arguments:   - uint8_t *msg: pointer to output message
*              - const poly *a: pointer to input polynomial
**************************************************/
void poly_tomsg(uint8_t msg[KYBER_INDCPA_MSGBYTES], const poly *a)
{
  unsigned int i;
  uint32_t base = 0;
  for(i=0;i<KYBER_N/8;i++) {
    msg[i] = 0;

    int16_t target[8] = {0};
    target[0] = a->coeffs[base];
    target[1] = a->coeffs[base + 1*KYBER_N/8];
    target[2] = a->coeffs[base + 2*KYBER_N/8];
    target[3] = a->coeffs[base + 3*KYBER_N/8];
    target[4] = a->coeffs[base + 4*KYBER_N/8];
    target[5] = a->coeffs[base + 5*KYBER_N/8];
    target[6] = a->coeffs[base + 6*KYBER_N/8];
    target[7] = a->coeffs[base + 7*KYBER_N/8];
    msg[i] |= decode_msg_8d(target);

    base += 1;
  }
}


#else
# error "KYBER_N_ENCODING_DIMENSIONS must be either 1, 2, 4 or 8"
#endif // !defined(KYBER_N_ENCODING_DIMENSIONS) || (KYBER_N_ENCODING_DIMENSIONS == 1)

/*************************************************
* Name:        poly_getnoise_eta1
*
* Description: Sample a polynomial deterministically from a seed and a nonce,
*              with output polynomial close to centered binomial distribution
*              with parameter KYBER_ETA1
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *seed: pointer to input seed
*                                     (of length KYBER_SYMBYTES bytes)
*              - uint8_t nonce: one-byte input nonce
**************************************************/
void poly_getnoise_eta1(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce)
{
  uint8_t buf[KYBER_ETA1*KYBER_N/4];
  prf(buf, sizeof(buf), seed, nonce);
  poly_cbd_eta1(r, buf);
}

/*************************************************
* Name:        poly_getnoise_eta2
*
* Description: Sample a polynomial deterministically from a seed and a nonce,
*              with output polynomial close to centered binomial distribution
*              with parameter KYBER_ETA2
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *seed: pointer to input seed
*                                     (of length KYBER_SYMBYTES bytes)
*              - uint8_t nonce: one-byte input nonce
**************************************************/
void poly_getnoise_eta2(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce)
{
  uint8_t buf[KYBER_ETA2*KYBER_N/4];
  prf(buf, sizeof(buf), seed, nonce);
  poly_cbd_eta2(r, buf);
}


/*************************************************
* Name:        poly_ntt
*
* Description: Computes negacyclic number-theoretic transform (NTT) of
*              a polynomial in place;
*              inputs assumed to be in normal order, output in bitreversed order
*
* Arguments:   - uint16_t *r: pointer to in/output polynomial
**************************************************/
void poly_ntt(poly *r)
{
  ntt(r->coeffs);
  poly_reduce(r);
}

/*************************************************
* Name:        poly_invntt_tomont
*
* Description: Computes inverse of negacyclic number-theoretic transform (NTT)
*              of a polynomial in place;
*              inputs assumed to be in bitreversed order, output in normal order
*
* Arguments:   - uint16_t *a: pointer to in/output polynomial
**************************************************/
void poly_invntt_tomont(poly *r)
{
  invntt(r->coeffs);
}

/*************************************************
* Name:        poly_basemul_montgomery
*
* Description: Multiplication of two polynomials in NTT domain
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_basemul_montgomery(poly *r, const poly *a, const poly *b)
{
  unsigned int i;
  for(i=0;i<KYBER_N/4;i++) {
    basemul(&r->coeffs[4*i], &a->coeffs[4*i], &b->coeffs[4*i], zetas[64+i]);
    basemul(&r->coeffs[4*i+2], &a->coeffs[4*i+2], &b->coeffs[4*i+2], -zetas[64+i]);
  }
}

/*************************************************
* Name:        poly_tomont
*
* Description: Inplace conversion of all coefficients of a polynomial
*              from normal domain to Montgomery domain
*
* Arguments:   - poly *r: pointer to input/output polynomial
**************************************************/
void poly_tomont(poly *r)
{
  unsigned int i;
  const int16_t f = (1ULL << 32) % KYBER_Q;
  for(i=0;i<KYBER_N;i++)
    r->coeffs[i] = montgomery_reduce((int32_t)r->coeffs[i]*f);
}

/*************************************************
* Name:        poly_reduce
*
* Description: Applies Barrett reduction to all coefficients of a polynomial
*              for details of the Barrett reduction see comments in reduce.c
*
* Arguments:   - poly *r: pointer to input/output polynomial
**************************************************/
void poly_reduce(poly *r)
{
  unsigned int i;
  for(i=0;i<KYBER_N;i++)
    r->coeffs[i] = barrett_reduce(r->coeffs[i]);
}

/*************************************************
* Name:        poly_add
*
* Description: Add two polynomials; no modular reduction is performed
*
* Arguments: - poly *r: pointer to output polynomial
*            - const poly *a: pointer to first input polynomial
*            - const poly *b: pointer to second input polynomial
**************************************************/
void poly_add(poly *r, const poly *a, const poly *b)
{
  unsigned int i;
  for(i=0;i<KYBER_N;i++)
    r->coeffs[i] = a->coeffs[i] + b->coeffs[i];
}

/*************************************************
* Name:        poly_sub
*
* Description: Subtract two polynomials; no modular reduction is performed
*
* Arguments: - poly *r:       pointer to output polynomial
*            - const poly *a: pointer to first input polynomial
*            - const poly *b: pointer to second input polynomial
**************************************************/
void poly_sub(poly *r, const poly *a, const poly *b)
{
  unsigned int i;
  for(i=0;i<KYBER_N;i++)
    r->coeffs[i] = a->coeffs[i] - b->coeffs[i];
}
