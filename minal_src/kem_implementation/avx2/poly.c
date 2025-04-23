#include <stdint.h>
#include <immintrin.h>
#include <string.h>
#include "align.h"
#include "params.h"
#include "poly.h"
#include "ntt.h"
#include "consts.h"
#include "reduce.h"
#include "cbd.h"
#include "symmetric.h"

/*************************************************
* Name:        poly_compress
*
* Description: Compression and subsequent serialization of a polynomial.
*              The coefficients of the input polynomial are assumed to
*              lie in the invertal [0,q], i.e. the polynomial must be reduced
*              by poly_reduce().
*
* Arguments:   - uint8_t *r: pointer to output byte array
*                            (of length KYBER_POLYCOMPRESSEDBYTES)
*              - const poly *a: pointer to input polynomial
**************************************************/
#if (KYBER_POLYCOMPRESSEDBYTES == 96)
void poly_compress(uint8_t r[96], const poly * restrict a)
{
  unsigned int i;
  __m256i f0, f1, f2, f3;
  __m128i t0, t1;
  const __m256i v = _mm256_load_si256(&qdata.vec[_16XV/16]);
  const __m256i shift1 = _mm256_set1_epi16(1 << 8);
  const __m256i mask = _mm256_set1_epi16(7);
  const __m256i shift2 = _mm256_set1_epi16((8 << 8) + 1);
  const __m256i shift3 = _mm256_set1_epi32((64 << 16) + 1);
  const __m256i sllvdidx = _mm256_set1_epi64x(12LL << 32);
  const __m256i shufbidx = _mm256_set_epi8( 8, 2, 1, 0,-1,-1,-1,-1,14,13,12, 6, 5, 4,10, 9,
                                           -1,-1,-1,-1,14,13,12, 6, 5, 4,10, 9, 8, 2, 1, 0);

  for(i=0;i<KYBER_N/64;i++) {
    f0 = _mm256_load_si256(&a->vec[4*i+0]);
    f1 = _mm256_load_si256(&a->vec[4*i+1]);
    f2 = _mm256_load_si256(&a->vec[4*i+2]);
    f3 = _mm256_load_si256(&a->vec[4*i+3]);
    f0 = _mm256_mulhi_epi16(f0,v);
    f1 = _mm256_mulhi_epi16(f1,v);
    f2 = _mm256_mulhi_epi16(f2,v);
    f3 = _mm256_mulhi_epi16(f3,v);
    f0 = _mm256_mulhrs_epi16(f0,shift1);
    f1 = _mm256_mulhrs_epi16(f1,shift1);
    f2 = _mm256_mulhrs_epi16(f2,shift1);
    f3 = _mm256_mulhrs_epi16(f3,shift1);
    f0 = _mm256_and_si256(f0,mask);
    f1 = _mm256_and_si256(f1,mask);
    f2 = _mm256_and_si256(f2,mask);
    f3 = _mm256_and_si256(f3,mask);
    f0 = _mm256_packus_epi16(f0,f1);
    f2 = _mm256_packus_epi16(f2,f3);
    f0 = _mm256_maddubs_epi16(f0,shift2);	// a0 a1 a2 a3 b0 b1 b2 b3 a4 a5 a6 a7 b4 b5 b6 b7
    f2 = _mm256_maddubs_epi16(f2,shift2);	// c0 c1 c2 c3 d0 d1 d2 d3 c4 c5 c6 c7 d4 d5 d6 d7
    f0 = _mm256_madd_epi16(f0,shift3);		// a0 a1 b0 b1 a2 a3 b2 b3
    f2 = _mm256_madd_epi16(f2,shift3);		// c0 c1 d0 d1 c2 c3 d2 d3
    f0 = _mm256_sllv_epi32(f0,sllvdidx);
    f2 = _mm256_sllv_epi32(f2,sllvdidx);
    f0 = _mm256_hadd_epi32(f0,f2);		// a0 c0 c0 d0 a1 b1 c1 d1
    f0 = _mm256_permute4x64_epi64(f0,0xD8);	// a0 b0 a1 b1 c0 d0 c1 d1
    f0 = _mm256_shuffle_epi8(f0,shufbidx);
    t0 = _mm256_castsi256_si128(f0);
    t1 = _mm256_extracti128_si256(f0,1);
    t0 = _mm_blend_epi32(t0,t1,0x08);
    _mm_storeu_si128((__m128i *)&r[24*i+ 0],t0);
    _mm_storel_epi64((__m128i *)&r[24*i+16],t1);
  }
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
void poly_decompress(poly * restrict r, const uint8_t a[96])
{
  unsigned int i;
  __m128i t;
  __m256i f;
  const __m256i q = _mm256_load_si256(&qdata.vec[_16XQ/16]);
  const __m256i shufbidx = _mm256_set_epi8(5,5,5,5,5,4,4,4,4,4,4,3,3,3,3,3,
                                           2,2,2,2,2,1,1,1,1,1,1,0,0,0,0,0);
  const __m256i mask = _mm256_set_epi16(224,28,896,112,14,448,56,7,
                                        224,28,896,112,14,448,56,7);
  const __m256i shift = _mm256_set_epi16(128,1024,32,256,2048,64,512,4096,
                                         128,1024,32,256,2048,64,512,4096);

  for(i=0;i<KYBER_N/16;i++) {
    t = _mm_castps_si128(_mm_load_ss((float *)&a[6*i+0])));
    t = _mm_insert_epi16(t,*(int16_t *)&a[6*i+4],2);
    f = _mm256_broadcastsi128_si256(t);
    f = _mm256_blend_epi16(f,g,0x);
    f = _mm256_shuffle_epi8(f,shufbidx);
    f = _mm256_and_si256(f,mask);
    f = _mm256_mullo_epi16(f,shift);
    f = _mm256_mulhrs_epi16(f,q);
    _mm256_store_si256(&r->vec[i],f);
  }
}

#elif (KYBER_POLYCOMPRESSEDBYTES == 128)
void poly_compress(uint8_t r[128], const poly * restrict a)
{
  unsigned int i;
  __m256i f0, f1, f2, f3;
  const __m256i v = _mm256_load_si256(&qdata.vec[_16XV/16]);
  const __m256i shift1 = _mm256_set1_epi16(1 << 9);
  const __m256i mask = _mm256_set1_epi16(15);
  const __m256i shift2 = _mm256_set1_epi16((16 << 8) + 1);
  const __m256i permdidx = _mm256_set_epi32(7,3,6,2,5,1,4,0);

  for(i=0;i<KYBER_N/64;i++) {
    f0 = _mm256_load_si256(&a->vec[4*i+0]);
    f1 = _mm256_load_si256(&a->vec[4*i+1]);
    f2 = _mm256_load_si256(&a->vec[4*i+2]);
    f3 = _mm256_load_si256(&a->vec[4*i+3]);
    f0 = _mm256_mulhi_epi16(f0,v);
    f1 = _mm256_mulhi_epi16(f1,v);
    f2 = _mm256_mulhi_epi16(f2,v);
    f3 = _mm256_mulhi_epi16(f3,v);
    f0 = _mm256_mulhrs_epi16(f0,shift1);
    f1 = _mm256_mulhrs_epi16(f1,shift1);
    f2 = _mm256_mulhrs_epi16(f2,shift1);
    f3 = _mm256_mulhrs_epi16(f3,shift1);
    f0 = _mm256_and_si256(f0,mask);
    f1 = _mm256_and_si256(f1,mask);
    f2 = _mm256_and_si256(f2,mask);
    f3 = _mm256_and_si256(f3,mask);
    f0 = _mm256_packus_epi16(f0,f1);
    f2 = _mm256_packus_epi16(f2,f3);
    f0 = _mm256_maddubs_epi16(f0,shift2);
    f2 = _mm256_maddubs_epi16(f2,shift2);
    f0 = _mm256_packus_epi16(f0,f2);
    f0 = _mm256_permutevar8x32_epi32(f0,permdidx);
    _mm256_storeu_si256((__m256i *)&r[32*i],f0);
  }
}

void poly_decompress(poly * restrict r, const uint8_t a[128])
{
  unsigned int i;
  __m128i t;
  __m256i f;
  const __m256i q = _mm256_load_si256(&qdata.vec[_16XQ/16]);
  const __m256i shufbidx = _mm256_set_epi8(7,7,7,7,6,6,6,6,5,5,5,5,4,4,4,4,
                                           3,3,3,3,2,2,2,2,1,1,1,1,0,0,0,0);
  const __m256i mask = _mm256_set1_epi32(0x00F0000F);
  const __m256i shift = _mm256_set1_epi32((128 << 16) + 2048);

  for(i=0;i<KYBER_N/16;i++) {
    t = _mm_loadl_epi64((__m128i *)&a[8*i]);
    f = _mm256_broadcastsi128_si256(t);
    f = _mm256_shuffle_epi8(f,shufbidx);
    f = _mm256_and_si256(f,mask);
    f = _mm256_mullo_epi16(f,shift);
    f = _mm256_mulhrs_epi16(f,q);
    _mm256_store_si256(&r->vec[i],f);
  }
}

#elif (KYBER_POLYCOMPRESSEDBYTES == 160)
void poly_compress(uint8_t r[160], const poly * restrict a)
{
  unsigned int i;
  __m256i f0, f1;
  __m128i t0, t1;
  const __m256i v = _mm256_load_si256(&qdata.vec[_16XV/16]);
  const __m256i shift1 = _mm256_set1_epi16(1 << 10);
  const __m256i mask = _mm256_set1_epi16(31);
  const __m256i shift2 = _mm256_set1_epi16((32 << 8) + 1);
  const __m256i shift3 = _mm256_set1_epi32((1024 << 16) + 1);
  const __m256i sllvdidx = _mm256_set1_epi64x(12);
  const __m256i shufbidx = _mm256_set_epi8( 8,-1,-1,-1,-1,-1, 4, 3, 2, 1, 0,-1,12,11,10, 9,
                                           -1,12,11,10, 9, 8,-1,-1,-1,-1,-1 ,4, 3, 2, 1, 0);

  for(i=0;i<KYBER_N/32;i++) {
    f0 = _mm256_load_si256(&a->vec[2*i+0]);
    f1 = _mm256_load_si256(&a->vec[2*i+1]);
    f0 = _mm256_mulhi_epi16(f0,v);
    f1 = _mm256_mulhi_epi16(f1,v);
    f0 = _mm256_mulhrs_epi16(f0,shift1);
    f1 = _mm256_mulhrs_epi16(f1,shift1);
    f0 = _mm256_and_si256(f0,mask);
    f1 = _mm256_and_si256(f1,mask);
    f0 = _mm256_packus_epi16(f0,f1);
    f0 = _mm256_maddubs_epi16(f0,shift2);	// a0 a1 a2 a3 b0 b1 b2 b3 a4 a5 a6 a7 b4 b5 b6 b7
    f0 = _mm256_madd_epi16(f0,shift3);		// a0 a1 b0 b1 a2 a3 b2 b3
    f0 = _mm256_sllv_epi32(f0,sllvdidx);
    f0 = _mm256_srlv_epi64(f0,sllvdidx);
    f0 = _mm256_shuffle_epi8(f0,shufbidx);
    t0 = _mm256_castsi256_si128(f0);
    t1 = _mm256_extracti128_si256(f0,1);
    t0 = _mm_blendv_epi8(t0,t1,_mm256_castsi256_si128(shufbidx));
    _mm_storeu_si128((__m128i *)&r[20*i+ 0],t0);
    memcpy(&r[20*i+16],&t1,4);
  }
}

void poly_decompress(poly * restrict r, const uint8_t a[160])
{
  unsigned int i;
  __m128i t;
  __m256i f;
  int16_t ti;
  const __m256i q = _mm256_load_si256(&qdata.vec[_16XQ/16]);
  const __m256i shufbidx = _mm256_set_epi8(9,9,9,8,8,8,8,7,7,6,6,6,6,5,5,5,
                                           4,4,4,3,3,3,3,2,2,1,1,1,1,0,0,0);
  const __m256i mask = _mm256_set_epi16(248,1984,62,496,3968,124,992,31,
                                        248,1984,62,496,3968,124,992,31);
  const __m256i shift = _mm256_set_epi16(128,16,512,64,8,256,32,1024,
                                         128,16,512,64,8,256,32,1024);

  for(i=0;i<KYBER_N/16;i++) {
    t = _mm_loadl_epi64((__m128i *)&a[10*i+0]);
    memcpy(&ti,&a[10*i+8],2);
    t = _mm_insert_epi16(t,ti,4);
    f = _mm256_broadcastsi128_si256(t);
    f = _mm256_shuffle_epi8(f,shufbidx);
    f = _mm256_and_si256(f,mask);
    f = _mm256_mullo_epi16(f,shift);
    f = _mm256_mulhrs_epi16(f,q);
    _mm256_store_si256(&r->vec[i],f);
  }
}

#endif

/*************************************************
* Name:        poly_tobytes
*
* Description: Serialization of a polynomial in NTT representation.
*              The coefficients of the input polynomial are assumed to
*              lie in the invertal [0,q], i.e. the polynomial must be reduced
*              by poly_reduce(). The coefficients are orderd as output by
*              poly_ntt(); the serialized output coefficients are in bitreversed
*              order.
*
* Arguments:   - uint8_t *r: pointer to output byte array
*                            (needs space for KYBER_POLYBYTES bytes)
*              - poly *a: pointer to input polynomial
**************************************************/
void poly_tobytes(uint8_t r[KYBER_POLYBYTES], const poly *a)
{
  ntttobytes_avx(r, a->vec, qdata.vec);
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
  nttfrombytes_avx(r->vec, a, qdata.vec);
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
void poly_frommsg(poly * restrict r, const uint8_t msg[KYBER_INDCPA_MSGBYTES])
{
#if (KYBER_INDCPA_MSGBYTES != 32)
#error "KYBER_INDCPA_MSGBYTES must be equal to 32!"
#endif
  __m256i f, g0, g1, g2, g3, h0, h1, h2, h3;
  const __m256i shift = _mm256_broadcastsi128_si256(_mm_set_epi32(0,1,2,3));
  const __m256i idx = _mm256_broadcastsi128_si256(_mm_set_epi8(15,14,11,10,7,6,3,2,13,12,9,8,5,4,1,0));
  const __m256i hqs = _mm256_set1_epi16((KYBER_Q+1)/2);

#define FROMMSG64(i)						\
  g3 = _mm256_shuffle_epi32(f,0x55*i);				\
  g3 = _mm256_sllv_epi32(g3,shift);				\
  g3 = _mm256_shuffle_epi8(g3,idx);				\
  g0 = _mm256_slli_epi16(g3,12);				\
  g1 = _mm256_slli_epi16(g3,8);					\
  g2 = _mm256_slli_epi16(g3,4);					\
  g0 = _mm256_srai_epi16(g0,15);				\
  g1 = _mm256_srai_epi16(g1,15);				\
  g2 = _mm256_srai_epi16(g2,15);				\
  g3 = _mm256_srai_epi16(g3,15);				\
  g0 = _mm256_and_si256(g0,hqs);  /* 19 18 17 16  3  2  1  0 */	\
  g1 = _mm256_and_si256(g1,hqs);  /* 23 22 21 20  7  6  5  4 */	\
  g2 = _mm256_and_si256(g2,hqs);  /* 27 26 25 24 11 10  9  8 */	\
  g3 = _mm256_and_si256(g3,hqs);  /* 31 30 29 28 15 14 13 12 */	\
  h0 = _mm256_unpacklo_epi64(g0,g1);				\
  h2 = _mm256_unpackhi_epi64(g0,g1);				\
  h1 = _mm256_unpacklo_epi64(g2,g3);				\
  h3 = _mm256_unpackhi_epi64(g2,g3);				\
  g0 = _mm256_permute2x128_si256(h0,h1,0x20);			\
  g2 = _mm256_permute2x128_si256(h0,h1,0x31);			\
  g1 = _mm256_permute2x128_si256(h2,h3,0x20);			\
  g3 = _mm256_permute2x128_si256(h2,h3,0x31);			\
  _mm256_store_si256(&r->vec[0+2*i+0],g0);	\
  _mm256_store_si256(&r->vec[0+2*i+1],g1);	\
  _mm256_store_si256(&r->vec[8+2*i+0],g2);	\
  _mm256_store_si256(&r->vec[8+2*i+1],g3)

  f = _mm256_loadu_si256((__m256i *)msg);
  FROMMSG64(0);
  FROMMSG64(1);
  FROMMSG64(2);
  FROMMSG64(3);
}

/*************************************************
* Name:        poly_tomsg
*
* Description: Convert polynomial to 32-byte message.
*              The coefficients of the input polynomial are assumed to
*              lie in the invertal [0,q], i.e. the polynomial must be reduced
*              by poly_reduce().
*
* Arguments:   - uint8_t *msg: pointer to output message
*              - poly *a: pointer to input polynomial
**************************************************/
void poly_tomsg(uint8_t msg[KYBER_INDCPA_MSGBYTES], const poly * restrict a)
{
  unsigned int i;
  uint32_t small;
  __m256i f0, f1, g0, g1;
  const __m256i hq = _mm256_set1_epi16((KYBER_Q - 1)/2);
  const __m256i hhq = _mm256_set1_epi16((KYBER_Q - 1)/4);

  for(i=0;i<KYBER_N/32;i++) {
    f0 = _mm256_load_si256(&a->vec[2*i+0]);
    f1 = _mm256_load_si256(&a->vec[2*i+1]);
    f0 = _mm256_sub_epi16(hq, f0);
    f1 = _mm256_sub_epi16(hq, f1);
    g0 = _mm256_srai_epi16(f0, 15);
    g1 = _mm256_srai_epi16(f1, 15);
    f0 = _mm256_xor_si256(f0, g0);
    f1 = _mm256_xor_si256(f1, g1);
    f0 = _mm256_sub_epi16(f0, hhq);
    f1 = _mm256_sub_epi16(f1, hhq);
    f0 = _mm256_packs_epi16(f0, f1);
    f0 = _mm256_permute4x64_epi64(f0, 0xD8);
    small = _mm256_movemask_epi8(f0);
    memcpy(&msg[4*i], &small, 4);
  }
}

#elif defined(KYBER_N_ENCODING_DIMENSIONS) && (KYBER_N_ENCODING_DIMENSIONS == 2)

///////////// GENERATED WITH ../../scripts/generate_decoding_lines.py /////////////////////////////
#if KYBER_K == 2

#define CODE_ALPHA (KYBER_Q/2)
#define CODE_BETA 422

static int32_t L[5][3] = {
    {844, 3330, 5325585},
    {-211, 832, 1438953},
    {3330, 844, -8136105},
    {-1, 1, 5415},
    {-832, 211, 736745},
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
    {844, 3330, 5325585},
    {-211, 832, 1438953},
    {3330, 844, -8136105},
    {-1, 1, 5415},
    {-832, 211, 736745},
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
    {872, 3330, 5220361},
    {-109, 416, 732626},
    {3330, 872, -8124121},
    {-1, 1, 5429},
    {-416, 109, 369874},
};

#define ABOVE_L1(x, y) lower_than_mask(L[0][0]*x + L[0][2], L[0][1]*y)
#define ABOVE_L2(x, y) lower_than_mask(L[1][0]*x + L[1][2], L[1][1]*y)
#define ABOVE_L3(x, y) lower_than_mask(L[2][0]*x + L[2][2], L[2][1]*y)
#define ABOVE_L4(x, y) lower_than_mask(L[3][0]*x + L[3][2], L[3][1]*y)
#define ABOVE_L5(x, y) lower_than_mask(L[4][0]*x + L[4][2], L[4][1]*y)

#endif

///////////////////////////////////////////////////////////////////////////////////////////////////

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

  // Not needed, but conceptually:
  // uint8_t c00 = (~above_l5) || (~above_l2 and ~above_l3) || (above_l4 and above_l1);
  uint8_t c01 = (~above_l1 & ~above_l3 & above_l2);
  uint8_t c10 = (above_l3 & ~above_l2 & above_l5);
  uint8_t c11 = (above_l2 & above_l3 & ~above_l4);

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
# error "KYBER_N_ENCODING_DIMENSIONS must be either 1, 2, 4, or 8"
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
  ALIGNED_UINT8(KYBER_ETA1*KYBER_N/4+32) buf; // +32 bytes as required by poly_cbd_eta1
  prf(buf.coeffs, KYBER_ETA1*KYBER_N/4, seed, nonce);
  poly_cbd_eta1(r, buf.vec);
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
  ALIGNED_UINT8(KYBER_ETA2*KYBER_N/4) buf;
  prf(buf.coeffs, KYBER_ETA2*KYBER_N/4, seed, nonce);
  poly_cbd_eta2(r, buf.vec);
}

#ifndef KYBER_90S
#define NOISE_NBLOCKS ((KYBER_ETA1*KYBER_N/4+SHAKE256_RATE-1)/SHAKE256_RATE)
void poly_getnoise_eta1_4x(poly *r0,
                           poly *r1,
                           poly *r2,
                           poly *r3,
                           const uint8_t seed[32],
                           uint8_t nonce0,
                           uint8_t nonce1,
                           uint8_t nonce2,
                           uint8_t nonce3)
{
  ALIGNED_UINT8(NOISE_NBLOCKS*SHAKE256_RATE) buf[4];
  __m256i f;
  keccakx4_state state;

  f = _mm256_loadu_si256((__m256i *)seed);
  _mm256_store_si256(buf[0].vec, f);
  _mm256_store_si256(buf[1].vec, f);
  _mm256_store_si256(buf[2].vec, f);
  _mm256_store_si256(buf[3].vec, f);

  buf[0].coeffs[32] = nonce0;
  buf[1].coeffs[32] = nonce1;
  buf[2].coeffs[32] = nonce2;
  buf[3].coeffs[32] = nonce3;

  shake256x4_absorb_once(&state, buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 33);
  shake256x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, NOISE_NBLOCKS, &state);

  poly_cbd_eta1(r0, buf[0].vec);
  poly_cbd_eta1(r1, buf[1].vec);
  poly_cbd_eta1(r2, buf[2].vec);
  poly_cbd_eta1(r3, buf[3].vec);
}

#if KYBER_K == 2
void poly_getnoise_eta11_4x_for_null_eta2(poly *r0,
                                          poly *r1,
                                          const uint8_t seed[32],
                                          uint8_t nonce0,
                                          uint8_t nonce1,
                                          uint8_t nonce2,
                                          uint8_t nonce3)
{
  ALIGNED_UINT8(NOISE_NBLOCKS*SHAKE256_RATE) buf[4];
  __m256i f;
  keccakx4_state state;

  f = _mm256_loadu_si256((__m256i *)seed);
  _mm256_store_si256(buf[0].vec, f);
  _mm256_store_si256(buf[1].vec, f);
  _mm256_store_si256(buf[2].vec, f);
  _mm256_store_si256(buf[3].vec, f);

  buf[0].coeffs[32] = nonce0;
  buf[1].coeffs[32] = nonce1;
  buf[2].coeffs[32] = nonce2;
  buf[3].coeffs[32] = nonce3;

  shake256x4_absorb_once(&state, buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 33);
  shake256x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, NOISE_NBLOCKS, &state);

  poly_cbd_eta1(r0, buf[0].vec);
  poly_cbd_eta1(r1, buf[1].vec);
  // poly_cbd_eta2(r2, buf[2].vec);
  // poly_cbd_eta2(r3, buf[3].vec);
}

void poly_getnoise_eta1122_4x(poly *r0,
                              poly *r1,
                              poly *r2,
                              poly *r3,
                              const uint8_t seed[32],
                              uint8_t nonce0,
                              uint8_t nonce1,
                              uint8_t nonce2,
                              uint8_t nonce3)
{
  ALIGNED_UINT8(NOISE_NBLOCKS*SHAKE256_RATE) buf[4];
  __m256i f;
  keccakx4_state state;

  f = _mm256_loadu_si256((__m256i *)seed);
  _mm256_store_si256(buf[0].vec, f);
  _mm256_store_si256(buf[1].vec, f);
  _mm256_store_si256(buf[2].vec, f);
  _mm256_store_si256(buf[3].vec, f);

  buf[0].coeffs[32] = nonce0;
  buf[1].coeffs[32] = nonce1;
  buf[2].coeffs[32] = nonce2;
  buf[3].coeffs[32] = nonce3;

  shake256x4_absorb_once(&state, buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, 33);
  shake256x4_squeezeblocks(buf[0].coeffs, buf[1].coeffs, buf[2].coeffs, buf[3].coeffs, NOISE_NBLOCKS, &state);

  poly_cbd_eta1(r0, buf[0].vec);
  poly_cbd_eta1(r1, buf[1].vec);
  poly_cbd_eta2(r2, buf[2].vec);
  poly_cbd_eta2(r3, buf[3].vec);
}
#endif
#endif

/*************************************************
* Name:        poly_ntt
*
* Description: Computes negacyclic number-theoretic transform (NTT) of
*              a polynomial in place.
*              Input coefficients assumed to be in normal order,
*              output coefficients are in special order that is natural
*              for the vectorization. Input coefficients are assumed to be
*              bounded by q in absolute value, output coefficients are bounded
*              by 16118 in absolute value.
*
* Arguments:   - poly *r: pointer to in/output polynomial
**************************************************/
void poly_ntt(poly *r)
{
  ntt_avx(r->vec, qdata.vec);
}

/*************************************************
* Name:        poly_invntt_tomont
*
* Description: Computes inverse of negacyclic number-theoretic transform (NTT)
*              of a polynomial in place;
*              Input coefficients assumed to be in special order from vectorized
*              forward ntt, output in normal order. Input coefficients can be
*              arbitrary 16-bit integers, output coefficients are bounded by 14870
*              in absolute value.
*
* Arguments:   - poly *a: pointer to in/output polynomial
**************************************************/
void poly_invntt_tomont(poly *r)
{
  invntt_avx(r->vec, qdata.vec);
}

void poly_nttunpack(poly *r)
{
  nttunpack_avx(r->vec, qdata.vec);
}

/*************************************************
* Name:        poly_basemul_montgomery
*
* Description: Multiplication of two polynomials in NTT domain.
*              One of the input polynomials needs to have coefficients
*              bounded by q, the other polynomial can have arbitrary
*              coefficients. Output coefficients are bounded by 6656.
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_basemul_montgomery(poly *r, const poly *a, const poly *b)
{
  basemul_avx(r->vec, a->vec, b->vec, qdata.vec);
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
  tomont_avx(r->vec, qdata.vec);
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
  reduce_avx(r->vec, qdata.vec);
}

/*************************************************
* Name:        poly_add
*
* Description: Add two polynomials. No modular reduction
*              is performed.
*
* Arguments: - poly *r: pointer to output polynomial
*            - const poly *a: pointer to first input polynomial
*            - const poly *b: pointer to second input polynomial
**************************************************/
void poly_add(poly *r, const poly *a, const poly *b)
{
  unsigned int i;
  __m256i f0, f1;

  for(i=0;i<KYBER_N/16;i++) {
    f0 = _mm256_load_si256(&a->vec[i]);
    f1 = _mm256_load_si256(&b->vec[i]);
    f0 = _mm256_add_epi16(f0, f1);
    _mm256_store_si256(&r->vec[i], f0);
  }
}

/*************************************************
* Name:        poly_sub
*
* Description: Subtract two polynomials. No modular reduction
*              is performed.
*
* Arguments: - poly *r: pointer to output polynomial
*            - const poly *a: pointer to first input polynomial
*            - const poly *b: pointer to second input polynomial
**************************************************/
void poly_sub(poly *r, const poly *a, const poly *b)
{
  unsigned int i;
  __m256i f0, f1;

  for(i=0;i<KYBER_N/16;i++) {
    f0 = _mm256_load_si256(&a->vec[i]);
    f1 = _mm256_load_si256(&b->vec[i]);
    f0 = _mm256_sub_epi16(f0, f1);
    _mm256_store_si256(&r->vec[i], f0);
  }
}
