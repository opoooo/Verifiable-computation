#include "md2.h"

/*! \struct pk_a
 * B is the polynomial defined as X^2+b0
 * R is the remainder polynomial
 */
typedef struct
{
  fq_poly_t B;
  fq_poly_t R;
} pk_a;

/*! \struct ek_a
 * d is the degree of polynomial
 * A is the evaluation polynomial
 */
typedef struct
{
  int d;
  fq_poly_t A;
  fq_poly_t Q;
} ek_a;

/*! \struct vk_x
 * B is the polynomial defined as X^2+b0
 * R is the remainder polynomial
 * pi is the proof
 */
typedef struct
{
  g1_t B;
  g2_t R;
  g2_t pi;
} vk_x;

/**
 * generate a random fq_t polynomial
 * @param[out] out       the out polynomial
 * @param[in] deg        the highest degree
 * @param[in] ctx        a parameter used in flint lib
 * @param[in] bound      the bound
 */
int fq_poly_rand(fq_poly_t pol, int deg, fq_ctx_t ctx, bn_t bound);

/**
 * multify two bn_t numbers and mod the result by an bn_t number
 * @param[out]out        the out number
 * @param[in] in1        a bn_t number
 * @param[in] in2        a bn_t number
 * @param[in] modulus    the bn_t number used to mod
 */
int bn_mul_modz(bn_t out, bn_t in1, bn_t in2, bn_t modulus);

/**
 * reduce the coefficients of a polynomial
 * @param[out] m         the out polynomial
 * @param[in] par        some parameter in flint lib
 */
int msg_modp(fq_poly_t m, lhe_par *par);

/**
 * generate the verification key
 * @param[out] vck       the publick key
 * @param[out] eka       the evaluation key
 * @param[in] A          the evaluation polynomial
 * @param[in] par        some parameter in flint lib   
 */
int vc_keygen(pk_a *pka, ek_a *eka, fq_poly_t A, lhe_par *par);

/**
 * generate the verification key
 * @param[out] vkx       the public verification key
 * @param[in] pka        the publick key
 * @param[in] x          the message x
 * @param[in] par        some parameter in flint lib   
 */
int vc_pgen(pk_a *pka, vk_x *vkx, bn_t x, lhe_par *par);

/**
 * generate the verification key
 * @param[out] y         the output y
 * @param[out] vkx       the public verification key
 * @param[in] x          the message
 * @param[in] eka        the evaluation key
 * @param[in] par        some parameter in flint lib   
 */
int vc_comp(g2_t y, vk_x *vkx, bn_t x, ek_a *eka, lhe_par *par);

/**
 * generate the verification key
 * @param[in] y          the output y
 * @param[in] vkx        the public verification key
 * @param[in] par        some parameter in flint lib   
 */
int vc_vrfy(g2_t y, vk_x *vkx, lhe_par *par);