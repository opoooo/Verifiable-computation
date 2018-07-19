#include "md2.h"

/*! \struct vc_k
 * m and d are the row and column of the F passing in and the vck->F going to generate.
 * F is the encrypted function
 * G and T is used to verify
 */
typedef struct
{
  int m;
  int d;
  fq_mat_t F;
  g1_t *G;
  g1_t *T;
} vc_k;


/*! \struct vc_p
 * C is the encrypted x
 * hk is the encrypted y
 */
typedef struct
{
  lhe_c **C;
  fq_poly_t s;
  hash_k *hk;
  gt_t pi;
} vc_p;

/**
 * multify two bn_t numbers and mod the result by an bn_t number
 * @param[out]out        the out number
 * @param[in] in1        a bn_t number
 * @param[in] in2        a bn_t number
 * @param[in] modulus    the bn_t number used to mod
 */
int bn_mul_modz(bn_t out, bn_t in1, bn_t in2, bn_t modulus);

/**
 * generate a random fq_t polynomial
 * @param[out] out       the out polynomial
 * @param[in] deg        the highest degree
 * @param[in] ctx        a parameter used in flint lib
 * @param[in] bound      the bound
 */
int fq_poly_rand(fq_poly_t pol, int deg, fq_ctx_t ctx, bn_t bound);

/**
 * reduce the coefficients of a polynomial
 * @param[out] m         the out polynomial
 * @param[in] par        some parameter in flint lib
 */
int msg_modp(fq_poly_t m, lhe_par *par);

/**
 * generate the verification key
 * @param[out] vck       output the encryped function vck->F and the verification key: vck->G ,vck ->T lhe key K
 * @param[in] F          the input function F
 * @param[in] par        some parameter in flint lib   
 */
int vc_keygen(vc_k *vck, fq_mat_t F, lhe_par *par);

/**
 * generate a random matrix F stored fq_t
 * @param[out]F        the random matrix
 * @param[in] m        the row of the random matrix
 * @param[in] d        the column of the random matrix
 * @param[in] par      some parameter in flint lib
 */
int fq_mat_randz(fq_mat_t F, int m, int d, lhe_par *par);

/**
 * generate the verification key
 * @param[out] vcp       output the encryped message vcp->C and the secret key
 * @param[in] vck        the publick key
 * @param[in] x          the message x
 * @param[in] par        some parameter in flint lib   
 */
int vc_pgen(vc_p *vcp, vc_k * vck, bn_t *x, lhe_par *par);

/**
 * compute the encrypted output y
 * @param[out] nv		 the encrypted output y
 * @param[in] vck        the vck->F is needed here to compute
 * @param[in] vcp        the vcp->C is needed here to compute
 * @param[in] par        some parameter in flint lib   
 */
int vc_comp(lhe_c **nv, vc_k *vck, vc_p *vcp, lhe_par *par);

/**
 * verify the outcome. output 0 if the outcome is true, otherwise 2
 * @param[out] y         the decoded output y
 * @param[in] vck   	 the vck->G and vck->T is needed here to verify
 * @param[in] vcp   	 the vcp->C is needed here to verify
 * @param[in] nv		 the encrypted output y
 * @param[in] par   	 some parameter in flint lib  
 */
int vc_vrfy(fq_poly_t *y, vc_k *vck, vc_p *vcp, lhe_c **nv, lhe_par *par);