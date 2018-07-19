#include "lhep.h"


/*! \struct lhe_f
 * the function F is encrypted as two parts p0 and p1
 */
typedef struct
{
  fq_poly_t p0;
  fq_poly_t p1;
} lhe_c;

/**
 * generate a LHE key
 * @param[out]s          the secret key get from this function
 * @param[in] par        some parameter in flint lib
 */
int lhe_keygen(fq_poly_t s, lhe_par *par);

/**
 * encrype the information
 * @param[out]c          the out polynomial
 * @param[in] s          the secret key
 * @param[in] m          the information wanted to be encryped
 * @param[in] par        some parameter in flint lib
 */
int lhe_enc(lhe_c *c, fq_poly_t s, fq_poly_t m,lhe_par *par);

/**
 * decode the information
 * @param[out]m1         the decoded result
 * @param[in] s          the secret key
 * @param[in] c          the encryped information
 * @param[in] par        some parameter in flint lib
 */
int lhe_dec(fq_poly_t m1, fq_poly_t s, lhe_c *c, lhe_par *par);

/**
 * compute
 * @param[out]c          the computation result
 * @param[in] L          one line of the function matrix with size 1*d
 * @param[in] C          the input message
 * @param[in] d          the column of the L and the row of the C
 * @param[in] par        some parameter in flint lib
 */
int lhe_eval(lhe_c *c, lhe_c **C, fq_t *L, int d, lhe_par *par);