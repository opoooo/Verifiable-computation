#include "poly_eval.h"

// this function generate a random polynomial that has small coefficients
int fq_poly_rand(fq_poly_t pol, int deg, fq_ctx_t ctx, bn_t bound)
{

  fq_poly_init(pol,ctx);

  bn_t z1;
  bn_null(z1);

  fmpz_t z2;
  fmpz_init(z2);

  fq_t z3;
  fq_init(z3,ctx);

  for (int n=0;n<deg;n++)
  {
      bn_rand_mod(z1,bound);
      bn2fmpz(z2,z1);
      fq_set_fmpz(z3,z2,ctx);  
      fq_poly_set_coeff(pol,n,z3,ctx);
  }


  // release the memory
  bn_free(z1);
  fmpz_clear(z2);
  fq_clear(z3,ctx);
  return 0;

}

int bn_mul_modz(bn_t out, bn_t in1, bn_t in2, bn_t modulus)
{
   bn_mul_basic(out,in1,in2);
   bn_mod_basic(out,out,modulus);
   return 0;
}

// this function reduce the coefficients of a polynomial
int msg_modp(fq_poly_t m, lhe_par *par)
{
 
   long len=fq_poly_length(m,par->ctx);

   // fq_t pr;
   // fq_init(pr,par->ctx);
   // fq_set_fmpz(pr,par->pf,par->ctx);
   // fq_print_pretty(pr,par->ctx);

   fq_t zq;
   fq_init(zq,par->ctx);

   fmpz_t zf;
   fmpz_init(zf);

   // fq_poly_print_pretty(m,"x",par->ctx);
   // printf("\n\n");

   for (int i=0;i<len;i++)
   {

     fq_poly_get_coeff(zq,m,i,par->ctx);
     // printf("\n\n");
     // fq_print_pretty(zq,par->ctx);
     fmpz_set_str(zf,fq_get_str_pretty(zq,par->ctx), 10);
     fmpz_mod(zf,zf,par->pf);
     fq_set_fmpz(zq,zf,par->ctx);
     // printf("\n\n");
     // fq_print_pretty(zq,par->ctx);
     fq_poly_set_coeff(m,i,zq,par->ctx);
   }

   // fq_poly_print_pretty(m,"x",par->ctx);
   // printf("\n\n");

   fq_clear(zq,par->ctx);
   fmpz_clear(zf);

   return 0;
 
}


// the vc_keygen algorithm
int vc_keygen(pk_a *pka, ek_a *eka, fq_poly_t A, lhe_par *par)
{
	// b0<-Fp*
	bn_t z1;
	bn_zero(z1);

	fmpz_t z2;
	fmpz_init(z2);

	fq_t b0;
	fq_init(b0,par->ctx);

	while (bn_is_zero(z1)) bn_rand_mod(z1,par->p);	
	bn2fmpz(z2,z1);
    fq_set_fmpz(b0,z2,par->ctx);  

	// B(X) = X^2+b0
	fq_t one;
	fq_init(one,par->ctx);

	fq_t zero;
	fq_init(zero,par->ctx);
	
	fq_poly_t B;	
	fq_poly_init(B,par->ctx);
	
	fq_zero(zero,par->ctx);
	fq_one(one,par->ctx);
	
	fq_poly_set_coeff(B,2,one,par->ctx);
	fq_poly_set_coeff(B,1,zero,par->ctx);
	fq_poly_set_coeff(B,0,b0,par->ctx);	

	// A = BQ+R
	fq_poly_t Q;
	fq_poly_init(Q,par->ctx);
	fq_poly_t R;
	fq_poly_init(R,par->ctx);
	fq_poly_divrem(Q,R,A,B,par->ctx);
	// fq_poly_print_pretty(A,"x",par->ctx);
	// printf("\n\n");
	// fq_poly_print_pretty(B,"x",par->ctx);
	// printf("\n\n");
	// fq_poly_print_pretty(Q,"x",par->ctx);
	// printf("\n\n");
	// fq_poly_print_pretty(R,"x",par->ctx);

	// pka
	fq_poly_init(pka->B,par->ctx);
	fq_poly_init(pka->R,par->ctx);
	fq_poly_set(pka->B,B,par->ctx);
	fq_poly_set(pka->R,R,par->ctx);

	// eka
	fq_poly_init(eka->A,par->ctx);
	fq_poly_init(eka->Q,par->ctx);
	fq_poly_set(eka->A,A,par->ctx);
	fq_poly_set(eka->Q,Q,par->ctx);

	return 0;
}

// the vc_pgen algorithm
int vc_pgen(pk_a *pka, vk_x *vkx, bn_t x, lhe_par *par)
{
	fmpz_t z2;
	fmpz_init(z2);

	fq_t z3;
	fq_init(z3,par->ctx);

	bn2fmpz(z2,x);
    fq_set_fmpz(z3,z2,par->ctx); 

	fq_t z;

	// vkx->B = x^2+b0
	g1_new(vkx->B);
	fq_init(z,par->ctx);
	// fq_poly_evaluate_fq(fq_t rop , const  fq_poly_t f, const fq_t a, const  fq_ctx_t  ctx)
	fq_poly_evaluate_fq(z,pka->B,z3,par->ctx);
	// fq_poly_print_pretty(pka->B,"x",par->ctx);
	// printf("\n\n");
	// fq_print(z,par->ctx);
	// printf("\n\n");
	bn_t temp;
	bn_new(temp);
  	fq2bn(temp,z,par->ctx);
	g1_mul(vkx->B,par->g,temp);
	// g1_print(vkx->B);
	// printf("\n\n");

	// vkx->R = r1*x+r0
	g2_new(vkx->R);
	fq_init(z,par->ctx);
	// fq_poly_evaluate_fq(fq_t rop , const  fq_poly_t f, const fq_t a, const  fq_ctx_t  ctx)
	fq_poly_evaluate_fq(z,pka->R,z3,par->ctx);
	bn_t temp1;
	bn_new(temp1);
  	fq2bn(temp1,z,par->ctx);
	g2_mul(vkx->R,par->h,temp1);
	// g1_print(vkx->R);
	// printf("\n\n");

	return 0;
}

// the vc_comp algorithm
int vc_comp(g2_t y, vk_x *vkx, bn_t x, ek_a *eka, lhe_par *par)
{
	fmpz_t z2;
	fmpz_init(z2);

	fq_t z3;
	fq_init(z3,par->ctx);

	bn2fmpz(z2,x);
    fq_set_fmpz(z3,z2,par->ctx); 

	fq_t z;
	fq_init(z,par->ctx);

	// y = A(x) mod p	
	fq_poly_t m;
	fq_poly_init(m,par->ctx);
	fq_poly_set(m,eka->A,par->ctx);
	// fq_poly_print_pretty(m,"x",par->ctx);
	msg_modp(m,par);
	// printf("\n\n");
	// fq_poly_print_pretty(m,"x",par->ctx);
	// fq_poly_evaluate_fq(fq_t rop , const  fq_poly_t f, const fq_t a, const  fq_ctx_t  ctx)
	fq_poly_evaluate_fq(z,m,z3,par->ctx);
	bn_t temp;
	bn_new(temp);
  	fq2bn(temp,z,par->ctx);
	g2_mul(y,par->h,temp);
	// gt_print(y);

	// vkx->pi = Q(x)
	g2_new(vkx->pi);
	fq_init(z,par->ctx);
	// fq_poly_evaluate_fq(fq_t rop , const  fq_poly_t f, const fq_t a, const  fq_ctx_t  ctx)
	fq_poly_evaluate_fq(z,eka->Q,z3,par->ctx);
	bn_t temp1;
	bn_new(temp1);
  	fq2bn(temp1,z,par->ctx);
	g2_mul(vkx->pi,par->h,temp1);
	// fq_poly_print_pretty(eka->Q,"x",par->ctx);
	// g2_print(par->h);
	// g1_print(vkx->pi);
	// printf("\n\n");

	return 0;
}

// the vc_vrfy algorithm
int vc_vrfy(g2_t y, vk_x *vkx, lhe_par *par)
{
	// y?=B(x)*vkx->pi+R(x)
 	gt_t right;
 	gt_new(right);
    gt_set_unity(right);

    gt_t left;
 	gt_new(left);
    gt_set_unity(left);

    // g2_print(y);
    // printf("\n\n");

    // gt_print(right);
    // printf("\n\n");
    // g1_print(vkx->B);
    // printf("\n\n");
    // g2_print(vkx->R);
    // printf("\n\n");
    // g2_print(vkx->pi);
    // printf("\n\n");

    gt_t temp1;
    gt_new(temp1);
    gt_set_unity(temp1);

    gt_t temp2;
    gt_new(temp2)
    gt_set_unity(temp2);

    pc_map(temp1,vkx->B,vkx->pi);
    pc_map(temp2,par->g,vkx->R);

    gt_mul(right,temp1,temp2);

    // gt_print(right);
    // printf("\n\n");
    
    // g1_print(vkx->B);
    // printf("\n\n");
    // g2_print(vkx->R);
    // printf("\n\n");
    // g2_print(vkx->pi);

    pc_map(left,par->g,y);

	int flag =gt_cmp(left,right);

	if (flag == 2) return 2;
	return 0;
}
