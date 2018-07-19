#include "lhe.h"

// the key generation algorithm of the linearly homomorphic encryption: LHE.KeyGen
int lhe_keygen(fq_poly_t s, lhe_par *par)
{  
  fq_poly_rand(s,par->ni,par->ctx,par->r);
  return 0;
}

// the encryption algorithm of the linearly homomorphic encryption: LHE.Enc
int lhe_enc(lhe_c *c, fq_poly_t s, fq_poly_t m,lhe_par *par)
{
   fq_poly_t a, e;
   fq_poly_rand(a,par->ni,par->ctx,par->q); 
   fq_poly_rand(e,par->ni,par->ctx,par->r);

   fq_poly_init(c->p0,par->ctx);
   fq_poly_init(c->p1,par->ctx);

  // set the value of p1
  fq_t u1,u2;
  fq_init(u1,par->ctx);
  fq_zero(u1,par->ctx);
  fq_init(u2,par->ctx);
  fq_one(u2,par->ctx);
  fq_sub(u1,u1,u2,par->ctx);

  // fq_print(u1,par->ctx);
  fq_poly_scalar_mul_fq(c->p1,a,u1,par->ctx); 

  // fq_poly_print(c->p1,par->ctx);

  // set the value of c0
  fq_poly_t v;
  fq_poly_init(v,par->ctx);
  // p0=p0+as
  fq_poly_mul(v,a,s,par->ctx);
  fq_poly_add(c->p0,c->p0,v,par->ctx);
  // p0=p0+pe
  fq_set_fmpz(u1,par->pf,par->ctx);
  // fmpz_print(par->pf);
  fq_poly_scalar_mul_fq(v,e,u1,par->ctx); 
  fq_poly_add(c->p0,c->p0,v,par->ctx);
  // p0=p0+m
  fq_poly_add(c->p0,c->p0,m,par->ctx);
  
  //reduction modulo x^n+1
  fq_poly_t Q;
  fq_poly_init(Q,par->ctx);
  fq_poly_divrem(Q,c->p0,c->p0,par->modf,par->ctx); 

  // release the memory
  fq_poly_clear(a,par->ctx);
  fq_poly_clear(e,par->ctx);
  fq_clear(u1,par->ctx);
  fq_clear(u2,par->ctx);
  fq_poly_clear(v,par->ctx);
  fq_poly_clear(Q,par->ctx);

  return 0;   

}

// the decryption function of the linearly homomorphic encryption scheme: LEH.Dec
int lhe_dec(fq_poly_t m1, fq_poly_t s, lhe_c *c, lhe_par *par)
{
    fq_poly_init(m1,par->ctx);
    fq_poly_mul(m1,s,c->p1,par->ctx);
    fq_poly_add(m1,m1,c->p0,par->ctx);

    // reduction modulo x^n+1
    fq_poly_t Q;
    fq_poly_init(Q,par->ctx);
    fq_poly_divrem(Q,m1,m1,par->modf,par->ctx); 
   
    // release the memory
    fq_poly_clear(Q,par->ctx);
    
    return 0;   
}


// the eval function of the linearly homomorphic encryption scheme: LHE.Eval
int lhe_eval(lhe_c *c, lhe_c **C, fq_t *L, int d, lhe_par *par)
{
   fq_poly_init(c->p0,par->ctx);
   fq_poly_init(c->p1,par->ctx);
  
   fq_poly_t z;
   fq_poly_init(z,par->ctx);

   fq_poly_t Q;
   fq_poly_init(Q,par->ctx);

   for (int i=0; i<d; i++)
   {
       fq_poly_scalar_mul_fq(z,C[i]->p0,L[i],par->ctx);
       fq_poly_add(c->p0,c->p0,z,par->ctx);
       fq_poly_divrem(Q,c->p0,c->p0,par->modf,par->ctx);
       fq_poly_scalar_mul_fq(z,C[i]->p1,L[i],par->ctx);
       fq_poly_add(c->p1,c->p1,z,par->ctx);
       fq_poly_divrem(Q,c->p1,c->p1,par->modf,par->ctx);
   }

   // release the memory
   fq_poly_clear(z,par->ctx);
   fq_poly_clear(Q,par->ctx);
   return 0;

}