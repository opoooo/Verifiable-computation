#include "lhep.h"

int lhep_new(lhe_par *par)
{	
	//get the order of the ecc groups
 	g1_get_ord(par->q);
 	bn2fmpz(par->qf,par->q);
	//define the n in the polynomial f(x)=x^n+1
	bn_set_dig(par->n,64);
	bn2int(&par->ni,par->n);
	bn2fmpz(par->nf,par->n);
	//define the discrete gaussian parameter r
	bn_set_dig(par->r,(1<<12));
 	bn2fmpz(par->rf,par->r);
	//define the prime p that guides the actual homomorphic computation
	bn_read_str(par->p,"1050031",20,10);
  	bn2fmpz(par->pf,par->p);

        //define the context, F_q
 	fq_ctx_init(par->ctx,par->qf,1,"ctx"); 

	//define the modf polynomial
	fq_t u; 
	fq_init(u,par->ctx);
	fq_one(u,par->ctx);
	fq_poly_init(par->modf,par->ctx);
	fq_poly_set_coeff(par->modf,par->ni,u,par->ctx);
	fq_poly_set_coeff(par->modf,0,u,par->ctx);
    
	//define the modp polynomial
	fq_set_fmpz(u,par->pf,par->ctx);
	fq_poly_init(par->modp,par->ctx);
	fq_poly_set_coeff(par->modp,0,u,par->ctx);

	// get the generators of the ECC groups G_1, G_2, and G_T
	//basic operations in the ecc groups
	g1_new(par->g);
	g1_get_gen(par->g);
	g2_new(par->h);
	g2_get_gen(par->h); 
	gt_new(par->gt);
	pc_map(par->gt, par->g, par->h); 
	
	return 0;
 
}

int lhep_free(lhe_par *par){
	bn_clean(par->p);
	fmpz_clear(par->pf);
	bn_clean(par->q);
	fmpz_clear(par->qf);
	bn_clean(par->n);
	fmpz_clear(par->nf);
	bn_clean(par->r);
	fmpz_clear(par->rf);
	fq_ctx_clear(par->ctx);
	fq_poly_realloc(par->modf,0,par->ctx);
	fq_poly_realloc(par->modp,0,par->ctx);
	g1_free(par->g);
	g2_free(par->h);
	gt_free(par->gt);
	return 0;
	
	
}

