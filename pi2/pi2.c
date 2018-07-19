#include "pi2.h"

int bn_mul_modz(bn_t out, bn_t in1, bn_t in2, bn_t modulus)
{
   bn_mul_basic(out,in1,in2);
   bn_mod_basic(out,out,modulus);
   return 0;
}

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


// this function reduce the coefficients of a polynomial
int msg_modp(fq_poly_t m, lhe_par *par)
{
 
   long len=fq_poly_length(m,par->ctx);



   fq_t zq;
   fq_init(zq,par->ctx);

   fmpz_t zf;
   fmpz_init(zf);

   for (int i=0;i<len;i++)
   {

     fq_poly_get_coeff(zq,m,i,par->ctx);
     fmpz_set_str(zf,fq_get_str_pretty(zq,par->ctx), 10);
     fmpz_mod(zf,zf,par->pf);
     fq_set_fmpz(zq,zf,par->ctx);
     fq_poly_set_coeff(m,i,zq,par->ctx);
   }

    fq_clear(zq,par->ctx);
    fmpz_clear(zf);

   return 0;
 
}


// the vc_keygen algorithm: generating the evaluation key and public verification key of the scheme
int vc_keygen(vc_k *vck, fq_mat_t F, lhe_par *par)
{
   vck->m=fq_mat_nrows(F,par->ctx);
   vck->d=fq_mat_ncols(F,par->ctx);

   fq_mat_init(vck->F,vck->m,vck->d,par->ctx);
   fq_mat_set(vck->F,F,par->ctx);

   vck->G=malloc(sizeof(g1_t)*vck->m);
   vck->T=malloc(sizeof(g1_t)*vck->d);

   fq_mat_t R;
   fq_mat_init(R,1,vck->m,par->ctx);
   bn_t zb;
   bn_new(zb);
   fmpz_t zf;
   fmpz_init(zf);
   fq_t zq;
   fq_init(zq,par->ctx);



   for (int i=0;i<vck->m;i++)
   {
      bn_rand_mod(zb,par->q);
      bn2fmpz(zf,zb);
   
      fq_set_fmpz(zq,zf,par->ctx);
      fq_mat_entry_set(R,0,i,zq,par->ctx);
      g1_new(vck->G[i]);
      g1_mul(vck->G[i],par->g,zb);      
   }
   

   printf("\n\n this is the random vector R:\n\n");
   fq_mat_print_pretty(R,par->ctx);


   fq_mat_t RF;
   fq_mat_init(RF,1,vck->d,par->ctx);
   fq_mat_mul(RF,R,F,par->ctx);

   printf("\n\n this is the random vector RF:\n\n");
   fq_mat_print_pretty(RF,par->ctx);
   
    
   for (int j=0;j<vck->d;j++)
   {
      bn_read_str(zb,fq_get_str_pretty(fq_mat_entry(RF,0,j),par->ctx),strlen(fq_get_str_pretty(fq_mat_entry(RF,0,j),par->ctx)),10);
      g1_new(vck->T[j]);
      g1_mul(vck->T[j],par->g,zb);   
   } 


   // release the memory
   fq_mat_clear(R,par->ctx);
   bn_free(zb);
   fmpz_clear(zf);
   fq_clear(zq,par->ctx);   
   return 0;
}


// generate a random matrix over fq
int fq_mat_randz(fq_mat_t F, int m, int d, lhe_par *par)
{
   bn_t zb;
   bn_new(zb);
   fmpz_t zf;
   fmpz_init(zf);
   fq_t zq;
   fq_init(zq,par->ctx);

   for (int i=0;i<m;i++)
   {
     for (int j=0;j<d;j++)
     {
         bn_rand_mod(zb,par->p);
         bn2fmpz(zf,zb);
         fq_set_fmpz(zq,zf,par->ctx);
         fq_mat_entry_set(F,i,j,zq,par->ctx);        
     }
   }
   return 0;
}



// the vc_pgen algorithm
int vc_pgen(vc_p *vcp, vc_k * vck, bn_t *x, lhe_par *par)
{
   // choose secret key for the lhe scheme
   lhe_keygen(vcp->s,par);

   // encrypting the message vector x
   fmpz_t zf;
   fmpz_init(zf);
   fq_t zq;
   fq_init(zq,par->ctx);

   fq_poly_t zqp;
   fq_poly_init(zqp,par->ctx);


   vcp->C=malloc(sizeof(lhe_c*)*vck->d);



   for (int j=0;j<vck->d;j++)
   {
      bn2fmpz(zf,x[j]);
      fq_set_fmpz(zq,zf,par->ctx);
      fq_poly_set_fq(zqp,zq,par->ctx);
      vcp->C[j]=malloc(sizeof(lhe_c));
      lhe_enc(vcp->C[j],vcp->s,zqp,par);
   }
   



   // generate the hash keys
   hash_k *zhk;
   vcp->hk=malloc(sizeof(*zhk));
   hash_keygen(vcp->hk,par);

   g2_t hv;
   g2_set_infty(hv);
   
   gt_set_unity(vcp->pi);
   gt_t zt;
   gt_set_unity(zt);

   for (int j=0;j<vck->d;j++)
   {
      hash_H(hv,vcp->C[j],vcp->hk,0,par);
      pc_map(zt,vck->T[j],hv);
      gt_mul(vcp->pi,vcp->pi,zt);
   }

   return 0;
}





// the vc_comp algorithm
int vc_comp(lhe_c **nv, vc_k *vck, vc_p *vcp, lhe_par *par)
{
  fq_t *L=malloc(sizeof(fq_t)*vck->d);


  for (int j=0;j<vck->d;j++)
  {
    fq_init(L[j],par->ctx);
  }



  for (int i=0;i<vck->m;i++)
  {
     for (int j=0;j<vck->d;j++)
     {
       fq_set(L[j],fq_mat_entry(vck->F,i,j),par->ctx);
     
     }
     nv[i]=malloc(sizeof(lhe_c));
     lhe_eval(nv[i],vcp->C,L,vck->d,par);     
  }    
  return 0;
}

// the vc_vrfy algorithm
int vc_vrfy(fq_poly_t *y, vc_k *vck, vc_p *vcp, lhe_c **nv, lhe_par *par)
{
  int flag=-1;
  /*
  g2_t hv;
  g2_new(hv);
  g2_set_infty(hv);
  */
  fq_t zq;
  g1_t temp;
  g1_t right;
  g1_t left;
  bn_t temp1;
  bn_t temp2;

  for (int l=0;l<2*par->ni;l++)
  {
     g1_new(left);
     g1_set_infty(left);

     for (int i=0;i<vck->m;i++)
     {
        g1_new(temp);
        g1_set_infty(temp);
        fq_init(zq,par->ctx);
        if (l < par->ni) fq_poly_get_coeff(zq,nv[i]->p0,l,par->ctx);
        else fq_poly_get_coeff(zq,nv[i]->p1,l-par->ni,par->ctx);
        bn_new(temp1);
        fq2bn(temp1,zq,par->ctx);
        g1_mul(temp,vck->G[i],temp1);
        g1_add(left,left,temp);
     }

     g1_new(right);
     g1_set_infty(right);

     for (int j=0;j<vck->d;j++)
     {
        g1_new(temp);
        g1_set_infty(temp);
        fq_init(zq,par->ctx);
        if (l < par->ni) fq_poly_get_coeff(zq,vcp->C[j]->p0,l,par->ctx);
        else fq_poly_get_coeff(zq,vcp->C[j]->p1,l-par->ni,par->ctx);
        bn_new(temp2);
        fq2bn(temp2,zq,par->ctx);
        g1_mul(temp,vck->T[j],temp2);
        g1_add(right,right,temp);
     }
     g1_print(left);
     printf("\n\n");
     g1_print(right);
     printf("\n\n");
     flag=g1_cmp(left,right);
     if (flag == 2) return 2;
  }
  /*
  for (int i=0;i<vck->m;i++)
  {
     g2_set_infty(hv); // this looks completely strange, why?
     hash_H(hv,nv[i],vcp->hk,1,par);
     g2_norm(hv,hv);
     printf("\n\n this is hv of nv[%d]\n\n",i);
     g2_print(hv);
     printf("\n\n");
     pc_map(zgt,vck->G[i],hv);
     gt_mul(left,left,zgt);
  }

  printf("\n\n the verification tag left1.1 is: \n\n");
  gt_print(left);


  gt_set_unity(left);
  for (int i=0;i<vck->m;i++)
  {
     hash_H(hv,nv[i],vcp->hk,0,par);
     printf("\n\n this is hv of nv[%d]\n\n",i);
     g2_print(hv);
     printf("\n\n");
     pc_map(zgt,vck->G[i],hv);
     gt_mul(left,left,zgt);
  }

  printf("\n\n the verification tag left0 is: \n\n");
  gt_print(left);

  flag=gt_cmp(left,vcp->pi);
   
  
  if (flag ==2)
  {
     return 2;
  }
*/
  for (int i=0;i<vck->m;i++)
  {
      lhe_dec(y[i],vcp->s,nv[i],par);
  }
   
  return 0;
}
