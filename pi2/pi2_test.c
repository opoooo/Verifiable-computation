#include "pi2.h"


int main(int argc, char **argv)
{

// initialize the project /////////////////////////////////////////////////        
        if (md2_init())
        {
                printf("Testing FAILED\n");
                printf("Problem initializing the library\n");
               return 1;
        }


  lhe_par *par=malloc(sizeof(*par));
  lhep_new(par);
printf("\n--------------------------------------begin test-------------------------\n\n\n");
 
        // test the rand matrix algorithm
        int mat_m=100;
        int mat_d=100;
        fq_mat_t F;
        fq_mat_init(F,mat_m,mat_d,par->ctx);
        fq_mat_randz(F,mat_m,mat_d,par);
   
        // print 
        printf("\n\n the matrix F is:\n\n");
        fq_mat_print_pretty(F,par->ctx);
      
              
        // test the vc_keygen algorithm
        vc_k *vck=malloc(sizeof(*vck));
        vc_keygen(vck,F,par);
        // print
        printf("\n\n the integer vck->m=%d\n\n",vck->m);
        printf("\n\n the integer vck->d=%d\n\n",vck->d);
        printf("\n\n the matrix vck->F is: \n\n");
        fq_mat_print_pretty(vck->F,par->ctx);
        printf("\n\n the vck->G vector is:\n\n");
        for (int i=0;i<vck->m;i++)
        {
         g1_print(vck->G[i]);
         printf("\n\n");
        }
        printf("\n\n the vck->T vector is:\n\n");
        for (int j=0;j<vck->d;j++)
        {
         g1_print(vck->T[j]);
         printf("\n\n");
  }
        
       

        // test the vc_pgen algorithm
        vc_p *vcp=malloc(sizeof(*vcp));
        // generate a random message vector of d elements
        bn_t x[vck->d];
        for (int j=0;j<vck->d;j++)
          bn_rand_mod(x[j],par->p);
        // compute the problem instance
        vc_pgen(vcp,vck,x,par);
        // print
        printf("\n\n the encryption secret key vcp->s is\n\n");
        fq_poly_print_pretty(vcp->s,"x",par->ctx);


        printf("\n\n the input x is:\n\n");
        for (int j=0;j<vck->d;j++)
        { 
           printf("\n\n this is x[%d]:\n\n",j);
           bn_print(x[j]);
           printf("\n\n");
        }


        fq_poly_t m1;
        fq_poly_init(m1,par->ctx);

        printf("\n\n the encoded input is:\n\n");
        for (int j=0;j<vck->d;j++)
        { 
           printf("\n\n this is mu[%d]->p0:\n\n",j);
           fq_poly_print_pretty((vcp->C[j])->p0,"x",par->ctx);
           printf("\n\n");
           printf("\n\n this is mu[%d]->p1:\n\n",j);
           fq_poly_print_pretty((vcp->C[j])->p1,"x",par->ctx);
           printf("\n\n");
           printf("\n\n the decryption of vcp->C[%d] is:\n\n",j);
           lhe_dec(m1,vcp->s,vcp->C[j],par);
           fq_poly_print_pretty(m1,"x",par->ctx);
           printf("\n\n\n");
        }

        printf("\n\n the verification key is vcp->pi:\n\n");
        gt_print(vcp->pi);

       

        // test the vc_comp algorithm
        lhe_c **nv=malloc(sizeof(lhe_c)*vck->m);
        vc_comp(nv,vck,vcp,par);
        // print
        printf("\n\n the computation result is:\n\n");
        for (int i=0;i<vck->m;i++)
        { 
           printf("\n\n this is nv[%d]->p0:\n\n",i);
           fq_poly_print_pretty(nv[i]->p0,"x",par->ctx);
           printf("\n\n");

           printf("\n\n this is nv[%d]->p1:\n\n",i);
           fq_poly_print_pretty(nv[i]->p1,"x",par->ctx);
           printf("\n\n");

           printf("\n\n the decryption of nv[%d] is:\n\n",i);
           lhe_dec(m1,vcp->s,nv[i],par);
           fq_poly_print_pretty(m1,"x",par->ctx);
           printf("\n\n\n");
        }

       // test the vc_vrfy algorithm
       fq_poly_t y[vck->m];
       int flag=vc_vrfy(y,vck,vcp,nv,par);
       printf("\n\n the verification result is flag=%d\n\n",flag);
       

      
 
printf("\n\n\n--------------------------------------end  test------------------------\n\n\n");







    
    

 
        





 

/*
        //release the memory
  bn_free(order);
  fq_clear(u,par->ctx);
  g2_free(zg);
  bn_free(zt);
  fq_poly_clear(s,par->ctx);
  fq_poly_clear(m,par->ctx);
  fq_poly_clear(m1,par->ctx);
  fq_poly_clear(m2,par->ctx);
*/


        return 0;
}

