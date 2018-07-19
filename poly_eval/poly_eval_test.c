#include "poly_eval.h"

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
 
        // rand polynomial algorithm
		    fq_poly_t A;
		    fq_poly_init(A,par->ctx);
        fq_poly_rand(A,64,par->ctx,par->p);

   
        // print 


              
        // test the vc_keygen algorithm
        pk_a *pka=malloc(sizeof(*pka));
        ek_a *eka=malloc(sizeof(*eka));
        vc_keygen(pka,eka,A,par);


        // print

  

        // test the vc_pgen algorithm
        vk_x *vkx=malloc(sizeof(*vkx));
        // generate a random message vector of d elements
        bn_t x;
        bn_new(x);
        bn_rand_mod(x,par->p);
        
        // compute the problem instance
        vc_pgen(pka,vkx,x,par);

        // print 


        // test the vc_comp algorithm
        g2_t y;
        g2_new(y);
        g2_set_infty(y);
        vc_comp(y,vkx,x,eka,par);
        // print
        
		// gt_print(y);



       // test the vc_vrfy algorithm     
       int flag=vc_vrfy(y,vkx,par);
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

