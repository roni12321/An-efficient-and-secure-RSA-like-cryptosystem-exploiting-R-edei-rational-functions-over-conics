
#include <stdlib.h> //for srand,rand
#include <time.h>   //for time
#include <gmp.h>    //for big numbers
#include <gmpxx.h>
#include <cstdlib>
#include <ctime>
#include <iostream> // I/O
#include <chrono>   //for measure time
#include <iomanip>

using namespace std::chrono;
using namespace std;

class Point{ /* class for point values */
public:
    mpz_t x;
    mpz_t y;
};

/* Function declaration */
void generatePublicKey(mpz_t e,mpz_t phi,gmp_randstate_t state);/* Define 1 < e < lcm(p+1,q+1) */
void nonQuadraticResidue(Point point,mpz_t N,mpz_t D); /* Define the parameter D */
void parametrization(Point point,mpz_t N, mpz_t M); /*parametrization from point with 2 values into one value */
void redeiSummation(mpz_t D,mpz_t e,mpz_t M,mpz_t N,mpz_t C); /* Encrypt function used by Rèdei rational function */
void inverseParametrization(mpz_t M,mpz_t D,mpz_t N,Point * DecryptedPoint);/* Inverse parametrization of the Decrypted message */

int main(int argc, const char * argv[]) {
    auto start = high_resolution_clock::now(); // start measure time for all program
    srand((unsigned int)time(NULL)); //cast for unsigned int for not loose values
    unsigned long long int seed;
    seed = (unsigned long long int)rand() % 100 + 1; /* Define the limit of the random number */
    gmp_randstate_t state;  /* Variable state for gmp_randinit, must be initialized */
    /* Create an object of Point and init the variables x,y*/
    Point newPoint; /* This point for using before decryption*/
    mpz_init(newPoint.x);
    mpz_init(newPoint.y);
    Point DecryptedPoint; /* The point with values that decrypted*/
    mpz_init(DecryptedPoint.x);
    mpz_init(DecryptedPoint.y);
    /* Define the parameters of the RSA algorithm */
    mp_bitcnt_t n;
    mpz_t p,q,N,e,d,D,M,C,pPlus1,qPlus1,lcmOfPplus1andQplus1,DecryptedMessage;
    mpz_init(p);
    mpz_init(q);
    mpz_init(N);
    mpz_init(pPlus1);
    mpz_init(qPlus1);
    mpz_init(lcmOfPplus1andQplus1);
    mpz_init(e);
    mpz_init(d);
    mpz_init(D);
    mpz_init(M);
    mpz_init(C);
    mpz_init(DecryptedMessage);
    n = 8; /* Number of bits for the number in range 0 to (2^n)-1, inclusive */
    gmp_randinit_default(state);    /* Initialize state with a default algorithm */
    gmp_randseed_ui ( state, seed ); /* Set an initial seed value into state */
    /* The size of a seed determines how many different sequences of random numbers that it’s possible to generate */
    while (true)
    {
        mpz_urandomb(p, state, n);  // Generate a uniformly distributed random integer
        if (mpz_probab_prime_p(p, 50) == 2) // Determine whether p is prime. Return 2 if p is definitely prime
            break;
    }
    while (true)
    {
        mpz_urandomb(q, state, n);  // Generate a uniformly distributed random integer
        if (mpz_probab_prime_p(q, 50) == 2) // Determine whether q is prime. Return 2 if q is definitely prime
            break;
    }
    mpz_mul(N, p, q);    /* Set N to p times q */
    gmp_printf("p =  %Z d\n", p);
    gmp_printf("q =  %Z d\n", q);
    gmp_printf("N =  %Z d\n", N);
    /* Calculate p+1 and q+1 */
    mpz_add_ui(pPlus1,p,1); /* Calculate p+1 and put the result into 'pPlus1'  */
    mpz_add_ui(qPlus1,q,1); /* Calculate q+1 and put the result into 'qPlus1'  */
    // now we should choose e such that gcd(e,lcm(p+1,q+1))=1
    mpz_lcm(lcmOfPplus1andQplus1, pPlus1, qPlus1); /* Calculate lcm(p+1,q+1) */
    generatePublicKey(e,lcmOfPplus1andQplus1,state); /* Calculate e by using the function generatePublicKey */
    gmp_printf("e = %Z d\n",e);
    if(mpz_invert(d, e, lcmOfPplus1andQplus1)==0){ /* choosing d such that it satisfies d*e = 1 mod lcm(p+1,q+1) */
        cout <<"There NO inverse parameter for e\n";
        return 0;
    }
    gmp_printf("d = %Z d\n",d);
    /* define the plain values */
    mpz_set_ui(newPoint.x, 83); /* Set the first part of message */
    mpz_set_ui(newPoint.y, 135); /* Set the second part of message */
    nonQuadraticResidue(newPoint,N,D); /* Calculate quadratic non-residue D */
    parametrization(newPoint,N,M); /* make parametrization for 2 values of point into one value */
    /* Encryption procces start */
    cout << "Encryption time is : ";
    redeiSummation(D,e,M,N,C); /* Encrypt by Rèdei rational function */
    /* Encryption procces end */
    gmp_printf("C =  %Z d\n", C);
    /* Decryption procces start */
     cout << "Decryption time is : ";
    redeiSummation(D,d,C,N,DecryptedMessage);
    /* Decryption procces end */
    gmp_printf("DecryptedMessage =  %Z d\n", DecryptedMessage);
    /* Inverse parametrization of the Decrypted message */
    inverseParametrization(M,D,N,&DecryptedPoint);
    gmp_printf("The Decrypted point value is:\nMx =  %Z d\nMy =  %Z d\n", DecryptedPoint.x,DecryptedPoint.y);
    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(N);
    mpz_clear(e);
    mpz_clear(d);
    mpz_clear(D);
    mpz_clear(M);
    mpz_clear(C);
    mpz_clear(pPlus1);
    mpz_clear(qPlus1);
    mpz_clear(lcmOfPplus1andQplus1);
    mpz_clear(DecryptedMessage);
    mpz_clear(newPoint.x);
    mpz_clear(newPoint.y);
    mpz_clear(DecryptedPoint.x);
    mpz_clear(DecryptedPoint.y);
    auto stop = high_resolution_clock::now(); // stop measure time
    auto duration = duration_cast<nanoseconds>(stop - start); // time measurement end - start
    cout <<"Total time for Fast-RSA : "<< fixed << setprecision(6) << (double)duration.count() /(1000000000)<<"\n";
    return 0;
}


void generatePublicKey(mpz_t e,mpz_t lcm,gmp_randstate_t state){/* Define 1 < e < lcm(p+1,q+1) */
    mpz_t one,zero,gcdResult; /* Declaration of function local parameters */
    /* Initialize the parameters */
    mpz_init(one);
    mpz_init(gcdResult);
    mpz_init(zero);/* Init zero to be 0 */
    mpz_add_ui(one, one, 1); /* define the parameter one to be with value 1 by adding */
    do{
        mpz_urandomm(e,state,lcm); /* Generate a uniformly distributed random integer */
        mpz_gcd (gcdResult, e, lcm); /* calculate the gcd(e,lcm(p+1,q+1)) and put the result into "gcdResult" */
    } /* The while checks if:
       * (1) "e" is not equals to zero.
       * (2) "e" is not equals to one.
       * (3) Gcd(e,lcm(p+1,q+1)) is not equals to zero.
       */
    while((mpz_cmp(e,zero) == 0)||(mpz_cmp(e,one) == 0)||(mpz_cmp(gcdResult,one) != 0));
    /* Clear the parameters that used in the function */
    mpz_clear(zero);
    mpz_clear(one);
}
void nonQuadraticResidue(Point point,mpz_t N,mpz_t D){
    mpz_t nom,denom;
    mpz_init(nom);
    mpz_init(denom);
    mpz_mul(nom, point.x, point.x);/* Calculate Mx^2 */
    mpz_sub_ui(nom, nom, 1); /* Calculate Mx^2 -1 */
    mpz_mul(denom, point.y, point.y);/* Calculate My^2 */
    if(mpz_invert(denom, denom, N)==0){ /* Calculate the inverse of My^2 mod N */
        cout <<"There NO inverse parameter for e\n";
        exit(1);
    }
    mpz_mul(D, nom, denom); /* D = nom * denom */
    mpz_mod(D,D,N); /* D = (nom * denom) mod N */
    mpz_clear(nom);
    mpz_clear(denom);
}
void parametrization(Point point,mpz_t N, mpz_t M){/*parametrization from point with 2 values into one value */
    mpz_t nom,denom;
    mpz_init(nom);
    mpz_init(denom);
    mpz_add_ui(nom, point.x, 1); /* Calculate Mx + 1 */
    if(mpz_invert(denom, point.y, N)==0){ /* Calculate the inverse of My mod N */
        cout <<"There NO inverse parameter for e\n";
        exit(1);
    }
    mpz_mul(M, nom, denom); /* D = nom * denom */
    mpz_mod(M,M,N); /* D = (nom * denom) mod N */
    mpz_clear(nom);
    mpz_clear(denom);
}
void redeiSummation(mpz_t D,mpz_t e,mpz_t M,mpz_t N,mpz_t C){/* Encrypt function used by Rèdei rational function */
    auto start2 = high_resolution_clock::now(); // start measure time for decryption
    unsigned long int k=0;
    mpz_t resultD,resultZ_An,index,sigmaLimit,e_fac,An,Bn,factorialResultAn1,factorialResultAn2,twiceIndex,nMinus2K,denom_nCk,cNk,cNkBn,term,twiceIndexPlus1,factorialResultBn1,nMinus2KMinus1,factorialResultBn2,resultZ_Bn,localD,localZ;
    mpz_init(resultD);
    mpz_init(resultZ_An);
    mpz_init(index);
    mpz_set_ui(index, 0); /* index = 0 */
    mpz_init(sigmaLimit);
    mpz_init(e_fac);
    mpz_init(An);
    mpz_init(Bn);
    mpz_init(factorialResultAn1); /* (2k)! */
    mpz_init(factorialResultAn2); /* (n - 2k)! */
    mpz_init(twiceIndex);
    mpz_init(nMinus2K);
    mpz_init(denom_nCk);
    mpz_init(cNk);
    mpz_init(cNkBn);
    mpz_init(term);
    mpz_init(twiceIndexPlus1);
    mpz_init(factorialResultBn1);
    mpz_init(nMinus2KMinus1);
    mpz_init(factorialResultBn2);
    mpz_init(resultZ_Bn);
    mpz_init(localD);
    mpz_init(localZ);
    mpz_set_ui(localD, 1);//D^0
    mpz_set_ui(twiceIndex,0);
    mpz_cdiv_q_ui(sigmaLimit, e, 2); /* Calculate the high limit of the sigma */
    /* put into D in case index 0*/
    mpz_set(resultD,localD);
    /* Z^(n - 2k) of An, z is the message M */
    mpz_powm(resultZ_An, M, nMinus2K, N);
    /* Z^(n - 2k - 1) of Bn */
    mpz_cdiv_q(resultZ_Bn, resultZ_An, M);
    while(mpz_cmp(index, sigmaLimit) != 0){
        /* An & Bn calculation */
        /* (2k)! (for An for An step)  */
        //        mpz_mul_ui(twiceIndex, index, 2); /* twiceIndex = 2 * index */
        //    k = (unsigned long int)twiceIndex;
        mpz_bin_ui(factorialResultAn1, e, k);
        k++;
        /* (2k + 1)! (for An for Bn step) */
        //        mpz_add_ui(twiceIndexPlus1,twiceIndex,1);
        // k = (unsigned long int)twiceIndexPlus1;
        mpz_bin_ui(factorialResultBn1, e, k);
        k++;
        if(mpz_get_ui(index)!=0){
           /* D^k -> the same result for An and Bn */
           mpz_mul(resultD, resultD, D);
           mpz_cdiv_q(resultZ_An, resultZ_Bn, M);
           mpz_cdiv_q(resultZ_Bn, resultZ_An, M);
        }
        /* multiply all factors for An */
        mpz_mul(term, resultZ_An, resultD);/* term = resultZ_An * resultD */
        mpz_mul(term, term, factorialResultAn1);/* term *= factorialResultAn1 */
        mpz_mod(term, term, N); /* term = term % N */
        /* An += term */
        mpz_add(An, An, term);
        /* multiply all factors for Bn */
        mpz_mul(term, resultZ_Bn, resultD);/* term = resultZ_Bn * resultD */
        mpz_mul(term, term, factorialResultBn1);/* term *= factorialResultBn1 */
        mpz_mod(term, term, N); /* term = term % N */
        /* Bn += term */
        mpz_add(Bn, Bn, term);
        mpz_add_ui(index, index, 1); /* index++ */
    }
    mpz_invert(Bn, Bn, N); /* inverse of Bn mod N */
    mpz_mul(C, An, Bn);/* c= An * Bn */
    mpz_mod(C, C, N);
    /* clear the parameters that used in the function */
    mpz_clear(resultD);
    mpz_clear(resultZ_An);
    mpz_clear(index);
    mpz_clear(sigmaLimit);
    mpz_clear(e_fac);
    mpz_clear(An);
    mpz_clear(Bn);
    mpz_clear(factorialResultAn1);
    mpz_clear(factorialResultAn2);
    mpz_clear(twiceIndex);
    mpz_clear(nMinus2K);
    mpz_clear(denom_nCk);
    mpz_clear(cNk);
    mpz_clear(cNkBn);
    mpz_clear(term);
    mpz_init(twiceIndexPlus1);
    mpz_clear(factorialResultBn1);
    mpz_clear(nMinus2KMinus1);
    mpz_clear(factorialResultBn2);
    mpz_clear(resultZ_Bn);
    auto stop2 = high_resolution_clock::now(); // stop measure time for encryption
    auto duration2 = duration_cast<nanoseconds>(stop2 - start2); // time measurement end - start
    cout << fixed << setprecision(6) << (double)duration2.count() / (1000000000)<<"\n";
}

void inverseParametrization(mpz_t M,mpz_t D,mpz_t N,Point * DecryptedPoint){
    mpz_t nom,denom;
    mpz_init(nom);
    mpz_init(denom);
    /*calculation for Mx - the x value in the point */
    /* nominator calculation*/
    mpz_powm_ui(nom, M, 2, N);/*     mpz_mul(nom, M, M) */
    mpz_add(nom, nom, D);
    /* denominator calculation */
    mpz_powm_ui(denom, M, 2, N);/*     mpz_mul(denom, M, M) */
    mpz_sub(denom, denom, D);
    /* calculation of the inverse of demominator */
    mpz_invert(denom, denom, N);
    /* calculation of the point mod N */
    mpz_mul(nom, nom, denom);
    mpz_mod(DecryptedPoint->x, nom, N);
    /* End of calculation for Mx */
    /*calculation for My - the y value in the point */
    /* nominator calculation*/
    mpz_mul_ui(nom, M, 2);
    /* The denominator is the same to the x value calculation */
    /* calculation of the point mod N */
    mpz_mul(nom, nom, denom);
    mpz_mod(DecryptedPoint->y, nom, N);
    mpz_clear(nom);
    mpz_clear(denom);
}



