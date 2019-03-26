#ifndef RANDOMWALKS_RANDOM_HPP
#define RANDOMWALKS_RANDOM_HPP

#include <cmath>

using namespace std;

/*
 *  The algorithm is a combination of a Fibonacci sequence (with lags of 97 and
 *  33, and operation "subtraction plus one, modulo one") and an "arithmetic
 *  sequence" (using subtraction).
 *
 *  It passes ALL of the tests for random number generators and has a period of
 *  2^144, is completely portable (gives bit identical results on all machines
 *  with at least 24-bit mantissas in the floating point representation).
 *
 *  Marsaglia, G., Zaman, A., & Tsang, W. W. (1990). Toward a universal random
 *  number generator. Statistics & Probability Letters, 9(1), 35-39.
 *  doi:10.1016/0167-7152(90)90092-l
 *
 *  James, F. (1990). A review of pseudorandom number generators. Computer
 *  Physics Communications, 60(3), 329-344.
 *  doi:10.1016/0010-4655(90)90032-v
 */

/* Globals */

static double u[97], c, cd, cm;

static int i97, j97;

static bool test = false;

/*
 *  void RandomInitialise(int ij, int kl)
 *
 *  This is the initialization routine for the random number generator.
 *  The seed variables can have values between:
 *  0 <= IJ <= 31328 and 0 <= KL <= 30081
 *
 *  The random number sequences created by these two seeds are of sufficient
 *  length to complete an entire calculation with. For example, if several
 *  different groups are working on different parts of the same calculation,
 *  each group could be assigned its own IJ seed. This would leave each group
 *  with 30000 choices for the second seed. That is to say, this random
 *  number generator can create 900 million different subsequences -- with
 *  each subsequence having a length of approximately 64^30.
 *
 *  Use IJ = 1802 & KL = 9373 to test the random number generator. The
 *  subroutine should be used to generate 20000 random numbers. Then, display
 *  the next six random numbers generated multiplied by 4096*4096. If the
 *  random number generator is working properly, the random numbers should be:
 *  6533892.0   14220222.0   7275067.0   6172232.0   8354498.0   64633180.0
 */

void RandomInitialise(int ij,int kl) {
    double s,t;
    int i,j,k,l,m;

    /*
         Handle the seed range errors
         First random number seed must be between 0 and 31328
         Second seed must have a value between 0 and 30081
     */

    if (ij < 0 || ij > 31328 || kl < 0 || kl > 30081) {
        ij = 1802;
        kl = 9373;
    }

    i = (ij / 177) % 177 + 2;
    j = (ij % 177)       + 2;
    k = (kl / 169) % 178 + 1;
    l = (kl % 169);

    for (int ii = 0; ii < 97; ii++) {
        s = 0.0;
        t = 0.5;

        for (int jj = 0; jj < 24; jj++) {
            m = (((i * j) % 179) * k) % 179;
            i = j;
            j = k;
            k = m;
            l = (53 * l + 1) % 169;

            if (((l * m % 64)) >= 32) {
                s += t;
            }
            t *= 0.5;
        }

        u[ii] = s;
    }

    c    = 362436.0 / 16777216.0;
    cd   = 7654321.0 / 16777216.0;
    cm   = 16777213.0 / 16777216.0;
    i97  = 97;
    j97  = 33;
    test = true;
}

/*
 *  double RandomUniform()
 *
 *  This is the call function for the random number generator. It returns a
 *  double between 0 and 1.
 */

double RandomUniform(void) {
    double uni;

    /* Make sure the initialisation routine has been called */
    if (!test) {
        RandomInitialise(1802, 9373);
    }

    uni = u[i97 - 1] - u[j97 - 1];

    if (uni <= 0.0) {
        uni++;
    }

    u[i97 - 1] = uni;

    i97--;

    if (i97 == 0) {
        i97 = 97;
    }

    j97--;

    if (j97 == 0) {
        j97 = 97;
    }

    c -= cd;

    if (c < 0.0) {
        c += cm;
    }

    uni -= c;

    if (uni < 0.0) {
        uni++;
    }

    return uni;
}

/*
 *  double RandomGaussian(double mean,double stddev)
 *
 *  This function returns a normally distributed pseudo-random number with a
 *  given mean and standard deviation. Calls are made to a function subprogram
 *  which must return independent random numbers uniform in the interval (0,1).
 *  The algorithm uses the ratio of uniforms method of A.J. Kinderman and J.F.
 *  Monahan augmented with quadratic bounding curves.
 *
 *  Leva, J. L. (1992). Algorithm 712; a normal random number generator. ACM
 *  Transactions on Mathematical Software, 18(4), 454-455.
 *  doi:10.1145/138351.138367
 */

double RandomGaussian(double mean,double stddev) {
    double  q, u, v, x, y;

    /*
         Generate P = (u,v) uniform in rect. enclosing acceptance region
         Make sure that any random numbers <= 0 are rejected, since
         gaussian() requires uniforms > 0, but RandomUniform() delivers >= 0.
     */

    do {
        u = RandomUniform();
        v = RandomUniform();

        if (u <= 0.0 || v <= 0.0) {
            u = 1.0;
            v = 1.0;
        }

        v = 1.7156 * (v - 0.5);

        /*  Evaluate the quadratic form */
        x = u - 0.449871;
        y = fabs(v) + 0.386595;
        q = x * x + y * (0.19600 * y - 0.25472 * x);

        /* Accept P if inside inner ellipse */
        if (q < 0.27597) {
            break;
        }

        /*  Reject P if outside outer ellipse, or outside acceptance region */
    } while ((q > 0.27846) || (v * v > -4.0 * log(u) * u * u));

    /*  Return ratio of P's coordinates as the normal deviate */
    return (mean + stddev * v / u);
}

/*
 *  int RandomInt(int lower,int upper)
 *
 *  This function returns a random integer between the upper and lower bounds,
 *  inclusive.
 */

int RandomInt(int lower,int upper) {
    return((int)(RandomUniform() * (upper - lower + 1)) + lower);
}

/*
 *  double RandomIDouble(double lower,double upper)
 *
 *  This function returns a random double between the upper and lower bounds.
 */

double RandomDouble(double lower,double upper) {
    return((upper - lower) * RandomUniform() + lower);
}

#endif //RANDOMWALKS_RANDOM_HPP
