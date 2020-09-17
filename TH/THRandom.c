#include "THGeneral.h"
#include "THRandom.h"

#ifndef _WIN32
#include <fcntl.h>
#include <unistd.h>
#endif

/* Code for the Mersenne Twister random generator.... */
#define n _MERSENNE_STATE_N
#define m _MERSENNE_STATE_M

/* Creates (unseeded) new generator*/
static THGenerator* THGenerator_newUnseeded()
{
  THGenerator *self = THAlloc(sizeof(THGenerator));
  memset(self, 0, sizeof(THGenerator));
  self->left = 1;
  self->seeded = 0;
  self->normal_is_valid = 0;
  return self;
}

/* Creates new generator and makes sure it is seeded*/
THGenerator* THGenerator_new()
{
  THGenerator *self = THGenerator_newUnseeded();
  THRandom_seed(self);
  return self;
}

THGenerator* THGenerator_copy(THGenerator *self, THGenerator *from)
{
    memcpy(self, from, sizeof(THGenerator));
    return self;
}

void THGenerator_free(THGenerator *self)
{
  THFree(self);
}

int THGenerator_isValid(THGenerator *_generator)
{
  if ((_generator->seeded == 1) &&
    (_generator->left > 0 && _generator->left <= n) && (_generator->next <= n))
    return 1;

  return 0;
}

#ifndef _WIN32
static unsigned long readURandomLong()
{
  int randDev = open("/dev/urandom", O_RDONLY);
  unsigned long randValue;
  if (randDev < 0) {
    THError("Unable to open /dev/urandom");
  }
  ssize_t readBytes = read(randDev, &randValue, sizeof(randValue));
  if (readBytes < sizeof(randValue)) {
    THError("Unable to read from /dev/urandom");
  }
  close(randDev);
  return randValue;
}
#endif // _WIN32

unsigned long THRandom_seed(THGenerator *_generator)
{
#ifdef _WIN32
  unsigned long s = (unsigned long)time(0);
#else
  unsigned long s = readURandomLong();
#endif
  THRandom_manualSeed(_generator, s);
  return s;
}

/* The next 4 methods are taken from http:www.math.keio.ac.jpmatumotoemt.html
   Here is the copyright:
   Some minor modifications have been made to adapt to "my" C... */

/*
   A C-program for MT19937, with initialization improved 2002/2/10.
   Coded by Takuji Nishimura and Makoto Matsumoto.
   This is a faster version by taking Shawn Cokus's optimization,
   Matthe Bellew's simplification, Isaku Wada's double version.

   Before using, initialize the state by using init_genrand(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.keio.ac.jp/matumoto/emt.html
   email: matumoto@math.keio.ac.jp
*/

/* Macros for the Mersenne Twister random generator... */
/* Period parameters */
/* #define n 624 */
/* #define m 397 */
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UMASK 0x80000000UL /* most significant w-r bits */
#define LMASK 0x7fffffffUL /* least significant r bits */
#define MIXBITS(u,v) ( ((u) & UMASK) | ((v) & LMASK) )
#define TWIST(u,v) ((MIXBITS(u,v) >> 1) ^ ((v)&1UL ? MATRIX_A : 0UL))
/*********************************************************** That's it. */

void THRandom_manualSeed(THGenerator *_generator, unsigned long the_seed_)
{
  int j;

  /* This ensures reseeding resets all of the state (i.e. state for Gaussian numbers) */
  THGenerator *blank = THGenerator_newUnseeded();
  THGenerator_copy(_generator, blank);
  THGenerator_free(blank);

  _generator->the_initial_seed = the_seed_;
  _generator->state[0] = _generator->the_initial_seed & 0xffffffffUL;
  for(j = 1; j < n; j++)
  {
    _generator->state[j] = (1812433253UL * (_generator->state[j-1] ^ (_generator->state[j-1] >> 30)) + j);
    /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
    /* In the previous versions, mSBs of the seed affect   */
    /* only mSBs of the array state[].                        */
    /* 2002/01/09 modified by makoto matsumoto             */
    _generator->state[j] &= 0xffffffffUL;  /* for >32 bit machines */
  }
  _generator->left = 1;
  _generator->seeded = 1;
}

unsigned long THRandom_initialSeed(THGenerator *_generator)
{
  return _generator->the_initial_seed;
}

void THRandom_nextState(THGenerator *_generator)
{
  unsigned long *p = _generator->state;
  int j;

  _generator->left = n;
  _generator->next = 0;

  for(j = n-m+1; --j; p++)
    *p = p[m] ^ TWIST(p[0], p[1]);

  for(j = m; --j; p++)
    *p = p[m-n] ^ TWIST(p[0], p[1]);

  *p = p[m-n] ^ TWIST(p[0], _generator->state[0]);
}

unsigned long THRandom_random(THGenerator *_generator)
{
  unsigned long y;

  if (--(_generator->left) == 0)
    THRandom_nextState(_generator);
  y = *(_generator->state + (_generator->next)++);

  /* Tempering */
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);

  return y;
}

/* generates a random number on [0,1)-double-interval */
static double __uniform__(THGenerator *_generator)
{
  /* divided by 2^32 */
  return (double)THRandom_random(_generator) * (1.0/4294967296.0);
}

/*********************************************************

 Thanks *a lot* Takuji Nishimura and Makoto Matsumoto!

 Now my own code...

*********************************************************/

double THRandom_uniform(THGenerator *_generator, double a, double b)
{
  return(__uniform__(_generator) * (b - a) + a);
}

double THRandom_normal(THGenerator *_generator, double mean, double stdv)
{
  THArgCheck(stdv > 0, 2, "standard deviation must be strictly positive");

  /* This is known as the Box-Muller method */
  if(!_generator->normal_is_valid)
  {
    _generator->normal_x = __uniform__(_generator);
    _generator->normal_y = __uniform__(_generator);
    _generator->normal_rho = sqrt(-2. * log(1.0-_generator->normal_y));
    _generator->normal_is_valid = 1;
  }
  else
    _generator->normal_is_valid = 0;

  if(_generator->normal_is_valid)
    return _generator->normal_rho*cos(2.*M_PI*_generator->normal_x)*stdv+mean;
  else
    return _generator->normal_rho*sin(2.*M_PI*_generator->normal_x)*stdv+mean;
}

double THRandom_exponential(THGenerator *_generator, double lambda)
{
  return(-1. / lambda * log(1-__uniform__(_generator)));
}

double THRandom_cauchy(THGenerator *_generator, double median, double sigma)
{
  return(median + sigma * tan(M_PI*(__uniform__(_generator)-0.5)));
}

/* Faut etre malade pour utiliser ca.
   M'enfin. */
double THRandom_logNormal(THGenerator *_generator, double mean, double stdv)
{
  THArgCheck(stdv > 0, 2, "standard deviation must be strictly positive");
  return(exp(THRandom_normal(_generator, mean, stdv)));
}

double THRandom_gamma(THGenerator *state, double shape)
{
  if (shape < 1.) {
    double u, v, x, y;
    while (1) {
      u = THRandom_uniform(state, 0, 1);
      v = THRandom_exponential(state, 1);
      if (u <= 1.0 - shape) {
        x = pow(u, 1./shape);
        if (x <= v) {
          return x;
        }
      }
      else {
        y = -log((1 - u)/shape);
        x = pow(1.0 - shape + shape*y, 1./shape);
        if (x <= (v + y)) {
          return x;
        }
      }
    }
  }
  else if (shape > 1.) {
    double d = shape - (1./3.);
    double c = 1./sqrt(9. * d);
    double u, v, x = 0;
    do {
      x = THRandom_normal(state, 0, 1);
      v = (1 + c * x) * (1 + c * x) * (1 + c * x);
      u = THRandom_uniform(state, 0, 1);
    } while (v <= 0. ||
             (((log(u) >= 0.5 * x * x + d * (1 - v + log(v)))) &&
              (u < 1.0 - 0.0331*(x*x)*(x*x))));
    return d * v;
  }
  else {
    return THRandom_exponential(state, 1);
  }
}

double THRandom_gamma2 (THGenerator *state, double a, double scale)
{
  const static double sqrt32 = 5.656854;
  const static double exp_m1 = 0.36787944117144232159;/* exp(-1) = 1/e */
  /* Coefficients q[k] - for q0 = sum(q[k]*a^(-k))
   * Coefficients a[k] - for q = q0+(t*t/2)*sum(a[k]*v^k)
   * Coefficients e[k] - for exp(q)-1 = sum(e[k]*q^k)
   */
  const static double q1 = 0.04166669;
  const static double q2 = 0.02083148;
  const static double q3 = 0.00801191;
  const static double q4 = 0.00144121;
  const static double q5 = -7.388e-5;
  const static double q6 = 2.4511e-4;
  const static double q7 = 2.424e-4;

  const static double a1 = 0.3333333;
  const static double a2 = -0.250003;
  const static double a3 = 0.2000062;
  const static double a4 = -0.1662921;
  const static double a5 = 0.1423657;
  const static double a6 = -0.1367177;
  const static double a7 = 0.1233795;

  /* State variables [FIXME for threading!] :*/
  static double aa = 0.;
  static double aaa = 0.;
  static double s, s2, d;    /* no. 1 (step 1) */
  static double q0, b, si, c;/* no. 2 (step 4) */

  double e, p, q, r, t, u, v, w, x, ret_val;

  if (isnan(a) || isnan(scale)) {
    THError("gamma2: shape and scale should be valid positive numbers.");
    return nan("");
  }
  if (a <= 0.0 || scale <= 0.0) {
    if(scale == 0. || a == 0.) return 0.;
    THError("gamma2: shape and scale should be positive values.");
    return nan("");
  }
  if(!isfinite(a) || !isfinite(scale)) return (1.0/0.0);

  if (a < 1.) { /* GS algorithm for parameters a < 1 */
    e = 1.0 + exp_m1 * a;
    for(;;) {
      p = e * THRandom_uniform(state, 0, 1);
      if (p >= 1.0) {
        x = -log((e - p) / a);
        if (THRandom_exponential(state, 1) >= (1.0 - a) * log(x))
          break;
      } else {
        x = exp(log(p) / a);
        if (THRandom_exponential(state, 1) >= x)
          break;
      }
    }
    return scale * x;
  }

  /* --- a >= 1 : GD algorithm --- */

  /* Step 1: Recalculations of s2, s, d if a has changed */
  if (a != aa) {
    aa = a;
    s2 = a - 0.5;
    s = sqrt(s2);
    d = sqrt32 - s * 12.0;
  }
  /* Step 2: t = standard normal deviate,
     x = (s,1/2) -normal deviate. */

  /* immediate acceptance (i) */
  t = THRandom_normal(state, 0, 1);
  x = s + 0.5 * t;
  ret_val = x * x;
  if (t >= 0.0)
    return scale * ret_val;

  /* Step 3: u = 0,1 - uniform sample. squeeze acceptance (s) */
  u = THRandom_uniform(state, 0, 1);
  if (d * u <= t * t * t)
    return scale * ret_val;

  /* Step 4: recalculations of q0, b, si, c if necessary */

  if (a != aaa) {
    aaa = a;
    r = 1.0 / a;
    q0 = ((((((q7 * r + q6) * r + q5) * r + q4) * r + q3) * r
           + q2) * r + q1) * r;

      /* Approximation depending on size of parameter a */
      /* The constants in the expressions for b, si and c */
      /* were established by numerical experiments */

      if (a <= 3.686) {
        b = 0.463 + s + 0.178 * s2;
        si = 1.235;
        c = 0.195 / s - 0.079 + 0.16 * s;
      } else if (a <= 13.022) {
        b = 1.654 + 0.0076 * s2;
        si = 1.68 / s + 0.275;
        c = 0.062 / s + 0.024;
      } else {
        b = 1.77;
        si = 0.75;
        c = 0.1515 / s;
      }
  }
  /* Step 5: no quotient test if x not positive */

  if (x > 0.0) {
    /* Step 6: calculation of v and quotient q */
    v = t / (s + s);
    if (fabs(v) <= 0.25)
      q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v
                                + a3) * v + a2) * v + a1) * v;
    else
      q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);


      /* Step 7: quotient acceptance (q) */
      if (log(1.0 - u) <= q)
        return scale * ret_val;
  }

  for(;;) {
    /* Step 8: e = standard exponential deviate
     *	u =  0,1 -uniform deviate
     *	t = (b,si)-double exponential (laplace) sample */
    e = THRandom_exponential(state, 1);
    u = THRandom_uniform(state, 0, 1);
    u = u + u - 1.0;
    if (u < 0.0)
      t = b - si * e;
    else
      t = b + si * e;
    /* Step	 9:  rejection if t < tau(1) = -0.71874483771719 */
    if (t >= -0.71874483771719) {
      /* Step 10:	 calculation of v and quotient q */
      v = t / (s + s);
      if (fabs(v) <= 0.25)
        q = q0 + 0.5 * t * t *
          ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v
            + a2) * v + a1) * v;
      else
        q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
        /* Step 11:	 hat acceptance (h) */
        /* (if q not positive go to step 8) */
        if (q > 0.0) {
          w = expm1(q);
          /*  ^^^^^ original code had approximation with rel.err < 2e-7 */
          /* if t is rejected sample again at step 8 */
          if (c * fabs(u) <= w * exp(e - 0.5 * t * t))
            break;
        }
    }
  } /* repeat .. until  `t' is accepted */
  x = s + 0.5 * t;
  return scale * x * x;
}

double THRandom_beta(THGenerator *_generator, double a, double b)
{
  double output;

  if (a <= 1 && b <= 1) {
    double U, V, X, Y;

    while (1) {
      U = THRandom_uniform(_generator, 0, 1);
      V = THRandom_uniform(_generator, 0, 1);
      X = pow(U, 1/a);
      Y = pow(V, 1/a);

      if ((X + Y) <= 1) {
        if ((X + Y) > 0) {
          output = X / (X + Y);
          break;
        }
        else {
          double logX = log(U)/a;
          double logY = log(V)/b;
          double logM = logX > logY ? logX : logY;
          logX -= logM;
          logY -= logM;

          output = exp(logX - log(exp(logX) + exp(logY)));
          break;
        }
      }
    }
  }
  else {
    double Ga = THRandom_gamma(_generator, a);
    double Gb = THRandom_gamma(_generator, b);
    output = Ga / (Ga + Gb);
  }

  return output;
}

int THRandom_geometric(THGenerator *_generator, double p)
{
  THArgCheck(p > 0 && p < 1, 1, "must be > 0 and < 1");
  return((int)(log(1-__uniform__(_generator)) / log(p)) + 1);
}

int THRandom_bernoulli(THGenerator *_generator, double p)
{
  THArgCheck(p >= 0 && p <= 1, 1, "must be >= 0 and <= 1");
  return(__uniform__(_generator) <= p);
}
