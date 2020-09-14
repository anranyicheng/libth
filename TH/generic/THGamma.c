#ifndef TH_GENERIC_FILE
#define TH_GENERIC_FILE "generic/THGamma.c"
#else

#if defined(TH_REAL_IS_DOUBLE)
static int
zeta_impl_series_dbl(double *a, double *b, double *s, const double x, const double machep) {
  int i = 0;
  while ((i < 9) || (*a <= 9.0)) {
    i += 1;
    *a += 1.0;
    *b = pow(*a, -x);
    *s += *b;
    if (fabs(*b / *s) < machep) return 1;
  }

  // Return whether we are done
  return 0;
}

static double
zeta_impl_dbl(double x, double q) {
  int i;
  double p, r, a, b, k, s, t, w;

  const double A[] = {
    12.0,
    -720.0,
    30240.0,
    -1209600.0,
    47900160.0,
    -1.8924375803183791606e9, /*1.307674368e12/691*/
    7.47242496e10,
    -2.950130727918164224e12,  /*1.067062284288e16/3617*/
    1.1646782814350067249e14,  /*5.109094217170944e18/43867*/
    -4.5979787224074726105e15, /*8.028576626982912e20/174611*/
    1.8152105401943546773e17,  /*1.5511210043330985984e23/854513*/
    -7.1661652561756670113e18  /*1.6938241367317436694528e27/236364091*/
  };

  const double maxnum = DBL_MAX;
  const double zero = 0.0, half = 0.5, one = 1.0;
  const double machep = 1e-15;

  if (x == one) return maxnum;

  if (x < one) {
    return zero;
  }

  if (q <= zero) {
    if (q == floor(q)) {
      return maxnum;
    }
    p = x;
    r = floor(p);
    if (p != r) return zero;
  }

  /* Permit negative q but continue sum until n+q > +9 .
   * This case should be handled by a reflection formula.
   * If q<0 and x is an integer, there is a relation to
   * the polygamma function.
   */
  s = pow(q, -x);
  a = q;
  b = zero;

  // Run the summation in a helper function that is specific to the doubleing
  // precision
  if (zeta_impl_series_dbl(&a, &b, &s, x, machep)) {
    return s;
  }

  w = a;
  s += b * w / (x - one);
  s -= half * b;
  a = one;
  k = zero;
  for (i = 0; i < 12; i++) {
    a *= x + k;
    b /= w;
    t = a * b / A[i];
    s = s + t;
    t = fabs(t / s);
    if (t < machep) return s;
    k += one;
    a *= x + k;
    b /= w;
    k += one;
  }
  return s;
};

static double
polynomial_evaluation_dbl(double x, const double *f, int n) {
  double result = 0.0;
  for (int i = 0; i < n; i++) {
    result *= x;
    result += f[i];
  }
  return result;
}

static double
digamma_impl_maybe_poly_dbl(const double s) {
  const double A[] = {8.33333333333333333333E-2, -2.10927960927960927961E-2,
    7.57575757575757575758E-3, -4.16666666666666666667E-3,
    3.96825396825396825397E-3, -8.33333333333333333333E-3,
    8.33333333333333333333E-2};

  double z;
  if (s < 1.0e17) {
    z = 1.0 / (s * s);
    return z * polynomial_evaluation_dbl(z, A, 7);
  } else
    return 0.0;
}

static double
digamma_impl_dbl(const double u) {
  double x = u;
  double p, q, nz, s, w, y;
  char negative;

  const double maxnum = FLT_MAX;
  const double m_pi = M_PI;

  negative = 0;
  nz = 0.0;

  const double zero = 0.0;
  const double one = 1.0;
  const double half = 0.5;

  if (x <= zero) {
    negative = one;
    q = x;
    p = floor(q);
    if (p == q) {
      return maxnum;
    }
    /* Remove the zeros of tan(m_pi x)
     * by subtracting the nearest integer from x
     */
    nz = q - p;
    if (nz != half) {
      if (nz > half) {
        p += one;
        nz = q - p;
      }
      nz = m_pi / tan(m_pi * nz);
    } else {
      nz = zero;
    }
    x = one - x;
  }

  /* use the recurrence psi(x+1) = psi(x) + 1/x. */
  s = x;
  w = zero;
  while (s < 10.0) {
    w += one / s;
    s += one;
  }

  y = digamma_impl_maybe_poly_dbl(s);

  y = log(s) - (half / s) - y - w;

  return (negative) ? y - nz : y;
}

static double
polygamma_impl_dbl(int n, double x) {
  if (n == 0) {
    return digamma_impl_dbl(x);
  }

  // dumb code to calculate factorials
  double factorial = 1.0;
  for (int i = 0; i < n; i++) {
    factorial *= (i + 1);
  }

  return pow(-1.0, n + 1) * factorial * zeta_impl_dbl(n + 1, x);
  /* return ((n + 1 % 2) ? -1.0 : 1.0) * factorial * zeta_impl_dbl(n + 1, x); */
}

static double
lbeta_impl_dbl(double a, double b) {
  a = fabs(a);
  b = fabs(b);
  return lgamma(a) + lgamma(b) - lgamma(a + b);
}

#else
static int
zeta_impl_series(float *a, float *b, float *s, const float x, const float machep) {
  int i = 0;
  while (i < 9) {
    i += 1;
    *a += 1.0f;
    *b = powf(*a, -x);
    *s += *b;
    if (fabsf(*b / *s) < machep) {
      return 1;
    }
  }

  // Return whether we are done
  return 0;
}

static float
zeta_impl(float x, float q) {
  int i;
  float p, r, a, b, k, s, t, w;

  const float A[] = {
    12.0,
    -720.0,
    30240.0,
    -1209600.0,
    47900160.0,
    -1.8924375803183791606e9, /*1.307674368e12/691*/
    7.47242496e10,
    -2.950130727918164224e12,  /*1.067062284288e16/3617*/
    1.1646782814350067249e14,  /*5.109094217170944e18/43867*/
    -4.5979787224074726105e15, /*8.028576626982912e20/174611*/
    1.8152105401943546773e17,  /*1.5511210043330985984e23/854513*/
    -7.1661652561756670113e18  /*1.6938241367317436694528e27/236364091*/
  };

  const float maxnum = FLT_MAX;
  const float zero = 0.0, half = 0.5, one = 1.0;
  const float machep = 1e-15;

  if (x == one) return maxnum;

  if (x < one) {
    return zero;
  }

  if (q <= zero) {
    if (q == floorf(q)) {
      return maxnum;
    }
    p = x;
    r = floorf(p);
    if (p != r) return zero;
  }

  /* Permit negative q but continue sum until n+q > +9 .
   * This case should be handled by a reflection formula.
   * If q<0 and x is an integer, there is a relation to
   * the polygamma function.
   */
  s = powf(q, -x);
  a = q;
  b = zero;

  // Run the summation in a helper function that is specific to the floating
  // precision
  if (zeta_impl_series(&a, &b, &s, x, machep)) {
    return s;
  }

  w = a;
  s += b * w / (x - one);
  s -= half * b;
  a = one;
  k = zero;
  for (i = 0; i < 12; i++) {
    a *= x + k;
    b /= w;
    t = a * b / A[i];
    s = s + t;
    t = fabsf(t / s);
    if (t < machep) return s;
    k += one;
    a *= x + k;
    b /= w;
    k += one;
  }
  return s;
};

static float
polynomial_evaluation(float x, const float *f, int n) {
  float result = 0.0;
  for (int i = 0; i < n; i++) {
    result *= x;
    result += f[i];
  }
  return result;
}

static float
digamma_impl_maybe_poly(const float s) {
  const float A[] = {-4.16666666666666666667E-3f, 3.96825396825396825397E-3f,
    -8.33333333333333333333E-3f, 8.33333333333333333333E-2f};
  float z;
  if (s < 1.0e8f) {
    z = 1.0f / (s * s);
    return z * polynomial_evaluation(z, A, 4);
  } else {
    return 0.0f;
  }
}

static float
digamma_impl(const float u) {
  float x = u;
  float p, q, nz, s, w, y;
  char negative;

  const float maxnum = FLT_MAX;
  const float m_pi = M_PI;

  negative = 0;
  nz = 0.0;

  const float zero = 0.0;
  const float one = 1.0;
  const float half = 0.5;

  if (x <= zero) {
    negative = one;
    q = x;
    p = floorf(q);
    if (p == q) {
      return maxnum;
    }
    /* Remove the zeros of tan(m_pi x)
     * by subtracting the nearest integer from x
     */
    nz = q - p;
    if (nz != half) {
      if (nz > half) {
        p += one;
        nz = q - p;
      }
      nz = m_pi / tanf(m_pi * nz);
    } else {
      nz = zero;
    }
    x = one - x;
  }

  /* use the recurrence psi(x+1) = psi(x) + 1/x. */
  s = x;
  w = zero;
  while (s < 10.0) {
    w += one / s;
    s += one;
  }

  y = digamma_impl_maybe_poly(s);

  y = logf(s) - (half / s) - y - w;

  return (negative) ? y - nz : y;
}

static float
polygamma_impl(int n, float x) {
  if (n == 0) {
    return digamma_impl(x);
  }

  // dumb code to calculate factorials
  float factorial = 1.0;
  for (int i = 0; i < n; i++) {
    factorial *= (i + 1);
  }

  return powf(-1.0, n + 1) * factorial * zeta_impl(n + 1, x);
}

static float
lbeta_impl(float a, float b) {
  a = fabsf(a);
  b = fabsf(b);
  return lgammaf(a) + lgammaf(b) - lgammaf(a + b);
}

#endif

void
THTensor_(polygamma)(long n, THTensor *input, THTensor *output) {
  real *inputData, *outputData;

  if (input->nDimension != 2) {
    THError("polygamma: input is supposed to be 2D.");
    return;
  }

  inputData = THTensor_(data)(input);
  outputData = THTensor_(data)(output);

  long inputStrideHeight = input->stride[0];
  long inputStrideWidth = input->stride[1];

  long outputStrideHeight = output->stride[0];
  long outputStrideWidth = output->stride[1];

  long height = input->size[0];
  long width = input->size[1];

  for (long x = 0; x < width; x++) {
    for (long y = 0; y < height; y++) {
      const long outAddress = x * outputStrideWidth + y * outputStrideHeight;
      const long inAddress = x * inputStrideWidth + y * inputStrideHeight;
#if defined(TH_REAL_IS_DOUBLE)
      outputData[outAddress] = polygamma_impl_dbl(n, inputData[inAddress]);
#else
      outputData[outAddress] = polygamma_impl(n, inputData[inAddress]);
#endif
    }
  }

  return;
}

void
THTensor_(lbeta)(THTensor *a, THTensor *b, THTensor *output) {
  real *aData, *bData, *outputData;

  if (a->nDimension != 2 || b->nDimension != 2) {
    THError("lbeta: inputs are supposed to be 2D.");
    return;
  }

  aData = THTensor_(data)(a);
  bData = THTensor_(data)(b);
  outputData = THTensor_(data)(output);

  long aStrideWidth = a->stride[0];
  long aStrideHeight = a->stride[1];

  long bStrideWidth = b->stride[0];
  long bStrideHeight = b->stride[1];

  long outputStrideWidth = output->stride[0];
  long outputStrideHeight = output->stride[1];

  long height = a->size[0];
  long width = a->size[1];

  for (int x = 0; x < width; x++) {
    for (int y = 0; y < height; y++) {
      const long outAddress = x * outputStrideWidth + y * outputStrideHeight;
      const long aAddress = x * aStrideWidth + y * aStrideHeight;
      const long bAddress = x * bStrideWidth + y * bStrideHeight;
      outputData[outAddress] =
#if defined(TH_REAL_IS_DOUBLE)
        lbeta_impl_dbl(aData[aAddress], bData[bAddress]);
#else
        lbeta_impl(aData[aAddress], bData[bAddress]);
#endif
    }
  }

  return;
}



#endif
