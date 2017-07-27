double griewankFunction(struct position x){

  double  f = 0;
  double  p = 1;
  double xd;
  int d;

  for (d = 0; d < x.size; d++){
    xd = x.x[d];
    f = f + xd * xd;
    p = p * cos (xd / sqrt ((double) (d + 1)));
  }

  f = f / 4000 - p + 1;

  return f;

}

double rastriginFunction(struct position x){

  int d;
  int k = 10;
  double f = 0;
  double xd;

  for (d = 0; d < x.size; d++){
    xd = x.x[d];
    f =f+ xd * xd - k * cos (2 * pi * xd);
  }
  f =f+ x.size * k;

  return f;

}
