---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.3
  kernelspec:
    display_name: Maxima
    language: maxima
    name: maxima
---

# Moment-based approximation to the BID process


This is the exact generating function solution to the birth-immigration-death process.

```maxima
genfunc (x, t, n0, kb, kd) := 
    ((kd * %e^(- (kb - kd) * t) - kd + 
        (kd - kb * %e^(- (kb - kd) * t)) * x) / 
        (kd * %e^(- (kb - kd) * t) - kb + 
        (kb - kb * %e^(- (kb - kd) * t)) * x))^n0;
```

```maxima tags=[]
probdist (m, t, n0, kb, kd) :=
block([gf],
    gf : genfunc(z, t, n0, kb, kd),
    makelist([j, subst (0, z, diff(gf, z, j) / (j!))],
        j, 0, m))$
```

The trial distribution.

```maxima
trialdist : (c0 + c1*x + c2*x^2 + c3*x^3) * %e^(-k*x);
```

The moments of the trial distribution.

```maxima
m0 : -subst(0, x, integrate(trialdist, x));
m1 : -subst(0, x, integrate(x * trialdist, x));
m2 : -subst(0, x, integrate(x^2 * trialdist, x));
m3 : -subst(0, x, integrate(x^3 * trialdist, x));
m4 : -subst(0, x, integrate(x^4 * trialdist, x));
```

```maxima tags=[]
simeqs: solve([m1 = v1, m2 = v2, m3 = v3, m4 = v4], [c0, c1, c2, c3])$
```

```maxima tags=[]
simeqs[1][1];
simeqs[1][2];
simeqs[1][3];
simeqs[1][4];
```

As an example, consider the parameters n0 = 3, t = 1, kb = 1.1, kd = 0.9

```maxima tags=[]
exdist : probdist(30, 1, 3, 1.1, 0.9)$
```

```maxima tags=[]
exdist[1];
exdist[2];
exdist[3];
```

Compute the moments for the exact distribution.

```maxima
numop(n, gf) := 
    if n=0 then gf
        else numop(n-1, z * diff(gf, z))$
mom(m, t, n0, kb, kd) := 
    factor(subst(1, z, numop(m, genfunc(z, t, n0, kb, kd))))$
```

```maxima
psubst([kb = 1.1, kd = 0.9], mom(1, 1, 3, kb, kd));
psubst([kb = 1.1, kd = 0.9], mom(2, 1, 3, kb, kd));
psubst([kb = 1.1, kd = 0.9], mom(3, 1, 3, kb, kd));
psubst([kb = 1.1, kd = 0.9], mom(4, 1, 3, kb, kd));
```

```maxima
sum(exdist[k][2], k, 1, 31);
sum((k - 1) * exdist[k][2], k, 1, 31);
sum((k - 1)^2 * exdist[k][2], k, 1, 31);
sum((k - 1)^3 * exdist[k][2], k, 1, 31);
sum((k - 1)^4 * exdist[k][2], k, 1, 31);
```

<!-- #region tags=[] -->
Solve for the parameters.  We frist solve for the c's because those equations are linear.
<!-- #endregion -->

```maxima tags=[]
csol : solve([m1 = 3.664, m2 = 21.54, m3 = 166.22, m4 = 1578.7], [c0, c1, c2, c3])$
```

```maxima tags=[]
csol[1][1];
csol[1][2];
csol[1][3];
csol[1][4];
```

We then solve a polynomial equation for k.

```maxima tags=[]
normer : factor(psubst(csol, m0) - 1);
allroots(normer);
```

```maxima tags=[]
ccsol : cons(k = 0.3313, subst(0.3313 , k, csol[1]));
```

<!-- #region -->


Check the answer.
<!-- #endregion -->

```maxima
psubst(ccsol, [m0, m1, m2, m3, m4]);
```

Substitute the numerical values into the trial distrbution.

```maxima
approxdist : psubst(ccsol, trialdist);
```

Plot against the exact solution.

```maxima tags=[]
set_plot_option([svg_file, "maxplot-moment-gen.svg"])$
```

```maxima tags=[]
plot2d([approxdist, [discrete, exdist]], [x, 1, 16],
[legend, "approx", "exact"], [style, lines, points]);
```
