---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.2
  kernelspec:
    display_name: Maxima
    language: maxima
    name: maxima
---

<center><font size="+4">Generating function analysis of the diffusion approximation to the birth-immigration-death process</font></center>

<!-- #region tags=[] -->
# Moment generating function
<!-- #endregion -->

```maxima tags=[]
genfunc (x, t, n0, kb, kd) := 
    ((kd * %e^(- (kb - kd) * t) - kd + 
        (kd - kb * %e^(- (kb - kd) * t)) * x) / 
        (kd * %e^(- (kb - kd) * t) - kb + 
        (kb - kb * %e^(- (kb - kd) * t)) * x))^n0;
```

Extract the probability distribution

```maxima tags=[]
probdist (m, t, n0, kb, kd) :=
block([gf],
    gf : genfunc(z, t, n0, kb, kd),
    makelist([j, subst (0, z, diff(gf, z, j) / (j!))],
        j, 0, m))$
```

Reparameterize it according to the diffusion approximation.

```maxima tags=[]

factor(genfunc(z / N, N * tau, N * x0, D + s/(2*N), D - s/(2*N)));

```

```maxima tags=[]

renoprob (m, N, tau, x0, D, s) :=
    map (lambda ([pa], [pa[1] / N, pa[2] * N]),
    probdist (m, N * tau, N * x0, D + s/(2*N), D - s/(2*N)))$

```

Plot with values tau = 1.0, x0 = 3.0, D = 1.0, s = 1.0

```maxima
set_plot_option([svg_file, "maxplot.svg"])$
```

```maxima
plot2d([[discrete, renoprob(30, 1, 1.0, 3, 1.0, 1.0)],
    [discrete, renoprob(60, 2, 1.0, 3, 1.0, 1.0)],
    [discrete, renoprob(90, 3, 1.0, 3, 1.0, 1.0)],
    [discrete, renoprob(120, 4, 1.0, 3, 1.0, 1.0)]],
    [legend, "N=1", "N=2", "N=3", "N=4"],
    [xlabel, "x"]);
```

Check that it is normalized and extract fixation probability.

```maxima tags=[]
genfunc(0, t, n0, kb, kd);
factor(genfunc(0, N * tau, N * x0, D + s/(2*N), D - s/(2*N)));
genfunc(1, t, n0, kb, kd);
```

<!-- #region tags=[] toc-hr-collapsed=true -->
# Functions to compute moments
<!-- #endregion -->

```maxima tags=[]
numop(n, gf) := 
    if n=0 then gf
        else numop(n-1, z * diff(gf, z))$
mom(m, t, n0, kb, kd) := 
    factor(subst(1, z, numop(m, genfunc(z, t, n0, kb, kd))))$
momren(m, t, n0, kb, kd) := 
                factor(psubst([t = N * tau, 
                        n0 = N * x0,
                        kb = D + s/(2*N), 
                        kd = D - s/(2*N)], 
                mom(m, t, n0, kb, kd) / N^m))$
cumu(m, t, n0, kb, kd) := 
    factor(subst(1, z, numop(m, 
            log(genfunc(z, t, n0, kb, kd)))))$
cumuren(m, zt, n0, kb, kd) := 
                factor(psubst([t = N * tau, 
                        n0 = N * x0,
                        kb = D + s/(2*N), 
                        kd = D - s/(2*N)], 
                cumu(m, t, n0, kb, kd) / N^m))$
```

<!-- #region tags=[] toc-hr-collapsed=true -->
# Compute cumulants
<!-- #endregion -->

```maxima tags=[]
cumu(1, t, n0, kb, kd);
cumu(2, t, n0, kb, kd);
cumu(3, t, n0, kb, kd);
cumu(4, t, n0, kb, kd);

```

Reparameterize cumulants and collect powers of N.

```maxima tags=[]

cumuren(1, t, n0, kb, kd);
cumuren(2, t, n0, kb, kd);

```

```maxima tags=[]
block([mo],
    mo : expand(cumuren(3, t, n0, kb, kd)),
    factor(coeff (mo, N, 0)) +
    factor(coeff (mo, N, -2)) / N^2);
block([mo],
    mo : expand(cumuren(4, t, n0, kb, kd)),
    factor(coeff (mo, N, 0)) +
    factor(coeff (mo, N, -2)) / N^2);

```

```maxima
plot2d(
    subst(1.0, s, 
        (%e^(s*tau)-1)*%e^(s*tau)*(%e^(s*tau)+1)/2),
[tau, 0, 1]);
```

<!-- #region tags=[] toc-hr-collapsed=true -->
# Appendix
<!-- #endregion -->

<!-- #region tags=[] -->
## Original maxima code
<!-- #endregion -->

```maxima
/* [wxMaxima batch file version 1] [ DO NOT EDIT BY HAND! ]*/
/* [ Created with wxMaxima version 18.02.0 ] */
/* [wxMaxima: comment start ]
Moment generating function.
   [wxMaxima: comment end   ] */


/* [wxMaxima: input   start ] */
genfunc (x, t, n0, kb, kd) := 
    ((kd * %e^(- (kb - kd) * t) - kd + 
        (kd - kb * %e^(- (kb - kd) * t)) * x) / 
        (kd * %e^(- (kb - kd) * t) - kb + 
        (kb - kb * %e^(- (kb - kd) * t)) * x))^n0;
/* [wxMaxima: input   end   ] */


/* [wxMaxima: comment start ]
Extract the probability distribution.
   [wxMaxima: comment end   ] */


/* [wxMaxima: input   start ] */
probdist (m, t, n0, kb, kd) :=
block([gf],
    gf : genfunc(z, t, n0, kb, kd),
    makelist([j, subst (0, z, diff(gf, z, j) / (j!))],
        j, 0, m))$
/* [wxMaxima: input   end   ] */


/* [wxMaxima: comment start ]
Reparameterize it according to diffusion approximation.
   [wxMaxima: comment end   ] */


/* [wxMaxima: input   start ] */
factor(genfunc(z / N, N * tau, N * x0, D + s/(2*N), D - s/(2*N)));
/* [wxMaxima: input   end   ] */


/* [wxMaxima: input   start ] */
renoprob (m, N, tau, x0, D, s) :=
    map (lambda ([pa], [pa[1] / N, pa[2] * N]),
    probdist (m, N * tau, N * x0, D + s/(2*N), D - s/(2*N)))$
/* [wxMaxima: input   end   ] */


/* [wxMaxima: comment start ]
Plot with values tau = 1.0, x0 = 3.0, D = 1.0, s = 1.0
   [wxMaxima: comment end   ] */


/* [wxMaxima: input   start ] */
wxplot2d([[discrete, renoprob(30, 1, 1.0, 3, 1.0, 1.0)],
    [discrete, renoprob(60, 2, 1.0, 3, 1.0, 1.0)],
    [discrete, renoprob(90, 3, 1.0, 3, 1.0, 1.0)],
    [discrete, renoprob(120, 4, 1.0, 3, 1.0, 1.0)]],
    [legend, "N=1", "N=2", "N=3", "N=4"],
    [xlabel, "x"]);
/* [wxMaxima: input   end   ] */


/* [wxMaxima: comment start ]
Check that it is normalized and extract
fixation probability.
   [wxMaxima: comment end   ] */


/* [wxMaxima: input   start ] */
genfunc(0, t, n0, kb, kd);
factor(genfunc(0, N * tau, N * x0, D + s/(2*N), D - s/(2*N)));
genfunc(1, t, n0, kb, kd);
/* [wxMaxima: input   end   ] */


/* [wxMaxima: comment start ]
Functions to compute moments.
   [wxMaxima: comment end   ] */


/* [wxMaxima: input   start ] */
numop(n, gf) := 
    if n=0 then gf
        else numop(n-1, z * diff(gf, z))$
mom(m, t, n0, kb, kd) := 
    factor(subst(1, z, numop(m, genfunc(z, t, n0, kb, kd))))$
momren(m, t, n0, kb, kd) := 
                factor(psubst([t = N * tau, 
                        n0 = N * x0,
                        kb = D + s/(2*N), 
                        kd = D - s/(2*N)], 
                mom(m, t, n0, kb, kd) / N^m))$
cumu(m, t, n0, kb, kd) := 
    factor(subst(1, z, numop(m, 
            log(genfunc(z, t, n0, kb, kd)))))$
cumuren(m, zt, n0, kb, kd) := 
                factor(psubst([t = N * tau, 
                        n0 = N * x0,
                        kb = D + s/(2*N), 
                        kd = D - s/(2*N)], 
                cumu(m, t, n0, kb, kd) / N^m))$
/* [wxMaxima: input   end   ] */


/* [wxMaxima: comment start ]
Have a look at the cumulants.
   [wxMaxima: comment end   ] */


/* [wxMaxima: input   start ] */
cumu(1, t, n0, kb, kd);
cumu(2, t, n0, kb, kd);
cumu(3, t, n0, kb, kd);
cumu(4, t, n0, kb, kd);
/* [wxMaxima: input   end   ] */


/* [wxMaxima: comment start ]
Reparameterize the cumulants and collect 
powers of N.
   [wxMaxima: comment end   ] */


/* [wxMaxima: input   start ] */
cumuren(1, t, n0, kb, kd);
cumuren(2, t, n0, kb, kd);
/* [wxMaxima: input   end   ] */


/* [wxMaxima: input   start ] */
block([mo],
    mo : expand(cumuren(3, t, n0, kb, kd)),
    factor(coeff (mo, N, 0)) +
    factor(coeff (mo, N, -2)) / N^2);
block([mo],
    mo : expand(cumuren(4, t, n0, kb, kd)),
    factor(coeff (mo, N, 0)) +
    factor(coeff (mo, N, -2)) / N^2);
/* [wxMaxima: input   end   ] */


/* [wxMaxima: input   start ] */
wxplot2d(
    subst(1.0, s, 
        (%e^(s*tau)-1)*%e^(s*tau)*(%e^(s*tau)+1)/2),
[tau, 0, 1]);
/* [wxMaxima: input   end   ] */

"Created with wxMaxima 18.02.0"$
```


```maxima

```
