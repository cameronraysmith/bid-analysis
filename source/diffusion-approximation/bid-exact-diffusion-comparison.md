---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.6
  kernelspec:
    display_name: Maxima
    language: maxima
    name: maxima
---

# Comparing the exact solution and diffusion approximation to the birth-immigration-death process

<!-- <center><font size="+4">Comparing the exact solution and diffusion approximation to the birth-immigration-death process</font></center> -->


## Diffusion limit via dimensionless generating function


Start with the generating function Z(z) for the BD process. 

```maxima
genfunc(z, tau, n0, sigma) :=
    (((sigma - 1)*(1 - %e^(-tau)) +
        (1 - sigma -(1 + sigma) * %e^(-tau)) * z) / 
    (-(1 + sigma) + (1 - sigma) * %e^(-tau) +
        (1 + sigma) * (1 - %e^(-tau)) * z))^(n0);
```

Reparameterize it in terms of s and epsilon = 1/N.

```maxima
geneps : genfunc(%e^(theta * epsilon), 
    tau, 1, s*epsilon);
```

Expand in powers of epsilon.

```maxima
cgfexpansion : (x0 / epsilon)*taylor(log(geneps), epsilon, 0, 3);
```

Collect the zeroth order term, which corresponds to the diffusion approximation and the correction term.

```maxima
cgfdiff : coeff(cgfexpansion, epsilon, 0);
cgfcorr : factor(coeff(cgfexpansion, epsilon, 2));
```

```maxima
tex(cgfdiff);
tex(cgfcorr);
```

Expand in powers of theta to generate the cumulants.

```maxima
cumser : taylor(log(genfunc(%e^(theta), tau, n0, sigma)),
    theta, 0, 4);
```

Extract the cumulants.

```maxima
kappa1 : factor(coeff(cumser, theta, 1));
kappa2 : factor(coeff(cumser, theta, 2));
kappa3 : factor(coeff(cumser, theta, 3));
kappa4 : factor(coeff(cumser, theta, 4));
```

Likewise obtain the cumulants from the rescaled generating function seperated into diffusion approximation and corrections.

```maxima
taylor(cgfdiff, theta, 0, 4);
taylor(cgfcorr, theta, 0, 4);
```

Form combinations of the cumulants.

```maxima tags=[]
block([R : expand(kappa1 * kappa3 / kappa2^2)],
    [factor(coeff(R, sigma, 0)),
    factor(coeff(R, sigma, 1)),
    factor(coeff(R, sigma, 2))]);
```

```maxima tags=[]
block([R : expand(kappa1^2 * kappa4 / kappa2^3)],
    [factor(coeff(R, sigma, 0)),
    factor(coeff(R, sigma, 1)),
    factor(coeff(R, sigma, 2))]);
```

```maxima tags=[]
factor(psubst(solve(y = (x+1)/(x-1), x), (3*x^2 - 2)/(3*(x - 1)^2)));
```

## Plot


```lisp
plot2d([[parametric, kap2(1, tau, 1), kap3(1, tau, 1), [tau, 0, 10]],
    [parametric, kap2(1, tau, 0.5), kap3(1, tau, 0.5), [tau, 0, 10]],
    [parametric, kap2(1, tau, 0), kap3(1, tau, 0), [tau, 0, 10]]],
    [legend, "N=1", "N=2", "N=inf"], [xlabel, "k_2/k_1^2"],
    [ylabel, "k_3/k_1^3"],
    [pdf_file, "fig/exact-cumulants-test.pdf"],
    [gnuplot_preamble, "set key left"],
    [gnuplot_pdf_term_command, "set term pdfcairo lw 3 size 17.2 cm, 12.9 cm font 'Latin Modern Roman,20'"]);
```

```maxima tags=[]
plot3d(cos(-x^2+y^3/4),[x,-4,4],[y,-4,4],
[gnuplot_pm3d,true],[grid,150,150],
[gnuplot_pdf_term_command, "set view map; unset surface; set term pdfcairo lw 3 size 17.2 cm, 12.9 cm font 'Latin Modern Roman,20'"]);
```

```maxima tags=[]
plot2d([contour, sigma^2*(3*%e^(2*tau)-2)/(3*(%e^tau-1)^2)], 
    [sigma, 0, 1], [tau, 0.02, 0.25], [levels, 100],
    [legend, false], [title, "fourth moment"],
    [pdf_file, "fig/example-contour.tex"],
    [gnuplot_pdf_term_command, "set term epslatex standalone lw 3 size 17.2 cm, 12.9 cm font 'Latin Modern Roman,20'"]);
```

<!-- #region tags=[] -->
### Absorbing boundary
<!-- #endregion -->

```maxima
genfunc(0, tau, x0/epsilon, s*epsilon);
```

```maxima
taylor(genfunc(0, tau, 1, s*epsilon), epsilon, 0, 2);
```

```maxima
genfunc(0, tau, n0, sigma);
```

```maxima
taylor((1 + s*epsilon)/(1 - s*epsilon), s, 0, 2);
```

```maxima tags=[]
plot2d([(1 - %e^(-tau))/(3 - %e^(-tau)),
        (1 - %e^(-tau))/(1.5 - %e^(-tau)),
        (1 - %e^(-tau))/(11/9 - %e^(-tau))], [tau, 0.1, 3],
    [legend, "sigma = 0.3", "sigma = 0.2", "sigma = 0.1"],
    [title, "P(0, tau) for n_0 = 0"],
    [pdf_file, "fig/p0n0.pdf"],
    [gnuplot_pdf_term_command, "set term pdfcairo lw 3 size 17.2 cm, 12.9 cm font 'Latin Modern Roman,20'"]);
```

```maxima tags=[]
plot2d([((1 - %e^(-tau))/(11/9 - %e^(-tau)))^10,
        %e^(-2/(1 - %e^(-tau)))], [tau, 0.1, 3],
    [legend, "sigma = 0.1, n_0 = 10", "diffusion"],
    [pdf_file, "fig/p0n10.pdf"],
    [gnuplot_pdf_term_command, "set term pdfcairo lw 3 size 17.2 cm, 12.9 cm font 'Latin Modern Roman,20'"]);
```

## Diffusion limit via generating function


Define the moment generating function.

```maxima tags=[]
genfunc (x, t, n0, kb, kd) := 
    ((kd * %e^(- (kb - kd) * t) - kd + 
        (kd - kb * %e^(- (kb - kd) * t)) * x) / 
        (kd * %e^(- (kb - kd) * t) - kb + 
        (kb - kb * %e^(- (kb - kd) * t)) * x))^n0;
```

Reparameterize it according to the diffusion approximation.

```maxima tags=[]
genreparam : factor(genfunc(%e^(-theta), tau*N, N * x0, D + s/(2*N), D - s/(2*N)));
```

Introduce epsilon = 1/N so as to make system size expansion.

```maxima tags=[]
geneps : genfunc(%e^(-theta * epsilon), tau / epsilon, 1, 
    D + s * epsilon/2, D - s* epsilon /2);
```

Expand the cumulant generating function.

```maxima tags=[]
cgfexpansion : (x0/epsilon)*taylor(log(geneps), epsilon, 0, 3);
```

```maxima
coeff(cgfexpansion, epsilon, 0);
```

```maxima tags=[]
texput(theta, "\\theta");
texput(tau, "(\\tau - \\tau_{0})");
```

```maxima tags=[]
tex(genreparam);
```

```maxima tags=[]
tex(factor(coeff(cgfexpansion, epsilon, 0)));
```

```maxima tags=[]
tex(factor(coeff(cgfexpansion, epsilon, 2)));
```

```maxima tags=[]
tex(factor(coeff(cgfexpansion, epsilon, 4)));
```

## Compute cumulants for exact and diffusion of the BD process

```maxima tags=[]
numop(n, gf) := 
    if n=0 then gf
        else numop(n-1, z * diff(gf, z))$
cumu(m, t, n0, kb, kd) := 
    factor(subst(1, z, numop(m, 
            log(genfunc(z, t, n0, kb, kd)))))$
```

```maxima tags=[]
cumu(0, t, n0, kb, kd);
```

```maxima tags=[]
cumu(1, t, n0, kb, kd);
cumu(2, t, n0, kb, kd);
cumu(3, t, n0, kb, kd);
cumu(4, t, n0, kb, kd);
```

### First order terms

```maxima tags=[]
difcum1 : taylor(epsilon * cumu(1, t / epsilon, x0 / epsilon, 
        D + s * epsilon / 2, D - s * epsilon / 2), epsilon, 0, 2);
difcum2 : taylor(epsilon^2 * cumu(2, t / epsilon, x0 / epsilon, 
        D + s * epsilon / 2, D - s * epsilon / 2), epsilon, 0, 2);
difcum3 : taylor(epsilon^3 * cumu(3, t / epsilon, x0 / epsilon, 
        D + s * epsilon / 2, D - s * epsilon / 2), epsilon, 0, 2);
difcum4 : taylor(epsilon^4 * cumu(4, t / epsilon, x0 / epsilon, 
        D + s * epsilon / 2, D - s * epsilon / 2), epsilon, 0, 2);
```

```maxima tags=[]
difcum1first : factor(taylor(epsilon * cumu(1, t / epsilon, x0 / epsilon, 
        D + s * epsilon / 2, D - s * epsilon / 2), epsilon, 0, 1));
difcum2first : factor(taylor(epsilon^2 * cumu(2, t / epsilon, x0 / epsilon, 
        D + s * epsilon / 2, D - s * epsilon / 2), epsilon, 0, 1));
difcum3first : factor(taylor(epsilon^3 * cumu(3, t / epsilon, x0 / epsilon, 
        D + s * epsilon / 2, D - s * epsilon / 2), epsilon, 0, 1));
difcum4first : factor(taylor(epsilon^4 * cumu(4, t / epsilon, x0 / epsilon, 
        D + s * epsilon / 2, D - s * epsilon / 2), epsilon, 0, 1));
```

```maxima tags=[]
difcum1 : subst(0, epsilon, epsilon * cumu(1, t / epsilon, x0 / epsilon, 
        D + s * epsilon / 2, D - s * epsilon / 2));
difcum2 : subst(0, epsilon, epsilon^2 * cumu(2, t / epsilon, x0 / epsilon, 
        D + s * epsilon / 2, D - s * epsilon / 2));
difcum3 : subst(0, epsilon, epsilon^3 * cumu(3, t / epsilon, x0 / epsilon, 
        D + s * epsilon / 2, D - s * epsilon / 2));
difcum4 : subst(0, epsilon, epsilon^4 * cumu(4, t / epsilon, x0 / epsilon, 
        D + s * epsilon / 2, D - s * epsilon / 2));
```

<!-- #region tags=[] -->
### Second order terms
<!-- #endregion -->

```maxima tags=[]
factor(coeff(difcum3, epsilon, 2));
```

```maxima tags=[]
factor(coeff(difcum4, epsilon, 2));
```

### Relations among terms

```maxima tags=[]
factor(2 * difcum3 * difcum1 - 3 * difcum2^2);
factor(difcum4 * difcum1^2 - 3 * difcum2^3 );
```

## Compute cumulants for inhomogeneous component of BID process


Cumulant generating function for inhomogenous equation.

```maxima
frac1 : (kb - kd * %e^(-(kb-kd)*t)) / ((kb - kd) * %e^(-(kb-kd)*t));
frac2 : (kb - kb * %e^(-(kb-kd)*t)) / (kb - kd * %e^(-(kb-kd)*t));
```

```maxima
inhfrac : frac1 * (1 - frac2 * %e^(theta * epsilon));
```

```maxima
psubst([kb = D + s/(2*N), kd = D - s/(2*N), t = N*tau], frac1);
psubst([kb = D + s/(2*N), kd = D - s/(2*N), t = N*tau], frac2);
```

```maxima
psubst([kb = D + s/(2*N), kd = D - s/(2*N), t = N*tau], inhfrac);
```

```maxima
inhgenf : psubst([kb = D + s * epsilon/2, kd = D - s * epsilon/2, 
         t = tau/epsilon], (M/kb) * log(inhfrac));
```

```maxima tags=[]
inhgenser : taylor(psubst([kb = D + s * epsilon/2, kd = D - s * epsilon/2, 
         t = tau/epsilon], (M/kb) * log(inhfrac)), epsilon, 0, 2);
```

```maxima
logcontract(coeff(inhgenser, epsilon, 0));
```

```maxima
logcontract(coeff(inhgenser, epsilon, 1));
```

```maxima
expand(logcontract(coeff(inhgenser, epsilon, 1)));
```

```maxima tags=[]
piece1 : (D*M*s*%e^(s*tau)*theta*log(-s/(D*%e^(s*tau)*theta-D*theta-s)))/(2*D^3*%e^(s*tau)*theta-2*D^3*theta-2*D^2*s)+(D*M*s*theta*log(-(D*%e^(s*tau)*theta)/s+(D*theta)/s+1))/(2*D^3*%e^(s*tau)*theta-2*D^3*theta-2*D^2*s)+(M*s^2*log(-(D*%e^(s*tau)*theta)/s+(D*theta)/s+1))/(2*D^3*%e^(s*tau)*theta-2*D^3*theta-2*D^2*s);
piece2 : (D^2*M*%e^(s*tau)*theta^2)/(2*D^3*%e^(s*tau)*theta-2*D^3*theta-2*D^2*s)-(D^2*M*theta^2)/(2*D^3*%e^(s*tau)*theta-2*D^3*theta-2*D^2*s)+(D*M*s*%e^(s*tau)*theta)/(2*D^3*%e^(s*tau)*theta-2*D^3*theta-2*D^2*s)-(D*M*s*theta)/(2*D^3*%e^(s*tau)*theta-2*D^3*theta-2*D^2*s);
```

```maxima
factor(piece1);
```

```maxima
factor(piece2);
```

```maxima
(M*s/(2*D^2*(D*%e^(s*tau)*theta-D*theta-s))) *
(- D * %e^(s*\tau) * theta + D*theta + s);
```

```maxima
factor((M*s/(2*D^2*(D*%e^(s*tau)*theta-D*theta-s))) *
(- D * %e^(s*\tau) * theta + D*theta + s));
```

```maxima
inhgenexpanded : expand(taylor(inhgenser, theta, 0, 4));
```

```maxima tags=[]
kap1 : factor(coeff(inhgenexpanded, theta, 1));
kap2 : factor(coeff(coeff(inhgenexpanded, theta, 2), epsilon, 0)) +
    (1/N) * factor(coeff(coeff(inhgenexpanded, theta, 2), epsilon, 1));
kap3 : factor(coeff(coeff(inhgenexpanded, theta, 3), epsilon, 0)) +
    (1/N) * factor(coeff(coeff(inhgenexpanded, theta, 3), epsilon, 1));
kap4 : factor(coeff(coeff(inhgenexpanded, theta, 4), epsilon, 0)) +
    (1/N) * factor(coeff(coeff(inhgenexpanded, theta, 4), epsilon, 1));
```

## Compare exact and approximate cumulants

```maxima tags=[]
kap2 (D, tau, epsilon) := 2* D * (1 - %e^(-tau));
kap3 (D, tau, epsilon) := 6* D^2 * (1 - %e^(-tau))^2 +
    (epsilon/2) * (1 - %e^(-2*tau));
```

```maxima tags=[]
plot2d([[parametric, kap2(1, tau, 1), kap3(1, tau, 1), [tau, 0, 10]],
    [parametric, kap2(1, tau, 0.5), kap3(1, tau, 0.5), [tau, 0, 10]],
    [parametric, kap2(1, tau, 0), kap3(1, tau, 0), [tau, 0, 10]]],
    [legend, "N=1", "N=2", "N=inf"], [xlabel, "k_2/k_1^2"],
    [ylabel, "k_3/k_1^3"],
    [pdf_file, "fig/exact-cumulants-test.pdf"],
    [gnuplot_preamble, "set key left"],
    [gnuplot_pdf_term_command, "set term pdfcairo lw 3 size 17.2 cm, 12.9 cm font 'Latin Modern Roman,20'"]);
```

<!-- #region tags=[] -->
```lisp
plot2d([[parametric, kap2(1, tau, 1), kap3(1, tau, 1), [tau, 0, 10]],
    [parametric, kap2(1, tau, 0.5), kap3(1, tau, 0.5), [tau, 0, 10]],
    [parametric, kap2(1, tau, 0), kap3(1, tau, 0), [tau, 0, 10]]],
    [legend, "$N=1$", "$N=2$", "$N=\\infty$"], [xlabel, "$\\frac{\\kappa_2}{\\kappa_1^2}$"],
    [ylabel, "$\\frac{\\kappa_3}{\\kappa_1^3}$"],
    [gnuplot_preamble, "set key top left"],
    [gnuplot_term, "cairolatex pdf lw 3 size 17.2 cm, 12.9 cm font 'Latin Modern Roman,20'"], 
    [gnuplot_out_file, "fig/exact-cumulants-test.tex"]);
```
<!-- #endregion -->
