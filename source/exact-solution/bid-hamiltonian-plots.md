---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.7
  kernelspec:
    display_name: Maxima
    language: maxima
    name: maxima
---

# Hamiltonian plots


## Characteristics

```maxima tags=[]
plot2d(map(lambda([theta], log((%e^(theta - t))/(1-(1-%e^(-t))*%e^(theta)))), 
        [-0.25, -0.375, -0.5, -0.625, -0.75, -0.875, -0.125, 0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875]),
        [t, 0.01, 1.0], [y, -1, 1], [legend, false],
        [xlabel, "t"], [ylabel, "𝜃"],
        [color, "#1F654C", "#226F54", "#559972", "#6EAE81", "#87C38F", "#BEDAA5"],
        [pdf_file, "fig/hamiltonian-characteristics.pdf"],
        [gnuplot_pdf_term_command, "set term pdfcairo enhanced lw 4.5 size 17.2 cm, 12.9 cm font 'Latin Modern Roman,30'"]);
```

<!-- #region tags=[] -->
## phase portrait
<!-- #endregion -->

```maxima tags=[]
plot2d(map(lambda([h], h/(%e^theta - 1)), 
        [0.1, 0.25, 0.5, 1.0, 1.75, 2.5, 3.5, 5.0]),
        [theta, 0.01, 1.0], [y, 0.0, 4.0], [legend, false],
        [xlabel, "𝜃"], [ylabel, "n"],
        [color, "#1F654C", "#226F54", "#559972", "#6EAE81", "#87C38F", "#BEDAA5"],
        [pdf_file, "fig/hamiltonian-phase-portrait.pdf"],
        [gnuplot_pdf_term_command, "set term pdfcairo enhanced lw 4.5 size 17.2 cm, 12.9 cm font 'Latin Modern Roman,30'"]);
```

<!-- #region tags=[] -->
## CGF
<!-- #endregion -->

```maxima
Gam(theta,t) := theta - t + log(1 - %e^theta + %e^(theta - t));
```

```maxima tags=[]
plot2d([Gam(theta, 0.5), Gam(theta, 1), Gam(theta, 2.5), Gam(theta, 5)], 
    [theta, -5,0], [legend, "0.5", "1.0", "2.5", "5.0"], 
    [xlabel, "𝜃"], [ylabel, "Γ(𝜃)"],
    [color, "#1F654C", "#226F54", "#559972", "#6EAE81", "#87C38F", "#BEDAA5"],
    [pdf_file, "fig/hamiltonian-cgf.pdf"],
    [gnuplot_preamble, "
set nokey
set label '0.5' at -3, -2.7
set label '1.0' at -1, -2.6
set label '2.5' at -1, -4.3
set label '5.0' at -1, -6.9
"],
    [gnuplot_pdf_term_command, "set term pdfcairo enhanced lw 4.5 size 17.2 cm, 12.9 cm font 'Latin Modern Roman,30'"]);
```

## MGF

```maxima
 Z(theta,t) := %e^(theta - t)/(1 - %e^theta + %e^(theta - t));
```

```maxima tags=[]
plot2d([Z(theta, 0.5), Z(theta, 1), Z(theta, 2.5), Z(theta, 5)], 
    [theta, -5,0], [legend, "0.5", "1.0", "2.5", "5.0"], 
    [xlabel, "𝜃"], [ylabel, "Z(𝜃)"],
    [color, "#1F654C", "#226F54", "#559972", "#6EAE81", "#87C38F", "#BEDAA5"],
    [pdf_file, "fig/hamiltonian-mgf.pdf"],
    [gnuplot_preamble, "
set key leftb top
"],
    [gnuplot_pdf_term_command, "set term pdfcairo enhanced lw 4.5 size 17.2 cm, 12.9 cm font 'Latin Modern Roman,30'"]
    );
```

## Distribution

```maxima
Prob(n,t) := if n=0 then 0 else %e^(-t) * (1 - %e^(-t))^(n - 1);
```

```maxima
plot2d(map(lambda([t], [discrete, makelist(Prob(n, t), n, 0, 15)]), 
        [0.5, 1.0, 2.5, 5.0]), [xlabel, "n"], [ylabel, "P(n)"],
        [legend, "0.5", "1.0", "2.5", "5.0"],
        [color, "#1F654C", "#226F54", "#559972", "#6EAE81", "#87C38F", "#BEDAA5"],
    [pdf_file, "fig/hamiltonian-pdf.pdf"],
    [gnuplot_preamble, "
set key right top
"],
    [gnuplot_pdf_term_command, "set term pdfcairo enhanced lw 4.5 size 17.2 cm, 12.9 cm font 'Latin Modern Roman,30'"]);
```

## Solution of the BID process using the method of characteristics

<!-- #region tags=[] -->
### Homogeneous equation
<!-- #endregion -->

Homogeneous PDE statisfied by the CGF.

```maxima
homcgfeqn (W) := diff(W, t) - 
    ((1 + sigma) / (2 * sigma)) * (z^2 - z) * diff(W, z) -
    ((1 - sigma) / (2 * sigma)) * (1 - z) * diff(W, z);
```

ODE for characteristics of the generating function.

```maxima
chareqn(z) := diff(z, t) + ((1 + sigma) / (2 * sigma)) * (z^2 - z) +
    ((1 - sigma) / (2 * sigma)) * (1 - z);
```

Expand the integrand for solving the characteristic equation as partial fractions.

```maxima
partfrac (-1 / (((1 + sigma) / (2 * sigma)) * (z^2 - z) +
    ((1 - sigma) / (2 * sigma)) * (1 - z)), z);
```

From the implicit solution of the characteristic equation, obtain z and z0 explicitly.

```maxima
zsol : part(solve ([(z - (1 - sigma) / (1 + sigma)) / (z - 1) =
    %e^t * (z0 - (1 - sigma) / (1 + sigma)) / (z0 - 1)], [z]), 1, 2);
z0sol : part(solve ([(z - (1 - sigma) / (1 + sigma)) / (z - 1) =
    %e^t * (z0 - (1 - sigma) / (1 + sigma)) / (z0 - 1)], [z0]), 1, 2);
```

Double check that the explicit solution indeed satisfied the ODE.

```maxima
factor(chareqn(zsol));
```

Reorder the solution for z for the publication.

```maxima
factor(coeff(expand(num(zsol)), z0, 1)) * z0;
factor(coeff(expand(num(zsol)), z0, 0));
factor(coeff(expand(denom(zsol)), z0, 1)) *z0;
factor(coeff(expand(denom(zsol)), z0, 0));
```

```maxima
newzsol :
((-(sigma + 1) + (1 - sigma) * %e^(-t)) * z0 +
    (1 - sigma) * (1 - %e^(-t))) /
(-(sigma + 1) * (1 - %e^(-t)) * z0 +
    1 - sigma - (1 + sigma) * %e^(-t));
```

```maxima
factor(zsol - newzsol);
```

Likewise, reorder the solution for z for typesetting.

```maxima
factor(coeff(expand(num(z0sol)), z, 1)) * z;
factor(coeff(expand(num(z0sol)), z, 0));
factor(coeff(expand(denom(z0sol)), z, 1)) *z;
factor(coeff(expand(denom(z0sol)), z, 0));
```

```maxima
newz0sol : (((sigma - 1) * (1 - %e^(-t))
   + (1 - sigma - (1 + sigma) * %e^(-t)) * z) / 
    ((-1 - sigma + (1 - sigma) * %e^(-t)) + 
        (1 + sigma) * (1 - %e^(-t)) * z));
```

```maxima
factor(z0sol - newz0sol);
```

Check that the solution satsfies the PDE.

```maxima
factor(homcgfeqn(z0sol));
```

<!-- #region tags=[] -->
### Inhomogeneous equation
<!-- #endregion -->

Find characteristic passing through point $(t, z)$ parameterized by $t_p$.

```maxima tags=[]
charz : factor(psubst([t = tp, z0 = z0sol], zsol));
```

```maxima
subst(t, tp, charz);
```

Pull off coefficients for typesetting.

```maxima
factor(coeff(num(charz), z, 0));
factor(coeff(num(charz), z, 1));
factor(coeff(denom(charz), z, 0));
factor(coeff(denom(charz), z, 1));
```

Carry out the integral.

```maxima
block([indef],
    indef : integrate(charz - 1, tp),
    charint : logcontract (subst(t, tp, indef) - subst(0, tp, indef)));
```

```maxima
part(charint, 1, 1, 1, 2, 1, 1, 2, 1);
```

Re-express it for typesetting.

```maxima
logcontract(ev(
((2 * sigma) / (1 + sigma)) * 
    log((2 * sigma * %e^(-t)) /
        ((1 + sigma) * (1 - %e^(-t)) * z
        - 1 - sigma + (1 - sigma) * %e^(-t)))
        - charint,
    logexpand = super));
```

Check that it satisfies the original CGF equation.

```maxima
factor(homcgfeqn(charint));
```

## Recovering the distribution

<!-- #region tags=[] -->
### Homogeneous equation
<!-- #endregion -->

Re-express the generating function in a form that makes it easy to expand as a geometric series.

```maxima tags=[]
apiece : factor(coeff(expand(num(z0sol)), z, 0) /
    coeff(expand(denom(z0sol)), z, 0));
bpiece : factor((coeff(expand(num(z0sol)), z, 1) * 
    coeff(expand(denom(z0sol)), z, 0) -
    coeff(expand(num(z0sol)), z, 0) * 
    coeff(expand(denom(z0sol)), z, 1)) / 
    coeff(expand(denom(z0sol)), z, 0)^2);
cpiece : factor(-coeff(expand(denom(z0sol)), z, 1) /
    coeff(expand(denom(z0sol)), z, 0));
```

```maxima tags=[]
factor(apiece + bpiece * z / (1 - cpiece * z)- z0sol);
```

Re-express the generating function in a form that makes it easy to expand as a binomial series.

```maxima tags=[]
psubst([t = 1.0, sigma = 0.2], [apiece, bpiece, cpiece]);
```

```maxima tags=[]
appi : factor(coeff(expand(num(z0sol)), z, 1) /
    coeff(expand(num(z0sol)), z, 0));
bppi : factor(-coeff(expand(denom(z0sol)), z, 1) /
    coeff(expand(denom(z0sol)), z, 0));
cppi : factor(coeff(expand(num(z0sol)), z, 0) /
    coeff(expand(denom(z0sol)), z, 0));
```

```maxima tags=[]
factor(cppi * (1 + appi * z) / (1 - bppi * z) - z0sol);
```

Compute example numerical values of the Taylor series expansion of $z_0$

```maxima tags=[]
psubst([sigma = 0.2, t = 1.0], taylor(newz0sol, z, 0, 6));
```

<!-- #region tags=[] -->
### Inhomogeneous equation
<!-- #endregion -->

```maxima
inhmgfrac : 
((1 + sigma) * (1 - %e^(-t)) * z
        - 1 - sigma + (1 - sigma) * %e^(-t)) / 
    (2 * sigma * %e^(-t));
```

```maxima tags=[]
factor(coeff(expand(inhmgfrac), z, 0));
factor(-coeff(expand(inhmgfrac), z, 1));
```
