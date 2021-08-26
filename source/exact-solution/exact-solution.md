---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.4
  kernelspec:
    display_name: Maxima
    language: maxima
    name: maxima
---

# Exact solulution to the BID process


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
