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

## Recovering the distribution


Re-express the generating function in a form that makes it easy to expand as a geometric series.

```maxima
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

```maxima
factor(apiece + bpiece * z / (1 - cpiece * z)- z0sol);
```

Re-express the generating function in a form that makes it easy to expand as a binomial series.

```maxima
psubst([t = 1.0, sigma = 0.2], [apiece, bpiece, cpiece]);
```

```maxima
appi : factor(coeff(expand(num(z0sol)), z, 1) /
    coeff(expand(num(z0sol)), z, 0));
bppi : factor(-coeff(expand(denom(z0sol)), z, 1) /
    coeff(expand(denom(z0sol)), z, 0));
cppi : factor(coeff(expand(num(z0sol)), z, 0) /
    coeff(expand(denom(z0sol)), z, 0));
```

```maxima
factor(cppi * (1 + appi * z) / (1 - bppi * z) - z0sol);
```

Check agreement numerically

```maxima tags=[]
psubst([sigma = 0.2, t = 1.0], taylor(newz0sol, z, 0, 6));
```
