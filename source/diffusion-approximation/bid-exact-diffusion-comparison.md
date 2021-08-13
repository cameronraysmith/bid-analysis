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

<center><font size="+4">Comparing the exact solution and diffusion approximation to the birth-immigration-death process</font></center>


# Diffusion limit via generating function


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

```maxima
texput(theta, "\\theta");
texput(tau, "(\\tau - \\tau_{0})");
```

```maxima
tex(genreparam);
```

```maxima
tex(factor(coeff(cgfexpansion, epsilon, 0)));
```

```maxima
tex(factor(coeff(cgfexpansion, epsilon, 2)));
```

```maxima
tex(factor(coeff(cgfexpansion, epsilon, 4)));
```

# Compute cumulants for exact and diffusion of the BD process

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

## First order terms

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
## Second order terms
<!-- #endregion -->

```maxima tags=[]
factor(coeff(difcum3, epsilon, 2));
```

```maxima tags=[]
factor(coeff(difcum4, epsilon, 2));
```

## Relations among terms

```maxima tags=[]
factor(2 * difcum3 * difcum1 - 3 * difcum2^2);
factor(difcum4 * difcum1^2 - 3 * difcum2^3 );
```
