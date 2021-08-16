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

# Solving the diffusion approximation to the birth-immigration-death process

<!-- <center><font size="+4">Solving the diffusion approximation to the birth-immigration-death process</font></center> -->


## Forward Laplace transform

```maxima tags=[]
diffusioneqn (rho, D, s, M) := 
    diff (rho, tau) -
    D * diff (x * rho, x, 2) +
    s * diff (x * rho, x, 1) +
    M * diff (rho, x);
```

```maxima tags=[]
cgfeqn (gam, D, s, M) := 
    diff (gam, tau) +
    (D * theta^2 - s* theta) * diff (gam, theta) -
    M * theta;
```

```maxima tags=[]
factor(diffusioneqn (rho0(x,tau), D, s, M));
```

## Characteristic equation


Compute the integral for the implicit solution.

```maxima tags=[]
logcontract(integrate(1/(D * thp^2 - s * thp), thp));
```

Solve the implicit equation.

```maxima tags=[]
thetasol : solve ([((D* theta - s) * theta0) / ((D* theta0 - s) * theta) =
    %e^(s*u)], [theta])[1];
theta0sol : solve ([((D* theta - s) * theta0) / ((D* theta0 - s) * theta) =
    %e^(s*u)], [theta0])[1];
```

Check the solution.

```maxima tags=[]
subst(0, u,  thetasol);
```

```maxima tags=[]
chareqn(theta) := diff(theta, u) - D * theta^2 + s * theta;
```

```maxima tags=[]
factor(chareqn(part(thetasol, 2)));
```

## Homogeneous solution


Re-expressing $u$ in terms of $\tau$ and including the initial condition $P(x,0) = \delta(x-x_0) \Rightarrow \Gamma_{\mathrm{hom}}(\theta, 0)=x_0 \, \theta$, where $\Rightarrow$ indicates taking the logarithm of the Laplace transform, yields

```maxima tags=[]
gamhom : x0 * s * theta * %e^(s*tau) / 
    (D * (%e^(s*tau) - 1) * theta + s);
```

```maxima tags=[]
factor(gamhom - x0 * subst(tau, u, part(theta0sol, 2)));
```

Noting $M=0$ in the homogeneous case

```maxima tags=[]
factor(cgfeqn (gamhom, D, s, 0));
```

## Inhomogeneous solutions


Find the characteristic which passes through a given point (theta_f, tau_f).

```maxima tags=[]
th0sol : solve([thetaf = 
        subst(uf, u, part(thetasol, 2))], [theta0])[1];
```

```maxima tags=[]
thfsol : factor(psubst(th0sol, thetasol));
```

```maxima tags=[]
subst(uf, u, thfsol);
```

Integrate along that characteristic.

```maxima tags=[]
ginhint : integrate((thetaf*s*%e^(uf*s))/
    (s*%e^(s*u)-D*thetaf*%e^(s*u)+D*thetaf*%e^(uf*s)), u);
```

```maxima tags=[]
logcontract(ev(factor(subst(uf, u, ginhint) - subst(0, u, ginhint)), 
        logexpand=super));
```

```maxima tags=[]
logcontract(ev(factor(subst(uf, u, ginhint) - subst(0, u, ginhint)), 
        logexpand=super));
```

Adding back the factor of $M$ provides the solution to the inhomogeneous CGF.

```maxima tags=[]
gaminh : (M/D) * log( (D/s) * (%e^(-s * tau) - 1) * theta + 1);
```

## Inverse Laplace transform


### Homogeneous solution


We proceed by constructing a change of variable to place the inverse Laplace transform integral into a form enabling the identification of a Bessel function solution. Combining the form of the Laplace transform with the form of `gamhom` we have

```maxima tags=[]
schemintgr : x * theta  - a * theta / (b * theta + 1);
```

The following change of variable $\theta \rightarrow \theta'$ will render `shemintgr` a function of $\theta' + \frac{1}{\theta'}$ as expected for a Bessel function

```maxima tags=[]
thetaChV: c * thp - 1/b;
```

```maxima tags=[]
schemintgrthp: expand(subst(thetaChV, theta, schemintgr));
```

Since we are looking for a form that includes a single factor of $\theta' + \frac{1}{\theta'}$ we look for $c x = \frac{a}{b^2 c}$

```maxima tags=[]
cbessel: a^(1/2)*x^(-1/2)/b$
thetaChVab: subst(cbessel, c, thetaChV);
schemintgrelimc: expand(subst(thetaChVab, theta, schemintgr));
```

To provide $\theta$ in terms of $\theta'$ we substitute the values of $a$ and $b$ given by the homogeneous solution `gamhom`. First we note that if we define $a = x_0  e^{s \tau}$ and $b = \frac{D}{s} (e^{s \tau} -1)$

```maxima tags=[]
avalue: x0 * %e^(s*tau);
bvalue: (D/s) * (%e^(s*tau)- 1);
```

that $\mathrm{gamhom} = \frac {a \theta}{b \theta +1}$

```maxima tags=[]
factor(gamhom - avalue*theta/(bvalue*theta +1));
```

Now we substitute values of $a$ and $b$ into the terms of `schemintgrelimc`


```maxima tags=[]
factor(psubst([a = avalue, 
               b = bvalue],
               thetaChVab));
```

```maxima tags=[]
factor(psubst([a = avalue, 
               b = bvalue],
               -a/b - x/b));
```

```maxima tags=[]
factor(psubst([a = avalue, 
               b = bvalue],
               sqrt(a*x)/b));
```

### Inhomogeneous solution

```maxima tags=[]
invinhvarchange : theta = 
    thetap/x - s / (D*(%e^(s * tau) - 1));
```

```maxima tags=[]
diff(part(invinhvarchange,2), thetap, 1);
```

```maxima tags=[]
factor(ratsimp(
  psubst(invinhvarchange, 
         diff(part(invinhvarchange,2), thetap, 1)
         *e^(x*theta)*((D/s) * (%e^(s * tau) - 1) * theta + 1)^(-M/D))
         ));
```

```maxima tags=[]
factor(part(%o116,1,2,2)+part(%o116,1,2,3));
```

```maxima tags=[]
inhdifsol : x^(M/D - 1) * 
    %e^(-s*x/(D * (%e^(s * tau) - 1))) *
    ((D*(%e^(s * tau) - 1)/s)^(-M/D));
```

```maxima tags=[]
inhdifsol : x^(M/D - 1) * 
    %e^(-s*x/(D * (%e^(s * tau) - 1))) *
    (((%e^(s * tau) - 1))^(-M/D));
```

```maxima tags=[]
factor(expand(diffusioneqn (inhdifsol, D, s, M)));
```
