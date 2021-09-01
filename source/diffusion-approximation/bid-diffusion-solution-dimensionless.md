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

# Solving the diffusion approximation to the dimensionless birth-immigration-death process

<!-- <center><font size="+4">Solving the diffusion approximation to the birth-immigration-death process</font></center> -->


## Forward Laplace transform

```maxima tags=[]
diffusioneqn (rho, s, M) := 
    diff (rho, tau) -
    (1/(2*s)) * diff (x * rho, x, 2) +
    diff (x * rho, x, 1) +
    M * diff (rho, x);
```

```maxima tags=[]
cgfeqn (gam, s, M) := 
    diff (gam, tau) +
    ((1/(2*s)) * theta^2 - theta) * diff (gam, theta) -
    M * theta;
```

```maxima tags=[]
factor(diffusioneqn (rho0(x,tau), s, M));
```

## Characteristic equation


Compute the integral for the implicit solution.

```maxima tags=[]
implicitint : logcontract(integrate(1/(thp^2/(2*s) -  thp), thp));
```

```maxima
factor(diff(implicitint, thp));
```

```maxima
logcontract(subst(theta, thp, implicitint) -
    subst(theta0, thp, implicitint));
```

```maxima tags=[]
:lisp (msetq $implicitsol '((MTIMES SIMP) ((MPLUS SIMP) ((MTIMES SIMP) -2 $S ) $THETA) $THETA0 ((MEXPT SIMP) ((MPLUS SIMP) ((MTIMES SIMP) -2 $S) $THETA0) -1) ((MEXPT SIMP) $THETA -1)))
```

```maxima tags=[]
logcontract(subst(theta, thp, implicitint) -
    subst(theta0, thp, implicitint) -
    log(implicitsol));
```

Solve the implicit equation.

```maxima tags=[]
thetasol : solve ([implicitsol = %e^(tau)], [theta])[1];
theta0sol : solve ([implicitsol = %e^(tau)], [theta0])[1];
```

Check the solution.

```maxima tags=[]
subst(0, tau,  thetasol);
```

```maxima tags=[]
chareqn(theta) := diff(theta, tau) - (1/(2*s)) * theta^2 + theta;
```

```maxima tags=[]
factor(chareqn(part(thetasol, 2)));
```

## Homogeneous solution


Re-expressing $u$ in terms of $\tau$ and including the initial condition $P(x,0) = \delta(x-x_0) \Rightarrow \Gamma_{\mathrm{hom}}(\theta, 0)=x_0 \, \theta$, where $\Rightarrow$ indicates taking the logarithm of the Laplace transform, yields

```maxima tags=[]
gamhom : x0 * part(theta0sol, 2);
```

Noting $M=0$ in the homogeneous case

```maxima tags=[]
factor(cgfeqn (gamhom, s, 0));
```

## Inhomogeneous solutions


Find the characteristic which passes through a given point (theta_f, tau_f).

```maxima tags=[]
th0sol : solve([thetaf = 
        subst(tauf, tau, part(thetasol, 2))], [theta0])[1];
```

```maxima tags=[]
thfsol : factor(psubst(th0sol, thetasol));
```

```maxima tags=[]
subst(tauf, tau, thfsol);
```

Integrate along that characteristic.

```maxima
ginhint : integrate(part(thfsol, 2), tau);
```

```maxima tags=[]
ginhint : integrate((thetaf*s*%e^(uf*s))/
    (s*%e^(s*u)-D*thetaf*%e^(s*u)+D*thetaf*%e^(uf*s)), u);
```

Check the solution.

```maxima
factor(cgfeqn (psubst([tauf = tau, thetaf = theta], gaminh), s, M));
```

```maxima tags=[]
gaminh : M *logcontract(ev(factor(subst(tauf, tau, ginhint) - 
            subst(0, tau, ginhint)), logexpand=super));
```

## Inverse Laplace transform


### Homogeneous solution


We proceed by constructing a change of variable to place the inverse Laplace transform integral into a form enabling the identification of a Bessel function solution. Combining the form of the Laplace transform with the form of `gamhom` we have

```maxima tags=[]
schemintgr : x * theta  - a * theta / (b * theta + 1);
```

The following change of variable $\theta \rightarrow \theta'$ will render `shemintgr` a function of $\theta' + \frac{1}{\theta'}$ as expected for a Bessel function

```maxima tags=[]
expand(subst(c * thp - 1/b, theta, schemintgr));
```

Since we are looking for a form that includes a single factor of $\theta' + \frac{1}{\theta'}$ we look for $c x = \frac{a}{b^2 c}$

```maxima tags=[]
expand(subst((a^(1/2) * x^(-1/2) * thp - 1)/b, theta, 
            schemintgr));
```

To provide $\theta$ in terms of $\theta'$ we substitute the values of $a$ and $b$ given by the homogeneous solution `gamhom`. First we note that if we define $a = x_0  e^{\tau}$ and $b = \frac{1}{2s} (e^{s \tau} -1)$

```maxima tags=[]
aquant : x0 * %e^(tau);
bquant : (%e^(tau) - 1) / (2 * s);
```

```maxima tags=[]
homthetasub: factor(psubst([a = aquant, b = bquant],
(a^(1/2) * x^(-1/2) * thp - 1)/b));
```

```maxima
thetarecipcoeff : factor(psubst([a = aquant, b = bquant],
(a^(1/2) * x^(-1/2) /b)));
```

```maxima
thetaexpconst : factor(psubst([a = aquant, b = bquant],
-a/b - x/b));
```

### Inhomogeneous solution

```maxima tags=[]
invinhvarchange : theta = 
    thetap/x - 2 * s / (%e^(tau) - 1);
```

```maxima
factor(psubst(invinhvarchange, 
  (%e^(tau) - 1) * theta / (2 * s) + 1));
```

```maxima tags=[]
inhdifsol : x^(2*M*s - 1) * 
    %e^(-2*s*x/(%e^(tau) - 1)) /
    ((%e^(tau) - 1)^(2*M*s));
```

```maxima tags=[]
factor(diffusioneqn (inhdifsol, s, M));
```
