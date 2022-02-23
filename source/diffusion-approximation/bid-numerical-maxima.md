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

# Numerical approximation of the birth-immigration-death process

<!-- <center><font size="+4">Comparing the exact solution and diffusion approximation to the birth-immigration-death process</font></center> -->


## Birth-death


### Exact solution

```maxima tags=[]
apiece(tau, sigma) := -((sigma-1)*(%e^tau-1)) / 
    (sigma*%e^tau+%e^tau+sigma-1);
bpiece(tau, sigma) := (4*sigma^2*%e^tau) / 
    (sigma*%e^tau+%e^tau+sigma-1)^2;
cpiece(tau, sigma) := ((sigma+1)*(%e^tau-1)) /   
    (sigma*%e^tau+%e^tau+sigma-1);
```

```maxima tags=[]
probdistarr(tau, sigma, dist, nmax) := block(
    [a, b, c, d, k],
    a : apiece(tau,sigma),
    b : bpiece(tau,sigma),
    c : cpiece(tau, sigma),
    
    /* Generate the spike at n = 0. */
    dist[0] : a,
    /* Generate the rest of the distribution. */
    for k : 1 thru nmax do dist[k] : b * c^(k-1))$
```

```maxima tags=[]
convolve(a, b, c, nmax) := block(
    [m, k],
    for m : 0 thru nmax do(
        c[m] : 0,
        for k : 0 thru m do
            c[m] : c[m] + a[k] * b[m - k]))$
```

```maxima tags=[]
convolpow(a, m, b, nmax) := block(
    [tmparr, /* Temporary working array. */
    doubler /* Store consecutive doublings of a. */ ,
    n /* Counter for for loops. */ ],
    tmparr : make_array(flonum, nmax + 1),
    doubler : make_array(flonum, nmax + 1),
        
    /* Initialize the output array suitably. */
    if mod(m, 2) = 1 then(
        /* Set b = a when m is odd. */ 
        for n : 0 thru nmax do b[n] : a[n],
        m : m - 1)
        /* Set b to identity element when m is even. */
        else(
            b[0] : 1.0,
            for n : 1 thru nmax do b[n] : 0),
    m : m/2,
    
    /* Initialize the doubler. */
    for n : 0 thru nmax do doubler[n] : a[n],
    
    /*  Compute the result. */
    while m > 0 do(
        /* print("entering loop m= ", m), */
        /* Double the doubler. */
        for n : 0 thru nmax do tmparr[n] : doubler[n],
        convolve(tmparr, tmparr, doubler, nmax),
        /* When m is odd we multiply the doubler into b. */
        if mod(m, 2) = 1 then(
            for n : 0 thru nmax do tmparr[n] : b[n],
            convolve(tmparr, doubler, b, nmax),
            m : m - 1),
        m : m/2))$
```

```maxima tags=[]
exact_dist(tau, sigma, n0, nmax) := block(
    [dist0, distn],
    dist0 : make_array(flonum, nmax + 1),
    distn : make_array(flonum, nmax + 1),
    
    probdistarr(tau, sigma, dist0, nmax),
    convolpow(dist0, n0, distn, nmax),
    listarray(distn))$
```

### Diffusion approximation

```maxima tags=[]
diff_dist (xi, tau, xi0) :=
    (%e^(tau/2)*sqrt(xi0))/((%e^tau-1)*sqrt(xi)) *
    %e^(-(%e^tau*xi0+xi)/(%e^tau-1)) *
    bessel_i(1, 2*%e^(tau/2)*sqrt(xi*xi0)/(%e^tau-1));
```

### Comparing exact and diffusion

```maxima tags=[]
rescaled_exact(N, tau, s, x0, xmax) := block(
    [dist0, distn, nmax],
    nmax : N * xmax,
    dist0 : make_array(flonum, nmax + 1),
    distn : make_array(flonum, nmax + 1),
    
    probdistarr(tau, s/N, dist0, nmax),
    convolpow(dist0, N*x0, distn, nmax),
    makelist([n/N, N*distn[n]], n, 0, nmax))$
```

```maxima tags=[]
plot2d([[discrete, rescaled_exact(1, 1.0, 0.33, 5, 50)],
    [discrete, rescaled_exact(2, 1.0, 0.33, 5, 50)],
    [discrete, rescaled_exact(3, 1.0, 0.33, 5, 50)],
    [discrete, rescaled_exact(4, 1.0, 0.33, 5, 50)],
    0.66 * diff_dist(0.66*x, 1.0, 0.66*5.0)], [x, 0.1, 50.0],
    [legend, "N=1", "N=2", "N=3", "N=4", "diffusion"],
    [title, "tau = 1.0, s = 0.33, x_0 = 5"],
    [pdf_file, "fig/rescaled-exact.pdf"],
    [gnuplot_preamble, "set key right"],
    [gnuplot_pdf_term_command, "set term pdfcairo lw 3 size 17.2 cm, 12.9 cm font 'Latin Modern Roman,20'"]);
```

## Fix first two moments and vary sigma


Numerically calculate the cumulants.

```maxima tags=[]
cumus(dist) := block([len, m1, m2, m3],
    len : length(dist) - 1,
    m1 : sum(j * dist[j + 1], j, 0, len),
    m2 : sum(j^2 * dist[j + 1], j, 0, len),
    m3 : sum(j^3 * dist[j + 1], j, 0, len),
    [m1, m2 - m1^2, m3 - 3*m1*m2 + 2*m1^3])$
```

Analytically calculate cumulants of the BD process.

```maxima tags=[]
kaps(tau, n0, sigma) := [
    n0 * %e^(tau),
    (n0 / sigma) * (%e^(2*tau) - %e^(tau)),
    (n0 / (2 * sigma^2)) * (%e^(2*tau) - %e^(tau)) *
        (3*(%e^(tau) - 1) + (%e^(tau) + 1) * sigma^2)]$
```

Given values tau1 and n01, determine tau and n0 such that the BD distribution with parameters (tau, n0, sigma) has the same values of the first two moments as the BD distribution wiwth parameters (tau1, n01, 1).

```maxima tags=[]
solve((%e^tau-1)/sigma =%e^tau1-1, tau);
```

```maxima tags=[]
psubst([tau=log(sigma*%e^tau1-sigma+1)], n01 * %e^(tau1 - tau));
```

```maxima tags=[]
momequiv (tau1, n01, sigma) :=
	    [log(sigma*%e^tau1-sigma+1), 
	    (n01*%e^tau1)/(sigma*%e^tau1-sigma+1)]$
```

Comparison for n01 = 12.

```maxima tags=[]
momequiv(0.1, 12, 0.2);
momequiv(0.1, 12, 0.8);
```

```maxima tags=[]
kaps(0.0208, 13, 0.2);
kaps(0.0808, 12, 0.8);
```

```maxima tags=[]
plot2d([[discrete, exact_dist(0.0208, 0.2, 13, 30)], 
    [discrete, exact_dist(0.0808, 0.8, 12, 30)]],
    [legend, "tau = 0.2", "tau = 0.8"],
    [pdf_file, "fig/n01-12.pdf"],
    [gnuplot_preamble, "set key right"],
    [gnuplot_pdf_term_command, "set term pdfcairo lw 3 size 17.2 cm, 12.9 cm font 'Latin Modern Roman,20'"]);
```

```maxima tags=[]
cumus(exact_dist(0.0208, 0.2, 13, 30));
cumus(exact_dist(0.0808, 0.8, 12, 30));
```

Comparison for n01 = 24

```maxima tags=[]
momequiv(0.1, 24, 0.2);
momequiv(0.1, 24, 0.8);
```

```maxima tags=[]
kaps(0.0208, 26, 0.2);
kaps(0.0808, 24, 0.8);
```

```maxima tags=[]
cumus(exact_dist(0.0208, 0.2, 26, 40));
cumus(exact_dist(0.0808, 0.8, 24, 40));
```

```maxima tags=[]
plot2d([[discrete, exact_dist(0.0208, 0.2, 26, 40)], 
    [discrete, exact_dist(0.0808, 0.8, 24, 40)]],
    [legend, "tau = 0.2", "tau = 0.8"],
    [pdf_file, "fig/n01-24.pdf"],
    [gnuplot_preamble, "set key right"],
    [gnuplot_pdf_term_command, "set term pdfcairo lw 3 size 17.2 cm, 12.9 cm font 'Latin Modern Roman,20'"]);
```

Comparison for n01 = 48

```maxima tags=[]
momequiv(0.1, 49, 0.2);
momequiv(0.1, 49, 0.8);
```

```maxima tags=[]
kaps(0.0208, 53, 0.2);
kaps(0.0808, 50, 0.8);
```

```maxima tags=[]
cumus(exact_dist(0.0208, 0.2, 53, 80));
cumus(exact_dist(0.0808, 0.8, 50, 80));
```

```maxima tags=[]
plot2d([[discrete, exact_dist(0.0208, 0.2, 53, 80)], 
    [discrete, exact_dist(0.0808, 0.8, 50, 80)]],
    [legend, "tau = 0.2", "tau = 0.8"],
    [pdf_file, "fig/n01-48.pdf"],
    [gnuplot_preamble, "set key right"],
    [gnuplot_pdf_term_command, "set term pdfcairo lw 3 size 17.2 cm, 12.9 cm font 'Latin Modern Roman,20'"]);
```

```maxima
load(implicit_plot);
```

```maxima
implicit_plot([(sigma^2/3) * coth(tau/2) = 0.25,
    (sigma^2/3) * coth(tau/2) = 0.5, 
    (sigma^2/3) * coth(tau/2) = 1.0,
    (sigma^2/3) * coth(tau/2) = 2.0,
    (sigma^2/3) * coth(tau/2) = 4.0,
    (sigma^2/3) * coth(tau/2) = 8.0], [sigma, 0.0, 1], [tau, 0.01, 0.5],
    [legend, false]
);


[(sigma^2/3) * coth(tau/2) = 0.25,
    (sigma^2/3) * coth(tau/2) = 0.5, 
    (sigma^2/3) * coth(tau/2) = 1.0,
    (sigma^2/3) * coth(tau/2) = 2.0,
    (sigma^2/3) * coth(tau/2) = 4.0,
    (sigma^2/3) * coth(tau/2) = 8.0], [sigma, 0.0, 1], [tau, 0.01, 0.5],
    [legend, false]
```

```maxima tags=[]
plot2d ([(sigma^2/3) * coth(tau/2) = 0.25,
         (sigma^2/3) * coth(tau/2) = 0.5, 
         (sigma^2/3) * coth(tau/2) = 1.0,
         (sigma^2/3) * coth(tau/2) = 2.0,
         (sigma^2/3) * coth(tau/2) = 4.0,
         (sigma^2/3) * coth(tau/2) = 8.0],
        [sigma, 0.0, 1], 
        [tau, 0.01, 0.5],
        [color, "#1F654C", "#226F54", "#559972", "#6EAE81", "#87C38F", "#BEDAA5"],
        [legend, false],
        [xlabel, "s"],
        [ylabel, "t"],
        [pdf_file, "fig/contour-plot-01.pdf"],
        [gnuplot_preamble, "
set label '0.25' at 0.3, 0.4 
set label '0.5' at 0.45, 0.37
set label '1.0' at 0.6, 0.3
set label '2.0' at 0.72, 0.21
set label '4.0' at 0.8, 0.14
set label '8.0' at 0.85, 0.08
set object circle at 0.2, 0.02 size 0.005 front fillstyle solid 1.0 fc rgb '#5B5B5B'
set label 'A' at 0.21, 0.02 front tc rgb '#5B5B5B'
set object circle at 0.8, 0.08 size 0.005 fillstyle solid 1.0 fc rgb '#5B5B5B'
set label 'B' at 0.81, 0.08 front tc rgb '#5B5B5B'
set object circle at 0.2, 0.2 size 0.005 fillstyle solid 1.0 fc rgb '#5B5B5B'
set label 'C' at 0.21, 0.2 front tc rgb '#5B5B5B'
set object circle at 0.8, 0.2 size 0.005 fillstyle solid 1.0 fc rgb '#5B5B5B'
set label 'D' at 0.81, 0.2 front tc rgb '#5B5B5B'"],
        [gnuplot_pdf_term_command, "set term pdfcairo enhanced lw 3.5 size 17.2 cm, 12.9 cm font 'Latin Modern Roman,23'"]);
```

## BID


Generate the distribution for the BD process with n_0 = 1 from the geometric series.  The answer goes into the array dist which runs from 0 to nmax.

```maxima tags=[]
BDdist(tau, sigma, dist, nmax) := block(
    [a, b, c, k],
    a : -((sigma-1)*(%e^tau-1)) / 
    (sigma*%e^tau+%e^tau+sigma-1),
    b : (4*sigma^2*%e^tau) / 
    (sigma*%e^tau+%e^tau+sigma-1)^2,
    c : ((sigma+1)*(%e^tau-1)) /   
    (sigma*%e^tau+%e^tau+sigma-1),
    
    /* Generate the spike at n = 0. */
    dist[0] : a,
    /* Generate the rest of the distribution. */
    for k : 1 thru nmax do dist[k] : b * c^(k-1))$
```

Given three arrays a, b, c which run form 0 to nmax, conolve a with b, writing the result into c. 

```maxima tags=[]
convolve(a, b, c, nmax) := block(
    [m, k],

    for m : 0 thru nmax do(
        c[m] : 0,
        for k : 0 thru m do
            c[m] : c[m] + a[k] * b[m - k]))$
```

Compute the repeated convolution of an array a with itself m times, writing the result into an array b.  Both arrays are assumed to run from 0 to nmax.

```maxima tags=[]
convolpow(a, m, b, nmax) := block(
    [tmparr, /* Temporary working array. */
    doubler /* Store consecutive doublings of a. */ ,
    n /* Counter for for loops. */ ],
    tmparr : make_array(flonum, nmax + 1),
    doubler : make_array(flonum, nmax + 1),
        
    /* Initialize the output array suitably. */
    if mod(m, 2) = 1 then(
        /* Set b = a when m is odd. */ 
        for n : 0 thru nmax do b[n] : a[n],
        m : m - 1)
        /* Set b to identity element when m is even. */
        else(
            b[0] : 1.0,
            for n : 1 thru nmax do b[n] : 0),
    m : m/2,
    
    /* Initialize the doubler. */
    for n : 0 thru nmax do doubler[n] : a[n],
    
    /*  Compute the result. */
    while m > 0 do(
        /* print("entering loop m= ", m), */
        /* Double the doubler. */
        for n : 0 thru nmax do tmparr[n] : doubler[n],
        convolve(tmparr, tmparr, doubler, nmax),
        /* When m is odd we multiply the doubler into b. */
        if mod(m, 2) = 1 then(
            for n : 0 thru nmax do tmparr[n] : b[n],
            convolve(tmparr, doubler, b, nmax),
            m : m - 1),
        m : m/2))$
```

Special case of the BID process with n_0 zero but m nonzero.

```maxima tags=[]
Idist(tau, sigma, m, dist, nmax) := block(
    [a, b, c, p, n],

    a : (1 + sigma - (1 - sigma) * %e^(-tau)) / 
         (2 * sigma * %e^(-tau)),
    b : (1 + sigma) * (1 - %e^(-tau)) /
         (1 + sigma - (1 - sigma) * %e^(-tau)),
    c : 2 * m * sigma / (1 + sigma),

    p : a^(-c),
    for n : 0 thru nmax do(
        dist[n] : p,
        p : p * b * (c + n) / (1 + n)))$
```

```maxima

```
