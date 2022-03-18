---
jupyter:
  jupytext:
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

```maxima tags=[]
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
