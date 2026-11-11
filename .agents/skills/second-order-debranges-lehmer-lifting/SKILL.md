---
name: second-order-debranges-lehmer-lifting
description: Scientific knowledge skill for Paper XIII on second-order de Branges-Nevanlinna kernels E_{alpha,beta}(z) = F(z) - i alpha F'(z) - beta F''(z). Contains the exact derivative-Wronskian splitting K_{alpha,beta}(w,w) = (alpha/pi) W_tilde_{0,1} + (alpha beta/pi) W_tilde_{1,2}, canonical Wigner phase-space normalization W_{1,2} = 16 int y sinh(2dy) I^{(1)}(y,t) dy with I^{(1)} = (y^2 + (1/4) d_t^2) I, radial 4D-Laplace dipole recursion, (1/40) (ln(t/2pi))^5 diagonal and O(t^{-5/2}) horizon-pair scaling, and the critical zero-theoretic blind spot on Davenport-Heilbronn (beta_coal = 0.04467500, beta_blind = 0.05349534) motivating the open differential threshold beta_star(F).
---

# Second-Order De Branges–Nevanlinna Kernels, Derivative-Wronskian Splitting, Davenport–Heilbronn Blind Spots, and the Differential Threshold $\beta^{\ast}$ (Paper XIII)

## 1. Scope and Critical Zero-Theoretic Guardrails

This skill encapsulates the verified analytical identities, canonical Wigner phase-space normalizations, interval-arithmetic enclosures, and Davenport–Heilbronn counterexample diagnostics of **Paper XIII** (author: Thijs Missiaen, ORCID: [0009-0007-3678-0810](https://orcid.org/0009-0007-3678-0810)).

### Mandatory Guardrails Against Overclaiming

1. **Hermite–Biehler Positivity of $E _ {\alpha,\beta}$ for Fixed $\beta > 0$ Is Strictly Weaker Than RH**: Writing $E _ {\alpha,\beta}(z) = A _ \beta(z) - \mathrm{i}\alpha B(z)$ with $A _ \beta(z) := F(z) - \beta F''(z)$ and $B(z) := F'(z)$, membership $E _ {\alpha,\beta} \in \overline{\mathcal{HB}}$ corresponds to all zeros of $A _ \beta(z)$ being real and interlacing with $F'(z)$. By Jensen (1913) and Pólya (1927), if $F \in \mathcal{LP}$ then $F - \beta F'' \in \mathcal{LP}$ for all $\beta \ge 0$, but the converse is **false** for any fixed $\beta > 0$ because $1 - \beta\partial_z^2$ pulls narrow conjugate off-line pairs onto $\mathbb{R}$ (note also that Ki 2006, *J. Number Theory* **120**, 120–131 proves asymptotic Farmer–Rhoades cosine universality $\lim _ {n\to\infty} A_n\Xi^{(2n)}(C_n z) = \cos z$ as $n \to \infty$, not real-rootedness of $F - \beta F''$ for fixed $\beta$).
2. **Explicit Davenport–Heilbronn Blind Spot (Theorem 5.1)**: For the self-dual Davenport–Heilbronn function $\Xi _ {\mathrm{DH}}(z)$ near its first off-line zero $z _ {\star} = 85.69934848537759 - 0.3085171824566374\mathrm{i}$ and negative local maximum $\tau _ {\mathrm{false}} = 85.74720082646130$ ($\Xi _ {\mathrm{DH}}(\tau _ {\mathrm{false}}) \approx -1.8610\times 10^{-29} < 0$, $\Xi _ {\mathrm{DH}}''(\tau _ {\mathrm{false}}) \approx -3.4788\times 10^{-28} < 0$):
   - At $\tau _ {\mathrm{false}}$, $c_1^{\mathrm{DH}}(\tau _ {\mathrm{false}}) = -6.474216\times 10^{-57} < 0$ whereas $c_1^{(1),\mathrm{DH}}(\tau _ {\mathrm{false}}) = +1.210239\times 10^{-55} > 0$, so $\pi K _ {1,\beta}^{\mathrm{DH}}(\tau _ {\mathrm{false}}, \tau _ {\mathrm{false}}) > 0$ for all $\beta > \beta _ {\mathrm{blind}}^{(1)}(\tau _ {\mathrm{false}}) = +\Xi _ {\mathrm{DH}}(\tau _ {\mathrm{false}})/\Xi _ {\mathrm{DH}}''(\tau _ {\mathrm{false}}) = 0.053495338542548$ (including $\beta = \frac{1}{4}$).
   - Directly at $z _ {\star}$, $\widetilde{W} _ {0,1}^{\mathrm{DH}}(d _ {\star}, \gamma _ {\star}) = 0$ while $\widetilde{W} _ {1,2}^{\mathrm{DH}}(d _ {\star}, \gamma _ {\star}) = +2.011047\times 10^{-55} > 0$, so $\pi K _ {1,\beta}^{\mathrm{DH}}(z _ {\star}, z _ {\star}) > 0$ for **every** $\beta > 0$.
   - Locally, the first off-line zero pair of $\Xi _ {\mathrm{DH}}(z) - \beta\Xi _ {\mathrm{DH}}''(z)$ coalesces onto $\mathbb{R}$ at $\beta _ {\mathrm{coal}}^{(1)} = 0.044675001258358$ ($t _ {\mathrm{coal}} = 85.629148445738$); globally, Rouché's theorem on the Bohr-almost-periodic zeros of $\mathcal{D} _ {\mathrm{DH}}$ in $\sigma > 1$ proves $\beta^{\star}(\Xi _ {\mathrm{DH}}) = +\infty$ unconditionally.

## 2. Canonical Identities and Theorems

### 2.1 Second-Order Kernel Splitting and Cauchy–Pick Decomposition (Theorem 2.1 & Proposition 2.2)

For any real-entire function $F(z)$, real $\alpha, \beta \in \mathbb{R}$, and $z_1 \neq \bar{z} _ 2$:

$$
K _ {\alpha,\beta}(z_1, z_2) := \frac{E _ {\alpha,\beta}(z_1)\overline{E _ {\alpha,\beta}(z_2)} - E _ {\alpha,\beta}^{\sharp}(z_1)\overline{E _ {\alpha,\beta}^{\sharp}(z_2)}}{2\pi\mathrm{i}(z_1 - \bar{z} _ 2)} \equiv \frac{\alpha}{\pi}\frac{F(z_1)\overline{F'(z_2)} - F'(z_1)\overline{F(z_2)}}{z_1 - \bar{z} _ 2} + \frac{\alpha\beta}{\pi}\frac{F'(z_1)\overline{F''(z_2)} - F''(z_1)\overline{F'(z_2)}}{z_1 - \bar{z} _ 2}.
$$

On the diagonal $w = t - \mathrm{i}d$ ($d \neq 0$):

$$
K _ {\alpha,\beta}(t - \mathrm{i}d, t - \mathrm{i}d) \equiv \frac{\alpha}{\pi}\widetilde{W} _ {0,1}(d,t) + \frac{\alpha\beta}{\pi}\widetilde{W} _ {1,2}(d,t), \qquad W _ {k,k+1}(d,t) := \frac{1}{2}\partial_d\left\lvert F^{(k)}(t - \mathrm{i}d)\right\rvert^2, \qquad \widetilde{W} _ {k,k+1}(d,t) := \frac{W _ {k,k+1}(d,t)}{d}.
$$

### 2.2 Critical-Line Lehmer Lifting at Order $d^0$ (Theorem 3.1)

For $F = \Xi$ as $d \to 0^+$:

$$
\pi K _ {\alpha,\beta}(t,t) \equiv \alpha\left(\left(\Xi'(t)\right)^2 - \Xi(t)\Xi''(t)\right) + \alpha\beta\left(\left(\Xi''(t)\right)^2 - \Xi'(t)\Xi'''(t)\right) = \alpha c_1(t) + \alpha\beta c_1^{(1)}(t).
$$

At an extremum $t _ {\mathrm{mid}}$ ($\Xi'(t _ {\mathrm{mid}}) = 0$) between consecutive simple real zeros $t _ \pm = t _ {\mathrm{mid}} \pm \delta$, $-\Xi'(t _ {\mathrm{mid}})\Xi'''(t _ {\mathrm{mid}}) = 0$ vanishes identically, giving $\pi K _ {\alpha,\beta}(t _ {\mathrm{mid}}, t _ {\mathrm{mid}}) = -\alpha\Xi(t _ {\mathrm{mid}})\Xi''(t _ {\mathrm{mid}}) + \alpha\beta(\Xi''(t _ {\mathrm{mid}}))^2 \ge +\alpha\beta(\Xi''(t _ {\mathrm{mid}}))^2 > 0$.

### 2.3 Canonical Wigner Normalization, Radial 4D-Laplace PDE, and Saddle Scaling (Theorems 4.1 & 4.2)

Under the canonical series normalization $\Xi(z) = 8\int_0^{\infty} \Phi(u)\cos(2zu)\mathrm{d}u$ and $I(y,t) = 2\int_0^{\infty} \Phi(\frac{y+r}{2})\Phi(\frac{y-r}{2})\cos(2tr)\mathrm{d}r$:

$$
W _ {0,1}(d,t) = 16\int_0^{\infty} y\sinh(2dy)I(y,t)\mathrm{d}y, \qquad W _ {1,2}(d,t) = 16\int_0^{\infty} y\sinh(2dy)I^{(1)}(y,t)\mathrm{d}y, \qquad I^{(1)}(y,t) := \left(y^2 + \frac{1}{4}\partial_t^2\right)I(y,t).
$$

The parity-reduced Wronskians satisfy

$$
\widetilde{W} _ {1,2}(d,t) \equiv \frac{1}{4}\left(\frac{\partial^2\widetilde{W} _ {0,1}}{\partial d^2} + \frac{2}{d}\frac{\partial\widetilde{W} _ {0,1}}{\partial d} + \frac{\partial^2\widetilde{W} _ {0,1}}{\partial t^2}\right), \qquad c _ {2m-1}^{(1)}(t) \equiv \frac{1}{4}\left(2m(2m+1)c _ {2m+1}(t) + c _ {2m-1}''(t)\right) \qquad (\forall m \ge 1).
$$

In the Riemann–Siegel main sum up to $N(t) = \lfloor\sqrt{t/2\pi}\rfloor$, the diagonal grows as $\mathcal{D} _ {1,2}(0,t) = 4\sum _ {n=1}^{N(t)} n^{-1}\omega_n(t)^4 \sim \frac{1}{40}(\ln\frac{t}{2\pi})^5$, and the exact $(m,n)$ cross-term $C _ {m,n}^{(1,2)}(t) = \frac{2\omega_m\omega_n}{\sqrt{mn}}[(\omega_m+\omega_n)^2\cos(\phi_m-\phi_n) - (\omega_m-\omega_n)^2\cos(\phi_m+\phi_n)]$ gives $\lvert C _ {N-1,N}^{(1,2)}\rvert = \mathcal{O}(N(t)^{-5}) = \mathcal{O}(t^{-5/2})$ at the horizon pair $(N-1, N)$ (without controlling low-index interior pairs).

### 2.4 Differential Reality Threshold $\beta^{\star}(F)$ (Definition 5.3)

$$
\beta^{\star}(F) := \inf\left\lbrace\beta \in \mathbb{R} : F(z) - \beta F''(z) \text{ has only real zeros}\right\rbrace, \qquad \mathrm{RH} \iff \beta^{\star}(\Xi) \le 0, \qquad \beta^{\star}(\Xi _ {\mathrm{DH}}) = +\infty.
$$

## 3. Reproducible Verification Workflow

```bash
python scripts/verify_paper13_second_order_debranges.py
```
