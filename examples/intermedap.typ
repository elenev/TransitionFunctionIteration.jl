// Document settings
#set page(paper: "us-letter", margin: 1.25in)
#set text(font: "Linux Libertine", size: 11pt)
#set par(justify: true)
// #set heading(numbering: "1.")

// Title Block
#align(center)[
  #text(16pt, weight: "bold")[Intermediary Asset Pricing Model] \
  #v(0.5em)
  *Vadim Elenev* \
  #datetime.today().display()
]

#v(1.5em)


// This is a typst document. Math is formatted using typst, NOT LATEX.

Time is discrete and infinite. There is a continuum of infinitely-lived households, each of which maximizes expected discounted utility from consumption. The representative household has preferences

$
E_0 sum_(t=0)^infinity beta^t u(c_t) \
u(c_t) = c_t^(1-gamma)/(1-gamma)
$

Households supply their unit labor endowment inelastically to firms, which produce output using a Cobb-Douglas production function with capital and labor. The technology growth rate follows an AR(1) process in logs.

Firms are owned by financial intermediaries. Intermediaries optimize their capital structure between equity and one-period debt, subject to a leverage constraint and costly equity issuance. Intermediaries maximize the expected discounted value of dividends/

Households own the equity of financial intermediaries, receiving dividends and buying their bonds.

== Household Problem

Let $a_t$ denote household wealth at current prices. We will use $cal(Z)_t^H$ to denote state variables exogenous to the household -- aggregate capital stock $K_t$, aggregate intermediary net worth $N_t$, the productivity trend $Z_t$, and its growth rate $g_t$. To simplify notation, we drop the explicit dependence of prices and dividends on these state variables.

 The household's problem can be written as:

$
  V^H (a_t, cal(Z)_t^H ) = max_(c_t, s_t, b_t^h) u(c_t) + beta E_t V^H (W^H_(t+1), cal(Z)_(t+1))
$

subject to:

$
  c_t + s_t Q_t + p_t b_t^h &= a_t + w_t \
  a_t &= (Q_t + D_t) s_(t-1) + b_(t-1)^h
$

=== Stationarized Problem

Let $hat(c)_t = c_t / Z_t$, $hat(a)_t = a_t / Z_t$, $hat(Q)_t = Q_t / Z_t$, $hat(D)_t = D_t / Z_t$, $hat(b)_(t-1)^h = b_(t-1)^h / Z_t$, $hat(hat(b))_t^h = b_t^h / Z_t$, and $hat(w)_t = w_t / Z_t$, and $hat(V)^H_t = V^H_t / Z_t^(1-gamma)$. Then the stationarized problem is:

$
  hat(V)^H (hat(a)_t, cal(hat(Z))_t^H ) = max_(hat(c)_t, s_t, b_t^h)  u(hat(c)_t) + beta  E_t [ hat(V)^H (hat(a)_(t+1), cal(hat(Z))_(t+1)) exp((1-gamma) g_(t+1)) ]
$
subject to:
$
  hat(c)_t + s_t hat(Q)_t + p_t hat(hat(b))_t^h &= hat(a)_t + hat(w)_t \
  hat(a)_t  &= (hat(Q)_t + hat(D)_t) s_(t-1) + exp(-g_t) hat(hat(b))_(t-1)^h
$
where $cal(hat(Z))_t^H = (hat(K)_(t-1), hat(N)_t, hat(g)_t)$.


=== Optimality Conditions

The first-order conditions for the household's problem are:

$
  u'(hat(c)_t) p_t &= beta E_t [ exp((1-gamma) g_(t+1)) u'(hat(c)_(t+1)) exp(-g_(t+1))  ] \
  u'(hat(c)_t) hat(Q)_t &= beta E_t [ exp((1-gamma) g_(t+1)) u'(hat(c)_(t+1)) (hat(Q)_(t+1) + hat(D)_(t+1)) ] \
$


The household optimality conditions define a stationarized stochastic discount factor 

$
hat(M)_(t,t+1) = beta exp((1-gamma) g_(t+1)) ( hat(c)_(t+1) / hat(c)_t )^(-gamma) 
$

The non-stationarized optimality condition, in contrast, imply the following SDF:

$
M_(t,t+1) = beta ( c_(t+1) / c_t )^(-gamma)
$

As a result, 

$
M_(t,t+1) =  beta ( hat(c)_(t+1) / hat(c)_t) ^(-gamma) exp(-gamma g_(t+1)) = 
hat(M)_(t,t+1) exp(-g_(t+1))
$


== Firm Problem

Let the static part of the firm problem be defined as:
$
  Pi(Z_t, K_(t-1)) = max_(L_t, I_t) ( K_t^alpha (Z_t L_t)^(1-alpha) - w_t L_t - I_t - phi /2 (I_t / K_t - tilde(delta) )^2 + q_t I_t )
$
where labor-augmenting productivity trend is given by:

$
z_t equiv log Z_t = z_(t-1) + g_t \
g_t = (1 - rho) macron(g) + rho g_(t-1) + sigma_(g,epsilon) epsilon_t^g
$

Define the intermediary's net worth as:

$
  N_t = Pi(Z_t, K_(t-1)) + (1 - delta) q_t K_(t-1) - B_(t-1)
$

The intermediary pays a dividend of $D_t = eta N_t - X_t$, where $eta$ is the target payout ratio and $X_t$ represents equity issuance, which incurs a cost $xi / 2 (X_t / N_t)^2 N_t $.

Use $cal(Z)_t^F$ to denote the state variables exogenous to the intermediary. The intermediary's problem is then:

$
  V^F (N_t, cal(Z)_t^F ) = max_(K_t, X_t, B_t) eta N_t - X_t + E_t [ M_(t,t+1) V^F (N_(t+1), cal(Z)_(t+1)) ]
$

subject to the budget constraint
$
  q_t K_t  = (1-eta) N_t + p_t B_t + X_t - xi / 2 (X_t / N_t)^2 N_t
$

and the leverage constraint:
$
  p_t B_t <= theta q_t K_t
$

=== Static Problem

Let $l_t = ( Z_t L_t ) / K_(t-1)$ and $i_t = I_t / K_(t-1)$. Write the stationarized wage $hat(w) = w_t / Z_t$. Then by dividing $Pi$ through by $K_(t-1)$ shows that the static problem is linear in capital:

$
  Pi(Z_t, K_(t-1)) = 
  K_(t-1) 
  underbrace( 
  max_(l_t, i_t) ( l_t^(1-alpha) - hat(w)_t l_t - i_t - phi /2 (i_t - tilde(delta) )^2 + q_t i_t )
  , r_t^k  )
$

Optimality conditions are:

$
  (1-alpha) l_t^(-alpha) &= hat(w)_t \
  1 + phi (i_t - tilde(delta)) &= q_t
$

Plugging in the wage optimality condition into $r^k_t$, we get:

$
  r^k_t equiv alpha l_t^(1-alpha) - phi/2 (i_t - tilde(delta))^2 + (q_t - 1) i_t
$

=== Stationarized Problem

Let $hat(K)_(t-1) = K_(t-1) / Z_t$ and $hat(hat(K))_t = K_t / Z_t$ denote stationarized start-of-period and end-of-period capital,respectively. The same convention applies to $hat(B)_(t-1) = B_(t-1) / Z_t$ and $hat(hat(B))_t = B_t / Z_t$. Let $hat(N)_t = N_t / Z_t$ be stationarized net worth and $hat(V)^F_t = V^F_t / Z_t$.

Plugging in the stationarized SDF $hat(M)_(t,t+1)$, we get the stationarized problem is:

$
  V^F (hat(N)_t, cal(hat(Z))_t^F ) = max_(hat(hat(K))_t,hat(X)_t, hat(hat(B))_t) eta hat(N)_t - hat(X)_t + E_t [hat(M)_(t,t+1)  hat(V)^F ( hat(N)_(t+1), cal(hat(Z))_(t+1)) ]
$
subject to:

$
  q_t hat(hat(K))_t &= p_t hat(hat(B))_t + (1-eta) hat(N)_t + hat(X)_t - xi / 2 (hat(X)_t / hat(N)_t)^2 hat(N)_t \
  hat(N)_(t+1) &= exp(-g_(t+1)) [R^K hat(hat(K))_t (r_t^k + (1-delta)) -  hat(hat(B))_t ] \
  hat(K)_t &= hat(hat(K))_t exp(-g_(t+1)) \
  p_t hat(hat(B))_t &<= theta q_t hat(hat(K))_t
$

where $cal(hat(Z))_t^F = (hat(K)_(t-1), hat(N)_t, hat(g)_t)$ denotes the state variables exogenous to the intermediary. Note that here $hat(K)_(t-1)$ and $hat(N)_t$ denote _aggregate_ variables, which affect equilibrium prices, but which the intermediary cannot optimize over.

Plugging in $hat(K)_(t-1) = K_(t-1) / Z_t$ into the definition for $l_t$, we can rewrite the return on capital as:

$
  r^k_t equiv alpha (L / hat(K)_(t-1))^(1-alpha) - phi/2 (i_t - tilde(delta))^2 + (q_t - 1) i_t
$

// === Optimality Conditions

// Assign $hat(mu)_t$ as the multiplier on the budget constraint and $hat(lambda)_t$ as the multiplier on the leverage constraint. The first-order conditions for the intermediary's problem are:

// $
//   hat(mu)_t q_t &= hat(lambda)_t theta + E_t [ hat(M)_(t,t+1) hat(mu)_(t+1) exp(-g_(t+1)) (r_(t+1)^k + (1-delta)) ] \
//   hat(mu)_t p_t &= hat(lambda)_t + E_t [ hat(M)_(t,t+1) hat(mu)_(t+1) exp(-g_(t+1)) ] \
//   1 &= hat(mu)_t (1 - xi hat(X)_t / hat(N)_t)
// $


=== In Units of Net Worth 
Define $k_t = hat(hat(K))_t / hat(N)_t$ and $b_t = hat(hat(B))_t / hat(N)_t$ as well as $v^F = V^F_t / hat(N)_t$. Then the intermediary's problem can be rewritten as:

$
  hat(N)_t v^F ( cal(Z)_t^F ) = max_(k_t, x_t, b_t) eta hat(N)_t - x_t hat(N)_t + E_t [ hat(M)_(t,t+1) v^F (cal(Z)_(t+1)^F) hat(N)_(t+1) ]
$

And net worth can be written as:

$
  hat(N)_(t+1) = exp(-g_(t+1))hat(N)_t [k_t (r_t^k + (1-delta)) - b_t] 
$

Dividing through by $hat(N)_t$ and defining $n_(t+1) = hat(N)_(t+1) / hat(N)_t$, we get:

$
  v^F (cal(Z)_t^F) = max_(k_t, x_t, b_t) eta - x_t + E_t [ hat(M)_(t,t+1) v^F (cal(Z)_(t+1)^F) n_(t+1) ]
$

subject to:

$
  q_t k_t  &= 1 - eta + p_t b_t + x_t - xi / 2 x_t^2 \
  n_t &= exp(-g_t) [k_(t-1) (r_t^k + (1-delta)) - b_(t-1)] \
  p_t b_t &<= theta q_t k_t
$

=== Optimality Conditions

Assign $hat(mu)_t$ as the multiplier on the budget constraint and $hat(lambda)_t$ as the multiplier on the leverage constraint. The first-order conditions for the intermediary's problem are:

$
  hat(mu)_t q_t &= hat(lambda)_t theta + E_t [ hat(M)_(t,t+1) v^F_(t+1) exp(-g_(t+1)) (r_(t+1)^k + (1-delta)) ] \
  hat(mu)_t p_t &= hat(lambda)_t + E_t [ hat(M)_(t,t+1) v^F_(t+1) exp(-g_(t+1)) ] \
  1 &= hat(mu)_t (1 - xi x_t) 
$

Finally, $v^F (cal(Z)_t^F)$ is given by the envelope condition of the problem in levels:

$
  (partial hat(V)^F_t) / (partial hat(N)_t) &= v^F (cal(Z)_t^F) = eta - x_t + hat(mu)_t (1 - eta - xi/2 x_t^2)
$

== Aggregation and Equilibrium

There is a measure 1 of identical households, so household aggregates are straightforward: $hat(C) = integral hat(c)_t = hat(c)_t$, $hat(A) = integral hat(a)_t = hat(a)_t$, $L = 1$, $hat(hat(B))^H_t = integral hat(hat(b))_t^h = hat(hat(b))_t^h$, and $S_t = integral s_t = s_t$

The following markets must clear in equilibrium:

+ Labor: $l_t hat(K)_(t-1) exp(-g_(t-1))  = L$
+ Capital: $k_t hat(N)_t = hat(hat(K))_t = (1-delta + i_t) hat(K)_(t-1) $
+ Bonds: $b_t hat(N)_t = hat(hat(B))^H_t$
+ Shares: $S_t = 1$
+ Goods: $hat(C)_t + (i_t + phi/2 i_t^2) hat(K)_(t-1) + xi/2 x_t^2 hat(N)_t = hat(K)_(t-1)^alpha L^(1-alpha)$

Since the value of the firm is $hat(N)_t hat(V)^F_t$, the value per share is $hat(Q)_t$, and the number of shares is $S_t = 1$, we have that $hat(Q) = hat(N)_t hat(V)^F_t$ and $hat(D)_t = hat(N)_t (eta + x_t)$.

== Deterministic Steady State

At steady state, it must be that $hat(c)_t = hat(c)$, $hat(K)_t = hat(K)$ for all $t$. Imposing the constant law of motion for capital yields investment:

$
  hat(K) = hat(hat(K)) exp(-g) &= (1-delta + i) hat(K) exp(-g) \
  i &= delta + exp(g) - 1 equiv tilde(delta)
$
This then implies that $q = 1$ from the investment optimality condition.

In steady state, as long as the optimal dividend is different from $eta$, intermediaries prefer to lever up to the constraint to avoid equity issuance costs, so it must be that:

$
  p hat(hat(B)) = theta hat(hat(K)), quad tilde(hat(lambda)) equiv hat(lambda)/hat(mu) > 0
$

We can then write $hat(hat(B)) = (theta / p) hat(hat(K))$, where the household optimality condition for bonds implies that $p = beta exp(-gamma g)$.

=== Equilibrium as a function of equity issuance $x$

We now express each equilibrium object as a function of equity issuance $x$, with the goal of finding the $x$ that satisfies all conditions.

*Step 1: $v^F/hat(mu)$ as a function of $x$.*
The equity FOC and envelope condition give:

$
  hat(mu)(x) &= 1/(1 - xi x) \
  v^F/hat(mu) &= (eta - x)(1 - xi x) + (1 - eta - xi/2 x^2)
$

Define $cal(v)(x) = (v^F (x))/(hat(mu)(x))$, i.e., the wedge between household and intermediary SDFs as a function of equity issuance. At $x = 0$, $cal(v)(0) = 1$. 

*Step 2: Normalized leverage multiplier $tilde(hat(lambda))(x)$.*
Define $tilde(hat(lambda)) = hat(lambda)/hat(mu)$ to write the Euler equations in standard asset-pricing form. Dividing the household bond Euler equation by $hat(mu)$:

$
  p = beta exp(-gamma g)
$

Plugging in for $beta exp(-gamma g)$ in the intermediary bond Euler equation and rearranging gives:

$
  tilde(hat(lambda))(x) = p( 1 - cal(v)(x) )
$

*Step 3: Return on capital $r^k (x)$.*

// $
//   1 = q = tilde(hat(lambda)) theta + p cal(v)(x) (r^k (x) + 1 - delta)
// $

Dividing the capital Euler equation by $hat(mu)$ and using $q = 1$:

$
  1 - tilde(hat(lambda))(x) theta = p cal(v)(x) (r^k (x) + 1 - delta)
$

Solving for $r^k$:

$
  r^k (x) = (1 - tilde(hat(lambda))(x) theta) / (p cal(v)(x)) - (1 - delta)
$

*Step 4: Capital-to-net-worth ratio $k(x)$.*
Setting $n = 1$ in the net worth evolution equation with $b = theta k / p$:

$
  // 1 &= ( (r^k (x) + 1 - delta) k(x) - theta k(x) / p )/G \
  // G &= [(r^k (x) + 1 - delta) - theta / p] k(x) \
  k(x) &= G / (r^k (x) + 1 - delta - theta/p)
$

=== Nonlinear equation for $x$

// The steady-state budget constraint is:

// $
//   k &= (1 - eta) + x - xi/2 x^2 + p b \
//   k (1 - theta) &= 1 - eta + x - xi/2 x^2
// $

Substituting $k(x)$ into the intermediary's budget constraint (with $q = 1$, $p b = theta k$):

$
  (1 - theta) k = 1 - eta + x - xi/2 x^2
$

This yields a single nonlinear equation in $x$:

$
  f(x) equiv (1 - theta) k(x) - (1 - eta + x - xi/2 x^2) = 0
$

=== Remaining steady-state values

Given $x^*$ solving $f(x^*) = 0$, all remaining objects follow directly. From the definition of $r^k$:

$
  r^k = alpha (L / hat(K))^(1-alpha) - phi/2(i - tilde(delta))^2 + (q-1) i
$

With $i = tilde(delta)$ and $q = 1$, this simplifies to just the MPK portion. At $L = 1$,

$
  r^k = alpha (L / hat(K))^(1-alpha) \
  hat(K) = ( r^k / alpha )^(1/(alpha-1))
$

From the capital law of motion, we solve for $hat(hat(K)) = hat(K) G$.

Since $k = hat(hat(K)) / hat(N)$, we can recover net worth:

$
  hat(N) &= hat(hat(K)) / k(x^*)
$

Finally, use $B = theta k / p$ and $hat(w) = (1-alpha) hat(K)^alpha$ to get consumption:
$
//  hat(C) &= (eta - x^*) hat(N) + (1-alpha) hat(K)^alpha + hat(B) - p hat(hat(B)) \
  hat(C) &= (eta - x^*) hat(N) + (1-alpha) hat(K)^alpha + theta hat(K) (1 / p - G) \
$

== With Epstein-Zin Preferences

Now, let's let household preferences be given by:

$
  V_t = [ (1-beta) c_t^(1-1/psi) + beta ( E_t V_(t+1)^(1-gamma) )^((1-1/psi)/(1-gamma)) ]^(1/(1-1/psi))
$

where $psi$ is the intertemporal elasticity of substitution and $gamma$ is the coefficient of relative risk aversion.

The household problem is now written as:

$
  V^H (a_t, cal(Z)_t^H ) = max_(c_t, s_t, b_t^h) [ (1-beta) c_t^(1-1/psi) + beta E_t [ ( V^H (a_(t+1),cal(Z)_(t+1)^H))^(1-gamma)]^((1-1/psi)/(1-gamma)) ]^(1/(1-1/psi))
$
s.t.
$
  c_t + s_t Q_t + p_t b_t^h &= a_t + w_t \
  a_t &= (Q_t + D_t) s_(t-1) + b_(t-1)^h
$

The stationarized value function then is:

$
  Z_t hat(V)^H (hat(a)_t, cal(hat(Z))_t^H ) = max_(hat(c)_t, s_t, b_t^h) [ (1-beta) Z_t^(1-1/psi) hat(c)_t^(1-1/psi) + beta ( E_t [Z_(t+1)^(1-gamma) ( V^H (a_(t+1),cal(Z)_(t+1)^H))^(1-gamma)] )^((1-1/psi)/(1-gamma)) ]^(1/(1-1/psi))
$

Writing $Z_(t+1) = exp(g_(t+1)) Z_t$ and pulling $Z_t$ out of the expectation, we can rewrite the stationarized problem as:

$
 hat(V)^H (hat(a)_t, cal(hat(Z))_t^H ) = max_(hat(c)_t, s_t, b_t^h) [ (1-beta) hat(c)_t^(1-1/psi) + beta ( E_t [exp((1-gamma) g_(t+1))) ( V^H (a_(t+1),cal(Z)_(t+1)^H))^(1-gamma))] )^((1-1/psi)/(1-gamma)) ]^(1/(1-1/psi))
$

#let CE = math.op(math.italic("CE"))
Let $hat(CE)_t (hat(X)_(t+1)) = E_t [exp((1-gamma) g_(t+1))) hat(X)_(t+1)^(1-gamma)]^(1/(1-gamma))$ denote the certainty equivalent over stationarized next-period quantities.

Then the stationarized problem can be rewritten as:

$
  hat(V)^H (hat(a)_t, cal(hat(Z))_t^H ) = max_(hat(c)_t, s_t, b_t^h) [ (1-beta) hat(c)_t^(1-1/psi) + beta hat(CE)_t ( 
   V^H (a_(t+1),cal(Z)_(t+1)^H))^(1-1/psi) ]^(1/(1-1/psi))
$
s.t.
$
  hat(c)_t + s_t hat(Q)_t + p_t hat(hat(b))_t^h &= hat(a)_t + hat(w)_t \
  hat(a)_t  &= (hat(Q)_t + hat(D)_t) s_(t-1) + exp(-g_t) hat(hat(b))_(t-1)^h
$

To faciliate the FOC derivations, recall that for the CES aggregator

$
  A(x,y) = (a x^rho + b y^rho)^(1/rho)
$
we have that 

$
(partial A) / (partial x) &= a (a x^rho + b y^rho)^(1/rho - 1) x^(rho - 1) = a rho A^(1-rho) x^(rho - 1) \
(partial A) / (partial y) &= b A^(1-rho) y^(rho - 1)
$ 

Let $hat(theta.alt_t)$ denote the multiplier on the time $t$ budget constraint. The FOCs are:

$
  hat(theta.alt_t) &= (1-beta) hat(V)_t^(1/psi) hat(c)_t^(-1/psi) \
  hat(theta.alt_t) p_t (x) &= E_t [ 
  (partial V^H_t)/(partial CE_t) 
  (partial CE_t)/(partial V^H_(t+1)) 
  (partial V^H_(t+1))/(partial hat(a)_(t+1)) 
  (partial hat(a)_(t+1))/(partial x) 
  ] \
$
where the second expression is a generic first order condition for asset $x$ with a price $p_t (x)$ and payoff $(partial hat(a)_(t+1) ) / (partial x)$.

The partial derivative of the value function with respect to the certainty equivalent is given by:
$
  (partial V^H_t)/(partial CE_t)  = beta hat(V)_t^(1/psi) hat(CE)_t^(-1/psi)
$

The envelope condition states that
$
  (partial hat(V)^H_(t+1)) / (partial hat(a)_(t+1)) = hat(theta.alt)_(t+1) = (1-beta) hat(V)_(t+1)^(-1/psi) hat(c)_(t+1)^(-1/psi) 
$

From the definition of the certainty equivalent, we get that:
$
  (partial CE_t)/(partial V^H_(t+1)) = exp((1-gamma) g_(t+1)) CE_t^(gamma) ( hat(V)^H_(t+1))^(-gamma)
$

Dividing the two asset FOCs by the consumption FOCs yields Euler Equation:

$
  p_t (x) = (
  E_t [ 
  beta hat(V)_t^(1/psi) hat(CE)_t^(-1/psi)

  (1-beta) hat(V)_(t+1)^(-1/psi) hat(c)_(t+1)^(-1/psi) 
   
  exp((1-gamma) g_(t+1)) CE_t^(gamma) ( hat(V)^H_(t+1))^(-gamma) 
  
  (partial hat(a)_(t+1))/(partial x) 
  ] ) / 
  ( (1-beta) hat(V)_t^(1/psi) hat(c)_t^(-1/psi) )
$

Collecting terms,

$
  p_t (x) =
  E_t [ 
  beta exp((1-gamma) g_(t+1))
  ( hat(c)_(t+1) / hat(c)_t )^(-1/psi)
  ( ( hat(CE)_t ) / ( hat(V)_(t+1)^H ) )^(gamma-1/psi)

  (partial hat(a)_(t+1))/(partial x) 
  ] 
$

This equation defines the stationarized SDF as:
$
  hat(M)_(t,t+1) = beta exp((1-gamma) g_(t+1))
  ( hat(c)_(t+1) / hat(c)_t )^(-1/psi)
  ( ( hat(CE)_t ) / ( hat(V)_(t+1)^H ) )^(gamma-1/psi)
$

As expected, at $gamma = 1/psi$, the SDF collapses to the standard power utility form.

The Euler Equations for bonds and stocks are then given by plugging in the appropriate payoffs for $partial hat(a)_(t+1) / partial x$.

$
  p_t &= E_t [ hat(M)_(t,t+1) exp(-g_(t+1)) ] \
  hat(Q)_t &= E_t [ hat(M)_(t,t+1) (hat(Q)_(t+1) + hat(D)_(t+1)) ]
$

=== Steady State SDF

In deterministic steady state, $hat(c)_t = hat(c)_{t+1}$ and 
$
hat(CE)_t (hat(V)^H) = [ exp((1-gamma) g) (hat(V)^H)^(1-gamma)]^(1/(1-gamma))
=
exp(g) hat(V)^H
$. The SDF then is:

$
  hat(M) = beta exp((1-gamma) g) exp(g (gamma - 1/psi) ) = beta exp((1-1/psi) g) 
$

The Euler equations for bonds are then:

$
  p = beta exp((1-1/psi) g) exp(-g) = beta exp(-g/psi)
$

The steady-state value function is:

$
  hat(V)^H = [ (1-beta) hat(c)^(1-1/psi) + beta exp((1-1/psi) g) (hat(V)^H)^(1-1/psi) ]^(1/(1-1/psi))
$
Solving for $hat(V)^H$:
$
  // ( hat(V)^H) ^(1-1/psi) = (1-beta) hat(c)^(1-1/psi) + beta exp((1-1/psi) g) (hat(V)^H)^(1-1/psi) \

  // (hat(V)^H)^(1-1/psi) (1 - beta exp((1-1/psi) g)) = (1-beta) hat(c)^(1-1/psi) \

  // (hat(V)^H)^(1-1/psi) 
  // = ( (1-beta) hat(c)^(1-1/psi) ) /
  // (1 - beta exp((1-1/psi) g))  \

  hat(V)^H
  = [ ( 1-beta ) /
  (1 - beta exp((1-1/psi) g)) ]^(1/(1-1/psi)) hat(c) \
$

In the case of $psi = 1$, the value function becomes Cobb-Douglas, which in steady state is:
$
  hat(V)^H = hat(c)^(1-beta) hat(CE)^(beta) \
$
Plugging in the certainty equivalent and re-arranging,
$
  hat(V)^H &= hat(c)^(1-beta) (exp(g) hat(V)^H)^(beta) \

  (hat(V)^H) &=  exp(beta/(1-beta) g) hat(c)
$