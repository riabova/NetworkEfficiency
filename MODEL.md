# Bilevel DEA-based Column Generation Method

The problem is defined over a network $G = (V, E)$. The master problem is presented by an arbitrary problem defined over the set of all possible (connected) subnetworks, $\mathcal{K}$. The subnetwork-generating subproblem is given by a bilevel problem that defines a subnetwork and features an embedded DEA model to ensure that a generated network has a negative reduced cost that depends on its DEA score.

We define the master problem as follows:

**Minimize:**

```math
    \min c^\intercal x \quad (1a)
```

**Subject to:**

```math
    Ax \geqslant b \quad (1b)
```
```math
    x \geqslant 0 \quad (1c)
```

The DEA-scores $c_k$ are used as the cost coefficients in the objective (1a) that minimizes the total score of the selected subnetworks. The constraints $Ax \geqslant b$ are application-specific and are given in 
a general form in this section. Examples of these constraints could be vertex cover constraints introduced in the application example, cardinality constraints, or packing constraints.

The dual of this master problem is:

**Maximize:**

```math
    \max b^\intercal \pi \quad (2a)
```

**Subject to:**

```math
    A^\intercal \pi \leqslant c \quad (2b)
```
```math
    \pi \geqslant 0 \quad (2c)
```

In the column generation framework, we start with a subset of subnetworks, $\mathcal{K}'\in \mathcal{K}$. A surnetwork here is defined as a subset of connected nodes and the respective arcs connecting them. Note that any single node can also be viewe as a triviall subnetwork of size 1. That will get us the restricted master problem and the corresponding restricted dual problem. Our goal then is to generate a new subnetwork $k'$ with the negative reduced cost $c_{k'} - \sum_{i \in \mathcal{D}}a_{ik'}\pi_i$. As was mentioned before, $c_k$'s are the DEA scores of the subnetworks. Replacing the DEA score of the network being constructed with the fractional DEA model, we get the following subproblem:

**Minimize:**

```math
	\min_a (\max_{u, \nu}\{\frac{\sum_r u_ry_r(a)}{\sum_l \nu_lx_l(a)}\} - \sum_{i \in \mathcal{D}} a_i\pi_i) \quad (3a)
```

**Subject to:**

```math
	\frac{\sum_r u_ry_{rk}}{\sum_l \nu_lx_{lk}} \leqslant 1 \quad \forall k\in \mathcal{K}'\quad (3b)
```
```math
  D^\intercal a \geqslant d \quad (3c)
```
```math
	u, \nu \geqslant 0, a_i \in \{0,1\} \quad (3e)
```

where the indices $l=1..\mathcal{L}$ and $r=1..\mathcal{R}$ are employed for numbering the DEA inputs and outputs.

Constraint set (3c) here represents application-specific constraints, such as connectivity, flow, distance etc.

This subproblem can be formulated as a bilevel LP. After the Charnes-Cooper transformation (fixing the denominator in the objective), that bilevel problem looks like:

**Minimize:**

```math
	\min_a (\sum_r \mu_ry_r(a_1, ..., a_n) - \sum_{i \in \mathcal{D}} a_i\pi_i) \quad (4a)
```

**Subject to:**

```math
	D^\intercal a \geqslant d  \quad (4b)
```
```math
	\mu, \nu \in argmax_{\mu, \nu}\{\sum_r \mu_ry_r(a_1, ..., a_n): \sum_r \mu_ry_{rk} - \sum_l \nu_lx_{lk} \leqslant 0 ~ \forall k\in K', \sum_l \nu_lx_l(a_1, ..., a_n) = 1\} \quad (4c)
```
```math
	a_i \in \{0,1\} \quad (4e)
```

This bilevel LP can be reformulated using strong duality. It has been pointed out in the bilevel literature, that the Karush-Kuhn-Tucker conditions are the preferable way of reformulating bilevel LPs, 
as using strong duality produces a need to use McCormick envelopes, which, in the case of multiplying two continuous variables, usually leads to worse performance than the existing approaches to handle the 
complementary slackness non-convexity. However, in our case, the multiplication happens between binary and continuous variables and the linearization is very straightforward, and can be handled by modern 
solvers, such as Gurobi.

After replacing the lower-level problem by its optimality conditions, the single-level reformulation of this problem then is:

**Minimize:**

```math
		\min_a (\sum_r \mu_ry_r(a_1, ..., a_n) - \sum_{i \in D} a_i\pi_i) \quad (5a)
```

**Subject to:**

```math
		D^\intercal a \geqslant d \quad (5b)
```
```math
		\sum_r \mu_ry_r(a_1, ..., a_n) = \theta \quad (5c)
```
```math
		\sum_r \mu_ry_{rk} - \sum_l \nu_lx_{lk} \leqslant 0 ~ \forall k\in K' \quad (5d)
```
```math
		\sum_l \nu_lx_l(a_1, ..., a_n) = 1 \quad (5e)
```
```math
		\sum_kx_{lk}\lambda_k \leqslant \theta x_l(a_1, ..., a_n) ~ \forall l \quad (5f)
```
```math
		\sum_ky_{rk}\lambda_k \geqslant y_r(a_1, ..., a_n) \quad (5g)
```
```math
		a_i \in \{0,1\}, ~ \lambda, \mu, \nu \geqslant 0 \quad (5h)
```

where $\theta$ and $\lambda$ are the dual variable of the DEA model in the lower level.

### Integer solutions

Just as with any other column generation framework, the method described in this section yields an optimal solution to a linear relaxation of the original problem. 
This solution is likely to be fractional and therefore we need to go one step further to obtain the integer solution to the original problem. The problem presented in this study is complicated by the fact that 
the DEA scores for all subnetworks in set $\mathcal{K}'$ are recalculated after each iteration. The scores tend to change as the comparison set is populated with new subnetworks, which in turn might affect 
the already processed node in a Branch-and-Price algorithm. Here, we settle for a heuristic method and resolve the final master problem as an integer program.

## Example: Distribution Network Efficiency

We present an example of applying the described method to a distribution network. The network ships a product that can be expressed using aggregated continuous flows between its nodes. Some nodes produce and 
store the product while others can only store it. Nodes can ship the product between each other and to the customers. In the context of DEA, we refer to both individual nodes and subnetworks of nodes as 
_decision-making units (DMUs)_. The summary of DEA features we consider is presented below:

| Inputs       | Description                                                                                                 |
|--------------|-------------------------------------------------------------------------------------------------------------|
| `tot_flow`   | Total outflow to customers                                                                                 |
| `ncust`      | Number of "closest customers"                                                                               |
| `cdfiff`     | Customers difficulty                                                                                        |
| `veh_TE`     | Vehicle capacity                                                                                             |
| Outputs      |                                                                                                             |
| `not_kept`   | Amount of product a DMU sent to other DMUs when it did not have enough product capacity                    |
| `help`       | Amount of product a DMU sent to DMUs in need when it had enough product capacity                           |

**Caption:** DEA inputs and outputs

The concept of ``closest customers" used to define the _ncust_ feature means the customers for whom a given DMU is the closest one; hence it would be logical if that DMU serves them. The concept of closest 
customers is illustrated below:

<p align="center">
<img src="https://github.com/riabova/NetworkEfficiency/blob/main/img/closest_custs.png" width="400">
</p>

We formulate the master problem as a minimal cover problem:

**Minimize:**

```math
	\min \sum_{k \in \mathcal{K}'} c_kz_k \quad (6a)
```

**Subject to:**

```math
	\sum_{k \in \mathcal{K}}\hat{a}_{ik}z_k \geqslant 1 \quad \forall i \in D\quad (6b)
```
```math
	z_k \in \{0, 1\} \quad (6c)
```

where $z_k$ are binary variables equal to 1 if subnetwork $k$ is used in the cover and 0 otherwise.

Binary parameters $\hat{a}_{ik}$ here indicate if DMU $i$ is included in subnetwork $k$ and parameter $c_k$ in the objective represents the DEA score of subnetwork $k$. Note that this problem is always feasible as there exists a solution in which all the $z_k$ are equal to 1 for $k \in \mathcal{D}$ which means that all individual DMUs are taken as their own subnetwork.

If we relax the integrality of $z$'s, the dual of this master problem is:

**Minimize:**

```math
	\max \sum_{i \in \mathcal{D}} \pi_i - \sum_{i \in \mathcal{K}'}\gamma_k\quad (7a)
```

**Subject to:**

```math
	\sum_{i \in \mathcal{D}}a_{ik}\pi_i  - \gamma_k \leq c_k \quad \forall k \in \mathcal{K}'\quad (7b)
```
```math
	\pi_i, \gamma_k >= 0 \quad (7c)
```

The variables $\pi_i$ here are induced by constraints (6b) and the variables $\gamma_k$ come from the upper bounds on $z_k$.

We take the set of individual DMUs as the initial subset of subnetworks, $\mathcal{K}'$. We then use the following subproblem to generate a new subnetwork $k'$ with the negative reduced cost $c_{k'} - \sum_{i \in \mathcal{D}}a_{ik'}\pi_i$, where the DEA score $c_{k'}$ is calculated as the objective function of the DEA model.

**Minimize:**

```math
	\min_a (\max_{u, \nu}\{\frac{\sum_r u_ry_r(a)}{\sum_l \nu_lx_l(a)}\} - \sum_{i \in \mathcal{D}} a_i\pi_i) + \gamma_k \quad (8a)
```

**Subject to:**

```math
	\frac{\sum_r u_ry_{rk}}{\sum_l \nu_lx_{lk}} \leqslant 1 \quad \forall k\in \mathcal{K}'
        \sum_{i \in \mathcal{S}} a_i <= |\mathcal{S}| - 1 ~ \forall \mathcal{S} \subseteq \mathcal{D}: \mathcal{S} \text{ is disconnected} \quad (8b)
```
```math
	\sum_{j \in \mathcal{D}}t_{ij}a_ia_j >= \varepsilon a_i * t_i^{out} - M(1 - s_i^1) \quad \forall i \in \mathcal{D} \quad (8c)
```
```math
	\sum_{j \in \mathcal{D}}t_{ji}a_ia_j >= \varepsilon a_i * t_i^{in} - M(1 - s_i^2) \quad \forall i \in \mathcal{D} \quad (8d)
```
```math
	\sum_{j \in \mathcal{D}}t_{ij}a_ia_j >= \varepsilon a_i * t_j^{out} - M(1 - s_j^1) \quad \forall j \in \mathcal{D} ~ \forall i \in \mathcal{D} \quad (8e)
```
```math
	\sum_{j \in \mathcal{D}}t_{ji}a_ia_j >= \varepsilon a_i * t_j^{in} - M(1 - s_j^2) \quad \forall j \in \mathcal{D} ~ \forall i \in \mathcal{D} \quad (8f)
```
```math
	\sum_{i\in \mathcal{D}}s_i^1 + \sum_{i\in \mathcal{D}}s_i^2 >= 1 \quad (8g)
```
```math
	a_i \in \{0,1\} \quad (8h)
```

Constraints (8b) are the original constraints of the CCR DEA model. The constraint set (8c) is needed to ensure the connectivity of the generated subnetworks. It is exponential in size and was implemented using a callback for delayed constraint generation.

We also added flow constraints (8d - 8h) that ensure that a DMU can be a part of a subnetwork only if it has a significant flow with at least one other DMU in that subnetwork. The parameters $t_{ij}$, $t_i^{in}$ and $t_i^{out}$ correspond to the total (undirected) flow between DMUs $i$ and $j$, total incoming flow of DMU $i$ and total outcoming flow of DMU $i$ respectively. Significance, in that case, is given as a percentage of the total (incoming or outcoming) flow, $\varepsilon$. Constraints (8d - 8g) present four different ways in which a flow can be significant. Particularly, constraint (8e) says that DMU $i$ can be a part of a subnetwork if it ships a significant amount of its outflow to the rest of the subnetwork. Constraint (8e) says that DMU $i$ can be a part of a subnetwork if it receives a significant amount of its inflow from the rest of the subnetwork. Constraint (8f) says that DMU $i$ can be a part of the subnetwork if the flow from $i$ to $j$ is significant for some DMU $j$ in the subnetwork. Constraint (8g) says that DMU $i$ can be a part of the subnetwork if the flow from $j$ to $i$ is significant for some DMU $j$ in the subnetwork. And finally, constraint (8h) says that at least one of these four cases should be true.

Note that the dual variable $\gamma_k$ present in the objective of the subproblem is defined for each subnetwork $k$ and we cannot access its value from the current restricted dual problem RDP solution. It is easy to see that the value of $\gamma_k$ should be 0 for a new subnetwork as one way of viewing it is as a subnetwork not currently in the basis of the restricted master problem, and thus the shadow price of its upper bound is 0. However, it needs careful consideration as ignoring $\gamma_k$ leads to the possibility of the subproblem generating the same subnetwork multiple times as it cannot account for the true reduced cost in that case. We suggest a way of handling that issue with the following block of constraints:

```math
	a_i + \hat{a}_{ik} \geqslant q_{ik} \quad \forall i \in \mathcal{D} ~ \forall k \in \mathcal{K}' \quad (9a)
```
```math
	a_i + \hat{a}_{ik} \leqslant 2 - q_{ik} \quad \forall i \in \mathcal{D} ~ \forall k \in \mathcal{K}' \quad (9b)
```
```math
	\sum_{i \in \mathcal{D}}q_{ik} \geqslant 1 \quad \forall k \in \mathcal{K}' \quad (9c)
```

These constraints are essentially saying that $\forall k \in \mathcal{K}' ~ \exists ~ i \in \mathcal{D}: ~ a_i + \hat{a}_{ik} = 1$ (subnetworks are different in at least one DMU).

The below relationships show how the DEA inputs and outputs for the subnetwork under evaluation can be calculated from the variables $a_i$, assuming that we have the data for all DMUs in $\mathcal{D}$ at hand.

```math
y_{ro} = \{not\_kept(a), help(a)\}
```
```math
not\_kept(a_1..a_n) = \sum_{i\in \mathcal{D}}not\_kept_ia_i
```
```math
help(a_1..a_n) = \sum_ihelp_ia_i
```

```math
x_{ro} = \{tot\_flow(a), ncust(a), cdiff(a), veh\_TE(a)\}
```
```math
tot\_flow(a_1..a_n) = \sum_{i\in \mathcal{D}}tot\_flow_ia_i
```
```math
ncust(a_1..a_n) = \sum_{i\in \mathcal{D}}ncust_ia_i
```
```math
cdiff(a_1..a_n) = \sum_{i\in \mathcal{D}} \frac{ncust_i}{ncust(a_1..a_n)}cdiff_i a_i
```
```math
veh\_TE(a_1..a_n) = \sum_{i\in \mathcal{D}}veh\_TE_ia_i
```

It can be observed that in the above definition, there is a non-linear relationship between $cdiff(a_1..a_n)$ and $ncust(a_1..a_n)$. We used the following constraints along with the built-in non-linearity support in Gurobi to model that:
```math
        f_1 * ncust(a_1..a_n) = \sum_{i\in\mathcal{D}}ncust_i*cdiff_i*a_i 
```
```math
        f_2 = - f_1 + \max_{i\in \mathcal{K}'}\{cdiff_i\} + \min_{i\in \mathcal{K}'}\{cdiff_i\}
```

where $f_1$ denotes the raw value of $cdiff$ (before the transformation for undesirable DEA features described in Section 3.2) and $f_2$ denotes the transformed $cdiff$. Constraint (10a) follows directly from the $cdiff$ definition given above. Constraint (10b) denotes the transformation for undesirable features from Section 3.2.
