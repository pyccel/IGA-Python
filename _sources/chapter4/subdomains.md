# Subdomains


In this section, we consider a domain $\Omega$ which a union of multiple subdomain, *i.e.*

$$
\Omega := \bigcup_{i=1}^n \Omega_i, \quad \Omega_i \bigcap \Omega_j = \emptyset, ~ i \neq j  
$$

We shall denote by $\mathcal{I}$ the set of all internal interfaces of $\Omega$.

We shall also need the following operators

- The jump of the function $u$, defined as $[\![ u ]\!] := u|\_{\Omega_{i_1}} -  u|\_{\Omega_{i_2}}$ for two adjacent subdomains $\Omega_{i_1}$ and $\Omega_{i_2}$
- The average of the function $u$, defined as $\{u\} := \frac{1}{2} \left( u|\_{\Omega_{i_1}} +  u|\_{\Omega_{i_2}} \right)$ for two adjacent subdomains $\Omega_{i_1}$ and $\Omega_{i_2}$

In general, these operators are defined based on the considered variational formulation. We shall see in the following examples how to define them.
