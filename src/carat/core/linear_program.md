### Linear Program Optimization Formulation

$$
\begin{alignedat}{4}
& \text{minimize} && \sum_{cbgpse}(z_{cbgpse}-q_{cbgpse})+\sum_{cpse}(z_{cpse}-q_{cpse}) \quad\quad\quad\quad\quad && &&(a) \\
& \text{s.t.} && && && \\
& && \beta_{cbgpsea} = \sum_{p's'}\psi_{p's'pse}\ \beta_{cp's'ea} \:,
&&\quad \forall c,b,g,p,s,e,a \quad && (b) \\
& && \beta_{cpsea} = \sum_{c'b'g'}\mu_{c'b'g'cp}\ \beta_{c'b'g'psea} \:,
&&\quad \forall c,p,s,e,a && (c) \\
& && \sum_a \beta_{cbgpsea}-z_{cbgpse}-q_{cbgpse}=1,
&&\quad \forall c,b,g,p,s,e && (d) \\
& && \sum_a \beta_{cpsea}-z_{cpse}-q_{cpse}=1,
&&\quad \forall c,p,s,e && (e) \\
& && \beta_{cbgpsea} \in [0,1],
&&\quad \forall c,b,g,p,s,e,a && (f) \\
& && \beta_{cpsea} \in [0,1],
&&\quad \forall c,p,s,e,a && (g) \\
& && z_{cbgpse} \in \mathbb{R}^+, \quad q_{cbgpse} \in \mathbb{R}^-,
&&\quad \forall c,b,g,p,s,e && (h) \\
& && z_{cpse} \in \mathbb{R}^+, \quad q_{cpse} \in \mathbb{R}^-,
&&\quad \forall c,p,s,e && (i)
\end{alignedat}
$$


*N.B. $c'$ denotes the inlet company code, whereas $c$ denotes the outlet company code.*

**Table 1. Value chain indices**

### Notation for indices used in the value chain model

| Index | Description |
|-------|-------------|
| $a$ | Elemental attribute (e.g., biogenic, fossil, etc.) |
| $b$ | Business process, anonymized coding: PLNTb |
| $c$ | Company code, anonymized coding: COMPc |
| $e$ | Element (e.g., carbon) |
| $g$ | Main product, same structure as $p$ |
| $p$ | Product, anonymized coding: PRODp |
| $s$ | Substance, represented by SMILES |


**Table 2. Decision variables, slack variables, and parameters**

| Notation | Description |
|----------|-------------|
| $\beta_{cbgpsea}$ | Fraction of elemental attribute $a$ of element $e$ in substance $s$, material $p$, at production node $(c, b, g) $. |
| $\beta_{cpsea}$ | Fraction of elemental attribute $a$ of element $e$ in substance $s$, material $p$, at mix node $(c, p) $. |
| $z_{cbgpse}$ | Positive slack variable for element $e$ in substance $s$, material $p$, at production node $(c, b, g)$. |
| $q_{cbgpse}$ | Negative slack variable for element $e$ in substance $s$, material $p$, at production node $(c, b, g)$. |
| $z_{cpse}$ | Positive slack variable for element $e$ in substance $s$, material $p$, at mix node $(c, p)$. |
| $q_{cpse}$ | Negative slack variable for element $e$ in substance $s$, material $p$, at mix node $(c, p) $. |
| $\mu_{c'b'g'cp}$ | Mix node share, i.e., the fraction of a virtual tank $(c,p)$ sourced from a production node $(c',b',g') $. |
| $\psi_{p's'pse}$ | Bill of atoms, i.e., the fraction of element $e$ in substance $s$ in product $p$, sourced from substance $s'$ in product $p'$. |
