---
title: Liouville's Theorem derivation
layout: post
author: wen
tags:
- physics
- derivation
---

In classical mechanics the liouville's theorem states that the intgeral element
\begin{align}
\partial p_0 \partial q_0 = \partial p_t \partial q_t
\end{align}
where $$p_0, q_0 $$ and $$p_t, q_t$$are the generalized coordinates at time 0 and t, respectively.

To prove this, we notice the Hamiltonian equation
\begin{align}
\dot p_0 = -\frac{\partial H}{\partial q_0},
\dot q_0 = \frac{\partial H}{\partial p_0}
\end{align}
and
\begin{align}
\dot p_t = -\frac{\partial H}{\partial q_t},
\dot q_t = \frac{\partial H}{\partial p_t}
\end{align}

Since $$q_t=q_t(p_0,q_0)$$ is a function of $$p_0$$, and $$q_0$$, we can write
\begin{align}
\dot q_t = \frac {\partial q_t}{\partial p_0} \cdot \dot p_0 + \frac{\partial q_t}{\partial q_0} \cdot  \dot q_0 = \frac{\partial q_t}{\partial p_0} \cdot (-\frac{\partial H}{\partial q_0}) + \frac{\partial q_t}{\partial q_0} \cdot \frac{\partial H}{\partial p_0} = -\frac{\partial q_t}{\partial p_0} \cdot (\frac{\partial p_t}{\partial q_0} \cdot \frac{\partial H}{\partial p_t} + \frac{\partial q_t}{\partial q_0} \cdot \frac{\partial H}{\partial q_t})+ \frac{\partial q_t}{\partial q_0} \cdot (\frac{\partial p_t}{\partial p_0} \cdot \frac{\partial H}{\partial p_t} + \frac{\partial q_t}{\partial p_0} \cdot \frac{\partial H}{\partial q_t}) \\\
\dot q_t = (-\frac{\partial q_t}{\partial p_0} \cdot \frac{\partial p_t}{\partial q_0} + \frac{\partial q_t}{\partial q_0} \cdot \frac{\partial p_t}{\partial p_0}) \cdot \frac{\partial H}{\partial p_t} + (-\frac{\partial q_t}{\partial p_0} \cdot \frac{\partial q_t}{\partial q_0} + \frac{\partial q_t}{\partial q_0} \cdot \frac{\partial q_t}{\partial p_0}) \cdot \frac{\partial H}{\partial q_t} = (-\frac{\partial q_t}{\partial p_0} \cdot \frac{\partial p_t}{\partial q_0} + \frac{\partial q_t}{\partial q_0} \cdot \frac{\partial p_t}{\partial p_0}) \cdot \frac{\partial H}{\partial p_t}
\end{align}

Combing this equation with the Hamiltonian equation for $$\dot q_t$$, we obtained
\begin{align}
\frac{\partial q_t}{\partial q_0} \cdot \frac{\partial p_t}{\partial p_0} - \frac{\partial q_t}{\partial p_0} \cdot \frac{\partial p_t}{\partial q_0} = 1
\end{align}

That is
\begin{align}
det (J(q_t, p_t, q_0, p_0)) = det(\frac{\partial (q_t, p_t)}{\partial (q_0, p_0)}) = \begin{vmatrix}
\frac{\partial q_t}{\partial q_0} & \frac{\partial q_t}{\partial p_0}\\\ \frac{\partial p_t}{\partial q_0} & \frac{\partial p_t}{\partial p_0}
\end{vmatrix} = 1
\end{align}
Where J is the Jacobian matrix.
So,
\begin{align}
\partial p_t \partial q_t = det\lvert J \rvert \cdot\partial p_0 \partial q_0 = \partial p_0 \partial q_0
\end{align}

Q.E.D