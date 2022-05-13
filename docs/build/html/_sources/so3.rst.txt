SO(3) - so(3) functions
==================================================

An element of SO(3) represent a matrix rotation in :math:`\mathbb{R}^3`.
Such an element can be generated from an element :math:`\phi` in so(3).

As so(3) is a vector space, we are able to do probabilities such as expected value, things that are not possible on a group (i.e., SO(3)).

Elements of so(3) are 3x3 skew-symetric matrices of the form 

.. math::
  M = \phi^\wedge = 
  \begin{bmatrix}
    \phi_1 \\ \phi_2 \\ \phi_3
  \end{bmatrix}^\wedge
  =
  \begin{bmatrix}
      0 & -\phi_3 & \phi_2 \\
      \phi_3 & 0 & -\phi_1 \\
      -\phi_2 & \phi_1 & 0
  \end{bmatrix}.

Do note that elements of so(3) are **matrices** and not vectors.
However, we prefer using the *coordinates* while working with these elements.
Refer to `wikipedia <https://en.wikipedia.org/wiki/Basis_(linear_algebra)>`_ for a precise definition of these terms.

The (linear) operation :math:`(\cdot)^\wedge` maps a coordinate to its associated element in so(3).

An element of so(3) is mapped on SO(3) using the exponential function, and an element of SO(3) can be sent to so(3) using the logarithm function.
For this mapping to be bijective, we need to restrict elements of so(3) in the ball of radius :math:`\pi`, that is

.. math::
  ||\phi|| < \pi.

**Careful:** The exponential mapping is not homomorphic, which translates to the inequality

.. math::
  \exp(\phi_1^\wedge)\exp(\phi_2^\wedge) \neq \exp(\phi_1^\wedge+\phi_2^\wedge).

The true result relies on the *Baker-Campbell-Hausdorff* formula, that is quite complicated and is an infinite series.
However, if one of the component is small, we have the following useful approximation:

.. math::
  \exp(\phi_1^\wedge)\exp(\phi_2^\wedge) \approx \exp( [\phi_1+J(\phi_1)^{-1}\phi_2]^\wedge) \quad\text{if } \phi_2 \text{ is small},

where :math:`J(\phi_1)` is the (left) Jacobian of SO(3).

Basic so(3) operations
-----------------------------------------
.. autofunction:: liblie.wedge_so3
.. autofunction:: liblie.vee_so3

SO(3) - so(3) mappings
-------------------------------------------------------
.. autofunction:: liblie.logm_so3
.. autofunction:: liblie.expm_so3

Left Jacobian of SO(3)
-----------------------------------------
.. autofunction:: liblie.jac_so3
.. autofunction:: liblie.jac_inv_so3

