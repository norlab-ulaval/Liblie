SE(3) - se(3) functions
==================================================
An element of SE(3) represent a transformation matrix in :math:`\mathbb{R}^3`.
Such an element can be generated from an element :math:`\xi` in se(3).

As se(3) is a vector space, we are able to do probabilities such as expected value, things that are not possible on a group (i.e., SE(3)).

Elements of se(3) are 4x4 matrices of the form 

.. math::
  M = \xi^\wedge = 
  \begin{bmatrix}
    \rho \\ \phi 
  \end{bmatrix}^\wedge
  =
  \begin{bmatrix}
  \phi^\wedge & \rho \\
  0^T & 0
  \end{bmatrix},

where :math:`\phi^\wedge` corresponds to the element of so(3) of coordinate :math:`\phi`.
See :ref:`SO(3) - so(3) functions` for details on this operation.

Do note that elements of se(3) are **matrices** and not vectors.
However, we prefer using the *coordinates* while working with these elements.
Refer to `wikipedia <https://en.wikipedia.org/wiki/Basis_(linear_algebra)>`_ for a precise definition of these terms.

The (linear) operation :math:`(\cdot)^\wedge` maps a coordinate to its associated element in se(3).

An element of se(3) is mapped on SE(3) using the exponential function, and an element of SE(3) can be sent to se(3) using the logarithm function.
For this mapping to be bijective, the restriction is the same as in `SO(3) - so(3) functions`.
We need to restrict the rotation element :math:`\phi` of :math:`\xi=[\rho^T\, \phi^T]^T` se(3) in the ball of radius :math:`\pi`, that is

.. math::
  ||\phi|| < \pi.

**Careful:** The exponential mapping is not homomorphic, which translates to the inequality

.. math::
  \exp(\xi_1^\wedge)\exp(\xi_2^\wedge) \neq \exp(\xi_1^\wedge+\xi_2^\wedge).

The true result relies on the *Baker-Campbell-Hausdorff* formula, that is quite complicated and is an infinite series.
However, if one of the component is small, we have the following useful approximation:

.. math::
  \exp(\xi_1^\wedge)\exp(\xi_2^\wedge) \approx \exp( [\xi_1+\mathcal{J}(\xi_1)^{-1}\xi_2]^\wedge) \quad\text{if } \xi_2 \text{ is small},

where :math:`\mathcal{J}(\xi_1)` is the (left) Jacobian of SE(3).

Basic se(3) operations
-----------------------------------------
.. autofunction:: liblie.wedge_se3
.. autofunction:: liblie.curly_wedge_se3
.. autofunction:: liblie.vee_se3

SE(3) - se(3) mappings
-------------------------------------------------------
.. autofunction:: liblie.logm_se3
.. autofunction:: liblie.expm_se3

Left Jacobian of se(3)
-----------------------------------------
.. autofunction:: liblie.jac_se3
.. autofunction:: liblie.jac_inv_se3

