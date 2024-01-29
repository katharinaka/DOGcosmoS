.. _rotation_curve_shapes:

Parameterisation of rotation curve shapes
------------------------------------------

In order to quantify the shapes of rotations curves, the rotation velocity is measured at an inner fiducial radius which is defined as

.. math::

  r_{fid} = 2\ v_{max} / (70\ km/s)\ kpc. 

The velocity at this fiducial radius, vfid, indicates if the inner rotation curve is more steeply rising (i.e. the galaxy has a cuspy density profile such as NFW) or more slowly rising (i.e. the galaxy has a more cored profile). Dividing by vmax gives us the normalised shape parameter

.. math::

  \eta_{rot} = v_{fid} / v_{max}.

With this parameterisation you can investigate the distribution of rotation curve shapes for your galaxy sample of interest.


.. automodule:: rotation_curve_shapes
   :members:
