#!/usr/bin/env python
# coding: utf-8

# (ch08)=
# # Strain

# Stresses acting through time within the Earth can lead to deformation. Deformation is more complicated than stress because it involves the comparison of states of the rock material at two different points in time. Therefore, one needs to establish both temporal and spatial reference frames. This chapter covers deformation and strain, including infinitesimal, finite, and progressive strain. This is a relatively long and complicated material, it comprises several chapters in classical books such as Ramsay (1967) and Means (1976). However, our purpose is not to discuss in detail the theory of strain, but rather introduce some interesting strain problems in geosciences, and their computation.

# (ch08-1)=
# ## Deformation and strain
# 
# Strictly speaking, deformation involves rigid body deformation (translation and rotation), and non-rigid body deformation (changes in shape and volume) or strain. In geology, rigid-body deformation is rather difficult to determine, with a few exceptions (e.g. paleolatitude of a continent from paleomagnetic pole). We often just focus on strain.
# 
# Consider the deformation shown in {numref}`Figure %s <ch08_fig01>`. The square and inscribed circle ({numref}`Figure %s <ch08_fig01>`a) are first translated, then rotated, and then sheared. The translation vector **t** ({numref}`Figure %s <ch08_fig01>`b) and rotation $\omega$ ({numref}`Figure %s <ch08_fig01>`c) describe the first stage of rigid body-deformation. Shear of the objects defines the second stage of non-rigid body deformation or strain ({numref}`Figure %s <ch08_fig01>`d).

# ```{figure} /figures/ch08_fig01.png
# :width: 400px
# :name: ch08_fig01
# 
# Deformation of a square of side $l_{i}$ and inscribed circle. **a.** Initial configuration. **b.** Translation $\textbf{t}$. **c.** Rotation $\omega$. **d.** Shear. $l_{f}$ and $\psi$ are the final length and angular shear of the long side of the parallelogram. Red and blue lines in ellipse are principal axes and lines of no finite elongation (LNFE), respectively.
# ```

# The changes in shape (distortion) of the bodies can be described by the changes in length and angle of lines. If the initial length of a line is $l_i$, and the deformed length of the line is $l_f$, the change in length of the line can be defined by either one of the following parameters:

# $$
# \begin{gathered}
#     e = \frac{l_f-l_i}{l_i} \\
#     S = \frac{l_f}{l_i} = 1 + e \\
#     \lambda = S^2 = (1+e)^2
# \end{gathered}
# $$ (ch08_eq01)

# where $e$ is elongation, $S$ is stretch, and $\lambda$ is quadratic elongation. For the long side of the parallelogram in {numref}`Figure %s <ch08_fig01>`d, $e$ is positive, and $S$ and $\lambda$ are larger than $1$. For the short side, $e$ is negative, and $S$ and $\lambda$ are lower than $1$. Lines that have their original length have $e = 0$, and $S$ and $\lambda = 1$ (blue lines, {numref}`Figure %s <ch08_fig01>`d). Notice than $e$ cannot be lower than $-1$, and $S$ and $\lambda$ cannot be lower than 0 (a line can’t be shortened more than its original length).
# 
# Changes in angle are measured by the angular shear $\psi$, which is the change in angle between two originally perpendicular lines ({numref}`Figure %s <ch08_fig01>`d). From the angular shear, one can calculate the shear strain $\gamma$:

# $$
# \gamma = \tan\psi
# $$ (ch08_eq02)

# Changes in area or volume (dilation) can be described using the areal or volumetric stretch:

# $$
# \begin{gathered}
#     S_A = \frac{A_f}{A_i} \\
#     S_V = \frac{V_f}{V_i}
# \end{gathered}
# $$ (ch08_eq03)

# where $A_i$ and $A_f$ are initial and final area, and $V_i$ and $V_f$ are initial and final volume. We can also determine areal or volumetric elongation.
# 
# In this chapter, we make two major assumptions about strain: Strain is continuous (it is distributed uniformly across the body), and strain is homogeneous (it is identical across the body). Clearly these assumptions are incorrect: Rocks are full of discontinuities, and geological structures (e.g. folds) exhibit heterogeneous strain. However at the appropriate scale, rocks can be described as continuum materials, and the heterogeneous strain of geological structures can be represented by domains of homogeneous (yet compatible) strain. This makes possible applying homogeneous strain to geological deformation.
# 
# For homogeneous strain, straight lines remain straight, parallel lines remain parallel, and circles become ellipses in 2D, or spheres become ellipsoids in 3D ({numref}`Figure %s <ch08_fig01>`d). The resultant ellipse (or ellipsoid) is called the strain ellipse (or ellipsoid). The stretches along the axes of this ellipse (or ellipsoid, red lines in {numref}`Figure %s <ch08_fig01>`d) are called the principal stretches, and they are denoted by the symbols $S_1$, $S_2$, $S_3$, for the maximum, intermediate, and minimum principal stretch, respectively. The volumetric stretch can also be defined as:

# $$
# S_V = S_1S_2S_3
# $$ (ch08_eq04)

# The deformation in {numref}`Figure %s <ch08_fig01>` is volume and area constant ($S_V$ and $S_A = 1$). Under this condition, there are two lines in the strain ellipse that have their original length (blue lines, {numref}`Figure %s <ch08_fig01>`d). These lines are known as the lines of no finite elongation (LNFE).

# (ch08-2)=
# ## Deformation and displacement gradients
# 
# Consider the deformation shown in {numref}`Figure %s <ch08_fig02>`. The cube is transformed into a brick-shaped body. It is shortened along the $\mathbf{X_2}$ axis by half ($e = -0.5$, $S = 0.5$), and it is stretched along the $\mathbf{X_3}$ axis twice ($e = 1$, $S = 2$).

# ```{figure} /figures/ch08_fig02.png
# :width: 400px
# :name: ch08_fig02
# 
# Deformation of a cube in to a brick-shaped body. The principal axes of strain are parallel to the coordinate axes before $\textbf{X}$ and after $\textbf{x}$ deformation. Modified from Means (1976).
# ```

# We can express the deformed coordinates $\mathbf{x}$ of any point in the cuboid in terms of the undeformed coordinates $\mathbf{X}$ of the same point in the cube:
# 
# :::{note}
# By convention, we use $\mathbf{X}$ to refer to the undeformed coordinates, and $\mathbf{x}$ to denote the deformed coordinates.
# :::

# $$
# \begin{gathered}
#     x_1 = 1 X_1 + 0 X_2 + 0 X_3 \\
#     x_2 = 0 X_1 + 0.5 X_2 + 0 X_3 \\
#     x_3 = 0 X_1 + 0 X_2 + 2 X_3
# \end{gathered}
# $$ (ch08_eq05)

# Likewise, we can express the undeformed coordinates $\mathbf{X}$ of any point in the cube in terms of the deformed coordinates $\mathbf{x}$ of the same point in the cuboid:

# $$
# \begin{gathered}
#     X_1 = 1 x_1 + 0 x_2 + 0 x_3 \\
#     X_2 = 0 x_1 + 2 x_2 + 0 x_3 \\
#     X_3 = 0 x_1 + 0 x_2 + 0.5 x_3
# \end{gathered}
# $$ (ch08_eq06)

# Equations {eq}`ch08_eq05` and {eq}`ch08_eq06` are called *coordinate transformations*, but they are fundamentally different from the transformation of coordinate axes in [Chapter 5](ch05). The coordinate transformations here are between two different states in time. Eq. {eq}`ch08_eq05` is called a *Green* transformation (new in terms of old), while Eq. {eq}`ch08_eq06` is a *Cauchy* transformation (old in terms of new).
# 
# We can take the partial derivatives of these equations, which are just the coefficients of the equations:

# $$
# \frac{\partial x_i}{\partial X_j}=\begin{bmatrix}1&0&0\\0&0.5&0\\0&0&2\end{bmatrix}
# $$ (ch08_eq07)

# and:

# $$
# \frac{\partial X_i}{\partial x_j}=\begin{bmatrix}1&0&0\\0&2&0\\0&0&0.5\end{bmatrix}
# $$ (ch08_eq08)

# These partial derivatives are known as the deformation gradients and unlike the coordinate transformations, they are homogeneous and independent of the coordinate axes (i.e. they are tensors). Eq. {eq}`ch08_eq07` is called the *Green* deformation gradient, while Eq. {eq}`ch08_eq08` is the *Cauchy* deformation gradient. Since in {numref}`Figure %s <ch08_fig02>` the principal axes of strain are parallel to the coordinate axes, we can easily see one thing: The components of the Green deformation gradient are the stretches along the coordinate axes, and the components of the Cauchy deformation gradient are the inverse of the stretches along the coordinate axes. For more complex situations this will not be the case, but in general the Green and Cauchy deformation gradient tensors are related to the stretch (Allmendinger et al., 2012).
# 
# Another way to study the deformation is to look at the displacement $\mathbf{u}$ between the undeformed and deformed states ({numref}`Figure %s <ch08_fig02>`). The displacement can either be expressed in the undeformed configuration:

# $$
# \begin{gathered}
#     u_1 = 0 X_1 + 0 X_2 + 0 X_3 \\
#     u_2 = 0 X_1 - 0.5 X_2 + 0 X_3 \\
#     u_3 = 0 X_1 + 0 X_2 + 1 X_3
# \end{gathered}
# $$ (ch08_eq09)

# or in the deformed configuration:

# $$
# \begin{gathered}
#     u_1 = 0 x_1 + 0 x_2 + 0 x_3 \\
#     u_2 = 0 x_1 - 1 x_2 + 0 x_3 \\
#     u_3 = 0 x_1 + 0 x_2 + 0.5 x_3
# \end{gathered}
# $$ (ch08_eq10)

# Eq. {eq}`ch08_eq09` is known as the *Lagrangian* displacement (in terms of old coordinates), while Eq. {eq}`ch08_eq09` is the *Eulerian* displacement (in terms of new coordinates). We can take the partial derivatives of these equations, which are just the coefficients of the equations:

# $$
# \frac{\partial u_i}{\partial X_j}=\begin{bmatrix}0&0&0\\0&-0.5&0\\0&0&1\end{bmatrix}
# $$ (ch08_eq11)

# and:

# $$
# \frac{\partial u_i}{\partial x_j}=\begin{bmatrix}0&0&0\\0&-1&0\\0&0&0.5\end{bmatrix}
# $$ (ch08_eq12)

# These partial derivatives are known as the displacement gradients and unlike the displacement, they are homogeneous and independent of the coordinate axes (i.e. they are tensors). Eq. {eq}`ch08_eq11` is called the *Lagrangian* displacement gradient, while Eq. {eq}`ch08_eq12` is the *Eulerian* displacement gradient. For the case of {numref}`Figure %s <ch08_fig02>`, we can see that the components of the Lagrangian displacement gradient are the elongations along the coordinate axes, while the components of the Eulerian displacement gradient are the elongations along the coordinate axes with respect to the deformed reference frame ($\tilde{e}=(l_f-l_i)/l_f$). This will not be the case for more complicated situations, but in general the Lagrangian and Eulerian displacement gradient tensors are related to the elongation (Allmendinger et al., 2012).
# 
# {numref}`ch08_tab01` summarizes this section.

# ```{list-table} Coordinate transformations and displacement via deformation and displacement gradient tensors. Rigid body deformation is not included.
# :header-rows: 1
# :name: ch08_tab01
# 
# * - 
#   - Old coordinates
#   - New coordinates
# * - Coordinate transformations
#   - Green $x_i=\frac{\partial x_i}{\partial X_j}X_j$
#   - Cauchy $X_i=\frac{\partial X_i}{\partial x_j}x_j$
# * - Displacement
#   - Lagrangian $u_i=\frac{\partial u_i}{\partial X_j}X_j$
#   - Eulerian $u_i=\frac{\partial u_i}{\partial x_j}x_j$
# ```

# (ch08-3)=
# ## Infinitesimal strain
# 
# Suppose that a line parallel to the $\mathbf{X_1}$ axis is elongated 1$\%$ of its initial length. In this case the displacement gradient is:

# $$
# \frac{\partial u_1}{\partial X_1} = \frac{0.01}{1.0} = 0.01 \quad\quad \text{and} \quad\quad \frac{\partial u_1}{\partial x_1} = \frac{0.01}{1.01} = 0.0099
# $$ (ch08_eq13)

# Thus, when strains are small ($e < 1\%$), $\frac{\partial u_i}{\partial X_i} \approx \frac{\partial u_i}{\partial x_i}$, and the difference between the displacement gradient in the initial or final state is insignificant. Small strains are called *infinitesimal strains*, and they are important in a number of fields in geosciences, particularly in geophysics.
# 
# For infinitesimal strain, the distinction between old and new coordinates is irrelevant, and therefore we can just use one column in {numref}`ch08_tab01`. The displacement gradient $\mathbf{e}$ is an asymmetric tensor, and it can be decomposed into a symmetric and an antisymmetric tensor (Eq. {eq}`ch08_eq14`):
# 
# :::{note}
# We use $\mathbf{e}$ to denote the displacement gradient. This is different than the non-bold italic letter $e$ which is used for elongation.
# :::

# $$
# e_{ij} = \varepsilon_{ij} + \omega_{ij}
# $$ (ch08_eq14)

# where:

# $$\label{eq8.14} \varepsilon_{ij}=\frac{1}{2}(e_{ij}+e_{ji})=\begin{bmatrix}e_{11}&\frac{e_{12}+e_{21}}{2}&\frac{e_{13}+e_{31}}{2}\\ \frac{e_{21}+e_{12}}{2}&e_{22}&\frac{e_{23}+e_{32}}{2}\\ \frac{e_{31}+e_{13}}{2}&\frac{e_{32}+e_{23}}{2}&e_{33}\end{bmatrix}
# $$ (ch08_eq15)

# and:

# $$
# \omega_{ij}=\frac{1}{2}(e_{ij}-e_{ji})=\begin{bmatrix}0&\frac{e_{12}-e_{21}}{2}&\frac{e_{13}-e_{31}}{2}\\ \frac{e_{21}-e_{12}}{2}&0&\frac{e_{23}-e_{32}}{2}\\ \frac{e_{31}-e_{13}}{2}&\frac{e_{32}-e_{23}}{2}&0\end{bmatrix}
# $$ (ch08_eq16)

# The symmetric tensor $\boldsymbol{\varepsilon}$ is called the strain tensor. The diagonal, $\varepsilon_{ii}$, terms of this tensor correspond to the elongations along the axes of the reference system, and the off-diagonal terms, $\varepsilon_{ij}$, correspond to half the shear strain ($\varepsilon_{ij} = \frac{\gamma_{ij}}{2}$). Like any symmetric tensor, the strain tensor has three principal axes which can be determined by computing the eigenvalues and eigenvectors of the tensor ([Section 6.2](ch06-2)). These principal axes define the infinitesimal strain ellipsoid.
# 
# The antisymmetric tensor $\boldsymbol{\omega}$ is also known as an axial vector. The Cartesian coordinates, $r_i$, of this vector give the orientation of the rotation axis:

# $$
# r_1=\frac{-(\omega_{23}-\omega_{32})}{2}\quad r_2=\frac{-(-\omega_{13}+\omega_{31})}{2}\quad \text{and}\quad r_3=\frac{-(\omega_{12}-\omega_{21})}{2}
# $$ (ch08_eq17)

# The amount of rotation in radians is just the length of the vector $\mathbf{r}$:

# $$
# \omega = \vert \mathbf{r}\vert = \sqrt{r_1^2+r_2^2+r_3^2}
# $$ (ch08_eq18)

# The function [inf_strain](https://github.com/nfcd/compGeo/blob/master/source/functions/inf_strain.py) computes the infinitesimal strain and rotation tensors from a displacement gradient. It also outputs the principal strains, and the components and amount of rotation:

# In[1]:


import numpy as np

from compgeo.cart_to_sph import cart_to_sph
from compgeo.zero_twopi import zero_twopi


def inf_strain(e):
    """
    inf_strain computes infinitesimal strain from an input
    displacement gradient tensor

    USE: eps,ome,pstrain,rotc,rot = inf_strain(e)

    e = 3 x 3 displacement gradient tensor
    eps = 3 x 3 strain tensor
    ome = 3 x 3 rotation tensor
    pstrain = 3 x 3 matrix with magnitude (column 1), trend
            (column 2) and plunge (column 3) of maximum
            (row 1), intermediate (row 2),and minimum 
            (row 3) principal strains
    rotc = 1 x 3 vector with rotation components
    rot = 1 x 3 vector with rotation magnitude and trend
        and plunge of rotation axis

    NOTE: Output trends and plunges of principal strains
        and rotation axes are in radians

    Python function translated from the Matlab function
    InfStrain in Allmendinger et al. (2012)
    """
    # Initialize variables
    eps = np.zeros((3,3))
    ome = np.zeros((3,3))
    pstrain = np.zeros((3,3))
    rotc = np.zeros(3)
    rot = np.zeros(3)

    # Compute strain and rotation tensors
    for i in range(3):
        for j in range(3):
            eps[i,j]=0.5*(e[i,j]+e[j,i])
            ome[i,j]=0.5*(e[i,j]-e[j,i])

    # Compute principal strains and orientations.
    # Here we use the function eigh. D is a vector
    # of eigenvalues (i.e. principal strains), and V is a
    # full matrix whose columns are the corresponding
    # eigenvectors (i.e. principal strain directions)
    D,V = np.linalg.eigh(eps)

    # Maximum principal strain
    pstrain[0,0] = D[2]
    pstrain[0,1],pstrain[0,2] = cart_to_sph(V[0,2],V[1,2],V[2,2])
    
    # Intermediate principal strain
    pstrain[1,0] = D[1] 
    pstrain[1,1],pstrain[1,2] = cart_to_sph(V[0,1],V[1,1],V[2,1])
    
    # Minimum principal strain
    pstrain[2,0] = D[0] 
    pstrain[2,1],pstrain[2,2] = cart_to_sph(V[0,0],V[1,0],V[2,0])

    # Calculate rotation components
    rotc[0]=(ome[1,2]-ome[2,1])*-0.5
    rotc[1]=(-ome[0,2]+ome[2,0])*-0.5
    rotc[2]=(ome[0,1]-ome[1,0])*-0.5

    # Compute rotation magnitude
    rot[0] = np.sqrt(rotc[0]**2+rotc[1]**2+rotc[2]**2)
    
    # Compute trend and plunge of rotation axis
    rot[1],rot[2] = cart_to_sph(rotc[0]/rot[0],
                                rotc[1]/rot[0],
                                rotc[2]/rot[0])
    
    # If plunge is negative
    if rot[2] < 0:
        rot[1] = zero_twopi(rot[1]+np.pi)
        rot[2] *= -1
        rot[0] *= -1

    return eps, ome, pstrain, rotc, rot


# (ch08-3-1)=
# ### Mohr circle for infinitesimal strain
# 
# As discussed in [Section 6.4.1](ch06-4-1), the rotation of the infinitesimal strain tensor about one principal axis can be represented by a Mohr circle. We start with the infinitesimal strain tensor in principal coordinates:

# $$
# \varepsilon_{ij}=\begin{bmatrix}\varepsilon_1&0&0\\ 0&\varepsilon_2&0\\ 0&0&\varepsilon_3\end{bmatrix}
# $$ (ch08_eq19)

# and perform a rotation $\theta$ about $\mathbf{X_2}$ ({numref}`Figure %s <ch08_fig03>`a), which is described by the transformation matrix:

# $$
# a_{ij} = \begin{pmatrix}\cos\theta&0&\sin\theta\\0&1&0\\-\sin\theta&0&\cos\theta\end{pmatrix}
# $$ (ch08_eq20)

# The tensor transformation equation is:

# $$
# \varepsilon_{ij}' = a_{ik}a_{jl}\varepsilon_{kl}
# $$ (ch08_eq21)

# which results in the new form of the strain tensor:

# $$
# \varepsilon'_{ij}=\begin{bmatrix}\varepsilon_1\cos^2\theta+\varepsilon_3\sin^2\theta&0&-(\varepsilon_1-\varepsilon_3)\sin\theta\cos\theta\\0&\varepsilon_2&0\\-(\varepsilon_1-\varepsilon_3)\sin\theta\cos\theta&0&\varepsilon_1\sin^2\theta+\varepsilon_3\cos^2\theta\end{bmatrix}
# $$ (ch08_eq22)

# ```{figure} /figures/ch08_fig03.png
# :width: 600px
# :name: ch08_fig03
# 
# Rotation of the infinitesimal strain tensor about principal axis $\mathbf{X_2}$. **a.** Physical plane representation, **b.** Mohr circle representation. Modified from Allmendinger et al. (2012).
# ```

# Upon rearranging, we get:

# $$
# \begin{gathered}
#     \varepsilon_{11}'=\frac{\varepsilon_1+\varepsilon_3}{2}+\frac{\varepsilon_1-\varepsilon_3}{2}\cos 2\theta \\
#     \varepsilon_{33}'=\frac{\varepsilon_1+\varepsilon_3}{2}-\frac{\varepsilon_1-\varepsilon_3}{2}\cos 2\theta \\
#     \varepsilon_{13}'=\varepsilon_{31}'=\frac{\gamma}{2}=-\frac{\varepsilon_1-\varepsilon_3}{2}\sin 2\theta
# \end{gathered}
# $$ (ch08_eq23)

# which are the equations of the Mohr circle for infinitesimal strain ({numref}`Figure %s <ch08_fig03>`b). This Mohr circle is not very useful, but perhaps the most important thing illustrated by it is that the two directions of maximum shear strain are at $\pm$<!-- -->45to the principal axes, $\varepsilon_1$ and $\varepsilon_3$ ({numref}`Figure %s <ch08_fig03>`b). Turning this around for the infinitesimal strain in a shear zone, the infinitesimal shortening and extension directions are always at 45to the shear zone. For example, veins and foliations at the edge of a shear zone, and P and T axes of earthquakes ([Section 8.3.2](ch08-3-2)) are all oriented 45to the shear (fault) zone ({numref}`Figure %s <ch08_fig04>`).

# ```{figure} /figures/ch08_fig04.png
# :width: 600px
# :name: ch08_fig04
# 
# The infinitesimal extension and shortening directions are always at 45$^o$ to the shear zone. **a.** Sigmoidal veins, **b.** Foliation, and **c.** P and T axes of earthquakes. Modified from Allmendinger et al. (2012).
# ```

# (ch08-3-2)=
# ### Applications of Infinitesimal Strain
# 
# (ch08-3-2-1)=
# #### **P** and **T** axes
# 
# Our first application is based on the observation we just made regarding the orientation of the principal axes of infinitesimal strain with respect to a shear zone ({numref}`Figure %s <ch08_fig04>`). In a fault, the infinitesimal shortening **P** and extension **T** axes are at 45$^o$ to the slip vector (e.g. fault striae), and they lie on the plane defined by the fault striae and the fault pole. This plane is known as the movement plane $M$ ({numref}`Figure %s <ch08_fig05>`).

# ```{figure} /figures/ch08_fig05.png
# :width: 600px
# :name: ch08_fig05
# 
# **a.** Reverse fault, striae and movement plane $M$. **b.** Representation of **a** on a lower hemisphere equal area stereonet. Modified from Marshak and Mitra (1988).
# ```

# The displacement gradient tensor, $e_{ij}$ of a population of $n$ small faults with poles vectors, $p_i$, and slip vectors, $u_i$, is given by (Allmendinger et al., 2012):

# $$
# e_{ij}=\frac{\sum_{i=1}^{n}\left(M_{g} u_{i} p_{j}\right)}{V}
# $$ (ch08_eq24)

# where $M_g$ is the geometric moment and $V$ is the volume of the region being deformed. As discussed in [Section 8.3](ch08.3), we can decompose this equation to yield the infinitesimal strain and rotation tensors:

# $$
# e_{i j}=\varepsilon_{i j}+\omega_{i j}=\frac{M_{g}\left(u_{i} p_{j}+u_{j} p_{i}\right)}{2 V}+\frac{M_{g}\left(u_{i} p_{j}-u_{j} p_{i}\right)}{2 V}
# $$ (ch08_eq25)

# Because $M_g/2V$ is a scalar, the orientations of the principal axes of infinitesimal strain are identical to the principal axes of the symmetric tensor $(u_i p_j+u_j p_i)$, which are just a function of the fault planes and striae orientations. The **P** and **T** axes are therefore a simple, direct representation of fault geometry and the orientation and sense of slip (Allmendinger et al., 1989).
# 
# The function [pt_axes](https://github.com/nfcd/compGeo/blob/master/source/functions/pt_axes.py) computes and plots the **P** and **T** axes from the orientation of several faults and their slip vectors:

# In[2]:


import numpy as np
import matplotlib.pyplot as plt

from compgeo.cart_to_sph import cart_to_sph
from compgeo.zero_twopi import zero_twopi
from compgeo.sph_to_cart import sph_to_cart
from compgeo.stereonet import stereonet
from compgeo.great_circle import great_circle
from compgeo.st_coord_line import st_coord_line
from compgeo.pole import pole_from_plane


def pt_axes(fault,slip,sense, fpsv):
    """
    pt_axes computes the P and T axes from the orientation
    of several fault planes and their slip vectors. Results
    are plotted in an equal area stereonet

    USE: P,T,senseC,fig,ax = pt_axes(fault,slip,sense)

    fault = nfaults x 2 vector with strikes and dips of
        faults
    slip = nfaults x 2 vector with trend and plunge of
        slip vectors
    sense = nfaults x 1 vector with sense of faults
    fpsv = A flag to tell wether the fault plane and
        slip vector are plotted (1) or not
    P = nfaults x 2 vector with trend and plunge of P axes
    T = nfaults x 2 vector with trend and plunge of T axes
    senseC = nfaults x 1 vector with corrected sense of slip

    fig and ax are handles to the figure and axes

    NOTE: Input/Output angles are in radians

    Python function based on the Matlab function
    PTAxes in Allmendinger et al. (2012)
    """
    pi = np.pi
    # Initialize some vectors
    p = np.zeros(3)
    u = np.zeros(3)
    eps = np.zeros((3,3))
    P = np.zeros((np.size(fault,0),2))
    T = np.zeros((np.size(fault,0),2))
    senseC = sense

    # For all faults
    for i in range(np.size(fault,0)):
        
        # Direction cosines of pole to fault and slip vector
        trd, plg = pole_from_plane(fault[i,0],fault[i,1])
        p[0],p[1],p[2] = sph_to_cart(trd, plg)
        u[0],u[1],u[2] = sph_to_cart(slip[i,0],slip[i,1])
        # Compute u(i)*p(j) + u(j)*p(i)
        for j in range(3):
            for k in range(3):
                eps[j,k]=u[j]*p[k]+u[k]*p[j]
        
        # Compute orientations of principal axes of strain
        # Here we use the function eigh
        _,V = np.linalg.eigh(eps)
        
        # P orientation
        P[i,0],P[i,1] = cart_to_sph(V[0,2],V[1,2],V[2,2])
        if P[i,1] < 0:
            P[i,0] = zero_twopi(P[i,0]+pi)
            P[i,1] *= -1
        
        # T orientation
        T[i,0],T[i,1] = cart_to_sph(V[0,0],V[1,0],V[2,0]) 
        if T[i,1] < 0.0:
            T[i,0] = zero_twopi(T[i,0]+pi)
            T[i,1] *= -1
        
        # Determine 3rd component of pole cross product slip
        cross = p[0] * u[1] - p[1] * u[0]
        
        # Use cross and first character in sense to
        # determine if kinematic axes should be switched
        s2 = "p"
        if sense[i][0] == "T" or sense[i][0] == "t": 
            s2 = "Y"
        
        if (sense[i][0]=="R" or sense[i][0]=="r") and cross>0.0:
            s2 = "Y"
        
        if (sense[i][0]=="L" or sense[i][0]=="l") and cross<0.0: 
            s2 = "Y"
        
        if s2 == "Y":
            temp1 = P[i,0]
            temp2 = P[i,1]
            P[i,0] = T[i,0]
            P[i,1] = T[i,1]
            T[i,0] = temp1
            T[i,1] = temp2
            
            if cross < 0.0: 
                senseC[i] = "TL"
            
            if cross > 0.0:
                senseC[i] = "TR"
        else:
            
            if cross < 0.0:
                senseC[i] = "NR"
            
            if cross > 0.0:
                senseC[i] = "NL"

    # Plot in equal area stereonet
    fig, ax = stereonet(0,90*pi/180,10*pi/180,1)
    
    # Plot P and T axes
    for i in range(np.size(fault,0)):
        
        if fpsv == 1:
            
            # Plot fault
            path = great_circle(fault[i,0],fault[i,1],1)
            ax.plot(path[:,0],path[:,1],"k")
            
            # Plot slip vector (black circle)
            xp,yp = st_coord_line(slip[i,0],slip[i,1],1)
            ax.plot(xp,yp,"ko","MarkerFaceColor","k")
        
        # Plot P axis (blue circle)
        xp,yp = st_coord_line(P[i,0],P[i,1],1)
        ax.plot(xp,yp,"bo","MarkerFaceColor","b")
        
        # Plot T axis (red circle)
        xp,yp = st_coord_line(T[i,0],T[i,1],1)
        ax.plot(xp,yp,"ro","MarkerFaceColor","r")

    return P, T, senseC, fig, ax


# Let’s use this function to plot the **P** and **T** axes for a group of faults from the Central Andes in Northern Argentina. The file [jujuy.txt](https://github.com/nfcd/compGeo/blob/master/source/data/ch8-1/jujuy.txt) contains the orientation (strike and dip) of the faults, and the orientation (trend and plunge) and sense of movement of the slip vectors (striae).
# 
# :::{note}
# The file [jujuy.txt](https://github.com/nfcd/compGeo/blob/master/source/data/ch8-1/jujuy.txt) is the same one included in the program [FaultKin](https://www.rickallmendinger.net/faultkin) by Richard Allmendinger.
# :::
# 
# The notebook [ch8-1](https://github.com/nfcd/compGeo/blob/master/source/notebooks/ch8-1.ipynb) shows how to do this:

# In[3]:


# Import libraries
import os
import numpy as np
pi = np.pi

# Import function pt_axes
from compgeo.pt_axes import pt_axes


# Read the faults
jujuy = np.loadtxt(os.path.abspath("data/ch8-1/jujuy.txt"),
                   usecols = [0,1,2,3])

fault = jujuy[:,0:2] * pi/180
slip = jujuy[:,2:4] * pi/180

sense = np.loadtxt(os.path.abspath("data/ch8-1/jujuy.txt"), 
                   usecols = 4, dtype = "str")

# Compute P and T axes and plot them
# Don't plot the faults and slip vectors
P, T, senseC, fig, ax = pt_axes(fault, slip, sense, 0)


# Here, the overall shortening direction given by the **P** axes (blue dots) is about E-W. From the **P** and **T** axes, it is possible to calculate a moment tensor summation and the kinematic axes, which define two nodal planes separating the regions of extension (**T** axes) and shortening (**P** axes). The regions of extension are typically shaded, and so the diagram looks like a beach ball ({numref}`Figure %s <ch08_fig04>`c). Beach ball diagrams are also used to display the focal mechanisms of earthquakes, where P-waves first arrivals pushing the ground up are marked as **T** axes, and those pushing the ground down are marked as **P** axes. We leave the moment tensor summation for the Exercises section.

# (ch08-3-2-2)=
# #### 2D strain from GPS data
# 
# The global positioning system (GPS) has revolutionized Earth sciences by providing geoscientists with real time monitoring of active deformation. Continuous GPS measurements provide mm resolution of the displacement of stations. Because the changes in distance between stations is very small (10s of mm) relative to the distance between stations (10s of km), the strains measured by GPS networks are infinitesimal. {numref}`Figure %s <ch08_fig06>` shows this problem in two dimensions:

# ```{figure} /figures/ch08_fig06.png
# :width: 450px
# :name: ch08_fig06
# 
# Three points in an initial configuration $\textbf{X}$ move along non-parallel vectors $\textbf{u}$ to a final configuration $\textbf{x}$, resulting in strain. Modified from Allmendinger et al. (2012).
# ```

# If the strain is homogeneous, the displacement of the stations is expressed by the following equation:

# $$
# u_{i}=t_{i}+e_{i j} X_{j}
# $$ (ch08_eq26)

# where $t_i$ is the translation vector and $e_{ij}$ is the displacement gradient tensor. In matrix form, this equation can be written as:

# $$
# \left[\begin{array}{c}{ }^{1} u_1 \\ { }^{1} u_{2} \\ { }^{2} u_{1} \\ { }^{2} u_{2} \\ \cdots \\ \cdots \\ { }^{n} u_{1} \\ { }^{n} u_{2}\end{array}\right]=\left[\begin{array}{cccccc}1 & 0 & { }^{1} X_{1} & { }^{1} X_{2} & 0 & 0 \\ 0 & 1 & 0 & 0 & { }^{1} X_{1} & { }^{1} X_{2} \\ 1 & 0 & { }^{2} X_{1} & { }^{2} X_{2} & 0 & 0 \\ 0 & 1 & 0 & 0 & { }^{2} X_{1} & { }^{2} X_{2} \\ \cdots & \cdots & \cdots & \cdots & \cdots & \cdots \\ \cdots & \cdots & \cdots & \cdots & \cdots & \cdots \\ 1 & 0 & { }^{n} X_{1} & { }^{n} X_{2} & 0 & 0 \\ 0 & 1 & 0 & 0 & { }^{n} X_{1} & { }^{n} X_{2}\end{array}\right]\left[\begin{array}{l}t_{1} \\ t_{2} \\ e_{11} \\ e_{12} \\ e_{21} \\ e_{22}\end{array}\right]
# $$ (ch08_eq27)

# where the superscripts 1 to $n$ refer to the stations. The column vector to the right of Eq. {eq}`ch08_eq27` contains the unknowns, which are the two components of the translation vector ($t_1$ and $t_2$) and the four components of the displacement gradient tensor ($e_{11}$, $e_{12}$, $e_{21}$ and $e_{22}$). Therefore, in 2D there are 6 unknowns, and since each station delivers 2 equations, we need a minimum of 3 non-colinear stations to determine the strain ellipse.
# 
# :::{note}
# In 3D there are 12 unknowns and each station delivers 3 equations. Therefore, we need a minimum of 4 non-colinear stations to determine the strain ellipsoid.
# :::
# 
# Notice that Eq. {eq}`ch08_eq27` is written not just for 3 stations but for $n$ stations. If there are more than 3 stations in 2D (or 4 stations in 3D), the system is overdetermined (more equations than unknowns), and we can use the extra information to assess the uncertainties of the model parameters.
# 
# The solution to Eq. {eq}`ch08_eq27` is a classical application of inverse theory, specifically the solution of the linear least-squares problem (Press et al., 1986). This problem has the form:

# $$
# \mathbf{y}=\mathbf{M} \mathbf{x}
# $$ (ch08_eq28)

# where **y** is the vector with known displacements, **M** is the matrix with the location of the stations (this matrix is known as the design matrix), and **x** is the vector with the unknowns. To solve for **x**, **y** is multiplied by the inverse of **M**:

# $$
# \mathbf{x}=\mathbf{M}^{-1} \mathbf{y}
# $$ (ch08_eq29)

# In Matlab or Python, there are routines specifically designed to solve this problem. We use the function [lscov](https://github.com/nfcd/compGeo/blob/master/source/functions/lscov.py) by Paul Müller to solve this problem.
# 
# There are three main strategies to compute the strain from a network of stations: i. To calculate the strain on triangular cells whose vertices are defined by the stations (Delaunay triangulation), ii. To calculate the strain on a regular grid of cells, using the stations within a radius, $r$, from the center of each cell (nearest neighbor method), or iii. To calculate the strain on a regular grid of cells, using all the stations weighted by their distance to the center of each cell (distance weighted method; Cardozo and Allmendinger; 2009). In this last case, the weighting factor is given by:

# $$
# W=\exp \left[\frac{-d^{2}}{2 \alpha^{2}}\right]
# $$ (ch08_eq30)

# where $d$ is the distance from the station to the center of the cell, and $\alpha$ is a constant that specifies how the impact of a station decays with distance. A larger value of $\alpha$ smooths out local variations.
# 
# The function [grid_strain](https://github.com/nfcd/compGeo/blob/master/source/functions/grid_strain.py) computes and plots the infinitesimal strain of a network of GPS stations. Notice that after computing the displacement gradient for each cell, we use the function `inf_strain` to compute the strain:

# In[4]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from scipy.spatial import Delaunay

from compgeo.lscov import lscov
from compgeo.inf_strain import inf_strain


def grid_strain(pos,disp,k,par,plotpar,plotst):
    """
    grid_strain computes the infinitesimal strain of a network
    of stations with displacements in x (east) and y (north).
    Strain in z is assumed to be zero (plane strain)

    USE: cent,eps,ome,pstrain,rotc,fig,ax = 
        grid_strain(pos,disp,k,par,plotpar,plotst)

    pos = nstations x 2 matrix with x (east) and y (north)
        positions of stations in meters
    disp = nstations x 2 matrix with x (east) and y (north)
        displacements of stations in meters
    k = Type of computation: Delaunay (k = 0), nearest
        neighbor (k = 1), or distance weighted (k = 2)
    par = Parameters for nearest neighbor or distance
        weighted computation. If Delaunay (k = 0), enter
        a scalar corresponding to the minimum internal
        angle of a triangle valid for computation.
        If nearest neighbor (k = 1), input a 1 x 3 vector
        with grid spacing, number of nearest neighbors,
        and maximum distance to neighbors. If distance
        weighted (k = 2), input a 1 x 2 vector with grid
        spacing and distance weighting factor alpha
    plotpar = Parameter to color the cells: Max elongation
        (plotpar = 0), minimum elongation
        (plotpar = 1), rotation (plotpar = 2),
        or dilatation (plotpar = 3)
    plotst = A flag to plot the stations (1) or not (0)
    cent = ncells x 2 matrix with x and y positions of 
        cells centroids
    eps = 3 x 3 x ncells array  with strain tensors of
        the cells
    ome = 3 x 3 x ncells array with rotation tensors of
        the cells
    pstrain = 3 x 3 x ncells array with magnitude and
        orientation of principal strains of
        the cells
    rotc = ncells x 3 matrix with rotation components
        of cells

    fig and ax are handles to the figure and axes

    NOTE: Input/Output angles are in radians. Output
        azimuths are given with respect to North
        pos, disp, grid spacing, max. distance to
        neighbors, and alpha should be in meters

    Python function translated from the Matlab function
    GridStrain in Allmendinger et al. (2012)
    """
    pi = np.pi
    
    # If Delaunay
    if k == 0:
        # Indexes of triangles vertices
        # Use function Delaunay
        tri = Delaunay(pos)
        inds = tri.simplices
        
        # Number of cells
        ncells = np.size(inds,0)
        
        # Number of stations per cell = 3
        nstat = 3
        
        # Centers of cells
        cent = np.zeros((ncells,2))
        for i in range(ncells):
            
            # Triangle vertices
            v1x=pos[inds[i,0],0]
            v2x=pos[inds[i,1],0]
            v3x=pos[inds[i,2],0]
            v1y=pos[inds[i,0],1]
            v2y=pos[inds[i,1],1]
            v3y=pos[inds[i,2],1]
            
            # Center of cell
            cent[i,0]=(v1x + v2x + v3x)/3.0
            cent[i,1]=(v1y + v2y + v3y)/3.0
            
            # Triangle internal angles
            s1 = np.sqrt((v3x-v2x)**2 + (v3y-v2y)**2)
            s2 = np.sqrt((v1x-v3x)**2 + (v1y-v3y)**2)
            s3 = np.sqrt((v2x-v1x)**2 + (v2y-v1y)**2)
            a1 = np.arccos((v2x-v1x)*(v3x-v1x)/(s3*s2)+\
                (v2y-v1y)*(v3y-v1y)/(s3*s2))
            a2 = np.arccos((v3x-v2x)*(v1x-v2x)/(s1*s3)+\
                (v3y-v2y)*(v1y-v2y)/(s1*s3))
            a3 = np.arccos((v2x-v3x)*(v1x-v3x)/(s1*s2)+\
                (v2y-v3y)*(v1y-v3y)/(s1*s2))
            
            # If any of the internal angles is less than
            # specified minimum, invalidate triangle
            if a1 < par or a2 < par or a3 < par:
                inds[i,:] = np.zeros(3)
    
    # If nearest neighbor or distance weighted
    else:
        
        # Construct grid
        xmin, xmax = min(pos[:,0]), max(pos[:,0])
        ymin, ymax = min(pos[:,1]), max(pos[:,1])
        cellsx = int(np.ceil((xmax-xmin)/par[0]))
        cellsy = int(np.ceil((ymax-ymin)/par[0]))
        xgrid = np.arange(xmin,(xmin+(cellsx+1)*par[0]),par[0])
        ygrid = np.arange(ymin,(ymin+(cellsy+1)*par[0]),par[0])
        XX,YY = np.meshgrid(xgrid,ygrid)
        
        # Number of cells
        ncells = cellsx * cellsy
        
        # Number of stations per cell (nstat) and
        # other parameters
        # If nearest neighbor
        if k == 1:
            nstat = par[1] # max neighbors
            sqmd = par[2]**2 # max squared distance
        
        # If distance weighted
        elif k == 2:
            nstat = np.size(pos,0) # all stations
            dalpha = 2.0*par[1]*par[1] # 2*alpha*alpha
        
        # Cells" centers
        cent = np.zeros((ncells,2))
        count = 0
        for i in range(cellsy):
            for j in range(cellsx):
                cent[count,0] = (XX[i,j]+XX[i,j+1])/2.0
                cent[count,1] = (YY[i,j]+YY[i+1,j])/2.0
                count += 1
        
        # Initialize stations indexes for cells to -1
        inds = np.ones((ncells,nstat), dtype=int)*-1
        
        # Initialize weight matrix for distance weighted
        wv = np.zeros((ncells,nstat*2))
        
        # For all cells set stations indexes
        for i in range(ncells):
            
            # Initialize sq distances to -1.0
            sds = np.ones(nstat)*-1.0
            
            # For all stations
            for j in range(np.size(pos,0)):
                
                # sq distance from cell center to station
                dx = cent[i,0] - pos[j,0]
                dy = cent[i,1] - pos[j,1]
                sd = dx**2 + dy**2
                
                # If nearest neighbor
                if k == 1:
                    
                    # If within the max sq distance
                    if sd <= sqmd:
                        minsd = min(sds)
                        mini = np.argmin(sds)
                        
                        # If less than max neighbors
                        if minsd == -1.0:
                            sds[mini] = sd
                            inds[i,mini] = j
                        
                        # If max neighbors
                        else:
                            
                            # If sq distance is less 
                            # than neighbors max sq distance
                            maxsd = max(sds)
                            maxi = np.argmax(sds)
                            if sd < maxsd:
                                sds[maxi] = sd
                                inds[i,maxi] = j
                
                # If distance weighted
                elif k == 2:
                    
                    # All stations indexes
                    inds[i,:] = np.arange(nstat)
                    
                    # Weight factor
                    weight = np.exp(-sd/dalpha)
                    wv[i,j*2] = weight
                    wv[i,j*2+1] = weight

    # Initialize arrays
    y = np.zeros(nstat*2)
    M = np.zeros((nstat*2,6))
    e = np.zeros((3,3)) 
    eps = np.zeros((3,3,ncells))
    ome = np.zeros((3,3,ncells)) 
    pstrain = np.zeros((3,3,ncells))
    rotc = np.zeros((ncells,3)) 

    # For each cell
    for i in range(ncells):
        
        # If required minimum number of stations
        if min(inds[i,:]) >= 0:
            
            # Displacements column vector y
            # and design matrix M. X1 = North, X2 = East
            for j in range(nstat):
                ic = inds[i,j]
                y[j*2] = disp[ic,1]
                y[j*2+1] = disp[ic,0]
                M[j*2,:] = [1.,0.,pos[ic,1],pos[ic,0],0.,0.]
                M[j*2+1,:] = [0.,1.,0.,0.,pos[ic,1],pos[ic,0]]
            
            # Find x using function lscov
            # If Delaunay or nearest neighbor
            if k == 0 or k == 1:
                x = lscov(M,y)
            
            # If distance weighted
            elif k == 2:
                x = lscov(M,y,wv[i,:])
            
            # Displacement gradient tensor
            for j in range(2):
                e[j,0] = x[j*2+2]
                e[j,1] = x[j*2+3]
            
            # Compute strain
            eps[:,:,i],ome[:,:,i],pstrain[:,:,i],\
                rotc[i,:],_ = inf_strain(e)

    # Variable to plot
    # If maximum principal strain
    if plotpar == 0:
        vp = pstrain[0,0,:]
        lcb = "emax"
    
    # If minimum principal strain
    elif plotpar == 1:
        vp = pstrain[2,0,:]
        lcb = "emin"
    
    # If rotation: 
    # For plane strain, rotation = rotc(3)
    elif plotpar == 2:
        vp = rotc[:,2]*180/pi
        lcb = r"Rotation ($\circ$)"
    
    # If dilatation
    elif plotpar == 3:
        vp = pstrain[0,0,:]+pstrain[1,0,:]+pstrain[2,0,:]
        lcb = "dilatation"

    # Make a figure
    fig, ax = plt.subplots(figsize=(15,7.5))

    # Patches and colors for cells
    patches = []
    colors = []

    # Fill cells patches and colors
    # If Delaunay
    if k == 0:
        for i in range(ncells):
            
            # If minimum number of stations
                if min(inds[i,:]) >= 0:
                    xpyp = [[pos[inds[i,0],0],pos[inds[i,0],1]],\
                        [pos[inds[i,1],0],pos[inds[i,1],1]],\
                        [pos[inds[i,2],0],pos[inds[i,2],1]]]
                    
                    # length in km
                    xpyp = np.divide(xpyp,1e3)
                    polygon = Polygon(xpyp, True)
                    patches.append(polygon)
                    colors.append(vp[i])
    
    # If nearest neighbor or distance weighted
    if k == 1 or k == 2:
        count = 0
        for i in range(cellsy):
            for j in range(cellsx):
                
                # If minimum number of stations
                if min(inds[count,:]) >= 0:
                    xpyp = [[XX[i,j],YY[i,j]],[XX[i,j+1],YY[i,j+1]],\
                        [XX[i+1,j+1],YY[i+1,j+1]],[XX[i+1,j],YY[i+1,j]]]
                    
                    # length in km
                    xpyp = np.divide(xpyp,1e3)
                    polygon = Polygon(xpyp, True)
                    patches.append(polygon)
                    colors.append(vp[count])
                count += 1

    # Collect cells patches
    pcoll = PatchCollection(patches)
    
    # Cells colors
    pcoll.set_array(np.array(colors))
    
    # Color map is blue to red
    pcoll.set_cmap("bwr")
    
    # Positive values are red, negative are
    # blue and zero is white
    vmin = min(vp) 
    vmax = max(vp)
    norm=mcolors.TwoSlopeNorm(vmin=vmin,vcenter=0.0,vmax=vmax)
    pcoll.set_norm(norm)

    # Draw cells
    ax.add_collection(pcoll)

    # Plot stations
    if plotst == 1:
        ax.plot(pos[:,0]*1e-3,pos[:,1]*1e-3,"k.",markersize=2) 

    # Axes
    ax.axis("equal")
    ax.set_xlabel("x (km)")
    ax.set_ylabel("y (km)")

    # Color bar with nice ticks
    intv = (vmax-vmin)*0.25
    ticks=[vmin,vmin+intv,vmin+2*intv,vmin+3*intv,vmax]
    lticks = ["{:.2e}".format(ticks[0]),\
        "{:.2e}".format(ticks[1]),"{:.2e}".format(ticks[2]),\
        "{:.2e}".format(ticks[3]),"{:.2e}".format(ticks[4])]
    cbar = fig.colorbar(pcoll, label=lcb, ticks=ticks)
    cbar.ax.set_yticklabels(lticks)

    return cent, eps, ome, pstrain, rotc, fig, ax


# Let’s use this function to compute the infinitesimal rotation in Tibet and the Himalaya. The file [tibet.txt](https://github.com/nfcd/compGeo/blob/master/source/data/ch8-2/tibet.txt) contains the UTM, east and west, coordinates of GPS stations in the Tibetan Plateau region and their displacements in meters (from Zhang et al., 2004). The notebook [ch8-2](https://github.com/nfcd/compGeo/blob/master/source/notebooks/ch8-2.ipynb) shows the solution to this problem using the Delaunay, nearest neighbor, and distance weighted methods.

# In[5]:


import numpy as np
pi = np.pi

# Import function grid_strain
from compgeo.grid_strain import grid_strain


# Load Zhang et al. GPS data from the Tibetan plateau
# load x, y coordinates and displacements
tibet = np.loadtxt(os.path.abspath("data/ch8-2/tibet.txt"))
pos = tibet[:,0:2]
disp = tibet[:,2:4] 

# Rotation from Delaunay triangulation, plot stations
par = 10 * pi/180 #Minimum internal angle of triangles
cent, eps, ome, pstrain, rotc, fig, ax = grid_strain(pos, disp, 0, par, 2, 1)


# In[6]:


# Rotation from nearest neighbor 
# Grid spacing = 75 km, neighbors 6, 
# max. distance = 150 km, plot stations
par = [75e3,6,150e3]
cent, eps, ome, pstrain, rotc, fig, ax = grid_strain(pos, disp, 1, par, 2, 1)


# In[7]:


# Rotation from distance Weighted 
# Grid spacing = 75 km, alpha = 150 km, plot stations
par = [75e3,150e3]
cent, eps, ome, pstrain, rotc, fig, ax = grid_strain(pos, disp, 2, par, 2, 1)


# Given the uneven distribution of the stations, the distance weighted method is perhaps the best representation of the infinitesimal strain. The rotation is about a downward vertical axis. The GPS-based rotation shows large coherent domains of clockwise (positive) and counterclockwise (negative) rotation. Interestingly enough, these domains are consistent with the permanent, long-term deformation of this region (Allmendinger et al., 2007).

# (ch08-4)=
# ## Finite strain
# 
# When deformations are large, the initial and final states are not identical:

# $$
# d X_{i} \neq d x_{i} \quad \text { and } \quad \frac{\partial u_{i}}{\partial X_{i}} \neq \frac{\partial u_{i}}{\partial x_{i}}
# $$ (ch08_eq31)

# and therefore, we need to use all the tensors in {numref}`ch08_tab01`.
# 
# The mathematics of finite strain is quite involved (e.g. Means, 1976; Allmendinger et al., 2012) and we will not cover it in detail here. Rather, we will focus on the main tensors required to compute finite strain. We start with the finite strain tensor in the undeformed configuration:

# $$
# E_{i j}=\frac{1}{2}\left[\frac{\partial u_{i}}{\partial X_{j}}+\frac{\partial u_{j}}{\partial X_{i}}+\frac{\partial u_{k}}{\partial X_{i}} \frac{\partial u_{k}}{\partial X_{j}}\right]=\frac{1}{2}\left[e_{i j}+e_{j i}+e_{k i} e_{k j}\right]
# $$ (ch08_eq32)

# where $E_{ij}$ is known as the *Lagrangian finite strain tensor*. Similarly, we can compute the finite strain tensor in the deformed configuration:

# $$
# \bar{E}_{i j}=\frac{1}{2}\left[\frac{\partial u_{i}}{\partial x_{j}}+\frac{\partial u_{j}}{\partial x_{i}}-\frac{\partial u_{k}}{\partial x_{i}} \frac{\partial u_{k}}{\partial x_{j}}\right]=\frac{1}{2}\left[\bar{e}_{i j}+\bar{e}_{j i}-\bar{e}_{k i} \bar{e}_{k j}\right]
# $$ (ch08_eq33)

# where $\bar{E}_{i j}$ is known as the *Eulerian finite strain tensor.
# 
# :::{note}
# We use an upper bar to denote tensors in the deformed configuration, e.g. $\bar{e}_{ij}$ and $\bar{E}_{ij}$.
# :::
# 
# Notice that Eqs. {eq}`ch08_eq29` and {eq}`ch08_eq30` are similar to Eq. {eq}`ch08_eq15` for the infinitesimal strain tensor, but with an extra term ($e_{k i} e_{k j}$ or $\bar{e}_{k i} \bar{e}_{k j}$).
# 
# From the deformation gradient $\mathbf{F}$, we can also compute deformation tensors. In the undeformed configuration:

# $$
# C_{i j}=\frac{\partial x_{k}}{\partial X_{i}} \frac{\partial x_{k}}{\partial X_{j}}=F_{k i} F_{k j}
# $$ (ch08_eq34)

# where $C_{ij}$ is known as the *Green deformation tensor*. Similarly, for the deformed configuration:

# $$
# \bar{C}_{i j}=\frac{\partial X_{k}}{\partial x_{i}} \frac{\partial X_{k}}{\partial x_{j}}=\bar{F}_{k i} \bar{F}_{k j}
# $$ (ch08_eq35)

# where $\bar{C}_{i j}$ is known as the *Cauchy deformation tensor*.
# 
# The finite strain and deformation tensors are related. In the undeformed configuration:

# $$
# E_{i j}=\frac{1}{2}\left(C_{i j}-\delta_{i j}\right) \quad \text { and } \quad C_{i j}=\delta_{i j}+2 E_{i j}
# $$ (ch08_eq36)

# and in the deformed configuration:

# $$
# \bar{E}_{i j}=\frac{1}{2}\left(\delta_{i j}-\bar{C}_{i j}\right) \quad \text { and } \quad \bar{C}_{i j}=\delta_{i j}-2 \bar{E}_{i j}
# $$ (ch08_eq37)

# where $\delta_{ij}$ is the Kronecker delta, a function that returns 1 if $i = j$, or 0 otherwise.
# 
# Thus, the finite strain tensors do not contain any more information than the deformation tensors, and viceversa. They are all symmetric tensors that have principal axes and can be represented by Mohr circles (see next section). From Eqs. {eq}`ch08_eq36` and {eq}`ch08_eq37`, it is clear that $E_{ij}$ and $C_{ij}$ have the same principal axes orientations. This is also the case for $\bar{E}_{ij}$ and $\bar{C}_{ij}$. In addition these tensors are related to the quadratic elongation, $\lambda_i$, along the coordinate axes:

# $$
# \lambda_{i}=C_{i i}=1+2 E_{i i}
# $$ (ch08_eq38)

# and:

# $$
# \frac{1}{\lambda_{i}}=\bar{C}_{i i}=1-2 \bar{E}_{i i}
# $$ (ch08_eq39)

# {numref}`ch08_tab02` shows the finite strain and deformation tensors for the deformation in {numref}`Figure %s <ch08_fig02>`.

# ```{list-table} Finite strain and deformation tensors for Figure 8.2.
# :header-rows: 0
# :name: ch08_tab02
# 
# * - Undeformed
#   - $E_{ij}=\left[\begin{array}{lll}0 & 0 & 0 \\ 0 & -0.375 & 0 \\ 0 & 0 & 1.5\end{array}\right]$
#   - $C_{ij} = \left[\begin{array}{lll}1 & 0 & 0 \\ 0 & 0.25 & 0 \\ 0 & 0 & 4\end{array}\right]$
# * - Deformed
#   - $\bar{E}_{ij}=\left[\begin{array}{lll}0 & 0 & 0 \\ 0 & -1.5 & 0 \\ 0 & 0 & 0.375\end{array}\right]$
#   - $\bar{C}_{ij} = \left[\begin{array}{lll}1 & 0 & 0 \\ 0 & 4 & 0 \\ 0 & 0 & 0.25\end{array}\right]$
# ```

# You can see that the diagonal components of the Green deformation tensor, $C_{ij}$, are the quadratic elongations, $\lambda_i$, along the coordinate axes, and the diagonal components of the Cauchy deformation tensor, $\bar{C}_{ij}$, are the inverse of the quadratic elongations, $1/\lambda_i$. The diagonal components of the Lagrangian and Eulerian finite strain tensors, $E_{ij}$ and $\bar{E}_{ij}$, can be found from Eqs. {eq}`ch08_eq38` and {eq}`ch08_eq39`.
# 
# The function [fin_strain](https://github.com/nfcd/compGeo/blob/master/source/functions/fin_strain.py) computes the finite strain from the displacement gradient tensor, in the undeformed (`frame` = 0) or deformed (`frame` = 1) configuration. Notice that for computing the maximum shear strain, the function assumes plane strain (Ramsay, 1967).

# In[8]:


import numpy as np

from compgeo.cart_to_sph import cart_to_sph


def fin_strain(e,frame):
    """
    fin_strain computes finite strain from an input
    displacement gradient tensor

    USE: eps,pstrain,dilat,maxsh = fin_strain(e,frame)

    e = 3 x 3 Lagrangian or Eulerian displacement gradient
        tensor
    frame = Reference frame. 0 = undeformed (Lagrangian)
        state, 1 = deformed (Eulerian) state
    eps = 3 x 3 Lagrangian or Eulerian strain tensor
    pstrain = 3 x 3 matrix with magnitude (column 1), trend
        (column 2) and plunge (column 3) of maximum
        (row 1), intermediate (row 2), and minimum
        (row 3) elongations
    dilat = dilatation
    maxsh = 1 x 2 vector with max. shear strain and
        orientation with respect to maximum principal
        strain direction. Only valid in 2D

    NOTE: Output angles are in radians

    Python function translated from the Matlab function
    FinStrain in Allmendinger et al. (2012)
    """
    # Initialize variables
    eps = np.zeros((3,3))
    pstrain = np.zeros((3,3))
    maxsh = np.zeros(2)

    # Compute strain tensor
    for i in range(3):
        for j in range(3):
            eps[i,j] = 0.5*(e[i][j]+e[j][i])
            for k in range(3):
                
                # If undeformed reference frame: 
                # Lagrangian strain tensor
                if frame == 0:
                    eps[i,j] = eps[i,j] + 0.5*(e[k][i]*e[k][j])
                
                # If deformed reference frame: 
                # Eulerian strain tensor
                elif frame == 1:
                    eps[i,j] = eps[i,j] - 0.5*(e[k][i]*e[k][j])

    # Compute principal elongations and orientations
    D, V = np.linalg.eigh(eps)

    # Principal elongations
    for i in range(3):
        ind = 2-i
        
        # Magnitude
        # If undeformed reference frame: 
        # Lagrangian strain tensor
        if frame == 0:
            pstrain[i,0] = np.sqrt(1.0+2.0*D[ind])-1.0
        
        # If deformed reference frame:
        # Eulerian strain tensor
        elif frame == 1:
            pstrain[i,0] = np.sqrt(1.0/(1.0-2.0*D[ind]))-1.0
        
        # Orientations
        pstrain[i,1],pstrain[i,2] = cart_to_sph(V[0,ind],
            V[1,ind],V[2,ind])

    # Dilatation
    dilat = (1.0+pstrain[0,0])*(1.0+pstrain[1,0])* \
        (1.0+pstrain[2,0]) - 1.0

    # Maximum shear strain: This only works if plane strain
    lmax = (1.0+pstrain[0,0])**2 # Maximum quad. elongation
    lmin = (1.0+pstrain[2,0])**2 # Minimum quad. elongation
    
    # Maximum shear strain: Ramsay (1967) Eq. 3.46
    maxsh[0] = (lmax-lmin)/(2.0*np.sqrt(lmax*lmin))
    
    # Angle of maximum shear strain with respect to maximum
    # principal strain. Ramsay (1967) Eq. 3.45
    # If undeformed reference frame
    if frame == 0:
        maxsh[1] = np.pi/4.0
    
    # If deformed reference frame
    elif frame == 1:
        maxsh[1] = np.arctan(np.sqrt(lmin/lmax))

    return eps, pstrain, dilat, maxsh


# (ch08-4-1)=
# ### Mohr circle for finite strain
# 
# The finite strain and deformation tensors are symmetric tensors. Therefore, a rotation about one of the principal axis of the tensor can be represented by a Mohr circle. As you may suspect, there are Mohr circles for finite strain in the undeformed and deformed configuration (Ramsay, 1967). For geoscientists, the Mohr circle for finite strain in the deformed configuration is the most important, since in nature we do observe deformed rocks.
# 
# To derive this Mohr circle, we can start with the Cauchy deformation tensor in a principal axes coordinate system:

# $$
# \bar{C}_{i j}=\left[\begin{array}{ccc}\bar{C}_{1} & 0 & 0 \\ 0 & \bar{C}_{2} & 0 \\ 0 & 0 & \bar{C}_{3}\end{array}\right]
# $$ (ch08_eq40)

# and perform a rotation $\theta$ about $\mathbf{X_2}$, which is described by the transformation matrix:

# $$
# a_{ij} = \begin{pmatrix}\cos\theta&0&\sin\theta\\0&1&0\\-\sin\theta&0&\cos\theta\end{pmatrix}
# $$ (ch08_eq41)

# The tensor transformation equation is:

# $$
# \bar{C}_{ij}' = a_{ik}a_{jl}\bar{C}_{kl}
# $$ (ch08_eq42)

# which results in the new form of the tensor:

# $$
# \bar{C}'_{ij}=\begin{bmatrix}\bar{C}_1\cos^2\theta+\bar{C}_3\sin^2\theta&0&-(\bar{C}_1-\bar{C}_3)\sin\theta\cos\theta\\0&\bar{C}_2&0\\-(\bar{C}_1-\bar{C}_3)\sin\theta\cos\theta&0&\bar{C}_1\sin^2\theta+\bar{C}_3\cos^2\theta\end{bmatrix}
# $$ (ch08_eq43)

# Upon rearranging, we get:

# $$
# \begin{gathered}
#     \bar{C}_{11}'=\frac{\bar{C}_1+\bar{C}_3}{2}+\frac{\bar{C}_1-\bar{C}_3}{2}\cos 2\theta  \\
#     \bar{C}_{33}'=\frac{\bar{C}_1+\bar{C}_3}{2}-\frac{\bar{C}_1-\bar{C}_3}{2}\cos 2\theta \\
#     \bar{C}_{13}'=\bar{C}_{31}'=-\frac{\bar{C}_1-\bar{C}_3}{2}\sin 2\theta
# \end{gathered}
# $$ (ch08_eq44)

# As stated in Eq. {eq}`ch08_eq39`, $\bar{C}_{ii} = 1/\lambda_i$. Also $\bar{C}_{ij}$ for $i \neq j$ is equal to $\gamma / \lambda$. Using $\lambda ' = 1/\lambda$ and $\gamma ' = \gamma / \lambda$, we get the equations for the Mohr circle for finite strain in the deformed configuration:

# $$
# \begin{gathered}
#     \lambda^{\prime}=\frac{\left(\lambda_{1}^{\prime}+\lambda_{3}^{\prime}\right)}{2}+\frac{\left(\lambda_{1}^{\prime}-\lambda_{3}^{\prime}\right)}{2} \cos 2 \theta  \\
#     \gamma^{\prime}=-\frac{\left(\lambda_{1}^{\prime}-\lambda_{3}^{\prime}\right)}{2} \sin 2 \theta
# \end{gathered}
# $$ (ch08_eq45)

# As an example, {numref}`Figure %s <ch08_fig07>`a shows a deformed circle and an inscribed triangle after 30$^o$ clockwise shear. {numref}`Figure %s <ch08_fig07>`b shows the Mohr circle for this deformation. In the Mohr circle, the horizontal axis is $\lambda'$, and the vertical axis is $\gamma'$. We follow the convention of Ragan (2009): Positive angular shear, $\psi$, corresponds to anticlockwise rotation of the original line’s normal ({numref}`Figure %s <ch08_fig07>`a), and the $\gamma'$ axis in the Mohr circle increases downwards ({numref}`Figure %s <ch08_fig07>`b). Notice that the angle between the horizontal axis and a line from the origin to any point in the Mohr circle is equal to $\psi$. Therefore, a line from the origin and tangent to the Mohr circle indicates $\psi_{\text{max}}$ ({numref}`Figure %s <ch08_fig07>`b).

# ```{figure} /figures/ch08_fig07.png
# :width: 600px
# :name: ch08_fig07
# 
# **a.** Physical plane, and **b.** Mohr circle for 30° clockwise shear of a circle and a triangle. The pole, the triangle bisectors a, b and c (red and blue), the principal strain axes $S_{1}$ and $S_{3}$ (black), the LNFE (green), and the lines of maximum shear strain (orange) are all drawn in the Mohr circle.
# ```

# From points in the Mohr circle, one can trace the corresponding lines with the same orientation than in the physical plane ({numref}`Figure %s <ch08_fig07>`b, triangle bisectors a, b and c). These lines will intersect at a point called the *pole* to the Mohr circle ({numref}`Figure %s <ch08_fig07>`b). From the pole, one can trace lines of any orientation; they will intersect the circle at points that represent the strain of lines with the same orientation in the physical plane. Thus, you can imagine the pole to be the center of the strain ellipse, and from it, it is easy to trace the principal axes of strain ($S_1$ and $S_3$), the lines of no finite elongation (LNFE, $\lambda'=1$), and the lines of maximum shear strain ({numref}`Figure %s <ch08_fig07>`b).

# (ch08-4-2)=
# ### 2D finite strain from displacement data
# 
# If we have a group of points or stations with displacement data, we can determine the finite strain following a strategy similar to the one we used for infinitesimal strain. The function [grid_fin_strain](https://github.com/nfcd/compGeo/blob/master/source/functions/grid_fin_strain.py) computes and plots the finite strain for a group of points with displacement data. This function is very similar to our previous function `grid_strain`, and therefore we only include its header here. Notice that after computing the displacement gradient in the undeformed (`frame` = 0) or deformed (`frame` = 1) configuration, we use our function `fin_strain` to compute the finite strain in each cell.

# In[9]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from scipy.spatial import Delaunay

from compgeo.lscov import lscov
from compgeo.fin_strain import fin_strain


def grid_fin_strain(pos,disp,frame,k,par,plotpar,plotst):
    """
    grid_fin_strain computes the finite strain of a group
    of points with displacements in x and y.
    Strain in z is assumed to be zero (plane strain)
    
    USE: cent,eps,pstrain,dilat,maxsh,fig,ax = 
        grid_fin_strain(pos,disp,frame,k,par,plotpar,plotst)
    
    pos = npoints x 2 matrix with x and y position
        of points
    disp = nstations x 2 matrix with x and y
        displacements of points
    frame = Reference frame. 0 = undeformed (Lagrangian)
        state, 1 = deformed (Eulerian) state
    k = Type of computation: Delaunay (k = 0), nearest
        neighbor (k = 1), or distance weighted (k = 2)
    par = Parameters for nearest neighbor or distance
        weighted computation. If Delaunay (k = 0), enter
        a scalar corresponding to the minimum internal
        angle of a triangle valid for computation.
        If nearest neighbor (k = 1), input a 1 x 3 vector
        with grid spacing, number of nearest neighbors,
        and maximum distance to neighbors. If distance
        weighted (k = 2), input a 1 x 2 vector with grid
        spacing and distance weighting factor alpha
    plotpar = Parameter to color the cells: Max elongation
        (plotpar = 0), minimum elongation
        (plotpar = 1), dilatation (plotpar = 2),
        or max. shear strain (plotpar = 3)
    plotst = A flag to plot the stations (1) or not (0)
    cent = ncells x 2 matrix with x and y positions of the
        cells centroids
    eps = 3 x 3 x ncells array    with strain tensors of
        the cells
    pstrain = 3 x 3 x ncells array with magnitude and
        orientation of principal strains of the cells
    dilat = ncells x 1 vector with dilatation of the cells
    maxsh = ncells x 2 matrix with max. shear strain and
        orientation with respect to maximum principal
        strain direction, of the cells.
        Only valid for plane strain
        
    fig and ax are handles to the figure and axes
    
    NOTE: Input/Output angles are in radians. Output
        azimuths are given with respect to y
        pos, disp, grid spacing, max. distance to
        neighbors, and alpha should be in the same
        length units
    """
    # If Delaunay
    if k == 0:
        # Indexes of triangles vertices
        # Use function Delaunay
        tri = Delaunay(pos)
        inds = tri.simplices
        # Number of cells
        ncells = np.size(inds,0)
        # Number of stations per cell = 3
        nstat = 3
        # Centers of cells
        cent = np.zeros((ncells,2))
        for i in range(ncells):
            # Triangle vertices
            v1x=pos[inds[i,0],0]
            v2x=pos[inds[i,1],0]
            v3x=pos[inds[i,2],0]
            v1y=pos[inds[i,0],1]
            v2y=pos[inds[i,1],1]
            v3y=pos[inds[i,2],1]
            # Center of cell
            cent[i,0]=(v1x + v2x + v3x)/3.0
            cent[i,1]=(v1y + v2y + v3y)/3.0
            # Triangle internal angles
            s1 = np.sqrt((v3x-v2x)**2 + (v3y-v2y)**2)
            s2 = np.sqrt((v1x-v3x)**2 + (v1y-v3y)**2)
            s3 = np.sqrt((v2x-v1x)**2 + (v2y-v1y)**2)
            a1 = np.arccos((v2x-v1x)*(v3x-v1x)/(s3*s2)+\
                (v2y-v1y)*(v3y-v1y)/(s3*s2))
            a2 = np.arccos((v3x-v2x)*(v1x-v2x)/(s1*s3)+\
                (v3y-v2y)*(v1y-v2y)/(s1*s3))
            a3 = np.arccos((v2x-v3x)*(v1x-v3x)/(s1*s2)+\
                (v2y-v3y)*(v1y-v3y)/(s1*s2))
            # If any of the internal angles is less than
            # specified minimum, invalidate triangle
            if a1 < par or a2 < par or a3 < par:
                inds[i,:] = np.zeros(3)
    # If nearest neighbor or distance weighted
    else:
        # Construct grid
        xmin = min(pos[:,0]); xmax = max(pos[:,0])
        ymin = min(pos[:,1]); ymax = max(pos[:,1])
        cellsx = int(np.ceil((xmax-xmin)/par[0]))
        cellsy = int(np.ceil((ymax-ymin)/par[0]))
        xgrid = np.arange(xmin,(xmin+(cellsx+1)*par[0]),par[0])
        ygrid = np.arange(ymin,(ymin+(cellsy+1)*par[0]),par[0])
        XX,YY = np.meshgrid(xgrid,ygrid)
        # Number of cells
        ncells = cellsx * cellsy
        # Number of stations per cell (nstat) and
        # other parameters
        # If nearest neighbor
        if k == 1:
            nstat = par[1] # max neighbors
            sqmd = par[2]**2 # max squared distance
        # If distance weighted
        elif k == 2:
            nstat = np.size(pos,0) # all stations
            dalpha = 2.0*par[1]*par[1] # 2*alpha*alpha
        # Cells" centers
        cent = np.zeros((ncells,2))
        count = 0
        for i in range(cellsy):
            for j in range(cellsx):
                cent[count,0] = (XX[i,j]+XX[i,j+1])/2.0
                cent[count,1] = (YY[i,j]+YY[i+1,j])/2.0
                count += 1
        # Initialize stations indexes for cells to -1
        inds = np.ones((ncells,nstat), dtype=int)*-1
        # Initialize weight matrix for distance weighted
        wv = np.zeros((ncells,nstat*2))
        # For all cells set stations indexes
        for i in range(ncells):
            # Initialize sq distances to -1.0
            sds = np.ones(nstat)*-1.0
            # For all stations
            for j in range(np.size(pos,0)):
                # Sq distance from cell center to station
                dx = cent[i,0] - pos[j,0]
                dy = cent[i,1] - pos[j,1]
                sd = dx**2 + dy**2
                # If nearest neighbor
                if k == 1:
                    # If within the max sq distance
                    if sd <= sqmd:
                        minsd = min(sds)
                        mini = np.argmin(sds)
                        # If less than max neighbors
                        if minsd == -1.0:
                            sds[mini] = sd
                            inds[i,mini] = j
                        # If max neighbors
                        else:
                            # If sq distance is less
                            # than neighbors max sq distance
                            maxsd = max(sds)
                            maxi = np.argmax(sds)
                            if sd < maxsd:
                                sds[maxi] = sd
                                inds[i,maxi] = j
                # If distance weighted
                elif k == 2:
                    # All stations indexes
                    inds[i,:] = np.arange(nstat)
                    # Eq. 8.27: Weight factor
                    weight = np.exp(-sd/dalpha)
                    wv[i,j*2] = weight
                    wv[i,j*2+1] = weight
    
    # Initialize arrays
    y = np.zeros(nstat*2)
    M = np.zeros((nstat*2,6))
    e = np.zeros((3,3)) 
    eps = np.zeros((3,3,ncells)) 
    pstrain = np.zeros((3,3,ncells))
    dilat = np.zeros((ncells,1))
    maxsh = np.zeros((ncells,2)) 
        
    # For each cell
    for i in range(ncells):
        # If required minimum number of stations
        if min(inds[i,:]) >= 0:
            # Eq. 8.24: Displacements column vector y
            # and design matrix M. X1 = y, X2 = x
            for j in range(nstat):
                ic = inds[i,j]
                y[j*2] = disp[ic,1]
                y[j*2+1] = disp[ic,0]
                M[j*2,:] = [1.,0.,pos[ic,1],pos[ic,0],0.,0.]
                M[j*2+1,:] = [0.,1.,0.,0.,pos[ic,1],pos[ic,0]]
            # Eqs. 8.25-8.26: Find x using function lscov
            # If Delaunay or nearest neighbor
            if k == 0 or k == 1:
                x = lscov(M,y)
            # If distance weighted
            elif k == 2:
                x = lscov(M,y,wv[i,:])
            # Displacement gradient tensor
            for j in range(2):
                e[j,0] = x[j*2+2]
                e[j,1] = x[j*2+3]
            # Compute strain
            eps[:,:,i],pstrain[:,:,i],dilat[i,:],\
                maxsh[i,:] = fin_strain(e,frame)

    # Variable to plot
    # If maximum principal strain
    if plotpar == 0:
        vp = pstrain[0,0,:]
        vmin = 0.0
        vmax = 2.0
        lcb = "emax"
    # If minimum principal strain
    elif plotpar == 1:
        vp = pstrain[2,0,:]
        vmin = -2.0
        vmax = 0.0
        lcb = "emin"
    # If dilatation:
    elif plotpar == 2:
        vp = dilat[:]
        vmin = -1.0
        vmax = 1.0
        lcb = "dilatation"
    # If max. shear strain
    elif plotpar == 3:
        vp = maxsh[:,0]
        vmin = 0.0
        vmax = 2.0
        lcb = "max. shear strain"

    # Make a figure
    fig, ax = plt.subplots(figsize=(15.0,7.5))
    
    # Patches and colors for cells
    patches = []
    colors = []
    
    # Fill cells patches and colors
    # If Delaunay
    if k == 0:
        for i in range(ncells):
            # If minimum number of stations
                if min(inds[i,:]) >= 0:
                    xpyp = [[pos[inds[i,0],0],pos[inds[i,0],1]],\
                            [pos[inds[i,1],0],pos[inds[i,1],1]],\
                            [pos[inds[i,2],0],pos[inds[i,2],1]]]
                    polygon = Polygon(xpyp, True)
                    patches.append(polygon)
                    colors.append(vp[i])
    # If nearest neighbor or distance weighted
    if k == 1 or k == 2:
        count = 0
        for i in range(cellsy):
            for j in range(cellsx):
                # If minimum number of stations
                if min(inds[count,:]) >= 0:
                    xpyp = [[XX[i,j],YY[i,j]],[XX[i,j+1],YY[i,j+1]],\
                        [XX[i+1,j+1],YY[i+1,j+1]],[XX[i+1,j],YY[i+1,j]]]
                    polygon = Polygon(xpyp, True)
                    patches.append(polygon)
                    colors.append(vp[count])
                count += 1

    # Collect cells patches
    pcoll = PatchCollection(patches)
    # Cells colors
    pcoll.set_array(np.array(colors))
    # Color map is blue to red
    pcoll.set_cmap("jet")
    norm = mcolors.Normalize(vmin=vmin,vmax=vmax)
    pcoll.set_norm(norm)
    
    # Draw cells
    ax.add_collection(pcoll)
    
    # Plot stations
    if plotst == 1:
        ax.plot(pos[:,0],pos[:,1],"k.",markersize=2) 
    
    # Axes
    ax.axis("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    
    # Color bar with nice ticks
    intv = (vmax-vmin)*0.25
    ticks=[vmin,vmin+intv,vmin+2*intv,vmin+3*intv,vmax]
    lticks = ["{:.2}".format(ticks[0]),\
            "{:.2}".format(ticks[1]),"{:.2}".format(ticks[2]),\
            "{:.2}".format(ticks[3]),"{:.2}".format(ticks[4])]
    cbar = fig.colorbar(pcoll, label=lcb, ticks=ticks)
    cbar.ax.set_yticklabels(lticks)
    
    return cent, eps, pstrain, dilat, maxsh, fig, ax


# Let’s use this function to compute the maximum shear strain of a discrete element model of a normal fault. The discrete element method is a mechanical method that simulates the rocks as an assembly of elements. The file [demfault.txt](https://github.com/nfcd/compGeo/blob/master/source/data/ch8-3/demfault.txt) contains the deformed $x$ and $y$ coordinates of the elements and their displacements in meters. The notebook [ch8-3](https://github.com/nfcd/compGeo/blob/master/source/notebooks/ch8-3.ipynb) shows the solution to this problem using the nearest neighbor method:

# In[10]:


import os
import numpy as np

# Import function grid_fin_strain
from compgeo.grid_fin_strain import grid_fin_strain

# load x, y deformed coordinates and displacements
demfault = np.loadtxt(os.path.abspath("data/ch8-3/demfault.txt"))
pos = demfault[:,0:2]
disp = demfault[:,2:4]

# Max. shear strain from nearest neighbor, Def. config. 
# Grid spacing = 0.2 m, neighbors 6, 
# max. distance = 0.4 m, plot stations
par = [0.2,6,0.4]
cent,eps,pstrain,dilat,maxsh,fig,ax = grid_fin_strain(pos,
                                        disp,1,1,par,3,1)

# Add units to the axes
ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")
plt.show()


# The maximum shear strain delineates the internal structure of the fault. You can find more information about computing strain from displacement or velocity data in Cardozo and Allmendinger (2009).

# (ch08-5)=
# ## Progressive strain
# 
# We have focused so far on the undeformed (initial) and deformed (final) states. However, finite strain is actually the cumulative result of a series of strain increments:

# $$
# \begin{aligned}{ }^{1} x_{i} &={ }^{1} F_{i j}{ }^{1} X_{j} \\{ }^{2} X_{i} &={ }^{1} x_{i} \\{ }^{2} x_{i} &={ }^{2} F_{i j}{ }^{2} X_{j} \\{ }^{3} X_{i} &={ }^{2} x_{i} \\{ }^{3} x_{i} &={ }^{3} F_{i j}{ }^{3} X_{j} \\ \cdots \\{ }^{n} x_{i} &={ }^{n} F_{i j}{ }^{n} X_{j} \end{aligned}
# $$ (ch08_eq46)

# where $F_{ij}$ is the Green deformation gradient, and $1$ to $n$ are the strain increments. This can also be written as:

# $$
# x_{i}={ }^{n} F_{i j} \ldots{ }^{3} F_{i j}{ }^{2} F_{i j}{ }^{1} F_{i j}{ }^{1} X_{j}
# $$ (ch08_eq47)

# Notice that since ${ }^{2}\mathbf{F}{ }^{1}\mathbf{F}\neq{ }^{1}\mathbf{F}{ }^{2}\mathbf{F}$, for finite strain we must know the order of deformation. For example, if we want to determine the finite strain of a group of faults in a region, we must know the order at which these faults formed. This is often very difficult to determine.
# 
# Let’s look at some simple deformations, which are characterized by the same incremental deformation gradient through time. For simplicity, we will assume that there is no strain along the $\mathbf{X_2}$ axis, and all strain is in the $\mathbf{X_1}\mathbf{X_3}$ plane. Thus, we will be dealing with plane strain.

# (ch08-5-1)=
# ### Pure shear
# 
# For pure shear, the principal stretches are along the coordinate axes (e.g. {numref}`Figure %s <ch08_fig02>`). Let’s assume the maximum principal stretch, $S_1$, is along $\mathbf{X_1}$, and the minimum principal strecth, $S_3$, is along $\mathbf{X_3}$ ({numref}`Figure %s <ch08_fig08>`a). $S_2 = 1$ (plane strain) and is parallel to $\mathbf{X_2}$. The Green deformation gradient for this case is:

# $$
# { }^{PS}F_{ij}=\left[\begin{array}{ccc}S_{1} & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 0 & S_{3}\end{array}\right]
# $$ (ch08_eq48)

# ```{figure} /figures/ch08_fig08.png
# :width: 600px
# :name: ch08_fig08
# 
# **a.** Pure shear, and **b.** Simple shear. Blue is undeformed, and red is deformed geometry. Large white and gray circles are initial and final positions, respectively. The displacement paths of the points are divided in 10 increments. Modified from Allmendinger et al. (2012).
# ```

# The function [pure_shear](https://github.com/nfcd/compGeo/blob/master/source/functions/pure_shear.py) deforms a collection of points using pure shear, a value of $S_1$ (`st1`), and assuming plane strain and area conservation ($S_1S_3 = 1$). The function displays the points’ displacement paths for a number of increments (`ninc`), and the progressive strain history as the value of $S_1$ versus the angle $S_1$ makes with the $\mathbf{X_1}$ axis ($\Theta$). The last two parameters are computed from the eigenvalues and eigenvectors of the Green deformation tensor:

# In[11]:


import numpy as np
import matplotlib.pyplot as plt


def pure_shear(pts,st1,ninc):
    """
    pure_shear computes and plots displacement paths and
    progressive finite strain history for pure shear with
    maximum stretching parallel to the X1 axis

    USE: paths,pfs,fig,ax = pure_shear(pts,st1,ninc)

    pts: npoints x 2 matrix with X1 and X3 coord. of points
    st1 = Maximum principal stretch
    ninc = number of strain increments
    paths = displacement paths of points
    pfs = progressive finite strain history. column 1 =
        orientation of maximum stretch with respect to X1
        in degrees, column 2 = maximum stretch magnitude

    fig and ax are handles to the figure and axes

    NOTE: Intermediate principal stretch is 1.0 (Plane 
        strain). Output orientations are in radians

    Python function based on the Matlab function
    PureShear in Allmendinger et al. (2012)
    """
    # Compute minimum principal stretch and incr. stretches
    st1inc=st1**(1.0/ninc)
    st3=1.0/st1
    st3inc=st3**(1.0/ninc)

    # Initialize displacement paths
    npts = np.size(pts,0) # Number of points
    paths = np.zeros((ninc+1,npts,2))
    paths[0,:,:] = pts # Initial points of paths

    # Calculate incr. deformation gradient tensor
    F = np.array([[st1inc, 0.0], [0.0, st3inc]])

    # Initialize figure
    fig, ax = plt.subplots(1, 2, figsize=(15,5)) # 1 x 2 figure

    # Compute displacement paths
    for i in range(npts): # for all points
        for j in range(ninc+1): # for all strain increments
            for k in range(2):
                for L in range(2):
                    paths[j,i,k] = F[k,L]*paths[j-1,i,L] + paths[j,i,k]
        
        # Plot displacement path of point
        ax[0].plot(paths[:,i,0], paths[:,i,1], "k.-")

    # Plot initial polygon
    inpol = np.zeros((npts+1,2))
    inpol[0:npts,]=paths[0,0:npts,:]
    inpol[npts,] = inpol[0,]
    ax[0].plot(inpol[:,0],inpol[:,1],"b-")
    
    # Plot final polygon
    finpol = np.zeros((npts+1,2))
    finpol[0:npts,]=paths[ninc,0:npts,:]
    finpol[npts,] = finpol[0,]
    ax[0].plot(finpol[:,0],finpol[:,1],"r-")

    # set axes
    ax[0].set_xlabel(r"$\mathbf{X_1}$")
    ax[0].set_ylabel(r"$\mathbf{X_3}$")
    ax[0].grid()
    ax[0].axis("equal")

    # Initalize progressive finite strain history
    pfs = np.zeros((ninc+1,2))
    pfs[0,:] = [0, 1] #Initial state

    # Calculate progressive finite strain history
    for i in range(1,ninc+1):
        
        # Determine the finite deformation gradient tensor
        finF = np.linalg.matrix_power(F, i)
        
        # Determine Green deformation tensor
        G = np.dot(finF,finF.conj().transpose())
        
        # Stretch magnitude and orientation: Maximum 
        # eigenvalue and their corresponding eigenvectors
        # of Green deformation tensor
        D, V = np.linalg.eigh(G)
        pfs[i,0] = np.arctan(V[1,1]/V[0,1])
        pfs[i,1] = np.sqrt(D[1])

    # Plot progressive finite strain history
    ax[1].plot(pfs[:,0]*180/np.pi,pfs[:,1],"k.-")
    ax[1].set_xlabel(r"$\Theta\;(\circ)$")
    ax[1].set_ylabel("Maximum finite stretch")
    ax[1].set_xlim(-90,90)
    ax[1].set_ylim(1,max(pfs[:,1])+0.5)
    ax[1].grid()

    return paths, pfs, fig, ax


# Let’s use this function to reproduce {numref}`Figure %s <ch08_fig08>`a. The notebook [ch8-4](https://github.com/nfcd/compGeo/blob/master/source/notebooks/ch8-4.ipynb) shows how to do this:

# In[12]:


# Import libraries
import numpy as np

# Import function pure_shear
from compgeo.pure_shear import pure_shear

# Initial points coordinates
pts = np.zeros((16,2))
pts[:,0]=[-1,-1,-1,-1,-1,-0.5,0,0.5,1,1,1,1,1,0.5,0,-0.5]
pts[:,1]=[-1,-0.5,0,0.5,1,1,1,1,1,0.5,0,-0.5,-1,-1,-1,-1]
st1 = 2.5
ninc = 10

paths, psf, fig, ax = pure_shear(pts, st1, ninc)


# As you can see the principal strain axes, $S_1$ and $S_3$, do not rotate throughout the deformation ($\Theta$ is 0 throughout the deformation). Pure shear is a *non-rotational* deformation.

# (ch08-5-2)=
# ### Simple shear
# 
# For plane strain and simple shear along the $\mathbf{X_1}$ axis ({numref}`Figure %s <ch08_fig08>`b), the Green deformation gradient is:

# $$
# { }^{SS}F_{i j}=\left[\begin{array}{lll}1 & 0 & \gamma \\ 0 & 1 & 0 \\ 0 & 0 & 1\end{array}\right]
# $$ (ch08_eq49)

# where $\gamma$ is the shear strain.
# 
# The function [simple_shear](https://github.com/nfcd/compGeo/blob/master/source/functions/simple_shear.py) deforms a collection of points using simple shear and a value of shear strain (`gamma`). The function displays the points’ displacement paths for a number of increments (`ninc`), and the progressive strain history as the value of $S_1$ versus $\Theta$:

# In[13]:


import numpy as np
import matplotlib.pyplot as plt


def simple_shear(pts,gamma,ninc):
    """
    simple_shear computes and plots displacement paths and
    progressive finite strain history for simple shear
    parallel to the X1 axis

    USE: paths,pfs,fig,ax = simple_shear(pts,gamma,ninc)

    pts: npoints x 2 matrix with X1 and X3 coord. of points
    gamma = Engineering shear strain
    ninc = number of strain increments
    paths = displacement paths of points
    pfs = progressive finite strain history. column 1 =
        orientation of maximum stretch with respect to X1 
        in degrees, column 2 = maximum stretch magnitude

    fig and ax are handles to the figure and axes

    NOTE: Intermediate principal stretch is 1.0 (Plane 
        strain). Output orientations are in radians

    Python function based on the Matlab function
    SimpleShear in Allmendinger et al. (2012)
    """
    # Incremental engineering shear strain
    gammainc = gamma/ninc

    # Initialize displacement paths
    npts = np.size(pts,0) # Number of points
    paths = np.zeros((ninc+1,npts,2))
    paths[0,:,:] = pts # Initial points of paths

    # Calculate incr. deformation gradient tensor Eq. 8.44
    F = np.array([[1.0, gammainc],[0.0, 1.0]])

    # Initialize figure
    fig, ax = plt.subplots(1, 2, figsize=(15,5)) # 1 x 2 figure

    # Compute displacement paths
    for i in range(npts): # for all points
        for j in range(ninc+1): # for all strain increments
            for k in range(2):
                for L in range(2):
                    paths[j,i,k] = F[k,L]*paths[j-1,i,L] + paths[j,i,k]
        
        # Plot displacement path of point
        ax[0].plot(paths[:,i,0], paths[:,i,1], "k.-")

    # Plot initial polygon
    inpol = np.zeros((npts+1,2))
    inpol[0:npts,]=paths[0,0:npts,:]
    inpol[npts,] = inpol[0,]
    ax[0].plot(inpol[:,0],inpol[:,1],"b-")
    
    # Plot final polygon
    finpol = np.zeros((npts+1,2))
    finpol[0:npts,]=paths[ninc,0:npts,:]
    finpol[npts,] = finpol[0,]
    ax[0].plot(finpol[:,0],finpol[:,1],"r-")

    # set axes
    ax[0].set_xlabel(r"$\mathbf{X_1}$")
    ax[0].set_ylabel(r"$\mathbf{X_3}$")
    ax[0].grid()
    ax[0].axis("equal")

    # Initalize progressive finite strain history
    pfs = np.zeros((ninc+1,2))
    
    # In. state: Max. extension is at 45 deg from shear zone
    pfs[0,:] = [np.pi/4.0, 1.0]

    # Calculate progressive finite strain history
    for i in range(1,ninc+1):
        
        # Determine the finite deformation gradient tensor
        finF = np.linalg.matrix_power(F, i)
        
        # Determine Green deformation tensor
        G = np.dot(finF,finF.conj().transpose())
        
        # Stretch magnitude and orientation: Maximum 
        # eigenvalue and their corresponding eigenvectors
        # of Green deformation tensor
        D, V = np.linalg.eigh(G)
        pfs[i,0] = np.arctan(V[1,1]/V[0,1])
        pfs[i,1] = np.sqrt(D[1])

    # Plot progressive finite strain history
    ax[1].plot(pfs[:,0]*180/np.pi,pfs[:,1],"k.-")
    ax[1].set_xlabel(r"$\Theta\;(\circ)$")
    ax[1].set_ylabel("Maximum finite stretch")
    ax[1].set_xlim(-90,90)
    ax[1].set_ylim(1,max(pfs[:,1])+0.5)
    ax[1].grid()

    return paths, pfs, fig, ax


# Let’s use this function to reproduce {numref}`Figure %s <ch08_fig08>`b. The notebook [ch8-5](https://github.com/nfcd/compGeo/blob/master/source/notebooks/ch8-5.ipynb) shows how to do this:

# In[14]:


# Import libraries
import numpy as np

# Import function simple_shear
from compgeo.simple_shear import simple_shear

# Initial points coordinates
pts = np.zeros((16,2))
pts[:,0]=[-1,-1,-1,-1,-1,-0.5,0,0.5,1,1,1,1,1,0.5,0,-0.5]
pts[:,1]=[-1,-0.5,0,0.5,1,1,1,1,1,0.5,0,-0.5,-1,-1,-1,-1]
gamma = 2.5
ninc = 10
paths, psf, fig, ax = simple_shear(pts, gamma, ninc)


# In this case, the principal strain axes, $S_1$ and $S_3$, rotate throughout the deformation. $\Theta$ is initially 45(as predicted by infinitesimal strain, {numref}`Figure %s <ch08_fig04>`), and it progressively decreases to a value of 20by the end of the deformation. Simple shear is a *rotational* deformation.

# (ch08-5-3)=
# ### General shear
# 
# Sub-simple shear (De Paor, 1983), or general shear, is a combination of pure shear and simple shear. If the shear direction and $S_1$ are parallel to $\mathbf{X_1}$, $S_2 = 1$ (plane strain), and area is constant, the Green deformation gradient is (Allmendinger et al., 2012):

# $$
# { }^{GS}F_{i j}=\left[\begin{array}{lll}S_1 & 0 & \frac{\gamma(S_1-S_3)}{2\ln S_1} \\ 0 & 1 & \quad\: 0 \\ 0 & 0 &\quad\: S_3\end{array}\right]
# $$ (ch08_eq50)

# If the shear direction is parallel to $\mathbf{X_1}$, but $S_1$ is parallel to $\mathbf{X_3}$, the Green deformation gradient is (Allmendinger et al., 2012):

# $$
# { }^{GS}F_{i j}=\left[\begin{array}{lll}S_3 & 0 & \frac{\gamma(S_3-S_1)}{2\ln S_3} \\ 0 & 1 & \quad\: 0 \\ 0 & 0 &\quad\: S_1\end{array}\right]
# $$ (ch08_eq51)

# A simple dimensionless measure of the ratio of simple to pure shear is the cosine of the acute angle between the eigenvectors of the Green deformation gradient $\mathbf{F}$. This measure is known as the *kinematic vorticity number*, $W_k$ (Truesdell, 1953). General shear is characterized by a $W_k$ between 0 (pure shear) and 1 (simple shear).
# 
# The function [general_shear](https://github.com/nfcd/compGeo/blob/master/source/functions/general_shear.py) deforms a collection of points using general shear and values of $S_1$ (`st1`) and $\gamma$ (`gamma`), for a shear direction parallel to $\mathbf{X_1}$, and $S_1$ parallel (`kk` = 0) or perpendicular to the shear direction (`kk` = 1). Similar to the two previous functions, `general_shear` displays the points’ displacement paths for a number of increments (`ninc`), and $S_1$ versus $\Theta$:

# In[15]:


import numpy as np
import matplotlib.pyplot as plt


def general_shear(pts,st1,gamma,kk,ninc):
    """
    general_shear computes displacement paths, kinematic
    vorticity numbers and progressive finite strain 
    history, for general shear with a pure shear stretch,
    no area change, and a single shear strain

    USE: paths,wk,pfs,fig,ax = 
        general_shear(pts,st1,gamma,kk,ninc)

    pts = npoints x 2 matrix with X1 and X3 coord. of points
    st1 = Pure shear stretch parallel to shear zone
    gamma = Engineering shear strain
    kk = An integer that indicates whether the maximum 
        finite stretch is parallel (kk = 0), or 
        perpendicular (kk=1) to the shear direction
    ninc = number of strain increments
    paths = displacement paths of points
    wk = Kinematic vorticity number
    pfs = progressive finite strain history. column 1 =
        orientation of maximum stretch with respect to 
        X1, column 2 = maximum stretch magnitude

    fig and ax are handles to the figure and axes

    NOTE: Intermediate principal stretch is 1.0 (Plane
        strain). Output orientations are in radians

    Python function translated from the Matlab function
    GeneralShear in Allmendinger et al. (2012)
    """
    # Compute minimum principal stretch and incr. stretches
    st1inc =st1**(1.0/ninc)
    st3 =1.0/st1
    st3inc =st3**(1.0/ninc)

    # Incremental engineering shear strain
    gammainc = gamma/ninc

    # Initialize displacement paths
    npts = np.size(pts,0) # Number of points
    paths = np.zeros((ninc+1,npts,2))
    paths[0,:,:] = pts # Initial points of paths

    # Initialize figure
    fig, ax = plt.subplots(1, 2, figsize=(15,5)) # 1 x 2 figure

    # Calculate incremental deformation gradient tensor
    # If max. stretch parallel to shear direction Eq. 8.45
    if kk == 0:
        F=np.zeros((2,2))
        F[0,]=[st1inc, (gammainc*(st1inc-st3inc))/
            (2.0*np.log(st1inc))]
        F[1,]=[0.0, st3inc]
    
    # If max. stretch perpendicular to shear direction Eq. 8.46
    elif kk == 1:
        F=np.zeros((2,2))
        F[0,]= [st3inc, (gammainc*(st3inc-st1inc))/
                            (2.0*np.log(st3inc))]
        F[1,]= [0.0, st1inc]

    # Compute displacement paths
    for i in range(npts): # for all points
        for j in range(ninc+1): # for all strain increments
            for k in range(2):
                for L in range(2):
                    paths[j,i,k] = F[k,L]*paths[j-1,i,L] + paths[j,i,k]
        
        # Plot displacement path of point
        xx = paths[:,i,0]
        yy = paths[:,i,1]
        ax[0].plot(xx,yy,"k.-")

    # Plot initial and final polygons
    inpol = np.zeros((npts+1,2))
    inpol[0:npts,]=paths[0,0:npts,:]
    inpol[npts,] = inpol[0,]
    ax[0].plot(inpol[:,0],inpol[:,1],"b-")
    finpol = np.zeros((npts+1,2))
    finpol[0:npts,]=paths[ninc,0:npts,:]
    finpol[npts,] = finpol[0,]
    ax[0].plot(finpol[:,0],finpol[:,1],"r-")

    # Set axes
    ax[0].set_xlabel(r"$\mathbf{X_1}$")
    ax[0].set_ylabel(r"$\mathbf{X_3}$")
    ax[0].grid()
    ax[0].axis("equal")

    # Determine the eigenvectors of the flow (apophyses)
    # Since F is not symmetrical, use function eig
    _,V = np.linalg.eig(F)
    theta2 = np.arctan(V[1,1]/V[0,1])
    wk = np.cos(theta2)

    # Initalize progressive finite strain history. 
    # We are not including the initial state
    pfs = np.zeros((ninc,ninc))

    # Calculate progressive finite strain history
    for i in range(1,ninc+1):
        
        # Determine the finite deformation gradient tensor
        finF = np.linalg.matrix_power(F, i)
        
        # Determine Green deformation tensor
        G = np.dot(finF,finF.conj().transpose())
        
        # Stretch magnitude and orientation: Maximum 
        # eigenvalue and their corresponding eigenvectors
        # of Green deformation tensor
        D, V = np.linalg.eigh(G)
        pfs[i-1,0] = np.arctan(V[1,1]/V[0,1])
        pfs[i-1,1] = np.sqrt(D[1])

    # Plot progressive finite strain history
    ax[1].plot(pfs[:,0]*180/np.pi,pfs[:,1],"k.-")
    ax[1].set_xlabel(r"$\Theta\;(\circ)$")
    ax[1].set_ylabel("Maximum finite stretch")
    ax[1].set_xlim(-90,90)
    ax[1].set_ylim(1,max(pfs[:,1])+0.5)
    ax[1].grid()

    return paths, wk, pfs, fig, ax


# Let’s use this function to simulate the deformation produced by $S_1 = 2.0$ and $\gamma = 0.5$, for a shear direction parallel to $\mathbf{X_1}$, and $S_1$ parallel or perpendicular to $\mathbf{X_1}$. The notebook [ch8-6](https://github.com/nfcd/compGeo/blob/master/source/notebooks/ch8-6.ipynb) shows the solution to this problem:

# In[16]:


# Import libraries
import numpy as np

# Import function general_shear
from compgeo.general_shear import general_shear

# Initial points coordinates
pts = np.zeros((16,2))
pts[:,0]=[-1,-1,-1,-1,-1,-0.5,0,0.5,1,1,1,1,1,0.5,0,-0.5]
pts[:,1]=[-1,-0.5,0,0.5,1,1,1,1,1,0.5,0,-0.5,-1,-1,-1,-1]
st1 = 2.0
gamma = 0.5
ninc = 10

# Max. finite stretch parallel to shear direction
kk = 0
paths1, wk1, psf1, fig, ax = general_shear(pts, st1, gamma, kk, ninc)
print("Wk = {:.4f}".format(wk1))


# In[17]:


# Max. finite stretch perpendicular to shear direction
kk = 1
paths2, wk2, pfs2, fig, ax = general_shear(pts, st1, gamma, kk, ninc)
print("Wk = {:.4f}".format(wk2))


# In this case, the contribution of simple shear is about one third that of pure shear ($W_k = 0.34$). Try running the notebook with different values of $S_1$ and $\gamma$, and check $W_k$ and $\Theta$.

# (ch08-6)=
# ## Exercises
# 
# 1. This exercise is from Rowland and Duebendorfer (1994). The file [bar.txt](https://github.com/nfcd/compGeo/blob/master/source/data/ch8-exercise1/bar.txt) contains the strike and dip (RHR), slip direction, and slip sense of faults within the Basin and Range province, specifically southern Nevada and the Lake Mead area. Plot the **P** and **T** axes for these data using the function `pt_axes`. Are these data consistent with the regional E-W extension direction in the Basin and Range province?
# 
# 2. From the **P** and **T** axes, one can compute an unweighted moment tensor summation as follows (Allmendinger et al., 1989):
# 
#     $$
#     \scriptstyle \mathbf{K} = \left[\begin{array}{ccc}\scriptstyle\sum (P_N)^{2} - (T_N)^{2}  & \scriptstyle\sum (P_N) (P_E) - (T_N) (T_E) & \scriptstyle\sum (P_N) (P_D) - (T_N) (T_D) \\ \scriptstyle\sum (P_E) (P_N) - (T_E) (T_N) & \scriptstyle\sum (P_E)^{2} - (T_E)^{2} & \scriptstyle\sum (P_E) (P_D) - (T_E) (T_D) \\ \scriptstyle\sum (P_D) (P_N) - (T_D) (T_N) & \scriptstyle\sum (P_D) (P_E) - (T_D) (T_E) & \scriptstyle\sum (P_D)^{2} - (T_D)^{2}\end{array}\right]
#     $$ (ch08_eq52)
# 
#     where the subscripts $N$, $E$, and $D$ refer to the east, north and  down direction cosines of the **P** and **T** axes. The eigenvalues and eigenvectors of $\mathbf{K}$ give the relative magnitudes and orientations of the kinematic axes. These kinematic axes define the orientation of two nodal planes, or possible faults, separating the areas of infinitesimal extension (**T** axes) from those of infinitesimal shortening (**P** axes). The faults intersect at the intermediate eigenvector, the minimum and maximum eigenvectors define the movement plane, and the faults slip vectors are on the movement plane at 45 from the minimum and maximum eigenvectors. Modify the function `pt_axes` to plot the nodal planes. This is not a simple problem, don’t give up.
# 
# 3. The file [andes.txt](https://github.com/nfcd/compGeo/blob/master/source/data/ch8-exercise2/andes.txt) contains GPS data (stations UTM coordinates and displacements in meters) from the Central Andes (Kendrick et al., 2001; Brooks et al., 2003). Compute the infinitesimal rotation in this region using the function `grid_strain`. What does the rotation tell you about the tectonics of this region? Check Allmendinger et al. (2005) to support your answer.
# 
# 4. The file [antithetic.txt](https://github.com/nfcd/compGeo/blob/master/source/data/ch8-exercise3/antithetic.txt) contains data from a discrete element model of a low angle normal fault (deformed elements’ coordinates and displacements). Use the function `grid_fin_strain` to compute the shear strain of this model.
# 
# 5. The file [trilobite.txt](https://github.com/nfcd/compGeo/blob/master/source/data/ch8-exercise4/trilobite.txt) contains the undeformed $x$ and $y$ coordinates of points delineating a trilobite specimen. Deform the trilobite using **a.** pure shear with $S_1 = 2.0$, **b.** simple shear with $\gamma = 1.0$, and **c.** general shear with $S_1 = 1.5$, $\gamma = 1.0$, and shear direction and $S_1$ parallel to $x$.
# 
# :::{hint}
# Use functions `pure_shear`, `simple_shear`, and `general_shear`.
# :::

# (ch08-7)=
# ## References
# 
# Allmendinger, R.W., Gephart, J.W. and Marrett, R.A. 1989. Notes on fault
# slip analysis. GSA short course.
# 
# Allmendinger, R.W., Smalley, R.J., Bevis, M. et al. 2005. Bending the
# Bolivian orocline in real time. Geology 33, 905-908.
# 
# Allmendinger, R.W., Reilinger, R. and Loveless, J. 2007. Strain and
# rotation rate from GPS in Tibet, Anatolia, and the Altiplano. Tectonics
# 26, TC3013.
# 
# Allmendinger, R.W., Cardozo, N. and Fisher, D.M. 2012. Structural
# Geology Algorithms: Vectors and Tensors. Cambridge University Press.
# 
# Brooks, B. A., Bevis, M., R. Smalley, J. et al. 2003. Crustal Motion in
# the Southern Andes (26-36S): do the Andes behave like a microplate?
# Geochemistry, Geophysics, Geosystems - G3 4, GC000505.
# 
# Cardozo, N. and Allmendinger, R.W. 2009. SSPX: A program to compute
# strain from displacement/velocity data. Computers and Geosciences 35,
# 1343-1357.
# 
# DePaor, D.G. 1983. Ortographic analysis of geological structures - I.
# Deformation theory. Journal of Structural Geology 5, 255-277.
# 
# Kendrick, E., Bevis, M., Smalley, R. and Brooks, B. 2001. An integrated
# crustal velocity field for the Central Andes. Geochemistry, Geophysics,
# Geosystems - G3 2, GC000191.
# 
# Marshak, S. and Mitra, G. 1988. Basic Methods of Structural Geology.
# Prentice Hall.
# 
# Means, W.D. 1976. Stress and Strain: Basic Concepts of Continuum
# Mechanics for Geologists. New York: Springer-Verlag.
# 
# Press, W.H., Flannery, B.P., Teoukolsky, S.A. and Vetterling, W.T. 1986.
# Numerical Recipes: The Art of Scientific Computing. Cambridge University
# Press.
# 
# Ragan, D.M. 2009. Structural Geology: An Introduction to Geometrical
# Techniques. Cambridge University Press.
# 
# Ramsay, J.G. 1967. Folding and Fracturing of Rocks. McGraw-Hill, New
# York.
# 
# Rowland, S.M. and Duebendorfer, E.M. 1994. Structural Analysis and
# Synthesis: A Laboratory Course in Structural Geology. Wiley.
# 
# Truesdell, C. 1953. Two measures of vorticity. Journal of Rational
# Mechanics and Analysis 2, 173-217.
# 
# Zhang, P., Zhengkang, S., Min, W., et al. 2004, Continuous deformation
# of the Tibetan Plateau from Global Positioning System data: Geology 32,
# 809-812.
