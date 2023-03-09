#!/usr/bin/env python
# coding: utf-8

# (ch09)=
# # Elasticity

# We have covered stress and strain separately, but now we should consider how stress and strain are related in rocks. This relationship depends on the rock mechanical properties, which themselves depend on physical conditions such as the state of stress, temperature, and strain rate. The mathematical models that describe the relationship between stress and strain are known as constitutive models. In this chapter, we will cover the simplest constitutive model, namely elasticity. Because rocks in the upper crust fracture (deform non-elastically) at even modest strains, elasticity is closely related to the concepts of infinitesimal strain ([Section 8.3](ch08-3)), and it is ideal to solve infinitesimal-strain problems such as the propagation of seismic waves through the Earth. However, elasticity is also useful to understand other geological problems such as the processes leading to rock fractures, the flexure of the lithosphere under crustal loads, and even fault related folding, which can be envisioned as the summation of small deformation, elastic increments through time.

# (ch09-1)=
# ## Theory
# 
# In or study of elasticity, we will make some basic assumptions (Gudmunsson, 2011):
# 
# - The rock is homogeneous; this means that its properties do not vary with location.
# 
# - The rock is isotropic; this means that its properties are the same regardless of direction.
# 
# - The relation between stress and strain in the rock is linear and it is described by an elastic modulus ({numref}`Figure %s <ch09_fig01>`a).
# 
# - The strains are infinitesimal ([Section 8.3](ch08-3)).
# 
# - When the rock is stressed it instantaneously becomes strained. When the stress is removed, the strain disappears immediately (double arrows in {numref}`Figure %s <ch09_fig01>`a).

# ```{figure} /figures/ch09_fig01.png
# :width: 550px
# :name: ch09_fig01
# 
# **a.** Linear-elastic, versus **b.** Actual strain-stress behavior of a rock sample under compression. Modified from Gudmunsson (2011) and Zoback (2010).
# ```

# Thus, elastic deformation is linear and independent of both, the loading path and time. This facilitates the solution of elastic problems, as we will see in [Section 9.2](ch09-2). Since elasticity is path-independent, it is also possible to run elastic models forwards and backwards, which is convenient for solving inverse problems ([Chapter 10](ch10)). However, rocks are not linear elastic materials, and their behavior is more complex than that ({numref}`Figure %s <ch09_fig01>`b). Upon compressional loading in the lab, a rock sample experiences initially inelastic crack closing, followed by elastic behavior at very small strains, non-elastic permanent deformation, and finally failure (Zoback, 2010; {numref}`Figure %s <ch09_fig01>`b). Elasticity captures just a small portion of that.
# 
# For the basic assumptions above, the elastic model is expressed by the following equation in indicial format:

# $$
# \label{eq9.1}
# \varepsilon_{ij}=\frac{1+\nu}{E}\sigma_{ij}-\frac{\nu}{E}\sigma_{kk}\delta_{ij}
# $$ (ch09_eq01)

# where $\boldsymbol{\varepsilon}$ and $\boldsymbol{\sigma}$ are the strain and stress tensors, $E$ is the Young’s modulus, $\nu$ is the Poisson’s ratio, and $\delta_{ij}$ is the Kronecker delta (1 if $i = j$, 0 otherwise). This equation is also known as the Hooke’s law.
# 
# $E$ and $\nu$ become more clear if we consider a cube of rock under uniaxial stress ($\sigma_{11}=\sigma_1$, all other tractions = 0; {numref}`Figure %s <ch09_fig02>`a). Under this condition and from Eq. {eq}`ch09_eq01`, the strains are:

# $$
# \begin{gathered}
#     \varepsilon_{11}=\frac{1}{E}\sigma_{11} \\
#     \varepsilon_{22} = \varepsilon_{33} = \frac{-\nu}{E}\sigma_{11} = -\nu\varepsilon_{11}
# \end{gathered}
# $$ (ch09_eq02)

# ```{figure} /figures/ch09_fig02.png
# :width: 600px
# :name: ch09_fig02
# 
# Cube of rock under **a.** Uniaxial stress, **b.** Simple shear, and **c.** Triaxial stress. Blue is initial and red is final configuration. Formulas show the elastic moduli $E$, $\nu$, $G$ and $K$. Modified from Zoback (2010).
# ```

# Therefore, $E$ is the ratio of the axial stress to the axial strain, $\sigma_{11}/\varepsilon_{11}$. Since strain is adimensional, $E$ has units of stress. In consolidated rocks, $E$ can vary greatly from $\sim$ 1 to 100 GPa (Gudmunsson, 2011, his Appendix D).
# 
# $\nu$ is the negative ratio of the lateral strain to the axial strain, $-\varepsilon_{33}/\varepsilon_{11}$. $\nu$ is adimensional, and it can have values between 0 and 0.5. A nearly incompressible material such as rubber has a $\nu \sim 0.5$. Rubber maintains a nearly constant volume by laterally expanding just enough to compensate for a given shortening. A perfectly compressible material such as cork has a $\nu \sim 0.0$. Cork can be shortened without much lateral expanding. Rocks are somewhat compressible, they have a $\nu \sim 0.25$ (Pollard and Fletcher, 2005).
# 
# :::{note}
# A nice analogy for $\nu$ is a bottle stopper. A cork stopper can be cylindrical, while a rubber stopper must be tapered.
# :::
# 
# Now, let’s look at the case of simple shear (all tractions 0 except $\sigma_{13}$, {numref}`Figure %s <ch09_fig02>`b). In this case, the strain is given by:

# $$
# \varepsilon_{13}=\frac{1+\nu}{E}\sigma_{13}
# $$ (ch09_eq03)

# which can be written in terms of the shear strain $\gamma_{13}$ as ([Section 8.3](ch08-3)):

# $$
# \label{eq9.4}
# \begin{gathered}
#     \gamma_{13}=\frac{2(1+\nu)}{E}\sigma_{13}\quad \\
#     \frac{\sigma_{13}}{\gamma_{13}}=\frac{E}{2(1+\nu)}=G
# \end{gathered}
# $$ (ch09_eq04)

# $G$ is called the shear modulus (the ratio of the shear stress to the shear strain). For perfectly compressible material, $\nu=0$ so $G=E/2$, and for incompressible material, $\nu=0.5$ so $G=E/3$. For a rock with $\nu=0.25$, $G=E/2.5$.
# 
# Finally, let’s look at the more general case of triaxial stress (principal stress coordinate system and non-zero principal stresses, {numref}`Figure %s <ch09_fig02>`c). In this case, the axial strains are:

# $$\begin{gathered}
#     \varepsilon_{11}=\frac{1}{E}\sigma_{11}-\frac{\nu}{E}\sigma_{22}-\frac{\nu}{E}\sigma_{33} \\
#     \varepsilon_{22}=-\frac{\nu}{E}\sigma_{11}+\frac{1}{E}\sigma_{22}-\frac{\nu}{E}\sigma_{33} \\
#     \varepsilon_{33}=-\frac{\nu}{E}\sigma_{11}-\frac{\nu}{E}\sigma_{22}+\frac{1}{E}\sigma_{33}
# \end{gathered}
# $$ (ch09_eq05)

# To visualize this deformation in terms of volumetric changes, we can compute the mean stress, $\sigma_m$ ([Section 7.3](ch07-3)), and the volumetric strain, $\Delta$ (also called dilatation), as shown in {numref}`Figure %s <ch09_fig02>`c. The ratio of the mean stress to the volumetric strain, $\sigma_m/\Delta$, is called the bulk modulus, $K$, and it is given by:

# $$
# K=\frac{E}{3(1-2\nu)}
# $$ (ch09_eq06)

# The bulk modulus, $K$, approaches an infinite value as $\nu$ approaches 0.5, i.e. as the material becomes incompressible. For a perfectly compressible material, $\nu=0$ and $K=E/3$. Liquids are nearly incompressible, whereas gases are highly compressible. High porosity rocks are somewhat compressible (low $K$), whereas low porosity rocks tend to be less compressible (high $K$).
# 
# $E$, $\nu$, $G$ and $K$ control the propagation of seismic waves in the Earth ([Section 9.2.4](ch09-2-4)). If one can measure the velocity of compressional, $V_p$, and shear, $V_s$, seismic waves at a given depth, the Poisson’s ratio can be calculated from them (Fossen, 2016).

# (ch09-1-1)=
# ### Horizontal stress in the crust
# 
# Elasticity provides a first estimate of the horizontal stress, $\sigma_h$, in the crust. Consider a cube of rock in the subsurface. The cube is subject to an overburden vertical stress, $\sigma_v$, given by Eq. {eq}`ch07_eq02`. Due to this stress, the cube will shorten vertically and tend to expand laterally ({numref}`Figure %s <ch09_fig02>`a), but because it is confined by the rocks around, the horizontal strains will be zero. Since strain is only along the vertical, this condition is known as uniaxial strain.
# 
# Let’s assume $\sigma_{11}=\sigma_v$ and $\sigma_{22}=\sigma_{33}=\sigma_h$. If the horizontal strains are zero, then we have $\varepsilon_h=\varepsilon_{22}=\varepsilon_{33}=0$. From Eq. {eq}`ch09_eq05` it follows that:

# $$
# \varepsilon_{22}=0=\frac{1}{E}[\sigma_h-\nu(\sigma_h+\sigma_v)]
# $$ (ch09_eq07)

# Since $1/E\neq 0$, then:

# $$
# \sigma_h-\nu\sigma_h-\nu\sigma_v=0
# $$ (ch09_eq08)

# or:

# $$
# \sigma_h=\frac{\nu}{1-\nu}\sigma_v
# $$ (ch09_eq09)

# Thus, according to elasticity, the ratio of the horizontal stress to the vertical stress, $\sigma_h/\sigma_v$, is a function of the Poisson’s ratio. This stress ratio is also known as $k$ in rock mechanics (Hoek, 2006). In incompressible materials, $\nu = 0.5$ and $k$ reaches its maximum value of 1.0. In rocks, $\nu\sim 0.25$ and $k\sim 0.33$. Eq. {eq}`ch09_eq09` may be applicable to young sedimentary rocks in large stable basins, or young volcanic areas with lava flows close to the surface. However, over geological time there will be tectonic stresses in active areas, which make the uniaxial assumption invalid (Gudmunsson, 2011). In fact, in the crust at a given depth instead of a single value of horizontal stress $\sigma_h$, we have a maximum, $\sigma_H$, and minimum, $\sigma_h$, horizontal stress (Ragan, 2009; Zoback, 2010). In crustal areas under tectonic compression, $\sigma_H>\sigma_v$ and $k>1.0$ (Hoek, 2006).

# (ch09-2)=
# ## Applications
# 
# In this section, we will discuss four applications of the theory of elasticity. These applications are varied and encompass the fields of geomechanics, basin analysis, structural geology, and geophysics. As in the rest of the book, we introduce briefly the theory and move directly to the code implementation. More details about the applications and the derivation of the equations can be found in the references provided.

# (ch09-2-1)=
# ### Stresses around a circular hole
# 
# The problem of the stresses around a circular hole in a material is one of the most important problems in rock mechanics (Jaeger et al., 2007, their section 8.5). Since most holes drilled through rock are of circular cross section (e.g. tunnels and boreholes), it is not hard to see why this is the case. This problem was first solved by the German engineer Kirsch in 1898. {numref}`Figure %s <ch09_fig03>`a shows the principal stress trajectories around a cylindrical hole in a biaxial stress field based on the Kirsch equations. Since the hole is a free surface, the principal stress trajectories are parallel and perpendicular to it. Where the trajectories of $\sigma_1$ converge, stresses are more compressive (at the azimuth of $\sigma_3$). Where the trajectories of $\sigma_1$ diverge, the stresses are less compressive (at the azimuth of $\sigma_1$).

# ```{figure} /figures/ch09_fig03.png
# :width: 600px
# :name: ch09_fig03
# 
# **a.** Principal stress trajectories around a cylindrical hole under far-field stresses. Continuous lines are $\sigma_{1}$ and dashed lines are $\sigma_{3}$ trajectories. **b.** Plate with a circular hole under far-field stresses. To solve this problem we use a polar coordinate system as indicated in the inset. Modified from Zoback (2010) and Allmendinger (2020).
# ```

# To solve this problem, we use a polar coordinate system as indicated in {numref}`Figure %s <ch09_fig03>`b. The hole has a radius, $a$, and the stresses are calculated at a distance, $r$, from the center of the hole. Any orientation is defined by an angle, $\theta$, measured counterclockwise from $\sigma_3$. For the case where $\sigma_{11} = \sigma_1$ and $\sigma_{33} = \sigma_3$ and there is fluid pressure, $P_f$, in the hole, the radial and tangential normal stresses are (Jaeger et al., 2007):

# $$\begin{gathered}
#     \sigma_{r r}=\frac{\left(\sigma_{1}+\sigma_{3}\right)}{2}\left[1-\left(\frac{a}{r}\right)^{2}\right]+\frac{\left(\sigma_{3}-\sigma_{1}\right)}{2}\left[1-4\left(\frac{a}{r}\right)^{2}+3\left(\frac{a}{r}\right)^{4}\right] \cos 2 \theta+P_{f}\left(\frac{a}{r}\right)^{2} \\ \\
#     \sigma_{\theta \theta}=\frac{\left(\sigma_{1}+\sigma_{3}\right)}{2}\left[1+\left(\frac{a}{r}\right)^{2}\right]-\frac{\left(\sigma_{3}-\sigma_{1}\right)}{2}\left[1+3\left(\frac{a}{r}\right)^{4}\right] \cos 2 \theta-P_{f}\left(\frac{a}{r}\right)^{2}
# \end{gathered}
# $$ (ch09_eq10)

# $\sigma_{\theta\theta}$ is commonly referred to as the *hoop stress*.
# 
# The function [hoop](https://github.com/nfcd/compGeo/blob/master/source/functions/hoop.py) computes the hoop and radial stresses near a circular hole, assuming the hole is on a principal plane, $\sigma_3$ is at $\theta=0\,\textrm{or}\, 180$, and $\sigma_1$ is at $\theta=90\,\textrm{or}\, 270$:

# In[1]:


import numpy as np
import matplotlib.pyplot as plt


def hoop(geom,stress):
    """
    hoop computes the hoop and radial stresses 
    around a circular hole. assuming that 
    the circle is on a principal plane,
    smax (s1) is N-S (theta = 90 or 270), 
    and smin (s3) is E-W (theta = 0 or 180)
    Based on Jaeger et al. (2007)

    USE: shm, srm, fig, ax = hoop(geom,stress)

    geom: A 1 x 2 vector with the number of points 
        along the radius, and the number of points 
        around the circle. These values are used
        to construct the grid around the circle
    stress: A 1 x 3 vector with the value of s1, s3,
        and fluid pressure (pf), all in MPa
    shm, srm: maximum hoop and radial stresses and
        their theta orientations

    fig and ax are handles to the figure and axes
    """
    pi = np.pi # pi 

    # Geometry
    M = geom[0] # points along radius
    N = geom[1] # points around circle
    R1 = 1.0 # radius of hole = 1.0 
    R2 = R1 * 5 # outer radius = 5.0
    nR = np.linspace(R1,R2,M)
    nT = np.linspace(0,2*pi,N)
    R, T = np.meshgrid(nR,nT)
    
    # Convert grid to cartesian coordintes
    X = R*np.cos(T)
    Y = R*np.sin(T)
    m,n = X.shape

    # Principal stresses and pore pressure (MPa)
    s1 = stress[0] 
    s3 = stress[1]
    pf = stress[2]

    # Initialize hoop and radial stresses
    sh = np.zeros(X.shape)
    sr = np.zeros(X.shape)

    # Initialize maximum hoop and radial stresses
    shm = np.zeros(2)
    srm = np.zeros(2)

    # Compute hoop and radial stresses
    for i in range(m):
        for j in range(n):
            
            # Hoop stress
            sh[i,j] = (s1+s3)/2*(1+(R1/R[i,j])**2) \
                -(s3-s1)/2*(1+3*(R1/R[i,j])**4) \
                *np.cos(2*T[i,j]) - pf*(R1/R[i,j])**2
            
            # Radial stress
            sr[i,j] = (s1+s3)/2*(1-(R1/R[i,j])**2) \
                +(s3-s1)/2*(1-4*(R1/R[i,j])**2+3*(R1/R[i,j])**4) \
                *np.cos(2*T[i,j]) + pf*(R1/R[i,j])**2
            
            # maximum hoop stress
            if sh[i,j] > shm[0]:
                shm[0] = sh[i,j]
                shm[1] = T[i,j]*180/pi
            
            # maximum radial stress
            if sr[i,j] > srm[0]:
                srm[0] = sr[i,j]
                srm[1] = T[i,j]*180/pi

    # Create figure
    fig, ax = plt.subplots(2,2,figsize=(10,10))

    # Plot hoop stress
    ax[0,0].axis("equal")
    ax[0,0].set_frame_on(False)
    ax[0,0].get_xaxis().set_visible(False)
    ax[0,0].get_yaxis().set_visible(False)  
    ax[0,0].plot(X[:,0],Y[:,0],"k",linewidth=1.5)
    ax[0,0].plot(X[:,n-1],Y[:,n-1],"k",linewidth=1.5)
    cbar = ax[0,0].contourf(X, Y, sh, cmap="jet")    
    cstr = ax[0,0].contour(X,Y,sh, colors = "k")
    fig.colorbar(cbar, ax=ax[0,0], label="MPa")
    ax[0,0].set_title("Hoop stress",fontweight="bold")

    # Plot radial stress
    ax[0,1].axis("equal")
    ax[0,1].set_frame_on(False)
    ax[0,1].get_xaxis().set_visible(False)
    ax[0,1].get_yaxis().set_visible(False)  
    ax[0,1].plot(X[:,0],Y[:,0],"k",linewidth=1.5)
    ax[0,1].plot(X[:,n-1],Y[:,n-1],"k",linewidth=1.5)
    cbar = ax[0,1].contourf(X, Y, sr, cmap="jet")
    cstr = ax[0,1].contour(X,Y,sr, colors = "k")
    fig.colorbar(cbar, ax=ax[0,1], label="MPa")
    ax[0,1].set_title("Radial stress",fontweight="bold")

    # Plot variation of hoop and radial stress along s3
    ax[1,0].plot(R[0,:]/R1,sh[0,:],"r.-", label="Hoop")
    ax[1,0].plot(R[0,:]/R1,sr[0,:],"b.-", label="Radial")
    ax[1,0].grid(visible=True)
    ax[1,0].set_xlabel("Normalized radial distance")
    ax[1,0].set_ylabel("Stress (MPa)")
    ax[1,0].legend(loc="upper right")
    ax[1,0].set_title("Stress variation along s3",
        fontweight="bold")

    # Plot variation of hoop and radial stress around circle
    ax[1,1].plot(T[:,0]*180/pi,sh[:,0],"r.-", label="Hoop")
    ax[1,1].plot(T[:,0]*180/pi,sr[:,0],"b.-", label="Radial")
    ax[1,1].grid(visible=True)
    ax[1,1].set_xlim([0, 360])
    ax[1,1].set_xticks([0, 90, 180, 270, 360])
    ax[1,1].set_xlabel("Angle around the hole (deg)")
    ax[1,1].set_ylabel("Stress (MPa)")
    ax[1,1].legend(loc="upper right")
    ax[1,1].set_title("Stress variation around circle",
        fontweight="bold")

    return shm, srm, fig, ax


# Let’s use this function to compute the hoop and radial stresses near a circular hole, for $\sigma_1 = 50$ and $\sigma_3 = 25$ MPa. The notebook [ch9-1](https://github.com/nfcd/compGeo/blob/master/source/notebooks/ch9-1.ipynb) shows the solution to this problem:

# In[2]:


# Import function hoop
import sys, os

from compgeo.hoop import hoop

# Define the geometry, 50 points along radius
# 100 points around circle
geom = [50, 100]

# Define stress: sigma1 = 50, sigma3 = 25, Pf = 0 MPa 
stress = [50, 25, 0]

# Compute and plot hoop and radial stresses
shm, srm, fig, ax = hoop(geom, stress)


# Around the hole, the hoop stress is maximum ($\sigma_{\theta\theta}=2.5\sigma_1$) at the azimuth of $\sigma_3$ ($\theta=0\,\textrm{or}\, 180$), and it is minimum ($\sigma_{\theta\theta}=\sigma_3$) at the azimuth of $\sigma_1$ ($\theta=90\,\textrm{or}\, 270$) (right graph). However, the hoop stress decreases quite rapidly to the far-field stress values. At $\theta=0\,\textrm{or}\, 180$, the hoop stress reaches the far-field value of $\sigma_1$ at a distance $5r$ (left graph). The radial stress, $\sigma_{rr}$, is zero around the hole, and reaches the far-field stress values at a distance $5r$. Try experimenting with different values of principal stresses and pore pressure.
# 
# Since in a circular hole, the hoop stress is maximum at the azimuth of $\sigma_3$, and minimum at the azimuth of $\sigma_1$, one would expect the rock to fail by compression near $\sigma_3$, and by tension near $\sigma_1$. This is exactly what happens in a borehole ({numref}`Figure %s <ch09_fig04>`a). In the regions of maximum hoop stress, the rock fails by compression and shear, forming shattered rock areas called *breakouts*. Near the regions of minimum hoop stress, the rock fails by tension, forming tension cracks. This allows determining from borehole images, the orientation of $\sigma_3$, and if the borehole is vertical, the orientation of $\sigma_H$ and $\sigma_h$ ({numref}`Figure %s <ch09_fig04>`b).

# ```{figure} /figures/ch09_fig04.png
# :width: 600px
# :name: ch09_fig04
# 
# **a.** Breakouts and tension cracks around a vertical borehole. **b.** Breakouts in ultrasonic borehole image appear as dark bands. Inset shows a cross-sectional view of the well made with televiewer data. Breakouts indicate a SE-NW $\sigma_{3}$ orientation. Modified from Allmendinger (2020) and Zoback (2010).
# ```

# (ch09-2-2)=
# ### Lithospheric flexure
# 
# The deflection of the surface of the Earth under crustal loads (e.g. a fold-thrust belt; {numref}`Figure %s <ch09_fig05>`a) can be reasonably estimated by an elastic or flexural model. In this model, the uppermost layer of the Earth (i.e. the elastic lithosphere) responds to crustal loads as an elastic beam floating in a weaker, fluid-like foundation, the astenospheric mantle ({numref}`Figure %s <ch09_fig05>`b) (Turcotte and Schubert, 1982; Watts, 2001).

# ```{figure} /figures/ch09_fig05.png
# :width: 600px
# :name: ch09_fig05
# 
# **a.** Fold-thrust belt and foreland basin. **b.** Conceptual elastic model for **a**. **c.** Uniformly distributed load over an infinite beam. The figures show the situation for a point $C$ under, to the left, and to the right of the load. **a** and **b** are from Van Der Pluijm and Marshak (2004) and **c** from Hetenyi (1946).
# ```

# This flexural model provides insight into how tectonic loads (e.g. fold-thrust belts) and sedimentary basins (e.g. foreland basins) are linked, and how the crust and mantle support loads. The flexural model has allowed geologists to understand the regional variations of strength of the lithosphere, and the implications these variations have for mountain building, sedimentary basin formation, and earthquakes (Watts, 2001).
# 
# The physics of this model is the same as that of an elastic beam on an elastic foundation. If the thickness of the beam is constant, the beam is infinite, and the height of the load is constant, the solution is relatively simple (Hetenyi, 1946). We first use the flexural ridigidity, $D$, to express the strength of the elastic lithosphere:

# $$
# D=\frac{ETe^3}{12(1-\nu^2)}
# $$ (ch09_eq11)

# where $Te$ is the thickness of the elastic lithosphere ({numref}`Figure %s <ch09_fig05>`b). We then invent a term $D_{(\alpha,x)}$:

# $$
# \label{eq9.12}
# D_{(\alpha,x)}=\exp(-x/\alpha)\cos(x/\alpha)
# $$ (ch09_eq12)

# where $x$ is the horizontal coordinate and $\alpha$ is a length parameter given by:

# $$
# \alpha=\left[\frac{4D}{\rho_f g}\right]^{\frac{1}{4}}
# $$ (ch09_eq13)

# where $\rho_f$ is the density of the foundation, which is the difference between the density of the mantle, $\rho_m$, and the density of the material filling the basin, $\rho_s$. $g$ is gravity. The deflection $u$ of any point $C$ along the beam can then be computed as ({numref}`Figure %s <ch09_fig05>`c):

# - If the point is under the load:
# 
#     $$
#     u = \frac{q}{2k}(2-D_{(\alpha,a)}-D_{(\alpha,b)})
#     $$ (ch09_eq14)
# 
# - If the point is to the left of the load:
#     $$
#     u = \frac{q}{2k}(D_{(\alpha,a)}-D_{(\alpha,b)})
#     $$ (ch09_eq15)
# 
# - If the point is to the right of the load:
#     $$
#     u = -\frac{q}{2k}(D_{(\alpha,a)}-D_{(\alpha,b)})
#     $$ (ch09_eq16)

# where $\frac{q}{k}=h\frac{\rho}{\rho_f}$, $h$ is the height of the load, and $\rho$ is the load density. $a$ and $b$ are measured as absolute distances from the point to the left and right borders of the load ({numref}`Figure %s <ch09_fig05>`c).
# 
# Although these equations allow us to compute only the deflection produced by a rectangular load column, we can approximate any load profile by discretizing it into narrower load columns. Since the deformation is elastic, the total deflection profile is the sum of the deflection profiles of the load columns (principle of superposition).
# 
# The function [flex2d](https://github.com/nfcd/compGeo/blob/master/source/functions/flex2d.py) computes the deflection profile produced by a group of load columns on an elastic lithosphere resting on a fluid-like foundation:

# In[3]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def flex2d(geom,elas,loads):
    """
    flex2d computes the deflection profile produced by 
    a group of load columns on an elastic lithosphere 
    resting on a fluid-like foundation 
    (i.e. astenosphere). The program is based on 
    Hetenyi (1946) solution for an infinite elastic 
    beam on a fluid-like foundation

    USE: w, wp, fig, ax = flex2d(geom,elas,loads)

    geom: A 1 x 2 vector with the lateral extent of the 
        domain in meters, and the distance between
        points (x interval) in meters where the
        deflection will be computed
    elas: A 1 x 4 vector with the Young Modulus (in Pa), 
        Poisson ratio, Elastic thickness in meters, 
        and density of the foundation in kg/m^3
    loads: A nloads x 4 vector with the left x coordinate,
        right x coordinate, height, and density of each
        load column. Lengths should be in meters and 
        density in kg/m^3
    w: The deflection of the lithosphere in meters at the 
        points specified by the extent and x interval
    wp: 3 x 2 matrix with key deflection parameters: 
        1st row is the maximum deflection (maximum basin
        depth) and x at this location
        2nd row is the zero deflection and x at this 
        location (basin width)
        3rd row is the minimum deflection (forebulge) 
        and x at this location

    fig and ax are handles to the figure and axes
    """
    # Geometry
    extent = geom[0] # Extent in x
    xint = geom[1] # Interval in x
    x = np.arange(0,extent+1,xint) # points in x
    w = np.zeros(len(x)) # initialize displacement

    # Elastic and flexural parameters
    E = elas[0] # Young Modulus
    v = elas[1] # Poisson ratio
    h = elas[2] # Elastic thickness
    
    # Flexural rigidity
    rigid = (E*h*h*h)/(12*(1.0-v*v)) 
    densup = elas[3] # Density of foundation
    g = 9.81 # Gravity
    k = densup*g # Support of foundation
    
    # Flexural parameter
    alpha = ((4*rigid)/(k))**(1/4) 

    # Loads
    lxmin = loads[:,0]
    lxmax = loads[:,1]
    lh = loads[:,2]
    ldens = loads[:,3]

    # Compute deflection profile
    # for all the loads columns
    for i in range(len(lxmin)):
        q = lh[i] * ldens[i] * 9.81
        
        # for all points in x
        for j in range(len(x)):
            tolf = abs(x[j] - lxmin[i])
            tort = abs(x[j] - lxmax[i])
            dA = np.exp(-tolf/alpha)*np.cos(tolf/alpha)
            dB = np.exp(-tort/alpha)*np.cos(tort/alpha)
            
            # If below the load
            if x[j] >= lxmin[i] and x[j] <= lxmax[i]:
                w[j] = w[j] + (q/(2.0*k))*(2.0-dA-dB)
            
            # If to the left of the load
            elif x[j] < lxmin[i]:
                w[j] = w[j] + (q/(2.0*k))*(dA-dB)
            
            # If to the right of the load
            elif x[j] > lxmax[i]:
                w[j] = w[j] - (q/(2.0*k))*(dA-dB)

    # Key deflection parameters
    wp = np.zeros((3,2))
    
    # maximum basin depth
    v = max(w)
    ind = w.argmax()
    wp[0,0] = v
    wp[0,1] = x[ind]
    
    # basin width
    ind = np.where(w <=0)[0]
    wp[1,0] = 0.0
    wp[1,1] = x[ind[0]]
    
    # forebulge
    v = min(w)
    ind = w.argmin()
    wp[2,0] = v
    wp[2,1] = x[ind]

    # Plot
    # Input loads profile
    fig, ax = plt.subplots(2, 1, figsize=(10,8))    

    for i in range(len(lxmin)):
        xcol = [lxmin[i], lxmin[i], lxmax[i], lxmax[i]]
        ycol = [0, lh[i], lh[i], 0]
        ax[0].plot(xcol,ycol,"k-")
    
    maxlh = max(lh)*1.25
    ax[0].axis([0, geom[0], -maxlh, maxlh])
    ax[0].grid()
    ax[0].set_xlabel("m")
    ax[0].set_ylabel("m")
    ax[0].xaxis.set_major_formatter(ticker.FormatStrFormatter
        ("%0.0e"))
    ax[0].set_title("Input loads",fontweight="bold")

    # Deflected loads plus deflection
    wl = np.interp(lxmin,x,w)
    wr = np.interp(lxmax,x,w)
    
    for i in range(len(lxmin)):
        xcol = [lxmin[i], lxmin[i], lxmax[i], lxmax[i]]
        ycol = [-wl[i], lh[i]-wl[i], lh[i]-wr[i], -wr[i]]
        ax[1].plot(xcol,ycol,"k-")
    
    ax[1].plot(x,-w,"b-")
    ax[1].axis([0, geom[0], -maxlh, maxlh])
    ax[1].grid()
    ax[1].set_xlabel("m")
    ax[1].set_ylabel("m")
    ax[1].xaxis.set_major_formatter(ticker.FormatStrFormatter
        ("%0.0e"))
    ax[1].set_title("Deflected loads",fontweight="bold")
    fig.subplots_adjust(hspace=0.6)

    return w, wp, fig, ax


# Let’s use this function to compute the deflection produced by a group of crustal loads. The file [loads.txt](https://github.com/nfcd/compGeo/blob/master/source/data/ch9-2/loads.txt) contains the load columns. Each column is specified by the left and right $x$ coordinates, height and density. The notebook [ch9-2](https://github.com/nfcd/compGeo/blob/master/source/notebooks/ch9-2.ipynb) shows the solution to this problem:

# In[4]:


# Import libraries
import numpy as np

# Import function flex2d
from compgeo.flex2d import flex2d

# Geometry
geom = np.zeros(2)
geom[0] = 500e3 # extent = 500 km
geom[1] = 5e3 # Interval in x =  5 km

# Elastic and flexural parameters
elas = np.zeros(4)
elas[0] = 70e9 # Young Modulus = 70 GPa
elas[1] = 0.25 # Poisson"s ratio = 0.25
elas[2] = 30e3 # Elastic thickness = 30 km
elas[3] = 3300 # Density of mantle in kg/m^3

# loads
loads=np.loadtxt(os.path.abspath("data/ch9-2/loads.txt"))

# Compute deflection profile
w, wp, fig, ax = flex2d(geom, elas, loads)


# (ch09-2-3)=
# ### Faults as elastic dislocations
# 
# One of the simplest ways to model a fault is as an elastic dislocation. A dislocation element consists of a line segment in 2D or a polygon (typically a triangle or rectangle) in 3D across which displacement is discontinuous. In the case of a fault, the value of the dislocation is the amount of slip on the fault. A dislocation element representing a fault is typically embedded in an elastic half space, which is a medium that deforms elastically and has infinite extent downward but a horizontal upper surface (the ground surface) along which no shear stress is allowed. By specifying the position and size of the dislocation element and the amount of slip on it, the displacements, strains, and stresses anywhere in the half space can be calculated.
# 
# {numref}`Figure %s <ch09_fig06>`a-b is a great example from a normal fault earthquake, the M = 7.0, October 28, 1983, Bora Peak, Idaho, earthquake on the Lost River fault (Stein and Barrientos, 1985). The dots in {numref}`Figure %s <ch09_fig06>`a are elevation changes measured at ground locations in a 50 years period. These can be considered coseismic displacements (Stein and Barrientos, 1985). The blue line is the elastic dislocation model for the elevation changes. As you can see, the model fit is quite good. {numref}`Figure %s <ch09_fig06>`b shows the modelled fault in cross section. The fault dips $\sim$<!-- -->45, extends from the surface to a depth of $\sim$<!-- -->14 km, and is consistent with aftershocks (circles) and the earthquake focal mechanism. {numref}`Figure %s <ch09_fig06>`c shows a geologic cross section along the same profile. Notice that the top basement profile is similar to the displacement profile ({numref}`Figure %s <ch09_fig06>`a) but enlarged $\sim$<!-- -->1000 times, as if it were the result of many similar earthquakes over geologic time!

# ```{figure} /figures/ch09_fig06.png
# :width: 600px
# :name: ch09_fig06
# 
# **a.** Elevation changes (dots) produced by the Bora Peak earthquake, and elastic dislocation model (blue). **b.** Cross section including the modelled fault (red), aftershocks (circles), and earthquake focal mechanism. **c.** Geologic cross section. **d.** Geometry of edge dislocation at depth in an elastic half space. Symbols are explained in the text. **a**-**c** are from Stein and Barrientos (1985) and d from Segall (2010).
# ```

# To model faults as elastic dislocation elements in 2D, we use the geometry illustrated in {numref}`Figure %s <ch09_fig06>`d. The dislocation is located at $\xi_1 = 0$, $\xi_2$. The displacements along $x_1$ and $x_2$ are (Segall, 2010):

# $$
# %\resizebox{1.05\hsize}{!}
# {\begin{aligned} u_{1}=& -\frac{s_{1}}{\pi(1-\nu)}\left\{\frac{(1-\nu)}{2}\left(\theta_{2}-\theta_{1}\right)+\frac{x_{1}\left(x_{2}-\xi_{2}\right)}{4 r_{1}^{2}}-\frac{x_{1}\left[x_{2}+(3-4 \nu) \xi_{2}\right]}{4 r_{2}^{2}}+\frac{\xi_{2} x_{2} x_{1}\left(x_{2}+\xi_{2}\right)}{r_{2}^{4}}\right\} \\ &+\frac{s_{2}}{\pi(1-\nu)}\left\{\frac{(1-2 \nu)}{4} \log \left(r_{2} / r_{1}\right)-\frac{\left(x_{2}-\xi_{2}\right)^{2}}{4 r_{1}^{2}}+\frac{\left[x_{2}^{2}+\xi_{2}^{2}-4(1-\nu) \xi_{2}\left(x_{2}+\xi_{2}\right)\right]}{4 r_{2}^{2}}\right.\\ &\left.+\frac{x_{2} \xi_{2}\left(x_{2}+\xi_{2}\right)^{2}}{r_{2}^{4}}\right\} \\ u_{2}=& -\frac{s_{1}}{\pi(1-\nu)}\left\{\frac{(1-2 \nu)}{4} \log \left(r_{2} / r_{1}\right)+\frac{\left(x_{2}-\xi_{2}\right)^{2}}{4 r_{1}^{2}}-\frac{\left[\left(x_{2}+\xi_{2}\right)^{2}-2 \xi_{2}^{2}-2(1-2 \nu) \xi_{2}\left(x_{2}+\xi_{2}\right)\right]}{4 r_{2}^{2}}\right.\\ &\left.+\frac{x_{2} \xi_{2}\left(x_{2}+\xi_{2}\right)^{2}}{r_{2}^{4}}\right\}+\frac{s_{2}}{\pi(1-\nu)}\left\{\frac{(1-\nu)}{2}\left(\theta_{1}-\theta_{2}\right)+\frac{x_{1}\left(x_{2}-\xi_{2}\right)}{4 r_{1}^{2}}-\frac{x_{1}\left[x_{2}+(3-4 \nu) \xi_{2}\right]}{4 r_{2}^{2}}\right.\\ &\left.-\frac{\xi_{2} x_{2} x_{1}\left(x_{2}+\xi_{2}\right)}{r_{2}^{4}}\right\} \end{aligned}$}
# $$ (ch09_eq17)

# where $s_1$ and $s_2$ are the slip components along $x_1$ and $x_2$, and:

# $$
# \begin{array}{l}r_{1}^{2}=x_{1}^{2}+\left(x_{2}-\xi_{2}\right)^{2} \\ \\ r_{2}^{2}=x_{1}^{2}+\left(x_{2}+\xi_{2}\right)^{2}\end{array}
# $$ (ch09_eq18)

# measure the squared distance from the observation point to the dislocation and the image dislocation (a mirror image of the dislocation above the ground surface), respectively ({numref}`Figure %s <ch09_fig06>`d). $\theta_1$ refers to the angle about the dislocation, and $\theta_2$ refers to the angle about the image dislocation ({numref}`Figure %s <ch09_fig06>`d):

# $$
# \begin{gathered}
#     \theta_{1}=\tan ^{-1}\left(\frac{x_{1}-\xi_{1}}{x_{2}-\xi_{2}}\right) \\
#     \theta_{2}=\tan ^{-1}\left(\frac{x_{1}-\xi_{1}}{x_{2}+\xi_{2}}\right)
# \end{gathered}
# $$ (ch09_eq19)

# Notice that in terms of elastic moduli, the displacements in Eq. {eq}`ch09_eq17` depend only on the Poisson’s ratio $\nu$. Also for a dislocation not at $\xi_1 = 0$, simply replace $x_1$ with $x_1 - \xi_1$ in Eq. {eq}`ch09_eq17`.
# 
# The function [disloc2d](https://github.com/nfcd/compGeo/blob/master/source/functions/disloc2d.py) computes the displacements on a 2D planar fault of finite extent, modeled by two edge dislocations in an elastic half space. Notice that `disloc2d` uses the function `displacement`, which computes Eqs. {eq}`ch09_eq17`-{eq}`ch09_eq19`.

# In[5]:


from numpy import arctan,arctan2,sin,cos,sign,sqrt,pi,mod,log


def disloc2d(tip,base,slip,nu,obsx,obsy):
    """
    This function calculates displacements on a 2D planar 
    fault of finite extent, modeled by two edge dislocations 
    in a homogeneous, isotropic elastic half space

    Arguments:
        tip = tuple of (x,y) coordinates of the fault tip
        base = tuple of (x,y) coordinates of the fault base
        slip = slip on the fault, + for reverse, - for normal
        nu = Poisson's ratio (scalar)
        obsx = x coordinates of observation points 
        obsy = y coordinates of observation points
    Returns:
        ux = x components of displ. vectors at obs. points 
        uy = y components of displ. vectors at obs. points

    disloc2d uses function displacement

    Written by David Oakley (david.o.oakley@uis.no)
    """
    dip = arctan2(tip[1]-base[1],tip[0]-base[0])
    s1 = slip*cos(dip)
    s2 = slip*sin(dip)
    
    [ux_part1,uy_part1] = displacement(tip[0],tip[1],s1,s2,
        nu,obsx,obsy)
    [ux_part2,uy_part2] = displacement(base[0],base[1],-s1,-s2,
        nu,obsx,obsy)
    
    ux = ux_part1+ux_part2
    uy = uy_part1+uy_part2
    
    return ux,uy


def displacement(xi1,xi2,s1,s2,nu,x1,x2):
    """
    Calculate displacements (u1,u2) in 2D in a half space at
    points (x1,x2) due to an edge dislocation at (xi1,xi2) with 
    slip vector (s1,s2). This uses Eqs from Segall (2010), 
    as corrected in the Errata to that book. Indices 1 and 2 
    correspond to the x and y directions respectively.
    Arguments:
        xi1 = x coordinate of the dislocation tip 
        xi2 = y coordinate of the dislocation tip 
        x1 = x coordinates of observation points 
        x2 = y coordinates of observation points 
        s1 = x component of slip vector 
        s2 = y component of slip vector 
        nu = Poisson"s ratio 
    Returns:
        ux = x components of displ. vectors at obs. points 
        uy = y components of displ. vectors at obs. points 
    """
    # This occurs if the fault dips to the left
    if sign(s1)==sign(s2): 
        # Flip to the sign convention that Segall"s equations use
        s1 = -s1
        s2 = -s2
    
    # These equations are written for xi1 = 0. If that's not 
    # the case, this makes it equivalent to that case
    x1 = x1-xi1 
    r1_sq = x1**2.+(x2-xi2)**2.
    r2_sq = x1**2.+(x2+xi2)**2.
    r1 = sqrt(r1_sq)
    r2 = sqrt(r2_sq)
    log_r2_r1 = log(r2/r1)
    
    # Calculate the angles relative to the vertical axis
    theta1 = arctan2(x1,(x2-xi2))
    theta2 = arctan2(x1,(x2+xi2))
    dip = arctan(s2/s1)
    
    # This puts dip in the range [0,pi]
    if dip<0:
        dip=dip+pi 
    
    # The following puts the atan branch cuts along the fault 
    # by rotating theta1 to point down the fault and theta2 to 
    # point up (above the half space) along it
    # Shift theta1 to be measured up from the fault
    theta1 = theta1+pi/2.+dip
    
    # Shift theta2 to be measured from a line pointing up 
    # opposite the fault from the image dislocation
    theta2 = theta2+dip-pi/2. 
    theta1 = mod(theta1,2.*pi)
    theta2 = mod(theta2,2.*pi)
    
    # Make a correction for rounding errors that can occur 
    # when theta1 or theta2 is very close to 0 or 2*pi.
    theta1[2.*pi-theta1<1e-10] = 0.
    theta2[2.*pi-theta2<1e-10] = 0.
    theta_diff = theta1-theta2
    
    # These are two very long equations....
    u1 = (-s1/(pi*(1.-nu)))*(((1.-nu)/2.)*(-theta_diff)+x1*(x2-xi2)/(4.*r1_sq) - x1*(x2+(3.-4.*nu)*xi2)/(4.*r2_sq)+xi2*x2*x1*(x2+xi2)/r2_sq**2.) + (s2/(pi*(1.-nu)))*(((1.-2.*nu)/4.)*log_r2_r1-(x2-xi2)**2./(4.*r1_sq) + (x2**2.+xi2**2.-4.*(1.-nu)*xi2*(x2+xi2))/(4.*r2_sq)+x2*xi2*(x2+xi2)**2./r2_sq**2.)
    u2 = (-s1/(pi*(1.-nu)))*(((1.-2.*nu)/4.)*log_r2_r1+(x2-xi2)**2./(4.*r1_sq) - ((x2+xi2)**2.-2.*xi2**2.-2.*(1.-2.*nu)*xi2*(x2+xi2))/(4.*r2_sq) + x2*xi2*(x2+xi2)**2./r2_sq**2.)+(s2/(pi*(1.-nu)))*(((1.-nu)/2.)*(theta_diff) + x1*(x2-xi2)/(4.*r1_sq)-x1*(x2+(3.-4.*nu)*xi2)/(4.*r2_sq)-xi2*x2*x1*(x2+xi2)/r2_sq**2.)

    return u1, u2


# Let’s use this function to compute the surface deformation, displacement vectors, and fold growth produced by a 30reverse fault. The notebook [ch9-3](https://github.com/nfcd/compGeo/blob/master/source/notebooks/ch9-3.ipynb) shows the solution to this problem. First, we import the required libraries and functions:

# In[6]:


# Import libraries
import numpy as np
from matplotlib import pyplot as plt

# Import function disloc2d
from compgeo.disloc2d import disloc2d


# Then, we define the geometry of the fault (upper and lower tips), Poisson’s ratio, and fault slip. The units of length are meters. In this case, we define a blind (non-surface breaking) 30dipping reverse fault, with its upper tip 500 m below the surface (and at $x = 0$ m) and its lower tip 1500 m below the surface. The fault slip is 2 m:

# In[7]:


# Define the fault:
dip = 30.0*np.pi/180.0 #Fault dip in radians
xt = 0.0
yt = -500.0
yb = -1500.0
xb = xt-(yt-yb)/np.tan(dip)
tip = (xt,yt)
base = (xb,yb)
nu = 0.25 #Poisson's ratio
slip = 2.0 #Positive for reverse, negative for normal


# To compute the deformation at the Earth’s surface produced by this fault, we create observation points at the surface on a transect perpendicular to the fault. We apply the elastic deformation, and plot the uplift of the points. We observe a maximum uplift of about 0.8 m above the fault as well as a small amount of subsidence away from the fault. Deformation of the Earth’s surface such as modeled here is often observed after an earthquake and can be used to make inferences about the geometry of and amount of slip on the fault that caused the earthquake.

# In[8]:


#Define the observation points at the surface
obsx1 = np.arange(-3000.0,3100.0,10.0,dtype=float)
obsy1 = np.zeros(obsx1.shape)

#Apply slip and deform the surface
ux,uy = disloc2d(tip,base,slip,nu,obsx1,obsy1)
obsx1 = obsx1+ux
obsy1 = obsy1+uy

#Plot the result
fig, ax = plt.subplots()
ax.plot(obsx1,obsy1,"k-")
ax.set_xlabel("Distance (m)")
ax.set_ylabel("Uplift (m)")
ax.set_title("Surface Deformation")
plt.show()


# We now create a two-dimensional grid of observation points around the fault and display the displacement vectors at these points. From the resulting figure, we can clearly see that material moves upward in the hanging wall of the fault and downward in its footwall with displacements greatest near the fault and continuous everywhere except across the fault.

# In[9]:


#Define observation points in a region around the fault.
obsx2, obsy2 = np.meshgrid(np.arange(-2500.0,1100.0,100.0,
    dtype=float),np.arange(-2500.0,100.0,100.0,dtype=float))

#Ignore divide by zero error for point right at fault tip.
np.seterr(divide='ignore', invalid='ignore') 
#Apply slip and calculate displacement vectors.
ux,uy = disloc2d(tip,base,slip,nu,obsx2,obsy2)

#Plot the result.
fig, ax = plt.subplots()
ax.plot([xt,xb],[yt,yb],'r-') #Plot the fault
q = ax.quiver(obsx2, obsy2, ux, uy)
ax.quiverkey(q, X=0.6, Y=1.05, U=1, 
             label='1 m displacement', labelpos='E')
ax.set_xlabel('Distance (m)')
ax.set_ylabel('Elevation (m)')
ax.set_title('Displacement Vectors',loc='left')
ax.axis('equal')
plt.show()


# While the plots above are from a single earthquake, repeated slip on a fault will cause much larger cumulative deformation ({numref}`Figure %s <ch09_fig06>`c). The rocks above blind faults are frequently folded, and elastic dislocation theory provides a way to model the formation of such folds in response to faulting. We create three horizontal rows of observation points, representing stratigraphic horizons, and we then apply 200 increments of slip (of 2 m each) to the fault. The result is that the stratigraphic layers above the fault are folded into an asymmetric anticline.

# In[10]:


#Define the observation points as horizontal beds
obsx3, obsy3 = np.meshgrid(np.arange(-2500.0,1100.0,
    100.0,dtype=float),[-400.0,-300.0,-200.0])

#Plot the initial state:
fig, ax = plt.subplots(2,1)
fig.subplots_adjust(hspace=0.5)
ax[0].plot([xt,xb],[yt,yb],'r-') #Plot the fault
for i in range(0,obsx3.shape[0]):
    ax[0].plot(obsx3[i,:],obsy3[i,:],'k-')
ax[0].set_ylabel('Elevation (m)')
ax[0].set_title('Initial State')
ax[0].axis('equal')

#Apply repeated slip increments and displace 
#the observation points.
n_inc = 200 #Total number of increments.
for i in range(n_inc):
    ux,uy = disloc2d(tip,base,slip,nu,obsx3,obsy3)
    obsx3 = obsx3+ux
    obsy3 = obsy3+uy
    
#Plot the final state:
ax[1].plot([xt,xb],[yt,yb],'r-') #Plot the fault
for i in range(obsx3.shape[0]):
    ax[1].plot(obsx3[i,:],obsy3[i,:],'k-')
ax[1].set_xlabel('Distance (m)')
ax[1].set_ylabel('Elevation (m)')
ax[1].set_title('Final State')
ax[1].axis('equal')
plt.show()


# Try a normal fault dipping 60.

# (ch09-2-4)=
# ### Wave propagation
# 
# The propagation of waves in elastic media can be derived from the equation of motion and Hooke’s law. The equation of motion relates the displacements, $\mathbf{u}$, to the stresses, $\boldsymbol{\sigma}$. Combining Hooke’s law with the equation of motion leads to a system of equations describing the propagation of stresses and displacements in a continuum. The equation of motion is given by (Pollard and Fletcher, 2005):

# $$
# \rho\frac{\partial^2 u_i}{\partial t^2} = \frac{\partial \sigma_{ij}}{\partial X_j}
# $$ (ch09_eq20)

# where $\rho$ is density. The Hooke’s law (Eq. {eq}`ch09_eq01`) can also be written as (Pollard and Fletcher, 2005):

# $$
# \sigma_{ij}= 2G \varepsilon_{ij} + \lambda \varepsilon_{kk}\delta_{ij}
# $$ (ch09_eq21)

# where $\lambda = \frac{E\nu}{(1 + \nu)(1-2\nu)}$ is the Lamé’s constant, $G$ is the shear modulus (Eq. {eq}`ch09_eq04`), and $\varepsilon_{ij} = \frac{1}{2}\left(\frac{\partial u_i}{\partial X_j}+\frac{\partial u_j}{\partial X_i}\right)$ is the infinitesimal strain tensor (Eq. {eq}`ch08_eq15`).
# 
# Let’s assume a 2-D isotropic and linearly elastic medium with a horizontal axis, $\mathbf{x}$, and a vertical axis, $\mathbf{z}$.
# 
# :::{note}
# For simplicity, we use in this section $\mathbf{x}$ and $\mathbf{z}$ instead of $\mathbf{X_1}$ and $\mathbf{X_2}$.
# :::
# 
# Eqs. {eq}`ch09_eq20` and {eq}`ch09_eq21` become:

# $$
# \rho\frac{\partial^2u_x}{\partial t^2}=\frac{\partial\sigma_{xx}}{\partial x}+\frac{\partial\sigma_{xz}}{\partial z}
# $$ (ch09_eq22)
# 
# $$
# \rho\frac{\partial^2u_z}{\partial t^2}=\frac{\partial\sigma_{xz}}{\partial x}+\frac{\partial\sigma_{zz}}{\partial z}
# $$ (ch09_eq23)
# 
# $$
# \sigma_{xx}=\left(\lambda+2G\right)\frac{\partial u_x}{\partial x}+\lambda\frac{\partial u_z}{\partial z}
# $$ (ch09_eq24)
# 
# $$
# \sigma_{zz}=\left(\lambda+2G\right)\frac{\partial u_z}{\partial z}+\lambda\frac{\partial u_x}{\partial x}
# $$ (ch09_eq25)

# and

# $$
# \sigma_{xz}=G\left(\frac{\partial u_x}{\partial z}+\frac{\partial u_z}{\partial x}\right)
# $$ (ch09_eq26)

# These equations are known as the *elastodynamic* equations. They can be further simplified into a first-order hyperbolic system of equations by using velocity instead of displacement (Virieux, 1986):

# $$
# \frac{\partial v_x}{\partial t}=b\left(\frac{\partial \sigma_{xx}}{\partial x}+\frac{\partial\sigma_{xz}}{\partial z}\right)
# $$ (ch09_eq27)
# 
# $$
# \frac{\partial v_z}{\partial t}=b\left(\frac{\partial \sigma_{xz}}{\partial x}+\frac{\partial\sigma_{zz}}{\partial z}\right)
# $$ (ch09_eq28)
# 
# $$
# \frac{\partial\sigma_{xx}}{\partial t}=\left(\lambda+2G\right)\frac{\partial v_x}{\partial x}+\lambda\frac{\partial v_z}{\partial z}
# $$ (ch09_eq29)
# 
# $$
# \frac{\partial\sigma_{zz}}{\partial t}=\left(\lambda+2G\right)\frac{\partial v_z}{\partial z}+\lambda\frac{\partial v_x}{\partial x}
# $$ (ch09_eq30)

# and

# $$
# \frac{\partial\sigma_{xz}}{\partial t}=G\left(\frac{\partial v_x}{\partial z}+\frac{\partial v_z}{\partial x}\right)
# $$ (ch09_eq31)

# where $v_x, \; v_z$ are the particle velocities, and $b$ is the buoyancy, the inverse of the density. Virieux (1986) solved Eqs. {eq}`ch09_eq27` to {eq}`ch09_eq31` using a finite-difference method and staggered grids. He also gave the numerical stability conditions of his scheme for different media.
# 
# :::{note}
# The finite difference method is a numerical method for solving differential equations as differences between points in a grid of regular cells. Staggered grids are superimposed grids that are slightly shifted in location. In our case, the stress and velocity grids are slightly shifted.
# :::
# 
# The program [CGeo$\_$elastic](https://github.com/nfcd/compGeo/blob/master/source/functions/CGeo_elastic.py) simulates the propagation of waves in elastic media. This code is implemented using Python classes, which are objects that can have both functions and data elements. There are four classes in the program: `Source` which implements the wave source, `Derivatives` which implements the finite differences that solve the derivatives, `Elastic`$\_$`model` which implements and plots the elastic model, and `Elastic`$\_$`waves` which simulates the elastic waves by solving Eqs. {eq}`ch09_eq27` to {eq}`ch09_eq31`.
# 
# The notebook [ch9-4](https://github.com/nfcd/compGeo/blob/master/source/notebooks/ch9-4.ipynb) illustrates the use of the `CGeo`$\_$`elastic` program. We first import the necessary libraries and classes:

# In[11]:


import numpy as np
import matplotlib.pyplot as plt
from IPython.display import clear_output

# Import CGeo_elastic classes
import compgeo.CGeo_elastic as ela


# We consider a region consisting of two materials with different seismic velocities but homogeneous densities separated by a 90corner. $V_p$ is 6000 m/s in the overburden and 9000 m/s inside the corner. The density, $\rho$, is 2500 kg/m$^3$. The S-wave velocities are computed using the equations $V_p = \sqrt{\frac{\lambda + 2G}{\rho}}$ and $V_s = \sqrt{\frac{G}{\rho}}$ (Sheriff and Geldart, 1995), and assuming $\lambda = G$, which corresponds to a Poisson ratio $\nu = 0.25$. The model is discretized using 100 m grid spacing in the $x$ and $z$ directions:

# In[12]:


# Model size in grid cells
nz = 249
nx = 381
dz = 100 
dx = 100

# Corner edge model with Poisson ratio = 0.25
vpdata = np.ones([nz,nx])*6000
vpdata[89:, 131:] = 9000
rhodata = np.ones([nz,nx])*2500
lam = rhodata*(vpdata**2)/3. 
G = lam 
vsdata = np.sqrt(G/rhodata)

# Initialize 2D elastic model class and plot model
Model = ela.Elastic_model(vpdata, vsdata, 
                          rhodata, dx, dz, 0.0, 0.0) 
fig, ax = Model.plot()


# The source wavelet consists of a Ricker wavelet with a dominant frequency of $\sqrt{\frac{40}{\pi^2}} \approx 2$ Hz. The time discretization and the number of time samples are $\Delta t = 0.005$, and $nt = 1201$, respectively. This gives a total modelling time of 6 seconds. The source is modelled as punctual explosion whose location is $x = 15000$ m and $z=5800$ m. Furthermore, seismic data is recorded right below the free surface at the top of the model:

# In[13]:


# Sampling and modeling length (6 seconds = (nt-1)*dt)
dt = 5e-3
nt = 1201

f0=(40/(np.pi**2))**.5  # Dominant frequency
t0=0.6 # Time delay

# Source position in grid points
sx = 150
sz = 58

# Receiver depth in grid points for the recordings
rz = 1

Source = ela.Source(nt, dt, sx, sz) # Source class
Source.Ricker(f0,t0,0) # Initializing a source wavelet
fig, ax = Source.plot()


# We then setup the finite difference scheme:

# In[14]:


# Initialize the classes to solve the Elastodynamic equations
Waves = ela.Elastic_waves(Model,nt,dt) 
Derivative = ela.Derivatives()

# Create containers for snapshots 
# These require lots of memory
P_snaps = np.zeros([nz,nx,nt]) # pressure
Vx_snaps = np.zeros([nz,nx,nt]) # vx
Vz_snaps = np.zeros([nz,nx,nt]) # vz

# Create containers for seismograms
# just below the free surface
P_record= np.zeros([nt,nx]) # pressure
Vx_record = np.zeros([nt,nx]) # vx 
Vz_record = np.zeros([nt,nx]) # vz

# Check stability of FD modeling scheme
dtstab = Waves.Courant_stability(np.max(Model.vp))
if(dt > dtstab):
    raise Exception("The value of dt should not exceed", 
                    "the stability limit of:", 
                    dtstab, "The value of dt was:", dt)


# Finally, the elastodynamic waves are solved using an explicit leap frog method where the wavefields at the next time step are solved from the wavefields at the previous time step. You will be presented with a progress indicator, please wait for the program to finish:

# In[15]:


# Loop over time
for it in range(0,nt):
    # Extrapolate waves one time step
    Waves.forwardStep(Derivative, Model)
    
    # Adding pressure source 
    Waves.insertPressure(Source, it)

    # Adding force source(s)
    # Clockwise angle gives force direction (180 is down)
    #angle = 180 
    #Waves.insertForce(Source, Model, angle, it)    
    
    # Recording a seismogram 
    P_record[it, :] = Waves.recordPressure(rz)
    Vx_record[it, :] = Waves.recordVelocity(rz, 'x')
    Vz_record[it, :] = Waves.recordVelocity(rz, 'z')
    
    # Record snapshots
    Vx_snaps[:,:,it] = Waves.Vx 
    Vz_snaps[:,:,it] = Waves.Vz 
    P_snaps[:,:,it] = 0.5*(Waves.Sxx + Waves.Szz) 
            
    if(it % np.floor(nt/20) == 0):
        clear_output(wait=True)
        print('Progress:', np.round(100*it/(nt-1)), '%')


# The results can be displayed in the form of numerical seismograms recorded just below the free surface. The left column corresponds to the $x$-component of the particle velocities ($v_x$), the middle column corresponds to the $z$-component of the particle velocities ($v_z$), and the right column corresponds to pressure ($\frac{\sigma_{xx} + \sigma_{zz}}{2}$). Both compressional and shear waves can be seen on the left and middle columns, whereas only compressional waves are visible on the right column:  

# In[16]:


fig, ax = plt.subplots(1,3, sharey=True, figsize=(16,24))
extents = [0, (nx-1)*dx, (nt-1)*dt, 0]

# Vx
vlim = 0.1*np.min(Vx_record)
ax[0].imshow(Vx_record,vmin=vlim, vmax=-vlim, 
             extent=extents, aspect=8000, cmap="gray")
ax[0].set_xlabel("Position (m)")
ax[0].set_ylabel("Time (s)")
ax[0].set_title("Vx")

# Vz
vlim = 0.1*np.min(Vz_record)
ax[1].imshow(Vz_record,vmin=vlim, vmax=-vlim, 
             extent=extents, aspect=8000, cmap="gray")
ax[1].set_xlabel("Position (m)")
ax[1].set_ylabel("Time (s)")
ax[1].set_title("Vz")

# Pressure
vlim = 0.1*np.min(P_record)
ax[2].imshow(P_record,vmin=vlim, vmax=-vlim, 
             extent=extents, aspect=8000, cmap="gray")
ax[2].set_xlabel("Position (m)")
ax[2].set_ylabel("Time (s)")
ax[2].set_title("Pressure")
plt.show()


# The results can also be shown as snapshots of the wave propagation at selected times. Again, the left column is $v_x$, the middle column is $v_z$, and the right column is pressure. Each row is a snapshot and time increases downwards. Compressional and shear waves are observed on the left and middle columns, whereas only compressional waves are visible on the right column:

# In[17]:


fig, ax = plt.subplots(4,3, sharex=True, 
                       sharey=True, figsize=(12,12))
extents = [0, (nx-1)*dx, (nz-1)*dz, 0]
snapit = int(np.floor(0.995/dt)) # initial time
delta_snap = int(np.floor(0.5/dt)) # delta time

for i in range(0,4): # 4 snaps
    
    # Vx
    vlim = 0.1*np.min(Vx_snaps[:,:,snapit])
    ax[i,0].imshow(Vx_snaps[:,:,snapit], vmin=vlim, 
                   vmax=-vlim, extent=extents, cmap="gray")
    ax[i,0].grid()
    ax[i,0].set_ylabel("Depth (m)")
    
    # Vz
    vlim = 0.1*np.min(Vz_snaps[:,:,snapit])
    ax[i,1].imshow(Vz_snaps[:,:,snapit], vmin=vlim, 
                   vmax=-vlim, extent=extents, cmap="gray")
    if i == 0:
        ax[i,1].set_title("Vz, " + str((snapit-1)*dt) + " s")
    
    else:
        ax[i,1].set_title(str((snapit-1)*dt) + " s")
    ax[i,1].grid()
    
    # Pressure
    vlim = 0.1*np.min(P_snaps[:,:,snapit])
    ax[i,2].imshow(P_snaps[:,:,snapit], vmin=vlim, 
                   vmax=-vlim, extent=extents, cmap="gray")
    ax[i,2].grid()
    snapit = snapit+delta_snap

ax[0,0].set_title("Vx");
ax[0,2].set_title("Pressure");
ax[3,0].set_xlabel("Position (m)")
ax[3,1].set_xlabel("Position (m)")
ax[3,2].set_xlabel("Position (m)")
plt.show()


# (ch09-3)=
# ## Exercises
# 
# 1. The tangential stress, $\sigma_{\eta\eta}$, at the surface of an elliptical hole under far field principal stresses, $\sigma_1$ and $\sigma_3$, is given by (Jaeger et al., 2007):
# 
#     $$
#     \sigma_{\eta \eta}=\frac{2 a b\left(\sigma_{1}+\sigma_{3}\right)+\left(\sigma_{1}-\sigma_{3}\right)\left[\left(a^{2}-b^{2}\right) \cos 2 \beta-(a+b)^{2} \cos 2(\beta-\eta)\right]}{\left(a^{2}+b^{2}\right)-\left(a^{2}-b^{2}\right) \cos 2 \eta}
#     $$ (ch09_eq32)
# 
#     where $a$ and $b$ are the major and minor axes of the hole, $\beta$  is the angle $\sigma_1$ makes with the major axis of the hole ({numref}`Figure %s <ch09_fig07>`), and $\eta$ is given by:
# 
#     $$
#     \tan\eta=\left(\frac{a}{b}\right)\tan\theta
#     $$ (ch09_eq33)
# 
#     where $\theta$ is the angle a vector connecting the center of the hole and the point of observation makes with the long axis of the hole ({numref}`Figure %s <ch09_fig07>`).

# ```{figure} /figures/ch09_fig07.png
# :width: 400px
# :name: ch09_fig07
# 
# Elliptical hole under far field stresses. Symbols are explained in the text.
# ```

#     Write a function that plots the tangential stress as function of $\theta$, on the surface of an elliptical hole under far field stresses, $\sigma_1$ and $\sigma_3$. Check your function with the three cases in Figure 8.6 of Jaeger et al. (2007).
# 
#     :::{hint}
#     Start with function `hoop`, and modify the `geom` and `stress` entries to fit this new problem. You will need to implement Eq. {eq}`ch09_eq32` for the tangential stress. You can use values of $\theta$ from 0 to 180 and 1increment.
#     :::
# 
# 2. The elastic lithosphere of the Earth can be broken by major structures, such as fault systems bounding terrain boundaries. In this case, a broken plate model is more appropriate than the continuous plate model. The implementation of the broken plate model with a free end at $x=0$, and for a rectangular load is explained in Hetenyi (1946). First, we introduce the following terms:
# 
#     $$
#     \begin{gathered}
#         A_{(\alpha,x)} = \exp(-x/\alpha)(\cos(x/\alpha) + \sin(x/\alpha)) \\ B_{(\alpha,x)} = \exp(-x/\alpha)\sin(x/\alpha) \\ 
#         C_{(\alpha,x)} = \exp(-x/\alpha)(\cos(x/\alpha) - \sin(x/\alpha))
#     \end{gathered}
#     $$ (ch09_eq34)
# 
#     and then these terms:
# 
#     $$
#     \begin{gathered}
#         aa = 2(B_{(\alpha,b')} - B_{(\alpha,a')}) + (C_{(\alpha,a')} - C_{(\alpha,b')}) \\ 
#         bb = (B_{(\alpha,b')} - B_{(\alpha,a')}) + (C_{(\alpha,a')} - C_{(\alpha,b')})
#     \end{gathered}
#     $$ (ch09_eq35)
# 
#     where $a'$ and $b'$ are the distances from the free end of the plate, $x=0$, to the left and right borders of the load, respectively. The deflection $u$ of any point along the beam ($x\geq 0$) is:
# 
#     - If the point is under the load:
#         $$
#         u = \frac{q}{2k}[(2-D_{(\alpha,a)}-D_{(\alpha,b)}) + bbA_{(\alpha,x)} - aaB_{(\alpha,x)}]
#         $$ (ch09_eq36)
# 
#     - If the point is to the left of the load:
#         $$
#         u = \frac{q}{2k}[(D_{(\alpha,a)}-D_{(\alpha,b)}) + bbA_{(\alpha,x)} - aaB_{(\alpha,x)}
#         $$ (ch09_eq37)
# 
#     - If the point is to the right of the load:
#         $$
#         u = \frac{q}{2k}[(D_{(\alpha,b)}-D_{(\alpha,a)}) + bbA_{(\alpha,x)} - aaB_{(\alpha,x)}
#         $$ (ch09_eq38)
# 
#     where $a$ and $b$ are the absolute distances from the point to the left and right borders of the load ({numref}`Figure %s <ch09_fig05>`c), and $D_{(\alpha,x)}$ is defined in Eq. {eq}`ch09_eq12`.
# 
#     Write a function that computes the deflection of a broken elastic lithosphere under a distribution of load columns. Check this function with the loads in [loads.txt](https://github.com/nfcd/compGeo/blob/master/source/data/ch9-2/loads.txt).
# 
#     :::{hint}
#     Modify function `flex2d` for the broken plate model. You will need to implement the new terms (Eqs. {eq}`ch09_eq34` and {eq}`ch09_eq35`), and modify the deflection according to Eqs. {eq}`ch09_eq36` to {eq}`ch09_eq38`.
#     :::
# 
# 3. The stress changes caused by an edge dislocation in an elastic half space in two dimensions are (Segall, 2010):
# 
#     $$
#     %\resizebox{1.05\hsize}{!}
#     {\begin{aligned}
#         \sigma_{11}=\frac{G s_{2}}{2 \pi(1-\nu)} &\left\{\frac{x_{1}\left[\left(x_{2}-\xi_{2}\right)^{2}-x_{1}^{2}\right]}{r_{1}^{4}}-\frac{x_{1}\left[\left(x_{2}+\xi_{2}\right)^{2}-x_{1}^{2}\right]}{r_{2}^{4}}\right.\\
#         +&\left.\frac{4 \xi_{2} x_{1}}{r_{2}^{6}}\left[\left(2 \xi_{2}-x_{2}\right)\left(x_{2}+\xi_{2}\right)^{2}+\left(3 x_{2}+2 \xi_{2}\right) x_{1}^{2}\right]\right\} \\
#         &+\frac{G s_{1}}{2 \pi(1-\nu)}\left\{\frac{\left(x_{2}-\xi_{2}\right)\left[\left(x_{2}-\xi_{2}\right)^{2}+3 x_{1}^{2}\right]}{r_{1}^{4}}-\frac{\left(x_{2}+\xi_{2}\right)\left[\left(x_{2}+\xi_{2}\right)^{2}+3 x_{1}^{2}\right]}{r_{2}^{4}}\right.\\
#         &\left.+\frac{2 \xi_{2}}{r_{2}^{6}}\left[6 x_{2}\left(x_{2}+\xi_{2}\right) x_{1}^{2}-\left(x_{2}-\xi_{2}\right)\left(x_{2}+\xi_{2}\right)^{3}-x_{1}^{4}\right]\right\} \\ \\
#         \sigma_{22}=\frac{-G s_{2}}{2 \pi(1-\nu)} &\left\{\frac{x_{1}\left[3\left(x_{2}-\xi_{2}\right)^{2}+x_{1}^{2}\right]}{r_{1}^{4}}-\frac{x_{1}\left[3\left(x_{2}+\xi_{2}\right)^{2}+x_{1}^{2}\right]}{r_{2}^{4}}-\frac{4 \xi_{2} x_{2} x_{1}}{r_{2}^{6}}\left[3\left(x_{2}+\xi_{2}\right)^{2}-x_{1}^{2}\right]\right\} \\
#         &+\frac{G s_{1}}{2 \pi(1-\nu)}\left\{\frac{\left(x_{2}-\xi_{2}\right)\left[\left(x_{2}-\xi_{2}\right)^{2}-x_{1}^{2}\right]}{r_{1}^{4}}-\frac{\left(x_{2}+\xi_{2}\right)\left[\left(x_{2}+\xi_{2}\right)^{2}-x_{1}^{2}\right]}{r_{2}^{4}}\right.\\
#         -&\left.\frac{2 \xi_{2}}{r_{2}^{6}}\left[6 x_{2}\left(x_{2}+\xi_{2}\right) x_{1}^{2}-\left(3 x_{2}+\xi_{2}\right)\left(x_{2}+\xi_{2}\right)^{3}+x_{1}^{4}\right]\right\} \\ \\
#         \sigma_{12}=\frac{G s_{2}}{2 \pi(1-\nu)} &\left\{\frac{\left(x_{2}-\xi_{2}\right)\left[\left(x_{2}-\xi_{2}\right)^{2}-x_{1}^{2}\right]}{r_{1}^{4}}-\frac{\left(x_{2}+\xi_{2}\right)\left[\left(x_{2}+\xi_{2}\right)^{2}-x_{1}^{2}\right]}{r_{2}^{4}}\right.\\
#         +&\left.\frac{2 \xi_{2}}{r_{2}^{6}}\left[6 x_{2}\left(x_{2}+\xi_{2}\right) x_{1}^{2}-x_{1}^{4}+\left(\xi_{2}-x_{2}\right)\left(x_{2}+\xi_{2}\right)^{3}\right]\right\} \\
#         &+\frac{G s_{1}}{2 \pi(1-\nu)}\left\{\frac{x_{1}\left[\left(x_{2}-\xi_{2}\right)^{2}-x_{1}^{2}\right]}{r_{1}^{4}}-\frac{x_{1}\left[\left(x_{2}+\xi_{2}\right)^{2}-x_{1}^{2}\right]}{r_{2}^{4}}\right.\\
#         &\left.+\frac{4 \xi_{2} x_{2} x_{1}}{r_{2}^{6}}\left[3\left(x_{2}+\xi_{2}\right)^{2}-x_{1}^{2}\right]\right\}
#     \end{aligned}\\$}
#     $$ (ch09_eq39)
# 
#     All symbols are like in [Section 9.2.3](ch09-2-3). Note that, unlike displacements (Eq. {eq}`ch09_eq17`), stresses do not depend on the angles about the dislocation and its image ($\theta_1$ and $\theta_2$) but do depend on the shear modulus ($G$).
# 
#     - Write a function that calculates the stresses at observation points due to an edge dislocation. Hint: Add to the `disloc2d.py` module a function to calculate stress. Call this stress function from function `disloc2d`
# 
#     - Calculate stresses for the same fault geometry and grid of observation points that you used to plot vectors of fault displacement in notebook [ch9-3](https://github.com/nfcd/compGeo/blob/master/source/notebooks/ch9-3.ipynb). Use 10 GPa for the shear modulus and 2 m of reverse-sense slip. Make plots of each of the three stress components: $\sigma_{11}$, $\sigma_{12}$, and $\sigma_{22}$. There are multiple ways you could make this plot. A simple way would be to use Pyplot `scatter` function with stress as the `c` (color) argument. A more advanced, but nicer looking way, would be to use `imshow`.
# 
#     - Calculate the three strain components $\varepsilon_{11}$, $\varepsilon_{12}$, and $\varepsilon_{22}$, and make plots of each. Eqs {eq}`ch09_eq17` and {eq}`ch09_eq39` are defined for a state of plane strain, in which $\varepsilon_{33}$ = 0 but $\sigma_{33} \neq 0$. Solving Eq. {eq}`ch09_eq01` with $\varepsilon_{33} = 0$ and writing the result in terms of $G$ and $\nu$ gives the following equations for plane strain:
# 
#         $$\begin{gathered}
#             \varepsilon_{11}=\frac{1}{2 G}\left[(1-\nu) \sigma_{11}-\nu \sigma_{22}\right] \\ 
#             \varepsilon_{22}=\frac{1}{2 G}\left[(1-\nu) \sigma_{22}-\nu \sigma_{11}\right] \\
#             \varepsilon_{12}=\frac{1}{2 G} \sigma_{12}
#         \end{gathered}
#         $$ (ch09_eq40)
# 
#     - Change the sense of slip to normal and repeat the stress and  strain calculations. Compare the results to those from the reverse slip case.
# 
# 4. The example in notebook [ch9-4](https://github.com/nfcd/compGeo/blob/master/source/notebooks/ch9-4.ipynb) considers elastic wave propagation in a region consisting of two materials with different seismic velocities but homogeneous densities.
# 
#     - Change the elastic model, such that the seismic velocities are homogeneous, but the density is heterogeneous and larger inside the 90 corner.
# 
#     - Change the source from a punctual explosion to a downward force source.
# 
#     - Change the source to a lateral force source.

# (ch09-4)=
# ## References
# 
# Allmendinger, R.W. 2020. Modern Structural Practice: A structural
# geology laboratory manual for the 21st century.
# \[[Online](https://www.rickallmendinger.net/download)\]. \[Accessed
# March, 2021\].
# 
# Fossen, H. 2016. Structural Geology. Cambridge University Press.
# 
# Gudmunsson, A. 2011. Rock fractures in geological processes. Cambridge
# University Press.
# 
# Hetenyi, M. 1946. Beams on Elastic Foundations. Theory with Applications
# in the Fields of Civil and Mechanical Engineering. The University of
# Michigan Press.
# 
# Hoek, E. 2006. Practical Rock Engineering. Available from this
# [link](https://www.rocscience.com/assets/resources/learning/hoek/Practical-Rock-Engineering-Full-Text.pdf).
# 
# Jaeger, J.C., Cook, N.G.W. and Zimmerman, R.W. 2007. Fundamentals of
# Rock Mechanics, fourth edition. Blackwell Publishing.
# 
# Pollard, D.D. and Fletcher, R.C. 2005. Fundamentals of Structural
# Geology. Cambridge University Press.
# 
# Ragan, D.M. 2009. Structural Geology: An Introduction to Geometrical
# Techniques. Cambridge University Press.
# 
# Segall, P. 2010. Earthquake and Volcano Deformation. Princeton
# University Press.
# 
# Sheriff, R. and Geldart, L.P. 1995. Exploration Seismology, second
# edition. Cambridge University Press.
# 
# Stein, R.S. and Barrientos, S.E. 1985. Planar High-Angle Faulting in the
# Basin and Range: Geodetic Analysis of the 1983 Borah Peak, Idaho,
# Earthquake. Journal of Geophysical Research 90, 11,355-11,366.
# 
# Turcotte, D.L. and Schubert, G. 1982. Geodynamics. Jon Wiley $\&$ Sons.
# 
# Van Der Pluijm, B.A. and Marshak, S. 2004. Earth Structure, second
# edition. Norton.
# 
# Virieux, J. 1986. P-SV wave propagation in heterogeneous media:
# Velocity-stress finite-difference method. Geophysics 51, 889-901.
# 
# Watts, A.B. 2001. Isostasy and Flexure of the Lithosphere. Cambridge
# University Press.
# 
# Zoback, M.D. 2010. Reservoir Geomechanics. Cambridge University Press.
