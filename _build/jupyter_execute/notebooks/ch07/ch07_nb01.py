#!/usr/bin/env python
# coding: utf-8

# (ch07_nb01)=
# # Principal axes of stress

# In[1]:


# Import libraries
import numpy as np
pi = np.pi

# Import Cauchy, TransformStress and PrincipalStress
import sys, os
sys.path.append(os.path.abspath("../../"))

from compgeo.cauchy import cauchy 
from compgeo.transform_stress import transform_stress
from compgeo.principal_stress import principal_stress


# In[2]:


# Stress tensor in principal stress coordinate system
stress = np.array([[40, 0, 0],[ 0, 30, 0],[ 0, 0, 20]])

# trend and plunge of X1, and trend of X3
tx1, px1, tx3 = np.radians([0, 90, 0])

# plane orientation
strike, dip = np.radians([40, 65])

# X1, X2 and X3 tractions on the plane
t,pt = cauchy(stress,tx1,px1,tx3,strike,dip)
print("X1, X2 and X3 tractions = ", t.round(3),"\n")

# Compute the normal and maximum shear tractions 
# on the plane: Eq. 7.6
l2 = pt[0]**2
m2 = pt[1]**2
n2 = pt[2]**2
s1 = stress[0,0]
s2 = stress[1,1]
s3 = stress[2,2]
s12 = s1 - s2
s23 = s2 - s3
s31 = s3 - s1
sigma = s1*l2 + s2*m2 + s3*n2
tau = np.sqrt(s12*s12*l2*m2 + s23*s23*m2*n2 + s31*s31*n2*l2)
print("Sigma = {:.3f}, Tau = {:.3f}\n".format(sigma,tau))

# New coordinate system
# trend and plunge of X"1,and trend of X"3
ntx1, npx1, ntx3 = np.radians([30, 45, 210])

# Transform stress to new coordinate system
nstress = transform_stress(stress,tx1,px1,tx3,ntx1,npx1,ntx3)
print("Stress in new coord. system = \n", 
      nstress.round(3),"\n")

# Principal stresses from new components
pstress, dcp = principal_stress(nstress,ntx1,npx1,ntx3)
pstress[:,1:3] = pstress[:,1:3]*180/pi
print("Sigma1 = {:.3f}, Trend = {:.1f}, Plunge = {:.1f}"
      .format(pstress[0,0],pstress[0,1],pstress[0,2]))
print("Sigma2 = {:.3f}, Trend = {:.1f}, Plunge = {:.1f}"
      .format(pstress[1,0],pstress[1,1],pstress[1,2]))
print("Sigma3 = {:.3f}, Trend = {:.1f}, Plunge = {:.1f}"
      .format(pstress[2,0],pstress[2,1],pstress[2,2]))


# In[ ]:




