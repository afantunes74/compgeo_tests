#!/usr/bin/env python
# coding: utf-8

# (ch01_nb01)=
# # Introduction to Python and Jupyter Notebooks

# This is a quick introduction to the Python programming language and the Jupyter Notebook framework. Notebooks let you run Python code interactively in a browser. We assume no prior programming experience, though if you have used MATLAB previously, you will find the Python syntax very similar. 
# 
# This notebook will only cover the basics of Python, with emphasis on features needed for this resource. However, Python is a popular general-purpose programming language used in a wide variety of projects and industries. Did you know that Dropbox and Spotify are mainly written in Python? 

# ## Basics

# A notebook is divided into computational units called "cells". Cells can contain text such as this one or Python code. Below is a cell with some common Python statements. Try changing the variables and re-run the cell. To run a cell, either click the "Run" button above, or type Ctrl+Enter.

# In[1]:


a = 2
b = 9.0
c = a + b
print("The sum is: ", c)

# This is just a comment

name = "Donald"
print("Hello, my name is", name)


# There are some other useful shortcuts you should know. To run a cell and move to the next cell, type Shift+Enter. To run a cell and insert a new cell below, type Alt+Enter. You can use the arrow keys to move quickly between cells. To run all the cells of a notebook, choose the `Cell` -> `Run All` menu.

# ### Conditionals
# A conditional is used to decide between different operations depending on a conditional statement. In Python, this is expressed in the following way:

# In[2]:


a = 3
b = 5

if a > b:
    print("a is bigger than b")
elif a < b:    
    print("a is smaller than b")
else:
    print("a is equal to b")


# Try changing the values of `a` and `b` to see how the output changes. Also, note that Python cares about the indentation of the line, so there must be a tab indent or 4 spaces for each operation in the `if` statement. If not you will get an error as shown here: 
# 

# In[3]:


if a == b:
    c = 3


# You can also use the logical operators `and`, `or` and `not` in the conditional statement:

# In[4]:


age = 30
if age > 18 and age < 34:
    print("You are a young adult")

if age < 18 or age > 80:
    print("You are not allowed to drive a car")


# ### Loops
# A loop is used to execute a group of statements multiple times. For instance, to print all numbers from 1 to 10 divisible by 3, we can use a `for` loop together with an `if` statement, and the modulus operator, `%`:

# In[5]:


print("Numbers divisible by three:")
for i in range(1, 11):
    if i % 3 == 0:
        print(i)


# `range` is a Python function that iterates from the given first number up to the second number, but without including this number. If we only give one number, the iteration will go from zero up to the given number - 1. We will give more examples of `for` loops later in this notebook.

# ### Functions and modules

# If we have written a useful piece of code, we often want to use it again without copying and pasting the code multiple times. To do this, we use functions and modules. For instance, if we want to convert an angle from degrees to radians, we can use the following formula, 
# \begin{equation}
#     \alpha_\text{radians} = \alpha_\text{degrees}\frac{\pi}{180}
# \end{equation}    
# To put this into a callable function, we use the `def` keyword:

# In[6]:


def deg_to_rad(angle_degrees):
    pi = 3.141592
    return pi*angle_degrees/180.0

angle_degrees = 45.0
print("Radians", deg_to_rad(angle_degrees))


# We can also include code from other places. This is useful to make your own library of functions that you can then use in any code.This is the modus operandi of this resource. We will implement and use functions to solve interesting problems in geosciences. Using a text editor, create a file called mylib.py and put it in the same folder the notebook is. In the file, write two functions to convert from degrees to radians, and from radians to degrees:
# ```
# def deg_to_rad(angle_degrees):
#     pi = 3.141592
#     return pi*angle_degrees/180
#     
# def rad_to_deg(angle_radians):
#     pi = 3.141592
#     return angle_radians*180/pi
# ```
# A file like this is called a `module`, and it can contain one or several functions. We can then import in the notebook the module and use its functions like this:

# In[7]:


try:
    import mylib

    # degrees to radians
    angle_degrees = 45
    angle_radians = mylib.deg_to_rad(angle_degrees)
    print(angle_degrees, "degrees is", 
          angle_radians, "radians")

    # radians to degrees
    print(angle_radians, "radians is", 
          mylib.rad_to_deg(angle_radians), "degrees")
    
except ModuleNotFoundError:
    print("Create a file called mylib.py")


# _Note_: If you make a change in mylib.py, the changes will not be immediately available in the notebook and it needs to be restarted . To circumvent this, we can use the following commands to always reload imported modules:

# In[8]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '')


# ## Mathematics

# To use Python as an environment for numerical mathematics, it is useful to use the `numpy` library for arrays and matrices, and the `matplotlib` library for plotting. See the links in the `Help` menu above for more information on these libraries. The following two lines import these libraries. 

# In[9]:


import numpy as np
import matplotlib.pyplot as plt


# To define an array, we use the `numpy` `array` function:

# In[10]:


a = np.array( [1, 2, 3, 4] )
print(a)


# To access an array element, we use brackets with the index of the element. A very important difference compared to Matlab is that in Python, the first element has index zero (like most other programming languages). We can also use negative indices to access values starting from the end of the array.

# In[11]:


print(a[0], a[2])
print(a[-1])


# Slicing is a very useful feature to extract subarrays. For instance, 

# In[12]:


print(a[2:])
print(a[1:3])


# Matrices are defined as multi-dimensional arrays.

# In[13]:


a_matrix = np.array( [[1, 2, 3], 
                   [4, 5, 6], 
                   [7, 8, 9]] )
b_matrix = np.array( [[2, 4],
                   [3, 5],
                   [5, 7]] )
print(a_matrix)
print(b_matrix)


# We can get the number of rows and columns from the `shape` method:

# In[14]:


nrow, ncol = b_matrix.shape
print("b has {} rows and {} columns".format(nrow, ncol))


# Let us make a function to multiply two matrices together. Considering an $n\times m$ matrix $\mathbf A$ and a $m\times p$ matrix $\mathbf B$. The formula to multiply these matrices can be written as:
# 
# \begin{equation}
#  C_{ij} = \sum_{k=1}^m A_{ik} B_{kj}
# \end{equation}
# for $i = 1, ..., n$ and $j = 1, ..., p$. Here, $\mathbf C$ will be a $n\times p$ matrix.
# 
# To implement this formula, we need to use a triple-nested loop, as shown in the function below:

# In[15]:


def matrix_multiply(A,B):
    n, m = A.shape
    nrow_B, p = B.shape
    
    # Check that the matrices are conformable
    if not nrow_B == m:
        raise ValueError("Error: Number of columns in A" 
                         "must equal number of rows in B")
    
    # Initialize C using the numpy zeros function
    C = np.zeros((n, p))
    
    for i in range(n):
        for j in range(p):
            for k in range(m):
                C[i,j] = C[i,j] + A[i,k]*B[k,j]
                
    return C

print(matrix_multiply(a_matrix, b_matrix))


# Verify by hand calculation that the above result is correct. Remember, the element in the first row and first
# column of $\mathbf C$ is equal to the sum of the product of the elements in the first row of $\mathbf A$ times the elements in the first column of $\mathbf B$, and so on.
# 
# What happens if you try the multiplication $\mathbf B \mathbf A$? Try it in the cell below:

# In[ ]:





# Although the above function is elegant, it is not very efficient. The `numpy` library contains super-optimized code for common operations such as matrix multiplication. The `numpy` `dot` function can be used for matrix multiplication. Let's repeat the matrix multiplication above using the `dot` function:

# In[16]:


C = np.dot(a_matrix, b_matrix)
print(C)


# When working with large matrices, there is a significant impact on the runtime. To illustrate this, let's generate two 100x100 matrices with random numbers. The `numpy` `random.rand` function generates the arrays and fill them with random numbers.

# In[17]:


N = 100
A = np.random.rand(N,N)
B = np.random.rand(N,N)


# Now, let's measure the difference in execution time between multiplying these matrices with our function, or the `numpy` `dot` function. The `time` function allows us to measure the time taken by each function in seconds:

# In[18]:


import time
start = time.time() # start time
C = matrix_multiply(A,B)
endt = time.time() # end time
print("Our function time = ",(endt-start))


# In[19]:


start = time.time() # start time
C = np.dot(A,B)
endt = time.time() # end time
print("NumPy dot function time = ",(endt-start))


# The `numpy` `dot` function is much faster than our function!

# ### Plotting
# Arrays can be easily plotted using the `matplotlib` `plot` command. Below we plot the sinusoidal function. We use the `numpy` `linspace` function to generate an array with equally spaced values between the start and end point, and the `numpy` `sin` function to take the sine of the array. The `plot` command plots the data and the semicolon after it removes any message output. With a low number of points, the curve is actually jagged. Increase the number of points n in the `linspace` command to get a smoother curve. Try values of n = 100, 1000, and 10000.

# In[20]:


# The linspace command gives us an equally spaced array, the syntax is
# linspace(start_point, end_point, number_of_points)
n = 10
x = np.linspace(0, 10, n)
y = np.sin(x)
plt.plot(x, y)
plt.show()


# We end with a slightly more advanced plot, showing how to change line style and markers, and add axis labels and a legend. The `numpy` `cos` function takes the cosine of the array. The `plt.subplots` command makes a new figure and returns handles to the figure `fig` and axes `ax`. We can then use `ax` to plot the functions (`plot`), set the axes labels (`set_xlabel` and `set_ylabel`), and add a legend (`legend`). The last line saves the figure as a png image. Try also different values of `n` to see how the plot changes:

# In[21]:


n = 25
x = np.linspace(0, 10,n)
y = np.cos(x)
# Make a figure
fig, ax = plt.subplots()
# plot cosine function as a red line
ax.plot(x, y, "r", label="Cosine")
x = np.linspace(0,10,5)
y = 0.01*x**2
# plot quadratic funcion as blue dashed line with dots
ax.plot(x,y,"bo--", label="Quadratic")
# label axes
ax.set_xlabel("Time")
ax.set_ylabel("Amplitude")
# Add legend
ax.legend()
# Save the figure as a png with resolution 300 dpi
fig.savefig("my_first_plot.png", dpi=300)
plt.show()


# In[ ]:




