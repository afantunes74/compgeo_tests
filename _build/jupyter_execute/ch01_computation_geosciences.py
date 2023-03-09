#!/usr/bin/env python
# coding: utf-8

# (ch01)=
# # Computation in Geosciences

# (ch01-1)=
# ## Solving problems by computation
# 
# Geology is an interpretive and historical science (Frodeman, 1995). We observe, collect, analyze, and interpret data (what) to tell a story (why). To collect data, we need to take measurements. All measurements have some uncertainty, and therefore uncertainty and error propagation are important in geosciences, and they are a recurring topic in this resource.
# 
# For the last 50 years or more, the methods geoscientists have used to visualize, analyze and interpret data are mostly graphical. For example, in structural geology, students typically learn two types of graphical constructions: orthographic and spherical projections (stereonets) (Ragan, 2009). Although these methods are great to visualize and solve geometrical problems in three-dimensions, they are not amenable to computation, and therefore applying these methods to large datasets with thousands of entries is impractical. Plane and spherical trigonometry allow deriving formulas (e.g. apparent dip formula) for computation (Ragan, 2009). However, these formulas give little insight about the problems. They are just formulas associated with complex geometric constructions, which bear no relation to each other, and which are difficult to combine to solve more complicated problems.
# 
# It turns out that many of the most interesting problems in geosciences can be solved using linear algebra, and vectors and tensors (Allmendinger et al., 2012). Linear algebra also happens to be the language of computation. The main purpose of this resource is to show how to solve problems in geosciences using computation. There are several advantages of following this approach. It will enhance your mathematical and computational skills, as well as promote your geological-mathematical problem solving disposition. In today’s digital age, these skills are very useful.

# (ch01-2)=
# ## Why Python?
# 
# The choice of programming language is important. While computer languages such as C or C++ are ideal to work with large datasets and computer-intensive operations, they involve a steep learning curve associated with their syntax, compilation, and execution (Jacobs et al., 2016). These coding details have little to do with the problem-solving approach of this resource. Interpretive languages such as Python, R or Matlab are a better choice because of their simpler syntax, and the interpretation and execution of commands as they are called (no need for compilation). In addition, these languages have access to an integrated development environment (IDE) that facilitates writing and debugging programs, and to many standard libraries that perform advanced tasks such as matrix operations and data visualization. Thus, Python, R or Matlab are “scientific packages” rather than just programming languages.
# 
# In this resource, the language of choice is Python. Besides the reasons above, Python has the following advantages:
# 
# -   Python can be learned quickly. It typically involves less code than other languages and its syntax is easier to read.
# -   Python comes with robust standard libraries for arrays and mathematical functions (NumPy), visualization (Pyplot), and scientific computing (SciPy).
# -   As of December 2022, Python is the most popular programming language followed by C, C++, and Java ([TIOBE index](https://tiobe.com/tiobe-index/)), with a large base of developers and users. It is used by every major technology company and it is almost a skill you must have in your CV to land a job as a geoscientist.
# -   Because of its large developers base, Python has access to a large amount of external libraries, including several libraries for geosciences. We make use of some of these libraries in this resource.
# -   Python can be installed easily through a single distribution that includes all the standard libraries and provides access to external libraries (see next section).
# -   Last but not least, Python is free and open source. This is probably why Python is more popular than its commercial counterparts.

# (ch01-3)=
# ## Installing Python
# 
# We recommend installing Python using the free Anaconda distribution. This distribution includes Python as well as many other useful applications, including Jupyter, which is the system we use to write the notebooks in this resource. Anaconda can be easily installed on any major operating system, including Windows, macOS, or Linux.
# 
# The installation process is quite straightforward. From the [Anaconda](https://www.anaconda.com) page, press the `Download` button. Windows and macOS users just need to download the Graphical installer, run it, and follow the steps to install Anaconda. Linux users need to download the installer and type a set of commands in a terminal window. Further instructions can be found in the [Anaconda installation](https://docs.anaconda.com/anaconda/install/) section.

# (ch01-4)=
# ## A first introduction to Python
# 
# In this section, we use our first Jupyter notebook to learn the basics of Python. Clone or download the resource [git repository](https://github.com/nfcd/compGeo). Open Anaconda and then launch Jupyter Notebook. This will open a browser with a list of files and folders in your home directory. Navigate to the folder `/notebooks` in the repository and open the notebook [ch1](https://github.com/nfcd/compGeo/blob/master/source/notebooks/ch1.ipynb). Alternatively, follow the notebook in the sections below. Surprisingly, few lines of code are required to introduce key topics such as conditionals, loops, functions, array mathematics, and plotting. This shows the power of Python.
# 
# :::{note}
# Notice that in the notebook we concentrate on features relevant to this resource. There is much more to learn about Python and there are very good online and book resources to do that.
# :::

# (ch01-4-1)=
# ### Basics
# 
# A notebook is divided into computational units called *cells*. Cells can contain text such as this one or Python code. Below is a cell with some typical Python statements.. Try changing the variables and re-run the cell. To run a cell, either click the `Run` button, or type `Ctrl+Enter`.
# 
# :::{note}
# In this book, a gray box with a green ribbon on the left indicates a code block; while a plain gray box indicates the output of a code block.
# :::

# In[1]:


a = 2
b = 9.0
c = a + b
print("The sum is: ", c)


# In[2]:


# this is just a comment
name = "Donald"
print("Hello, my name is ", name)


# There are some other useful shortcuts you should know. To run a cell and move to the next cell, type `Shift+Enter`. To run a cell and insert a new cell below, type `Alt+Enter`. You can use the arrow keys to move quickly between cells. To run all the cells of a notebook, choose the `Cell ` $\rightarrow$ ` Run All` menu.

# (ch01-4-2)=
# ### Conditionals
# 
# A conditional is used to decide between different operations depending on a conditional statement. In Python, this is expressed in the following way:

# In[3]:


a = 3
b = 5

if a > b:
    print("a is bigger than b")

elif a < b:    
    print("a is smaller than b")

else:
    print("a is equal to b")


# Try changing the values of `a` and `b` to see how the output changes. Also, note that Python cares about the indentation of the line, so there must be a tab indent or 4 spaces for each operation in the if statement. You can also use the logical operators `and`, `or`, and `not` in the conditional statement:

# In[4]:


age = 30

if age > 18 and age < 34:
    print("You are a young adult")

if age < 18 or age > 80:
    print("You are not allowed to drive a car")


# (ch01-4-3)=
# ### Loops
# 
# A loop is used to execute a group of statements multiple times. For instance, to print all numbers from 1 to 10 divisible by 3, we can use a `for` loop together with an `if` statement, and the modulus operator `%`:

# In[5]:


print("Numbers divisible by three:")

for i in range(1, 11):
    if i % 3 == 0:
        print(i)


# `range` is a Python function that iterates from the given first number up to the second number, but without including this number. If we only give one number, the iteration will go from zero up to the given number $-1$. We give more examples of `for` loops later in this notebook.

# (ch01-4-4)=
# ### Functions and modules
# 
# If we have written a useful piece of code, we often want to use it again without copying and pasting the code multiple times. To do this, we use functions and modules. For instance, if we want to convert an angle from degrees to radians, we can use the following formula:
# 
# $$
# \alpha_\text{radians} = \alpha_\text{degrees}\frac{\pi}{180}
# $$ (ch01_eq01)
# 
# :::{note}
# This example is just for demonstration. Python has specialized functions to convert from degrees to radians (`math.radians`) and vice-versa (`math.degrees`).
# :::
# 
# To put this into a callable function, we use the `def` keyword:

# In[6]:


def deg_to_rad(angle_degrees):
    pi = 3.141592
    return pi * angle_degrees/180.0


# In[7]:


angle_degrees = 45.0
print("Radians", deg_to_rad(angle_degrees))


# We can also include code from other places. This is useful to make your own library of functions that you can then use in any code. This is the modus operandi of this resource. We will implement and use functions to solve problems in geosciences. Using a text editor, create a file called `mylib.py` and put it in the same folder the notebook is. In the file, write two functions to convert from degrees to radians, and from radians to degrees:

# In[8]:


def deg_to_rad(angle_degrees):
    pi = 3.141592
    return pi * angle_degrees/180

def rad_to_deg(angle_radians):
    pi = 3.141592
    return angle_radians * 180/pi


# A file like this is called a *module*, and it can contain one or several functions. We can then import in the notebook the module and use its functions like this:

# In[9]:


try:
    import mylib

    # degrees to radians
    angle_degrees = 45
    angle_radians = mylib.deg_to_rad(angle_degrees)
    print(angle_degrees, "degrees is",  angle_radians, "radians")

    # radians to degrees
    print(angle_radians, "radians is", mylib.rad_to_deg(angle_radians), "degrees")
    
except ModuleNotFoundError:
    print("Create a file called mylib.py")


# :::{note}
# If you make a change in `mylib.py`, the changes will not be immediately available in the notebook and it needs to be restarted. To circumvent this, we can use the following commands to always reload imported modules:
# 
# ``` bash
# %autoreload
# ```
# :::

# (ch01-4-5)=
# ### Mathematics
# 
# To use Python as an environment for numerical mathematics, it is useful to use the NumPy library for arrays and matrices, and the Matplotlib library for plotting. See the links in the Jupyter Notebook `Help` menu for more information on these libraries. The following two lines import these libraries:

# In[10]:


import numpy as np
import matplotlib.pyplot as plt


# To define an array, we use the NumPy `array` function:

# In[11]:


a = np.array([1, 2, 3, 4])
print(a)


# To access an array element, we use brackets with the index of the element. A very important difference compared to Matlab is that in Python, the first element has index zero (like most other programming languages). We can also use negative indices to access values starting from the end of the array.

# In[12]:


print(a[0], a[2])
print(a[-1])


# Slicing is a very useful feature to extract subarrays. For instance:

# In[13]:


print(a[2:])
print(a[1:3])


# Matrices are defined as multi-dimensional arrays:

# In[14]:


a_matrix = np.array([[1, 2, 3], 
                     [4, 5, 6], 
                     [7, 8, 9]])

b_matrix = np.array([[2, 4],
                     [3, 5],
                     [5, 7]])

print(a_matrix)
print(b_matrix)


# :::{note}
# In Python, a long line within parenthesis can be broken into several lines.
# :::
# 
# We can get the number of rows and columns of the matrix from the `shape` method:

# In[15]:


nrow, ncol = b_matrix.shape
print("b has {} rows and {} columns".format(nrow, ncol))


# Let us make a function to multiply two matrices. Consider a $n\times m$ matrix **A** and a $m\times p$ matrix **B**. The formula to multiply these matrices can be written as:
# 
# $$
# C_{ij} = \sum_{k=1}^m A_{ik}B_{kj}
# $$ (ch01_eq02)
# 
# for $i = 1,...,n$ and $j = 1,...,p$. Here **C** will be a $n\times p$ matrix. To implement this formula, we need to use a triple-nested loop, as shown in the function below:

# In[16]:


def matrix_multiply(A, B):
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


# :::{note}
# In Python, a long string within parenthesis can be split as shown in the code above.
# :::

# In[17]:


print(matrix_multiply(a_matrix, b_matrix))


# Verify by hand calculation that the above result is correct. Remember, the element in the first row and first column of **C** is equal to the sum of the product of the elements in the first row of **A** times the elements in the first column of **B**, and so on. Try the multiplication **BA**. What happens?
# 
# Although the function above is elegant, it is not very efficient. The NumPy library contains super-optimized code for common operations such as matrix multiplication. The NumPy `dot` function can be used for matrix multiplication. Let’s repeat the matrix multiplication above using this function:

# In[18]:


C = np.dot(a_matrix, b_matrix)
print(C)


# When working with large matrices, there is a significant impact on the runtime. To illustrate this, let’s generate two 100 $\times$ 100 matrices with random numbers. The NumPy `random.rand` function generates the arrays and fill them with random numbers.

# In[19]:


N = 100
A = np.random.rand(N, N)
B = np.random.rand(N, N)


# Now, let’s measure the difference in execution time between multiplying these matrices with our function, or the NumPy `dot` function. The `time` function allows us to determine the time taken by each function in seconds:

# In[20]:


import time


start = time.time() # start time
C = matrix_multiply(A, B)
endt = time.time() # end time
print("Our function time = ", (endt - start))


# In[21]:


start = time.time() # start time
C = np.dot(A, B)
endt = time.time() # end time
print("NumPy dot function time = ",(endt - start))


# The NumPy `dot` function is much faster than our function!

# (ch01-4-6)=
# ### Plotting
# 
# Arrays can be easily plotted using the Matplotlib `plot` command. Below we plot the sine function. We use the NumPy `linspace` function to generate an array with equally spaced values between the start and end point, and the NumPy `sin` function to take the sine of the array. The `plot` command plots the data and the semicolon after it removes any message output. With a low number of points, the curve looks jagged. Increase the number of points `n` to get a smoother curve. Try values of `n` = 100, 1000, and 10000.

# In[22]:


# The linspace command gives us an equally spaced array
# The syntax is: 
# linspace(start_point, end_point, number_of_points)
n = 10
x = np.linspace(0, 10, n)
y = np.sin(x)
plt.plot(x, y);


# We end with a slightly more advanced plot, showing how to change line style and markers, and add axis labels and a legend. The NumPy `cos` function takes the cosine of the array. The `plt.subplots` command makes a new figure and returns handles to the figure `fig` and axes `ax`. We can then use `ax` to plot the functions (`plot`), set the axes labels (`set_xlabel` and `set_ylabel`), and add a legend (`legend`). The last line saves the figure as a png image. Try also different values of `n` to see how the plot changes:

# In[23]:


n = 25
x = np.linspace(0, 10, n)
y = np.cos(x)

# Make a figure
fig, ax = plt.subplots()

# plot cosine function as a red line
ax.plot(x, y, "r", label="Cosine")
x = np.linspace(0,10,5)
y = 0.01 * x**2

# plot quadratic funcion as blue dashed line with dots
ax.plot(x, y, "bo--", label="Quadratic")

# label axes
ax.set_xlabel("Time")
ax.set_ylabel("Amplitude")

# Add legend
ax.legend()

# Save the figure as a png with resolution 300 dpi
fig.savefig("data/output/my_first_plot.png", dpi=300);


# (ch01-5)=
# ## Exercises
# 
# 1. Write a program that prints each number from 1 to 20 on a new line. For each multiple of 3, print "Fizz" instead of the number. For each multiple of 5, print "Buzz" instead of the number. For numbers which are multiples of both 3 and 5, print "FizzBuzz" instead of the number. The correct answer is: 1 2 Fizz 4 Buzz Fizz 7 8 Fizz Buzz 11 Fizz 13 14 FizzBuzz 16 17 Fizz 19 Buzz.
# 
# :::{hint}
# You will need to use a loop and conditionals to solve this problem.
# :::
# 
# 2. Use the module `mylib.py` to convert the following angles in degrees to radians and vice versa: 0, 45, 90, 135, 180, 225, 270, 315, 360. Print the results in a neatly way.
# 
# 3. Given two 3 $\times$ 3 matrices **A** = [[1, 2, 3], [4, 5, 6], [7, 8, 9]] and **B** = [[5, 7, 2], [3, 5, 1], [2, 4, 3]], compute:
# 
#     1. The sum of the matrices (**A** + **B**),
#     
#     2. The difference of the matrices (**A** - **B**),
#     
#     3. The product of the matrices (**AB**),
#     
#     4. The sum of all elements of matrix **B**,
#     
#     5. The column sum of matrix **A**,
#     
#     6. The row sum of matrix **A**,
#     
#     7. The transpose $(\textbf{A}^T)$ of matrix **A**,
#     
#     8. The product $\textbf{AA}^T$. What is this product equal to?
# 
# 	:::{tip}
# 	Check the functions `add`, `subtract`, `dot`, `sum`, and `transpose` in the NumPy library.
# 	:::
# 
# 4.  The apparent dip $\alpha$ of a plane is given by the equation $\tan\alpha=\tan\delta\sin\beta$, where $\delta$ is the true dip of the plane, and $\beta$ is the orientation of the vertical profile along which the dip is measured ({numref}`Figure %s <ch03_fig01>`b).
# 
#     1. Make a function to compute the apparent dip $\alpha$ from the true dip $\delta$ and the orientation of the profile $\beta$. Angles should be entered in radians.
# 
#     2. Use this function to make a graph of profile orientation $\beta$ (0 to 90) versus apparent dip $\alpha$ (0 to 90), for values of true dip $\delta$ of 10, 20, 30, 40, 50, 60, 70, and 80. The graph should look like {numref}`Figure %s <ch01_fig01>` below.
# 	
# 	:::{tip}
# 	You need to use the NumPy and Matplotlib libraries. This problem is hard, don’t give up.
# 	:::

# ```{figure} /figures/ch01_fig01.png
# :width: 450px
# :name: ch01_fig01
# 
# Apparent dip $\alpha$ as function of the profile orientation $\beta$ and true dip $\delta$.
# ```

# (ch01-6)=
# ## References
# 
# Allmendinger, R.W., Cardozo, N. and Fisher, D.W. 2012. Structural
# Geology Algorithms: Vectors and Tensors. Cambridge University Press.
# 
# Frodeman, R. 1995. Geological reasoning: Geology as an interpretive and
# historical science. GSA Bulletin 107, 960-968.
# 
# Jacobs, C.T., Gorman, G.J., Rees, H.E. and Craig, L.E. 2016. Experiences
# with efficient methodologies for teaching computer programming to
# geoscientists. Journal of Geological Education 64, 183-198.
# 
# Ragan, D.M. 2009. Structural Geology: An Introduction to Geometrical
# Techniques. Cambridge University Press.

# In[ ]:




