import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib.patches as patches
import csv
import os
import time
from scipy.linalg import solveh_banded
from scipy.sparse import csr_matrix
from scipy.sparse import csc_matrix 
from scipy.sparse.linalg import cg

penal = 3

volfrac = 0.5

np.set_printoptions(precision=4)
fig, ax = plt.subplots()

def timer(func):
    def wrapper(*args, **kwargs):
        t1 = time.time()
        mat = func(*args, **kwargs)
        t2 = time.time()
        print(f'{func.__name__} finished in: {t2-t1}')
        return mat
    return wrapper

script_dir = os.path.dirname(os.path.abspath(__file__))
coord_path = os.path.join(script_dir, 'Plane_coord.csv')
con_path = os.path.join(script_dir, 'Plane_connect.csv')
prop_path = os.path.join(script_dir, 'Plane_properties.csv')
bc_path = os.path.join(script_dir, 'Plane_bound_con.csv')
force_path = os.path.join(script_dir, 'Plane_force_list.csv')
sens_path = os.path.join(script_dir, 'SensitivityFilter.csv')

with open(prop_path,'r')as fl:
    props = []
    for i in csv.reader(fl):
        for j in i:
            props.append(float(j))
    E = props[0]
    y = props[1]
    Yield = props[2]
    t =props[3]

with open(coord_path, 'r') as file:
    try:
        val_list = []
        for i in csv.reader(file):
            val_list.append([float(j) for j in i])
        val_array= np.array(val_list)
        with open(con_path,'r') as con:
            con_list= []
            for c in csv.reader(con):
                con_list.append([int(b) for b in c])
            con_array= np.array(con_list)

    except Exception:
        print('error')

with open(bc_path, 'r') as bc:
        bc_lst = []
        c = 0
        for c1 in csv.reader(bc):
            pos = int(c1[0])*2
            if c1[1]=='x':
                if c1[0]==c:
                    continue
                bc_lst.append(pos)
            elif c1[1] == 'y':
                if c1[0]==c:
                    continue
                bc_lst.append(pos+1)

with open(force_path, 'r') as force:
        force_list =[]
        lst = csv.reader(force)
        for row in lst:
            if row == []: continue
            for i in row[2:]:
                force_list.append([int(i),[float(j) for j in row[:2]]])


t1 = time.time()
sumHf = []
Hf = []
findices = []

h1 = []
with open(sens_path, 'r') as sens:
    reader = csv.reader(sens)
    for i, sen in enumerate(reader):
        if i%2==0:
            f = np.array([[float(i) for i in sen]])
            sumHf.extend([np.sum(f)])
            Hf.extend(f)
        else:
            ind = np.array([[int(i) for i in sen]])
            findices.extend(ind)

t2 = time.time()

def bitmapPlotyy(cord, con, x, indice):
    for i, j in enumerate(con):
        triangle = np.array([cord[j[0]], cord[j[1]], cord[j[2]]])
        plt.text(np.sum(triangle[:,0])/3, np.sum(triangle[:,1])/3, s = f"{indice[i]}")

        gray_shade = 1-x[i]  # Change this value (0 to 1) to adjust grayscale intensity
        #h = float(i/len(con))*np.array([1.0,1.0,1.0])
        triangle = patches.Polygon(triangle, closed=True, color = (gray_shade, gray_shade, gray_shade))
        ax.add_patch(triangle)
        #ax.set_xlim(-2, 2)
        #ax.set_ylim(-2, 2)
        plt.axis("equal")
#ll=456
#bitmapPlotyy(val_array, con_array[findices[ll], :], Hf[ll], findices[ll])
#print(len(findices[ll]))
#plt.show()



#print(Hf[0:5], findices[0:5], t2-t1)

#print(force_list)
""" with open(force_path, 'r') as force:
        force_list =[]
        lst = csv.reader(force) 
        for row in lst:
            if row == []: continue
            for i in range(0, len(row[2:]), 2):
                force_list.append([[int(row[i+2]), int(row[i+3])],[float(j) for j in row[:2]]]) """
#print(np.array(force_list)[:])

#E, y, t =  1, 0.2, 1
#E, y, t =  30e6, 0.25, 1

E,y,t=1,0.3,1
D = (E/(1-y**2))*np.array([ [1, y,       0],
                            [y, 1,       0],
                            [0, 0, (1-y)/2]])


""" Ds = (E/((1+y)*(1-2*y)))*np.array([ [1-y,   y,         0],
                                    [  y, 1-y,         0],
                                    [  0,   0, (1-2*y)/2]]) """



def force_matrix(node,force,coord):
    f_matrix= np.zeros([2*len(coord),1])
    k = 2*node
    for i in range(2):
        f_matrix[k+i] += force[i]
    return f_matrix



#@timer

def cstPrime():
    cord=val_array
    con = con_array
    force_mat=force_matrix(0,[0.0,0.0], cord)

    for ac in force_list:
        force_mat += force_matrix(*ac, cord)


    def area(element):
        area = 0.5*np.linalg.det(np.array([[1, *cord[element[0]]],
                                            [1, *cord[element[1]]],
                                            [1, *cord[element[2]]]]))
        return abs(area)
    local_matrices = np.zeros([len(con),6,6])
    def local_stiffness_matrix(e):
        a = area(e)
        x1, x2, x3 = cord[e[0]][0], cord[e[1]][0], cord[e[2]][0]
        y1, y2, y3 = cord[e[0]][1], cord[e[1]][1], cord[e[2]][1]
        q1, q2, q3, r1, r2, r3 = y2-y3, y3-y1, y1-y2, x3-x2, x1-x3, x2-x1

        B = (0.5/a)*np.array([[q1,  0, q2,  0, q3,  0],
                                [ 0, r1,  0, r2,  0, r3],
                                [r1, q1, r2, q2, r3, q3]])
        
        #print(B)
        BtD = np.transpose(B)@D
        return a*BtD@B, B

@timer
def cst(densityMat):
    cord=val_array
    con = con_array
    force_mat=force_matrix(0,[0.0,0.0], cord)

    for ac in force_list:
        force_mat += force_matrix(*ac, cord)


    def area(element):
        area = 0.5*np.linalg.det(np.array([[1, *cord[element[0]]],
                                            [1, *cord[element[1]]],
                                            [1, *cord[element[2]]]]))
        return abs(area)
    local_matrices = np.zeros([len(con),6,6])
    def local_stiffness_matrix(e):
        a = area(e)
        x1, x2, x3 = cord[e[0]][0], cord[e[1]][0], cord[e[2]][0]
        y1, y2, y3 = cord[e[0]][1], cord[e[1]][1], cord[e[2]][1]
        q1, q2, q3, r1, r2, r3 = y2-y3, y3-y1, y1-y2, x3-x2, x1-x3, x2-x1

        B = (0.5/a)*np.array([[q1,  0, q2,  0, q3,  0],
                                [ 0, r1,  0, r2,  0, r3],
                                [r1, q1, r2, q2, r3, q3]])
        
        #print(B)
        BtD = np.transpose(B)@D
        return a*BtD@B, B

    #print(local_stiffness_matrix(con[1]))


    def global_stiff_mat(coord_list):
        matrix= np.zeros([len(cord)*2,len(cord)*2])
        for p,i in enumerate(coord_list):
            t1 = np.array([[i[0],i[0]], [i[0],i[1]], [i[0],i[2]], [i[1],i[0]], [i[1],i[1]], [i[1],i[2]], [i[2],i[0]], [i[2],i[1]], [i[2],i[2]]])
            t2= np.array([[0,0], [0,2], [0,4], [2,0], [2,2], [2,4], [4,0], [4,2], [4,4]])
            ki = local_stiffness_matrix(i)[0]
            local_matrices[p] = ki
            k = (densityMat[0,p]**penal)*ki
            
            for j in range(9):
                matrix[2*t1[j,0]:2*t1[j,0]+2 , 2*t1[j,1]:2*t1[j,1]+2] += k[t2[j,0]:t2[j,0]+2 , t2[j,1]:t2[j,1]+2]
        return matrix
    

    def boundary_conditions(matrix,bc_list):
        matrix= np.delete(matrix,bc_list,0)
        if len(matrix[0])!=1:
            matrix= np.delete(matrix,bc_list,1)
        return matrix


    global_stiff_matrix= global_stiff_mat(con)
    #print(global_stiff_matrix)
    #APPLYING BOUNDARY CONDITIONS
    bc_stiff_mat= boundary_conditions(global_stiff_matrix,bc_lst)
    #print(bc_stiff_mat)
    bc_force_matrix= boundary_conditions(force_mat,bc_lst)
    #print(bc_force_matrix)
    t1  =time.time()
    #bc_stiff_mat_inv= np.linalg.inv(bc_stiff_mat) #matrix_inverter(bc_stiff_mat)
    #MULTIPLYING INVERSE TO FORCE MATRIX
    
    def rnge():
        max_arr = [abs(i[0]-i[1]) for i in con_array ]
        return max(max_arr)*3+3


    bc_stiff_mat= boundary_conditions(global_stiff_matrix,bc_lst)
    bc_force_matrix= boundary_conditions(force_mat,bc_lst)
    #bandwidth = rnge()
    bc_stiff_mat_inv= np.linalg.inv(bc_stiff_mat) 

    #diagonals = np.array([np.append(np.diag(bc_stiff_mat, -i), np.zeros([1,i])[0]) for i in range(bandwidth)])

    #nodal_disp = solveh_banded(diagonals, bc_force_matrix, lower= True)

    nodal_disp= bc_stiff_mat_inv@bc_force_matrix

    t2 =time.time()

    for ij in bc_lst:
        nodal_disp= np.insert(nodal_disp,ij,0)
    return nodal_disp, local_matrices


#OPTIMALITY CRITERIA METHOD
def OC(conn, x, volfrac, dc):
    l1, l2, move = 0, 100000, 0.2
    while l2 - l1 > 1e-4:
        lmid = 0.5 * (l1 + l2)
        #xnew1 = x * np.sqrt(np.maximum(-dc / lmid,0))
        xnew1 = x * np.sqrt(abs(-dc / lmid))

        xnew = np.maximum(0.001, np.maximum(x - move, np.minimum(1.0, np.minimum(x + move, xnew1))))

        if np.sum(xnew) - volfrac*len(conn) > 0:
            #print("l1")
            l1 = lmid
        else:
            #print("l2")
            l2 = lmid
    
    return xnew

def bitmapPlot(cord, con, x):
    for i, j in enumerate(con):
        triangle = np.array([cord[j[0]], cord[j[1]], cord[j[2]]])
        gray_shade = 1-x[0,i]  # Change this value (0 to 1) to adjust grayscale intensity
        #h = float(i/len(con))*np.array([1.0,1.0,1.0])
        triangle = patches.Polygon(triangle, closed=True, color = (gray_shade, gray_shade, gray_shade))
        
        ax.add_patch(triangle)
        plt.axis("equal")
        #ax.set_xlim(-5, 5)
        #ax.set_ylim(-6, 6)



plt.ion()  # Turn on interactive mode
#fig, ax = plt.subplots()

#FINITE ELEMENT ANALYSIS, OBJECTIVE FUNCTION AND FINDING THE SENSITIVITIES 
X = volfrac*np.ones([1, len(con_array)])
dc = np.zeros([1,len(con_array)])

change=1
iter = 1
#for i in range(20):
while change>0.01:
    xold = X
    #break
    print(iter)
    U,K = cst(X)
    c=0
    
    for ind, k in enumerate(con_array):
        Ke = K[ind]
        xe = X[0,ind]
        Ue = np.array([U[2*i: 2*i+2]  for i in k]).ravel()
        UeTKeUe = (Ue@Ke@Ue.T)
        c = c + (xe**penal)*UeTKeUe
        dc[0,ind] = -penal*(xe**(penal-1))*UeTKeUe

    # SENSITIVITY FILTER UPDATION

    def filterE(indices, X, xe, Hf, sumHf, dc): #for the element e 
        sumHxdcf = 0
        for f, valf in enumerate(indices):
            sumHxdcf += np.maximum(0,Hf[f])*X[0,valf]*dc[0,valf]
        dcdxe = sumHxdcf/(xe*sumHf)
        return dcdxe
    
    dcnew = np.zeros([1,len(con_array)])
    for e, val in enumerate(dcnew[0,:]):
        dcnew[0,e] = filterE(findices[e], X, X[0,e], Hf[e], sumHf[e], dc)
    dc = dcnew


    # DESIGN UPDATE BY OPTIMALITY CRITERIA METHOD
    X = OC(con_array, X, volfrac, dc)
    change = np.max(np.abs(X-xold))
    if change<0.02:
        print("convergence")
    
    iter+=1

    ax.clear()
    cord=val_array
    bitmapPlot(cord, con_array, X)
    plt.pause(0.01)  # Pause for a short duration to update the figure

plt.ioff()  # Turn off interactive mode
plt.show()

#print(dc[0:10])
#bitmapPlot(cord, con_array, X)

plt.show()
