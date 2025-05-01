import bpy
import bmesh
import numpy as np
import time
import os
import csv
script_dir = os.path.dirname(os.path.abspath(__file__))
dir_path = os.path.dirname(script_dir)
sens_path = os.path.join(dir_path, 'SensitivityFilter.csv')
t1 = time.time()

ob = bpy.context.edit_object
me = ob.data
bm = bmesh.from_edit_mesh(me)
print("NEW")
rminFactor =3
recursionDepth = 2

it = 0
with open(sens_path, 'w') as clr:
   clr.write('') 
   
with open(sens_path, "a+", newline= '') as file:
    writer = csv.writer(file)
    for faceI in bm.faces:

        faceS = bm.select_history.active

        active_median = faceI.calc_center_median()
        rmin = rminFactor*(active_median - faceI.verts[0].co).length

        #face.edges takes less time but need more recursion
        #face.verts takes more time but need less recursion

        def surroundFaces(surround, faces, count, level):
            if count==level:
                return list(set(surround))
            else:
                for face in faces:
                    H = np.array([ f
                        for e in face.verts # and face.edges #face.edges is short but lengthy
                        for f in e.link_faces if f is not face 
                        #if (rmin - (f.calc_center_median() - active_median).length)<=0
                        ])  

                    surround = np.concatenate((surround, H))
                count +=1
                return surroundFaces(list(set(surround)), list(set(surround)), count, level)

        levelI = surroundFaces(np.array([]), [faceI], 0,recursionDepth)
        His = []

        for fj in levelI:
            His.extend([max(0.0, rmin - (fj.calc_center_median() - active_median).length)])
        

        indices = [i.index for i in levelI]
        if len(indices)!=len(His):

            print("NOT EQUAL")
            break
        
        writer.writerow(His)
        writer.writerow(indices)
        if it==faceS.index:
            for f in levelI:
                #continue
                f.select = True
            #break
        #print("no of faces;",len(indices))
        it+=1

 

t2 = time.time()
print("time taken", t2-t1)

    
    


surround1 = []
def sample():
    faceI = bm.select_history.active
    print(faceI.index, "Index")
    #print("NEW")
    active_median = faceI.calc_center_median()
    rmin = rminFactor*(active_median - faceI.verts[0].co).length
    Hi = []


    def surroundFaces(faces, count, level):
        if count==level:
            return list(set(surround1))
        else:
            for face in faces:

                link_faces_edges = [f
                    for e in face.verts 
                    for f in e.link_faces if f is not face
                    if (rmin - (f.calc_center_median() - active_median).length)>=0]
                surround1.extend(link_faces_edges)
            count +=1
            return surroundFaces(list(set(surround1)), count, level)

    levelN = []
    level = recursionDepth
    levelI = surroundFaces([faceI], 0,level)
    print('sample hi',len(levelI))

#sample()


bmesh.update_edit_mesh(me)
me.update()
for f in surround1:
    #continue
    f.select = True
    
    
