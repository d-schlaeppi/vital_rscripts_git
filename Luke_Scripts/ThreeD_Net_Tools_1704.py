#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import numpy as np
import networkx as nx
import matplotlib.cm as cm
import matplotlib
import argparse
import matplotlib.pyplot as plt
import itertools
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import math


# In[ ]:


def nodez_1 (df): #this gets node from df - does not exclude by length
    list_node = []
    for index, row in df.iterrows():
        x = row['x']
        y = row['y']
        z = row['z']
        node = [x, y, z]
        list_node.append(node)
    array_node = np.array(list_node)
    return list_node, array_node


# In[ ]:


def nodez (list_edge):#use this after short edge exclusion
    list_node = []
    for e in list_edge:
        node1 = tuple(e[0])
        node2 = tuple(e[1])
        list_node.append(node1)
        list_node.append(node2)
    list_node = [*set(list_node)] #remove duplicate
    array_node = np.array(list_node)
    return list_node, array_node


# In[3]:


def edgez_1 (df):
    list_edge = []
    for index, row in df.iterrows():
        x1 = row['V1 x']
        y1 = row['V1 y']
        z1 = row['V1 z']
        x2 = row['V2 x']
        y2 = row['V2 y']
        z2 = row['V2 z']
       # end1=np.array([x1, y1, z1])
        length = row['Branch length']
        end1= [x1, y1, z1]
        end2= [x2 ,y2, z2]
        length=[length]
        edge = [end1, end2, length]
        list_edge.append(edge)
    array_edge = np.array(list_edge)
    return list_edge, array_edge


# In[4]:
def convert_cm_to_px (img_stack, x_cm, y_cm, z_cm):
    z_pix=len(img_stack)
    x_pix=len(img_stack[0][0])
    y_pix=len(img_stack[0])
    x_scale=x_cm/x_pix
    y_scale=y_cm/y_pix
    z_scale=z_cm/z_pix
    z_xy_scaling=z_scale/x_scale 
    return x_scale, y_scale, z_scale, z_xy_scaling

def tuple_convert (edge_list): #this converts a list of list format edge list to a nested tuple list
    tuple_edge_list = []
    for i in edge_list:
        end1=tuple(i[0])
        end2=tuple(i[1])
        length=tuple(i[2])
        edge=tuple([end1, end2, length])
        tuple_edge_list.append(edge)
    return tuple_edge_list


# In[5]:


def list_convert (tuple_edge_list): #this converts from tuple to list
    edge_list = []
    for i in tuple_edge_list:
        end1=list(i[0])
        end2=list(i[1])
        length=list(i[2])
        edge=[end1, end2, length]
        edge_list.append(edge)
    array_edge_list=np.array(edge_list)
    return edge_list, array_edge_list


# In[6]:


#this calculates the distance between to xyz nodes
def calculate_distance (node1, node2):
    np_node1=np.array(node1)
    np_node2=np.array(node2)
    squared_dist = np.sum((np_node1-np_node2)**2, axis=0) #this calculates distance between two points
    dist = np.sqrt(squared_dist)
    return dist


# In[7]:


#this finds a new midpoint node, between two old nodes
def midpoint_node (node1, node2):
    node1_x = node1[0] 
    node1_y = node1[1] 
    node1_z = node1[2] 
    node2_x = node2[0] 
    node2_y = node2[1] 
    node2_z = node2[2] 

    x = np.sum((node1_x+node2_x)/2) #c[0][0] -c[1][0]
    y = np.sum((node1_y+node2_y)/2)
    z = np.sum((node1_z+node2_z)/2)
    midpoint_node = (x, y, z)
    return midpoint_node


# In[8]:


#this prunes the nework/skeleton based on euclidean length
def prune_euclidean (edge_list, prune_length):
    node_list=[]
    pruned_edges=[]
    edge_list_pruned=[]

    for edge in edge_list:
        node1=edge[0]
        node2=edge[1]
        node_list.append(node1)
        node_list.append(node2)
    for edge in edge_list:
        branch_length=edge[2][0]
        node1=edge[0]
        node2=edge[1]
        dist1=calculate_distance(node1, node2)
        for c in range(1):
            node=edge[c]
            if node_list.count(node)==1:
                if dist1<prune_length:
                    pruned_edges.append(edge)
    for edge in edge_list:
        if edge not in pruned_edges:
            edge_list_pruned.append(edge)

    return edge_list_pruned


# In[9]:


#this prunes the nework/skeleton based on branch length (not euclidean dist)
def prune_branch (edge_list, prune_length):
    node_list=[]
    pruned_edges=[]
    edge_list_pruned=[]

    for edge in edge_list:
        node1=edge[0]
        node2=edge[1]
        node_list.append(node1)
        node_list.append(node2)
    for edge in edge_list:
        branch_length=edge[2][0]
        for c in range(1):
            node=edge[c]
            if node_list.count(node)==1:
                if branch_length<prune_length:
                    pruned_edges.append(edge)
    for edge in edge_list:
        if edge not in pruned_edges:
            edge_list_pruned.append(edge)

    return edge_list_pruned


# In[10]:


#removes edges below a certain euclidean (straight line length)
def reduce_network_euclidean (tuple_edge_list, length, its): 
    l=0
    edge_list2=tuple_edge_list
    edge_list3=edge_list2
    
    while l<its:
        l+=1
        new_edge_list=[]
        iterated_nodes=[]
        for edge1 in edge_list2:
            new_edge_list=[] #may even want in second for loop
            edge1_node1=edge1[0]
            edge1_node2=edge1[1]
            edge1_length=edge1[2][0]
            edge1_length_half=(edge1_length/2)
            dist1=calculate_distance(edge1_node1, edge1_node2)
            if dist1==0:
                edge_list2.remove(edge1)
            if dist1<length:
                deleted_edge=edge1
                iterated_nodes.append(edge1_node1)
                iterated_nodes.append(edge1_node2)
                new_node=midpoint_node(edge1_node1, edge1_node2)    
                break
        for edge2 in edge_list2:
            if edge2 != edge1:
                if edge2[0] in iterated_nodes:
                    kept_node=edge2[1]
                    #distance calculation
                    edge2_length=edge2[2][0]
                    #print(edge2_length)
                    new_edge_length=tuple([edge2_length+edge1_length_half])
                    #print(new_edge_length)
                    new_edge=(kept_node,new_node,new_edge_length)
                    new_edge_list.append(new_edge)                
                elif edge2[1] in iterated_nodes: 
                    kept_node=edge2[0]
                    edge2_length=edge2[2][0]
                    new_edge_length=tuple([edge2_length+edge1_length_half])
                    new_edge=(kept_node,new_node,new_edge_length)
                    new_edge_list.append(new_edge)
        for edge3 in edge_list2:
            if edge3[0] not in iterated_nodes:
                if edge3[1] not in iterated_nodes:
                    new_edge_list.append(edge3)
        edge_list2=new_edge_list
    edge_list_duplicates_removed=[]
    for edge in edge_list2:
        if edge not in edge_list_duplicates_removed:
            edge_list_duplicates_removed.append(edge)

    return edge_list_duplicates_removed
            
                


# In[11]:


#removes edges below a certain branch (true) length
def reduce_network_branch (tuple_edge_list, length, its): 
    l=0
    edge_list2=tuple_edge_list
    edge_list3=edge_list2
    
    while l<its:
        l+=1
        new_edge_list=[]
        iterated_nodes=[]
        for edge1 in edge_list2:
            new_edge_list=[] #may even want in second for loop
            edge1_node1=edge1[0]
            edge1_node2=edge1[1]
            edge1_length=edge1[2][0]
            edge1_length_half=(edge1_length/2)
            dist1=calculate_distance(edge1_node1, edge1_node2)
            if dist1==0:
                edge_list2.remove(edge1)
            if edge1_length<length:
                deleted_edge=edge1
                iterated_nodes.append(edge1_node1)
                iterated_nodes.append(edge1_node2)
                new_node=midpoint_node(edge1_node1, edge1_node2)    
                break
        for edge2 in edge_list2:
            if edge2 != edge1:
                if edge2[0] in iterated_nodes:
                    kept_node=edge2[1]
                    #distance calculation
                    edge2_length=edge2[2][0]
                    #print(edge2_length)
                    new_edge_length=tuple([edge2_length+edge1_length_half])
                    #print(new_edge_length)
                    new_edge=(kept_node,new_node,new_edge_length)
                    new_edge_list.append(new_edge)                
                elif edge2[1] in iterated_nodes: 
                    kept_node=edge2[0]
                    edge2_length=edge2[2][0]
                    new_edge_length=tuple([edge2_length+edge1_length_half])
                    new_edge=(kept_node,new_node,new_edge_length)
                    new_edge_list.append(new_edge)
        for edge3 in edge_list2:
            if edge3[0] not in iterated_nodes:
                if edge3[1] not in iterated_nodes:
                    new_edge_list.append(edge3)
        edge_list2=new_edge_list
    edge_list_duplicates_removed=[]
    for edge in edge_list2:
        if edge not in edge_list_duplicates_removed:
            edge_list_duplicates_removed.append(edge)

    return edge_list_duplicates_removed
            
                


# In[12]:


def tuple_convert_nest (edge_list): #this nests the edge nodes in another tuple, so the whole edge can be referenced by
    #[i][0]
    tuple_edge_list = []
    for i in edge_list:
        end1=tuple(i[0])
        end2=tuple(i[1])
        length=tuple(i[2])
        edge=tuple([end1, end2])
        #edge_dist=tuple([length])
        edge_complete=tuple([edge, length])
        tuple_edge_list.append(edge_complete)
    return tuple_edge_list


# In[13]:


def tuple_denest (edge_list): #this denestsit a layer
    #[i][0]
    tuple_edge_list = []
    for i in edge_list:
        end1=tuple(i[0][0])
        end2=tuple(i[0][1])
        length=tuple(i[1])
        #edge_dist=tuple([length])
        edge_complete=tuple([end1,end2, length])
        tuple_edge_list.append(edge_complete)
    return tuple_edge_list


# In[ ]:


#for removing duplicated with different distances
def clean_duplicate_node (edge_list_nest):
    cleaned_edges=[]
    dirty_edges=[]
    edge_list_nest2=edge_list_nest
    for edge in edge_list_nest:
        used_edge=[]
        used_edge.append(edge)
        #edge_list_nest2=set[edge_list_nest-used_edge]
        just_edge=edge[0]
        node1=edge[0][0]
        node2=edge[0][1]
        #print(edge[0])
        path=np.array(edge[1])
        dist1=calculate_distance(node1,node2)
        if dist1==0:
            dirty_edges.append(edge)

               # print(edge, edge2,clean_edge)
    no_dup_edges=list(set(edge_list_nest) - set(dirty_edges))
    no_dup_edges.extend(cleaned_edges)
    
    return(no_dup_edges)

#
def find_double_edge (edge_list):
    edge_list_nest=tuple_convert_nest(edge_list)
    double_edges=[]
    edge_list_nest2=edge_list_nest
    for edge in edge_list_nest:
        used_edge=[]
        used_edge.append(edge)
        #edge_list_nest2=set[edge_list_nest-used_edge]
        just_edge=edge[0]
        node1=edge[0][0]
        node2=edge[0][1]
        #print(edge[0])
        path=np.array(edge[1])
        for edge2 in edge_list_nest2:
            node3=edge2[0][0]
            node4=edge2[0][1]
            path2=np.array(edge2[1])
            dist1=calculate_distance(node1,node3)
            dist2=calculate_distance(node2,node4)
            if dist1==0 and dist2==0 and path==path2:
                continue
            if dist1==0 and dist2==0:
                double_edges.append(edge)
    
               # print(edge, edge2,clean_edge)
    return double_edges, len(double_edges)


# In[21]:


def remove_intermediate_nodes_clean (edge_list, iterations): 
    #first make a node list to count node occurence
    l=0
    edge_list=tuple_convert_nest(edge_list)
    edge_list=clean_duplicates(edge_list)
    old_edge_list=[]
    iterated_node_list=[]
    previous_edge=[]
    while l<iterations:
        node_list=[]
        edge_list_no_double=[]
        old_node_list=[]
        edge_list_duplicates_removed=[]
        edge_list_no_dist=[]
        for edge in edge_list:
            if edge not in edge_list_duplicates_removed:
                edge_list_duplicates_removed.append(edge)
              #  print(l, edge)
        edge_list=edge_list_duplicates_removed
        for edge in edge_list:
            node1=edge[0][0]
            node2=edge[0][1]
            node_list.append(node1)
            node_list.append(node2)
            edge_length=edge[1][0]  
            edge_mid=midpoint_node(node1, node2)
            if node1!=node2:
                if edge not in edge_list_no_double:
                    edge_list_no_double.append(edge)
                    edge_list_no_dist.append(edge_mid)
        edge_list=edge_list_no_double

        for node in node_list:
            node_count=node_list.count(node)
            test_edge=[]
            if node_count==2:
               # print('yerp:',node)
                old_node=node
                old_node_list.append(old_node)
                new_edge=[]
                new_length=[]
                for edge1 in edge_list:
                    
                    node1=edge1[0][0]
                    node2=edge1[0][1]
                    edge1_length=edge1[1][0]
                    for c in range(2):
                        if edge1[0][c]==node:
                            #print(edge1)
                            test_edge.append(edge1)
                            if edge1[0][c]==node1:
                                if node2 in old_node_list:
                                    new_edge.append(node1)
                                else:
                                    new_edge.append(node2)
                                new_length.append(edge1_length)
                            if edge1[0][c]==node2:
                                if node1 in old_node_list:
                                    new_edge.append(node2)
                                else:
                                    new_edge.append(node1)
                                new_length.append(edge1_length)
                break
        #print(edge1, 'test', new_edge, 'node:',)
        if l>iterations:
            edge_list.append(combined_edge)
            break
        else:
            new_edge_list=[]
            if new_edge not in old_edge_list: 
                node1=new_edge[0]
                node2=new_edge[1]
                length1=new_length[0]
                length2=new_length[1]
                combined_length=tuple([length1+length2])
                new_edge_no_dist=tuple([node1, node2])
                new_edge_mid=midpoint_node(node1,node2)

                if node1 == node2:
                    old_edge_list.append(new_edge)
                    old_node_list.append(old_node) 
                combined_edge=(new_edge_no_dist,combined_length)
                previous_edge=combined_edge
                if edge[0][0]!=edge[0][1]:
                    if new_edge_mid not in edge_list_no_dist:
                        new_edge_list.append(combined_edge)
                       # print(combined_edge,len(edge_list), l)
                for edge in edge_list:
                    edge_no_dist=edge[0]
                    if edge[0][0] not in old_node_list:
                        if edge[0][1] not in old_node_list:
                            if edge[0][0]!=edge[0][1]:
                                new_edge_list.append(edge)
            edge_list=new_edge_list
            
        l+=1
    edge_list_duplicates_removed=[]
    for edge in edge_list:
        if edge not in edge_list_duplicates_removed:
            edge_list_duplicates_removed.append(edge)  
    edge_list_duplicates_removed=tuple_denest(edge_list_duplicates_removed)     
    return edge_list_duplicates_removed
    #This has no cleaning step, so multiple edges are allowed between nodes.  
          
def remove_intermediate_nodes (rm_loop, edge_list, iterations): 
    #first make a node list to count node occurence
    if not isinstance(rm_loop, bool):
        raise TypeError("rm_loops must be a boolean value (True or False) if you want loop removal")
    l=0
    edge_list=tuple_convert_nest(edge_list)
    edge_list=clean_duplicate_node(edge_list)
    #print(edge_list)
    old_edge_list=[]
    iterated_node_list=[]
    previous_edge=[]
    while l<iterations:
        node_list=[]
        edge_list_no_double=[]
        old_node_list=[]
        edge_list_duplicates_removed=[]
        edge_list_no_dist=[]
        for edge in edge_list:
            if edge not in edge_list_duplicates_removed:
                edge_list_duplicates_removed.append(edge)
              #  print(l, edge)
        edge_list=edge_list_duplicates_removed
        
        for edge in edge_list:
            node1=edge[0][0]
            node2=edge[0][1]
            node_list.append(node1)
            node_list.append(node2)
            edge_length=edge[1][0]  
            edge_mid=midpoint_node(node1, node2)
            if rm_loop==True:
                if node1!=node2:
                    if edge not in edge_list_no_double:
                        edge_list_no_double.append(edge)
                        edge_list_no_dist.append(edge_mid)
            elif rm_loop==False:
                if edge not in edge_list_no_double:
                    edge_list_no_double.append(edge)
                    edge_list_no_dist.append(edge_mid)

        edge_list=edge_list_no_double
        for node in node_list:
            node_count=node_list.count(node)
            test_edge=[]
            if node_count==2:
               # print('yerp:',node)
                old_node=node
                old_node_list.append(old_node)
                new_edge=[]
                new_length=[]
                for edge1 in edge_list:
                    
                    node1=edge1[0][0]
                    node2=edge1[0][1]
                    edge1_length=edge1[1][0]
                    for c in range(2):
                        if edge1[0][c]==node:
                            #print(edge1)
                            test_edge.append(edge1)
                            if edge1[0][c]==node1:
                                if node2 in old_node_list:
                                    new_edge.append(node1)
                                else:
                                    new_edge.append(node2)
                                new_length.append(edge1_length)
                            if edge1[0][c]==node2:
                                if node1 in old_node_list:
                                    new_edge.append(node2)
                                else:
                                    new_edge.append(node1)
                                new_length.append(edge1_length)
                break
        #print(edge1, 'test', new_edge, 'node:',)
        if l>iterations:
            edge_list.append(combined_edge)
            break
        else:
            new_edge_list=[]
            if new_edge not in old_edge_list: 
                node1=new_edge[0]
                node2=new_edge[1]
                length1=new_length[0]
                length2=new_length[1]
                combined_length=tuple([length1+length2])
                new_edge_no_dist=tuple([node1, node2])
                new_edge_mid=midpoint_node(node1,node2)

                if node1 == node2:
                    old_edge_list.append(new_edge)
                    old_node_list.append(old_node) 
                combined_edge=(new_edge_no_dist,combined_length)
                previous_edge=combined_edge
                if edge[0][0]!=edge[0][1]:
                    if new_edge_mid not in edge_list_no_dist:
                        new_edge_list.append(combined_edge)
                       # print(combined_edge,len(edge_list), l)
                for edge in edge_list:
                    edge_no_dist=edge[0]
                    if edge[0][0] not in old_node_list:
                        if edge[0][1] not in old_node_list:
                            if edge[0][0]!=edge[0][1]:
                                new_edge_list.append(edge)
            if len(new_edge_list)>0:
                edge_list=new_edge_list
            else:
                break
            
        l+=1
    edge_list_duplicates_removed=[]
    for edge in edge_list:
        if edge not in edge_list_duplicates_removed:
            edge_list_duplicates_removed.append(edge)  
    edge_list_duplicates_removed=tuple_denest(edge_list_duplicates_removed)     
    return edge_list_duplicates_removed


# In[22]:


def remove_length (concat_edge): #use this to make into an array without distances (for 3d visualization)
    edge_list=[]
    for e in concat_edge:
        end1 = e[0]
        end2 = e[1]
        edge = [end1, end2]
        edge_list.append(edge)
    edge_array = np.array(edge_list)
    
    return(edge_array)


# In[23]:


def remove_tuple_length (concat_edge): #use this to make into an array without distances (for 3d visualization)
    edge_list=[]
    for e in concat_edge:
        end1 = e[0]
        end2 = e[1]
        edge = tuple([end1, end2])
        edge_list.append(edge)
    
    return(edge_list)



# In[ ]:


def test (edge_list_intermediate_removed,l):
    node_list=[]
    edge_list2=edge_list_intermediate_removed
    

    for edge in edge_list_intermediate_removed:
        print('edge_count:', edge_list_intermediate_removed.count(edge))
        node1=edge[0]
        node2=edge[1]
        node_list.append(node1)
        node_list.append(node2)
        for edge1 in edge_list2:
            new_edge_list=[] #may even want in second for loop
            edge_node1=edge1[0]
            edge_node2=edge1[1]            
            edge1_node1=edge1[0]
            edge1_node2=edge1[1]
            dist1=calculate_distance(edge_node1, edge1_node2)
            dist2=calculate_distance(edge1_node1, edge_node2)
            dist3=calculate_distance(edge1_node2, edge1_node1)
            if dist1<l:
                print('shit')
            if dist2<l:
                print('shit')
            if dist3<l:
                print('shit')
            
    for node in node_list:
        print(node_list.count(node))
        if node_list.count(node)==2:
            print('shit')


# In[25]:


def nodes_by_degree (edge_list_intermediate_removed):
    terminal_node_list=[]
    internal_node_list=[]
    node_list=[]
    node_list_dup_rm=[]
    edge_list2=edge_list_intermediate_removed
    

    for edge in edge_list_intermediate_removed:
        node1=edge[0]
        node2=edge[1]
        node_list.append(node1)
        node_list.append(node2)

    for node in node_list:
        if node_list.count(node)==1:
            terminal_node_list.append(node)
        else:
            internal_node_list.append(node)
        if node not in node_list_dup_rm:
            node_list_dup_rm.append(node)
    terminal_node_list_dup_rm=[]
    internal_node_list_dup_rm=[]
    for node in terminal_node_list:
        if node not in terminal_node_list_dup_rm:
            terminal_node_list_dup_rm.append(node)
    for node in internal_node_list:
        if node not in internal_node_list_dup_rm:
            internal_node_list_dup_rm.append(node)
    return terminal_node_list_dup_rm, internal_node_list_dup_rm, node_list_dup_rm
    
    

# In[2212]:


def parse_xyz (df_att):
   # df_att["xyz"] = pd.to_numeric(df_att["xyz"])
    for index, row in df_att.iterrows():
    #for xyz in df_att['xyz']:
        xyz=row['xyz']
        value=[]
        x=[]
        y=[]
        z=[]
        l=0
        separator=''
        for i in xyz:
            if i == ' ':
                continue
            if i =='(':
                continue
            if i ==')':
                value_str=(separator.join(value))
                z=value_str
                break
            if i == 0 or 1 or 2 or 3 or 4 or 5 or 6 or 7 or 8 or 9 or '.' or ',':
                if i ==',':
                    if l==0:
                        value_str=(separator.join(value))
                        x=value_str
                        value=[]
                        l+=1
                        continue
                    if l==1:
                        value_str=(separator.join(value))
                        y=value_str
                        value=[]
                        l+=1
                        continue
                else:
                    value.append(i)
        tup_x=tuple([x])
        tup_y=tuple([y])
        tup_z=tuple([z])
        float_x=float(tup_x[0])
        float_y=float(tup_y[0])
        float_z=float(tup_z[0])

        tup_xyz=tuple([float_x,float_y,float_z])
        df_att.at[index,'xyz']=tup_xyz
        df_att.at[index,'x']=tup_x
        df_att.at[index,'y']=tup_y
        df_att.at[index,'z']=tup_z
        
    return df_att


# In[ ]:


#requires pixels coords and radius
def create_sphere(cx,cy,cz, r_z, r_x, r_y, z_xy_scaling, resolution=360):
    '''
    create sphere with center (cx, cy, cz) and radius r_z
    '''
    phi = np.linspace(0, 2*np.pi, 2*resolution)
    theta = np.linspace(0, np.pi, resolution)

    theta, phi = np.meshgrid(theta, phi)
    #to account for differences in different xy to z pixel scaling
    
    #r_z=round(r*z_xy_scaling)
    
    r_xy=r_z*z_xy_scaling
    r_xy = r_xy*np.sin(theta)
    #r_z=r_xy/z_xy_scaling
    x = cx + np.cos(phi) * r_xy
    y = cy + np.sin(phi) * r_xy
    z = cz + r_z * np.cos(theta)
    points  = np.vstack( (x, y, z)).T
    return np.stack([x,y,z])



# In[ ]:


def combine_sphere_xyz (sphere, x_scale, y_scale, z_scale): #misses px
    x_stack = sphere[0]
    y_stack= sphere[1]
    z_stack = sphere[2]
    xyz_list=[]
    xyz_px_list=[]
    for i in range(len(x_stack)):
        xs=x_stack[i]
        ys=y_stack[i]
        zs=z_stack[i]
        for ii in range(len(xs)):
            x=xs[ii]
            y=ys[ii]
            z=zs[ii]
            x_px=round(x/x_scale)
            y_px=round(y/y_scale)
            z_px=round(z/z_scale)
            xyz=tuple([x, y, z])
            xyz_px=tuple([x_px, y_px, z_px])
            xyz_list.append(xyz)
            xyz_px_list.append(xyz_px)
    return xyz_list, xyz_px_list


# In[ ]:


def combine_sphere_xyz_px (sphere): #for asphere made with px
    x_stack = sphere[0]
    y_stack= sphere[1]
    z_stack = sphere[2]
    xyz_list=[]
    xyz_px_list=[]
    for i in range(len(x_stack)):
        xs=x_stack[i]
        ys=y_stack[i]
        zs=z_stack[i]
        for ii in range(len(xs)):
            x=xs[ii]
            y=ys[ii]
            z=zs[ii]
            x_px=round(x)
            y_px=round(y)
            z_px=round(z)
            xyz=tuple([x, y, z])
            xyz_px=tuple([x_px, y_px, z_px])
            xyz_list.append(xyz)
            xyz_px_list.append(xyz_px)
    xyz_px_list = [*set(xyz_px_list)]
        #rounded_xyz_list.append(rounded_xyz)
#     xyz_px_list_duplicates_removed=[]
#     for xyz_pv in rxyz_list:
#         if rounded_xyz not in rounded_xyz_list_duplicates_removed:
#             rounded_xyz_list_duplicates_removed.append(rounded_xyz)
    return xyz_list, xyz_px_list


# In[ ]:


def create_shell_sphere (cx, cy, cz, r_inner, r_outer, z_xy_scaling, down_sampling):
    spheres=[]
    for radii in range(r_inner, r_outer, down_sampling):
        sphere= create_sphere(cx,cy,cz,radii, z_xy_scaling, resolution=360)
        xyz_list, xyz_px_list = combine_sphere_xyz_px(sphere)
        spheres.extend(xyz_px_list)
        
    array_shell=np.array(spheres)
    return spheres, array_shell


# In[ ]:


#give radii and cemters in voxels
def create_ellipse(cx,cy,cz, r_x, r_y, r_z):
    '''
    create sphere with center (cx, cy, cz) and radius r_z
    '''
    phi = np.linspace(0, 2*np.pi, 256).reshape(256, 1)
    theta = np.linspace(0, np.pi, 256).reshape(-1,256)

    theta, phi = np.meshgrid(theta, phi)
    x = cx + r_x*np.sin(theta)*np.cos(phi)
    y = cy + r_y*np.sin(theta)*np.sin(phi)
    z = cz + r_z*np.cos(theta)
    points  = np.vstack( (x, y, z)).T
    
#     fig = plt.figure()  # Square figure
#     ax = fig.add_subplot(111, projection='3d')
#     ax.plot_surface(x, y, z, color='b')
    return np.stack([x,y,z])


# In[ ]:


#give radii and cemters in voxels
def create_ellipse_plot(cx,cy,cz, r_x, r_y, r_z):
    '''
    create sphere with center (cx, cy, cz) and radius r_z
    '''
    phi = np.linspace(0, 2*np.pi, 256).reshape(256, 1)
    theta = np.linspace(0, np.pi, 256).reshape(-1,256)

    theta, phi = np.meshgrid(theta, phi)
    x = cx + r_x*np.sin(theta)*np.cos(phi)
    y = cy + r_y*np.sin(theta)*np.sin(phi)
    z = cz + r_z*np.cos(theta)
    points  = np.vstack( (x, y, z)).T
    
    return x, y, z


# In[ ]:


#create shell of ellipses from first radius to max
#down sampling does nothing.
#can hopefully use a negative scaling factor to fill
def create_shell (cx, cy, cz, r_x, r_y, r_z, scaling, step_factor):
    spheres=[]
    try:
    	r_outer_x = round(r_x*scaling)
    	r_outer_y =  round(r_y*scaling)
    	r_outer_z = round(r_z*scaling)
    	r_diff_x=r_outer_x-r_x
    	r_diff_max=r_diff_x
    	r_diff_y=r_outer_y-r_y
    except ValueError:
    	spheres=[0] 
    	array_shell=np.array(spheres)
    	return spheres, array_shell
	 
    if r_diff_y>r_diff_max:
        r_diff_max=r_diff_y
    r_diff_z=r_outer_z-r_z
    if r_diff_z>r_diff_max:
        r_diff_max=r_diff_z

    while r_x<r_outer_y and r_y<r_outer_y and r_z<r_outer_z and r_x>2 and r_y>2 and r_z>2:
        r_x=round(r_x)
        r_y=round(r_y)
        r_z=round(r_z)
        sphere= create_ellipse(cx,cy,cz,r_x, r_y, r_z)
        xyz_list, xyz_px_list = combine_sphere_xyz_px(sphere)
        spheres.extend(xyz_px_list)
        r_x=r_x*step_factor
        r_y=r_y*step_factor
        r_z=r_z*step_factor
        
    spheres=[*set(spheres)] 
    array_shell=np.array(spheres)
    return spheres, array_shell


# In[ ]:


def measure_dims (x_px,y_px,z_px, img_stack):
#     x_px=round(node[0]/x_scale)
#     y_px=round(node[1]/y_scale)
#     z_px=round(node[2]/z_scale)
    z_slice=img_stack[z_px]
    y_slice=z_slice[y_px]
    px=y_slice[x_px]
    #first for x
    isna=np.isnan(px)
    if isna==True:
        x_diam, x_rad, y_diam, y_rad, z_diam, z_rad, x_mid, y_mid, z_mid, mean_diam = np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan
        return x_diam, x_rad, y_diam, y_rad, z_diam, z_rad, x_mid, y_mid, z_mid, mean_diam
    n=1
    x_pos_vals=[]
    while isna == False:
        i=x_px+n
        if i in range(0,len(img_stack[0][0])):
            try:
                x_pos_px=img_stack[z_px][y_px][i]
                isna=np.isnan(x_pos_px)
                x_max=i-1
                n+=1
                x_pos_vals.append(x_pos_px)
            except IndexError:
                break
        else:
            break
    isna=np.isnan(px)
    n=1
    x_neg_vals=[]
    while isna == False:
        i=x_px-n
        if i in range(0,len(img_stack[0][0])):
            try:
                x_neg_px=img_stack[z_px][y_px][i]
                isna=np.isnan(x_neg_px)
                x_min=i+1
                n+=1
                x_pos_vals.append(x_neg_px)
            except IndexError:
                break
        else:
            break
    x_diam=len(x_pos_vals)+len(x_neg_vals)
    x_rad=x_diam/2
    
    #for y
    isna=np.isnan(px)
    n=1
    y_pos_vals=[]
    while isna == False:
        i=y_px+n
        if i in range(0,len(img_stack[0][0])):
            try:
                y_pos_px=img_stack[z_px][i][x_px]
                isna=np.isnan(y_pos_px)
                y_max=i-1
                n+=1
                y_pos_vals.append(y_pos_px)
            except IndexError:
                break
        else:
            break
    isna=np.isnan(px)
    n=1
    y_neg_vals=[]
    while isna == False:
        i=y_px-n
        if i in range(0,len(img_stack[0][0])):
            try:
                y_neg_px=img_stack[z_px][i][x_px]
                isna=np.isnan(y_neg_px)
                y_min=i+1
                n+=1
                y_neg_vals.append(y_neg_px)
            except IndexError:
                break
        else:
            break
    y_diam=len(y_pos_vals)+len(y_neg_vals)
    y_rad=y_diam/2
    #now z
    isna=np.isnan(px)
    n=1
    z_pos_vals=[]
    while isna == False:
        i=z_px+n
        if i in range(len(img_stack[0])):
            try:
                z_pos_px=img_stack[i][y_px][x_px]
                isna=np.isnan(z_pos_px)
                z_max=i-1
                n+=1
                z_pos_vals.append(z_pos_px)
            except IndexError:
                break
        else:
            break
    isna=np.isnan(px)
    n=1
    z_neg_vals=[]
    while isna == False:
        i=z_px-n
        if i in range(0,len(img_stack[0])):
            try:
                z_neg_px=img_stack[i][y_px][x_px]
                isna=np.isnan(z_neg_px)
                z_min=i+1
                n+=1
                z_neg_vals.append(z_neg_px)
            except IndexError:
                break
        else:
            break
    z_diam=len(z_pos_vals)+len(z_neg_vals)
    z_rad=round(z_diam/2)
    
    x_mid=round(x_min+((x_max-x_min)/2))
    y_mid=round(y_min+((y_max-y_min)/2))
    z_mid=round(z_min+((z_max-z_min)/2))
    mean_diam=((x_diam+y_diam+z_diam)/3)
    return x_diam, x_rad, y_diam, y_rad, z_diam, z_rad, x_mid, y_mid, z_mid, mean_diam


# In[ ]:


def measure_chamber_slope (x_px,y_px,z_px, img_stack, x_scale, y_scale, z_scale):
    x_diam, x_rad, y_diam, y_rad, z_diam, z_rad, x_mid, y_mid, z_mid, mean_diam=measure_dims(x_px,y_px,z_px, img_stack)
    z_slice=img_stack[z_px]
    y_slice=z_slice[y_px]
    node_origin=y_slice[x_px]
    try:
        x_lower_bound=x_px-x_rad
        x_upper_bound=x_px+x_rad
        y_lower_bound=y_px-y_rad
        y_upper_bound=y_px+y_rad
    except TypeError or ValueError:
        x_mean_slope_cm=np.nan
        y_mean_slope_cm=np.nan
        return x_mean_slope_cm, y_mean_slope_cm
    x_gradient_lengths=[]
    y_gradient_lengths=[]
    x_gradient_lengths_neg=[]
    y_gradient_lengths_neg=[]
    x_slopes=[]
    y_slopes=[]
    #for gradient calc - takes a list down along an x axes at equidistant point from x_px
    #and calculates gradient between them.Repeats for all x_px along axis and gets mean.
    try:
    	x_px=int(x_px)
    	x_upper_bound=int(x_upper_bound)
    	y_px=int(y_px)
    	y_upper_bound=int(y_upper_bound)
    	z_px=int(z_px)
    except ValueError:
        x_mean_slope_cm=np.nan
        y_mean_slope_cm=np.nan
        return x_mean_slope_cm, y_mean_slope_cm
    i=1
    ii=1
    for px in range(x_px, x_upper_bound):
        px=px+1
        if px in range(0,len(img_stack[0][0])): #to stop sampling out of img bounds
            x_point=img_stack[z_px][y_px][px]#xpixel to sample down z from
            negative_px=x_px-i
            i+=1
            #print('px=', px, 'neg_px=', negative_px)
            x_point_neg=img_stack[z_px][y_px][negative_px]
            px_distance=(px*2)+1 #x_point-x_point_neg
            cm_distance=px_distance*x_scale
            #print(img_stack[z_px][y_px][px])
            x_gradient_points=[]
            x_gradient_points_neg=[]
            if x_point>0:
                for i in range(z_px-z_diam,z_px):
                    if i in range(0, len(img_stack[0])):#stops sampling outside z bounds
                        #print('z_px=', i, img_stack[i][y_px][px])
                        #print('z_neg_px=',img_stack[i][y_px][negative_px])
                        if img_stack[i][y_px][px]>0:
                            x_gradient_points.append(img_stack[i][y_px][px])
                            #print('zs=', img_stack[i][y_px][px])
                        else:
                            continue
                        if img_stack[i][y_px][negative_px]>0:
                            x_gradient_points_neg.append(img_stack[i][y_px][negative_px])
                            #print('zneg=',img_stack[i][y_px][negative_px])
                        else:
                            continue
                    else:
                        continue
                x_gradient_len=len(x_gradient_points)*z_scale
                x_gradient_lengths.append(x_gradient_len)
                x_gradient_neg_len=len(x_gradient_points_neg)*z_scale
                x_gradient_lengths_neg.append(x_gradient_neg_len)

                x_difference=+x_gradient_len-x_gradient_neg_len
                x_slope=x_difference/cm_distance
                x_slopes.append(x_slope)
                #print('x_slope=',x_slope)
                #print('x_len=',x_gradient_len, 'x_neg_len=',x_gradient_neg_len)
        else:
            continue
    for px2 in range(y_px, y_upper_bound):
        px2=px2+1
        if px2 in range(0,len(img_stack[0][0])):
            y_point=img_stack[z_px][px2][x_px]#xpixel to sample down z from
            negative_px2=y_px-ii
            ii+=1
            y_point_neg=img_stack[z_px][negative_px2][x_px]
            px2_distance=(px2*2)+1 #x_point-x_point_neg
            cm_distance2=px2_distance*y_scale
            y_gradient_points=[]
            y_gradient_points_neg=[]
            if y_point>0:
                for i in range(z_px-z_diam,z_px):
                    if i in range(0, len(img_stack[0])):#stops sampling outside z bounds
                        if img_stack[i][px2][x_px]>0:
                            x_gradient_points.append(img_stack[i][px2][y_px])
                        else:
                            continue
                        if img_stack[i][negative_px2][x_px]>0:
                            y_gradient_points_neg.append(img_stack[i][negative_px2][y_px])
                        else:
                            break
                    else:
                        continue
                y_gradient_len=len(y_gradient_points)*z_scale
                y_gradient_lengths.append(y_gradient_len)
                y_gradient_neg_len=len(y_gradient_points_neg)*z_scale
                y_gradient_lengths_neg.append(y_gradient_neg_len)
                y_difference=+y_gradient_len-y_gradient_neg_len
                y_slope=y_difference/cm_distance2
                y_slopes.append(y_slope)
                #print('y_slope=',y_slope)
        else:
            continue
    x_slope_array=np.array(x_slopes)
    y_slope_array=np.array(y_slopes)
    x_mean_slope_cm=np.mean(x_slope_array)
    y_mean_slope_cm=np.mean(y_slope_array)
    
    print('x_mean_slope=',x_mean_slope_cm, 'y_mean_slope=',y_mean_slope_cm)
    
    return x_mean_slope_cm, y_mean_slope_cm


# In[ ]:


#may want to only iterate over a mon 3 degree node_list. Provide the xy radius to sample over. This is scaled for differences in z scaling
def sample_intensity_sphere_px (node, img_stack, r_inner, r_outer, down_sampling, x_scale, y_scale, z_scale, z_xy_scale):
    r_inner=int(r_inner/z_scale)
    r_inner=round(r_inner)
    #print("r_inner=",r_inner)
    r_outer=int(r_outer/z_scale)
    r_outer=round(r_outer)
    #print("r_outer", r_outer)
    x_point=node[0]
    y_point=node[1]
    z_point=node[2] 
    x_px=int(round(node[0]/x_scale))
    y_px=int(round(node[1]/y_scale))
    z_px=int(round(node[2]/z_scale))
    node_px=img_stack[z_px][y_px][x_px]
    origin=int(0)
    inner_px_xyz_list, array_inner_sphere = create_shell(x_px,y_px,z_px, 0, r_inner, z_xy_scale, down_sampling)
    outer_px_xyz_list, array_outer_sphere= create_shell(x_px,y_px,z_px, r_inner, r_outer, z_xy_scale, down_sampling)
    print('check size of ellipsoids (px) inner:', len(inner_px_xyz_list),'outer:', len(array_outer_sphere))
    inner_intensities=[]
    for px in inner_px_xyz_list:
        x=px[0]
        y=px[1]
        z=px[2]
        px_intensity=img_stack[z][y][x]
        inner_intensities.append(px_intensity)

    shell_intensities=[]
    for px in outer_px_xyz_list:
        x=px[0]
        y=px[1]
        z=px[2]
        #to prevent error being thrown if pixels outside of image stack area are sampled
        if x and y in range(0,len(img_stack[0][0])) and z in range(0,len(img_stack[0])):
            px_intensity=img_stack[z][y][x]
            shell_intensities.append(px_intensity)
        else:
            continue
    inner_intensities[inner_intensities==0.]=np.nan
    #print('inner pixels sampled=',len(inner_intensities))
    inner_mean=np.nanmean(inner_intensities)
    shell_intensities[shell_intensities==0.]=np.nan
    #print('outer pixels sampled=', len(shell_intensities))
    shell_mean=np.nanmean(shell_intensities)
    node_keyed=tuple([node, inner_mean, shell_mean])
    #print('inner mean intensity=')
    return node_keyed, inner_px_xyz_list, outer_px_xyz_list



# In[ ]:


def sample_intensity_ellipse_px (x_px,y_px,z_px, img_stack, r_x, r_y, r_z, scale_factor, increment, decrement, x_scale, y_scale, z_scale, z_xy_scale):
    node_px=img_stack[z_px][y_px][x_px]
    origin=int(0)
    if np.isnan(x_px) or np.isnan(y_px) or np.isnan(z_px):
    	inner_mean=0
    	shell_mean=0
    	node_keyed=tuple([node_px, inner_mean, shell_mean])
    	print('inner mean intensity=', inner_mean, 'outer_mean_intensity=', shell_mean)
    	return node_keyed, inner_mean, shell_mean
    inner_px_xyz_list, array_inner_sphere = create_shell(x_px,y_px,z_px,r_x,r_y,r_z,scale_factor,decrement)
    outer_px_xyz_list, array_outer_sphere= create_shell(x_px,y_px,z_px,r_x,r_y,r_z,scale_factor,increment)
    outer_px_xyz_list = list(set(outer_px_xyz_list) - set(inner_px_xyz_list))
    array_outer_sphere=np.array(outer_px_xyz_list)
    isna=np.isnan(node_px)
    if len(array_outer_sphere)==0:
        inner_mean, shell_mean=0,0
        node_keyed=tuple([node_px, inner_mean, shell_mean])
        return node_keyed, inner_mean, shell_mean
    #print('check size of inner:', len(inner_px_xyz_list),'outer:', len(array_outer_sphere),'shell:')
    inner_intensities=[]
    for px in inner_px_xyz_list:
        x=px[0]
        y=px[1]
        z=px[2]
        try:
            px_intensity=img_stack[z][y][x]
            inner_intensities.append(px_intensity)
        except IndexError:
            continue

    shell_intensities=[]
    for px in outer_px_xyz_list:
        x=px[0]
        y=px[1]
        z=px[2]
        #to prevent error being thrown if pixels outside of image stack area are sampled
        if x and y in range(0,len(img_stack[0][0])) and z in range(0,len(img_stack[0])):
            try:
                px_intensity=img_stack[z][y][x]
                shell_intensities.append(px_intensity)
            except IndexError:
                continue
        else:
            continue
    inner_intensities[inner_intensities==0.]=np.nan
    #print('inner pixels sampled=',len(inner_intensities))
    inner_mean=np.nanmean(inner_intensities)
    try:
        shell_intensities[shell_intensities==0.]=np.nan
    except IndexError:
                pass
    #print('outer pixels sampled=', len(shell_intensities))
    shell_mean=np.nanmean(shell_intensities)
    node_keyed=tuple([node_px, inner_mean, shell_mean])
    print('inner mean intensity=', inner_mean, 'outer_mean_intensity=', shell_mean)
    return node_keyed, inner_mean, shell_mean


# In[ ]:


def find_maxdims (x_px,y_px,z_px, img_stack, allowance, x_scale, y_scale, z_scale):
    x_diam, x_rad, y_diam, y_rad, z_diam, z_rad, x_mid, y_mid, z_mid, mean_diam=measure_dims(x_px,y_px,z_px, img_stack)
    old_node = tuple([x_px*x_scale,y_px*y_scale,z_px*z_scale])
    best_node=old_node
    new_node=False
    if not np.isnan(x_mid) and not np.isnan(y_mid) and not np.isnan(z_mid):
        z_slice=img_stack[z_mid]
        y_slice=z_slice[y_mid]
        mid_px=y_slice[x_mid]
        isna=np.isnan(mid_px)
        new_node=True
    else:
        new_node=False
        best_node=old_node
        isna=True
    if isna==False:
        mid_node = tuple([x_mid*x_scale,y_mid*y_scale,z_mid*z_scale])
        best_node=mid_node
        best_node_px=tuple([round(x_px), round(y_px), round(z_px)])
        distance=calculate_distance(old_node, mid_node)
        new_node=True
    else:
        distance=0
        isna=True
    max_diam=mean_diam
    

    while isna == False and distance<allowance:
        new_node=img_stack[z_mid][y_mid][x_mid]
        distance=calculate_distance(old_node, best_node)
        isna=np.isnan(new_node)
        #print('mids=',x_mid, y_mid, z_mid)
        x_diam1, x_rad1, y_diam1, y_rad1, z_diam1, z_rad1, x_mid1, y_mid1, z_mid1, mean_diam1=measure_dims(x_mid,y_mid,z_mid, img_stack)
        x_diam2, x_rad2, y_diam2, y_rad2, z_diam2, z_rad2, x_mid2, y_mid2, z_mid2, mean_diam2=measure_dims(x_mid,y_mid,z_px, img_stack)
        x_diam3, x_rad3, y_diam3, y_rad3, z_diam3, z_rad3, x_mid3, y_mid3, z_mid3, mean_diam3=measure_dims(x_mid,y_px,z_mid, img_stack)
        x_diam4, x_rad4, y_diam4, y_rad4, z_diam4, z_rad4, x_mid4, y_mid4, z_mid4, mean_diam4=measure_dims(x_px,y_mid,z_mid, img_stack)
        try:
            if mean_diam1>max_diam:
                max_diam=mean_diam1
                x_mid, y_mid, z_mid = x_mid1, y_mid1, z_mid1
                x_rad, y_rad, z_rad = x_rad1, y_rad1, z_rad1
                best_node=tuple([x_mid*x_scale,y_mid*y_scale,z_mid*z_scale])
            elif mean_diam2>max_diam:
                max_diam=mean_diam2
                x_mid, y_mid, z_mid = x_mid2, y_mid2, z_mid2
                x_rad, y_rad, z_rad = x_rad2, y_rad2, z_rad2
                best_node=tuple([x_mid*x_scale,y_mid*y_scale,z_mid*z_scale])
            elif mean_diam3>max_diam:
                max_diam=mean_diam3
                x_mid, y_mid, z_mid = x_mid3, y_mid3, z_mid3
                x_rad, y_rad, z_rad = x_rad3, y_rad3, z_rad3
                best_node=tuple([x_mid*x_scale,y_mid*y_scale,z_mid*z_scale])
            elif mean_diam4>max_diam:
                max_diam=mean_diam4
                x_mid, y_mid, z_mid = x_mid4, y_mid4, z_mid4
                x_rad, y_rad, z_rad = x_rad4, y_rad4, z_rad4
                best_node=tuple([x_mid*x_scale,y_mid*y_scale,z_mid*z_scale])
        except TypeError or ValueError:
            break
        else:
            break
    print('Finding center of chamber - old_node=',old_node, 'max_dim_node=',best_node)
    if new_node==True:
        best_node_px=tuple([round(x_mid), round(y_mid), round(z_mid)])
        x_rad, y_rad, z_rad=round(x_rad), round(y_rad), round(z_rad)
    else:
        best_node_px=tuple([round(x_px), round(y_px), round(z_px)])
    
    return best_node, best_node_px, max_diam, x_rad, y_rad, z_rad
            



# In[ ]:


def define_nest_entrance(D1NE, edge_list, threshold_from_max):
    if not isinstance(D1NE, bool):
        raise TypeError("D1NE must be a boolean value (True or False) if you want nest entrances to only be D1")
    terminal_nodes, internal_nodes, node_list=nodes_by_degree(edge_list)
    if D1NE==True:
        node_list=terminal_nodes
    zmax=0
    terminal_nodes_defined=[]
    for node in node_list:
        z=node[2]
        if z>zmax:
            zmax=z
        else:
            continue
    not_possible_nes=[]       
    for edge in edge_list:
        node1=edge[0]
        node2=edge[1]
        node1_z=node1[2]
        node2_z=node2[2]
        if node1_z>node2_z:
            not_possible_nes.append(node2)
        elif node2_z>node1_z:
            not_possible_nes.append(node1)
        elif node2_z==node1_z:
            not_possible_nes.append(node1)
            not_possible_nes.append(node2)
    possible_nes=list(set(node_list)-set(not_possible_nes))        
    z_threshold=zmax-threshold_from_max
    for node in possible_nes:
        z=node[2]
        if z>z_threshold:
            terminal_nodes_defined.append(node)
    return terminal_nodes_defined, zmax


# In[ ]:


def assign_chambers (D1NE, edge_list, img_stack, NE_threshold, slope_max_thresh, intensity_thresh, allowance, scale_factor, increment, decrement, x_cm=None, y_cm=None, z_cm=None, voxel=None):
    if not isinstance(D1NE, bool):
        raise TypeError("D1NE must be a boolean value (True or False) if you want nest entrances to only be D1")
    if voxel:
        x_scale = y_scale = z_scale = voxel
        z_xy_scaling = 1
    else:
    	x_scale, y_scale, z_scale, z_xy_scaling = convert_cm_to_px(img_stack, x_cm, y_cm, z_cm)
    terminal_nodes, internal_nodes, node_list=nodes_by_degree(edge_list)
    nest_entrances=define_nest_entrance(D1NE, edge_list, NE_threshold)
    #print('nest entrances=',nest_entrances)
    internal_nodes=list(set(internal_nodes)-set(nest_entrances))
    tunnel_ends=list(set(terminal_nodes)-set(nest_entrances))
    chambers=[]
    junctions=[]
    ellipsoid_vols=[]
    for node in internal_nodes:
        or_x_px=round(node[0]/x_scale)
        or_y_px=round(node[1]/y_scale)
        or_z_px=round(node[2]/z_scale)
        best_node, best_node_px, mean_diam, x_rad, y_rad, z_rad=find_maxdims(or_x_px,or_y_px,or_z_px, img_stack, allowance, x_scale, y_scale, z_scale)
        #print('best node=', best_node_px)
        x_px, y_px, z_px=best_node_px[0], best_node_px[1], best_node_px[2]
        if np.isnan(x_px) or np.isnan(y_px) or np.isnan(z_px):
            junctions.append(node)
            continue
        new_node=img_stack[z_px][y_px][x_px]
        isna=np.isnan(new_node)
        if isna==True:
            #print(isna)
            x_px, y_px,z_px=or_x_px, or_y_px,or_z_px
            best_node=node
        ##is the potential chamber flat? measure slope
        #print(x_px, y_px, z_px, "pxls")
        x_mean_slope_cm, y_mean_slope_cm=measure_chamber_slope(x_px,y_px,z_px, img_stack, x_scale, y_scale, z_scale)
        print('slopes=', x_mean_slope_cm, y_mean_slope_cm)
        ##is the potential chamber enlarge? sample px intensity
        try:
            node_keyed, inner_mean, outer_mean = sample_intensity_ellipse_px(x_px,y_px,z_px, img_stack, x_rad, y_rad, z_rad, scale_factor, increment, decrement, x_scale, y_scale, z_scale, z_xy_scaling)
        except TypeError or ValueError:
            junctions.append(node)
            continue
        ##test slope and intenstisty
        if inner_mean==0 or outer_mean==0:
            intensity_diff=0
        else:
            intensity_diff=inner_mean/outer_mean
        x_mean_slope_cm, y_mean_slope_cm=abs(x_mean_slope_cm), abs(y_mean_slope_cm)
        mean_slope=np.nanmean(x_mean_slope_cm, y_mean_slope_cm)
        print('intensity diff=', intensity_diff)
        if mean_slope<slope_max_thresh and intensity_diff>intensity_thresh:
            ellipsoid_vol= 1.33 * math.pi * x_rad * y_rad * z_rad
            ellipsoid_vols.append(ellipsoid_vol)
            chambers.append(node)
            print('chamber=', node)
        else:
            junctions.append(node)
            #print('node=', node)
    return nest_entrances, tunnel_ends, chambers, junctions, ellipsoid_vols



#THIS TAKES TUN DIAM AS CRITERIA ALSO
def assign_chambers_2 (D1NE, MEAN_TUN_DIAMS, ch_tun_thrsh, edge_list, img_stack, NE_threshold, slope_max_thresh, slope_max_thresh2, intensity_thresh, intensity_thresh2, slope_max_thresh3, intensity_thresh3, slope_max_thresh4, allowance, scale_factor, increment, decrement, x_cm=None, y_cm=None, z_cm=None, voxel=None):
    if not isinstance(D1NE, bool):
        raise TypeError("D1NE must be a boolean value (True or False) if you want nest entrances to only be D1")
    if voxel:
        x_scale = y_scale = z_scale = voxel
        z_xy_scaling = 1
    else:
    	x_scale, y_scale, z_scale, z_xy_scaling = convert_cm_to_px(img_stack, x_cm, y_cm, z_cm)
    terminal_nodes, internal_nodes, node_list=nodes_by_degree(edge_list)
    nest_entrances=define_nest_entrance(D1NE, edge_list, NE_threshold)
    #print('nest entrances=',nest_entrances)
    tunnel_ends=list(set(terminal_nodes)-set(nest_entrances))
    
    internal_nodes=list(set(internal_nodes)-set(nest_entrances))
    
    chambers=[]
    junctions=[]
    ellipsoid_vols=[]
    for node in internal_nodes:
        #connected_ne=False
        #for edge in edge_list:
            #if node==edge[0] and edge[1] in nest_entrances:
             #   connected_ne=True
            #if node==edge[1] and edge[0] in nest_entrances:
             #   connected_ne=True
        #if connected_ne==True:   
            #junctions.append(node)
           # continue

        or_x_px=round(node[0]/x_scale)
        or_y_px=round(node[1]/y_scale)
        or_z_px=round(node[2]/z_scale)
        best_node, best_node_px, mean_diam, x_rad, y_rad, z_rad=find_maxdims(or_x_px,or_y_px,or_z_px, img_stack, allowance, x_scale, y_scale, z_scale)
        #print('best node=', best_node_px)
        x_px, y_px, z_px=best_node_px[0], best_node_px[1], best_node_px[2]
        if np.isnan(x_px) or np.isnan(y_px) or np.isnan(z_px):
            junctions.append(node)
            continue
        new_node=img_stack[z_px][y_px][x_px]
        isna=np.isnan(new_node)
        if isna==True:
            #print(isna)
            x_px, y_px,z_px=or_x_px, or_y_px,or_z_px
            best_node=node
        ##is the potential chamber flat? measure slope
        #print(x_px, y_px, z_px, "pxls")
        x_mean_slope_cm, y_mean_slope_cm=measure_chamber_slope(x_px,y_px,z_px, img_stack, x_scale, y_scale, z_scale)
        #print('slopes=', x_mean_slope_cm, y_mean_slope_cm)
        
        ##is the potential chamber enlarge? sample px intensity
        try:
            node_keyed, inner_mean, outer_mean = sample_intensity_ellipse_px(x_px,y_px,z_px, img_stack, x_rad, y_rad, z_rad, scale_factor, increment, decrement, x_scale, y_scale, z_scale, z_xy_scaling)
        except TypeError or ValueError:
            junctions.append(node)
            continue
        ##test slope and intenstisty
        if inner_mean==0 or outer_mean==0:
            intensity_diff=0
        else:
            intensity_diff=inner_mean/outer_mean
        x_mean_slope_cm, y_mean_slope_cm=abs(x_mean_slope_cm), abs(y_mean_slope_cm)
        if x_mean_slope_cm==0 and y_mean_slope_cm==0:
            mean_slope=0
        else:
            mean_slope=(x_mean_slope_cm+y_mean_slope_cm)/2
        #print('intensity diff=', intensity_diff)
        if mean_slope<slope_max_thresh and intensity_diff>intensity_thresh or mean_slope<slope_max_thresh2 and (mean_diam*x_scale)>(MEAN_TUN_DIAMS*ch_tun_thrsh) or mean_slope<slope_max_thresh3 and intensity_diff>intensity_thresh2 or mean_slope<slope_max_thresh4 and intensity_diff>intensity_thresh3:
            ellipsoid_vol= 1.33 * math.pi * x_rad * y_rad * z_rad
            ellipsoid_vols.append(ellipsoid_vol)
            chambers.append(node)
            print('CHAMBER=', node)
        else:
            junctions.append(node)
            #print('node=', node)
    return nest_entrances, tunnel_ends, chambers, junctions, ellipsoid_vols
# In[ ]:


def make_key(dlist_edge):
    conn_list=[]
    weight_list=[]
    for e in dlist_edge:
        con1 = tuple(e[0][0])
        key1= e[0][1]
        con2 = tuple(e[1][0])
        key2=e[1][1]
        d=e[2][0]
#             d =int(d)
            #d=d
        connection = [con1, key1, con2, key2, d]
        conn_list.append(connection)
        weight_list.append(d)
    return conn_list, weight_list


# In[ ]:


def temporal_linkage (ALL_NE, ALL_TE,ALL_CH, ALL_JU, list_of_edge_lists, threshold):
    keyed_list_of_list_nodes=[]
    l=1
    all_keyed_nodes=[]
    used_nodes=[]
    keyed_edge_list_dist_all=[]
    array_edge_list=remove_length(list_of_edge_lists[0])
    list_node1, array_node = nodez(array_edge_list)
    node_list_keyed1=[]
    for node in list_node1:
        keyed_node=tuple([node, l])
        l+=1
        node_list_keyed1.append(keyed_node)
        all_keyed_nodes.append(keyed_node)
    keyed_list_of_list_nodes.append(node_list_keyed1)
    for edge_list in range(1, len(list_of_edge_lists)):
        node_list_keyed1=[]
        node_list_keyed2=[]
        used_nodes=[]
        used_nodes2=[]
        array_edge_list=remove_length(list_of_edge_lists[edge_list])
        #print(len(array_edge_list))
        list_node, array_node = nodez(array_edge_list)
        list_node_CH=ALL_CH[edge_list]
        prev_list_node_CH=ALL_CH[edge_list-1]
        list_node_JU=ALL_JU[edge_list]
        prev_list_node_JU=ALL_JU[edge_list-1]
        list_node_NE=ALL_NE[edge_list]
        prev_list_node_NE=ALL_NE[edge_list-1]
        list_node_TE=ALL_TE[edge_list]
        prev_list_node_TE=ALL_TE[edge_list-1]
        
        list_node2=list_node
        list_node3=list_node
        for node in list_node:
            issamelocation=False
            previous_list_node=keyed_list_of_list_nodes[edge_list-1]
            for keyed_node in previous_list_node:
                unkeyed_node=keyed_node[0]
                issameloc=calculate_distance(unkeyed_node,node)
                if issameloc==0:
                    #print('sameloc', keyed_node)
                    keyed_node1=keyed_node
                    node_list_keyed1.append(keyed_node1)
                    issamelocation=True
            if issamelocation==True:
                continue
            #if node is a chamber, it could only have been a junction or chamber beforehand
            if node in list_node_CH:
                list_CH=[]
                for keyed_node in previous_list_node:
                    unkeyed_node=keyed_node[0]
                    if unkeyed_node in prev_list_node_CH or prev_list_node_JU:
                        list_CH.append(keyed_node)
                        #print('cham', keyed_node)
                previous_list_node=list_CH
            #if it was a tunnel enf, it could only have been a tunnel end pr
            elif node in list_node_TE:
                list_TE=[]
                for keyed_node in previous_list_node:
                    unkeyed_node=keyed_node[0]
                    if unkeyed_node in prev_list_node_TE:
                        list_TE.append(keyed_node)
                previous_list_node=list_TE
            #JUNCTION COULD NOT BE A CHAMBER BEFOREHAND
            elif node in list_node_JU:
                list_JU=[]
                for keyed_node in previous_list_node:
                    unkeyed_node=keyed_node[0]
                    if unkeyed_node not in prev_list_node_CH:
                        list_JU.append(keyed_node)
                previous_list_node=list_JU
                
                
            further_node=False
            even_further_node=False
            min_distance=threshold
            
            for keyed_node in previous_list_node:
                if keyed_node not in used_nodes:
                    distance = calculate_distance(node, keyed_node[0])
                    if distance<min_distance:
                        assigned_node=node
                        min_distance=distance
                        closest_node=keyed_node
                        min_distance=distance
                else:
                    continue
        
            if min_distance<threshold:
                key=closest_node[1]
                #print(edge_list, key)
                node2=closest_node[0]
                used_nodes.append(closest_node)
                keyed_node1=tuple([node,key])
                keyed_node2=tuple([node2,key])
                 ##refinement_step - could be a closer node
                for node3 in list_node2:
                    dist=calculate_distance(node3, node)
                    if dist!=0:
                        dist2=calculate_distance(node3, node2)
                        if dist2<min_distance and node3 not in used_nodes2:
                            further_node=True
                            #print(key, 'fn', node3, node2)
                
                if further_node==True:
                    for node4 in list_node3:
                        dist3=calculate_distance(node3, node4)
                        if dist3!=0:
                            if dist3<min_distance and node4 not in used_nodes2:
                                even_further_node=True
                    if even_further_node==True:
                        node_list_keyed1.append(keyed_node1)
                        used_nodes2.append(node)
                        all_keyed_nodes.append(keyed_node1)
                        node_list_keyed2.append(keyed_node2)
                        all_keyed_nodes.append(keyed_node2)
                        continue
                    else:
                        keyed_node1=tuple([node,l])
                        node_list_keyed1.append(keyed_node1)
                        all_keyed_nodes.append(keyed_node1)
                        l+=1
                        continue
               
                if keyed_node1 not in all_keyed_nodes:
                    #print(key, 'fn', node3, node2)
                    node_list_keyed1.append(keyed_node1)
                    used_nodes2.append(node)
                    all_keyed_nodes.append(keyed_node1)
                    node_list_keyed2.append(keyed_node2)
                    all_keyed_nodes.append(keyed_node2)
                else:
                    keyed_node1=tuple([node,l])
                    node_list_keyed1.append(keyed_node1)
                    all_keyed_nodes.append(keyed_node1)
                    used_nodes2.append(node)
                    l+=1
                if keyed_node2 not in all_keyed_nodes:
                    node_list_keyed2.append(keyed_node2)
                    all_keyed_nodes.append(keyed_node2)
            else:
                keyed_node1=tuple([node,l])
                node_list_keyed1.append(keyed_node1)
                all_keyed_nodes.append(keyed_node1)
                used_nodes2.append(node)
                l+=1
        keyed_list_of_list_nodes.append(node_list_keyed1)
    all_keyed_nodes = [x for x in all_keyed_nodes if not isinstance(x, int)]
    #print(all_keyed_nodes)
    #now to pair indices with edges and length
    ##probelm
    for edge_list in range(0, len(list_of_edge_lists)):
        keyed_edge_list_dist=[]
        for edge in list_of_edge_lists[edge_list]:
            node1=edge[0]
            node2=edge[1]
            dist=edge[2]
            for keyed_node in keyed_list_of_list_nodes[edge_list]: 
                isint=isinstance(keyed_node, int)
                if isint==True:
                    continue
                distance1=calculate_distance(node1, keyed_node[0])
                if distance1==0:
                    key_node1=keyed_node[1]
                    keyed_dist_node1=tuple([node1, key_node1])
                    all_keyed_nodes.append(key_node1)
                distance2=calculate_distance(node2, keyed_node[0])
                if distance2==0:
                    key_node2=keyed_node[1]
                    keyed_dist_node2=tuple([node2, key_node2])
                    all_keyed_nodes.append(key_node2)

            keyed_dist_edge=tuple([keyed_dist_node1, keyed_dist_node2, dist])
            keyed_edge_list_dist.append(keyed_dist_edge)
        keyed_edge_list_dist_all.append(keyed_edge_list_dist)
    
    return keyed_list_of_list_nodes, keyed_edge_list_dist_all


# In[ ]:


def calculate_angle (point1, point2, point3):
    distA=calculate_distance(point1, point2)
    distB=calculate_distance(point1, point3)
    #dist3es not represent an edge, and is imaginary
    distC=calculate_distance(point2, point3)
    cosC=((distA**2)+(distB**2)-(distC**2))/(2*(distA*distB))
    angle_rad=math.acos(cosC)
    angle_deg=angle_rad*180/math.pi
    print("angle between edges rad, deg=", angle_rad, angle_deg)
    ##now for vertical angle
    vector1=np.array(point3)-np.array(point1)
    vector2=np.array(point2)-np.array(point1)
    vector1=[abs(x) for x in vector1]
    vector2=[abs(x) for x in vector2]
    v1_x, v1_y, v1_z=vector1[0], vector1[1], vector1[2]
    v2_x, v2_y, v2_z=vector2[0], vector2[1], vector2[2]
    dot_product = (v1_x * v2_x) + (v1_y * v2_y) + (v1_z * v2_z)
    magnitude_v1 = np.linalg.norm(vector1)
    magnitude_v2 = np.linalg.norm(vector2)
#     cos_theta = dot_product / (magnitude_v1 * magnitude_v2)
#     angle_rad = math.acos(cos_theta)
#     print("angle between edges rad, deg=", angle_rad, angle_deg)
    #now vertical angle 
    x1, y1, z1 =point1[0], point1[1], point1[2]
    x2, y2, z2 =point2[0], point2[1], point2[2]
    x3, y3, z3 =point3[0], point3[1], point3[2]
    #this is in the same xy plane as the second point, but has point3s z, so just for vertical
    point4=tuple([x2,y2,z3])
    distD=calculate_distance(point1, point4)
    distE=calculate_distance(point2, point4)
    cosE=((distA**2)+(distD**2)-(distE**2))/(2*(distA*distD))
    vertical_angle = math.acos(cosE)
    #print("vertical angle=", vertical_angle)
    #now vertical
    distF=calculate_distance(point3, point4)
    cosF=((distD**2)+(distB**2)-(distF**2))/(2*(distD*distB))
    horizontal_angle = math.acos(cosF)
    #print("horizontal angle=", horizontal_angle)
    return angle_rad, vertical_angle, horizontal_angle


# In[ ]:
#calculate angles
def edge_angles (list_edge):
    list_edge1=list_edge
    list_3D_angles=[]
    list_vertical_angles=[]
    list_horizontal_angles=[]
    for edge in list_edge:
        node1=edge[0]
        node2=edge[1]
        for edge1 in list_edge1:
            node3=edge1[0]
            node4=edge1[1]
            dist1=calculate_distance(node1, node3)
            dist2=calculate_distance(node2, node4)
            if dist1==0 and dist2==0:
                continue
            else:
                if dist1==0:
                    D3a, va,ha=calculate_angle(node1, node2, node4)
                    list_3D_angles.append(D3a)
                    list_vertical_angles.append(va)
                    list_horizontal_angles.append(ha)
                elif dist2==0:
                    D3a, va,ha=calculate_angle(node2, node1, node3)
                    list_3D_angles.append(D3a)
                    list_vertical_angles.append(va)
                    list_horizontal_angles.append(ha)
    list_3D_angles=np.array([*set(list_3D_angles)])
    list_vertical_angles=np.array([*set(list_vertical_angles)])
    list_horizontal_angles=np.array([*set(list_horizontal_angles)])

    return list_3D_angles, list_vertical_angles, list_horizontal_angles


def chamber_depth (node_list, chamber_list):
    zmax=0
    depth_list=[]
    for node in node_list:
        z=node[2]
        if z>zmax:
            zmax=z
        else:
            continue
    for node2 in chamber_list:
        cham=node2[0]
        z=cham[2]
        depth=zmax-z
        depth_list.append(depth)
    depth_list=np.array(depth_list)
    mean_depth=np.mean(depth_list)
    depth_dev=np.std(depth_list)
    return depth_list, mean_depth, depth_dev
    
def calculate_depth (node_list, chamber_list):
    zmax=0
    zmin=100
    depth_list=[]
    for node in node_list:
        z=node[2]
        if z>zmax:
            zmax=z
        if z<zmin:
            zmin=z
        else:
            continue
        
    for node in node_list:
        z=node[2]
        depth=zmax-z
        depth_list.append(depth)
    max_depth=zmax-zmin
    depth_list=np.array(depth_list)
    mean_depth=np.mean(depth_list)
    depth_dev=np.std(depth_list)
    return depth_list, mean_depth, depth_dev, max_depth

def tunnel_dims (edge, img_stack, x_cm=None, y_cm=None, z_cm=None, voxel=None):
    n=500  
        #for edge in edge_list:
    if voxel:
        x_scale = y_scale = z_scale = voxel
        z_xy_scaling = 1
    else:
    	x_scale, y_scale, z_scale, z_xy_scaling = convert_cm_to_px(img_stack, x_cm, y_cm, z_cm)

    node1=np.array(edge[0])
    node2=np.array(edge[1])
    node1_px=np.array([node1[0]/x_scale, node1[1]/y_scale, node1[0]/z_scale])
    node2_px=np.array([node2[0]/x_scale, node2[1]/y_scale, node2[0]/z_scale])
    vector=node1-node2
    vector_px=node1_px-node2_px
    vector=vector_px
    midpoint=midpoint_node(node1,node2)

    mid_x_px=round(midpoint[0]/x_scale)
    mid_y_px=round(midpoint[1]/y_scale)
    mid_z_px=round(midpoint[2]/z_scale)
    midpoint_px=[mid_x_px, mid_y_px, mid_z_px]
    mp=midpoint_px
    isna=np.isnan(img_stack[mp[2], mp[1], mp[0]])
    if isna==True:
        dim1=np.nan
        dim2=np.nan
        mean_dim=np.nan
        area=np.nan
        return dim1, dim2, mean_dim, midpoint, area
    else:
        # Calculate the first orthogonal vector
        v1 = np.array([-vector[1], vector[0], 0])
        v1 = v1 / np.linalg.norm(v1)
        t = np.linspace(-n, n, 2*n+1)
        index_v1 = midpoint_px + t[:, np.newaxis] * v1[np.newaxis, :]
        index_v1 = np.round(index_v1).astype(int)
        #generate the second orthogonal vector
        v2 = np.array([-vector[2], 0, vector[0]])
        v2 = v2 / np.linalg.norm(v2)
        t = np.linspace(-n, n, 2*n+1)
        index_v2 = midpoint_px + t[:, np.newaxis] * v2[np.newaxis, :]
        index_v2 = np.round(index_v2).astype(int)

        pixels_v1 = []
        ind=0
        for i in index_v1:
            if i[0] >= 0 and i[0] < img_stack.shape[2] and i[1] >= 0 and i[1] < img_stack.shape[1] and i[2] >= 0 and i[2] < img_stack.shape[0]:
                pixels_v1.append(img_stack[i[2], i[1], i[0]])
                ind+=1
                try:
                    if tuple(i)==tuple(midpoint_px):
                        ind1=ind
                except TypeError or ValueError:
                    continue

        pixels_v2 = []
        ind=0
        for i in index_v2:
            if i[0] >= 0 and i[0] < img_stack.shape[2] and i[1] >= 0 and i[1] < img_stack.shape[1] and i[2] >= 0 and i[2] < img_stack.shape[0]:
                pixels_v2.append(img_stack[i[2], i[1], i[0]])
                ind+=1
                try:
                    if tuple(i)==tuple(midpoint_px):
                        ind2=ind
                except TypeError or ValueError:
                    continue

        all_pixels_v1=[]
        isna=False
        try:
            for i in range(ind1,len(pixels_v1)):
                try:
                    while not isna:
                        isna=np.isnan(pixels_v1[i])
                        all_pixels_v1.append(pixels_v1[i])
                        i+=1
                except IndexError:
                    break
            isna=False
            for i in range(ind1,0,-1):
                try:
                    while isna==False:
                        isna=np.isnan(pixels_v1[i])
                        all_pixels_v1.append(pixels_v1[i])
                        i-=1
                except IndexError:
                    break
            all_pixels_v2=[]
            isna=False
            for i in range(ind2,len(pixels_v2)):
                try:
                    while isna==False:
                        isna=np.isnan(pixels_v2[i])
                        all_pixels_v2.append(pixels_v2[i])
                        i+=1
                except IndexError:
                    break
            isna=False
            for i in range(ind2,0,-1):
                try:
                    while isna==False:
                        isna=np.isnan(pixels_v2[i])
                        all_pixels_v2.append(pixels_v1[i])
                        i-=1
                except IndexError:
                    break
        except UnboundLocalError:
            dim1=np.nan
            dim2=np.nan
            mean_dim=np.nan
            area=np.nan
            return dim1, dim2, mean_dim, midpoint, area	     
        
        dim1_px=len(all_pixels_v1)
        dim2_px=len(all_pixels_v2)
        #need to make a combined scale as tunnels may be diagonal
        combined_scale=(x_scale+y_scale)/2
        dim1=dim1_px*combined_scale
        dim2=dim2_px*combined_scale
        mean_dim=(dim1+dim2)/2
        #calculate cross-sectional area
        r1=dim1/2
        r2=dim2/2
        cs_area=math.pi * (r1**2 + r2**2 + r1 * r2) / 4
    return dim1, dim2, mean_dim, midpoint, cs_area

def measure_dims_xy (node, img_stack, x_cm=None, y_cm=None, z_cm=None, voxel=None):
    
    if voxel:
        x_scale = y_scale = z_scale = voxel
        z_xy_scaling = 1
    else:
    	x_scale, y_scale, z_scale, z_xy_scaling = convert_cm_to_px(img_stack, x_cm, y_cm, z_cm)
    x_px=round(node[0]/x_scale)
    y_px=round(node[1]/y_scale)
    z_px=round(node[2]/z_scale)
    z_slice=img_stack[z_px]
    y_slice=z_slice[y_px]
    px=y_slice[x_px]
    #first for x
    isna=np.isnan(px)
    if isna==True:
        x_diam,y_diam,mean_diam = np.nan,np.nan,np.nan
        return x_diam, y_diam, mean_diam
    n=1
    x_pos_vals=[]
    while isna == False:
        i=x_px+n
        if i in range(0,len(img_stack[0][0])):
            try:
                x_pos_px=img_stack[z_px][y_px][i]
                isna=np.isnan(x_pos_px)
                x_max=i-1
                n+=1
                x_pos_vals.append(x_pos_px)
            except IndexError:
                break
        else:
            break
    isna=np.isnan(px)
    n=1
    x_neg_vals=[]
    while isna == False:
        i=x_px-n
        if i in range(0,len(img_stack[0][0])):
            try:
                x_neg_px=img_stack[z_px][y_px][i]
                isna=np.isnan(x_neg_px)
                x_min=i+1
                n+=1
                x_pos_vals.append(x_neg_px)
            except IndexError:
                break
        else:
            break
    x_diam=len(x_pos_vals)+len(x_neg_vals)
    x_diam=x_diam*x_scale
    x_rad=x_diam/2
    
    #for y
    isna=np.isnan(px)
    n=1
    y_pos_vals=[]
    while isna == False:
        i=y_px+n
        if i in range(0,len(img_stack[0][0])):
            try:
                y_pos_px=img_stack[z_px][i][x_px]
                isna=np.isnan(y_pos_px)
                y_max=i-1
                n+=1
                y_pos_vals.append(y_pos_px)
            except IndexError:
                break
        else:
            break
    isna=np.isnan(px)
    n=1
    y_neg_vals=[]
    while isna == False:
        i=y_px-n
        if i in range(0,len(img_stack[0][0])):
            try:
                y_neg_px=img_stack[z_px][i][x_px]
                isna=np.isnan(y_neg_px)
                y_min=i+1
                n+=1
                y_neg_vals.append(y_neg_px)
            except IndexError:
                break
        else:
            break
    y_diam=len(y_pos_vals)+len(y_neg_vals)
    y_diam=y_diam*y_scale
    y_rad=y_diam/2
    #now z
    isna=np.isnan(px)
    n=1

    x_mid=round(x_min+((x_max-x_min)/2))
    y_mid=round(y_min+((y_max-y_min)/2))
    mean_diam=((x_diam+y_diam)/2)
    return x_diam, y_diam, mean_diam

def count_boxes(image, box_size):
    # Divide the image into boxes of size `box_size`
    box_rows = int(np.ceil(image.shape[2] / box_size))
    box_cols = int(np.ceil(image.shape[1] / box_size))
    box_depths = int(np.ceil(image.shape[1] / box_size))
    boxes = np.zeros((box_rows, box_cols, box_depths), dtype=bool)
    for i in range(box_rows):
        for j in range(box_cols):
            for k in range(box_depths):
                row_start = int(i * box_size)
                row_end = int((i + 1) * box_size)
                col_start = int(j * box_size)
                col_end = int((j + 1) * box_size)
                depth_start = int(k * box_size)
                depth_end = int((k + 1) * box_size)
                box = image[row_start:row_end, col_start:col_end, depth_start:depth_end]
                boxes[i, j, k] = np.sum(box) > 0
    return np.sum(boxes)
  
    
def fractal_dimension(image, box_size):
    # Binarize the image
    image = image > 0
    
    # Compute N(r) for a range of box sizes
    log_r = []
    log_N = []
    for i in range(1, 20):
        r = box_size * 2**(-i)
        box_count = count_boxes(image, r)
        log_r.append(np.log(r))
        log_N.append(np.log(box_count))
        # Fit a line to the log-log plot
    fit = np.polyfit(log_r, log_N, 1)
    return fit[0]
   
 #this looks confusing. It calculates the shortest paths between all pairs of nodes in a g.. It does everything twice, using edge number and dijkstra distances. It returns the diamete in edges first, then the associated path, then for dijkstra it gives
 #edge number between each pair then summed edge wights  then for network measure it does the  same.   
def calculate_diameter(G):
    #Gs = sorted(nx.connected_components(G), key=len, reverse=True)
    #Gmax = G.subgraph(Gs[0])
    
    # Calculate spatial diameter
    max_spatial_diameter = 0
    spatial_edges = []
    spatial_lengths=[]
    for node1, node2 in itertools.combinations(G.nodes, 2):
        
        path = nx.dijkstra_path(G, node1,node2, weight='weight')
        edge_weights = [G.get_edge_data(path[i], path[i+1])['weight'] for i in range(len(path) - 1)]
        edges = [(path[i], path[i+1]) for i in range(len(path)-1)]
        max_spatial_path=np.sum(edge_weights)
        #shortest_paths = list(shortest_paths)
        path_length=len(edges)
        spatial_edges.append(path_length)
        spatial_lengths.append(max_spatial_path)
        if max_spatial_path > max_spatial_diameter:
            max_spatial_diameter = max_spatial_path
            max_spatial_edges=path_length


    # Calculate edge diameter
    edge_paths=[]
    max_edge_diameter = 0
    edge_lengths=[]
    for node1, node2 in itertools.combinations(G.nodes, 2):
        path = nx.shortest_path(G, node1, node2)
        edge_weights = [G.get_edge_data(path[i], path[i+1])['weight'] for i in range(len(path) - 1)]
        edges = [(path[i], path[i+1]) for i in range(len(path)-1)]
        path_length=len(edges)
        edge_paths.append(path_length)
        max_edge_path=np.sum(edge_weights)
        edge_lengths.append(max_edge_path)
        if path_length > max_edge_diameter:
            max_edge_diameter = path_length
            max_edge_spatial_path=np.sum(edge_weights)
    
    return max_spatial_edges, max_edge_diameter,  max_spatial_diameter, max_edge_spatial_path, spatial_edges, spatial_lengths, edge_paths, edge_lengths
    
def cs_area(dim1, dim2):
	r1=dim1/2
	r2=dim2/2
	area=math.pi * (r1**2 + r2**2 + r1 * r2) / 4
	return area

def key_width (edge, tun_dims):
    mid=midpoint_node(edge[0], edge[1])
    i=0
    for dim in tun_dims:
        width=dim[1]
        dist=calculate_distance(mid,dim[0])
        if i==0:
            min_dist=dist
            keyed_width=width
        elif dist<min_dist:
            min_dist=dist
            keyed_width=width
        i+=1
    return keyed_width, min_dist
    
    
def chamber_vols (ch_list, img_stack, scale_factor, increment, decrement,allowance, x_cm=None, y_cm=None, z_cm=None, voxel=None):
    if voxel:
        x_scale = y_scale = z_scale = voxel
        z_xy_scaling = 1
    else:
    	x_scale, y_scale, z_scale, z_xy_scaling = convert_cm_to_px(img_stack, x_cm, y_cm, z_cm)
    chambers=[]
    junctions=[]
    ellipsoid_vols=[]
    for node in ch_list:
        xyz=node[0]
        or_x_px=round(xyz[0]/x_scale)
        or_y_px=round(xyz[1]/y_scale)
        or_z_px=round(xyz[2]/z_scale)
        best_node, best_node_px, mean_diam, x_rad, y_rad, z_rad=find_maxdims(or_x_px,or_y_px,or_z_px, img_stack, allowance, x_scale, y_scale, z_scale)
        ellipsoid_vol= 1.33 * math.pi * (x_rad*x_scale) * (y_rad*y_scale) * (z_rad*z_scale)
        print('Chamber vol=', ellipsoid_vol)
        ellipsoid_vols.append(ellipsoid_vol)
    return ellipsoid_vols
    
    ##TO BE RUN AFTER FIRST TEMPORAL LINKAGE, TO RETAIN CHAMS
    ##TO BE RUN AFTER FIRST TEMPORAL LINKAGE, TO RETAIN CHAMS
def temporal_linkage_2 (keyed_list_of_list_nodes, ALL_CH, ALL_JU, ALL_NE, ALL_TE, list_of_edge_lists, THRSH, CH_THRSH):
    LL_CHAMS=[]
    LL_CHAM_NUMS=[]
    LL_JU=[]
    LL_JU_NUMS=[]
    LL_NE=[]
    LL_NE_NUMS=[]
    LL_TE=[]
    LL_TE_NUMS=[]
    keyed_list_of_list_nodes2=[]
    L=0
    node_len=len(keyed_list_of_list_nodes)-1
    for node in range(0,len(keyed_list_of_list_nodes[node_len])):
        l=keyed_list_of_list_nodes[node_len][node][1]
        if l>L:
            L=l
        else:
            continue
    L+=1
    
    for i in range(0,len(keyed_list_of_list_nodes)):
        edge_list=list_of_edge_lists[i]
        
        NE=ALL_NE[i]
        TE=ALL_TE[i]
        CH=ALL_CH[i]
        JU=ALL_JU[i]

        chambers=[]
        cham_nums=[]
        jus=[]
        ju_nums=[]
        nes=[]
        ne_nums=[]
        tes=[]
        te_nums=[]

        if i == 0:
            for node1 in keyed_list_of_list_nodes[i]:
                number=node1[1]
                node=node1[0]
                paired=0
                while paired==0:
                    for ch in CH:
                        dist=calculate_distance(node, ch)
                        if dist==0:
                            chambers.append(node1)
                            cham_nums.append(number)
                        
                            paired=1
                    for ne in NE:
                        dist=calculate_distance(node, ne)
                        if dist==0:
                            nes.append(node1)
                            ne_nums.append(number)
                            paired=1
                    for ju in JU:
                        dist=calculate_distance(node, ju)
                        if dist==0:
                            jus.append(node1)
                            ju_nums.append(number)
                            paired=1
                    for te in TE:
                        dist=calculate_distance(node, te)
                        if dist==0:
                            tes.append(node1)
                            te_nums.append(number)
                            paired=1
            LL_CHAMS.append(chambers)
            LL_CHAM_NUMS.append(cham_nums)
            LL_JU.append(jus)
            LL_JU_NUMS.append(ju_nums)
            LL_NE.append(nes)
            LL_NE_NUMS.append(ne_nums)
            LL_TE.append(tes)
            LL_TE_NUMS.append(te_nums)
            ALL_NODES=[]
            ALL_NODES.extend(chambers)
            ALL_NODES.extend(jus)
            ALL_NODES.extend(nes)
            ALL_NODES.extend(tes)
            keyed_list_of_list_nodes2.append(ALL_NODES)
            
                
        else:
            #finds nodes with paired nodes higher than them. these cant be nes   
            not_possible_nes=[]
            possible_nes=[]
            just_nodes=[]
            z_max=0
            for edge in edge_list:
                node1=edge[0]
                node2=edge[1]
                node1_z=node1[2]
                node2_z=node2[2]
                if node1_z>node2_z:
                    not_possible_nes.append(node2)
                elif node2_z>node1_z:
                    not_possible_nes.append(node1)
                elif node2_z==node1_z:
                    not_possible_nes.append(node1)
                    not_possible_nes.append(node2)
                if node1_z>z_max:
                    z_max=node1_z
                if node2_z>z_max:
                    z_max=node1_z
                just_nodes.append(node1)
                just_nodes.append(node2)
            z_thresh=z_max-THRSH
            poss_nes1=list(set(just_nodes)-set(not_possible_nes))
            for j in poss_nes1:
                if j[2]>z_thresh:
                    possible_nes.append(j)
            #CHAMBERS CANNOT BE IN POSSIBLE NES
            
            deleted_nodes=[]
            node_list=keyed_list_of_list_nodes[i]
            for node2 in node_list:
                number=node2[1]
                node=node2[0]
                paired=0
                while paired==0:
                    #this ensures that if it was a chamber or previously it will continue to be
                    for ch in CH:
                        dist=calculate_distance(node, ch)
                        if dist==0:
                            k_node=tuple([node, number])
                            chambers.append(k_node)
                            cham_nums.append(number)
                            paired=1
                        
                    for ne in NE:
                        dist=calculate_distance(node, ne)
                        if dist==0:
                            k_node=tuple([node, number])
                            nes.append(k_node)
                            ne_nums.append(number)
                            paired=1
                    paired=1
                     ###Now check lens.
            unpaired_ch=list(set(LL_CHAM_NUMS[i-1])-set(cham_nums))
            #print(unpaired_ch)
            if len(unpaired_ch)>0:
                for ch in LL_CHAMS[i-1]:
                    number = ch[1]
                    if number not in unpaired_ch:
                        continue
                    node = ch[0]
                    linked=False
                    pair_ch_main=False
                    pair_ch_te=False
                    #this checks against the main node list first for the exact number
                    if number in unpaired_ch:
                        #print(number, ch, 'YEE')
                        for n in node_list:
                            if n in deleted_nodes:
                                continue
                            node2=n[0]
                            num2=n[1]
        
                            if number==num2:
                                if node2 in possible_nes:
                                    continue
                                if node2 in JU:
                                    k_node=[]
                                    k_node.append(n)
                                    node_list=list(set(node_list)-set(k_node))
                                    JU=list(set(JU)-set(k_node))
                                    TE=list(set(TE)-set(k_node))
                                    chambers.append(n)
                                    cham_nums.append(num2)
                                    pair_ch_main=True
                                else: #this if same number is a te - it 
                                    k_node=[]
                                    k_node.append(n)
                                    preserved_no=num2
                                    node_list=list(set(node_list)-set(k_node))
                                    TE=list(set(TE)-set(k_node))
                                    new_node=tuple([node2,L])
                                    print(num2)

                                    L+=1
                                    pair_ch_te=True
 
                        if pair_ch_te==True:
                            deleted_nodes.append(new_node)
                            number=preserved_no
                            node_list.append(new_node)
                            #print('NEW_NODE', new_node)
                            
                        if pair_ch_main==True:
                            continue
                        
                                
                        
                    min_dist=CH_THRSH
                    for ju in JU:
                        if ju in possible_nes:
                            continue
                        dist=calculate_distance(node, ju)
                        if dist<min_dist:
                            dist=min_dist
                            new_ch_node=ju
                            new_ch_number=number
                            new_ch=tuple([ju, number])
                            #print('hHH', new_ch)
                            linked=True
                    if linked==True:
                        #print('NEW CH=', number, new_ch_node)
                        for node in node_list:
                            dist1=calculate_distance(node[0], new_ch_node)
                            if dist1==0:
                                rm_node=[]
                                rm_node.append(node)
                                node_list=list(set(node_list)-set(rm_node))
                            
                               # print('REMOVED',node)
                                break
                        rm_node2=[]
                        rm_node2.append(new_ch_node)
                        #print(new_ch_node, number)
                        JU=list(set(JU)-set(rm_node2))
                        node_list=list(set(node_list)-set(rm_node))
                        chambers.append(new_ch)
                        cham_nums.append(new_ch_number)
                        

                        
            unpaired_ne=list(set(LL_NE_NUMS[i-1])-set(ne_nums))
            if len(unpaired_ne)>0:
                for ne in LL_NE[i-1]:
                    number = ne[1]
                    if number not in unpaired_ne:
                        continue
                    node = ne[0]
                    linked_ju=False
                    linked_te=False
                    pair_ne_main=False
                    #this checks against the main node list first for the exact number
                    if number in unpaired_ne:
                        for n in node_list:
                            node2=n[0]
                            num2=n[1]
                            if node2 in not_possible_nes:
                                continue
                            if number==num2:
                                k_node=[]
                                k_node.append(n)
                                node_list=list(set(node_list)-set(k_node))
                                JU=list(set(JU)-set(k_node))
                                TE=list(set(TE)-set(k_node))
                                nes.append(n)
                                ne_nums.append(num2)
                                pair_ne_main=True
                        if pair_ne_main==True:
                            continue
                    
                                
                    min_dist=THRSH
                    for ju in JU:
                        if ju in not_possible_nes:
                                continue
                        dist=calculate_distance(node, ju)
                        if dist<min_dist:
                            dist=min_dist
                            new_ne_node=ju
                            new_ne_number=number
                            new_ne=tuple([ju, number])
                            linked_ju=True
                    for te in TE:
                        if te in not_possible_nes:
                                continue
                        dist=calculate_distance(node, te)
                        
                        if dist<min_dist:
                            dist=min_dist
                            new_ne_node=te
                            new_ne_number=number
                            new_ne=tuple([te, number])
                            linked_te=True
                            linked_ju=False

                    if linked_ju==True:
                        for node in node_list:
                            dist1=calculate_distance(node[0], new_ne_node[0])
                            if dist1==0:
                                rm_node=[]
                                rm_node.append(node)
                                node_list=list(set(node_list)-set(rm_node))
                                break
                        rm_node=[]
                        rm_node.append(new_ne_node)
                        JU=list(set(JU)-set(rm_node))
                        #node_list=list(set(node_list)-set(new_ne_node))
                        nes.append(new_ne)
                        ne_nums.append(new_ne_number)
                    if linked_te==True:
                        for node in node_list:
                            dist1=calculate_distance(node[0], new_ne_node[0])
                            if dist1==0:
                                rm_node=[]
                                rm_node.append(node)
                                node_list=list(set(node_list)-set(rm_node))
                                #print('REMOVED',rm_node)
                                break
                        rm_node=[]
                        rm_node.append(new_ne_node)
                        TE=list(set(TE)-set(rm_node))
                        nes.append(new_ne)
                        ne_nums.append(new_ne_number)

            for node2 in node_list:
                
                number=node2[1]
                node=node2[0]
                paired=0
                while paired==0:
                    for ju in JU:
                        dist=calculate_distance(node, ju)
                        if dist==0:
                            k_node=tuple([node, number])
                            jus.append(k_node)
                            ju_nums.append(number)
                            paired=1

                    for te in TE:
                        dist=calculate_distance(node, te)
                        if dist==0:
                            k_node=tuple([node, number])
                            tes.append(k_node)
                            te_nums.append(number)
                            paired=1
                    paired=1

            LL_CHAMS.append(chambers)
            LL_CHAM_NUMS.append(cham_nums)
            LL_JU.append(jus)
            LL_JU_NUMS.append(ju_nums)
            LL_NE.append(nes)
            LL_NE_NUMS.append(ne_nums)
            LL_TE.append(tes)
            LL_TE_NUMS.append(te_nums)
            ALL_NODES=[]
            ALL_NODES.extend(chambers)
            ALL_NODES.extend(jus)
            ALL_NODES.extend(nes)
            ALL_NODES.extend(tes)
            keyed_list_of_list_nodes2.append(ALL_NODES)          
    return keyed_list_of_list_nodes2, LL_CHAMS, LL_CHAM_NUMS, LL_JU, LL_JU_NUMS, LL_NE, LL_NE_NUMS, LL_TE, LL_TE_NUMS
     
        
def make_edge_list_keyed (keyed_edge_list_dist_all, keyed_list_of_list_nodes2):
    keyed_edge_list_dist_all2=[]
    for el in range(0,len(keyed_edge_list_dist_all)):
        edge_list=keyed_edge_list_dist_all[el]
        edge_list2=[]
        for edge in edge_list:
            node1=edge[0][0]
            node2=edge[1][0]
            dist=edge[2]
            for node in keyed_list_of_list_nodes2[el]:
                dist1=calculate_distance(node[0], node1)
                if dist1==0:
                    keyed_node1=node

                dist2=calculate_distance(node[0], node2)
                if dist2==0:
                    keyed_node2=node
            edge_list2.append(tuple([keyed_node1, keyed_node2, dist]))
        keyed_edge_list_dist_all2.append(edge_list2)
    return keyed_edge_list_dist_all2     
        #if node not in LL_CHAMS[i-1]:
