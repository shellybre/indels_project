import math
from scipy.optimize import fsolve
from random import randint
from numpy import *
import warnings
import os

path = "/home/user/git_indels/"
class Node:
  def __init__(self):
    self.name = ""
    self.children = []
    self.lengthS = ""
    self.length = 0
    self.parent = None
    self.genom = []
    self.distances = []
    self.genomLength = []
class Tree:
    def __init__(self):
        self.leaves = []
        self.root = None
        self.total = []
    def getLeaves (self):
        for x in range(len(self.leaves)):
            print(self.leaves[x].name)
    def getRandomLeaf (self):
        r1 = random.randint(0,len(self.leaves) -1)
        return r1


#calculate si
def SyntenyIndex_indels(Genom1, Genom2, k=6):
    total = (len(Genom1)+len(Genom2))/2
    t = 0
    for item in range(k, len(Genom1) - k):
        if (Genom1[item] in Genom2 and Genom1[item]!=-1):
          set1 = set(Genom1[item - k:item])
          set2 = set(Genom1[item + 1:item + k + 1])
          set3 = set1.union(set2)
          item1 = Genom2.index(Genom1[item])
          set12 = set(Genom2[item1 - k:item1])
          set22 = set(Genom2[item1 + 1:item1 + k + 1])
          set32 = set12.union(set22)
          set3_filtered = {x for x in set3 if x != -1}
          set32_filtered = {x for x in set32 if x != -1}
	# Calculate the intersection of the filtered sets
          intersection = set3_filtered.intersection(set32_filtered)
          len_intersection = len(intersection)
          b = (1 / (2 * k))
          JSI = b * len_intersection
          t += JSI
    if t==0:
      t=1      
    #total = (len(Genom1)+len(Genom2))/2
    syntenyindex = (1 / total) * t
    return syntenyindex, total

def Xi_indels(Genom1, Genom2):  
    sharedPairs=0
    for x in range (len(Genom1)-1):
      if (Genom1[x] in Genom2 and Genom1[x]!=-1 and Genom1[x+1]!=-1):
         indx = Genom2.index(Genom1[x])
         if(indx+1!=len(Genom2)):
            if(Genom2[indx+1]==Genom1[x+1]):
               sharedPairs+=1
    L0 = (len(Genom1) + len(Genom2))/ 2
    if (sharedPairs==0):
      sharedPairs=1          
    P00=sharedPairs/(L0-1)
    '''
    try:
      xi=0.5*((1-1/P00**0.5)*log(P00))**0.5
    except:
      xi=LCS(Genom1, Genom2)
      return xi 
    '''
    xi=0.5*((1-1/P00**0.5)*log(P00))**0.5   
    itr=100
    for itern in range(itr):
       xi=-(log(1+xi)+log(P00))/2            
    return xi, sharedPairs

#for k = 6 in SI
def calcExactDistFromSi_6(t, synIn):
  if (synIn == 0.0):
               print ("fgfgg")
               return(1000.0)
  #if (synIn == 1):
  #        return 0 
  else:
    return(math.exp(2*t)*synIn - (6*t**10 + 30*t**9 + 115*t**8 + 230*t**7 + 362*t**6 + 362*t**5 + 280*t**4 + 140*t**3 + 50*t**2 + 10*t + 1)/(t + 1)**11)    
    
#get = array of nodes
#calculate for each two nodes: fsolve (SI). creates : matrix of the results, and string of the matrix for neighbor, and writes to file called "infile"
#returns the matrix and the string 
def calculate_len_for_SI (arr, k):
    siLength = [[0]*len(arr) for i in range(len(arr))]  
    commonGenes = 0
    SI = 0
    genomeLen = 0
    compares = 0
    for x in range (len(arr)):
      genomeLen+= len(arr[x].genom)      
    for x in range (len(arr)):
       for y in range (x, len(arr)):
          if (x == y):
            siLength[x][x] = 0
          else: 
           compares+=1  
           result, total = SyntenyIndex_indels(arr[x].genom, arr[y].genom, k)
           SI += result
           commonGenes += total
           if (result == 0.0):
                dist = [] 
                dist.append(1000)
                #dist.append(1)
           elif (result == 1):
                 dist = []
                 dist.append(0)
           else: 
             dist = [] 
             if (k==10):  
                dist = fsolve(calcExactDistFromSi, 1.0, args = result) #1.0 is the initial value for t
             elif (k == 6):
                dist = fsolve(calcExactDistFromSi_6, 1.0, args = result)  #1.0 is the initial value for t
           siLength[x][y] = round(dist[0],6) 
           siLength[y][x] = round(dist[0],6)      
    newStr = str(len(arr)) + "\n"      
    for x in range (len(arr)):
       newStr += str(arr[x].name)
       for y in range (len(arr)):
            newStr+= " " + str(siLength[x][y])  
       newStr+="\n"     
    f = open("infile", "w")
    f.write(newStr)
    f.close()  
    return SI / compares, genomeLen / len(arr) , commonGenes/ compares           


#get = array of nodes
#calculate for each two nodes: fsolve (SI). creates : matrix of the results, and string of the matrix for neighbor, and writes to file called "infile"
#returns the matrix and the string 
def calculate_len_for_newForm (arr):
    siLength = [[0]*len(arr) for i in range(len(arr))]  
    commonGenes = 0
    SI = 0
    genomeLen = 0
    compares = 0
    for x in range (len(arr)):
      genomeLen+= len(arr[x].genom)      
    for x in range (len(arr)):
       for y in range (x, len(arr)):
          if (x == y):
            siLength[x][x] = 0
          else:  
           compares+=1  
           if(arr[x].genom == arr[y].genom):
             result=0
             total=len(arr[x].genom)
           else:  
             result, total = Xi_indels(arr[x].genom, arr[y].genom)
           SI += result
           commonGenes += total
           siLength[x][y] = round(result,6) 
           siLength[y][x] = round(result,6)      
    newStr = str(len(arr)) + "\n"      
    for x in range (len(arr)):
       newStr += str(arr[x].name)
       for y in range (len(arr)):
            newStr+= " " + str(siLength[x][y])  
       newStr+="\n"     
    f = open("infile", "w")
    f.write(newStr)
    f.close()  
    return SI / compares, genomeLen / len(arr) , commonGenes/ compares        

#func to create string of newick format from tree 
def createNewickFormat (root,newick):
    if (len(root.children) == 0):
        if (root.length != 0):        
            newick += str(root.name) + ":" + str("{:.3f}".format(root.length))
            return newick
    else:
        newick += "("
        newick = createNewickFormat(root.children[0], newick)
        newick += ","
        newick = createNewickFormat(root.children[1], newick)
        newick += ")"
        if (root.length != 0):
           newick+= ":" + str("{:.3f}".format(root.length))
    return newick   


#function to add length(float) to nodes 
def createLength (root):
  if root == None :
    return
  if (root.parent != None):
      root.length = float(root.lengthS)
  x = 0     
  while (len(root.children) > x):
        createLength(root.children[x])
        x+=1

#function to add genoms to leaves in tree
#get : Tree object that has leaves, name= the name of the file we want to read from
#change path to the correct one.
def addTotree (tree, name):
    global path
    path1 = path + name + '/cogs/' 
    files = []
    for (root, dirs, file) in os.walk(path1):
      #open all files inside our file  
      for f in file:
            file_path = f"{path1}\{f}"
            print('addTotree    ', file_path)
            file_path = file_path.replace('\\', '')
            with open(file_path, 'r') as ft:
              name = ft.read()
              print('addTotree, open ', file_path)
            x = 0
            #value = the name of the leaf
            value = ''
            #find the name of the leaf (strip unnecessary chars as "cogs_") 
            while x < len(f):
                while x < len(f) and not f[x].isdigit():
                  x+=1
                while x < len(f) and f[x].isdigit():
                  value+=str(f[x])
                  x+=1
                break
            count = 0
            #find the leaf that contains the same name 
            for x in range (len(tree.leaves)):
                  if value in tree.leaves[x].name:
                      tempstr = ""
                      y = 0
                      while y < len(name):
                          #if we see "," then we encounter new number, we save the previous number
                          if (name[y] == ","):
                              tree.leaves[x].genom.append(tempstr)
                              y+=1
                              tempstr = ""
                          else:
                              tempstr+= name[y]
                              y+=1
                      break

                            
#create new tree structure
#get : nwstr = newick format string, tree = instance of Tree object, count = pointer to char in the string 'nwstr', prnt = Node's parent to be, in the first call its None (because root doesnt have parent)
#returns : at each iteration returns the last node we finished creating, last iteration : the root of the tree;  count = pointer to char in the string 'nwstr'                    
def createTree (nwstr,tree, count, prnt):
    while (len(nwstr) > count and (nwstr[count] == '\n' or nwstr[count] == ' ')):
         count+=1
    #if we at the end of the string, return      
    if (len(nwstr) == count):
      return None, count
    #if we see "(" , its a new node
    if (nwstr[count] == "("):
       root = Node()
       root.parent= prnt
       while (len(nwstr) > count and (nwstr[count] == '\n' or nwstr[count] == ' ')):
         count+=1
       count+=1
       #call function again, the func will return the child of that node 
       child, count = createTree (nwstr,tree, count, root)
       root.children.append(child)
       #while wee see "," , there are more children  
       while (nwstr[count] == ','):
         while (len(nwstr) > count and (nwstr[count] == '\n' or nwstr[count] == ' ')):
           count+=1
         count+=1
        #call function again, the func will return the child of that node 
         child, count = createTree (nwstr,tree, count, root)
         root.children.append(child)
       count+=1
       while (len(nwstr) > count and (nwstr[count] == '\n' or nwstr[count] == ' ')):
         count+=1
       #if we at the end of the string, return      
       if (len(nwstr) <= count):
        return root,count
       else :
        if (nwstr[count]  == ":"):
          count +=1
        #we want to read the name, we search till we see  ":", knowing after that comes the length   
        else :
          while (count < len (nwstr) and nwstr[count]  != ":"):
            while (len(nwstr) > count and (nwstr[count] == '\n' or nwstr[count] == ' ')):
               count+=1
            if (len(nwstr) <= count ):
              break
            root.name += nwstr[count]
            count+=1
          count+=1
        while (len(nwstr) > count and (nwstr[count] == '\n' or nwstr[count] == ' ')):
           count+=1
        #we read the length, stoped when we see "," or ")"   
        while (len(nwstr) > count and (nwstr[count] != "," and nwstr[count] != ")")):
          while (len(nwstr) > count and (nwstr[count] == '\n' or nwstr[count] == ' ')):
               count+=1
          root.lengthS +=nwstr[count]
          count+=1
       return root,count
    #if we dont see "(", the new node is a leaf
    else:
        leaf = Node()
        leaf.parent = prnt
        while (len(nwstr) > count and (nwstr[count] == '\n' or nwstr[count] == ' ')):
           count+=1
        if (nwstr[count]  == ":"):
          count +=1
        else :
        #we want to read the name, we search till we see  ":", knowing after that comes the length   
          while (nwstr[count]  != ":"):
            while (nwstr[count] == '\n' or nwstr[count] == ' '):
                count+=1
            leaf.name += nwstr[count]
            count+=1
          count+=1
       #we read the length, stoped when we see "," or ")"     
        while (len(nwstr) > count and (nwstr[count] != "," and nwstr[count] != ")")):
          while (nwstr[count] == '\n' or nwstr[count] == ' '):
                count+=1
          leaf.lengthS +=nwstr[count]
          count+=1
        #add the new leaf to tree.leaves  
        tree.leaves.append(leaf)
        return leaf, count

#function to create tree with genoms to each leaf
#get : the name of the tree (ATGC..)
#return : Tree object, we can access by it the root of the tree        
def templateTree(name):
  try :
    #open txt contains newick format
    f = open(path + str(name) + '/' + 'listCogs/' + str(name) + ".nhx", "r")                        
  except FileNotFoundError:
    print('FileNotFoundError ', path + str(name))
    return None
  else:
    strTree = f.read()
    #removes unnecessary gaps 
    strTree = strTree.strip('\n')
    strTree = strTree.strip('\t')
    strTree = strTree.strip(' ')
    #create new Tree object
    tree = Tree()
    #call for func that retruns the root of the tree, and pointer to the end of the string
    root, count = createTree(strTree, tree, 0, None)
    tree.root = root    
    createLength(tree.root)
    #add genoms to leaves in tree
    addTotree(tree, name)
    count=0
    #counts how many leaves have genom 
    for y in range (len(tree.leaves)):
        if len(tree.leaves[y].genom) > 0:
          count += 1
    print('templateTree returns tree ', name, ' with ', str(count) , "leaves")    
    return tree
    

def create_new_Tree_from_real(name,k):
   tree = templateTree(name)
   for x in range (len(tree.leaves)):
      print(tree.leaves[x].name)
      tree.leaves[x].name = "species_" + str(x)
      print(tree.leaves[x].name)
   newstr = createNewickFormat(tree.root, "") + ";" 
   f = open("intree", "w")
   f.write(newstr)
   f.close()     
   if(k==-1):
   	SI, genomeLen, commonGens = calculate_len_for_newForm(tree.leaves)
   else:
        SI, genomeLen, commonGens = calculate_len_for_SI(tree.leaves,k)
   #call neighbor
   os.system("echo y| phylip neighbor")
   #rename files 
   old_name = r"outfile"
   new_name = r"new_tree" + str(name) + ".txt"
   if os.path.exists(new_name):
    os.remove(new_name)
   os.rename(old_name, new_name)
   old_name = r"outtree"
   new_name =  r"intree2"
   os.rename(old_name, new_name)
   #call treedist
   stat = os.system("phylip treedist < phylip-param.txt > dist.out.txt"); 
   #rename files 
   old_name = r"outfile"
   if(k==-1):
     new_name = r"treedistCGA" + str(name) + ".txt"
   else:  
     new_name = r"treedist" + str(name) + ".txt"
   if os.path.exists(new_name):
    os.remove(new_name)
   os.rename(old_name, new_name)
   old_name = r"intree"
   new_name = r"first" + str(name) + ".txt"
   if os.path.exists(new_name):
    os.remove(new_name)
   os.rename(old_name, new_name)
   old_name = r"intree2"
   if(k==-1):
        new_name = r"secondCGA" + str(name) + ".txt"
   else:
      new_name = r"second" + str(name) + ".txt"
   if os.path.exists(new_name):
    os.remove(new_name)
   os.rename(old_name, new_name)
   return tree, SI, genomeLen, commonGens


#use with name of tree, -1 => cga, 6 => si 
create_new_Tree_from_real('ATGC007',6)
create_new_Tree_from_real('ATGC007',-1)
