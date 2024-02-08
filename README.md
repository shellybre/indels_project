# indels_project

the code runs only in linux. 
in order to run it in windows you should remove the 'phylip' ('phylip neighbor' change to 'neighbor',  'phylip treedist' change  to 'treedist'). 
you should have phylip installed in order to run it. 
in simulations : in order to run hp or dcj you should download UniMoG (https://bibiserv.cebitec.uni-bielefeld.de/dcj). 
put UniMoG in same folder as the code. 
you should have the input files in order to run the scripts. 

#indles_simulation : 
the main function is : createExcel_for_tree_multiple (methods,name, indels, l, add,n, Gsize). 
the function gets array of methods you want to simulate (0=hp, 1=cga, 2=LCS, 3=g_ab, 4=jaccard, 5=dcj, 6=si & 1-si (go together to preserve time as they do the same calculation)). 
name = the name you want for the files. 
indels = from 0-1.1, when 0 = indels only, 1.1=jump only. any value between 0-1 is combination of the two (as closer to 0 - more indels, as closer to 1 - more jumps). 
l = first mean branch length. 
add = by how much the branch length grow (for example, l = 0.1, add=0.2. run for 3 times : simulate for mean branch length of 0.1, 0.3, 0.5).
n = number of leaves in the simulated tree. 
Gsize = mean genome size of nodes in the tree). 
the function defualt are : l = 0.001, add = 0.001,n=20, Gsize=5000.
the function creates text file and excel file with all the resuls. 
this function calles for function : createTreeSimulationsMultiple (methods, name,indels, l=0.001, add=0.001,n=20, Gsize=5000). this function runs for 10 mean brnach length and 50 repeats for each (for each method). 
you can change the number of branch lengths and the number of repeates in the first 2 lines in function (itr=50, limit=10).


