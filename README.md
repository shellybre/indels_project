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
Gsize = mean genome size of nodes in the tree. 
the function defualt are : l = 0.001, add = 0.001,n=20, Gsize=5000.
the function creates text file and excel file with all the resuls. 
this function calles for function : createTreeSimulationsMultiple (methods, name,indels, l=0.001, add=0.001,n=20, Gsize=5000). this function runs for 10 mean brnach length and 50 repeats for each length (and for each method. so if you chose 2 method it will simulate for 2 X 10 X 50 = 1,000 times.). 
you can change the number of branch lengths and the number of repeates in the first 2 lines in function (itr=50, limit=10).

input file : phylip-param.txt (essential for the run). this file is essential for running phylip treedist (it answers the quaestion of treedist according to what we need). 

output files : 8 or 10 files (10 if you run dcj/hp, 8 otherwise). 2 main files are interesting : excel (with the name you gave) and file with the name : results_ + name you gave. 
excel file has the results for every method and every mean branch length checked : mean rf results, mean normalized rf results, std of rf and std of normalized rf. 
the other file have all the results from the run (all the 50 values of the repeats). in the output_simulations example the code ran for 2 repeats and 2 branch lengths. 
the other files you see : output files from phylip or dcj/hp runs (for the very last simulation that ran,as each simulation those files overwritten).

#realDataIndels : 
the main function is create_new_Tree_from_real(name,k). 
this function gets the name of the tree you want to check and which method you want (-1=cga, 6=si). 

In order to create the real tree you must have folder (with the name you gave for the function) that contains two files : one called 'cogs', that contains the list of cogs for each leaf in the tree. 
the other one called 'listCogs' that contains the newick file (the name have to be the same as the name of the main file, with .nhx ending). see example in input_real_data for this kind of file (in the example, listCogs have more files, you only need the .nhx one).  also replace the path in the beggining of the script to the path where this file is loacted. 
you also must have the file phylip-param.txt (essential for the run). this file is essential for running phylip treedist (it answers the quaestion of treedist according to what we need). this file is in the input files.
this function creates new tree from the distance matrix created given the cogs and the method. it also calculates rf between the given and the new tree. 

output : you have 6 output files. 'treedist'+name : contains the treedist output (rf result). 'second'+name : contains newick of the new tree. 'new_tree'+name contains the structure of this tree. 'infile': distance matrix. 'first'+name:the original newick. 'dist.out' : treedist description of the run. 








