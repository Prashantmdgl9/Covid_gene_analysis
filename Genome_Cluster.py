#Imports
import numpy as np
import os.path
import itertools
from sklearn import preprocessing
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
%matplotlib inline
import matplotlib
import matplotlib.pyplot as plt
import turicreate as tc

#Files can be downloaded from:
# https://www.ncbi.nlm.nih.gov/nuccore/NC_019843.3?report=fasta
# https://www.ncbi.nlm.nih.gov/nuccore/MN908947.3?report=fasta

path="Linear_Algebra/Covid Gene Analysis/COVID-19.fa"

def loadGenome(path):
    genome=""
    with open(path, "r") as f:
        _ = f.readline().lower() # ignore the header
        l = f.readline().lower()
        while l:
            genome += l[:-1] # ignore new line
            l = f.readline().lower()
    genome_length=len(genome)
    return genome, genome_length

genome, genome_length = loadGenome(path)

block_length=300
combination_length=[1,2,3,4]
genome_substrings=genome_length//block_length
genome_substrings
#Divide the genome
genome_blocks = []
start_idx = 0
end_idx = block_length
for i in range(genome_substrings):
    chunk = genome[start_idx:end_idx]
    genome_blocks.append(chunk)
    start_idx = start_idx + block_length
    end_idx = start_idx + block_length
    #print(chunk)

permutations_dicts = []
def gene_permutation_generation():

    for ws in range(len(combination_length)):
        #print(ws)
        #print(type(ws))
        gene_permutations = [''.join(ele) for ele in itertools.product('actg', repeat=(ws+1))]
        permutations_dict = {key : val for val,key in enumerate(gene_permutations)}
        permutations_dicts.append(permutations_dict)
    return permutations_dicts

permutations_dict = gene_permutation_generation()
permutations_dict

def num_table(genome,genome_substrings,combination_length):
    arrays = [np.zeros((genome_substrings, 4**x)) for x in combination_length]
    for wl in combination_length:  #1,2,3,4
        k = 0
        for block in genome_blocks: # 1018 blocks
        #for b in range(genome_substrings):
            #block = genome[b*block_length:(b+1)*block_length]
            for i in range(block_length//wl): #block length is 300. 300, 150, 100, 75
                word = block[i*wl:(i+1)*wl]
                col = permutations_dict[wl-1][word]
                arrays[wl-1][k,col] += 1
            k = k + 1
    return arrays

arr = num_table(genome,genome_substrings,combination_length)
arr[3].shape
#Standardize the data
def standardize(arr, combination_length):
    for wl in combination_length:
        std_scale = preprocessing.StandardScaler().fit(arr[wl-1])
        arr[wl-1] = std_scale.transform(arr[wl-1])
    return arr

arr = standardize(arr, combination_length)

# Principal Component Analysis

def computePrincipalComponents(arr, combination_length, n_components=2):
    arrPCA={}
    for wl in combination_length:
        pca = PCA(n_components=n_components).fit(arr[wl-1])
        arrPCA[wl-1] = pca.transform(arr[wl-1])
    return arrPCA

arrPCA = computePrincipalComponents(arr, combination_length)

len(arrPCA)
#fig_size = [8, 8]
#plt.rcParams["figure.figsize"] = fig_size
def plot_pca():
    plt.figure(figsize=(16,16))
    fig, axes = plt.subplots(2,2)
    #colors = np.random.rand(99)
    axes[0,0].scatter(arrPCA[0].T[0], arrPCA[0].T[1], s=50, alpha = 0.5, c = 'tab:green')
    axes[0,1].scatter(arrPCA[1].T[0], arrPCA[1].T[1], s=50, alpha = 0.5)
    axes[1,0].scatter(arrPCA[2].T[0], arrPCA[2].T[1], s=50, alpha = 0.5, c = 'tab:orange')
    axes[1,1].scatter(arrPCA[3].T[0], arrPCA[3].T[1], s=50, alpha = 0.5, c = 'tab:gray')


plot_pca()


kmeans = KMeans(n_clusters=3, random_state=0).fit(arr[2])

def cluster_plt():
    fig, ax = plt.subplots(1,1)
    ax.scatter(arrPCA[2].T[0], arrPCA[2].T[1], s=50, c=kmeans.labels_, alpha = 0.5,)

cluster_plt()
