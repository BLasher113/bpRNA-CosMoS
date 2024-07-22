from numpy.linalg import norm
from numpy import dot
import collections

def get_kmers(seq_1, seq_2, k):
    kmers = [] 
    kmers_list_1 = [] 
    kmers_list_2 = []
    n_kmers_1 = len(seq_1) - k + 1
    n_kmers_2 = len(seq_2) - k + 1
    for n_kmer in range(n_kmers_1):
        kmer = seq_1[n_kmer: n_kmer + k]
        kmers_list_1.append(kmer)
    for n_kmer in range(n_kmers_2):
        kmer = seq_2[n_kmer: n_kmer + k]
        kmers_list_2.append(kmer)
    kmers = set(kmers_list_1 + kmers_list_2)
    return kmers, kmers_list_1, kmers_list_2

def get_kmer_vectors(kmers_list_1, kmers_list_2, kmers, kmer_distances):
    V1 = []
    V2 = []
    if kmer_distances:
        kmer_counts_1 = dict(collections.Counter(kmers_list_1))
        kmer_counts_2 = dict(collections.Counter(kmers_list_2))
        for kmer in kmers:
            if kmer in kmer_counts_1.keys():
                counts_1 = kmer_counts_1[kmer]
            else: 
                counts_1 = 0
            if kmer in kmer_counts_2.keys():
                counts_2 = kmer_counts_2[kmer]
            else: 
                counts_2 = 0
            if kmer in kmer_distances.keys():
                pseudo_options = kmer_distances[kmer].keys()
                for p_kmer in pseudo_options:
                    e_dist = kmer_distances[kmer][p_kmer]
                    if p_kmer in kmer_counts_1.keys():
                        p_count_1 = kmer_counts_1[p_kmer]
                        counts_1 += p_count_1*e_dist
                    if p_kmer in kmer_counts_2.keys():
                        p_count_2 = kmer_counts_2[p_kmer]
                        counts_2 += p_count_2*e_dist
            V1.append(counts_1)
            V2.append(counts_2)
    else:
        for i, kmer in enumerate(kmers):
            counts_1 = kmers_list_1.count(kmer)
            counts_2 = kmers_list_2.count(kmer)
            V1.append(counts_1)
            V2.append(counts_2)
    return V1, V2

def get_score(V1, V2, w):
    cos_sim = (dot(V1, V2)/(norm(V1)*norm(V2)))*(1-w)
    return cos_sim

def get_similarity_score(ss_1, ss_2, k, kmer_distances = False):
    kmers, kmers_1, kmers_2 = get_kmers(ss_1, ss_2, k)
    #print(kmers, kmers_1, kmers_2)
    V1, V2 = get_kmer_vectors(kmers_1, kmers_2, kmers, kmer_distances)
    #print(V1, V2)
    len_1 = len(ss_1)
    len_2 = len(ss_2)
    l_dif = float(abs(len_1 - len_2))
    weight = 2*l_dif/(len_1 + len_2)
    sim_score = get_score(V1, V2, weight)
    return sim_score
