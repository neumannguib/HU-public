#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd
import allel
import biotite.sequence.phylo as phylo
import scipy.spatial
import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


breeds={'DSN':[]}
breed_names={'DSN':[]}
bulls_project=pd.read_excel("/home/guilherme/NAS/DNA/1000bulls/Run9-TAUraw-animalDistributionList-20210630.xlsx")
bulls_project.index=bulls_project['SampleID']
callset = allel.read_vcf("/halde/guilherme/Run9/paper/ChrY.vcf.gz",alt_number=1)
for i,sample in enumerate(callset['samples'][:]):
    try:
        sample_breed = bulls_project.loc[sample,'Breed']
        if sample_breed=='FjällCattle':
            sample_breed='Fjäll'
        if sample_breed=='TraditionalAngler':
            sample_breed='GermanRedAngler'
        if sample_breed=='HolsteinFriesian':
            sample_breed='Holstein'
        if  sample_breed not in list(breeds.keys())+['DeutschesSchwarzbuntesNiederungsrind']:
            breeds[sample_breed]=[]
            breed_names[sample_breed]=[]
        if sample_breed=='DeutschesSchwarzbuntesNiederungsrind':
            breeds['DSN'].append(i)
            breed_names['DSN'].append(sample)
        else:
            breeds[sample_breed].append(i)
            breed_names[sample_breed].append(sample)
    except:
        print(sample)
for key in breeds.keys():
    print(key,len(breeds[key]))


# In[6]:


df=pd.DataFrame()
df['ID']=breed_names['MeuseRhineYssel'] + breed_names['RedWhiteDualPurpose'] + breed_names['Rotbunt']+breed_names['DSN']
df.to_csv('/halde/guilherme/rotbunt_dsn_ids.txt',index=False,sep='\t',header=False)
df


# In[54]:


for chro in list(range(1,30))+['X']:
    command="vcftools --gzvcf /halde/guilherme/Run9/paper/Chr"+str(chro)+".vcf.gz --recode --mac 1 --max-missing 0.90 --keep /halde/guilherme/rotbunt_dsn_ids.txt --stdout | bgzip -c > /halde/guilherme/Run9/Rotbunt_DSN/Chr"+str(chro)+".vcf.gz &"
    os.system(command)


# In[3]:


breeds={'DSN':[],'Red White Dual Purpose':[]}
breed_names={'DSN':[],'Red White Dual Purpose':[]}
bulls_project=pd.read_excel("/home/guilherme/NAS/DNA/1000bulls/Run9-TAUraw-animalDistributionList-20210630.xlsx")
bulls_project.index=bulls_project['SampleID']
callset = allel.read_vcf("/halde/guilherme/Run9/Rotbunt_DSN/Chr1.vcf.gz",alt_number=1)
samples=[]
for i,sample in enumerate(callset['samples'][:]):
    sample_breed = bulls_project.loc[sample,'Breed']
    if  sample_breed not in list(breeds.keys())+['DeutschesSchwarzbuntesNiederungsrind']:
        breeds[sample_breed]=[]
        breed_names[sample_breed]=[]
    if sample_breed=='DeutschesSchwarzbuntesNiederungsrind':
        breeds['DSN'].append(i)
        breed_names['DSN'].append(sample+'_DSN')
        sample_breed='DSN'
    else:
        breeds[sample_breed].append(i)
        breed_names[sample_breed].append(sample+"_"+sample_breed)
    samples.append(sample+'_'+sample_breed)
    if sample_breed=='MeuseRhineYssel' or sample_breed=='Rotbunt' or sample_breed=='RedWhiteDualPurpose':
            sample_breed='Red White Dual Purpose'
            breeds[sample_breed].append(i)
            breed_names[sample_breed].append(sample+'_Red White Dual Purpose')
    
for key in breeds.keys():
    print(key,len(breeds[key]))


# In[4]:


gt = allel.GenotypeChunkedArray(callset['calldata/GT'])
chroms=list(callset['variants/CHROM'])
positions=list(callset['variants/POS'])
soma=len(gt)
calls=pd.DataFrame(gt.to_n_alt(),columns=samples)
dist =scipy.spatial.distance.pdist(calls.T, lambda u, v: np.nansum([np.abs(v-u)]))
dist_all=scipy.spatial.distance.squareform(dist)  
for chro in list(range(2,30)):    
    callset = allel.read_vcf("/halde/guilherme/Run9/Rotbunt_DSN/Chr"+str(chro)+".vcf.gz",alt_number=1)
    gt = gt.concatenate([allel.GenotypeChunkedArray(callset['calldata/GT'])], axis=0)
    soma+=len(gt)
    chroms+=list(callset['variants/CHROM'])
    positions+=list(callset['variants/POS'])
    dist=scipy.spatial.distance.pdist(calls.T, lambda u, v: np.nansum([np.abs(v-u)]))
    dist=scipy.spatial.distance.squareform(dist)
    dist_all=dist_all+dist


# In[62]:


gt


# In[5]:


tree = phylo.upgma(dist_all)
file=open('/halde/guilherme/tree_run9_upgma_dsn_rotbunt.txt','w')
file.write(tree.to_newick(labels=samples,include_distance=True))
file.close()


# In[8]:


import seaborn as sns
sns.clustermap(pd.DataFrame(dist_all,columns=samples,index=samples), standard_scale=1,figsize=(27, 29))
plt.show()


# In[7]:


coords1, model1 = allel.pca(gt_all.to_n_alt(), n_components=4, scaler='patterson')
def plot_pca_coords(coords, model, pc1, pc2, ax, sample_population):
    pop_colours={"MeuseRhineYssel":'r',"RedWhiteDualPurpose":'g',"DSN":'b',"Rotbunt":'orange'}
    sns.despine(ax=ax, offset=5)
    x = coords[:, pc1]
    y = coords[:, pc2]
    sample_population=np.array(sample_population)
    for pop in ['DSN','MeuseRhineYssel',"RedWhiteDualPurpose","Rotbunt"]:
        flt = (sample_population == pop)
        ax.plot(x[flt], y[flt], marker='o', linestyle=' ', color=pop_colours[pop], label=pop, markersize=6, mec='k', mew=.5)
    ax.set_xlabel('PC%s (%.1f%%)' % (pc1+1, model.explained_variance_ratio_[pc1]*100))
    ax.set_ylabel('PC%s (%.1f%%)' % (pc2+1, model.explained_variance_ratio_[pc2]*100))
    

def fig_pca(coords, model, title, sample_population=None):
    if sample_population is None:
        sample_population = df_samples.population.values
    # plot coords for PCs 1 vs 2, 3 vs 4
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(1, 2, 1)
    plot_pca_coords(coords, model, 0, 1, ax, sample_population)
    ax = fig.add_subplot(1, 2, 2)
    plot_pca_coords(coords, model, 2, 3, ax, sample_population)
    ax.legend(bbox_to_anchor=(1, 1), loc='upper left')
    fig.suptitle(title, y=1.02)
    fig.tight_layout()
pop = ["DSN"]*302+ ['MeuseRhineYssel']*25 + ['RedWhiteDualPurpose']*17 + ['Rotbunt']*10
fig_pca(coords1, model1, 'PCA',pop)


# In[10]:


for chro in list(range(1,30))+['MT']:
    command='bgzip /halde/guilherme/Run9/paper/prunned_'+str(chro)+'_raw.vcf &'
    os.system(command)


# In[ ]:




