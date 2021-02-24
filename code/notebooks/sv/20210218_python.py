#!/usr/bin/env python
# coding: utf-8

# # Load libraries

# In[2]:


import os
import glob
import pandas as pd
import numpy as np


# # Read data 

# In[3]:


#proj_dir = "/hps/research1/birney/users/ian/mikk_paper/mikk_genome"
proj_dir = "/Users/brettell/Documents/Repositories/mikk_genome"


# In[4]:


data_dir = os.path.join(proj_dir, "data/sv_analysis/20210217_sviper_filter_pass")
target_files = glob.glob(data_dir + "/*.csv")
target_files


# In[5]:


# Read files
dfs = list()
for f in target_files:
  sample = os.path.basename(f).strip('.csv')
  df = pd.read_csv(f,
                   header = None,
                   names = ['CHROM', 'POS', 'ALT', 'SVLEN', 'SVTYPE', 'CHR2', 'END','GT','LN','ST'],
                   dtype = {'CHROM' : str,
                            'POS' : int,
                            'ALT' : str, 
                            'SVLEN' : int, 
                            'SVTYPE' : str,
                            'CHR2' : str, 
                            'END' : int, 
                            'GT' : str, 
                            'LN' : int, 
                            'ST' : str})
  df['SAMPLE'] = sample
  dfs.append(df)
    
# bind together
full_df = pd.concat(dfs)


# In[6]:


full_df.head()


# In[7]:


full_df.describe(include='all')


# # Plot

# In[9]:


import plotly.express as px


# In[10]:


fig = px.violin(full_df, x="SAMPLE", y="SVLEN", color="SAMPLE", box=True,facet_row = "SVTYPE")
fig.show()


# In[12]:


full_df.loc[full_df["CHROM"] == "2", :]


# ### Types of SV

# In[24]:


fig = px.histogram(full_df,
                   x="SVTYPE", color="SVTYPE",
                   template = "simple_white",
                   labels=dict(SVTYPE="STRUCTURAL VARIANT TYPE"))
fig.update_layout(showlegend=False)
fig.show()


# In[ ]:





# In[23]:


fig = px.histogram(full_df,
                   x="SVTYPE", color="SVTYPE",
                   template = "simple_white",
                   labels=dict(SVTYPE="STRUCTURAL VARIANT TYPE",
                               count = "COUNT"),
                   facet_row = "SAMPLE")
fig.update_layout(showlegend=False)
fig.show()


# In[19]:


fig = px.histogram(full_df, x="SVLEN", color="SAMPLE",facet_row = "SVTYPE")
fig.show()


# In[ ]:




