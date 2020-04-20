## importing pymol
from ipymol import viewer as pymol
import os

## importing numpy and pandas
import numpy as np
import pandas as pd

# setting wd
os.chdir('/Users/alan/gdrive/documentos/doutorado/isb/20200210-ribosomeHighlighting/')

##
##RUN 'pymol -R' from main console before starting server
##
pymol.start()   # Start PyMOL RPC server

## loading pdb ribosome models
pymol.do('load 4v6u.cif')
pymol.do('load 4v9f.cif')

## style/view preferences
pymol.do('show cartoon, all')

## printing info about chains
chains = pd.read_table('chainIDs.txt')
chains[0:]

## label proteins/chains
for i in np.arange(len(chains)):
        pymol.do('sele ' + chains.geneID[i] + ', ' + chains.PDB[i] + ' and chain ' + chains.Chains[i])

## align halo LSU to pyro LSU
pymol.do('sele pyrorRNA, 4v6u and chain B*')
pymol.do('align 4v9f, pyrorRNA')

## hide pyro LSU except two subunits (L1, L40e) and SSU L7ae, unknown SX chain A9 and halo.arcLX (chain 6)
pymol.do('hide everything, pyrorRNA or chain A3 or chain A9 or chain 6')
pymol.do('show spheres, rpl1')
pymol.do('show spheres, rpl40e')

# color subunits
pymol.do('color lightblue, 4v6u')
pymol.do('color deepsalmon, 4v9f')

## color RNA grey
pymol.do('color grey, chain A2 or chain 0 or chain 9') ##rRNA
pymol.do('color grey, chain A0 or chain A1') ##tRNA

## show proteins as spheres
pymol.do('show spheres, rp*')

## getting entries corresponding to ribosome proteins
## in our network co-ip model
chains.loc[chains['ID'].isin(['VNG1689G', 'VNG1692G', 'VNG0177G', 'VNG1433G', 'VNG1668G',])]

##Color chains corresponding to outlier and unClustered subunits (red and blue):
def color_chains(gene):
    pymol.do('color tv_red, chain ' +
    chains.Chains[gene] +
    '; as spheres, chain ' +
    chains.Chains[gene])

def color_chains_blue(gene):
    pymol.do('color marine, chain ' +
    chains.Chains[gene] +
    '; as spheres, chain ' +
    chains.Chains[gene])

color_chains_blue(11)
color_chains_blue(27)
color_chains(36)
color_chains(48)
color_chains(55)

## create image and rotated; save both
#pymol.do('reset')
#pymol.do('png ribosomeNetworkHighlights.png, dpi=900')
#pymol.do('rotate y,270')
#pymol.do('png ribosomeNetworkHighlights-270.png, dpi=900')
