import pandas as pd
import warnings
workingRoot ='/home/ubuntu/PTM/data/'

class Uniprot2PDB:
    def __init__(self,workingRoot=workingRoot):
        self.workingRoot = workingRoot
        self.pdbChainUniprot = pd.read_csv(workingRoot+'pdb_chain_uniprot.csv')
        self.validRow = (self.pdbChainUniprot["PDB_BEG"].str.isnumeric()) &\
                        (self.pdbChainUniprot["PDB_END"].str.isnumeric())
        self.pdbChainIntegrity = self.pdbChainUniprot[self.validRow]
    def is_pdb_info_available(self,uniprotAccession,position):
        rowWithAccession = self.pdbChainIntegrity[self.pdbChainIntegrity["SP_PRIMARY"].isin([uniprotAccession])]

        if (rowWithAccession.shape[0] < 1):
            warnings.warn(f'No mapping of Uniprot {uniprotAccession}')
            return []
        else:
            rowWithAccession = rowWithAccession[(rowWithAccession["SP_BEG"] < position)
                                                & (rowWithAccession["SP_END"] > position)]
            rowWithAccession.loc[:, "Ksites"] = position + (rowWithAccession["PDB_BEG"].astype(int)
                                                            - rowWithAccession["SP_BEG"].astype(int))
            return [(x[1][0],x[1][1],x[1][9]) for x in rowWithAccession.iterrows()]