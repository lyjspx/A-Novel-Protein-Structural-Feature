from src import PSSM, Spatial_Feature
from src import others_feature
from src import utils
from Bio.PDB import *
from src import order_neighbors
import pandas as pd
import time

if __name__ == '__main__':
    start_time = time.time()
    #PSSM.PSSM_to_DB(pdbID='1apy',chain='B')
    # p = PDBParser()
    # structure = p.get_structure("X",'/home/ubuntu/PTM/data/1apy/pdb/pdb1apy.pdb')
    # print(structure[0]["A"].get_unpacked_list()[0].get_id()[1])
    # print(structure[0]["A"].get_unpacked_list()[0].get_resname())
    # x = structure[0]["A"].get_residues()
    # print(next(x).get_full_id())
    # print(next(x).get_full_id())
    # f = open("/home/ubuntu/download/mmcif_pdbx/1APY.cif",'r')
    # from pdbx import PdbxSyntaxError
    # from pdbx.reader import PdbxReader
    # from pdbx import loads as read_cifstr
    # reader = PdbxReader(f)
    # print(reader)
    # data_list = []
    # reader.read(data_list)
    # block = data_list[0]
    # print(block.__dict__)
    # struct_sheet = block["struct_sheet"]
    # print(struct_sheet)
    x = Spatial_Feature.SpatialFeature()
    print(x.pdbNum_mapping_resNum('4l3t','A',resNumToPdbNum=True))
    # print(x.get_sequence_segment("2bkr","B",15,N=12))
    # print(len(x.get_sequence_segment("2bkr","B",15,N=12)))

    #print(x.all_neighbor_fetch('2bkr','B',6,sort=True ,sortByDis=False))
    #print(x.pdbNum_mapping_resNum('1apy','B'))
    #print(x.get_N_nearest_neighbor('1apy','B',4,12,False))
    # y = PSSM.Pssm()
    #
    # print(y.get_seq_neighbor_PSSM('2bkr','A',210,N=12))
    # for i in range(10):
    #     print(others.get_AAindex("ALA"))
    #     print(others.get_AAindex("GLY"))

    # uniport = utils.Uniprot2PDB()
    # print(uniport.is_pdb_info_available("P24182",27))

    # z = order_neighbors.NeighborProcessing('PTM_k_NN24.csv')
    # x = z.get_similar_neighbor_in_each_sample()
    # print(x.shape)
    # print(x.head())
    # nearDf = x
    # nearDf["total_dis"] = z.get_group_dis(x)
    # nearDfCopy = z.get_K_neighbor(nearDf,K=2)
    # print(pd.DataFrame(z.get_final_feature(nearDfCopy,M=6)))

    print("--- %s seconds ---" % (time.time() - start_time))