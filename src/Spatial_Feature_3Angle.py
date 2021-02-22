import sqlite3
from Bio.PDB import *
import os
import pdbx
import pandas as pd
import numpy as np
from src.utils import workingRoot
from src.PSSM import threeAndOne

backBone = ["N", "CA", "C", "O"]
aminoAcidCodes = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLY", "GLU", "HIS", "ILE", "LEU", "LYS",
                  "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]


class SpatialFeature:
    def __init__(self, workingRoot=workingRoot, defaultDB='Spatial_feature_3Angle.db'):
        self.workingRoot = workingRoot
        self.defaultDB = workingRoot + defaultDB

    def sql_connection(self):
        try:
            con = sqlite3.connect(self.defaultDB)
            # print(f"Connection is established: Database is created/connected in {self.defaultDB}")
            return con
        except Exception as e:
            print(e)

    def PDB_downloader(self, pdbID):
        PDB = PDBList()
        if os.path.exists(self.workingRoot + pdbID + '.cif'):
            print(f'{pdbID} mmcif file already downloaded')
            return True
        retrieve_status = PDB.retrieve_pdb_file(pdb_code=pdbID, pdir=self.workingRoot, file_format="mmCif")
        if os.path.exists(retrieve_status):
            print(f'{pdbID} mmcif file just downloaded')
            return True
        else:
            print(f'{pdbID} mmcif file fails to be downloaded')
            return False

    def sql_table_exist(self, tableName) -> bool:
        cursor = self.sql_connection().cursor()
        cursor.execute(f''' SELECT count(name) FROM sqlite_master
                        WHERE type='table' AND name='{tableName}' ''')
        if cursor.fetchone()[0] == 1:
            cursor.close()
            return True
        else:
            cursor.close()
            return False

    def sql_create_table(self, tableName):
        con = self.sql_connection()

        '''
        some definitions
        distance: R group - R group
        angle: dihedral angle between Calpha(1)-Rcenter(1) and Calpha(1)-Rcenter(2)
        angleCa: dihedral angle between Calpha(1)-Rcenter(1) and Calpha(1)-Calpha(2)
        angleC: dihedral angle between Calpha(1)-Rcenter(1) and Calpha(1)-C(2)
        '''
        sqlCommand = f'''CREATE TABLE '{tableName}' (
                        pdbID CHAR(10),
                        chain CHAR(4),
                        pdbNum1 INT,
                        authNum1 INT,
                        resType1 CHAR(5),
                        pdbNum2 INT,
                        authNum2 INT,
                        resType2 CHAR(5),
                        distance FLOAT,
                        angle FLOAT,
                        angleCa FLOAT,
                        angleC FLOAT
                        ); '''
        try:
            con.cursor().execute(sqlCommand)
            con.commit()
            print(f'{tableName} created')
        except Exception as e:
            print(e)
            print(f'errors in create table {tableName}')
        finally:
            con.close()

    def spatial_feature_to_DB(self, pdbID, chain):
        tableName = '_'.join([pdbID, chain])
        if self.sql_table_exist(tableName):
            print(f"{tableName} is already in the DB")
            return
        if not self.PDB_downloader(pdbID):
            return
        with open(self.workingRoot + pdbID + '.cif', 'r') as f:
            cifData = pdbx.load(f)[0]
        atomList = pd.DataFrame(cifData.get_object('atom_site')._row_list
                                , columns=cifData.get_object('atom_site')._attribute_name_list)
        ############
        # Resolution not found in mmcif
        # Leave space here
        #####
        atomList = atomList[atomList['auth_asym_id'] == chain]
        if atomList.shape[0] < 1:
            print(f'chain {chain} not found in {pdbID} cif file')
            return
        res_info_chain = np.empty(2000, dtype=[('pdbID', 'S5'), ('chain', 'S2'),
                                               ('type', 'S4'), ('pdbNum', 'int16'),
                                               ('authNum', 'int16'), ('center', 'float16,float16,float16'),
                                               ('direction', 'float16,float16,float16'),
                                               ('caPosition', 'float16,float16,float16'),
                                               ('cPosition', 'float16,float16,float16')],)
        allPDBNum = atomList['label_seq_id'].unique()
        i = 0
        for PDBNum in allPDBNum:
            atomList = atomList[atomList["pdbx_PDB_model_num"] == '1']
            oneRes = atomList[atomList["label_seq_id"] == PDBNum]
            if oneRes['label_comp_id'].iloc[0] not in aminoAcidCodes:
                continue
            if oneRes['label_comp_id'].iloc[0] != 'GLY':
                try:  # in case some components missing
                    rGroup = oneRes[oneRes["label_atom_id"].isin(backBone)]
                    rGroupX = sum(pd.to_numeric(rGroup['Cartn_x'])) / rGroup.shape[0]
                    rGroupY = sum(pd.to_numeric(rGroup['Cartn_y'])) / rGroup.shape[0]
                    rGroupZ = sum(pd.to_numeric(rGroup['Cartn_z'])) / rGroup.shape[0]
                    cToRGroup = np.subtract((float(oneRes[oneRes['label_atom_id'] == 'CA']['Cartn_x']),
                                             float(oneRes[oneRes['label_atom_id'] == 'CA']['Cartn_y']),
                                             float(oneRes[oneRes['label_atom_id'] == 'CA']['Cartn_z'])),
                                            (rGroupX, rGroupY, rGroupZ))

                    res_info_chain[i] = (pdbID, chain, oneRes['label_comp_id'].iloc[0],
                                         PDBNum, oneRes["auth_seq_id"].iloc[0], (rGroupX, rGroupY, rGroupZ),
                                         tuple(cToRGroup),
                                         (float(oneRes[oneRes['label_atom_id'] == 'CA']['Cartn_x']),
                                          float(oneRes[oneRes['label_atom_id'] == 'CA']['Cartn_y']),
                                          float(oneRes[oneRes['label_atom_id'] == 'CA']['Cartn_z'])),
                                         (float(oneRes[oneRes['label_atom_id'] == 'C']['Cartn_x']),
                                          float(oneRes[oneRes['label_atom_id'] == 'C']['Cartn_y']),
                                          float(oneRes[oneRes['label_atom_id'] == 'C']['Cartn_z']))                                         )
                    i += 1
                except Exception as e:
                    print(e)
            else:
                try:
                    rGroupX = float(oneRes[oneRes['label_atom_id'] == 'C']['Cartn_x']) + \
                              float(oneRes[oneRes['label_atom_id'] == 'N']['Cartn_x']) + \
                              float(oneRes[oneRes['label_atom_id'] == 'O']['Cartn_x'])
                    rGroupY = float(oneRes[oneRes['label_atom_id'] == 'C']['Cartn_y']) + \
                              float(oneRes[oneRes['label_atom_id'] == 'N']['Cartn_y']) + \
                              float(oneRes[oneRes['label_atom_id'] == 'O']['Cartn_y'])
                    rGroupZ = float(oneRes[oneRes['label_atom_id'] == 'C']['Cartn_z']) + \
                              float(oneRes[oneRes['label_atom_id'] == 'N']['Cartn_z']) + \
                              float(oneRes[oneRes['label_atom_id'] == 'O']['Cartn_z'])
                    cToRGroup = np.subtract((float(oneRes[oneRes['label_atom_id'] == 'CA']['Cartn_x']),
                                             float(oneRes[oneRes['label_atom_id'] == 'CA']['Cartn_y']),
                                             float(oneRes[oneRes['label_atom_id'] == 'CA']['Cartn_z'])),
                                            (rGroupX, rGroupY, rGroupZ))
                    res_info_chain[i] = (pdbID, chain, oneRes['label_comp_id'].iloc[0],
                                         PDBNum, oneRes["auth_seq_id"].iloc[0], (rGroupX / 3, rGroupY / 3, rGroupZ / 3),
                                         tuple(cToRGroup),
                                         (float(oneRes[oneRes['label_atom_id'] == 'CA']['Cartn_x']),
                                          float(oneRes[oneRes['label_atom_id'] == 'CA']['Cartn_y']),
                                          float(oneRes[oneRes['label_atom_id'] == 'CA']['Cartn_z'])),
                                         (float(oneRes[oneRes['label_atom_id'] == 'C']['Cartn_x']),
                                          float(oneRes[oneRes['label_atom_id'] == 'C']['Cartn_y']),
                                          float(oneRes[oneRes['label_atom_id'] == 'C']['Cartn_z'])))
                    i += 1
                except Exception as e:
                    print(e)
        res_info_chain = res_info_chain[0:i]
        self.sql_create_table(tableName)
        conn = self.sql_connection()
        cursor = conn.cursor()
        try:
            for res1 in range(i):
                for res2 in range(i):
                    corrdinatesSubstract = np.subtract(tuple(res_info_chain[res1]["center"]),
                                                       tuple(res_info_chain[res2]["center"]))
                    distance = np.sqrt(np.sum(corrdinatesSubstract ** 2))
                    vector1 = np.array(tuple(res_info_chain[res1]["direction"]))
                    vector2 = np.subtract(tuple(res_info_chain[res2]['center']),
                                          tuple(res_info_chain[res1]['caPosition']))
                    vector3 = np.subtract(tuple(res_info_chain[res2]['caPosition']),
                                          tuple(res_info_chain[res1]['caPosition']))
                    vector4 = np.subtract(tuple(res_info_chain[res2]['cPosition']),
                                          tuple(res_info_chain[res1]['caPosition']))

                    if abs(distance) > 0.0001:  # Exclue 0 distance, e.g. self to self
                        angle = 180 * np.arccos(0.99 * np.dot(vector1, vector2) / (
                                np.linalg.norm(vector1) * np.linalg.norm(vector2))) / np.pi
                        angleCa = 180 * np.arccos(0.99 * np.dot(vector1, vector3) / (
                                np.linalg.norm(vector1) * np.linalg.norm(vector3))) / np.pi
                        angleC = 180 * np.arccos(0.99 * np.dot(vector1, vector4) / (
                                np.linalg.norm(vector1) * np.linalg.norm(vector4))) / np.pi
                    else:
                        angle = 0
                        angleCa = 0
                        angleC = 0
                    sqlCommand = f''' INSERT INTO `{tableName}`
                            (pdbID , chain,  pdbNum1 ,authNum1,resType1 ,
                                pdbNum2 ,authNum2 ,resType2 ,distance,angle,angleCa,angleC)
                    VALUES ('{pdbID}','{chain}',
                            '{res_info_chain[res1]['pdbNum']}',
                            '{res_info_chain[res1]['authNum']}',
                            '{res_info_chain[res1]['type'].decode('UTF-8')}',
                            '{res_info_chain[res2]['pdbNum']}',
                            '{res_info_chain[res2]['authNum']}',
                            '{res_info_chain[res2]['type'].decode('UTF-8')}',
                            {distance},{angle},{angleCa},{angleC})
                                   '''
                    cursor.execute(sqlCommand)
        except Exception as e:
            print(e)

        conn.commit()
        cursor.close()
        conn.close()
        print(f'{pdbID} {chain} done')
        return

    def all_neighbor_fetch(self, pdbID, chain, pdbNum, sort=True,
                           sortByDis=True):  # sortByDis: True - by dis| False - by Seq
        conn = self.sql_connection()
        tableName = '_'.join([pdbID, chain])
        if not self.sql_table_exist(tableName):
            print(f"{tableName} is not yet in the DB")
            return
        sqlCommand = f''' SELECT * From `{tableName}` WHERE `pdbNum1`={pdbNum}; '''
        cursor = conn.cursor()
        cursor.execute(sqlCommand)
        result = cursor.fetchall()
        cursor.close()
        conn.close()
        if sort:
            if sortByDis:
                return sorted(result, key=lambda x: x[-2])
            else:
                return sorted(result, key=lambda x: x[5])
        else:
            return result

    def pdbNum_mapping_resNum(self, pdbID, chain, resNumToPdbNum=True):
        conn = self.sql_connection()
        tableName = '_'.join([pdbID, chain])
        if not self.sql_table_exist(tableName):
            print(f"{tableName} is not yet in the DB")
            return
        sqlCommand = f''' SELECT distinct `pdbNum2` ,`authNum2` From `{tableName}`; '''
        cursor = conn.cursor()
        cursor.execute(sqlCommand)
        result = cursor.fetchall()
        cursor.close()
        conn.close()
        if resNumToPdbNum:
            return {v: k for k, v in dict(result).items()}
        else:
            return dict(result)

    def get_N_nearest_neighbor(self, pdbID, chain, pdbNum, N=12, byDistance=True):
        if byDistance:  # N defines total neighbors got
            allNeighbor = self.all_neighbor_fetch(pdbID=pdbID,
                                                  chain=chain,
                                                  pdbNum=pdbNum,
                                                  sort=True, sortByDis=True)
            return [(x[5], x[6], x[7], x[8], x[9], x[10], x[11]) for x in allNeighbor[1:N + 1]]
        else:  # N define number of upstream or downstream neighbor; thus 2*N neighbor returned
            allNeighbor = self.all_neighbor_fetch(pdbID=pdbID,
                                                  chain=chain,
                                                  pdbNum=pdbNum,
                                                  sort=True, sortByDis=False)
            totalRes = len(allNeighbor)
            pdbNumList = [x[5] for x in allNeighbor]
            currentPosInTable = pdbNumList.index(pdbNum)
            return [(-9, -9, "O", 0, 0, 0, 0) for j in range(max(0, N - currentPosInTable))] + \
                   [(allNeighbor[i][5], allNeighbor[i][6], allNeighbor[i][7],
                     allNeighbor[i][8], allNeighbor[i][9],allNeighbor[i][10], allNeighbor[i][11]) \
                    for i in range(max(0, currentPosInTable - N), currentPosInTable)] + \
                   [(allNeighbor[i][5], allNeighbor[i][6], allNeighbor[i][7],
                     allNeighbor[i][8], allNeighbor[i][9], allNeighbor[i][10], allNeighbor[i][11]) \
                    for i in range(currentPosInTable + 1, min(currentPosInTable + 1 + N, totalRes))] + \
                   [(-1, -1, "O", 0, 0, 0, 0) for j in range(max(0, currentPosInTable + 1 + N - totalRes))]

    def get_sequence_segment(self, pdbID, chain, pdbNum, N=12, byDistance=False, asString=True):
        if not byDistance:
            allNeighbor = self.all_neighbor_fetch(pdbID=pdbID, chain=chain, pdbNum=pdbNum, sort=True, sortByDis=False)
            totalRes = len(allNeighbor)
            pdbNumList = [x[5] for x in allNeighbor]
            currentPosInTable = pdbNumList.index(pdbNum)
            resList = [('O') for j in range(max(0, N - currentPosInTable))] + \
                      [(allNeighbor[i][7]) for i in range(max(0, currentPosInTable - N),
                                                          min(currentPosInTable + 1 + N, totalRes))] + \
                      [('O') for j in range(max(0, currentPosInTable + 1 + N - totalRes))]
            if asString:
                return ''.join([threeAndOne[x] for x in resList])
