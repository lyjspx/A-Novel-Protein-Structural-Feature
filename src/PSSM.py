import sqlite3
from pssmgen import PSSM
from Bio.PDB import *
import os
import pandas as pd
from src.utils import workingRoot

threeAndOne = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q',
                'GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K',
                'MET':'M','PHE':'F','PRO':'P','PYL':'O','SER':'S','SEC':'U',
                'THR':'T','TRP':'W','TYR':'Y','VAL':'V','ASX':'B','GLX':'Z',
                'XAA':'X','XLE':'J',
               'A':'ALA','R':'ARG','N':'ASN','D':'ASP','C':'CYS','Q':'GLN',
                'E':'GLU','G':'GLY','H':'HIS','I':'ILE','L':'LEU','K':'LYS',
                'M':'MET','F':'PHE','P':'PRO','O':'PYL','S':'SER','U':'SEC',
                'T':'THR','W':'TRP','Y':'TYR','V':'VAL','B':'ASX','Z':'GLX',
                'X':'XAA','J':'XLE'}

class Pssm:
    def __init__(self,workingRoot = workingRoot,defaultDB = 'PSSM.db',blastExe = '/usr/bin/psiblast',
                 blastDatabase=workingRoot+'uniref50.fasta'):
        self.workingRoot = workingRoot
        self.defaultDB = workingRoot + defaultDB
        self.blastExe = blastExe
        self.blastDatabase = blastDatabase

    def sql_connection(self):
        try:
            con = sqlite3.connect(self.defaultDB)
            #print(f"Connection is established: Database is created/connected in {self.defaultDB}")
            return con
        except Exception as e:
            print(e)

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

    def sql_create_table(self,tableName):
        con = self.sql_connection()
        sqlCommand = f'''CREATE TABLE '{tableName}' (
                        pdbID CHAR(10),
                        chain CHAR(4),
                        pdbNum INT,
                        resType CHAR(5),
                        A FLOAT, R FLOAT, N FLOAT, D FLOAT,C FLOAT,Q FLOAT,E FLOAT,G FLOAT,H FLOAT,I FLOAT,
                        L FLOAT,K FLOAT,M FLOAT,F FLOAT,P FLOAT,S FLOAT,T FLOAT,W FLOAT,Y FLOAT,V FLOAT,IC FLOAT
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

    def PSSM_to_DB(self,pdbID,chain,numThreads=8):
        tableName = '_'.join([pdbID,chain])
        if self.sql_table_exist(tableName):
            print(f"{tableName} is already in the DB")
            return
        '''
        check if pssm database contains the chain of PDB
        if yes, pass
        if not, use pssmgen to generate mapped pssm and upload
        :param pdbID:
        :return:
        '''
        PDB = PDBList()
        retrieve_status = PDB.retrieve_pdb_file(pdb_code=pdbID, pdir=self.workingRoot + pdbID + "/pdb",
                                                file_format="pdb")
        os.rename(retrieve_status, retrieve_status[0:(len(retrieve_status) - 3)] + 'pdb')
        gen = PSSM(work_dir=self.workingRoot + pdbID)
        gen.configure(blast_exe = self.blastExe,database=self.blastDatabase,num_threads=numThreads,
                      evalue=0.0001, comp_based_stats='T',
                max_target_seqs=2000, num_iterations=3, outfmt=7,
                save_each_pssm=True, save_pssm_after_last_round=True)
        gen.get_fasta(pdb_dir='pdb', chain=[chain], out_dir='fasta')
        gen.get_pssm(fasta_dir='fasta', out_dir='pssm_raw', run=True)
        gen.map_pssm(pssm_dir='pssm_raw', pdb_dir='pdb', out_dir='pssm', chain=[chain])
        gen.get_mapped_pdb(pdbpssm_dir='pssm', pdb_dir='pdb', pdbnonmatch_dir='pdb_nonmatch')
        finalPssmPath = f'{self.workingRoot}/{pdbID}/pssm/pdb{pdbID}.{chain}.pdb.pssm'
        if os.path.isfile(finalPssmPath):
            conn = self.sql_connection()
            cursor = conn.cursor()
            pssmMatrix = pd.read_table(filepath_or_buffer=finalPssmPath,
                                       comment="#",header= 0,sep='\s+')
            self.sql_create_table(tableName)
            try:
                for row in pssmMatrix.iterrows():
                    current_row = list(row)
                    sqlCommand = f'''INSERT INTO `{tableName}`(pdbID,chain,pdbNum,resType,A,R,N,D,C,Q,E,G,H,I,
                                                    L,K,M,F,P,S,T,W,Y,V,IC) 
                                                    Values(?,?,?,?,?,?,?,?,?,?,
                                                    ?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);'''
                    value = [pdbID, chain] + [list(current_row[1])[0]] + [threeAndOne[list(current_row[1])[1]]] + list(current_row[1])[4:]
                    cursor.execute(sqlCommand, value)
                conn.commit()
                print(f'{pdbID} {chain}\'s PSSM matrix was stored in DB')
            except Exception as e:
                print(e)
                print(f'error in {pdbID} {chain}')
            finally:
                cursor.close()
                conn.close()
        else:
            print(f'Failed to get PSSM matrix for {pdbID} {chain}')

    def get_PSSM(self,pdbID,chain,pdbNum):
        tableName = '_'.join([pdbID,chain])
        conn = self.sql_connection()
        cursor = conn.cursor()
        result = []
        try:
            sqlCommand = f''' SELECT * From `{tableName}` WHERE `pdbNum`={pdbNum}; '''
            cursor = conn.cursor()
            cursor.execute(sqlCommand)
            result = cursor.fetchall()
        except Exception as e:
            print(e)
        cursor.close()
        conn.close()
        return result

    def get_seq_neighbor_PSSM(self,pdbID,chain,pdbNum,N=12,includeSelf=False):
        tableName = '_'.join([pdbID,chain])
        conn = self.sql_connection()
        cursor = conn.cursor()
        result = []
        try:
            sqlCommand = f''' SELECT * From `{tableName}`; '''
            cursor = conn.cursor()
            cursor.execute(sqlCommand)
            result = cursor.fetchall()
        except Exception as e:
            print(e)
        cursor.close()
        conn.close()
        result = sorted(result,key= lambda x:x[2]) # sort rows by pdbNum
        # for these have less N neighbors, fill 0

        pdbNumList = [x[2] for x in result]
        totalRes = len(pdbNumList)
        currentPosInTable = pdbNumList.index(pdbNum)
        if not includeSelf:
            neighPSSM = [tuple([0 for i in range(20)]) for j in range(max(0,N - currentPosInTable))] + \
                        [result[i][4:24] for i in range(max(0,currentPosInTable - N),currentPosInTable )] + \
                        [result[i][4:24] for i in range(currentPosInTable + 1, min(currentPosInTable + 1 + N, totalRes))] + \
                        [tuple([0 for i in range(20)]) for j in range(max(0,currentPosInTable + 1 +  N - totalRes))]
        return neighPSSM




