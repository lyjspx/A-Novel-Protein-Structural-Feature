#import external modules
import pandas as pd
import numpy as np
from sklearn import preprocessing
from sklearn.cluster import KMeans
from multiprocessing import Pool

#import internal modules
from src.utils import workingRoot

class NeighborProcessing:
    def __init__(self,fileName,workingRoot=workingRoot):
        self.filename = fileName
        self.workingRoot = workingRoot
        self.neighbors = pd.read_csv(workingRoot+fileName).iloc[:,1:]
        self.kClusterDic = {}
        self.importantClusters = []
        self.v_one_column = np.vectorize(self._one_column)
        self.v_norm_distance = np.vectorize(self._norm_distance)


    def get_similar_neighbor_in_each_sample(self,numNeigh=18,numberProcess=18):
        '''
        Find similar neighors for each neighbor within each site

        :param numNeigh:
        :param numberProcess:
        :return:
        '''
        totalDis = []
        with Pool(processes=numberProcess) as pool:
            for col in [3 * x for x in range(numNeigh)]:
                x = pool.starmap(self._wrap_novectorize, [(row,
                                                        col) for row in range(self.neighbors.shape[0])])
                temp = [[z for y in x[i] for z in y] for i in range(len(x))]
                tempList = pd.DataFrame(temp)
                tempList['order'] = list(range(tempList.shape[0]))
                totalDis.append(tempList)
            # for i in pool.map(self._wrap, [3 * x for x in range(numNeigh)]):
            #     tempList = []
            #     for x in i:
            #         temp = [z for y in x for z in y]
            #         tempList.append(temp)
            #     tempList = pd.DataFrame(tempList)
            #     tempList['order'] = list(range(tempList.shape[0]))
            #     totalDis.append(tempList)

        totalDis = pd.concat(totalDis)
        totalDis['numNum'] = [i for i in range(numNeigh) for j in range(tempList.shape[0])]
        totalDis['resName'] = np.ravel(self.neighbors.iloc[:,[3 * i     for i in range(numNeigh)]], order='F')
        totalDis['resDis'] = np.ravel(self.neighbors.iloc[:, [3 * i + 1 for i in range(numNeigh)]], order='F')
        totalDis['resAng'] = np.ravel(self.neighbors.iloc[:, [3 * i + 2 for i in range(numNeigh)]], order='F')
        return totalDis

    def get_group_dis(self,similarNeigh):
        similarNeigh["resName"] = similarNeigh["resName"].fillna('21')
        similarNeigh["resDis"] = similarNeigh["resDis"].fillna(999)
        similarNeigh["resAng"] = similarNeigh["resAng"].fillna(999)
        numNeigh = int((similarNeigh.shape[1]-5)/3)
        disMatrix = [self.v_norm_distance(similarNeigh.iloc[:,3*i+1],
                                            similarNeigh.iloc[:,3*i+2],
                                            similarNeigh.iloc[:,-2],
                                            similarNeigh.iloc[:,-1])
                                            for i in range(numNeigh)]
        return [sum(i) for i in np.column_stack(disMatrix)]

    def get_K_neighbor(self,neighInfo,K=2):
        numNeigh = int((neighInfo.shape[1] - 6) / 3)
        colIndex = [j for i in range(numNeigh) for j in (3*i+1,3*i+2)]
        resSegment = []
        for res in neighInfo['resName'].unique():
            oneRes = neighInfo[neighInfo['resName'] == res].copy(deep=True)
            kmeans = KMeans(n_clusters=min(oneRes.shape[0]-1,K)).fit(oneRes.iloc[:,colIndex])
            self.kClusterDic[res] = kmeans
            oneRes['Kmeans_cluster'] = kmeans.labels_
            resSegment.append(oneRes)
        result = pd.concat(resSegment)
        result["total_dis"][result['total_dis'] < 1] = max(result['total_dis']) * 10
        return result

    def get_final_feature(self,clusteredNearDf:pd.DataFrame,M=4):
        allClusters = clusteredNearDf.sort_values(by=['total_dis']).loc[:,
                      ['resName','Kmeans_cluster']].drop_duplicates().iloc[0:M,]
        self.importantClusters = [str(allClusters.iloc[i,0]) +\
                             str(allClusters.iloc[i,1]) for i in range(M)]
        numNeigh = int((clusteredNearDf.shape[1] - 7) / 3)
        finalFeature = []
        for i in range(numNeigh):
            oneSample = []
            oneData = clusteredNearDf[clusteredNearDf["order"]==i].copy(deep=True)
            oneData["cluster"] = oneData["resName"] + oneData["Kmeans_cluster"].astype('str')
            for keyCluster in self.importantClusters:
                if oneData[oneData['cluster'] == keyCluster].shape[0] == 0:
                    oneSample.extend(['unknown',999,999])
                else:
                    oneSample.extend(list(oneData[oneData['cluster']==keyCluster].sample(1).iloc[:,[-6,-5,-4]].values[0]))
            finalFeature.append(oneSample)
        return finalFeature
 
    def _one_column(self,res_name, res_dis, res_angle):
        return self.neighbors.apply(lambda x: self._process_one_row(x, res_name, res_dis, res_angle), axis=1)

    def _wrap(self,i):
        '''
        deprecated, using too much memory
        :param i:
        :return:
        '''
        return self.v_one_column(self.neighbors.iloc[:,i],
                                 self.neighbors.iloc[:,i+1],
                                 self.neighbors.iloc[:,i+2])

    def _wrap_novectorize(self, row, column):
        return self._one_column(self.neighbors.iloc[row,column],
                                 self.neighbors.iloc[row,column+1],
                                 self.neighbors.iloc[row,column+2])

    def _process_one_row(self,one_row_protein, res_name, res_dis, res_angle):
        one_row = []
        length = sum(one_row_protein == res_name)
        if length < 1:
            one_row.append(21)
            one_row.append(999)
            one_row.append(999)
        elif length == 1:
            for index in one_row_protein[one_row_protein == res_name].index:
                index = int(index)
                one_row.append(one_row_protein[index - 1])
                one_row.append(one_row_protein[index])
                one_row.append(one_row_protein[index + 1])
        else:
            candidate_index = []
            compare = []
            compare_candidate = []
            for index in one_row_protein[one_row_protein == res_name].index:
                index = int(index)
                candidate_index.append(index)
                compare.append(np.sqrt((one_row_protein[index] - res_dis) ** 2 \
                                       + (one_row_protein[index + 1] - res_angle) ** 2))
                compare_candidate.append(one_row_protein[index - 1])
                compare_candidate.append(one_row_protein[index])
                compare_candidate.append(one_row_protein[index + 1])
            min_index = compare.index(min(compare))
            one_row.append(compare_candidate[min_index * 3])
            one_row.append(compare_candidate[min_index * 3 + 1])
            one_row.append(compare_candidate[min_index * 3 + 2])
        return one_row

    def _norm_distance(self,ang1,dis1,ang2,dis2):
        normArray = preprocessing.normalize([[ang1, ang2], [dis1, dis2]], axis=0)
        totalDis = np.sqrt((normArray[0][0] - normArray[0][1]) ** 2 + \
                              (normArray[1][0] - normArray[1][1]) ** 2)
        return (totalDis)