import src.Spatial_Feature
import pandas as pd

cdhit70ListPath = '/home/ubuntu/PTM/data/cdhit0.7_N_link.csv'
if __name__ == '__main__':
    print(cdhit70ListPath)
    cdhit70List = pd.read_csv(cdhit70ListPath)

    spatialFeature = src.Spatial_Feature.SpatialFeature()
    for row in cdhit70List.iterrows():
        print(row[1]["PDB"],row[1]["CHAIN"])
        spatialFeature.spatial_feature_to_DB(row[1]["PDB"],row[1]["CHAIN"])
