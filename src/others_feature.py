from propy.AAIndex import GetAAIndex1, get
from sklearn.preprocessing import LabelBinarizer
from .PSSM import threeAndOne
from .Spatial_Feature import aminoAcidCodes

aminoAcidPropertyByRes = {"ALA":0,"ILE":0,"LEU":0,"MET":0,"VAL":0,
                      "PHE":1, "TRP":1,"TYR":1,
                       "ASN":2,"CYS":2,"GLN":2,"SER":2,"THR":2,
                      "ASP":3,"GLU":3,"ARG":4,"HIS":4,"LYS":4,"GLY":5,"PRO":5,"unknown":6}

def get_AAindex(resType):
    #aaindex.grep('physicochemical')
    num1AAindex = get('WILM950101')[threeAndOne[resType]]
    #aaindex.grep('hydropho')
    num2AAindex = get('ARGP820101')[threeAndOne[resType]]
    #aaindex.grep('polarity')
    num3AAindex = get('GRAR740102')[threeAndOne[resType]]
    #polarizability
    num4AAindex = get('CHAM820101')[threeAndOne[resType]]
    #Hydration potential
    num5AAindex = get('WOLR810101')[threeAndOne[resType]]
    #Accessibility reduction ratio
    num6AAindex = get('PONP800107')[threeAndOne[resType]]
    #Net charge
    num7AAindex = get('KLEP840101')[threeAndOne[resType]]
    #Molecular weight
    num8AAindex = get('FASG760101')[threeAndOne[resType]]
    #PK-N
    num9AAindex = get('FASG760104')[threeAndOne[resType]]
    #PK-C
    num10AAindex = get('FASG760105')[threeAndOne[resType]]
    #melting point
    num11AAindex = get('FASG760102')[threeAndOne[resType]]
    #optical rotation
    num12AAindex = get('FASG760103')[threeAndOne[resType]]
    #entropy of formation
    num13AAindex = get('HUTJ700103')[threeAndOne[resType]]
    #heat capacity
    num14AAindex = get('HUTJ700101')[threeAndOne[resType]]
    #absolute entropy
    num15AAindex = get('HUTJ700102')[threeAndOne[resType]]
    return (num1AAindex,num2AAindex,num3AAindex,num4AAindex,num5AAindex,num6AAindex,\
           num7AAindex,num8AAindex,num9AAindex,num10AAindex,num11AAindex,num12AAindex,\
           num13AAindex,num14AAindex,num15AAindex)

def onehot_encoder_res(resType, numClass=20): #support 20 or 6
    if numClass not in [6, 20]:
        print(f'{numClass} is not supported')
        return
    elif numClass == 20:
        return LabelBinarizer().fit_transform(aminoAcidCodes)[aminoAcidCodes.index(resType)]
    else:
        return LabelBinarizer().fit_transform(range(7))[aminoAcidPropertyByRes[resType]]


