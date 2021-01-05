# -*- coding: utf-8 -*-
"""
Created on Sat Dec 26 17:08:20 2020

@author: DIY
"""
import argparse
import pandas as pd

from get_sequence import ExtractSeq
from extract import ExtractFeatures
from selection import select


from sklearn.externals import joblib

###### main function ######
def main():
    parser = argparse.ArgumentParser(description ='PREDAIP: Computational Prediction and Analysis for Anti-inflammatory Peptide via a Hybrid Feature Selection Technique.')
    parser.add_argument('-input',  dest ='inputfile', type=str, help='Query peptide sequences to be predicted in txt format.', required=True)
    parser.add_argument('-threshold', dest='threshold_value', type=float, help='Please input a value between 0 and 1', required=True)
    parser.add_argument('-output',  dest='outputfile', type=str, help='Saving the prediction results in csv format.', required=False)
    args = parser.parse_args()
    inputfile = args.inputfile;
    threshold = args.threshold_value;
    outputfile = args.outputfile;
    
    print ('Peptide sequences are loading...')
    sequence,namelist = ExtractSeq(inputfile)
    
    print ('FusedFeatures are extracting...')
    data = ExtractFeatures(sequence,1,2,3,4,5,6,7,8,9,10,11,12)
    X_test = select(data)
    print ('Loading model...')
    model = joblib.load('best.model')
    prediction = model.predict_proba(X_test)
    n = int(len(sequence))
    peptide_name=[]
    Sequences=[]
    Probability=[]
    print ('AIPs were predicted as follows.')
    print ('-------------------------------')
    for i in range(n):
        if prediction[i][1] > threshold:
            peptide_name.append(namelist[i])
            Sequences.append(sequence[i])
            print (namelist[i])
            print (sequence[i])
            var = '%.11f' % prediction[i][1]
            Probability.append(var)
            print ('probability value:' + str(var))
            print ('-------------------------------')
        else:
            print (namelist[i])
            print (sequence[i])
            print ('The peptide is not AIP')
            print ('-------------------------------')

    
    AA = {'a':peptide_name,'b':Sequences,'c':Probability}
    AIP = pd.DataFrame(AA)
    AIP.to_csv(outputfile, index=False, header=['peptide name','Sequences','Probability value'])
            
if __name__ == "__main__":
    main() 