import sys
sys.path.append("asap/py")
import pandas as pd
import re
from asap import *
from Bio import SeqIO
from scripts.extract_peptide_from_sequence import Seq
from scripts.add_custom_features import filter_pep, extend_df
from sklearn.externals import joblib
from sklearn.ensemble import RandomForestClassifier
import os

class ColumnSelector(object):
    """
    A feature selector for scikit-learn's Pipeline class that returns
    specified columns from a numpy array.
    """
    def __init__(self, cols):
        self.cols = cols

    def transform(self, X, y=None):
        return X[:, self.cols]

    def fit(self, X, y=None):
        return self

def get_prediction_input(in_file='ES_linked_VSGs.fasta'):
    res = []
    with open(in_file, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            record_id = record.id
            record_seq = str(record.seq)
            new_obj = Seq(record_seq)
            temp_res = new_obj.format_prediction(extend, pattern, record_id)
            if len(temp_res) == 0:
                print 'warning', record_id, 'is empty'
            else:
                res.append('\n'.join(temp_res))
    return '\n'.join(res)


 
#for the way asap work the 
extend = 10
pattern = 'N.[ST]'
#clf = joblib.load('models/voting_calibrated_classifier_rf_etc_svm.joblib.pkl')

def predict(temp__id, seq, clf):
    seq = str(seq).replace('*','A')
    window_extend = 10
    def windows_filter(window): 
        return window.get_aa_seq()[window_extraction_params.window_hot_index] == 'N'

    window_extraction_params = WindowExtractionParams(
                                window_prefix = window_extend, 
                                window_suffix = window_extend,
                                windows_filter = windows_filter)
    
    
    res = extract_windows_from_seq(seq, 
        window_extraction_params = window_extraction_params) 

    df = pd.read_csv(res)
    df = filter_pep(df)

    res_df=pd.DataFrame()
    if df.shape[0]>0:
        df = extend_df(in_df=df, hot_index=window_extend)
        
        #print df
        df_selection = df.iloc[:,4:]
        #print df_selection.head()
        # print df_selection.shape
        #print df_selection.shape()
        scores = clf.predict_proba(df_selection.values)
        #print scores
        res_df = pd.DataFrame()
        res_df['Prot']=[temp__id for n in df.loc[df_selection.index.values]['window_neighbourhood'] ]
        res_df['Seq']= [n for n in df.loc[df_selection.index.values]['window_seq']]
        res_df['Pos']= [n for n in df.loc[df_selection.index.values]['window_hot_index']]
        res_df['Score'] = scores[:,1]
        res_df['Prediction'] = ['Complex/Paucimannose' if n > 0.5 else 'Oligomannose' for n in scores[:,1]]
        new_pep = []
        for a,b in zip(res_df['Seq'],res_df['Score']):
            if b > 0.5:
                pep = a[0:10]+'['+a[10]+']'+a[11:]
            else:
                pep = a[0:10]+'['+a[10]+']'+a[11:]
            new_pep.append(pep)
        res_df['Seq']=new_pep


    return res_df


if __name__ == '__main__':
    pass




    

            

