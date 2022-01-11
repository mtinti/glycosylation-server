

from jinja2 import Environment, FileSystemLoader
import os
import re

# import Tornado
import tornado.ioloop
import tornado.web
import pandas as pd
import numpy as np
import pickle
from string import strip
import json
import itertools
#from PIL import Image
import io
#import pylab
import socket
import urlparse
import pandas as pd
from protein import Seq,Prot
from HTMLParser import HTMLParser
from Bio import SeqIO
from StringIO import StringIO  
from sklearn.externals import joblib
from predict_seq import predict


#Complex/Paucimannose , not biantennary
#Oligomannose, not triantennary 

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

clf = joblib.load('models/voting_calibrated_classifier_rf_etc_svm.joblib.pkl')



_PATH = os.path.dirname(os.path.abspath(__file__))
templateLoader = FileSystemLoader(searchpath="templates")
templateEnv = Environment(loader=templateLoader)
settings = {'debug': True,}
static_path=os.path.join(_PATH, 'static')


gene_table = pd.DataFrame.from_csv('data/gene_table_28.txt', sep='\t')
print gene_table.tail()


mapping_table = pd.DataFrame.from_csv('data/conv_28_54_GeneByLocusTag_Summary.txt', sep='\t')
#conv_ids = {}
conv_desc = dict(zip(mapping_table['Input ID'],mapping_table['Product Description']))
new_desc = []
for prot_id, desc  in zip(gene_table.index.values, gene_table['description']):
    if prot_id in conv_desc:
        new_desc.append(conv_desc[prot_id])
    else:
        new_desc.append(desc)
gene_table['description']=new_desc



pkl_file = open('data/Proteins2.pkl', 'rb')
prot_dict = pickle.load(pkl_file)
pkl_file.close()


table_predicted_protein = []
protein_selection  = [prot_dict[n].prot_id for n in prot_dict]# if prot_dict[n].hasPrediction()]
for n in  protein_selection:
    #print(n)
    temp_res = gene_table.loc[n][['description','comments','go']]
    temp_res = temp_res.fillna('')
    temp_res = [n]+list(temp_res.values)
    table_predicted_protein.append(temp_res)
table_predicted_protein = {'dataTable':table_predicted_protein}   


class MainHandler(tornado.web.RequestHandler):
    def get(self):
        TEMPLATE_FILE = "prot.html"
        template = templateEnv.get_template(TEMPLATE_FILE)
        index_html_output = template.render(title="welcome to the complex navigator")
        self.write(index_html_output)

class FetchData(tornado.web.RequestHandler):
    def get(self):
        data_tag = str(strip(self.request.uri.split('/')[-1]))
        if data_tag == 'trip_ids':
            self.write(json.dumps(prot_dict.keys()))

class DumpProtein(tornado.web.RequestHandler):

    def set_default_headers(self):
        print("setting headers!!!")
        self.set_header("access-control-allow-origin", "*")
        self.set_header("Access-Control-Allow-Headers", "x-requested-with")
        self.set_header('Access-Control-Allow-Methods', 'GET, PUT, DELETE, OPTIONS')
        # HEADERS!
        self.set_header("Access-Control-Allow-Headers", "access-control-allow-origin,authorization,content-type") 

    def options(self):
        # no body
        self.set_status(204)
        self.finish()

    def get(self):
        protein_id = str(strip(self.request.uri.split('/')[-1])).split('?_')[0]
        if protein_id in prot_dict:
            descriptions = gene_table.loc[protein_id]
            temp_prot = prot_dict[protein_id]          
            temp_dict = temp_prot.dump()
            temp_dict['go']=str(descriptions['go'])
            temp_dict['description']=descriptions['description']          
        else:
            temp_dict = {'id':protein_id,
                         'go':[],
                        'description':'not found',
                         'seq':'',
                         'Deamidated':[],
                         'HexNAc':[],
                        'Prediction':[],
                        'SignalPeptide':[],
                        'ExtraCellular':[],
                         'TransMembrane':[],
                          'Coverage':[],
                         'dataTable':[],

                         }        
        #print temp_prot.dump()
        #for key in temp_dict:
            #print key
            #print json.dumps(temp_dict[key],allow_nan=False)
        self.write(json.dumps(temp_dict,allow_nan=False))


class DumpPredictedProteins(tornado.web.RequestHandler):
    def get(self):
        self.write(json.dumps(table_predicted_protein,allow_nan=False))

class Prediction(tornado.web.RequestHandler):
    def get(self):
        TEMPLATE_FILE = "pred.html"
        template = templateEnv.get_template(TEMPLATE_FILE)
        index_html_output = template.render(title="welcome to prediction")
        self.write(index_html_output)


class PredictProtein(tornado.web.RequestHandler):

    def set_default_headers(self):
        print("setting headers!!!")
        self.set_header("access-control-allow-origin", "*")
        self.set_header("Access-Control-Allow-Headers", "x-requested-with")
        self.set_header('Access-Control-Allow-Methods', 'GET, PUT, DELETE, OPTIONS')
        # HEADERS!
        self.set_header("Access-Control-Allow-Headers", "access-control-allow-origin,authorization,content-type") 
    
    def options(self):
        # no body
        self.set_status(204)
        self.finish()  

    def get(self):
        #https://stackoverflow.com/questions/66923711/how-to-parse-fasta-string-using-seq-io-from-biopython
        #h = HTMLParser()
        #fasta_seq = h.unescape(str(self.request.uri))

        
        fasta_seq = str(self.request.uri)
        #print(fasta_seq)

        if len(fasta_seq) > 2000:
            self.write(json.dumps({'table': 'seq too long'}, allow_nan=False))



        fasta_seq = fasta_seq.replace("%3Cbr/%3E",'\n')
        fasta_seq = fasta_seq.replace("%3E",'>')
        fasta_seq = fasta_seq.replace("%20",' ')
        fasta_seq = fasta_seq.split('/')[-1]



        if fasta_seq[0] != '>':
            self.write(json.dumps({'error': 'sfasta seq'}, allow_nan=False))        


        fasta_io = StringIO(fasta_seq) 
        records = SeqIO.parse(fasta_io, "fasta") 
        for record in records:
            break

        #print(record)
        predictions_df = predict(record.id, record.seq, clf)
        #print(predictions_df)
        html_table = predictions_df.to_html()
        html_table = html_table.replace('table border', 
        'table id=\"html_table\" border')


        self.write(json.dumps({'id':record.id, 'seq':str(record.seq), 'table': html_table}, allow_nan=False))

        
handlers=[
    #(r"/prediction", Prediction),
    (r"/home", MainHandler),
    (r"/prediction", Prediction),
    (r"/data/trip_ids", FetchData),
    (r"/dumpProtein/.*", DumpProtein),
    (r"/dumpPredictedProteins", DumpPredictedProteins),
    (r"/predictProtein/.*", PredictProtein),
    (r"/static/(.*)", tornado.web.StaticFileHandler, {'path':static_path} ),
]

application = tornado.web.Application(handlers , **settings)   
  

if __name__ == "__main__":
    # Setup the server
    print 'start'
    PORT=8070
    application.listen(PORT)
    tornado.ioloop.IOLoop.instance().start()
