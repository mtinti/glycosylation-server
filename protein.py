from Bio import SeqIO
import numpy as np
import pickle
from string import strip
import json
import io
import pandas as pd
import re
import csv

#oligomannose/biantennary (D)
#paucimannose/triantennary 

class Seq():
    def __init__(self, seq):
        self.seq = seq
        
    def extract_seq(self, position, extend):
        seq = str(self.seq)
        if position < extend and position >= len(seq)-extend:
            add_front = extend-position
            add_back = extend - (len(seq) - position-1)
            res = add_front*'-'+seq+'-'*add_back
        
        elif position >= extend and position < len(seq)-extend:
            res = seq[position-extend:position+extend+1]
        
        elif position < extend:
            res = seq[:position+extend+1]
            add = extend-position
            res = '-'*add + res
        
        elif position >= len(seq)-extend:
            res = seq[position-extend:]
            add = extend - (len(seq) - position-1)
            res = res +'-'*add
        
        if len(res) != extend*2+1:
            #print position, extend, res,  len(seq)
            raise Exception('bed', 'bed')  
        return res
    
    
    
    def dump(self):
        res = str(self.seq)
        return res
    

class Prot():
    def __init__(self,prot_id='', seq = ''):
        self.prot_id = prot_id
        self.seq = seq
        self.features = {
        'Deamidated':[],
        'HexNAc':[],
        'Prediction':[],
        'SignalPeptide':[],
        'ExtraCellular':[],
        'TransMembrane':[],      
        'Domains':[],
        'Coverage':[],
        'Coverage2':[],
        'GPI':[], 
        }
        
    def load_feature(self,name='',start=0,end=0,score=0,desc=''):
        if name in self.features:
            temp = self.features[name]
            temp.append([start,end,score,desc])
            self.features[name] = temp
        else:
            raise Exception('feature', 'not inplemented')
    
    def sum_feature(self, name='', start=0, end=0, score=0, desc=''):
        if name in self.features:
            temp = self.features[name]
            new_list = []
            a=0
            for item in temp:
                if item[0] == start and item[1] == end:
                    a=1
                    new_list.append([item[0],item[1],item[2]+score,desc])
                else:
                    new_list.append(item)
            if a == 0:
                new_list.append([start,end,score,desc])
            
            self.features[name] = new_list
        else:
            raise Exception('feature', 'not inplemented')        
    
    
    def dump(self, extend=6):
        res = self.features
        #print(res)
        res['id']=self.prot_id
        res['seq']=self.seq.dump()
        '''
        res['Deamidated']=self.features['Deamidated']
        res['HexNAc']=self.features['HexNAc']
        res['Deamidated2']=self.features['Deamidated']
        res['HexNAc2']=self.features['HexNAc']        
        res['Prediction']=self.features['Prediction']   
        res['SignalPeptide']=self.features['SignalPeptide']
        res['ExtraCellular']=self.features['ExtraCellular']  
        res['TransMembrane']=self.features['TransMembrane']
        res['Domains']=self.features['Domains']
        res['Coverage']=self.features['Coverage']
        '''
        data_table = []
        
        for prediction in self.features['Prediction']:
            #print prediction
            site = int(prediction[0])
            score = round(float(prediction[2]),2)
            if score > 0.5:
                temp_pred = 'Complex/Paucimannose'
            else:
                temp_pred = 'Oligomannose'
            pep = self.seq.extract_seq(site, extend)
            Deamidated = 0
            HexNAc = 0
            for item in self.features['Deamidated']:
                if item[0] == site:
                    Deamidated = item[2]
            for item in self.features['HexNAc']:
                if item[0] == site:
                    HexNAc = item[2]                                      
            data_table.append( [pep,site+1,Deamidated,HexNAc,temp_pred] )                    
                       
        res['dataTable']=data_table

        print(res['GPI'])
        return res
        
    def hasPrediction(self):        
        if len(self.features['SignalPeptide'])>0 and len(self.features['ExtraCellular'])>0:
            for prediction in self.features['Prediction']:
                site = prediction[0]
                for extc in self.features['ExtraCellular']:
                    test = set(range(extc[0],extc[1]+1))
                    if site in test:
                        return True            
            return False     
        else:
            return False 
            
    def hasExperimentalData(self):
        if len(self.features['HexNAc'])>0 or len(self.features['Deamidated'])>0 :            
            return True
        else:
            return False
            
    
    def __str__(self):
        temp_string = '\n-------------------\n'
        temp_string +='its class prot\n'
        temp_string+=self.prot_id+'\n'
        temp_string+= '-------------------\n'
        return temp_string

'''
#not used
def update_coverage(my_dict, prot):
    new_coverage_1 =[]
    new_coverage_2 =[]
    for a,b in zip(my_dict[prot].features['Coverage'], my_dict[prot].features['Coverage2']) :
        print ( a,b, a['y']+b['y']) 
        if a['y']+b['y'] == 0:
           new_coverage_1.append(a) 
           new_coverage_2.append(b) 
        else:
            new_a = {}
            new_a['y']=round(float(a['y'])/(a['y']+b['y']),2)
            new_a['x']=a['x']

            new_b = {}
            new_b['y']=round(float(b['y'])/(a['y']+b['y']),2)
            new_b['x']=b['x']

            new_coverage_1.append(new_a) 
            new_coverage_2.append(new_b)
    my_dict[prot].features['Coverage']=new_coverage_1
    my_dict[prot].features['Coverage2']=new_coverage_2
'''

def update_coverage(my_dict, prot):
    de = my_dict[prot].features['Deamidated']
    hx = my_dict[prot].features['HexNAc']

    de_dict = {}
    hx_dict = {}
    for site_1 in de:
        de_dict[site_1[0]]=site_1
    for site_2 in hx:
        hx_dict[site_2[0]]=site_2

    sum_site = {}
    for site_1 in de:
        for site_2 in hx:
            if site_1[0] == site_2[0]:
                #print(site_1,site_2)
                temp = []
                for n in site_1:
                    temp.append(n)
                temp[2] = site_1[2]+site_2[2]
                sum_site[temp[0]]=temp
                #print('found', site_1, site_2, temp)
    
    new_coverage_hx = []
    new_coverage_de = []

    for n in my_dict[prot].features['Coverage']:
        pos = n['x']
        #print(pos)
        if pos in sum_site:
            new_coverage_hx.append({'x':pos+1,'y':  round(float(hx_dict[pos][2]) / sum_site[pos][2], 2)   })
            new_coverage_de.append({'x':pos+1,'y':  round(float(de_dict[pos][2]) / sum_site[pos][2] ,2)   })
        else:
            if pos in de_dict:
                new_coverage_de.append({'x':pos+1, 'y':  1.0  })
            else:
                new_coverage_de.append({'x':pos+1, 'y':  0.0  })   

            if pos in hx_dict:
                new_coverage_hx.append({'x':pos+1, 'y':  1.0  })
            else:
                new_coverage_hx.append({'x':pos+1, 'y':  0.0  })

    print(len(my_dict[prot].features['Coverage']),len(new_coverage_hx),len(new_coverage_de))
    
    #this to fix a viz issue, when ther is only one site with y=1 the viz brakes
    to_mod = 0
    for n in new_coverage_de:
        if n['y'] > 0:
            to_mod=1
            break
    if to_mod == 1:
        new_coverage_de = [{'y': 0.1, 'x': -100},{'y': 0, 'x': 1}] + new_coverage_de
        new_coverage_hx = [{'y': 0.1, 'x': -100},{'y': 0, 'x': 1}] + new_coverage_hx
    
    to_mod = 0
    for n in new_coverage_hx:
        if n['y'] > 0:
            to_mod=1
            break
    if to_mod == 1:
        new_coverage_hx = [{'y': 0.1, 'x': -100},{'y': 0, 'x': 1}] + new_coverage_hx
        new_coverage_de = [{'y': 0.1, 'x': -100},{'y': 0, 'x': 1}] + new_coverage_de

    #new_coverage_hx = [{'y': 0, 'x': 0},{'y': 0, 'x': 1}] + new_coverage_hx

    my_dict[prot].features['Coverage']=new_coverage_de
    my_dict[prot].features['Coverage2']=new_coverage_hx


def run_update_coverage():
    file = open("data/Proteins.pkl",'r')
    my_dict = pickle.load(file)

    
    for prot in my_dict:
        update_coverage(my_dict, prot)
    
    output = open('data/Proteins2.pkl', 'wb')
    pickle.dump(my_dict, output)
    output.close()   



def make_dict():
        dict_prot = {}
        
        handle = open("data/TriTrypDB-28_TbruceiTREU927_AnnotatedProteins.fasta", "rU")
        for record in SeqIO.parse(handle, "fasta"):
            seq = Seq(str(record.seq))
            prot = Prot(record.id, seq)
            dict_prot[record.id]=prot
        
        
        for item in  open('data/phobius_prediction.txt').read().split('//'):
            temp_res = {}
            for l in item.split('\n'):
                if l.startswith('ID'):
                    prot_id = strip(l.split()[-1])
                    temp_res['id']=prot_id
                if 'NON CYTOPLASMIC.' in l:
                    #print('NON CYTOPLASMIC.',l)
                    temp_data = l.split()
                    start = int(temp_data[2])
                    end = int(temp_data[3])
                    score = 0
                    temp_prot = dict_prot[prot_id]
                    temp_prot.load_feature('ExtraCellular',start,end,score)
                    dict_prot[prot_id] = temp_prot
                
                if l.startswith('FT   TRANSMEM'):
                    #print('TRANSMEM',l)
                    temp_data = l.split()
                    start = int(temp_data[2])
                    end = int(temp_data[3])
                    score = 0
                    temp_prot = dict_prot[prot_id]
                    temp_prot.load_feature('TransMembrane',start,end,score)
                    dict_prot[prot_id] = temp_prot            
                
                if l.startswith('FT   SIGNAL'):
                    temp_data = l.split()
                    start = int(temp_data[2])
                    end = int(temp_data[3])
                    score = 0
                    temp_prot = dict_prot[prot_id]
                    temp_prot.load_feature('SignalPeptide',start,end,score)
                    dict_prot[prot_id] = temp_prot          
                    
                    
        for l in  open('data/all_peptide_predicted.csv').readlines()[1:]:
            item_list = l.split(',')
            prot_id = item_list[0].split('_')[0]
            site = int(item_list[0].split('_')[1])
            score = float(strip(item_list[2])) 
            temp_prot = dict_prot[prot_id]
            temp_prot.load_feature('Prediction',site,site,round(score,2))
            dict_prot[prot_id] = temp_prot 

        for l in  open('data/Table_gpi.csv').readlines()[1:]:
            item_list = l.split(',')
            prot_id = item_list[1]
            site = int(item_list[2])
            score = float(strip(item_list[3]))
            prot_len = int(item_list[5])
            if prot_id in dict_prot:
                temp_prot = dict_prot[prot_id]
                temp_prot.load_feature('GPI', site, prot_len, round(score,2))
                dict_prot[prot_id] = temp_prot




        
        for l in  open('data/F235663_4.txt').readlines()[1:]:
            item_list = l.split('\t') 
            Deamidated = int(item_list[1])
            HexNAc = int(item_list[2])
            prots =  strip(item_list[-1]).split(',')   
            for x in prots:
                prot_id = x.split('_')[0]
                site = int(x.split('_')[1])-1
                temp_prot = dict_prot[prot_id]
                
                if HexNAc > 0:
                    temp_prot.load_feature('HexNAc',site,site,HexNAc)
                    dict_prot[prot_id] = temp_prot                   
                if Deamidated > 0:
                    temp_prot.load_feature('Deamidated',site,site,Deamidated)
                    dict_prot[prot_id] = temp_prot
        
        for l in  open('data/F238646_4.txt').readlines()[1:]:
            item_list = l.split('\t') 
            Deamidated = int(item_list[1])
            HexNAc = int(item_list[2])
            prots =  strip(item_list[-1]).split(',')   
            for x in prots:
                prot_id = x.split('_')[0]
                site = int(x.split('_')[1])-1
                temp_prot = dict_prot[prot_id]
                if HexNAc > 0:
                    temp_prot.sum_feature('HexNAc',site,site,HexNAc)
                    dict_prot[prot_id] = temp_prot                   
                if Deamidated > 0:
                    temp_prot.sum_feature('Deamidated',site,site,Deamidated)
                    dict_prot[prot_id] = temp_prot    
        
        
        
        for l in open('data/batch_all_hitdata.txt'):
            item_list = l.split('\t') 
            prot_id =  strip(item_list[0].split('>')[1])
            start = int(item_list[3])
            end = int(item_list[4]) 
            if float(item_list[5]) == 0.0:
                score = 150
            else:
                score = round(-np.log10(float(item_list[5])),2)
                
            description = item_list[7] +' '+item_list[8]
            temp_prot = dict_prot[prot_id]
            temp_prot.load_feature('Domains',start,end,score,description)
            dict_prot[prot_id] = temp_prot 

        
        
        temp_coverage = []
        temp_prot = 'Tb927.5.1810'
        table = open('data/F235663_1.txt').readlines()[1:]
        for line in table:
            row = line.split(',')
            #print row
            prot = row[2]
            start = int(row[-15])
            end = int(row[-14])
            #print prot,start,end

            if prot == temp_prot:
                temp_coverage+=range(start,end+1)
            else:
                sequence_length= len(str(dict_prot[temp_prot].seq.dump()))
                sequence_range = range(1,sequence_length+1)
                temp_list = []
                for pos in sequence_range:
                    temp_list.append({'x':pos,'y':temp_coverage.count(pos)})


                starting = dict_prot[temp_prot]
                features = starting.features

                features['Coverage'] = temp_list
                starting.features = features
                dict_prot[temp_prot]=starting
                
                temp_coverage = []
                temp_prot = prot
                temp_coverage+=range(start,end+1)
            
        
        temp_coverage = []
        temp_prot = 'Tb927.5.1810'
        table = open('data/F238646_1.txt').readlines()[1:]
        for line in table:
            row = line.split(',')
            #print row
            prot = row[3]
            start = int(row[-18])
            end = int(row[-17])
            #print prot,start,end

            if prot == temp_prot:
                temp_coverage+=range(start,end+1)
            else:
                sequence_length= len(str(dict_prot[temp_prot].seq.dump()))
                sequence_range = range(1,sequence_length+1)
                temp_list = []
                for pos in sequence_range:
                    temp_list.append({'x':pos,'y':temp_coverage.count(pos)})


                starting = dict_prot[temp_prot]
                features = starting.features

                features['Coverage2'] = temp_list
                starting.features = features
                dict_prot[temp_prot]=starting
                
                temp_coverage = []
                temp_prot = prot
                temp_coverage+=range(start,end+1)            
                

        output = open('data/Proteins.pkl', 'wb')
        pickle.dump(dict_prot, output)
        output.close()

if __name__ == '__main__':
    make_dict()
    run_update_coverage()

    file = open("data/Proteins2.pkl",'r')
    my_dict = pickle.load(file)
    print(my_dict['Tb927.5.4570'].features['Coverage'])
  




    #run_update_coverage()

    
    
    
    #print('---------')
    #for s in hx:
    #    print(s)
    



    #my_dict2 = update_coverage(my_dict, 'Tb927.5.4570')
    #print(my_dict['Tb927.5.4570'])
    #print(my_dict['Tb927.5.4570'].features['Coverage'][134:149])
    #print(my_dict['Tb927.5.4570'].features['Coverage2'][134:149])

    #update_coverage(my_dict, 'Tb927.5.4570')
    #print(my_dict['Tb927.5.4570'].features['Coverage'][134:149])
    #print(my_dict['Tb927.5.4570'].features['Coverage2'][134:149])

    #output = open('data/Proteins2.pkl', 'wb')
    #pickle.dump(my_dict, output)
    #output.close()


    #print(my_dict2['Tb927.5.4570'])
    #print(len(my_dict2['Tb927.5.4570'].features['Coverage']))
    #print(len(my_dict2['Tb927.5.4570'].features['Coverage2']))
        






    
        #'''
           
    #print  dict_prot['Tb927.5.1830'].dump()
    #print  dict_prot['Tb927.5.1830'].features['Coverage']




    
    



