
import re

class Seq():
    def __init__(self, seq):
        self.seq = seq
    
    def find(self, query):
        res = []
        for item in re.finditer(query, self.seq):
            res.append( item.start() )
        return res
    
    def extract_seq(self, position, extend):
        
        if position < extend and position >= len(self.seq)-extend:
            add_front = extend-position
            add_back = extend - (len(self.seq) - position-1)
            res = add_front*'-'+self.seq+'-'*add_back
        
        elif position >= extend and position < len(self.seq)-extend:
            res = self.seq[position-extend:position+extend+1]
        
        elif position < extend:
            res = self.seq[:position+extend+1]
            add = extend-position
            res = '-'*add + res
        
        elif position >= len(self.seq)-extend:
            res = self.seq[position-extend:]
            add = extend - (len(self.seq) - position-1)
            res = res +'-'*add
        
        if len(res) != extend*2+1:
            print position, extend, res,  len(self.seq)
            raise Exception('spam', 'eggs')
        #print res   
        return res
        
    def format_prediction(self, extend, pattern, record_id):
        res = []
        positions = self.find(pattern)
        for position in positions:
            pep = self.extract_seq(position, extend)
            code = '0'*extend +'1'+'0'*extend
            temp = '>'+record_id+'_'+str(position)+'_'+pep+'\n'+pep+'\n'+code
            res.append(temp) 
        return res
        
        
        

        
if __name__ == '__main__':
    
    seq ='''MPSNQEARLFLAVLVLAQVLPILVDSAAEKGFKQAFWQPLCQVSEELDDQPKGALFTLQA
    AASKIQKMRDAALRASIYAEINHGTNRAKAAVIVANHYAMKADSGLEALKQTLSSQEVTA
    TATASYLKGRIDEYLNLLLQTKESGTSGCMMDTSGTNTVTKAGGTIGGVPCKLQLSPIQP
    KRPAATYLGKAGYVGLTRQADAANNFHDNDAECRLASGHNTNGLGKSGQLSAAVTMAAGY
    VTVANSQTAVTVQALDALQEASGAAHQPWIDAWKAKKALTGAETAEFRNETAGIAGKTGV
    TKLVEEALLKKKDSEASEIQTELKKYFSGHENEQWTAIEKLISEQPVAQNLVGDNQPTKL
    GELEGNAKLTTILAYYRMETAGKFEVLTQKHKPAESQQQAAETEGSCNKKDQNECKSPCK
    WHNDAENKKCTLDKEEAKKVADETAKDGKTGNTNTTGSSNSFVISKTPLWLAVLLF
    '''
    
    record_id = 'pippo'
    extend = 10
    pattern = 'N.[ST]'
    
    seq = seq.replace('\n','')
    seq = Seq(seq)
    positions = seq.find(pattern)
    for position in positions:
        pep = seq.extract_seq(position,extend)
        print pep
        code = '0'*extend +'1'+'0'*extend
        temp = '>'+record_id+'_'+str(position)+'_'+pep+'\n'+pep+'\n'+code+'\n'
        print temp
    
    res = seq.format_prediction(extend, pattern, record_id)
    print '\n'.join(res)
        

    
    
    
    
    
    
    