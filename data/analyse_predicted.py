from scipy import linspace, polyval, polyfit, sqrt, stats, randn
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use('ggplot')
import pprint
import numpy as np
from sklearn.metrics import r2_score
from scipy.stats import skew, boxcox
from sklearn.linear_model import LinearRegression
from scipy import polyfit
from scipy import stats

all_prediction = pd.DataFrame.from_csv('all_peptide_predicted.csv')
print all_prediction.shape
all_prediction = all_prediction.drop_duplicates('window_seq_win-10')


fig, ax = plt.subplots()
all_prediction['rf'].plot(kind='kde', ax = ax)
plt.show()


fig, ax = plt.subplots()
all_prediction['rf'].plot(kind='box', ax = ax)
plt.show()

median= all_prediction['rf'].median()
std = all_prediction['rf'].std()
limit_up = round(median+std,2)
limit_down = round(median-std,2)
print 'limit up',limit_up
print 'limit down',limit_down


up = all_prediction[all_prediction['rf']>limit_up]
down = all_prediction[all_prediction['rf']<limit_up]

up['window_seq_win-10'].to_csv('pep_deamidated_ml.txt', index=False)
down['window_seq_win-10'].to_csv('pep_hexnac_ml.txt', index=False,)




start_data = pd.DataFrame.from_csv('F235663_4.txt',sep='\t')
start_data['filter']=[1  if '*' in n or '-' in n else 0 for n in start_data.index.values ]
start_data = start_data[start_data['filter']==0]
del start_data['filter']
test = list(all_prediction['window_seq_win-10'])
test_trim = [n[3:-3] for n in test]

#print test_trim.index('VRRGSGNHTIIYF')

start_data['extended_pep'] = [  test[test_trim.index(pep)] for pep in start_data.index.values]
print start_data.head()


new_data = pd.DataFrame.from_csv('F238646_4.txt',sep='\t')
new_data['filter']=[1  if '*' in n or '-' in n else 0 for n in new_data.index.values ]
new_data = new_data[new_data['filter']==0]
del new_data['filter']

merge = start_data.join(new_data, how='inner',  rsuffix='_joao')
merge['sum_ali'] = merge['Deamidated']+merge['HexNAc']
merge['fold_deam_ali'] = np.log2((merge['Deamidated']+1)  / (merge['HexNAc']+1))


merge['sum_joao'] = merge['Deamidated_joao']+merge['HexNAc_joao']
merge['fold_deam_joao'] = np.log2((merge['Deamidated_joao']+1)  / (merge['HexNAc_joao']+1))


regr = LinearRegression()
print merge.head()
fig,ax = plt.subplots()
temp = merge[ (merge['sum_ali']>1) & (merge['sum_joao'] >1)]
temp = temp[['sum_ali','sum_joao','fold_deam_joao','fold_deam_ali','extended_pep']]
temp = temp.dropna()
temp.plot(kind='scatter', x='fold_deam_joao', y='fold_deam_ali',ax=ax )

print temp['fold_deam_joao'].shape
print temp['fold_deam_ali'].shape

x = temp['fold_deam_joao'].values
y = temp['fold_deam_ali'].values
gradient, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print "R-squared", round(r_value**2,2)
abline_values = [gradient * i + intercept for i in x]
ax.plot(x, abline_values, 'b')
temp['abline_values']=abline_values
temp[['extended_pep','abline_values']].to_csv('pep_to_predict.csv')
#(ar,br)=polyfit(temp['fold_deam_joao'],temp['fold_deam_ali'],1)
#xr=polyval([ar,br],temp['fold_deam_ali'])
#compute the mean square error
#err=sqrt(sum((xr-temp['fold_deam_ali'])**2)/n)
#print('Coefficients: \n', regr.coef_)
#ax.plot(temp['fold_deam_joao'].values, regr.predict(temp['fold_deam_joao'].values), color='blue',
         #linewidth=3)

#ax.plot(temp['fold_deam_ali'],xr,'r--')         
ax.set_xlim(-8,8)
ax.set_ylim(-8,8)
plt.show()
print temp.shape

temp['abline_values'].plot(kind='hist', bins=15, normed=1)
temp['abline_values'].plot(kind='kde')

'''
for limit in range(0,40):
    temp = merge[ (merge['sum_ali']>limit) & (merge['sum_joao'] >limit)]
    x = temp['fold_deam_joao'].values
    y = temp['fold_deam_ali'].values
    gradient, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    print limit+1, temp.shape[0], "R-squared", r_value**2
'''


'''
merge['avg_deamidated'] = (merge['Deamidated']+merge['Deamidated_joao']+1)/2
merge['avg_exnac'] = (merge['HexNAc']+merge['HexNAc_joao']+1)/2



merge['tot_fold'] = np.log2(merge['avg_deamidated']/merge['avg_exnac'])
merge['tot_fold'] .plot(kind='hist')
plt.show()

#merge['tot_fold_transform'], lam = boxcox(merge['tot_fold'])
#merge['tot_fold_transform'] .plot(kind='hist')
'''

'''
fig,ax = plt.subplots()
merge.plot(kind='scatter', x='perc_deam_joao', y='perc_deam',ax=ax )

fig,ax = plt.subplots()
temp = merge[ (merge['sum']>2) & (merge['sum_joao'] >2)]
temp.plot(kind='scatter', x='perc_deam_joao', y='perc_deam',ax=ax )
plt.show()
print temp.shape
fig,ax = plt.subplots()
temp = merge[ (merge['sum']>2) & (merge['sum_joao'] >2)]
temp.plot(kind='scatter', x='fold_deam_joao', y='fold_deam',ax=ax )
plt.show()
print temp.shape
'''




