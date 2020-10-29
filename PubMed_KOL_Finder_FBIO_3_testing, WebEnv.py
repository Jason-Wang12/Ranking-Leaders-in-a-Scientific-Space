from Bio import Entrez
from Bio import Medline
from tqdm import tqdm
import pandas as pd
pd.set_option('display.max_colwidth', -1)
import numpy as np

import re
#Graphing modules
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter


# Change this email to your email address
Entrez.email = "jwang@fortressbiotech.com"

people = input('please enter a search: ')

def lenth_results(x):
    Entrez.email=Entrez.email
    keyword = x
    handle = Entrez.esearch(db ='pubmed',
                            retmax=10,
                            retmode ='text',
                            term = keyword)
    results= Entrez.read(handle)
    len_results = results['Count']
    print('Total number of publications that contain the term {}: {}'.format(keyword, results['Count']))    
    return int(len_results)

def search(x):
    Entrez.email=Entrez.email
    pubs = []
    keyword = x
    handle = Entrez.esearch(db ='pubmed',
                            retmax=int(len_results),
                            retmode ='text',
                            mindate = 2005,
                            maxdate = 2020,
                            usehistory='y',
                            idtype = 'acc',
                            term = keyword)
    results= Entrez.read(handle)
    pubs.append(results)
    print('Total number of publications that contain the term {}: {}'.format(keyword, results['Count']))    
    return pubs

def history(x):
    webenv = []
    for x in pubs:
        for i in [x['WebEnv']]:
            webenv.append(i)
    return webenv
def key(x):    
    query_key = []
    for x in pubs:
        for i in [x['QueryKey']]:
            query_key.append(i)
    return query_key

def id_list(x):
    id_list = []
    for x in pubs:
        for i in [x['IdList']]:
            id_list.append(i)
    id_list = id_list[0] #Right now it's a numpy array, turn it into a 1d list
    return id_list

def count(x):    
    count = []
    for x in pubs:
        for i in [x['Count']]:
            count.append(int(i))
    return count

def search_result(history, key, count):
    search_results = {'History':[], 'QueryKey':[], 'Count':[]}
    search_results['History'].append(webenv)
    search_results['QueryKey'].append(query_key)
    search_results['Count'].append(count)
    #Turn to list to more easily iterate over
    search_results = [search_results]
    search_results2 = pd.DataFrame(pubs)
    lst_col = 'IdList' #Set the column you want to expand
    search_results2 = pd.DataFrame({
      col:np.repeat(search_results2[col].values, search_results2[lst_col].str.len())
      for col in search_results2.columns.drop(lst_col)}
    ).assign(**{lst_col:np.concatenate(search_results2[lst_col].values)})[search_results2.columns]
    return search_results2

def batch(x):
    ids = list(search_results2['IdList'])
    batch_size = 100
    batches = [ids[x: x+100] for x in range(0, len(ids), batch_size)]
    return batches

def fetch(x):
    record_list = []
    batch_size = 100
    for x in tqdm(batches):
        handle = Entrez.efetch(db = 'pubmed',
                               rettype = 'medline',
                               retmode = 'text',
                               id = x,
                               retmax = batch_size,
                               webenv = webenv,
                               query_key = query_key,
                               idtype = 'acc')
        records = Medline.parse(handle)
        record_list.extend(list(records))
        print('Complete')  
    return record_list
 
def citations(x):
    ids = list(search_results2['IdList'])
    batch_size = 100
    batches = [ids[x: x
                   + 100] for x in range(0, len(ids), batch_size)]
    citation_list = []
    for i in tqdm(batches):
        handle = Entrez.elink(dbfrom="pubmed",
                              db="pmc",
                              LinkName="pubmed_pmc_refs",
                              id=i,
                              retmax = batch_size)
        citations = Entrez.read(handle)
        citation_list.extend(list(citations))
        print('Complete')
    return citation_list

def publication_volume(x):
    #Count number of authors, and weigh by impact factor
    df = pd.DataFrame(record_list, columns = ['FAU', 'TA', 'AD', 'TI'])
    #This splits up the list into searate rows
    df['FAU'] = df['FAU'].apply(pd.Series) \
        .merge(df, left_index = True, right_index = True
           )
    df['TA'] = df['TA'].str.lower()

    #Create weighting for influential journals
    mask = {'nature':40,
            'lancet':40,
            'jama': 80,
            'n engl j med': 80,
            'nat commun': 5,
            r'(?s).*':1
            }

    #Basically replicates the journal names column as a new column
    df['weight'] = df['TA'].replace(mask, regex = True)
    df['weight'] = df['weight'].fillna(1)
    df = df.reindex(df.index.repeat(df.weight))
    authors_flat = [i for i in list(df['FAU'].values.tolist())]

    n = 100
    top_authors = pd.DataFrame.from_records(
        Counter(authors_flat).most_common(n), columns=["Name", "Count"]
        )


    sns.barplot(x = 'Count', y = 'Name', data = top_authors[:20], palette = 'RdBu_r')
    plt.show()
    return plt

def clean_citations(citation_list):
    citations = pd.DataFrame(citation_list, columns = ['LinkSetDb'])
    citation2 = citations['LinkSetDb'].to_list()
    citation3 = pd.DataFrame(citation2, columns = ['Link'])
    res = pd.concat([citation3,citation3.Link.apply(pd.Series)],axis=1) #separating out the ids nexted within the dictionary link
    res.columns = ['original', 'list', 'dbname', 'pub'] #need to rename columns since some of them are the same now
    res['list']= res['list'].replace(np.nan, str('a')) #change nans to strings, to differentiate from the lists of citation ids
    citations4 = res['list'].to_list()
    lst = []
    for x in citations4:
        if type(x) == list:
            lst.append(len(x))
        else:
            lst.append(0)
    return lst

def senior_author(x):
    senior_authors = data_set['FAU'].to_list()
    sr_au = []
    for x in senior_authors:
        if type(x) == list:
            sr_au.append(x[-1:])
        else:
            sr_au.append(str(x))
    sr_df = pd.DataFrame(sr_au, columns = ['Authors'])#Turn list of lists
    sr_au = [i for i in sr_df['Authors']] #..into a list
    return sr_au

def put_together():
    final_data = pd.DataFrame(sr_au, columns = ['Authors'])
    final_data['publications'] = data_set['PMID']
    final_data['citations'] = lst
    final_data['location'] =data_set['PL']
    publication_cnts = []
    for x in data_set['PMID']:
        publication_cnts.append(1)
    final_data['publication counts'] = publication_cnts
    data_set['EDAT'].replace(np.nan, '2020-01-01', inplace = True)
    final_data['Year'] = data_set['EDAT'].astype(str).str[0:4].astype(int)
    return final_data

def graph_citations(x):
    graph_cits = final_data.groupby('Authors')['citations'].sum(inplace = True).reset_index()
    graph_cits2 = final_data.groupby('Authors')['publication counts'].sum(inplace = True).reset_index()
    graph_cits['publication volume'] = graph_cits2['publication counts']
    country = final_data.groupby('Authors')['location'].first(inplace = True).reset_index()
    graph_cits['location'] = country['location']
    #Since the data is still in the same order as when you parsed it, append
    #the data about institutions/contact info so you can extract the emails
    graph_cits = graph_cits.sort_values(by = 'citations', ascending = False)
    plt.figure(figsize = (16,8))
    sns.barplot(x = 'citations', y = 'Authors', data = graph_cits[:20], palette='coolwarm')
    #graph_cits3 = graph_cits[~(graph_cits['publication volume']<=20)] #Remove authors with too few publications
    graph_cits.to_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/KOLS/{} authors, location.xlsx'.format(people))    
   
    # use the scatter function
    plt.figure(figsize=(30,25))
    ax = sns.scatterplot(x = 'publication volume', y = 'citations', data = graph_cits,
                     size = 'publication volume', sizes = (40,400), alpha = 0.5, hue = 'citations', palette = 'muted', legend = False)
    def label_point(x, y, val, ax):
        a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
        for i, point in a.iterrows():
            ax.text(point['x']+.02, point['y'], str(point['val']))

    label_point(graph_cits['publication volume'], graph_cits['citations'], graph_cits['Authors'], ax)
   
    return graph_cits

def graph_publication_volume():
     
#Number of publications on the subject over the years
    data_set.dropna(subset=['EDAT'], inplace=True)
    data_set["Year"] = (
        data_set["EDAT"].astype(str).str[0:4].astype(int)
        )
    yearly = pd.DataFrame(data_set["Year"].value_counts().reset_index())
    yearly.columns = ["Year", "Count"]
    pal = sns.color_palette("Greens_d", len(yearly))
    plt.figure(figsize = (16, 8))
    sns.barplot(x="Year", y="Count", data=yearly, palette = np.array(pal[::-1]))
    plt.title("Publications over Time")
    plt.show()

def contact_info(x):
    #GEt list of contact info for the authors
    authors = list(graph_cits['Authors']) #Only x top authors in list

    emails = pd.DataFrame(data_set, columns = ['FAU', 'AD', 'TI']).dropna() #Pull out authors and, contact info
    #This will sepearate out the names in the list into individual rows
    emails['FAU'] = emails.FAU.apply(pd.Series) \
        .merge(emails, left_index = True, right_index = True)
        #Checks to see if the name of the top 30 are in the author column
    mask = emails.FAU.apply(lambda x: any(item for item in authors if item in x))
    #Returns only rows that contain at least one of the authors
    emails = emails[mask]
    emails = emails.set_index('FAU')
    emails = emails.loc[authors] #Reorder the authors according to their rankings
    #Now reset the index so that FAU can be called again for splitting the first/last name
    emails = emails.reset_index()
    emails['AD'] = emails['AD'].apply(lambda x: re.findall(r'[\w\.-]+@[\w\.-]+', str(x))) #removes everything except the email(s)
    emails['AD'] = emails['AD'].apply(', '.join) #Turns the list of emails into strings
    #Rename the columns
    emails = emails.rename(columns = {'FAU':'Name', 'AD':'email', 'TI':'biblio', 'AB':'abstract'})
    #Split first and last name in separate columns
    splitted = emails['Name'].str.split()
    emails['First'] = splitted.str[1]
    emails['Last'] = splitted.str[0]
    #Swap email/last name location
    column_titles = ['First', 'Last', 'email', 'Name', 'biblio', 'abstract']
    emails = emails.reindex(columns = column_titles)
    #Check if last name in email address
    emails['Last'] = emails['Last'].str.lower()
    emails['email'].replace('', np.nan, inplace=True)
    emails.dropna(subset=['email'], inplace=True)
    #Drop duplicates
    emails2 = emails.drop_duplicates(subset = 'Name',
                                     keep = 'first',
                                     inplace = False)
       
    print(emails2[['Name', 'email']])

    emails2[['Name', 'email']].to_excel('C:/Users/jason/Google Drive/Python/Biopython/CSV/KOLS/ {} emails.xlsx'.format(people))
    return emails2
   
       
if __name__ == '__main__':
    len_results = lenth_results(people)
    pubs = search(people) #search
    webenv = history(pubs)#Create a list of WebEnv historys & Query Keys
    query_key = key(pubs)
    id_list = id_list(pubs)
    count = count(pubs)
    search_results2 = search_result(webenv, query_key, count) #separate out query key/webenv to iterate over for parsing and downloading results
    batches = batch(search_results2) #break into batches of 100
    record_list = fetch(batches) #fetch the records
    data_set = pd.DataFrame(record_list)   #put into a dataframe
    #data_set.to_excel('C:/Users/jason/Google Drive/Python/Biopython/CSV/KOLSdataset.xlsx')
    data_set['FAU'].replace(np.nan, 'n', regex = True, inplace = True)
    publication_volume(record_list)   #graph out the authors by publication volume
    citation_list = citations(batches) #get the citations using the same batches
    lst = clean_citations(citation_list)
    sr_au = senior_author(data_set)
    final_data = put_together()
    final_data.to_excel('C:/Users/jwang/Desktop/Python/Biopython/CSV/KOLS/{} years.xlsx'.format(people))
    graph_cits = graph_citations(final_data)
    graph_publication_volume()
    emails = contact_info(graph_cits)
   




'''
    TESTING
    '''

   
'''
#Show publications over time
top_authors = list(asdf['Authors'].loc[:5])
time_set = pd.DataFrame(record_list, columns = ['FAU', 'EDAT'])
time_set['Year'] = time_set['EDAT'].astype(str).str[0:4].astype(int)
senior_authors = time_set['FAU'].to_list()
sr_au = []
for x in senior_authors:
    if type(x) == list:
        sr_au.append(x[-1:])
    else:
        sr_au.append(str(x))
time_set['Names'] = sr_au
time_set['Names'] = [','.join(map(str, l)) for l in time_set['Names']] #convert list of author names to string
time_set[time_set['Names'].isin(top_authors)] #filter for only authors with lots of citations
time_set = time_set.drop(columns = ['EDAT', 'FAU'])
time_set['count'] = 1
time_set[:10].pivot_table('count', 'Year', 'Names', aggfunc = 'count').plot(
    kind = 'line', marker = 'o', xticks = time_set.Year.unique(), legend = False)


count = time_set.groupby(['Year', 'Names'])['count'].size().reset_index() #Group by years, and count the occurance of number of publications
plt.figure(figsize = (30,20))
df1 = pd.pivot_table(count, values='count', index='Year', columns='Names')
ax = df1.plot(kind='line', legend = True)
'''