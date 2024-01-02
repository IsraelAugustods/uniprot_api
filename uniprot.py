# Imports 
import pandas as pd
import requests
import difflib
    
# Definiçao da lista de proteinas e ID dos organismos que queremos.
data = pd.read_excel('datasets/T1xT4_14DAA - RNA-seq Cana_Tupaciguara_2022.xlsx')
protein_names = data['protein_name'].head(5).to_list()



print(protein_names)

# Arabidopsis -> Arroz -> Milho -> Sorgo
organism_id_c3 = ['3702','4530','4577','4558']

# Sorgo -> Milho -> Setaria
organism_id_c4 = ['4558', '4577', '4554']

blacklist = ['type=', 'hypothetical', 'Hypothetical']

# Definição de alguns tipos de arquivo vazios, que vamos precisar pro codigo rodar.
df3 = {'Organism': [], 'Protein names': [], 'Gene Ontology (biological process)': [], 'Gene Ontology (cellular component)':[], 'Gene Ontology (molecular function)':[],'Function [CC]': []}

df4 = {'Organism': [], 'Protein names': [], 'Gene Ontology (biological process)': [], 'Gene Ontology (cellular component)':[], 'Gene Ontology (molecular function)':[], 'Function [CC]': []}

# Função de similaridade
def similaridade(alvo, nomes_prot_tabela):
        global pontuacoes, index_mais_similar
        pontuacoes =[]
        for item in nomes_prot_tabela:
            item = item.replace("[...]", "").strip()
            pontuacao = difflib.SequenceMatcher(None, alvo, item).ratio()*100
            pontuacoes.append(pontuacao)
        index_mais_similar = pontuacoes.index(max(pontuacoes))
        return index_mais_similar

# Função para organismos c3
def get_c3(protein_names, organism_id_c3):

    
#leitura das proteinas da lista
    global prot_lower  
    for protein in protein_names:
        invalid_prot = False    
        prot_lower = protein.lower()

#black list 
        for item in blacklist:
            if item in prot_lower:
                invalid_prot = True
                
        if invalid_prot == True:
            continue            
        
        global r 
        for x in range(len(organism_id_c3)):

            url = f'https://rest.uniprot.org//uniprotkb/search?query=(cc_function:*){protein} AND taxonomy_id={int(organism_id_c3[x])}&fields=organism_name&format=tsv' 
            r = requests.get(url).text
            
            if len(r) == 9:
                continue

            url1 = f'https://rest.uniprot.org//uniprotkb/search?query=(cc_function:*){protein} AND taxonomy_id={int(organism_id_c3[x])}&fields=protein_name&format=tsv' 
            r1 = requests.get(url1).text

            url2 = f'https://rest.uniprot.org//uniprotkb/search?query=(cc_function:*){protein} AND taxonomy_id={int(organism_id_c3[x])}&fields=go_p&format=tsv' 
            r2 = requests.get(url2).text

            url3 = f'https://rest.uniprot.org//uniprotkb/search?query=(cc_function:*){protein} AND taxonomy_id={int(organism_id_c3[x])}&fields=go_c&format=tsv' 
            r3 = requests.get(url3).text

            url4 = f'https://rest.uniprot.org//uniprotkb/search?query=(cc_function:*){protein} AND taxonomy_id={int(organism_id_c3[x])}&fields=go_f&format=tsv' 
            r4 = requests.get(url4).text

            url5 = f'https://rest.uniprot.org//uniprotkb/search?query=(cc_function:*){protein} AND taxonomy_id={int(organism_id_c3[x])}&fields=cc_function&format=tsv' 
            r5 = requests.get(url5).text


            a = r, r1, r2, r3, r4, r5
                
            org =    a[0].split(sep = '\n')
            prot_n = a[1].split(sep = '\n')
            go_p =   a[2].split(sep = '\n')
            go_c =   a[3].split(sep = '\n')
            go_f =   a[4].split(sep = '\n')
            func =   a[5].split(sep = '\n')

            
            df = pd.DataFrame( {org[0]:org, prot_n[0]:prot_n, go_p[0]:go_p, go_c[0]:go_c, go_f[0]:go_f, func[0]:func})
            df = df.drop(0)
            df = df.reset_index()

            df2 = df.iloc[similaridade(protein, df['Protein names'])]

            if pontuacoes[index_mais_similar] > 50 and len(df2['Function [CC]']) != 0:

                df3['Organism'].append(df2['Organism'])
                df3['Protein names'].append(df2['Protein names'])
                df3['Gene Ontology (biological process)'].append(df2['Gene Ontology (biological process)'])
                df3['Gene Ontology (cellular component)'].append(df2['Gene Ontology (cellular component)'])
                df3['Gene Ontology (molecular function)'].append(df2['Gene Ontology (molecular function)'])
                df3['Function [CC]'].append(df2['Function [CC]'])
                break 
  
        if len(r) == 9:
                df3['Protein names'].append(prot_lower)
                df3['Organism'].append('NA')
                df3['Gene Ontology (biological process)'].append('NA')
                df3['Gene Ontology (cellular component)'].append('NA')
                df3['Function [CC]'].append('NA')
                df3['Gene Ontology (molecular function)'].append('NA')
        
        
    
    
    global df_c3, df_final
    df_c3 = pd.DataFrame(df3)
    
    df_c3 = df_c3.rename(columns={'Organism': 'organism_c3', 'Gene Ontology (molecular function)': 'MF_c3', 'Gene Ontology (cellular component)': 'CC_c3', 'Gene Ontology (biological process)': 'BP_c3'})
    
    
    df_final = data.copy()
    df_final['organism_c3'] = df_c3['organism_c3']
    df_final['function_c3'] = df_c3['Function [CC]']
    df_final['MF_c3'] = df_c3['MF_c3']
    df_final['CC_c3'] = df_c3['CC_c3']
    df_final['BP_c3'] = df_c3['BP_c3']

    return  df_final.to_excel('t1.xlsx')
