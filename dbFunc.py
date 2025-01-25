import pandas
import sqlite3

def bringTable(route: str, is_excel=False):
        if(is_excel):
            table = pandas.read_excel(route, sheet_name=1, header=1)
        else:
            table = pandas.read_table(route,sep='\t',header=1)    
        return(table)

#Queries for the database

def nameSearch(name: str, signif: bool = False, signifThr: float = 0.05):
    if signif:
        sQuery = """SELECT expression.ORF, expression.short_name, expression.logFC, expression.FDR, pannzerAn.description, goAn.description, EggNogAn.long_desc 
        FROM ((expression LEFT JOIN goAn ON expression.ORF = goAn.ORF) LEFT JOIN pannzerAn ON expression.ORF = pannzerAn.ORF) LEFT JOIN EggNogAn ON expression.ORF = EggNogAn.ORF 
        WHERE expression.short_name LIKE ?
        GROUP BY expression.logFC, expression.FDR
        HAVING expression.FDR < ?
        ORDER BY expression.logFC DESC;"""
        extraParam = [name.join(['%','%']), signifThr]
    else:
        sQuery = """SELECT expression.ORF, expression.short_name, expression.logFC, expression.FDR, pannzerAn.description, goAn.description, EggNogAn.long_desc 
        FROM ((expression LEFT JOIN goAn ON expression.ORF = goAn.ORF) LEFT JOIN pannzerAn ON expression.ORF = pannzerAn.ORF) LEFT JOIN EggNogAn ON expression.ORF = EggNogAn.ORF 
        WHERE expression.short_name LIKE ?
        GROUP BY expression.logFC, expression.FDR
        ORDER BY expression.logFC DESC;"""
        extraParam = [name.join(['%','%'])]
    
    with sqlite3.connect('transcription.db') as con:
        resultsTab = pandas.read_sql_query(sQuery,con = con, params = extraParam)
        
    return resultsTab

def ecSearch(ec: str, signif: bool = False, signifThr: float = 0.05):
    if signif:
        sQuery = """SELECT expression.ORF, expression.short_name, expression.logFC, expression.FDR, pannzerAn.description, goAn.description, EggNogAn.long_desc 
        FROM ((expression LEFT JOIN goAn ON expression.ORF = goAn.ORF) LEFT JOIN pannzerAn ON expression.ORF = pannzerAn.ORF) LEFT JOIN EggNogAn ON expression.ORF = EggNogAn.ORF 
        WHERE EggNogAn.EC LIKE ?
        GROUP BY expression.logFC, expression.FDR
        HAVING expression.FDR < ?
        ORDER BY expression.logFC DESC;"""
        extraParam = [ec.join(['%','%']), signifThr]
    else:
        sQuery = """SELECT expression.ORF, expression.short_name, expression.logFC, expression.FDR, pannzerAn.description, goAn.description, EggNogAn.long_desc 
        FROM ((expression LEFT JOIN goAn ON expression.ORF = goAn.ORF) LEFT JOIN pannzerAn ON expression.ORF = pannzerAn.ORF) LEFT JOIN EggNogAn ON expression.ORF = EggNogAn.ORF 
        WHERE EggNogAn.EC LIKE ?
        GROUP BY expression.logFC, expression.FDR
        ORDER BY expression.logFC DESC;"""
        extraParam = [ec.join(['%','%'])]
    
    with sqlite3.connect('transcription.db') as con:
        resultsTab = pandas.read_sql_query(sQuery, con = con, params = extraParam)
    
    return resultsTab

def pfamSearch(pfam: str, signif: bool = False, signifThr: float = 0.05):
    if signif:
        sQuery = """SELECT expression.ORF, expression.short_name, expression.logFC, expression.FDR, pannzerAn.description, goAn.description, EggNogAn.long_desc 
        FROM ((expression LEFT JOIN goAn ON expression.ORF = goAn.ORF) LEFT JOIN pannzerAn ON expression.ORF = pannzerAn.ORF) LEFT JOIN EggNogAn ON expression.ORF = EggNogAn.ORF 
        WHERE EggNogAn.pfam LIKE ?
        GROUP BY expression.logFC, expression.FDR
        HAVING expression.FDR < ?
        ORDER BY expression.logFC DESC;"""
        extraParam = [pfam.join(['%','%']), signifThr]
    else:
        sQuery = """SELECT expression.ORF, expression.short_name, expression.logFC, expression.FDR, pannzerAn.description, goAn.description, EggNogAn.long_desc 
        FROM ((expression LEFT JOIN goAn ON expression.ORF = goAn.ORF) LEFT JOIN pannzerAn ON expression.ORF = pannzerAn.ORF) LEFT JOIN EggNogAn ON expression.ORF = EggNogAn.ORF 
        WHERE EggNogAn.pfam LIKE ?
        GROUP BY expression.logFC, expression.FDR
        ORDER BY expression.logFC DESC;"""
        extraParam = [pfam.join(['%','%'])]
    
    with sqlite3.connect('transcription.db') as con:
        resultsTab = pandas.read_sql_query(sQuery, con = con, params = extraParam)
    
    return resultsTab

def keggReactSearch(kegg: str, signif: bool = False, signifThr: float = 0.05):
    if signif:
        sQuery = """SELECT expression.ORF, expression.short_name, expression.logFC, expression.FDR, pannzerAn.description, goAn.description, EggNogAn.long_desc 
        FROM ((expression LEFT JOIN goAn ON expression.ORF = goAn.ORF) LEFT JOIN pannzerAn ON expression.ORF = pannzerAn.ORF) LEFT JOIN EggNogAn ON expression.ORF = EggNogAn.ORF 
        WHERE EggNogAn.kegg_reaction LIKE ?
        GROUP BY expression.logFC, expression.FDR
        HAVING expression.FDR < ?
        ORDER BY expression.logFC DESC;"""
        extraParam = [pfam.join(['%','%']), signifThr]
    else:
        sQuery = """SELECT expression.ORF, expression.short_name, expression.logFC, expression.FDR, pannzerAn.description, goAn.description, EggNogAn.long_desc 
        FROM ((expression LEFT JOIN goAn ON expression.ORF = goAn.ORF) LEFT JOIN pannzerAn ON expression.ORF = pannzerAn.ORF) LEFT JOIN EggNogAn ON expression.ORF = EggNogAn.ORF 
        WHERE EggNogAn.kegg_reaction LIKE ?
        GROUP BY expression.logFC, expression.FDR
        ORDER BY expression.logFC DESC;"""
        extraParam = [pfam.join(['%','%'])]

    with sqlite3.connect('transcription.db') as con:
        resultsTab = pandas.read_sql_query(sQuery, con = con, params = extraParam)
    
    return resultsTab

def descriptionSearch(desc: str, signif: bool = False, signifThr: float = 0.05):
    if signif:
        sQuery = """SELECT expression.ORF, expression.short_name, expression.logFC, expression.FDR, pannzerAn.description, goAn.description, EggNogAn.long_desc 
        FROM ((expression LEFT JOIN goAn ON expression.ORF = goAn.ORF) LEFT JOIN pannzerAn ON expression.ORF = pannzerAn.ORF) LEFT JOIN EggNogAn ON expression.ORF = EggNogAn.ORF 
        WHERE pannzerAn LIKE ? OR goAn.description LIKE ? OR EggNogAn.long_desc LIKE ?
        GROUP BY expression.logFC, expression.FDR
        HAVING expression.FDR < ?
        ORDER BY expression.logFC DESC;"""
        extraParam = [pfam.join(['%','%']), pfam.join(['%','%']), pfam.join(['%','%']), signifThr]
    else:
        sQuery = """SELECT expression.ORF, expression.short_name, expression.logFC, expression.FDR, pannzerAn.description, goAn.description, EggNogAn.long_desc 
        FROM ((expression LEFT JOIN goAn ON expression.ORF = goAn.ORF) LEFT JOIN pannzerAn ON expression.ORF = pannzerAn.ORF) LEFT JOIN EggNogAn ON expression.ORF = EggNogAn.ORF 
        WHERE EggNogAn.kegg_reaction LIKE ?
        GROUP BY expression.logFC, expression.FDR
        ORDER BY expression.logFC DESC;"""
        extraParam = [pfam.join(['%','%']), pfam.join(['%','%']), pfam.join(['%','%'])]
    
    return resultsTab

if __name__ == '__main__':
    
    conexion = sqlite3.connect('transcription.db')
    cursor = conexion.cursor()
    #Table creation
    cursor.execute("""CREATE TABLE IF NOT EXISTS expression(
    ORF VARCHAR(15) PRIMARY KEY,
    short_name VARCHAR,
    logFC DOUBLE,
    logCPM DOUBLE,
    PValue DOUBLE,
    FDR DOUBLE);""")    
    
    cursor.execute("""CREATE TABLE IF NOT EXISTS goAn(
    id INT PRIMARY KEY,
    ORF VARCHAR(15),
    ontology VARCHAR,
    goID INT,
    description VARCHAR,
    argotScore DOUBLE,
    argotPPV DOUBLE);""")

    cursor.execute("""CREATE TABLE IF NOT EXISTS pannzerAn(
    id INT PRIMARY KEY,
    ORF VARCHAR(15),
    type VARCHAR(10),
    score DOUBLE,
    PPV DOUBLE,
    description VARCHAR);""")

    cursor.execute("""CREATE TABLE IF NOT EXISTS EggNogAn(
    id INT PRIMARY KEY,
    ORF VARCHAR(15),
    long_desc VARCHAR,
    evalue DOUBLE,
    eggNOGOG VARCHAR,
    GO VARCHAR,
    EC VARCHAR,
    kegg_ko VARCHAR,
    kegg_path VARCHAR,
    kegg_module VARCHAR,
    kegg_reaction VARCHAR,
    kegg_rclass VARCHAR,
    pfam VARCHAR,
    eggnog_score DOUBLE,
    brite VARCHAR);""")
    
    #Filling the newly created tables
    expresion = bringTable('expresion_diferencial_cds.txt', False)
    expresion.to_sql('expression', conexion, if_exists='fail', index=False)
    pannzer = bringTable('Pannzer_para_db.xlsx',True)
    pannzer.to_sql('pannzerAn', if_exists="fail", index=True, index_label='id')
    go = bringTable('panzerGO_para_db.xlsx',True)
    go.to_sql('goAn', conexion, if_exists='fail', index=True, index_label='id')
    eggnog = bringTable('Pannzer-Eggnog annotation_para_db.xlsx', True)
    eggnog.to_sql('EggNogAn', conexion, if_exists='fail', index=True, index_label='id')
    conexion.commit()
    conexion.close()