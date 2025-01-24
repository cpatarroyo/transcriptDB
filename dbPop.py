import pandas
import sqlite3

if __name__ == '__main__':
    
    def bringTable(route: str, is_excel=False):
        if(is_excel):
            table = pandas.read_excel(route, sheet_name=1, header=1)
        else:
            table = pandas.read_table(route,sep='\t',header=1)
        
        return(table)
    
    conexion = sqlite3.connect('transcription.db')

