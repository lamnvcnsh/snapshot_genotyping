import pickle
table_file = open('resource/tables.pdata', 'rb')

tables = pickle.load(table_file)

maker_table = tables['marker_table']

for name in tables:
    file_name=name+".xlsx"
    tables[name].to_excel("resource/"+file_name, index=False)
