import os
import csv
import sys
import tempfile
import fileinput
import subprocess
from Bio import SeqIO
from pathlib import Path
import matplotlib.pyplot as plt
from subprocess import DEVNULL, PIPE
from dna_features_viewer import GraphicFeature, GraphicRecord

import pandas as pd

ZH_EXECUTABLE = Path("/home/leonid/zhunt3")
assert ZH_EXECUTABLE.is_file()

# Функция, предсказывающая участки Z-ДНК
def zhunt(query: str, windowsize: int = 6, minsize: int = 3, maxsize: int = 6):
    assert set(query).issubset({"A", "C", "G", "T", "N"})
    fd, temp = tempfile.mkstemp()
    os.close(fd)
    with open(temp, 'w') as stream:
        stream.write(query)

    subprocess.run(
        [ZH_EXECUTABLE,
         str(windowsize), str(minsize), str(maxsize), temp],
        check=True, stdout=PIPE, stderr=DEVNULL,
        input=query, encoding='ascii'
    )
    with open(temp + ".Z-SCORE", 'r') as stream:
        df = pd.read_csv(stream,
                         names=['Start', 'End', 'nu-1', 'nu-2', 'nu-3',
                                'ZH-Score', 'Sequence', 'Conformation'],
                         skiprows=1, sep='\s+')
    os.remove(temp)
    os.remove(temp + ".Z-SCORE")
    return df[['Start', 'End', 'ZH-Score', 'Sequence', 'Conformation']]

# Функция, создающая .bed файл
def make_bed(data, fileName, genome_file):
    to_write = data.copy()
    to_write.insert(0, 'chr', str(genome_file).split('.')[0])
    to_write['name'] = to_write.chr.astype(str) + '_' + to_write.Start.astype(str)
    to_write[['chr', 'Start', 'End', 'name', 'ZH-Score']].to_csv(fileName, sep='\t', header=False, index=False)

# Локальный адрес файлов имеет вид "/home/leonid/PycharmProjects/python/minor_2022/project/falciparum/1.fasta"
# species = "/yoelii/"
'''
# Этот блок заменяет все символы в .fasta файлах на N (из тех, что не A, C, G, T)
for i in range(1, 15):
    path_to_fasta = "/home/leonid/PycharmProjects/python/minor_2022/project" + species + str(i) + ".fasta"
    #current_fasta = "/home/leonid/PycharmProjects/python/minor_2022/project" + species + "current_file.fasta"
    #os.rename(path_to_fasta, current_fasta)
    start_changing = True
    for line in fileinput.input(path_to_fasta, inplace=1):
        if start_changing is True:
            start_changing = False
            new_line = line
        else:
            new_line = ''
            for i in line.rstrip():
                if i in ["A", "C", "G", "T", "N", "a", "c", "g", "t", "n"]:
                    new_line += i
                else:
                    new_line += "N"
            new_line += '\n'
        sys.stdout.write(new_line)
'''
'''
# Этот блок создаёт .bed файлы с Z-SCORE для хромосом
for i in range(1, 15):
    path_to_fasta = "/home/leonid/PycharmProjects/python/minor_2022/project" + species + str(i) + ".fasta"
    with open(path_to_fasta, 'r') as stream:
        sequences = list(SeqIO.parse(stream, format='fasta'))
    
    #for s in sequences:
    #    print(f"{s.description} \n\t=> {len(s)}")
    
    seq = sequences[0]
    df = zhunt(str(seq.seq))
    z = df[df['ZH-Score'] >= 500]
    print(z)
    name_for_bed = "Z_" + str(i) + "_yoelii" + ".bed"
    make_bed(z, name_for_bed, path_to_fasta)
'''
'''
# Этот блок создаёт .bed файлы с Z-SCORE для митохиндриального генома и апикопласта 
# (нашёлся лишь один участок в митохондрии falciparum)
for i in ["mit", "api"]:
    for j in ["/falciparum/", "/knowlesi/", "/yoelii/"]:
        path_to_fasta = "/home/leonid/PycharmProjects/python/minor_2022/project" + j + i + ".fasta"
        with open(path_to_fasta, 'r') as stream:
            sequences = list(SeqIO.parse(stream, format='fasta'))

        # for s in sequences:
        #    print(f"{s.description} \n\t=> {len(s)}")

        seq = sequences[0]
        df = zhunt(str(seq.seq))
        z = df[df['ZH-Score'] >= 500]
        print(z)

        name_for_bed = "Z_" + i + "_" + j[1:-1] + ".bed"
        make_bed(z, name_for_bed, path_to_fasta)
# Обрезаем путь до файлов
for i in range(1, 15):
    for j in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
        path_to_bed = "/home/leonid/PycharmProjects/python/minor_2022/project/final_Z_" + str(i) + "_" + j + ".bed"
        curr_z = pd.read_csv(path_to_bed, names=["chromosome", "start", "end", "chr+start", "Z-SCORE"], sep=' ')
        curr_z["chromosome"] = curr_z["chromosome"].str[55:]
        curr_z["chr+start"] = curr_z["chr+start"].str[55:]
        curr_z.to_csv(path_to_bed, sep=' ', index=False, header=None)
'''

'''
# Этот блок извлекает координаты промоторв из файлов feature table
feature_tables = {}
for i in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
    path_to_ft = "/home/leonid/PycharmProjects/python/minor_2022/project/" + i + ".txt"
    feature_tables[i] = pd.read_csv(path_to_ft, sep='\t', dtype=str)

tss_tables = {}
for i in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
    chromosome_list = []
    start_list = []
    end_list = []
    ID_list = []
    for index, row in feature_tables[i].iterrows():
        if row['# feature'] == 'CDS':
            chromosome_list.append(i + "/" + str(feature_tables[i]["chromosome"][index]))
            if type(feature_tables[i]["product_accession"][index]) is not str:
                ID_list.append("Unknown")
            else:
                ID_list.append(feature_tables[i]["product_accession"][index])
            if row['strand'] == '+':
                if (int(feature_tables[i]["start"][index - 1]) - 300) > 0:
                    start_list.append((int(feature_tables[i]["start"][index - 1]) - 300))
                    end_list.append((int(feature_tables[i]["start"][index - 1]) + 300))
                else:
                    start_list.append((int(feature_tables[i]["start"][index - 1])))
                    end_list.append((int(feature_tables[i]["start"][index - 1]) + 300))
            else:
                if (int(feature_tables[i]["end"][index - 1]) - 300) > 0:
                    start_list.append((int(feature_tables[i]["end"][index - 1]) - 300))
                    end_list.append((int(feature_tables[i]["end"][index - 1]) + 300))
                else:
                    start_list.append((int(feature_tables[i]["end"][index - 1])))
                    end_list.append((int(feature_tables[i]["end"][index - 1]) + 300))
    df_len = len(ID_list)
    tss_tables[i] = pd.DataFrame({'chromosome': chromosome_list, 'start': start_list, 'end': end_list, 'ID': ID_list})
    print(tss_tables)
for i in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
    tss_tables[i].to_csv((i + "_tss" ".bed"), sep='\t', index=False, header=None)
'''

'''
# Этот блок объединяет final .bed файлы
for i in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
    path_to_bed = ("/home/leonid/PycharmProjects/python/minor_2022/project/final_Z_1_" + i + ".bed")
    curr_df = pd.read_csv(path_to_bed, names=["chromosome", "start", "end", "chr+start", "Z-SCORE"], sep=' ')
    for j in range(2, 15):
        path_to_bed = "/home/leonid/PycharmProjects/python/minor_2022/project/final_Z_" + str(j) + "_" + i + ".bed"
        new_df = pd.read_csv(path_to_bed, names=["chromosome", "start", "end", "chr+start", "Z-SCORE"], sep=' ')
        curr_df = curr_df.append(new_df, ignore_index=True)
    name_for_final = "final_Z_" + i + ".bed"
    curr_df.to_csv(name_for_final, sep='\t', index=False, header=None)
'''


'''
# Дополняем таблицу Z-SCORE
z_data = {}
for i in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
    path_to_z_data = "/home/leonid/PycharmProjects/python/minor_2022/project/" + i + "_intersection.bed"
    z_data[i] = pd.read_csv(path_to_z_data, sep='\t', names=["chromosome", "start", "end", "ID", "chromosome_2",
                                                             "start_2", "end_2", "chr+start", "Z-SCORE"])
    z_data[i] = z_data[i].drop(["chromosome_2", "start_2", "end_2", "chr+start"], axis=1).sort_values(by=["Z-SCORE"],
                                                                                                      ascending=False)
ortho_df = pd.read_csv("myproject.proteinortho.tsv", sep='\t')
print(z_data)
print(ortho_df)

score_dict = {"falciparum": [], "gaboni": [], "knowlesi": [], "vivax": [], "yoelii": [], "z_sum": []}
for index, row in ortho_df.iterrows():
    current_score = {"falciparum": 0, "gaboni": 0, "knowlesi": 0, "vivax": 0, "yoelii": 0}
    for i in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
        label_name = i + "_protein.faa"
        for this_ID in row[label_name].split(','):
            number_of_ID = len(row[label_name].split(','))
            if this_ID in list(z_data[i]["ID"]):
                z_index = z_data[i].index[z_data[i]["ID"] == this_ID].tolist()[0]
                current_score[i] += z_data[i]["Z-SCORE"][z_index]
            else:
                if number_of_ID >= 2:
                    number_of_ID -= 1
        score_dict[i].append(current_score[i] / number_of_ID)
    score_dict["z_sum"].append(sum(current_score.values()))
ortho_df = ortho_df.join(pd.DataFrame(score_dict))
ortho_df.sort_values(by=["z_sum"], ascending=False).to_csv("clusters_by_z_sum.tsv", sep='\t', index=False)
'''
'''
# Ищем кластеры
z_data = {}
for i in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
    path_to_z_data = "/home/leonid/PycharmProjects/python/minor_2022/project/" + i + "_intersection.bed"
    z_data[i] = pd.read_csv(path_to_z_data, sep='\t', names=["chromosome", "start", "end", "ID", "chromosome_2",
                                                             "start_2", "end_2", "chr+start", "Z-SCORE"])
    z_data[i] = z_data[i].drop(["chromosome_2", "start_2", "end_2", "chr+start"], axis=1).sort_values(by=["Z-SCORE"],
                                                                                                      ascending=False)
print(z_data)
data = pd.read_csv("clusters_by_z_sum.tsv", sep='\t')
print(data.drop(["# Species", "Genes", "falciparum_protein.faa", "gaboni_protein.faa",
                 "knowlesi_protein.faa", "yoelii_protein.faa"], axis=1).loc[(data["gaboni"] > 0) & (data["yoelii"] > 0)])
print(data.drop(["# Species", "Genes", "Alg.-Conn.", "z_sum"], axis=1).iloc[[1388]])
'''

fig, ((ax1, ax2 ,ax3, ax4)) = plt.subplots(nrows=4, ncols=1,figsize=(12, 6))

features=[
    GraphicFeature(start=20, end=40, strand=+1, color="#ffd700",
                   label="Z-DNA"),
    GraphicFeature(start=20, end=500, strand=+1, color="#ffcccc",
                   label="Gene 1 with a very long name"),

]
record = GraphicRecord(sequence_length=1000, features=features)
record.plot(ax=ax1)

features2=[
    GraphicFeature(start=22, end=45, strand=+1, color="#ffd700",
                   label="Z_DNA"),
           GraphicFeature(start=100, end=110, strand=+1, color="#ffd700",
                   label="Z_DNA"),
    GraphicFeature(start=20, end=500, strand=+1, color="#ffcccc",
                   label="Gene 1 with a very long name"),
]
record2 = GraphicRecord(sequence_length=1000, features=features2)
record2.plot(ax2)
features3=[
    GraphicFeature(start=10, end=600, strand=+1, color="#ffcccc",
                   label="GENE_X"),
    GraphicFeature(start=20, end=38, strand=+1, color="#ffd700",
                   label="z-dna"),
]
record3 = GraphicRecord(sequence_length=600, features=features3)
plt.savefig(("project.pdf"))