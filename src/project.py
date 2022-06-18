import os
import csv
import sys
import tempfile
import fileinput
import subprocess
from Bio import SeqIO
import seaborn as sns
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
'''
# Этот блок заменяет все символы в .fasta файлах на N (из тех, что не A, C, G, T)
# species = "/yoelii/"
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
# Этот блок работает с статистикой
length_dict = {"falciparum": 0, "vivax": 0, "gaboni": 0, "knowlesi": 0, "yoelii": 0}
for i in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
    for j in range(1, 15):
        path_to_gb = "/home/leonid/PycharmProjects/python/minor_2022/project/" + i + "/" + str(j) + ".gb"
        with open(path_to_gb) as f:
            contents = f.readlines()
            length = int(contents[0].rstrip().split()[2])
            length_dict[i] += length
GC_dict = {"falciparum": 0, "vivax": 0, "gaboni": 0, "knowlesi": 0, "yoelii": 0}
for i in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
    GC_count = 0
    for j in range(1, 15):
        path_to_fasta = "/home/leonid/PycharmProjects/python/minor_2022/project/" + i + "/" + str(j) + ".fasta"
        start_reading = True
        with open(path_to_fasta) as f:
            contents = f.readlines()
            for line in contents:
                if start_reading is True:
                    start_reading = False
                else:
                    for letter in line.rstrip():
                        if letter == "G" or letter == "C" or letter == "g" or letter == "c":
                            GC_count += 1
    GC_dict[i] += (GC_count / length_dict[i]) * 100
print(length_dict)
print(GC_dict)
ann_count = {"falciparum": 0, "vivax": 0, "gaboni": 0, "knowlesi": 0, "yoelii": 0}
ann_length = {"falciparum": 0, "vivax": 0, "gaboni": 0, "knowlesi": 0, "yoelii": 0}
ann_percentage = {"falciparum": 0, "vivax": 0, "gaboni": 0, "knowlesi": 0, "yoelii": 0}
for i in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
    path_to_ft = "/home/leonid/PycharmProjects/python/minor_2022/project/" + i + ".txt"
    ft_df = pd.read_csv(path_to_ft, sep='\t', dtype=str)
    for index, row in ft_df.iterrows():
        if row["# feature"] == "gene" and row["class"] == "protein_coding":
            ann_count[i] += 1
            ann_length[i] += int(row["end"]) - int(row["start"])
for i in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
    ann_percentage[i] += ann_length[i] / length_dict[i]
print(ann_count)
print(ann_length)
print(ann_percentage)
'''
'''
# Этот блок создаёт .bed файлы с Z-SCORE для хромосом
# species = "/yoelii/"
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
# Этот блок извлекает координаты промоторов из файлов feature table
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
# Статистика по Z-ДНК
z_count = {"falciparum": 0, "vivax": 0, "gaboni": 0, "knowlesi": 0, "yoelii": 0}
z_length = {"falciparum": 0, "vivax": 0, "gaboni": 0, "knowlesi": 0, "yoelii": 0}
for i in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
    path_to_bed = ("/home/leonid/PycharmProjects/python/minor_2022/project/final_Z_" + i + ".bed")
    z_df = pd.read_csv(path_to_bed, names=["chromosome", "start", "end", "chr+start", "Z-SCORE"], sep='\t')
    z_count[i] += len(z_df)
    for index, row in z_df.iterrows():
        z_length[i] += int(row["end"]) - int(row["start"])
print(z_count)
print(z_length)
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
# Гистограммы для кластеров
clust = pd.read_csv("clusters_by_z_sum.tsv", sep='\t')
number_of_genoms = []
number_of_proteins = []
for index, row in clust.iterrows():
    curr_gene_num = 0
    curr_prot_num = 0
    if row["falciparum_protein.faa"] != '*':
        curr_gene_num += 1
        curr_prot_num += len(row["falciparum_protein.faa"].split(','))
    if row["vivax_protein.faa"] != '*':
        curr_gene_num += 1
        curr_prot_num += len(row["vivax_protein.faa"].split(','))
    if row["gaboni_protein.faa"] != '*':
        curr_gene_num += 1
        curr_prot_num += len(row["gaboni_protein.faa"].split(','))
    if row["knowlesi_protein.faa"] != '*':
        curr_gene_num += 1
        curr_prot_num += len(row["knowlesi_protein.faa"].split(','))
    if row["yoelii_protein.faa"] != '*':
        curr_gene_num += 1
        curr_prot_num += len(row["yoelii_protein.faa"].split(','))
    number_of_genoms.append(curr_gene_num)
    number_of_proteins.append(curr_prot_num)
clust["number_of_genoms"] = number_of_genoms
clust["number_of_proteins"] = number_of_proteins
# fig, ax = plt.subplots()
# sns.histplot(data=clust, x="number_of_genoms", discrete=True)
# ax.locator_params(axis='x', integer=True)
# plt.tight_layout()
# plt.savefig("clust_by_number_of_genoms.png", dpi=300)

# fig, ax = plt.subplots()
# sns.histplot(data=clust, x="number_of_proteins", discrete=True)
# ax.locator_params(axis='x', integer=True)
# plt.tight_layout()
# plt.savefig("clust_by_number_of_proteins.png", dpi=300)
print(clust)
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
# for i in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
#     print(z_data[i].head(10))

data = pd.read_csv("clusters_by_z_sum.tsv", sep='\t')
print(data.drop(["# Species", "Genes", "falciparum_protein.faa", "gaboni_protein.faa", "knowlesi_protein.faa",
                 "yoelii_protein.faa"], axis=1).loc[((data["falciparum"] > 0) & (data["vivax"] > 0) & (data["gaboni"] > 0) & (data["knowlesi"] > 0)) |
                                                    ((data["falciparum"] > 0) & (data["vivax"] > 0) & (data["gaboni"] > 0) & (data["yoelii"] > 0)) |
                                                    ((data["falciparum"] > 0) & (data["vivax"] > 0) & (data["yoelii"] > 0) & (data["knowlesi"] > 0)) |
                                                    ((data["falciparum"] > 0) & (data["yoelii"] > 0) & (data["gaboni"] > 0) & (data["knowlesi"] > 0)) |
                                                    ((data["yoelii"] > 0) & (data["vivax"] > 0) & (data["gaboni"] > 0) & (data["knowlesi"] > 0))])

# print(data.drop(["# Species", "Genes", "Alg.-Conn.", "z_sum"], axis=1).loc[(data["vivax_protein.faa"] == "EDL44826.1")])
# print(data.drop(["# Species", "Genes", "Alg.-Conn.", "z_sum"], axis=1).iloc[[1388]])
# data.drop(["# Species", "Genes", "Alg.-Conn.", "z_sum"], axis=1).iloc[[21, 33, 45, 258, 339, 394, 552, 566, 607, 780]].to_csv("top_clusters_by_z_sum.tsv", sep='\t', index=False)
# Кластеры с индексами 21, 33, 45, 258, 339, 394, 552, 566, 780
'''
'''
feature_tables = {}
for i in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
    path_to_ft = "/home/leonid/PycharmProjects/python/minor_2022/project/" + i + ".txt"
    feature_tables[i] = pd.read_csv(path_to_ft, sep='\t', dtype=str)

view_clusters = pd.read_csv("top_clusters_by_z_sum.tsv", sep='\t')
for index, row in view_clusters.iterrows():
    name_for_image = "cluster" + str(index + 1) + ".png"
    data_for_features = {"protein_ID": [], "z_dna_start": [], "z_dna_end": [], "length": [],
                         "gene_start": [], "gene_end": [], "strand": [], "function": []}
    for i in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
        name_for_label = i + "_protein.faa"
        data_for_features["protein_ID"].append(row[name_for_label])

        if row[name_for_label] in z_data[i]["ID"].tolist():
            z_data_index = z_data[i].index[z_data[i]["ID"] == row[name_for_label]].tolist()[0]
            data_for_features["z_dna_start"].append(int(z_data[i]["start"][z_data_index]))
            data_for_features["z_dna_end"].append(int(z_data[i]["end"][z_data_index]))
        else:
            data_for_features["z_dna_start"].append(19000000000)
            data_for_features["z_dna_end"].append(57000000000)

        feature_index = feature_tables[i].index[feature_tables[i]["product_accession"] == row[name_for_label]].tolist()[0]
        data_for_features["gene_start"].append(int(feature_tables[i]["start"][feature_index]))
        data_for_features["gene_end"].append(int(feature_tables[i]["end"][feature_index]))
        data_for_features["length"].append(int(feature_tables[i]["end"][feature_index]) - int(feature_tables[i]["start"][feature_index]) + 600)
        data_for_features["function"].append(feature_tables[i]["name"][feature_index])
        if feature_tables[i]["strand"][feature_index] == '+':
            data_for_features["strand"].append(+1)
        else:
            data_for_features["strand"].append(-1)
    for i in range(5):
        switch = data_for_features["gene_start"][i] - 300
        data_for_features["gene_start"][i] -= switch
        data_for_features["gene_end"][i] -= switch
        data_for_features["z_dna_start"][i] -= switch
        data_for_features["z_dna_end"][i] -= switch
    fig, ((ax1, ax2, ax3, ax4, ax5)) = plt.subplots(nrows=5, ncols=1, figsize=(12, 6))

    feature1 = [GraphicFeature(start=data_for_features["z_dna_start"][0], end=data_for_features["z_dna_end"][0],
                               strand=data_for_features["strand"][0], color="#ffd700", label="Z-DNA"),
                GraphicFeature(start=data_for_features["gene_start"][0], end=data_for_features["gene_end"][0],
                               strand=data_for_features["strand"][0], color="#ffcccc", label=data_for_features["function"][0])]
    record1 = GraphicRecord(sequence_length=data_for_features["length"][0], features=feature1)
    record1.plot(ax=ax1)

    feature2 = [GraphicFeature(start=data_for_features["z_dna_start"][1], end=data_for_features["z_dna_end"][1],
                               strand=data_for_features["strand"][1], color="#ffd700", label="Z-DNA"),
                GraphicFeature(start=data_for_features["gene_start"][1], end=data_for_features["gene_end"][1],
                               strand=data_for_features["strand"][1], color="#ffcccc",
                               label=data_for_features["function"][1])]
    record2 = GraphicRecord(sequence_length=data_for_features["length"][1], features=feature2)
    record2.plot(ax=ax2)

    feature3 = [GraphicFeature(start=data_for_features["z_dna_start"][2], end=data_for_features["z_dna_end"][2],
                               strand=data_for_features["strand"][2], color="#ffd700", label="Z-DNA"),
                GraphicFeature(start=data_for_features["gene_start"][2], end=data_for_features["gene_end"][2],
                               strand=data_for_features["strand"][2], color="#ffcccc",
                               label=data_for_features["function"][2])]
    record3 = GraphicRecord(sequence_length=data_for_features["length"][2], features=feature3)
    record3.plot(ax=ax3)

    feature4 = [GraphicFeature(start=data_for_features["z_dna_start"][3], end=data_for_features["z_dna_end"][3],
                               strand=data_for_features["strand"][3], color="#ffd700", label="Z-DNA"),
                GraphicFeature(start=data_for_features["gene_start"][3], end=data_for_features["gene_end"][3],
                               strand=data_for_features["strand"][3], color="#ffcccc",
                               label=data_for_features["function"][3])]
    record4 = GraphicRecord(sequence_length=data_for_features["length"][3], features=feature4)
    record4.plot(ax=ax4)

    feature5 = [GraphicFeature(start=data_for_features["z_dna_start"][4], end=data_for_features["z_dna_end"][4],
                               strand=data_for_features["strand"][4], color="#ffd700", label="Z-DNA"),
                GraphicFeature(start=data_for_features["gene_start"][4], end=data_for_features["gene_end"][4],
                               strand=data_for_features["strand"][4], color="#ffcccc",
                               label=data_for_features["function"][4])]
    record5 = GraphicRecord(sequence_length=data_for_features["length"][4], features=feature5)
    record5.plot(ax=ax5)
    plt.savefig(name_for_image)
'''

# Работа с квадруплексами
'''
# Этот блок объединяет final .bed файлы
for i in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
    path_to_bed = ("/home/leonid/PycharmProjects/python/minor_2022/project/q/" + i + "_1_q.bed")
    curr_df = pd.read_csv(path_to_bed, names=["chromosome", "start", "end", "chr+start", "Q-SCORE", "strand"], sep='\t')
    curr_df = curr_df.drop(["chr+start"], axis=1)
    curr_df["chromosome"] = [(str(i) + "/1")] * len(curr_df)
    for j in range(2, 15):
        path_to_bed = "/home/leonid/PycharmProjects/python/minor_2022/project/q/" + i + "_" + str(j) + "_q.bed"
        new_df = pd.read_csv(path_to_bed, names=["chromosome", "start", "end", "chr+start", "Q-SCORE", "strand"], sep='\t')
        new_df = new_df.drop(["chr+start"], axis=1)
        new_df["chromosome"] = [(str(i) + "/" + str(j))] * len(new_df)
        curr_df = curr_df.append(new_df, ignore_index=True)
    name_for_final = i + "_q.bed"
    curr_df.to_csv(name_for_final, sep='\t', index=False, header=None)
'''

'''
# Дополняем таблицу Q-SCORE
q_data = {}
for i in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
    path_to_q_data = "/home/leonid/PycharmProjects/python/minor_2022/project/" + i + "_q_intersection.bed"
    q_data[i] = pd.read_csv(path_to_q_data, sep='\t', names=["chromosome", "start", "end", "ID", "chromosome_2",
                                                             "start_2", "end_2", "Q-SCORE", "strand"])
    q_data[i] = q_data[i].drop(["chromosome_2", "start_2", "end_2"], axis=1).sort_values(by=["Q-SCORE"], ascending=False)
print(q_data)
ortho_df = pd.read_csv("myproject.proteinortho.tsv", sep='\t')
print(ortho_df)
score_dict = {"falciparum": [], "gaboni": [], "knowlesi": [], "vivax": [], "yoelii": [], "q_sum": []}
for index, row in ortho_df.iterrows():
    current_score = {"falciparum": 0, "gaboni": 0, "knowlesi": 0, "vivax": 0, "yoelii": 0}
    for i in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
        label_name = i + "_protein.faa"
        for this_ID in row[label_name].split(','):
            number_of_ID = len(row[label_name].split(','))
            if this_ID in list(q_data[i]["ID"]):
                q_index = q_data[i].index[q_data[i]["ID"] == this_ID].tolist()[0]
                current_score[i] += q_data[i]["Q-SCORE"][q_index]
            else:
                if number_of_ID >= 2:
                    number_of_ID -= 1
        score_dict[i].append(current_score[i] / number_of_ID)
    score_dict["q_sum"].append(sum(current_score.values()))
ortho_df = ortho_df.join(pd.DataFrame(score_dict))
ortho_df.sort_values(by=["q_sum"], ascending=False).to_csv("clusters_by_q_sum.tsv", sep='\t', index=False)

'''

'''
data = pd.read_csv("clusters_by_q_sum.tsv", sep='\t')
non_zero_number = []
for index, row in data.iterrows():
    curr_non_zero_number = 0
    for i in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
        if row[i] != 0.0:
            curr_non_zero_number += 1
    non_zero_number.append(curr_non_zero_number)
data["non_zero_number"] = non_zero_number
data.drop(["# Species", "Genes", "Alg.-Conn."], axis=1).sort_values(by=["non_zero_number", "q_sum"],
                                                                    ascending=False).head(11).to_csv("top_clusters_by_q_sum.tsv", sep='\t', index=False)
'''


'''
feature_tables = {}
for i in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
    path_to_ft = "/home/leonid/PycharmProjects/python/minor_2022/project/" + i + ".txt"
    feature_tables[i] = pd.read_csv(path_to_ft, sep='\t', dtype=str)

q_data = {}
for i in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
    path_to_q_data = "/home/leonid/PycharmProjects/python/minor_2022/project/" + i + "_q_intersection.bed"
    q_data[i] = pd.read_csv(path_to_q_data, sep='\t', names=["chromosome", "start", "end", "ID", "chromosome_2",
                                                             "start_2", "end_2", "Q-SCORE", "strand"])
    q_data[i] = q_data[i].drop(["chromosome_2", "start_2", "end_2"], axis=1).sort_values(by=["Q-SCORE"], ascending=False)
print(q_data)
print(feature_tables)

view_clusters = pd.read_csv("top_clusters_by_q_sum.tsv", sep='\t')
for index, row in view_clusters.iterrows():
    name_for_image = "cluster" + str(index + 1) + ".png"
    data_for_features = {"protein_ID": [], "q_dna_start": [], "q_dna_end": [], "length": [],
                         "gene_start": [], "gene_end": [], "strand": [], "function": []}
    for i in ["falciparum", "vivax", "gaboni", "knowlesi", "yoelii"]:
        name_for_label = i + "_protein.faa"
        data_for_features["protein_ID"].append(row[name_for_label])

        if row[name_for_label] in q_data[i]["ID"].tolist():
            q_data_index = q_data[i].index[q_data[i]["ID"] == row[name_for_label]].tolist()[0]
            data_for_features["q_dna_start"].append(int(q_data[i]["start"][q_data_index]))
            data_for_features["q_dna_end"].append(int(q_data[i]["end"][q_data_index]))
        else:
            data_for_features["q_dna_start"].append(19000000000)
            data_for_features["q_dna_end"].append(57000000000)

        feature_index = feature_tables[i].index[feature_tables[i]["product_accession"] == row[name_for_label]].tolist()[0]
        data_for_features["gene_start"].append(int(feature_tables[i]["start"][feature_index]))
        data_for_features["gene_end"].append(int(feature_tables[i]["end"][feature_index]))
        data_for_features["length"].append(int(feature_tables[i]["end"][feature_index]) - int(feature_tables[i]["start"][feature_index]) + 600)
        data_for_features["function"].append(feature_tables[i]["name"][feature_index])
        if feature_tables[i]["strand"][feature_index] == '+':
            data_for_features["strand"].append(+1)
        else:
            data_for_features["strand"].append(-1)
    for i in range(5):
        switch = data_for_features["gene_start"][i] - 300
        data_for_features["gene_start"][i] -= switch
        data_for_features["gene_end"][i] -= switch
        data_for_features["q_dna_start"][i] -= switch
        data_for_features["q_dna_end"][i] -= switch
    fig, ((ax1, ax2, ax3, ax4, ax5)) = plt.subplots(nrows=5, ncols=1)

    feature1 = [GraphicFeature(start=data_for_features["q_dna_start"][0], end=data_for_features["q_dna_end"][0],
                               strand=data_for_features["strand"][0], color="#ffd700", label="q-DNA"),
                GraphicFeature(start=data_for_features["gene_start"][0], end=data_for_features["gene_end"][0],
                               strand=data_for_features["strand"][0], color="#ffcccc", label=data_for_features["function"][0])]
    record1 = GraphicRecord(sequence_length=data_for_features["length"][0], features=feature1)
    record1.plot(ax=ax1)

    feature2 = [GraphicFeature(start=data_for_features["q_dna_start"][1], end=data_for_features["q_dna_end"][1],
                               strand=data_for_features["strand"][1], color="#ffd700", label="q-DNA"),
                GraphicFeature(start=data_for_features["gene_start"][1], end=data_for_features["gene_end"][1],
                               strand=data_for_features["strand"][1], color="#ffcccc",
                               label=data_for_features["function"][1])]
    record2 = GraphicRecord(sequence_length=data_for_features["length"][1], features=feature2)
    record2.plot(ax=ax2)

    feature3 = [GraphicFeature(start=data_for_features["q_dna_start"][2], end=data_for_features["q_dna_end"][2],
                               strand=data_for_features["strand"][2], color="#ffd700", label="q-DNA"),
                GraphicFeature(start=data_for_features["gene_start"][2], end=data_for_features["gene_end"][2],
                               strand=data_for_features["strand"][2], color="#ffcccc",
                               label=data_for_features["function"][2])]
    record3 = GraphicRecord(sequence_length=data_for_features["length"][2], features=feature3)
    record3.plot(ax=ax3)

    feature4 = [GraphicFeature(start=data_for_features["q_dna_start"][3], end=data_for_features["q_dna_end"][3],
                               strand=data_for_features["strand"][3], color="#ffd700", label="q-DNA"),
                GraphicFeature(start=data_for_features["gene_start"][3], end=data_for_features["gene_end"][3],
                               strand=data_for_features["strand"][3], color="#ffcccc",
                               label=data_for_features["function"][3])]
    record4 = GraphicRecord(sequence_length=data_for_features["length"][3], features=feature4)
    record4.plot(ax=ax4)

    feature5 = [GraphicFeature(start=data_for_features["q_dna_start"][4], end=data_for_features["q_dna_end"][4],
                               strand=data_for_features["strand"][4], color="#ffd700", label="q-DNA"),
                GraphicFeature(start=data_for_features["gene_start"][4], end=data_for_features["gene_end"][4],
                               strand=data_for_features["strand"][4], color="#ffcccc",
                               label=data_for_features["function"][4])]
    record5 = GraphicRecord(sequence_length=data_for_features["length"][4], features=feature5)
    record5.plot(ax=ax5)
    plt.savefig(name_for_image, dpi=600)
'''
