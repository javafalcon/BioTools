import os
import re
import subprocess
import requests
import json
import pandas as pd
import numpy as np
from Bio.Blast.Applications import NcbiblastpCommandline

def csv2json():
    df = pd.read_csv('accession_ids.txt', header=None, index_col=False).values.tolist()
    prot_pdb_dict = {}
    for p in df:
        prot_pdb_dict[p[0]] = p[1:]

    with open('accession_ids.json', 'w') as fp:
        json.dump(prot_pdb_dict, fp, indent=4)

class AlphaFoldDB():
    def __init__(self, protpdb_dir='protpdb', AlphaFoldPSDUrl='https://alphafold.ebi.ac.uk/files/'):
        with open('accession_ids.json', 'r') as fp:
            self.prot_pdb_dict = json.load(fp)
        self.protpdb_dir = protpdb_dir # 本地的文件夹名，用来保存下载的pdb文件
        self.url = AlphaFoldPSDUrl # AlphaFold Protein Structure Database下载pdb文件网址的前缀
        self.protpdb_file_list = os.listdir(self.protpdb_dir) # 在本地文件夹protpdb_dir下所有pdb文件名列表

    def getPDBFile(self, upid):
        """
        从AlphaFold Protein Structure Database下载蛋白质的三维结构预测结果，预测结果被保留在pdb文件；
        如果某个蛋白质（Uniprot accession）的pdb已经下载到本地，则直接返回文件名
        :param upid: 蛋白质Uniprot accession
        :return: PDB文件名或None（如果获取网页资源失败）
        """
        pf = self.prot_pdb_dict[upid]
        pdbfile = "{}-model_v{}.pdb".format(pf[2], pf[3])

        if pdbfile not in self.protpdb_file_list:
            try:
                r = requests.get(self.url + pdbfile, timeout=30)  # eg. 'AF-X1WFM8-F1-model_v2.pdb'
                r.raise_for_status()
                r.encoding = 'utf-8'
                with open(os.path.join(self.protpdb_dir, pdbfile), 'wb') as fw:
                    fw.write(r.content)
            except:
                return None

        return pdbfile

    def getAtomsCoordinate(self,pdbfile):
        """
        读取pdb文件，返回每一个原子的坐标
        :param upid: 蛋白质Uniprot accession
        :return: 3-D列表, 列表的长度是蛋白质序列中氨基酸的长度，列表的元素是2-D列表，其长度是每个氨基酸所包含的原子数，每个原子表示为3元向量
        """
        atoms, amino_acid= [], []
        p = '1'
        with open(os.path.join(self.protpdb_dir, pdbfile),'r', encoding='utf-8') as fr:
            for line in fr:
                if line.startswith('ATOM'):
                    ls = line.replace("\n","")
                    ls = re.sub(r'\s+', ' ', ls)
                    ls = ls.split(" ")
                    if ls[5] == p:
                        amino_acid.append([float(a) for a in ls[6:9]])
                    else:
                        atoms.append(amino_acid)
                        amino_acid = []
                        p = ls[5]
                        amino_acid.append([float(a) for a in ls[6:9]])
            atoms.append(amino_acid)
        return atoms

    def getAminoAcidsCoordinate(self, upid):
        """
        计算每一个氨基酸的坐标。定义氨基酸坐标为组成氨基酸的所有原子的质心（重心）坐标
        :param upid: 蛋白质序列的Uniport accession
        :return: np.ndarray with shape L*3，L是蛋白质序列的长度（氨基酸个数）
        """
        pdbfile = self.getPDBFile(upid)

        if pdbfile is not None:
            atoms = self.getAtomsCoordinate(pdbfile)
            n_aa = len(atoms)
            points = np.ndarray((n_aa, 3), dtype=float)
            for i in range(n_aa):
                points[i] = np.average(np.array(atoms[i]), 0)#每个氨基酸的坐标以其包含所有原子的质心表示
            return points
        else:
            return None

    def getAminoAcidsCoordinateFromHomoprot(self, protseq):
        # 执行blastp命令，找到序列的同源蛋白。同源蛋白名称输出到output.txt文件
        with open('input.fasta', 'w', encoding='utf-8') as fw:
            fw.write(protseq)
        cmd = 'blastp -db uniprot -query input.fasta -out output.txt -outfmt 7'
        try:
            subprocess.run(cmd, shell=False)
        except:
            return ""

        # 分析output.txt文件，找到得分最高的同源蛋白，以其序列片段作为最后的输出结果
        pdbfile, start, end = "", 0, 0
        with open('output.txt', 'r', encoding='utf-8') as fr:
            for line in fr:
                if line.startswith('#'):
                    continue
                else:
                    ls = line.replace("\n","")
                    ls = re.sub(r'\s+', ' ', ls)
                    ls = ls.split(" ")
                    pdbfile, start, end = ls[1], int(ls[8]), int(ls[9])
                    pdbfile = pdbfile[5:]+"-model_v2.pdb"
                    break
        # 读取pdb文件获取原子坐标
        atoms = self.getAtomsCoordinate(pdbfile)
        points = np.ndarray((end-start+1, 3), dtype=float)
        for i in range(start, end+1):
            points[i-start] = np.average(np.array(atoms[i]), 0)

        return points

if __name__ == "__main__":
    # upid = 'Q9JIL3'
    protseq = 'MNALMRLNQLKPGLQYKLISQTGPVHAPIFTMSVEVDGSNFEASGPSKKTAKLHVAVKVLQDMGL\
    PILTKHGKNPVMELNEKRRGLKYELISETGGSHDKRFVMEVEVDGQKFQGAGSNKKVAKAYAALAALEKLFP'
    alphafold = AlphaFoldDB()
    points = alphafold.getAminoAcidsCoordinateFromHomoprot(protseq)
