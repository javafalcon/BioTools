import os
import re
import subprocess
import requests
import json
import pandas as pd
import numpy as np
from Bio import PDB

abbr_AA_dict = {'GLY':'G', 'ALA':'A', 'VAL':'V', 'LEU':'L', 'ILE':'I',
                   'PRO':'P', 'PHE':'F', 'TYR':'Y', 'TRP':'W', 'SER':'S',
                   'THR':'T', 'CYS':'C', 'MET':'M', 'ASN':'N', 'GLN':'Q',
                   'ASP':'D', 'GLU':'E', 'LYS':'K', 'ARG':'R', 'HIS':'H'
                   }

def csv2json():
    df = pd.read_csv('accession_ids.txt', header=None, index_col=False).values.tolist()
    prot_pdb_dict = {}
    for p in df:
        prot_pdb_dict[p[0]] = p[1:]

    with open('accession_ids.json', 'w') as fp:
        json.dump(prot_pdb_dict, fp, indent=4)

class AlphaFoldDB():
    def __init__(self, protpdb_dir=r'E:\Repoes\pytorch\BioTools\protpdb', AlphaFoldPSDUrl='https://alphafold.ebi.ac.uk/files/'):
        with open(r'E:\Repoes\pytorch\BioTools\accession_ids.json', 'r') as fp:
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
        # pf = self.prot_pdb_dict[upid]
        # pdbfile = "{}-model_v{}.pdb".format(pf[2], pf[3])
        pdbfile = 'AF-{}-F1-model_v2.pdb'.format(upid)
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

    def getHomoProtPDBFiles(self, protseq, num=1):
        pdbfiles = []
        # 执行blastp命令，找到序列的同源蛋白。同源蛋白名称输出到output.txt文件
        with open('input.fasta', 'w', encoding='utf-8') as fw:
            fw.write(protseq)
        cmd = 'blastp -db uniprot -query input.fasta -out output.txt -outfmt 7'
        try:
            subprocess.run(cmd, shell=False)
        except:
            return pdbfiles

        # 分析output.txt文件，找到得分最高的同源蛋白，以其序列片段作为最后的输出结果
        pdbfile, start, end = "", 0, 0
        count = 0
        with open('output.txt', 'r', encoding='utf-8') as fr:
            for line in fr:
                if line.startswith('#'):
                    continue
                else:
                    ls = line.replace("\n", "")
                    ls = re.sub(r'\s+', ' ', ls)
                    ls = ls.split(" ")
                    pid = ls[1].split("-")[1] # e.g.: ls[1]='AFDB:AF-Q9Z1X4-F1'
                    self.getPDBFile(pid)
                    start, end = int(ls[8]), int(ls[9])
                    pdbfile = ls[1][5:] + "-model_v2.pdb"
                    pdbfiles.append([pdbfile, start, end])
                    count += 1
                    if count == num:
                        break

        return pdbfiles

class PDBFileParser():
    def __init__(self, pdbfile, protpdb_dir=r'E:\Repoes\pytorch\BioTools\protpdb', start=1, end=None):
        """
        :param protpdb_dir: 保存pdb文件的目录（路径）
        :param pdbfile: pdb文件名称（不含路径）
        :param start: pdb文件中氨基酸序列的起始位置。注：在pdb文件中氨基酸序号是从1开始计数的。
        :param end: pdb文件中氨基酸序列的结束位置（含），如果是None，表示到最后
        """
        self.pdbfile = os.path.join(protpdb_dir, pdbfile)
        self.parser = PDB.PDBParser()
        self.start = start
        self.end = end
        self.structure = self.parser.get_structure('X', self.pdbfile)
        self.model = self.structure[0]
        self.residues = self.model.get_residues()
        seq = ""
        for residue in self.residues:
            seq = seq + abbr_AA_dict[residue.resname]
        self.seq = seq[start-1: end]

    def getAtomsCoordinate(self):
        """
        读取pdb文件，返回每一个原子的坐标
        :param upid: 蛋白质Uniprot accession
        :return: 3-D列表, 列表的长度是蛋白质序列中氨基酸的长度，列表的元素是2-D列表，其长度是每个氨基酸所包含的原子数，每个原子表示为3元向量
        """
        atoms, amino_acid = [], []
        p = '1'
        with open(self.pdbfile, 'r', encoding='utf-8') as fr:
            for line in fr:
                if line.startswith('ATOM'):
                    ls = line.replace("\n", "")
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

    def getAminoAcidsCoordinate(self):
        """
        计算每一个氨基酸的坐标。定义氨基酸坐标为组成氨基酸的所有原子的质心（重心）坐标
        :return: np.ndarray with shape L*3，L是蛋白质序列的长度（氨基酸个数）
        """
        if self.pdbfile is not None:
            atoms = self.getAtomsCoordinate()
            if self.end is None:
                self.end = len(atoms)
            n_aa = self.end-self.start+1
            points = np.ndarray((n_aa, 3), dtype=float)
            for i in range(self.start-1, self.end):
                points[i - self.start + 1] = np.average(np.array(atoms[i]), 0)#每个氨基酸的坐标以其包含所有原子的质心表示
            return points
        else:
            return None

    def secondaryStructure(self):
        structure = self.parser.get_structure("X", self.pdbfile)
        dssp = PDB.DSSP(structure[0], self.pdbfile, dssp='mkdssp')
        ss = []
        if self.end is None:
            self.end = len(dssp)
        for a_key in list(dssp.keys()[self.start-1: self.end]):
            ss.append(dssp[a_key][2])
        return ss


    def dsspData(self):
        structure = self.parser.get_structure("X", self.pdbfile)
        dssp = PDB.DSSP(structure[0], self.pdbfile, dssp='mkdssp')
        res_ls = []
        for key in dssp.keys():
            t = dssp[key]
            #(index,aminoacid):(secondary structure, relative asa, phi, psi)
            res_ls.append((t[2], t[3], t[4], t[5]))
        return res_ls

if __name__ == "__main__":
    # upid = 'Q9JIL3'
    protseq = 'MNALMRLNQLKPGLQYKLISQTGPVHAPIFTMSVEVDGSNFEASGPSKKTAKLHVAVKVLQDMGLPILTKHGKNPVMELNEKRRGLKYELISETGGSHDKRFVMEVEVDGQKFQGAGSNKKVAKAYAALAALEKLFP'
    alphafold = AlphaFoldDB()
    pdbfiles = alphafold.getHomoProtPDBFiles(protseq)
    pdbfile, start, end = pdbfiles[0]
    parser = PDBFileParser(pdbfile, r'E:\Repoes\pytorch\BioTools\protpdb',  start, end)
    ss = parser.secondaryStructure()
    points = parser.getAminoAcidsCoordinate()
