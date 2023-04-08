import sys
import os
import argparse
import random
import re
import math


def getargs():
    parser = argparse.ArgumentParser(prog="generate plink data")
    parser.add_argument("indi_num", help="number of individual", type=int)
    parser.add_argument("vari_num", help="number of variants", type=int)
    parser.add_argument("--indi_file", help="fam file, if given, will use it instead generate", default=None)
    parser.add_argument("--vari_file", help="bim file, if given, will use it instead of generate", default=None)
    parser.add_argument("-o", "--output", help="name of output file name, default is out", default="out")
    parser.add_argument("--indi_prefix", help="prefix of individual id, default is ind", default="indi")
    parser.add_argument("--vari_prefix", help="prefix of variants id, defualt if var", default="vari")
    parser.add_argument("--na_rate", help="NA rate, defaule is 0.05", default=0.05)
    parser.add_argument("--max_family_size", help="NA rate, defaule is 0.05", default=5, type=int)
    parser.add_argument("--species", help="species", default=None, choices=["human", "mouse"])
    parser.add_argument("--chrom", help="chromosome used when generate variant list, for example '--chrom=1-15,17,19'",
                        type=str, default=None)
    parser.add_argument("--ref_fasta", help="fasta file for reference", type=str, default=None)

    args = parser.parse_args()
    
    return (args.indi_num, args.vari_num, args.indi_file, args.vari_file,
            args.output, args.indi_prefix, args.vari_prefix, args.na_rate,
            args.max_family_size, args.species, args.chrom, args.ref_fasta)


def generate_fam(indi_num, indi_file, output, indi_prefix, max_family_size):
    if (indi_file):
        if not os.path.exists(indi_file):
            print("indi fam file not exists", file=sys.stderr)
            exit(1)
        fout_fam = open(output + ".fam", "w")
        line_count = 0
        for l in open(indi_file):
            line_count += 1
            if line_count > indi_num:
                break
            l = l.rstrip().split()
            print("\t".join(l), file=fout_fam)
        fout_fam.close()

    else:
        fam_size = []
        fam_size_count = 0
        while (True):
            fam_size_rd = random.randint(1, max_family_size)
            if fam_size_count + fam_size_rd >= indi_num:
                fam_size.append(indi_num - fam_size_count)    
                break
            fam_size.append(fam_size_rd)
            fam_size_count += fam_size_rd
            
        fout_fam = open(output + ".fam", "w")
        str_format = "{0:0>" + str(len(str(indi_num))) + "}"
        fam_num = len(fam_size)
        for i in range(fam_num):
            family_id = indi_prefix + str_format.format(i + 1)
            for  j in range(fam_size[i]):
                indi_id = family_id + "_" + str(j + 1)
                father_id = family_id + "_F"
                mother_id = family_id + "_M"
                sex = str(random.randint(1, 2))
                phenotype_value = str(random.randint(1, 2))
                print("\t".join([family_id, indi_id, father_id, mother_id, sex, phenotype_value]), file=fout_fam)
        fout_fam.close()

    return


def parse_chrom_str(chrom_str):
    chroms = chrom_str.split(",")
    dtout = []
    re_sigle_digit = re.compile(r"^\d+$")
    re_digit_range = re.compile(r"^\d+-\d+$")
    for ele in chroms:
        res1 = re_sigle_digit.match(ele)
        res2 = re_digit_range.match(ele)
        if res1:
            dtout.append(int(ele))
        elif res2:
            ele = ele.split("-")
            for i in range(int(ele[0]), int(ele[1]) + 1):
                dtout.append(i)
        else:
            dtout.append(ele)
    
    return dtout
            

def generate_bim(vari_num, vari_file, output, vari_prefix, species, chrom, ref_fasta, max_vari_len):
    if (vari_file):
        if not os.path.exists(vari_file):
            print("variant bim file not exists", file=sys.stderr)
            exit(1)
        fout_bim = open(output + ".bim", "w")
        line_count = 0
        for l in open(vari_file):
            line_count += 1
            if line_count > vari_num:
                break
            l = l.rstrip().split()
            print("\t".join(l), file=fout_bim)
        fout_bim.close()

    else:
        fout_bim = open(output + ".bim", "w")
        if chrom:
            chrom = parse_chrom_str(chrom)
        else:
            chrom = [1]

        if species and not ref_fasta:
            if species == "human":
                chrom_range = {1: 240000000, 2: 240000000, 3: 200000000, 4: 190000000, 5: 180000000,
                    6: 170000000, 7: 150000000, 8: 140000000, 9:13000000, 10: 130000000, 11: 130000000,
                    12: 130000000, 13: 110000000, 14: 100000000, 15: 100000000, 16: 90000000, 17: 80000000,
                    18: 70000000, 19: 60000000, 20: 60000000, 21: 40000000, 22: 40000000, "X": 150000000,
                    "Y": 50000000, "MT": 16000}
            else:
                print("speciese not recgnized", file=sys.sdterr)
                exit(1)
        elif ref_fasta:
            print("reference genome not implimented yet.")
            exit(0)
        else:
            chrom_range = {1: 240000000, 2: 240000000, 3: 200000000, 4: 190000000, 5: 180000000,
                6: 170000000, 7: 150000000, 8: 140000000, 9:13000000, 10: 130000000, 11: 130000000,
                12: 130000000, 13: 110000000, 14: 100000000, 15: 100000000, 16: 90000000, 17: 80000000,
                18: 70000000, 19: 60000000, 20: 60000000, 21: 40000000, 22: 40000000, "X": 150000000,
                "Y": 50000000, "MT": 16000}

        chrom_dic = {}
        chrom_len_sum = 0
        for e in chrom:
            if e in chrom_range:
                chrom_len_sum += chrom_range[e]
                chrom_dic[e] = [chrom_range[e]]
            else:
                print("error, chrom not found in condidate species")
                exit(1)
        vari_sum_2 = 0
        for e in chrom_dic:
            vari_num_per_chrom = math.ceil(chrom_dic[e][0] / chrom_len_sum * vari_num)
            chrom_dic[e].append(vari_num_per_chrom)
            vari_sum_2 += vari_num_per_chrom

        chrom_dic[chrom[-1]][1] = chrom_dic[chrom[-1]][1] - (vari_sum_2 - vari_num) 

        for e in chrom_dic:
            chosed = random.choices(range(1, chrom_dic[e][0] + 1), k=chrom_dic[e][1])        
            chosed.sort()
            chrom_dic[e].append(chosed)
        #   print(chrom_dic)

        chrom.sort()
        for e in chrom:
            i = 0
            for pos in chrom_dic[e][2]:
                i += 1
                allel1 = ""
                allel2 = ""
                while(allel1 == allel2):
                    allel1 = ""
                    allel2 = ""
                    allel1_len = random.choices([1, 2], [90, 10], k=1)
                    allel1_len = allel1_len[0]
                    if allel1_len == 1:
                        allel1 = random.choice(["A", "T", "G", "C"])
                    else:
                        allel1_len = random.choice(range(1, max_vari_len + 1))
                        for i in range(allel1_len):
                            allel1 += random.choice(["A", "T", "G", "C"])
                    allel2_len = random.choices([1, 2], [90, 10], k=1)
                    allel2_len = allel2_len[0]
                    if allel2_len == 1:
                        allel2 = random.choice(["A", "T", "G", "C"])
                    else:
                        allel2_len = random.choice(range(1, max_vari_len + 1))
                        for i in range(allel2_len):
                            allel2 += random.choice(["A", "T", "G", "C"])
                print("\t".join([str(e), vari_prefix + "_" + str(e) + "_" + str(i), "0", str(pos), allel1, allel2]), file=fout_bim)
        fout_bim.close()
    return


def generate_bed(indi_num, vari_num, output, na_rate):
    # 0 homozygous of second allele, 11
    # 1 heterozygouse, 10
    # 2 homozygouse of first allele, 00
    # 3 missing, 01
    dt_per_vari = []
    dt_len = math.ceil(indi_num / 4) * 4
    choice_rate = [1 - na_rate, na_rate]
    fout_bed = open(output + ".bed", "wb")
    for i in range(vari_num):
        dt_per_vari = [0] * dt_len
        for j in range(indi_num):
            if random.choices([1, 2], choice_rate, k=1)[0] == 2:
                dt_per_vari[j] = 3
            else:
                dt_per_vari[j] = random.choice([0, 1, 2])
            
        for j in range(indi_num):
            if dt_per_vari[j] == 0:
                dt_per_vari[j] = '11'
            elif dt_per_vari[j] == 1:
                dt_per_vari[j] = '10'
            elif dt_per_vari[j] == 2:
                dt_per_vari[j] = '00'
            else:
                dt_per_vari[j] = '01'
        dt_new = b''
        for i in list(range(dt_len))[::4]:
            tmp = dt_per_vari[i: i + 4]
            tmp.reverse()
            dt_new += int("".join(tmp), base=2).to_bytes()
        fout_bed.write(dt_new)

    fout_bed.close()        
    return


def main():
    (indi_num, vari_num, indi_file, vari_file, output,
        indi_prefix, vari_prefix, na_rate, max_family_size,
        species, chrom, ref_fasta) = getargs()    

    #generate family data 
    generate_fam(indi_num, indi_file, output, indi_prefix, max_family_size)
    
    #generate variants data
    generate_bim(vari_num, vari_file, output, vari_prefix, species, chrom, ref_fasta, 10)

    #generate bed data
    generate_bed(indi_num, vari_num, output, na_rate)

    return


if __name__ == "__main__":
    main()
