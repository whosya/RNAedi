
in_f = "/home/whosy/workspace/data/三百草原始数据/三白草/to-upload/Sauch_all_gene_delerrogene.gff3"
out_f = "/home/whosy/workspace/data/三百草原始数据/三白草/to-upload/Sauch_all_gene_delerrogene_recds.gff3"

dir = {}
with open(in_f,"r") as f:
    for lines in f:
        if lines.startswith("#"):
            pass
        else:
            type = lines.strip().split("\t")[2]
            if type == "gene":
                gene_key = lines.strip().split("\t")[8].split(";")[0].replace("ID=","")
                if gene_key not in dir:
                    i=0
                dir[gene_key] = [lines]
            elif type == "mRNA":
                gene_key = lines.strip().split("\t")[8].split(";")[1].replace("Parent=","")
                dir[gene_key].append(lines)
            elif type == "exon" :
                gene_key = lines.strip().split("\t")[8].split(";")[1].replace("Parent=","").replace(".1","")
                dir[gene_key].append(lines)
            elif type == "CDS":
                gene_key = lines.strip().split("\t")[8].split(";")[1].replace("Parent=","").replace(".1","")
                i+=1
                dir[gene_key].append(lines.replace("ID=CDS","ID=CDS."+str(i)))
            else:
                pass
with open(out_f,"w") as f2:
    for v in dir.values():
        for line in v:
            f2.write(line)

                




