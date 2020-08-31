import sys, os
from shutil import copyfile


def main(argv):
    if len(argv) < 1:
        sys.stderr.write('''Usage: %s <ncbi_genome_metadata>
        <ncbi_genome_metadata>\t ncbi genome metadata from ncbi genome download \n''' % os.path.basename(sys.argv[0]))
        sys.exit(1)

    header = 0
    taxa = {}

    rank = {'Chromosome': 1,
            'Complete Genome' : 2,
            'Contig' : 4,
            'Scaffold' : 3}

    outfile = 'ncbi_genome_metadata_subset.txt'
    
    with open(sys.argv[1], 'r') as fh:
        for line in fh:
            if header == 0:
                header = 1
                continue
            line = line.strip().split("\t")
            keep = [line[13], line[22], line[0]]
            key = " ".join(line[9].split(" ")[0:2])
            rank_key = rank[line[13]]
            if key not in taxa:
                taxa[key] = {}
                taxa[key][rank_key] = [keep]
            else:
                try:
                    taxa[key][rank_key].append(keep)
                except:
                    taxa[key][rank_key] = [keep]




    genomes_out = []
    genera = {}
    for i in taxa:
        genus = i.split()[0]
        if not genus in genera:
            genera[genus] = 1
        else:
            genera[genus] += 1
        if genus == 'Trichoderma':
            pass
        elif genera[genus] > 2:
            continue
        ranks = [x for x in taxa[i].keys()]
        ranks.sort()
        genome_keep = []
        for rank in ranks:
            genomes = taxa[i][rank]
            for genome in genomes:
                if len(genome_keep) < 1: ## only keep one genome per species
                    genome_keep = genome
                    genome_keep.append(i)
                else:
                    break
        genomes_out.append(genome_keep)

    if not os.path.exists('genomes'):
        os.mkdir('genomes')

    out_f = open(outfile, 'w')
    for i in genomes_out:
        dst = os.path.join('genomes', i[2] + '.fna.gz')
        #print(i)
        copyfile(i[1], dst)
        out_f.write(i[1] + "\t" + i[3] + "\n")
    out_f.close()

        
if __name__ == '__main__':
    main(sys.argv[1:])
    
