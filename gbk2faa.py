#!/usr/bin/env python3

from Bio import SeqIO

import sys





def gbk2faa(gbk_filename, faa_filename):
    input_handle  = open(gbk_filename, "r")
    output_handle = open(faa_filename, "w")
    for seq_record in SeqIO.parse(input_handle, "genbank") :
        #print("Dealing with GenBank record " + seq_record.id)
        for seq_feature in seq_record.features :
            if seq_feature.type=="CDS" :
                assert len(seq_feature.qualifiers['translation'])==1
                output_handle.write(">{} from {}\n{}\n".format(seq_feature.qualifiers['locus_tag'][0],
                                                               seq_record.name,
                                                               seq_feature.qualifiers['translation'][0]))

    output_handle.close()
    input_handle.close()

    
def main(argv):
    if len(argv) != 2:
        sys.stderr.write("Usage: %s <fastafile> <n_seq_file>\n" % os.path.basename(sys.argv[0]))
        sys.exit(1)
    gbk_filename = sys.argv[1]
    faa_filename = sys.argv[2]
    gbk2faa(gbk_filename = gbk_filename, faa_filename = faa_filename)


if __name__ == '__main__':
    main(sys.argv[1:])
