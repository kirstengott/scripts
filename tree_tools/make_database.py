


import os, sys, math, shutil

def find_all(name, path):
    result = []
    for root, dirs, files in os.walk(path):
        if name in files:
            result.append(os.path.join(root, name))
    return result


def main(argv):

    if len(argv) < 2:
        sys.stderr.write('''Usage: %s <input_busco_dir> <outputdir> <cutoff float> <seq_type> 
        input_busco_dir: where the busco input is located
        outputdir:       where to write the output to
        cutoff_float:    an optional decimal percentage value of busco presence to include in database
        seq_type:        one of fna or faa, defaults to fna''' % os.path.basename(sys.argv[0]))
    else:
        dirs          = os.listdir(sys.argv[1])
        if 'fna' in sys.argv:
            seq_type      = 'fna'
        elif 'faa' in sys.argv:
            seq_type      = 'faa'
        else:
            seq_type      = 'fna'

            
        busco_out_dir = sys.argv[2]
        cutoff = 0
        if len(sys.argv) > 3:
            cutoff_in     = float(sys.argv[3])
            cutoff = math.trunc(len(dirs) * cutoff_in)
        
        if os.path.exists(busco_out_dir):
            print(busco_out_dir, 'exists, removing..')
            shutil.rmtree(busco_out_dir)
        os.mkdir(busco_out_dir)




        c = {}
        for d in dirs:
            path = '{}/{}/run_ascomycota_odb10/'.format(sys.argv[1], d)
            single_copies =  path  + "busco_sequences/single_copy_busco_sequences/"
            all_sc       = os.listdir(single_copies)
            all_sc_paths = [os.path.join(single_copies, x) for x in all_sc]
            
            for i in range(0, len(all_sc)):
                x    = all_sc[i]
                path = all_sc_paths[i]
                if seq_type not in x:
                    continue
                else:
                    try:
                        c[x][0] += 1
                        c[x][1].append(path)
                    except:
                        c[x] = [1, [path]]


        



        for busco in c:
            if c[busco][0] >= cutoff:
                outfile = os.path.join(busco_out_dir, busco)
                all_files = c[busco][1]
                for seq_file in all_files:
                    genome = seq_file.split("/")[1]
                    #print(genome)
                    #sys.exit()
                    command = "cat {} | sed -e 's/>.*$/>{}/' >>{}".format(seq_file, genome, outfile)
                    os.system(command)
                        


if __name__ == '__main__':
    main(sys.argv[1:])
