#!/usr/bin/python

import os
import sys
import re

re_ns = re.compile('(n+)')

source_name = "fasta"
source_type = "gap"
id_prefix   = "fastagap_"
score, orientation, phase = [ '.', '.', '.' ]

def parse_seq(ofh, seq_name, seq_seq, min_size=1):
    if len(seq_seq) == 0:
        return

    ofh.write("##sequence-region %s 0 %d\n" % (seq_name, len(seq_seq)))

    print "saving chromosome", seq_name, "len", len(seq_seq)
    hit_num = 1
    for m in re_ns.finditer(seq_seq.lower()):
        start_pos = m.start() + 1
        end_pos   = m.end()
        match_seq = m.group()
        diff_pos  = end_pos - start_pos
        match_len = len(match_seq)
        #print seq_name, start_pos, end_pos, diff_pos, match_len, match_seq
        
        if match_len < min_size:
            continue
        
        row_id     = id_prefix + seq_name + '_' + str(hit_num)
        attributes = "ID=%s;Name=%s;length=%d" % ( row_id, row_id, diff_pos )
        
        cols       = [ seq_name, source_name, source_type, start_pos, end_pos, score, orientation, phase, attributes ]
        #print cols
        
        ofh.write("\t".join( [ str(x) for x in cols ] ) + "\n")
        
        hit_num += 1

def main(args):
    infasta  = args[0]
    min_size = 1
    if len(args) > 1:
        min_size = int(args[1])
    
    
    outgff   = infasta + '.gff3'
    
    with open(infasta, 'r') as ifh:
        with open(outgff, 'w') as ofh:
            ofh.write("##gff-version 3\n")
            ofh.write("#infile  : %s\n" % infasta)
            ofh.write("#min_size: %d\n" % min_size)

            seq_name = None
            seq_seq  = ""
            for line in ifh:
                line = line.strip()
            
                if len(line) == 0:
                    continue
                
                if line[0] == ">":
                    if seq_name is not None:
                        parse_seq(ofh, seq_name, seq_seq, min_size=min_size)

                    seq_seq  = ""
                    seq_name = line[1:]
                
                else:
                    seq_seq += line
            
            if seq_name is not None:
                parse_seq(ofh, seq_name, seq_seq, min_size=min_size)
                pass
        
    

if __name__ == '__main__':
    main(sys.argv[1:])